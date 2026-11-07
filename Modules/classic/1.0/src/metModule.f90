!> \file
!> Performs disaggregation of input meteorological forcing arrays to the model physics timestep
!! @author V. Arora, J. Melton
module metDisaggModule

  use modelStateDrivers, only : metInputTimeStep, metTime, metFss, metFdl, metPre, metTa, metQa, metUv, metPres
  use generalUtils,      only : closeEnough
  use classicParams,     only : pi

  implicit none

  public :: disaggMet
  public :: makebig
  public :: stepInterpolation
  public :: linearInterpolation
  public :: precipDistribution
  public :: diurnalDistribution
  public :: zenithAngles
  public :: daylightIndices
  public :: distributeDiurnally
  public :: timeZone
  public :: timeShift

  integer :: numberPhysInMet      !< Number of physics timesteps that fit into one MET input file timestep
  integer :: numberMetInputinDay  !< Number of MET input file timesteps that fit into one day
  integer :: tslongCount          !< Number of physics timesteps in the original met arrays
  integer :: shortSteps           !< Number of physics timesteps in one day
  
  logical, parameter :: allLocalTime = .false.  !< If your gridded meteorology is relative to Greenwich then set this 
                                                !! to false. If the meteorology is all in local time set it to true. To
                                                !! determine which, look at your shortwave radiation. As you move through time 
                                                !! do you see the sun move across longitudes (allLocalTime = .false.) or across 
                                                !! latitudes (allLocalTime = .true.). It seems reanalysis will generally be false 
                                                !! while climate model outputs are generally true.

contains

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup metdisaggmodule_disaggMet
  !! @{
  !> Main subroutine to disaggregate input meteorology to that of the physics timestep

  subroutine disaggMet (longitude, latitude) ! longitude, latitude

    use classicParams, only : delt

    implicit none

    real, intent(in) :: longitude, latitude  ! in degrees
    integer          :: vcount, vcountPlus

    !> First check that we should be doing the disaggregation. If we are already
    !! at the needed timestep (delt) then we can return to the main.
    if (closeEnough(metInputTimeStep,delt,0.001)) then
      return
    end if

    !> Determine how many physics timesteps fit into a metInputTimeStep
    numberPhysInMet = nint(metInputTimeStep / delt)

    !> Determine how many metInputTimeStep fit into a day
    numberMetInputinDay = nint(86400./metInputTimeStep)


    !> Determine the number of physics timesteps that fit in a day
    shortSteps = numberPhysInMet * numberMetInputinDay

    !> Compute the timestep count (both original timestep and the needed one)
    tslongCount = size(metTime)              ! original timesteps

    ! needed number of timesteps at the physics timestep
    vcount = tslongCount * numberPhysInMet

    !> Add two extra days to the vcount since we need them for interpolation of the
    !! first and last days of our time period
    vcountPlus = vcount + 2 * numberPhysInMet * numberMetInputinDay

    !> Expand original MET data to physics timestep. This increases the array size
    !! and puts the old values within the larger array. It also duplicates the first
    !! and last days for the interpolation. This also expands the time array (metTime)
    call makebig(vcount, vcountPlus)

    !> Perform the interpolations that are either step (for longwave), linear (for
    !! for air temperature, specific humidity, wind, pressure), randomly distributed
    !! over a number of wet timesteps (precip), or diurnally distributed (shortwave)
    call stepInterpolation(metFdl)
    call linearInterpolation(metTa)
    call linearInterpolation(metQa)
    call linearInterpolation(metUv)
    call linearInterpolation(metPres)
    call precipDistribution(metPre)
    call diurnalDistribution(metFss, latitude)

    !> Adjust the met arrays for the timezone relative to Greenwich and also
    !! trim off the added two days.
    call timeShift(timeZone(longitude), vcount, vcountPlus)

  end subroutine disaggMet
  !! @}

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup metdisaggmodule_makebig
  !! @{
  !> Expands the meteorological variable arrays so they can accomodate the amount of
  !! timesteps on the physics timestep. Also copies the first day and last day values
  !! into the array to give a startday -1 and endday + 1 values for the interpolations.
  !! While we are at it we also expand the time array
  subroutine makebig(vcount, vcountPlus)

    use classicParams, only : delt

    implicit none

    integer, intent(in) :: vcount
    integer, intent(in) :: vcountPlus
    real, allocatable, dimension(:) :: tmpFss, tmpFdl, tmpPre, tmpTa, tmpQa, tmpUv, tmpPres, tmpTime
    integer :: i, j, timePlusTwoDays

    !> Transfer the present contents of the met*** arrays to temp ones.
    call move_alloc(metFss, tmpFss)
    call move_alloc(metFdl, tmpFdl)
    call move_alloc(metPre, tmpPre)
    call move_alloc(metTa, tmpTa)
    call move_alloc(metQa, tmpQa)
    call move_alloc(metUv, tmpUv)
    call move_alloc(metPres, tmpPres)
    call move_alloc(metTime, tmpTime)

    !> Allocate the met*** arrays in preparation to fill with old contents
    !! the allocated size is now including the met on the physics timestep.
    !! The shortwave and precipitation arrays don't require the extra two
    !! padding days.
    allocate(metFss(vcountPlus), metFdl(vcountPlus), metPre(vcountPlus), metTa(vcountPlus), &
             metQa(vcountPlus), metUv(vcountPlus), metPres(vcountPlus), metTime(vcount))

    ! initialize to zero so they are not filled in with random value by the compiler
    metFdl = 0. ; metFss = 0. ; metPre = 0. ; metTa = 0. ; metQa = 0. ; metUv = 0. ; metPres = 0.

    !> The remainder of the meteorological variables need to have a couple extra
    !! days pinned on,one at the start and one at the end.
    timePlusTwoDays = tslongCount + 2 * numberMetInputinDay

    !> Now fill the newly expanded arrays with the values you have in the position
    !! of their timesteps
    do j = 1,timePlusTwoDays
      if (j < numberMetInputinDay + 1) then ! Then copy the first day into that time
        i = j
      else if (j > (timePlusTwoDays - numberMetInputinDay)) then ! Copy the last day into that time
        i = j - (2 * numberMetInputinDay)
      else ! use the value for that time
        i = j - numberMetInputinDay
      end if
      metFdl ((j - 1) * numberPhysInMet + 1) = tmpFdl(i)
      metFss ((j - 1) * numberPhysInMet + 1) = tmpFss(i)
      metPre ((j - 1) * numberPhysInMet + 1) = tmpPre(i)
      metTa  ((j - 1) * numberPhysInMet + 1) = tmpTa(i)
      metQa  ((j - 1) * numberPhysInMet + 1) = tmpQa(i)
      metUv  ((j - 1) * numberPhysInMet + 1) = tmpUv(i)
      metPres((j - 1) * numberPhysInMet + 1) = tmpPres(i)
    end do

    !> Now expand the time dimension and fill in the values.
    do j = 1,tslongCount ! Don't need the buffer days so use tslongCount
      metTime((j - 1) * numberPhysInMet + 1) = tmpTime(j)
      do i = 1,numberPhysInMet - 1
        metTime((j - 1) * numberPhysInMet + 1 + i) = tmpTime(j) + i * delt / 86400.
      end do
    end do

    deallocate(tmpFss,tmpFdl,tmpPre,tmpTa,tmpQa,tmpUv,tmpPres,tmpTime)

  end subroutine makebig
  !! @}

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup metdisaggmodule_stepInterpolation
  !! @{
  !> Step interpolation applies the same value across all physics timesteps
  subroutine stepInterpolation (var)

    implicit none

    real, intent(inout) :: var(:)
    integer            :: i, j, countr

    countr = size(var)

    do i = 1,countr / numberPhysInMet
      do j = 2,numberPhysInMet
        var((i - 1) * numberPhysInMet + j) = var((i - 1) * numberPhysInMet + 1)
      end do
    end do

  end subroutine stepInterpolation
  !! @}

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup metdisaggmodule_linearInterpolation
  !! @{
  !> This method splits the dataset into intervals and performs linear interpolation of each individual interval
  subroutine linearInterpolation (var)

    implicit none

    real, intent(inout) :: var(:)
    integer             :: i, j, shortIntervals, start, endpt, countr
    real                :: t

    countr = size(var)

    ! Determine how many short intervals we'll process
    shortIntervals = countr / numberPhysInMet
    ! For every interval except the last one, do:
    do i = 1,shortIntervals - 1
      ! Copy the first value of the next interval to the last position of the current interval
      var(i * numberPhysInMet) = var(i * numberPhysInMet + 1)
    end do
    ! For the last interval, copy the first value of the interval to the last value
    var(countr) = var(countr - numberPhysInMet)

    ! For every interval
    do i = 1,shortIntervals
      ! Determine the start and end indices of the interval
      do j = 0,numberPhysInMet - 1
        start = (i - 1) * numberPhysInMet + 1
        endpt = i * numberPhysInMet
        ! Compute the interval offset
        t = real(j) / numberPhysInMet
        ! Compute the interpolated values based on the difference betweent the start and end as interval offset
        var(start + j) = var(start) + (var(endpt) - var(start)) * t
      end do
    end do

  end subroutine linearInterpolation
  !! @}

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup metdisaggmodule_precipDistribution
  !! @{
  !> Precipitation distribution occurs randomly, but conservatively, over the number of wet timesteps.
  subroutine precipDistribution (var)

    use classicParams, only : delt, zero

    implicit none

    real, intent(inout)                :: var(:)
    integer                            :: i, j, k, T, start, endpt, countr, wetpds, attempts
    real                               :: temp, startpre
    real, allocatable, dimension(:)    :: random
    integer, allocatable, dimension(:) :: sort_ind
    real, allocatable, dimension(:)    :: incomingPre, tmpvar
    logical                            :: needDistrib

    ! Set some initial conditions
    needDistrib = .true.
    attempts = 1

    allocate(random(numberPhysInMet), sort_ind(numberPhysInMet))

    countr = size(var)

    ! Loop through the metInputTimeStep timesteps (commonly 6 hr)

    ! Adjust the precip from mm/s to mm/6h (or mm/3h). The relationship below is
    ! derived for mm/6h originally.
    var = var * metInputTimeStep

    ! Save the incoming precip for balance check later
    allocate(incomingPre(countr), tmpvar(countr))
    incomingPre = var

    do while (needDistrib)
      do i = 1,countr/numberPhysInMet

        start = (i - 1) * numberPhysInMet + 1
        endpt = i * numberPhysInMet

        !> If precipitation for this metInputTimeStep is greater than 0, use formula
        !! to produce number of wet half hours

        if (var(start) > 0.) then

          startpre = var(start) * 1.E6 ! bump this up so we don't have problems due to truncation later.
          ! the numbers can be small so get truncated which leads to no balance.

          ! We expect precipitation for this relation to be in mm/6h
          wetpds = nint( &
                   max( &
                   min( &
                   abs(2.6 * log10(6.93 * var(start) ) ) &
                   , real(numberPhysInMet)) &
                   , 1.0) &
                   )

          ! Create an array of random numbers
          call random_number(random(:))

          ! Create an array of integers from 1 to numberPhysInMet
          do k = 1,numberPhysInMet
            sort_ind(k) = k
          end do

          ! Now use the random numbers to create a list of randomly sorted
          ! integers that will be used as indices later. Here the integers are
          ! sorted such that the highest integer :: ends up in the index of the largest
          ! random number. The sorting continues until integer :: 1 is in the index position
          ! of the smallest random number generated.

          do k = 1,numberPhysInMet
            do j = k,numberPhysInMet
              if (random(k) < random(j)) then
                temp = random(k)
                random(k) = random(j)
                random(j) = temp
                T = sort_ind(k)
                sort_ind(k) = sort_ind(j)
                sort_ind(j) = T
              end if
            end do
          end do

          ! Set all values in random to 0 in preparation for reassignment
          random(:) = 0.

          ! Produces random list of 1s and 0s
          do j = 1,wetpds
            call random_number(random(sort_ind(j)))
          end do

          random(:) = random(:) / sum(random(:))

          ! Check if random now sums to 1, if not readjust.
          if (sum(random(:)) /= 1.0) then
            random(:) = random(:) / sum(random(:))
          end if

          ! Disperse precipitiation randomly across the random indice
          k = 1
          do j = start,endpt
            ! Assign this time period its precipitation
            tmpvar(j) = random(k) * startpre
            k = k + 1
          end do

        else ! No precip, move on.
          wetpds = 0
          tmpvar(start:endpt) = 0.
        end if

      end do

      tmpvar = tmpvar * 1E-6 ! bump back down so it back in expected units.

      ! Balance check that we have conserved our precip
      if ((sum(tmpvar) - sum(incomingPre)) > 1.0e-5) then
        if (attempts > 3) then
          print * ,'Warning: In precipDistribution,precip is not being conserved',sum(var),sum(incomingPre)
          call errorHandler('metModule', - 1)
          return ! this is needed here as it ensures we don't get caught in a loop where it keeps failing but not moving on.
        else ! retry, could have just been a bad draw.
          needDistrib = .true.
          attempts = attempts + 1
        end if
      else
        ! All is well, finish up.
        needDistrib = .false.
      end if

    end do ! needDistrib

    ! So we now have the amount of precip in each physics timestep (mm/delt)
    ! we now need to then convert back to mm/s
    var = tmpvar / delt

    deallocate(random, sort_ind, incomingPre, tmpvar)

  end subroutine precipDistribution
  !! @}

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup metdisaggmodule_diurnalDistribution
  !! @{
  !> Diurnal distribution over the entire timespan
  subroutine diurnalDistribution (shortWave, latitude)

    use classicParams, only : delt

    implicit none

    real, intent(inout)                         :: shortWave(:)
    real, intent(in)                            :: latitude
    integer                                     :: i, d, start, endpt, midday
    integer, allocatable                        :: daylightIndYear(:,:)
    real, allocatable                           :: zenithAngYear(:,:)
    real, allocatable                           :: zenithNoon(:)
    real                                        :: latRad
    real                                        :: swMean
    real, allocatable, dimension(:)             :: timesteps

    midday = shortSteps / 2 + 1

    allocate(timesteps(shortSteps))

    do i = 1,shortSteps       ! Creates the time steps in physics timesteps (in seconds)
      timesteps(i) = real(i - 1) * delt
    end do

    allocate(daylightIndYear(365,shortSteps),zenithAngYear(365,shortSteps), &
                    zenithNoon(365))

    latRad = latitude * pi / 180.00

    ! Assume we are going to run the model for at least one full year so just find the
    ! zenith angles and daylight timesteps for one year and store the values.
    do d = 1,365 ! FLAG not set up for leap yet !
      zenithAngYear(d,:) = zenithAngles(timesteps,d,latRad,365,shortSteps)
      ! Find the zenith angle at noon
      zenithNoon(d) = zenithAngYear(d,midday)
      daylightIndYear(d,:) = daylightIndices(zenithAngYear(d,:),shortSteps)
    end do

    d = 365 ! We start on day 365 since we have the extra day added on the start
    do i = 1,size(shortWave) / shortSteps      ! Once, every day
      start = (i - 1) * shortSteps + 1
      endpt = i * shortSteps

      swMean = sum(shortWave(start:endpt)) / numberMetInputinDay ! Finds the mean of this day's values

      shortWave(start:endpt) = distributeDiurnally(zenithAngYear(d,:),daylightIndYear(d,:),zenithNoon(d),swMean,shortSteps)
      d = d + 1
      if (d > 365) d = 1
    end do

    deallocate(daylightIndYear,zenithAngYear,zenithNoon)

  end subroutine diurnalDistribution
  !! @}

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup metdisaggmodule_zenithAngles
  !! @{
  !> Find zenith angles depending on day of year and latitude
  function zenithAngles (timesteps, day, latRad, lastDOY, countr)

    implicit none

    integer, intent(in)          :: day
    real, intent(in)          :: timesteps(:)
    integer, intent(in)          :: countr
    integer, intent(in)          :: lastDOY
    real, intent(in)          :: latRad     !< In radians
    real, allocatable               :: zenithAngles(:)
    integer                         :: i
    real                            :: radHourAngle, degHourAngle, psi, dec

    real, parameter                 :: secondsInHour = 3600.
    real, parameter, dimension(4)   :: a = (/0.006918, - 0.399912, - 0.006758, - 0.002697/)
    real, parameter, dimension(4)   :: b = (/0.0,0.070257,0.000907,0.001480/)
    real, parameter, dimension(4)   :: n = (/0., 1., 2., 3./)

    allocate(zenithAngles(countr))

    psi = 2. * pi * real(day - 1) / real(lastDOY)
    dec = sum((a * cos(n * psi)) + (b * sin(n * psi)))

    ! Find the hour angle, convert it to radians then find the zenith angle(s).
    do i = 1,countr
      degHourAngle = 15. * (12. - (timesteps(i) / secondsInHour))
      radHourAngle = (degHourAngle / 360.) * 2. * pi
      zenithAngles(i) = acos((sin(latRad) * sin(dec)) + &
                        (cos(latRad) * cos(dec) * cos(radHourAngle)))
    end do

  end function zenithAngles
  !! @}

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup metdisaggmodule_daylightIndices
  !! @{
  !> Find day lengths depending on zenith angles
  function daylightIndices (zenithAngles, countr)

    use classicParams, only : delt

    implicit none

    integer                         :: i, daylightCount
    integer, intent(in)             :: countr
    real, intent(in)                :: zenithAngles(:)
    integer, allocatable            :: daylightIndices(:)
    real                            :: zenithCos
    real                            :: dayLength

    allocate(daylightIndices(countr))
    daylightCount = 0
    do i = 1,countr
      zenithCos = cos(zenithAngles(i))
      if (zenithCos >= 0) then
        daylightCount = daylightCount + 1
        daylightIndices(i) = i
      else
        daylightIndices(i) = - 1
      end if
    end do
    ! dayLength = daylightCount * delt ! purely diagnostic.
    ! print*, 'daylength in seconds', dayLength

  end function daylightIndices
  !! @}

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup metdisaggmodule_distributeDiurnally
  !! @{
  !> Determines correction for input values based on daylight indices and zenith angles
  function distributeDiurnally (zenithAngles, daylightIndices, zenithNoon, swMean, countr)

    implicit none

    integer, intent(in)          :: countr
    real, intent(in)          :: zenithAngles(:), zenithNoon, swMean
    integer, intent(in)          :: daylightIndices(:)
    real, allocatable         :: distributeDiurnally(:)
    real, allocatable         :: diurnalDistrib(:)
    integer                         :: i, daylightIndex
    real                            :: vsum, correction, diurnalMean

    allocate(diurnalDistrib(countr), distributeDiurnally(countr))

    vsum = 0.0
    do i = 1,countr
      daylightIndex = daylightIndices(i)              ! check with daylight_indices to see if it there is daylight at the time,
      if (daylightIndex > 0) then                     ! if there is find the value for s
        diurnalDistrib(i) = swMean * pi / 2. * cos((zenithAngles(daylightIndex) - zenithNoon) &
                            / (pi / 2. - zenithNoon) * pi / 2.)
      else ! else set it to 0
        diurnalDistrib(i) = 0
      end if
      vsum = vsum + diurnalDistrib(i)
    end do

    diurnalMean = vsum / countr

    if (diurnalMean > 0.) then
      correction = swMean / diurnalMean
    else
      correction = 0.
    end if
    distributeDiurnally = diurnalDistrib * correction

    deallocate(diurnalDistrib)

  end function distributeDiurnally
  !! @}

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup metdisaggmodule_timeZone
  !! @{
  !> Find the local timezone relative to Greenwich.
  real function timeZone (longitude)

    use classicParams, only : delt

    implicit none

    real, intent(in)      :: longitude

    !> Each timezone is 15 deg of longitude per hour. So convert our physics
    !! timestep into a minutes equivalent.
    timeZone = real(nint(longitude / 15. * 60.)) / (delt / 60.)

    if (longitude >= 180. .or. longitude < 0) then ! The longitudes run from 0 to 360 or from -180 to 180
      timeZone = shortSteps - timezone
    else ! Because of how cshift works,this needs to be negative.
      timeZone = - timeZone
    end if

  end function timeZone
  !! @}

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup metdisaggmodule_timeShift
  !! @{
  !> Perform a circular shift on the met arrays to account for the timezone the present cell is in.
  !! Also trims the met arrays to remove the two extra days added for the interpolations
  subroutine timeShift (timeZoneOffset, vcount, vcountPlus)

    implicit none

    real, intent(in)                :: timeZoneOffset
    integer, intent(in)             :: vcount
    integer, intent(in)             :: vcountPlus
    real, allocatable, dimension(:) :: tmpFss, tmpFdl, tmpPre, tmpTa, tmpQa, tmpUv, tmpPres

    !> Do a circular shift to the arrays (except Fss !) to adjust for the timezone
    !! but only if you have the meteorology relative to Greenwich.
    if (.not. allLocalTime) then 
      metFdl = cshift(metFdl,int(timeZoneOffset))
      metPre = cshift(metPre,int(timeZoneOffset))
      metTa = cshift(metTa,int(timeZoneOffset))
      metQa = cshift(metQa,int(timeZoneOffset))
      metUv = cshift(metUv,int(timeZoneOffset))
      metPres = cshift(metPres,int(timeZoneOffset))
    end if 

    !> Transfer the present contents of the met*** arrays to temp ones.
    call move_alloc(metFss,tmpFss)
    call move_alloc(metFdl,tmpFdl)
    call move_alloc(metPre,tmpPre)
    call move_alloc(metTa,tmpTa)
    call move_alloc(metQa,tmpQa)
    call move_alloc(metUv,tmpUv)
    call move_alloc(metPres,tmpPres)

    allocate(metFss(vcount),metFdl(vcount),metPre(vcount),metTa(vcount), &
                 metQa(vcount),metUv(vcount),metPres(vcount))

    metFss = tmpFss(shortSteps:vcountPlus - shortSteps)
    metFdl = tmpFdl(shortSteps:vcountPlus - shortSteps)
    metPre = tmpPre(shortSteps:vcountPlus - shortSteps)
    metTa = tmpTa(shortSteps:vcountPlus - shortSteps)
    metQa = tmpQa(shortSteps:vcountPlus - shortSteps)
    metUv = tmpUv(shortSteps:vcountPlus - shortSteps)
    metPres = tmpPres(shortSteps:vcountPlus - shortSteps)

    deallocate(tmpFss,tmpFdl,tmpPre,tmpTa,tmpQa,tmpUv,tmpPres)

  end subroutine timeShift
  !! @}
  !> \namespace metdisaggmodule
  !> Performs disaggregation of input meteorological forcing arrays to the model
  !! physics timestep (commonly half-hour)

end module metDisaggModule
