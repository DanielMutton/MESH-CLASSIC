!> \file
!> Central module for all general utilities
!!
module generalUtilsGAT

  implicit none

  public :: abandonCell
  public :: findDaylength
  public :: calcEsat
  public :: findCloudinessGAT
  public :: findLeapYears
  public :: findPermafrostVarsGAT
  public :: parseTimeStamp
  public :: closeEnough
  public :: initRandomSeed
  public :: checksumCalc

  logical :: run_model           !< Simple logical switch to either keep run going or finish

contains

  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_abandonCell
  !! @{
  !> Used to stop running a model grid cell. For errors that need to be caught early in a run,
  !! the fortran intrinsic 'stop' is preferred but for errors later in a run or simple fails on
  !! single grid cells, abandonCell is best since it allows the netcdf files to continue to
  !! written to and won't disrupt the MPI processes (as stop does)
  !! @author Joe Melton and Ed Wisernig
  !!
  subroutine abandonCell (errmsg)

    use classStateVars, only : class_rot

    implicit none

    character( * ), intent(in), optional :: errmsg

    real, pointer, dimension(:) :: DLONROW !<
    real, pointer, dimension(:) :: DLATROW !<

    DLATROW => class_rot%DLATROW
    DLONROW => class_rot%DLONROW

    run_model = .false.
    if (present(errmsg)) then
      print * ,errmsg
!MESH      print * ,'exiting cell: ',DLONROW,DLATROW
    end if
!MESH    print * ,'died on',DLONROW,DLATROW
    return

  end subroutine abandonCell
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_findDaylength
  !! @{
  !> Calculate the daylength based on the latitude and day of year
  !! @author Joe Melton
  !!
  real function findDaylength (solday, radl)

    ! Joe Melton Dec 18 2015 (taken from phenlogy.f)

    use classicParams, only : pi

    implicit none

    real, intent(in) :: solday  ! day of year
    real, intent(in) :: radl    ! latitude
    real :: theta               ! temp var
    real :: decli               ! temp var
    real :: term                ! temp var

    theta = 0.2163108 + 2.0 * atan(0.9671396 * tan(0.0086 * (solday - 186.0)))
    decli = asin(0.39795 * cos(theta))      ! declination ! note I see that CLASS does this also but with different formula...
    term = (sin(radl) * sin(decli))  /(cos(radl) * cos(decli))
    term = max( - 1.0,min(term,1.0))
    findDaylength = 24.0 - (24.0/pi) * acos(term)

  end function findDaylength
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_calcEsat
  !! @{
  !> Calculate the saturated vapour pressure in Pa. Based upon 
  !! the parameterization of Emanuel, 1994 \cite Emanuel1994-dt. 
  !! @author Joe Melton
  !!
  real function calcEsat(ta)

    use classicParams, only : TFREZ

    implicit none

    real, intent(in) :: ta  ! air/canopy temperature (K)
    
    if (ta >= tfrez) then
      calcEsat = exp(53.67957 - 6743.769/ta - 4.8451 * log(TA)) * 100. !100 converts from hPa to Pa.
    else !
      calcEsat = exp(23.33086 - 6111.72784/ta + 0.15215 * log(TA)) * 100. !100 converts from hPa to Pa.
    end if
      
  end function calcEsat
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_findLeapYears
  !! @{
  !> Check if this year is a leap year
  !! @author Joe Melton
  !!
  subroutine findLeapYears (iyear, leapnow, lastDOY)

    use classicParams, only : monthend, mmday, monthdays

    implicit none

    logical, intent(out)   :: leapnow
    integer, intent(in)    :: iyear
    integer, intent(inout) :: lastDOY

    if (mod(iyear,4) /= 0) then ! it is a common year
      leapnow = .false.
    else if (mod(iyear,100) /= 0) then ! it is a leap year
      leapnow = .true.
    else if (mod(iyear,400) /= 0) then ! it is a common year
      leapnow = .false.
    else ! it is a leap year
      leapnow = .true.
    end if

    if (leapnow) then
      lastDOY = 366
    else
      lastDOY = 365
    end if

    ! We do not check the MET files to make sure the incoming MET is in fact
    ! 366 days if leapnow. You must verify this in your own input files. Later
    ! in the code it will fail and print an error message to screen warning you
    ! that your file is not correct.
    if (leapnow) then ! adjust the calendar and set the error check.
      monthdays = (/ 31,29,31,30,31,30,31,31,30,31,30,31 /)
      monthend = (/ 0,31,60,91,121,152,182,213,244,274,305,335,366 /)
      mmday = (/ 16,46,76,107,137,168,198,229,260,290,321,351 /)

    else
      monthdays = [ 31,28,31,30,31,30,31,31,30,31,30,31 ] !< days in each month
      monthend  = [ 0,31,59,90,120,151,181,212,243,273,304,334,365 ] !< calender day at end of each month
      mmday     = [ 16,46,75,106,136,167,197,228,259,289,320,350 ] !< mid-month day
    end if

  end subroutine findLeapYears
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_findCloudiness
  !! @{
  !> The cosine of the solar zenith angle COSZ is calculated from the day of
  !> the year, the hour, the minute and the latitude using basic radiation geometry,
  !> and (avoiding vanishingly small numbers) is assigned to CSZROW.  The fractional
  !> cloud cover FCLOROW is commonly not available so a rough estimate is
  !> obtained by setting it to 1 when precipitation is occurring, and to the fraction
  !> of incoming diffuse radiation XDIFFUS otherwise (assumed to be 1 when the sun
  !> is at the horizon, and 0.10 when it is at the zenith).
  !! @author Diana Verseghy
  !!
  subroutine findCloudinessGAT (radjgat, cszgat, pregat, xdiffusgat, fclogat, ilg, il1, il2, imin, ihour, iday, lastDOY)

    use classicParams, only : pi
!MESH    use classStateVars, only : class_rot

    implicit none

!MESH    integer, intent(in) :: nltest
    integer, intent(in) :: ilg
    integer, intent(in) :: il1
    integer, intent(in) :: il2
    integer, intent(in) :: imin
    integer, intent(in) :: ihour
    integer, intent(in) :: iday
    integer, intent(in) :: lastDOY

    integer :: i
    real :: day
    real :: decl    !< Declination
    real :: hour
    real, dimension(ilg) :: coszgat    !< Cosine of the zenith angle

    real, dimension(ilg) :: radjgat !< Latitude of grid cell (positive north of equator) [rad]
    real, dimension(ilg) :: cszgat  !< Cosine of solar zenith angle [ ]
    real, dimension(ilg) :: pregat  !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg) :: xdiffusgat !< Fraction of diffused radiation
    real, dimension(ilg) :: fclogat !< Fractional cloud cover [ ]

!MESH    RADJROW => class_rot%RADJROW
!MESH    CSZROW => class_rot%CSZROW
!MESH    PREROW => class_rot%PREROW
!MESH    XDIFFUS => class_rot%XDIFFUS
!MESH    FCLOROW => class_rot%FCLOROW

    day = real(iday) + (real(ihour) + real(imin)/60.)/24.
    decl = sin(2. * pi * (284. + day)/real(lastDOY)) * 23.45 * pi/180.
    hour = (real(ihour) + real(imin)/60.) * pi/12. - pi

    do i = il1,il2

      coszgat(i) = sin(radjgat(i)) * sin(decl) + cos(radjgat(i)) * cos(decl) * cos(hour)

      cszgat(i) = sign(max(abs(coszgat(i)),1.0e-3),coszgat(i))
      if (pregat(i) > 0.) then
        xdiffusgat(i) = 1.0
      else
        xdiffusgat(i) = max(0.0,min(1.0 - 0.9 * coszgat(i),1.0))
      end if
      fclogat(i) = xdiffusgat(i)
    end do

  end subroutine findCloudinessGAT
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_parseTimeStamp
  !! @{
  !> Parses a time stamp in the expected form "day as %Y%m%d.%f"
  !! Returns an array with 1) year, 2) month, 3) day, 4) fraction of day
  !! 5) day of year
  !! @author Joe Melton, Ed Wisernig
  !!
  function parseTimeStamp (timeStamp)

    use classicParams, only : monthdays
	use iso_fortran_env
	
    implicit none

    real, dimension(5) :: parseTimeStamp
    !real, intent(in)  :: timeStamp
    !real              :: date, moment
    integer            :: intdate, day, month, year, totdays, t
	real, intent(in)   :: timeStamp
	character(250)   :: timeStampTemp
	real(real64)       :: dayfraction
	real(real64)       :: temp

	!print *, timeStamp
    !date = floor(timeStamp) ! remove the part days
    !parseTimeStamp(4) = timeStamp - date ! save the part days
    !intdate = int(date)
    !day = mod(intdate,100);   intdate = intdate / 100
    !month = mod(intdate,100); intdate = intdate / 100
    !year = intdate
    !parseTimeStamp(1) = real(year)
    !parseTimeStamp(2) = real(month)
    !parseTimeStamp(3) = real(day)
    !totdays = 0
    !if (month > 1) then
    !  do t = 1,month - 1
    !    totdays = totdays + monthdays(t)
    !  end do
    !end if
    !parseTimeStamp(5) = real(totdays + day)
	
	WRITE(timeStampTemp, *) timeStamp
	
	timeStampTemp = timeStampTemp(3:)
	!print *, timeStampTemp(2:5)
	read( timeStampTemp(2:5), * ) year
	!print *, timeStampTemp(6:7)
	read( timeStampTemp(6:7), * ) month
	!print *, timeStampTemp(8:9)
	read( timeStampTemp(8:9), * ) temp
	!print *, timeStampTemp(10:)
	read( timeStampTemp(10:), * ) dayfraction
	day = temp
	
	
	parseTimeStamp(1) = real(year)
    parseTimeStamp(2) = real(month)
    parseTimeStamp(3) = real(day)
	parseTimeStamp(4) = dayfraction
	totdays = 0
	if (month > 1) then
      do t = 1,month - 1
        totdays = totdays + monthdays(t)
      end do
    end if
	parseTimeStamp(5) = real(totdays + day)

  end function parseTimeStamp
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_findPermafrostVars
  !! @{
  !> Finds the active layer depth and depth to the frozen water table.
  !! @author Joe Melton
  !!
  subroutine findPermafrostVarsGAT ( &
    ftablegat, actlyrgat, maxAnnualActLyrgat, actLyrThisYrgat, &
    tbargat, thlqgat, thicgat, isndgat, dlzwgat, delz, thmgat, dlatgat, &
    ilg, il1, il2, ig, iday)

    use classicParams, only : ignd, tfrez, eftime, efoldfact
    use classStateVars, only : class_rot, class_gat

    implicit none

!MESH    integer, intent(in) :: nmtest
!MESH    integer, intent(in) :: nltest
    integer, intent(in) :: ilg
    integer, intent(in) :: il1
    integer, intent(in) :: il2
    integer, intent(in) :: ig
    integer, intent(in) :: iday
    real, dimension(ilg)  :: ftablegat      !< Depth to frozen water table (m)
    real, dimension(ilg)  :: actlyrgat      !< Active layer depth (m)
    real, dimension(ilg) :: dlatgat        !< Latitude (degrees)
    real, dimension(ilg,ig) :: tbargat    !< Temperature of soil layers [K]
    integer, dimension(ilg,ig) :: isndgat !< Sand content flag, used to delineate non-soils.
    real, dimension(ilg,ig) :: thicgat    !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, dimension(ilg,ig) :: thlqgat    !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, dimension(ilg,ig) :: dlzwgat    !< Permeable thickness of soil layer [m]
    real, dimension(ilg) :: delz           !< Overall thickness of soil layer [m]
    real, dimension(ilg,ig) :: thmgat     !< Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$
    real, dimension(ilg) :: actLyrThisYrgat !< Annual active layer depth maximum starting from summer solstice for the present year (m)
    real, dimension(ilg) :: maxAnnualActLyrgat !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).

    integer :: i, j

!MESH    ftable  => class_rot%ftable
!MESH    actlyr  => class_rot%actlyr
!MESH    maxAnnualActLyr => class_rot%maxAnnualActLyrROT
!MESH    actLyrThisYr => class_rot%actLyrThisYrROT
!MESH    tbarrot => class_rot%tbarrot
!MESH    thlqrot => class_rot%thlqrot
!MESH    thicrot => class_rot%thicrot
!MESH    isndrot => class_rot%isndrot
!MESH    dlzwrot => class_rot%dlzwrot
!MESH    delz    => class_gat%delz
!MESH    thmrot => class_rot%thmrot
!MESH    dlatrow => class_rot%dlatrow
    !---

    actlyrgat = 0.0
    ftablegat = 0.0
    do j = 1,ig
      do i = il1,il2
!MESH        do m = 1,nmtest
          if (abs(tbargat(i,j) - tfrez) < 0.0001) then
            if (isndgat(i,j) > - 3) then
              actlyrgat(i) = actlyrgat(i) + (thlqgat(i,j) / (thlqgat(i,j) &
                            + thicgat(i,j))) * dlzwgat(i,j)
              ftablegat(i) = ftablegat(i) + (thicgat(i,j) / (thlqgat(i,j) &
                            + thicgat(i,j) - thmgat(i,j))) * dlzwgat(i,j)
              ! else if (isndgat(1,j)==-3) then
              !    actlyr=actlyr+delz(j)
              !    ftable=ftable+delz(j)
            end if
          else if (tbargat(i,j) > tfrez) then
            actlyrgat(i) = actlyrgat(i) + delz(j)
            ftablegat(i) = ftablegat(i) + delz(j)
          end if

!MESH        end do
      end do
    end do

    do i = il1,il2
!MESH      do m = 1,nmtest

        ! Once a year we adjust the maximum annual active layer depth
        ! in an e-folding sense with the
        ! present maxmimum active layer depth for the year ending on the summer
        ! solstice. The maximum values are used in bio2str
        ! to ensure roots are not placed into frozen soil layers.

        if ((dlatgat(i) > 0. .and. iday == 355) & ! Boreal :: winter solstice.
            .or. (dlatgat(i) < 0. .and. iday == 172)) then  ! Austral winter solstice.
          maxAnnualActLyrgat(i) = maxAnnualActLyrgat(i) * efoldfact &
                                 + actLyrThisYrgat(i) * (1.0 - efoldfact)
          actLyrThisYrgat(i) = 0.0
        else
          ! Compare the present active layer depth against actLyrThisYr
          actLyrThisYrgat(i) = max(actLyrThisYrgat(i),actlyrgat(i))
        end if
!MESH      end do
    end do

  end subroutine findPermafrostVarsGAT
  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup generalutils_findPermafrostVars
  !! @{
  !> As real :: numbers are not precise, this is a simple way to compare two reals
  !! @author Joe Melton
  !!
  logical function closeEnough (num1, num2, error)

    implicit none

    real, intent(in)     :: num1, num2
    real, intent(in)     :: error
    if (abs(num1 - num2) < error) then
      closeEnough = .true.
    else
      closeEnough = .false.
    end if
  end function closeEnough
  !! @}
  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup generalutils_initRandomSeed
  !! @{
  !! This subroutine sets a repeatable seed for the random number generator.
  !! @author J. Melton
  !!
  subroutine initRandomSeed

    implicit none

    ! NOTE: this subroutine will eventually be replaced by "call random_init()" which
    ! is an intrinsic in the  Fortran 2018 standard.
    integer :: n
    integer, allocatable :: seed(:)

    call random_seed(size = n)

    allocate(seed(n))

    seed = 589389089

    call random_seed(put = seed)

  end subroutine initRandomSeed
  !! @}

  !> \ingroup generalutils_checksumCalc
  !! @{
  !! This subroutine takes the lonIndex and latIndex of a cell, accesses many
  !! attributes of the cell after the run, and creates a checksum from those attributes.
  !! This checksum is written to a .csv file in the output directory, which is then
  !! compared against the checksum from a previous run.
  !! @author M. Fortier
  !!
  subroutine checksumCalc (lonIndex, latIndex)

    use ctemStateVars,  only : c_switch, vrot
    use classStateVars, only : class_rot
    use classicParams,  only : icc, nmos, ignd, icp1, modelpft, iccp2, TFREZ

    implicit none

    ! arguments
    integer, intent(in) :: lonIndex       !< Index of the longitude of this cell
    integer, intent(in) :: latIndex       !< Index of the latitude of this cell
    ! pointers:
    real, pointer, dimension(:,:,:) :: FCANROT
    real, pointer, dimension(:,:)   :: FAREROT
    real, pointer, dimension(:,:,:) :: TBARROT
    real, pointer, dimension(:,:,:) :: THLQROT
    real, pointer, dimension(:,:,:) :: THICROT
    real, pointer, dimension(:,:)   :: TCANROT
    real, pointer, dimension(:,:)   :: TSNOROT
    real, pointer, dimension(:,:)   :: TPNDROT
    real, pointer, dimension(:,:)   :: ZPNDROT
    real, pointer, dimension(:,:)   :: RCANROT
    real, pointer, dimension(:,:)   :: SCANROT
    real, pointer, dimension(:,:)   :: SNOROT
    real, pointer, dimension(:,:)   :: ALBSROT
    real, pointer, dimension(:,:)   :: RHOSROT
    real, pointer, dimension(:,:)   :: GROROT

    logical, pointer :: ctem_on
    logical, pointer :: PFTCompetition
    logical, pointer :: lnduseon
    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:,:) :: gleafmasrow
    real, pointer, dimension(:,:,:) :: bleafmasrow
    real, pointer, dimension(:,:,:) :: stemmassrow
    real, pointer, dimension(:,:,:) :: rootmassrow
    real, pointer, dimension(:,:) :: twarmm
    real, pointer, dimension(:,:) :: tcoldm
    real, pointer, dimension(:,:) :: gdd5
    real, pointer, dimension(:,:) :: aridity           !< aridity index, ratio of potential evaporation to precipitation
    real, pointer, dimension(:,:) :: srplsmon          !< number of months in a year with surplus water i.e.precipitation more than potential evaporation
    real, pointer, dimension(:,:) :: defctmon          !< number of months in a year with water deficit i.e.precipitation less than potential evaporation
    real, pointer, dimension(:,:) :: anndefct          !< annual water deficit (mm)
    real, pointer, dimension(:,:) :: annsrpls          !< annual water surplus (mm)
    real, pointer, dimension(:,:) :: annpcp            !< annual precipitation (mm)
    real, pointer, dimension(:,:) :: dry_season_length !< length of dry season (months)
    real, pointer, dimension(:,:,:) :: litrmassrow
    real, pointer, dimension(:,:,:) :: soilcmasrow
    integer, pointer, dimension(:,:,:) :: lfstatusrow
    integer, pointer, dimension(:,:,:) :: pandaysrow
    real, pointer, dimension(:,:) :: Cmossmas          !< C in moss biomass, \f$kg C/m^2\f$
    real, pointer, dimension(:,:) :: litrmsmoss        !< moss litter mass, \f$kg C/m^2\f$
    real, pointer, dimension(:,:) :: dmoss             !< depth of living moss (m)


    ! local variables for this subroutine
    integer                             :: checksum, k !, PFTCover_size
    real, allocatable, dimension(:)     :: PFTCover, Soil, Canopy_real, Snow, CPools, Peatlands, CompetClimate !< Checksum values to track
    integer, allocatable, dimension(:)  :: Canopy_int
    character(len = 10)                 :: lonchar, latchar, checksumchar !< string representations
    character(len = 500)                :: buffer, filename

    ! point pointers:
    ctem_on           => c_switch%ctem_on
    PFTCompetition    => c_switch%PFTCompetition
    lnduseon          => c_switch%lnduseon
    fcancmxrow        => vrot%fcancmx
    gleafmasrow       => vrot%gleafmas
    bleafmasrow       => vrot%bleafmas
    stemmassrow       => vrot%stemmass
    rootmassrow       => vrot%rootmass
    twarmm            => vrot%twarmm
    tcoldm            => vrot%tcoldm
    gdd5              => vrot%gdd5
    aridity           => vrot%aridity
    srplsmon          => vrot%srplsmon
    defctmon          => vrot%defctmon
    anndefct          => vrot%anndefct
    annsrpls          => vrot%annsrpls
    annpcp            => vrot%annpcp
    dry_season_length => vrot%dry_season_length
    litrmassrow       => vrot%litrmass
    soilcmasrow       => vrot%soilcmas
    lfstatusrow       => vrot%lfstatus
    pandaysrow        => vrot%pandays
    Cmossmas          => vrot%Cmossmas
    litrmsmoss        => vrot%litrmsmoss
    dmoss             => vrot%dmoss
    FCANROT           => class_rot%FCANROT
    FAREROT           => class_rot%FAREROT
    TBARROT           => class_rot%TBARROT
    THLQROT           => class_rot%THLQROT
    THICROT           => class_rot%THICROT
    TCANROT           => class_rot%TCANROT
    TSNOROT           => class_rot%TSNOROT
    TPNDROT           => class_rot%TPNDROT
    ZPNDROT           => class_rot%ZPNDROT
    RCANROT           => class_rot%RCANROT
    SCANROT           => class_rot%SCANROT
    SNOROT            => class_rot%SNOROT
    ALBSROT           => class_rot%ALBSROT
    RHOSROT           => class_rot%RHOSROT
    GROROT            => class_rot%GROROT


    write(lonchar, '(I4)')lonIndex
    write(latchar, '(I4)')latIndex
    ! generate the proper, formatted filename
    !ignoreLint(3) (it messes with the file path)
    write(filename, "(A,A,A,'_',A,A)") TRIM(adjustl(c_switch%output_directory)), '/checksums/', &
                                       TRIM(adjustl(lonchar)), TRIM(adjustl(latchar)), '.csv'
    open(unit = 500, file = TRIM(adjustl(filename)))


    checksum = 0
    allocate(PFTCover(SIZE(FAREROT) + SIZE(FCANROT) + SIZE(GROROT) + SIZE(fcancmxrow))) !< allocate the array
    PFTCover = (/ pack(FAREROT, .true.),pack(FCANROT, .true.),pack(GROROT, .true.),pack(fcancmxrow, .true.) /)
    do k = 1,SIZE(PFTCover)
      checksum = checksum + bitcount(PFTCover(k))
    end do
    write(checksumchar, '(I4)')checksum   !< transfer to a string variable
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "PFTCover",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Soil(SIZE(THLQROT) + SIZE(THICROT) + SIZE(TPNDROT) + SIZE(ZPNDROT)))
    Soil = (/ pack(THLQROT, .true.),pack(THICROT, .true.),pack(TPNDROT, .true.),pack(ZPNDROT, .true.) /)
    do k = 1,SIZE(Soil)
      checksum = checksum + bitcount(Soil(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "Soil",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Canopy_real(SIZE(TCANROT) + SIZE(RCANROT) + SIZE(SCANROT)))
    Canopy_real = (/ pack(TCANROT, .true.),pack(RCANROT, .true.),pack(SCANROT, .true.) /)
    allocate(Canopy_int(SIZE(lfstatusrow) + SIZE(pandaysrow)))
    Canopy_int = (/ pack(lfstatusrow, .true.),pack(pandaysrow, .true.) /)
    do k = 1,SIZE(Canopy_real)
      checksum = checksum + bitcount(Canopy_real(k))
    end do
    do k = 1,SIZE(Canopy_int)
      checksum = checksum + bitcount_int(Canopy_int(k))
    end do

    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "Canopy",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Snow(SIZE(TSNOROT) + SIZE(SNOROT) +  SIZE(ALBSROT) + SIZE(RHOSROT)))
    Snow = (/ pack(TSNOROT, .true.),pack(SNOROT, .true.),pack(ALBSROT, .true.),pack(RHOSROT, .true.) /)
    do k = 1,SIZE(Snow)
      checksum = checksum + bitcount(Snow(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "Snow",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(CPools(SIZE(gleafmasrow) + SIZE(bleafmasrow) + SIZE(stemmassrow) + SIZE(rootmassrow) + SIZE(litrmassrow) + SIZE(soilcmasrow)))
    CPools = (/ pack(gleafmasrow, .true.),pack(bleafmasrow, .true.),pack(stemmassrow, .true.),pack(rootmassrow, .true.),pack(litrmassrow, .true.),pack(soilcmasrow, .true.) /)
    do k = 1,SIZE(Cpools)
      checksum = checksum + bitcount(Cpools(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "Cpools",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(Peatlands(SIZE(Cmossmas) + SIZE(litrmsmoss) + SIZE (dmoss)))
    Peatlands = (/ pack(Cmossmas, .true.),pack(litrmsmoss, .true.),pack(dmoss, .true.) /)
    do k = 1,SIZE(Peatlands)
      checksum = checksum + bitcount(Peatlands(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "Peatlands",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    checksum = 0
    allocate(CompetClimate(SIZE(twarmm) + SIZE(tcoldm) + SIZE(gdd5) + SIZE(aridity) + SIZE(srplsmon) + SIZE(defctmon) + SIZE(anndefct) + SIZE(annsrpls) + SIZE(annpcp) + SIZE(dry_season_length)))
    CompetClimate = (/ pack(twarmm, .true.),pack(tcoldm, .true.),pack(gdd5, .true.),pack(aridity, .true.),pack(srplsmon, .true.),pack(defctmon, .true.),pack(anndefct, .true.),pack(annsrpls, .true.),pack(annpcp, .true.),pack(dry_season_length, .true.) /)
    do k = 1,SIZE(CompetClimate)
      checksum = checksum + bitcount(CompetClimate(k))
    end do
    write(checksumchar, '(I4)')checksum
    write(buffer,"(A,',',A,',',A,',',A)") TRIM(adjustl(lonchar)),TRIM(adjustl(latchar)), &
                                          "CompetClimate",TRIM(adjustl(checksumchar))
    write(500, "(A)") TRIM(adjustl(buffer))


    close(500)

  end subroutine checksumCalc
  !! @}

  !> \ingroup generalutils_bitcount
  !! @{
  !! This function generates the bitcount of the specified variable
  !! @author M. Fortier
  !!
  integer function bitcount (scalar)

    implicit none
    real    :: scalar     !< Scalar to be bit-counted
    integer(SELECTED_REAL_KIND(15,307)) :: scalar_int !< Integer with memory representation of 'scalar'

    integer :: bit
    bitcount = 0
    scalar_int = TRANSFER(scalar,scalar_int)
    do bit = 0,BIT_SIZE(scalar_int) - 1
      if (BTEST(scalar_int,bit) ) bitcount = bitcount + 1
    end do

  end function bitcount
  !! @}


  !> \ingroup generalutils_bitcount_int
  !! @{
  !! This function generates the bitcount of the specified variable, from an integer
  !! @author M. Fortier
  !!
  integer function bitcount_int (scalar)

    implicit none
    integer    :: scalar     !< Scalar to be bit-counted
    integer    :: bit

    bitcount_int = 0
    do bit = 0,BIT_SIZE(scalar) - 1
      if (BTEST(scalar,bit) ) bitcount_int = bitcount_int + 1
    end do

  end function bitcount_int
  !! @}
  ! -----------------------------------------------------------------------------------------------

  !> \namespace generalUtilsGAT
  !> Central module for all general utilities
  !!
  !! The checksum subroutine computes content-based checksums of several groups of variables after
  !! a run has completed. These checksums are used when making non-logical changes to
  !! the code-base. This is accomplished by comparing the computed checksums, to the
  !! checksums of an identical run (same parameters, met files, etc.) before the changes
  !! were made. This is a modified form of regression testing specific to non-logical
  !! software changes.
  !!
  !! Checksums are not an infallible means of ensuring no logical changes, as two numbers
  !! may have different binary representations with the same number of flipped bits. However,
  !! with the number of variables we are checking, it is extremely unlikely to render a false-
  !! positive. If we take a data item of bit-length \f$n\f$, having \f$b\f$ flipped bits in
  !! its representation, the number of same-length data items with the same checksum can
  !! be expressed through the binary coefficient
  !!
  !! \f[ {n\choose b} \f]
  !!
  !! If we assume the worst case where \f$b=\frac{n}{2}\f$, then we are left with
  !!
  !! \f[ {n\choose n/2} \f]
  !!
  !! Dividing by the total number of possible values for an n-digit binary number,
  !! we find the probability of a false positive checksum to be, in the worst case:
  !!
  !! \f[ \large \frac{{n\choose n/2}}{2^n} \f]
  !!
  !! Due to the cascading effects of any logical changes to the model, all 7 groups
  !! of variables (6 if peatlands is diabled) will be affected. The smallest of these
  !! groups consists of 3 arrays of numbers for which checksums are computed. With the
  !! conservative assumption that each array contains only a single 32-bit number, the
  !! probability of a false positive (in the absolute worst case) becomes:
  !!
  !! \f[ \large \left({\frac{96\choose 48}{2^{96}}}\right)^{6} ~= 2.87105123*10^{-7}       \f]
  !!

end module generalUtilsGAT
