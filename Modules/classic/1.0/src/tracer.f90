!> \file
!> Module containing all relevant subroutines for the model
!! tracer.
!! @author Joe Melton
!!

module tracerModule

  implicit none

  public :: prepTracer
  public :: decay14C
  public :: fractionate13C
  public :: convertTracerUnits
  public :: checkTracerBalance

  ! Set to true if you want to check for tracer balance.
  ! NOTE: See the comments in the prepTracer header for how to set
  ! up run for this check.
  logical :: doTracerBalance  = .true. !< Logical to determine if the tracer pool balance check is performed.

contains

  !> \ingroup tracermodule
  !! @{
  !> Tracks C flow through the system.
  !! No fractionation effects are applied in this subroutine.
  !! The tracer's value depends on how
  !! the model is initialized and the input file used.
  !!
  !! The tracer trackes the C movement through the green leaf,
  !! brown leaf, root, stem, litter and soil C. Carbon that is
  !! incorporated into the plants are given a tracer value of
  !! tracerValue, which corresponds to that read in from the
  !! tracerCO2 file for the year simulated. As the simulation
  !! runs and C is transferred from the living pools to the detrital
  !! pools, the tracer also is transfered.
  !!
  !! The subroutine contains a tracer pools balance check that should be
  !! run before using the tracer subroutines, in case any other model
  !! developments have not be added to the fluxes for the tracer pools but
  !! have been applied elsewhere. To do the balance check, set useTracer
  !! to 1 (for simple tracer). This ensure not fractionation or decay is
  !! performed on the tracer. Ensure that the initFile tracer pool sizes
  !! are the same as the model C pools otherwise it will fail immediately.
  !! Next set checkTracerBalance to True in the prepTracer code. This
  !! replaces the read in tracerCO2 value with a value of 1, which means
  !! the tracer is given the same amount of C as the normal C pools. Then
  !! at the end of prepTracer a small section of code is called
  !! that checks that the tracer C pools are the same size as the model C pools.
  !! If they differ it indicates a missing flux term in prepTracer or
  !! that the initFile was not set up with identical tracer and normal C pools
  !! at the start of the balance check.
  !!
  !> @author Joe Melton
  subroutine prepTracer (il1, il2, ilg, tracerCO2, & ! In
                         tracerValue) ! Out

    use classicParams, only : icc, deltat, iccp2, zero, grass, ignd, nlat
    use ctemStateVars, only : tracer, c_switch, vgat

    implicit none

    integer, intent(in) :: il1 !< other variables: il1=1
    integer, intent(in) :: il2 !< other variables: il2=ilg
    integer, intent(in) :: ilg !< number of grid cells in this latitude band.

    real, intent(out) :: tracerValue(ilg) !< Temporary variable containing the tracer CO2 value to be used. Units vary.

    real, pointer :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer :: tracerMossCMass(:)      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, pointer :: tracerMossLitrMass(:)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, intent(in) :: tracerCO2(:)           !< Tracer CO2 value read in from tracerCO2File, units vary (simple: ppm, 14C \f$\Delta ^{14}C\f$)
    real, pointer :: gleafmas(:,:)     !< Green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, pointer :: bleafmas(:,:)     !< Brown leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, pointer :: stemmass(:,:)     !< Stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, pointer :: rootmass(:,:)     !< Root mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    ! COMBAK PERLAY
    ! real,pointer :: litrmass(:,:,:)   !< Litter mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    ! real,pointer :: soilcmas(:,:,:)   !< Soil carbon mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    ! COMBAK PERLAY
    real, pointer :: Cmossmas(:)       !< C in moss biomass, \f$kg C/m^2\f$
    real, pointer :: litrmsmoss(:)     !< moss litter mass, \f$kg C/m^2\f$
    integer, pointer :: ipeatland(:) !< Peatland flag: 0 = not a peatland, 1= bog, 2 = fen
    !
    ! Local
    integer :: i, j, k

    ! Point pointers
    tracerGLeafMass   => tracer%gLeafMassgat
    tracerBLeafMass   => tracer%bLeafMassgat
    tracerStemMass    => tracer%stemMassgat
    tracerRootMass    => tracer%rootMassgat
    tracerLitrMass    => tracer%litrMassgat
    tracerSoilCMass   => tracer%soilCMassgat
    tracerMossCMass   => tracer%mossCMassgat
    tracerMossLitrMass => tracer%mossLitrMassgat
    gleafmas          => vgat%gleafmas
    bleafmas          => vgat%bleafmas
    stemmass          => vgat%stemmass
    rootmass          => vgat%rootmass
    ! COMBAK PERLAY
    ! litrmass          => vgat%litrmass
    ! soilcmas          => vgat%soilcmas
    Cmossmas          => vgat%Cmossmas
    litrmsmoss        => vgat%litrmsmoss
    ipeatland        => vgat%ipeatland

    ! ---------

    !> If doTracerBalance is true, Check for mass balance for the tracer.
    !!  First set the tracer CO2 value to 1 so it gets the same inputs as the
    !! model CO2 pools. Otherwise just use the read in tracerValue.
    if (doTracerBalance) then
      tracerValue = 1.
    else
      tracerValue = tracerCO2
    end if

    !> We don't need to do any conversion of the carbon that is uptaked
    !! in this timestep. We assume that the tracerValue should just be applied as is.

    !> If the normal pools are empty, make the tracer pools the same.
    do i = il1,il2
      do j = 1,iccp2
        if (j <= icc) then ! these are just icc sized arrays.
          if (rootmass(i,j) < zero) tracerRootMass(i,j) = 0.
          if (stemmass(i,j) < zero) tracerStemMass(i,j) = 0.
          if (bleafmas(i,j) < zero) tracerBLeafMass(i,j) = 0.
          if (gleafmas(i,j) < zero) tracerGLeafMass(i,j) = 0.
        end if
        ! COMBAK PERLAY
        ! if (sum(litrmass(i,j,:)) < zero) tracerLitrMass(i,j,:) = 0.
        ! if (sum(soilcmas(i,j,:)) < zero) tracerSoilCMass(i,j,:) = 0.
        ! COMBAK PERLAY
      end do ! j

      if (ipeatland(i) > 0) print * ,'Tracer not set up yet for peatlands.'

    end do ! i

  end subroutine prepTracer
  !! @}
  ! -------------------------------------------------------
  !> \ingroup tracermodule
  !! @{
  !> Calculates the decay of \f$^{14}C\f$ in the tracer pools.
  !!
  !! Once a year we calculate the decay of \f$^{14}C\f$ in the tracer pools.
  !! This calculation is only called if useTracer == 2.
  !! If using spinfast to equilibrate the model we need to
  !! adjust the decay of 14C in the model soil C pool. Since
  !! spinfast increases the turnover time of soil C by the factor
  !! 1/spinfast, the 14C in the soils generated with a spinfast > 1
  !! will be too young by the same factor. Following Koven et al .2013
  !! \cite Koven2013-dd we adjust for that faster equilbration by
  !! increasing the decay in the soil C tracer pool by spinfast.
  !> @author Joe Melton
  subroutine decay14C (il1, il2)

    use ctemStateVars, only : c_switch, iccp2, icc, tracer, ignd
    use classicParams, only : lambda14C

    implicit none

    integer, intent(in) :: il1 !< other variables: il1=1
    integer, intent(in) :: il2 !< other variables: il2=ilg

    real, pointer :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer :: tracerMossCMass(:)        !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, pointer :: tracerMossLitrMass(:)     !< Tracer mass in moss litter, \f$kg C/m^2\f$
    integer, pointer :: spinfast               !< Set this to a higher number up to 10 to spin up
    !< soil carbon pool faster

    integer  :: i, j, k    ! counters
    real :: dfac    !< decay constant. \f$y^{-1}\f$

    ! point pointers
    tracerGLeafMass   => tracer%gLeafMassgat
    tracerBLeafMass   => tracer%bLeafMassgat
    tracerStemMass    => tracer%stemMassgat
    tracerRootMass    => tracer%rootMassgat
    tracerLitrMass    => tracer%litrMassgat
    tracerSoilCMass   => tracer%soilCMassgat
    tracerMossCMass   => tracer%mossCMassgat
    tracerMossLitrMass => tracer%mossLitrMassgat
    spinfast          => c_switch%spinfast
    !----------
    ! begin calculations

    dfac = exp( - 1. / lambda14C)

    do i = il1,il2
      tracerMossCMass = ((tracerMossCMass(i) + 1.) * dfac) - 1.
      tracerMossLitrMass = ((tracerMossLitrMass(i) + 1.) * dfac) - 1.
      do j = 1,iccp2
        if (j <= icc) then
          tracerGLeafMass(i,j) = ((tracerGLeafMass(i,j) + 1.) * dfac) - 1.
          tracerBLeafMass(i,j) = ((tracerBLeafMass(i,j) + 1.) * dfac) - 1.
          tracerStemMass(i,j) = ((tracerStemMass(i,j) + 1.) * dfac) - 1.
          tracerRootMass(i,j) = ((tracerRootMass(i,j) + 1.) * dfac) - 1.
        end if
        do k = 1,ignd
          tracerLitrMass(i,j,k) = ((tracerLitrMass(i,j,k) + 1.) * dfac) - 1.
          tracerSoilCMass(i,j,k) = ((tracerSoilCMass(i,j,k) + 1.) * dfac * spinfast) - 1.
        end do
      end do
    end do

  end subroutine decay14C
  !! @}
  ! -------------------------------------------------------
  !> \ingroup tracermodule
  !! @{ 13C tracer -- NOT IMPLEMENTED YET.
  subroutine fractionate13C

    use ctemStateVars,only : c_switch

    implicit none

    print * ,'fractionate13C: Not implemented yet. Usetracer cannot == 3 !'
    call errorHandler('tracer', - 1)
  end subroutine fractionate13C
  !! @}
  ! -------------------------------------------------------
  !> \ingroup tracermodule
  !! @{
  !> Converts the units of the tracers, depending on the tracer
  !! being simulated.
  !!
  !! If the tracer is a simple tracer, no conversion of units is needed.
  !! If the tracer is \f$^{14}C\f$ then we apply a conversion as follows.
  !! Incoming units expected are \f$\Delta ^{14}C\f$ reported
  !! relative to the Modern standard, including corrections for
  !! age and fractionation following Stuiver and Polach 1977 \cite Stuiver1977-yj
  !!
  !! \f$\Delta ^{14}C = 1000 \left( \left[ 1 + \frac{\delta ^{14}C}{1000} \right]
  !! \frac{0.975^2}{1 + \frac{\delta ^{13}C}{1000}} - 1 \right)\f$
  !!
  !! Here \f$\delta ^{14}C\f$ is  the  measured value and \f$\Delta ^{14}C\f$ corrects
  !! for isotopic fractionation by mass-dependent processes. The 0.975 term is the
  !! fractionation of \f$^{13}C\f$ by photosynthesis. Since we have only one tracer
  !! we make a few simplifying assumptions. First we assume that no fractionation of
  !! \f$^{14}C\f$ occurs, thus 0.975 becomes 1 and the \f$\delta ^{13}C\f$ becomes 0.
  !! \f$\Delta ^{14}C\f$ then simplifies to \f$\delta ^{14}C\f$ and becomes:
  !!
  !! \f$\Delta ^{14}C = \delta ^{14}C = 1000 \left(\frac{A_s}{A_{abs}} - 1 \right)\f$
  !!
  !! where the \f$A_s\f$ is the \f$^{14}C/C\f$ ratio in a given sample. We follow Koven
  !! et al. (2013) \cite Koven2013-dd in assuming a background preindustrial atmospheric
  !! \f$^{14}C/C\f$ ratio (\f$A_{abs}\f$) of \f$10^{-12}\f$. Technically, the standard value of this is
  !! 0.2260 \f$\pm\f$ 0.0012 Bq/gC where a Bq is 433.2 x \f$10^{-15}\f$ mole (\f$^{14}C\f$).
  !! The \f$^{14}C/C\f$ ratio, or \f$A_s\f$ can then be solved to be:
  !!
  !! \f$ A_s = 10^{-12}\left( \frac{\Delta ^{14}C}{1000} + 1 \right)\f$
  !!
  !! @author Joe Melton
  function convertTracerUnits (tracerco2conc)

    use ctemStateVars, only : c_switch, tracer, nlat, nmos

    implicit none

    real, dimension(nlat,nmos)       :: convertTracerUnits
    real, dimension(:,:), intent(in) :: tracerco2conc
    integer, pointer                 :: useTracer !< useTracer = 0, the tracer code is not used.
    ! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the
    !               tracer values in the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    ! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    ! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.

    useTracer         => c_switch%useTracer

    ! ------------------------

    if (useTracer == 0) then

      print * ,'Error,entered convertTracerUnits with a useTracer value of 0'
      call errorHandler('tracer',1)

    else if (useTracer == 1) then

      ! Simple tracer, no conversion of units needed.
      convertTracerUnits = tracerco2conc

    else if (useTracer == 2) then

      ! Incoming units expected are \Delta ^{14}C so need to convert to 14C/C ratio
      convertTracerUnits = (tracerco2conc / 1000. + 1.) * 1E-12

    else if (useTracer == 3) then

      print * ,'convertTracerUnits: Error, 13C not implemented (useTracer == 3 used)'
      call errorHandler('tracer', 2)

    end if

  end function convertTracerUnits
  !! @}
  ! -------------------------------------------------------
  !> \ingroup tracermodule
  !! @{
  !> Checks for balance between the tracer pools and the
  !! model normal C pools. The subroutine is called when
  !! the doTracerBalance logical is set to true in updateSimpleTracer
  !! and the model initFile is set up as described in the
  !! notes of updateSimpleTracer. checkTracerBalance uses the
  !! same tolerance for comparisions as balcar.
  !! @author Joe Melton
  subroutine checkTracerBalance (il1, il2)

    use classicParams, only : icc, iccp2, ignd, tolrance
    use ctemStateVars, only : tracer, vgat, c_switch

    implicit none

    integer, intent(in) :: il1 !< other variables: il1=1
    integer, intent(in) :: il2 !< other variables: il2=ilg

    real, pointer :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer :: tracerMossCMass(:)      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, pointer :: tracerMossLitrMass(:)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, pointer :: gleafmas(:,:)     !< Green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, pointer :: bleafmas(:,:)     !< Brown leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, pointer :: stemmass(:,:)     !< Stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, pointer :: rootmass(:,:)     !< Root mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    ! COMBAK PERLAY
    ! real, pointer :: litrmass(:,:,:)   !< Litter mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    ! real, pointer :: soilcmas(:,:,:)   !< Soil carbon mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    integer, pointer :: useTracer !< useTracer = 0, the tracer code is not used.
    ! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the
    !               tracer values in the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    ! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    ! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.

    integer :: i, j, k

    tracerGLeafMass   => tracer%gLeafMassgat
    tracerBLeafMass   => tracer%bLeafMassgat
    tracerStemMass    => tracer%stemMassgat
    tracerRootMass    => tracer%rootMassgat
    tracerLitrMass    => tracer%litrMassgat
    tracerSoilCMass   => tracer%soilCMassgat
    tracerMossCMass   => tracer%mossCMassgat
    tracerMossLitrMass => tracer%mossLitrMassgat
    gleafmas          => vgat%gleafmas
    bleafmas          => vgat%bleafmas
    stemmass          => vgat%stemmass
    rootmass          => vgat%rootmass
    ! COMBAK PERLAY
    ! litrmass          => vgat%litrmass
    ! soilcmas          => vgat%soilcmas
    useTracer         => c_switch%useTracer

    ! ---------

    if (useTracer /= 1) then
      print * ,'ERROR ! checkTracerBalance: useTracer must == 1 for this check, you have:',useTracer
      print * ,'Either turn off checkTracerBalance (doTracerBalance == .false.) or set useTracer == 1'
      print * ,'Ending run.'
      print * ,' *** ^^^ *** ^^^ *** ^^^ *** ^^^ *** ^^^ *** ^^^ '
      call errorHandler('checkTracerBalance', 1)
    end if

    do i = il1,il2
      do j = 1,iccp2
        if (j <= icc) then

          ! green leaf mass
          if (abs(tracerGLeafMass(i,j) - gleafmas(i,j)) > tolrance) then
            print * ,'checkTracerBalance: Tracer balance fail for green leaf mass'
            print * ,i,'PFT',j,'tracerGLeafMass',tracerGLeafMass(i,j), &
                    'gleafmas',gleafmas(i,j)
            call errorHandler('checkTracerBalance', - 1)
          end if

          ! brown leaf mass
          if (abs(tracerBLeafMass(i,j) - bleafmas(i,j)) > tolrance) then
            print * ,'checkTracerBalance: Tracer balance fail for brown leaf mass'
            print * ,i,'PFT',j,'tracerBLeafMass',tracerBLeafMass(i,j), &
                    'bleafmas',bleafmas(i,j)
            call errorHandler('checkTracerBalance', - 2)
          end if

          ! stem mass
          if (abs(tracerStemMass(i,j) - stemmass(i,j)) > tolrance) then
            print * ,'checkTracerBalance: Tracer balance fail for stem mass'
            print * ,i,'PFT',j,'tracerStemMass',tracerStemMass(i,j), &
                    'stemmass',stemmass(i,j)
            call errorHandler('checkTracerBalance', - 1)
          end if

          ! root mass
          if (abs(tracerRootMass(i,j) - rootmass(i,j)) > tolrance) then
            print * ,'checkTracerBalance: Tracer balance fail for root mass'
            print * ,i,'PFT',j,'tracerRootMass',tracerRootMass(i,j), &
                    'rootmass',rootmass(i,j)
            call errorHandler('checkTracerBalance', - 1)
          end if
        end if

        ! COMBAK PERLAY
        ! do k = 1,ignd
        !   ! green leaf mass
        !   if (abs(tracerLitrMass(i,j,k) - litrmass(i,j,k)) > tolrance) then
        !     print*,'checkTracerBalance: Tracer balance fail for litter mass'
        !     print*,i,'PFT',j,'layer',k,'tracerLitrMass',tracerLitrMass(i,j,k), &
        !               'litrmass',litrmass(i,j,k)
        !     call errorHandler('checkTracerBalance',-1)
        !   end if
        !   ! green leaf mass
        !   if (abs(tracerSoilCMass(i,j,k) - soilcmas(i,j,k)) > tolrance) then
        !     print*,'checkTracerBalance: Tracer balance fail for soil C mass'
        !     print*,i,'PFT',j,'layer',k,'tracerSoilCMass',tracerSoilCMass(i,j,k), &
        !               'soilcmas',soilcmas(i,j,k)
        !     call errorHandler('checkTracerBalance',-1)
        !   end if
        ! end do
      end do
    end do

  end subroutine checkTracerBalance
  !! @}
  ! -------------------------------------------------------
  !> \namespace tracermodule
  !!
  !! Contains all relevant subroutines for the model tracer.
  !!
  !> \file
end module tracerModule
