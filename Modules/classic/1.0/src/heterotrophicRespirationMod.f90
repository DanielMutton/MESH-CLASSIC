!> \file
!> Central module for all heterotrophic respiration-related operations
module heterotrophicRespirationMod

  implicit none

  ! Subroutines contained in this module:
  public  :: heterotrophicRespiration
  public  :: hetresg
  public  :: hetresv
  public  :: updatePoolsHetResp

contains

  !> \ingroup heterotrophicrespirationmod_heterotrophicRespiration
  !! @{
  !> Heterotrophic respiration subroutine that calls hetresg and hetresv.
  !> @author Vivek Arora and Joe Melton

  subroutine heterotrophicRespiration (il1, il2, ilg, ipeatland, fcancmx, fc, & ! In
                                       litrmass, soilcmas, delzw, thpor, tbar, & ! In
                                       psisat, thliq, sort, bi, isand, thice, & ! In
                                       fg, litrmsmoss, peatdep, wtable, zbotw, roottemp, & ! In
                                       useTracer, tracerLitrMass, tracerSoilCMass, tracerMossLitrMass, & ! In
                                       ltresveg, scresveg, litresmoss, socres_peat, & ! Out
                                       resoxic, resanoxic, & ! Out
                                       ltResTracer, sCResTracer, litResMossTracer, soCResPeatTracer) ! Out

    use classicParams, only : iccp2, ignd, icc, iccp1
    use peatlandsMod,  only : hetresPeat

    implicit none

    integer, intent(in) :: il1
    integer, intent(in) :: il2
    integer, intent(in) :: ilg
    integer, intent(in) :: ipeatland(:)
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    real, intent(in) :: fcancmx(:,:)      !< max. fractional coverage of CTEM's pfts, but this can be
    !< modified by land-use change, and competition between pfts
    real, intent(in) :: fc(:)             !<
    ! COMBAK PERLAY
    real, dimension(ilg,icc), intent(in) :: roottemp       !< root temperature, k
    real, intent(in) :: litrmass(:,:)   !< litter mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, intent(in) :: soilcmas(:,:)   !< soil carbon mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    ! real, intent(in) :: litrmass(:,:,:)   !< litter mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    ! real, intent(in) :: soilcmas(:,:,:)   !< soil carbon mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    ! COMBAK PERLAY
    real, intent(in) :: delzw(:,:)         !< thicknesses of the soil layers
    real, intent(in) :: thpor(:,:)         !< Soil total porosity \f$(cm^3 cm^{-3})\f$ -
    real, intent(in) :: tbar(:,:)          !< Soil temperature, K
    real, intent(in) :: psisat(:,:)        !< Saturated soil matric potential (m)
    real, intent(in) :: thliq(:,:)         !< liquid mois. content of soil layers
    integer, intent(in) :: sort(:)         !<
    real, intent(in) :: bi(:,:)            !< Brooks and Corey/Clapp and Hornberger b term
    integer, intent(in) :: isand(:,:)      !<
    real, intent(in) :: thice(:,:)         !< frozen mois. content of soil layers
    real, intent(in) :: fg(:)              !<
    real, intent(inout) :: litrmsmoss(:)   !< moss litter C pool (kgC/m2)
    real, intent(inout) :: peatdep(:)      !< peat depth (m)
    real, intent(in) :: wtable(:)          !< water table (m)
    real, intent(in) :: zbotw(:,:)         !< bottom of soil layers

    ! COMBAK PERLAY
    real, intent(out) :: ltresveg(ilg,iccp2)     !< fluxes for each pft: litter respiration for each pft + bare fraction
    real, intent(out) :: scresveg(ilg,iccp2)     !< soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
    ! real, intent(out) :: ltresveg(ilg,iccp2,ignd)     !< fluxes for each pft: litter respiration for each pft + bare fraction
    ! real, intent(out) :: scresveg(ilg,iccp2,ignd)     !< soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
    ! COMBAK PERLAY
    real, intent(out) :: litresmoss(ilg)    !< moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: socres_peat(ilg)   !< heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: resoxic(ilg)       !< oxic respiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: resanoxic(ilg)     !< anoxic respiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real, intent(in) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(in) :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(in) :: tracerMossLitrMass(:)   !< Tracer mass in moss litter, \f$kg C/m^2\f$


    real :: ltResTracer(ilg,iccp2,ignd)  !< Tracer fluxes for each pft: litter respiration for each pft + bare fraction TODO: Units (intent(out))
    real :: sCResTracer(ilg,iccp2,ignd)  !< Tracer soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts TODO: Units (intent(out))
    real, intent(out) :: litResMossTracer(ilg)    !< Tracer moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) TODO: Units
    real, intent(out) :: soCResPeatTracer(ilg)    !< Tracer heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) TODO: Units

    ! ------

    ! COMBAK PERLAY
    call hetresv(fcancmx, fc, litrmass(:,1:icc), soilcmas(:,1:icc), & ! In
                 delzw, thpor, il1, il2, & ! In
                 ilg, tbar, psisat,thliq, & ! In
                 roottemp, zbotw, sort, bi, & ! In
                 isand, thice, ipeatland, & ! In
                 ltresveg, scresveg) ! Out
    ! call    hetresv ( fcancmx,fc,litrmass(:,1:icc,:),soilcmas(:,1:icc,:), & ! In
    ! delzw,thpor,il1,il2, & ! In
    ! ilg,tbar,psisat,thliq, & ! In
    ! sort,bi,isand,thice,ipeatland, & ! In
    ! useTracer,tracerLitrMass(:,1:icc,:),tracerSoilCMass(:,1:icc,:), & ! In
    ! ltresveg(:,1:icc,:),scresveg(:,1:icc,:), & ! Out
    ! ltResTracer(:,1:icc,:),sCResTracer(:,1:icc,:)) ! Out
    ! COMBAK PERLAY

    !! Find heterotrophic respiration rates from bare ground subarea
    ! COMBAK PERLAY
    call hetresg(litrmass(:,iccp1), soilcmas(:,iccp1), delzw, thpor, & ! In
                 il1, il2, ilg, tbar, & ! In
                 psisat, bi, thliq, zbotw, & ! In
                 thice, fg, isand, & ! In
                 ltresveg(:,iccp1), scresveg(:,iccp1)) ! Out
    ! call  hetresg  (litrmass(:,iccp1,:),soilcmas(:,iccp1,:),delzw,thpor, & ! In
    !    il1,il2,ilg,tbar,    & ! In
    ! psisat,bi,thliq,      & ! In
    ! thice,fg,isand,    & ! In
    ! useTracer,tracerLitrMass(:,iccp1,:),tracerSoilCMass(:,iccp1,:), & ! In
    ! ltresveg(:,iccp1,:),scresveg(:,iccp1,:), & ! Out
    ! ltResTracer(:,iccp1,:),sCResTracer(:,iccp1,:)) ! Out
    ! COMBAK PERLAY

    !>  Find heterotrophic respiration rates for the LUC litter and soilc pools
    !!  The LUC litter and soil C respiration rates are assumed to
    !!  be applied over the entire tile but kept in layer 1
    ! COMBAK PERLAY
    call hetresg(litrmass(:,iccp2), soilcmas(:,iccp2), delzw, thpor, & ! In
                 il1, il2, ilg, tbar, & ! In
                 psisat, bi, thliq, zbotw, & ! In
                 thice, fg, isand, & ! In
                 ltresveg(:,iccp2), scresveg(:,iccp2)) ! Out
    ! call  hetresg  (litrmass(:,iccp2,:),soilcmas(:,iccp2,:),delzw,thpor, & ! In
    !    il1,il2,ilg,tbar,    & ! In
    ! psisat,bi,thliq,   & ! In
    ! thice,fg,isand,    & ! In
    ! useTracer,tracerLitrMass(:,iccp2,:),tracerSoilCMass(:,iccp2,:), & ! In
    ! ltresveg(:,iccp2,:),scresveg(:,iccp2,:), & ! Out
    ! ltResTracer(:,iccp2,:),sCResTracer(:,iccp2,:)) ! Out
    ! COMBAK PERLAY

    !! When peatlands are simulated find the peat heterotrophic respiration.
    ! Called regardless to set outputs to 0.
    call hetresPeat(il1, il2, ilg, ipeatland, & ! In
                    isand, litrmsmoss, peatdep, wtable, & ! In
                    tbar, thliq, thice, thpor, & ! In
                    bi, zbotw, delzw, psisat, & ! In
                    useTracer, tracerMossLitrMass, & ! In
                    litresmoss, socres_peat, resoxic, resanoxic, & ! Out
                    litResMossTracer, soCResPeatTracer) ! Out


  end subroutine heterotrophicRespiration
  !! @}
  !!
  !> \ingroup heterotrophicrespirationmod_hetresg
  !! @{
  !> Heterotrophic respiration subroutine for bare ground fraction
  !> @author Vivek Arora and Joe Melton
  ! COMBAK PERLAY
  subroutine hetresg (litrmass, soilcmas, delzw, thpor, &! In
                      il1, il2, ilg, tbar, &! In
                      psisat, b, thliq, zbotw, &! In
                      thice, frac, isand, &! In
                      litres, socres)! Out
    ! subroutine hetresg (litrmass, soilcmas, delzw, thpor,    & ! In
    !                          il1, il2, ilg, tbar,    & ! In
    !                        psisat, b, thliq,   & ! In
    !                        thice, frac, isand,    & ! In
    !                        useTracer, tracerLitrMass, tracerSoilCMass, & ! In
    !                        litres, socres, & ! Out
    !                        ltResTracer, sCResTracer) ! Out
    ! COMBAK PERLAY

    !               Canadian Terrestrial Ecosystem Model (CTEM)
    !           Heterotrophic Respiration Subroutine For Bare Fraction
    !
    !     11  Apr. 2003 - this subroutine calculates heterotrophic respiration
    !     V. Arora        over the bare subarea of a grid cell (i.e. ground only
    !                     and snow over ground subareas).
    !
    !     change history:

    !      7  Jun 2016  - Bring in step function for reduction in resp when soil freezes.
    !     J. Melton
    !
    !      6  Jun 2016  - Het resp reduces at depth in soil column.
    !     J. Melton
    !
    !      8  Feb 2016  - Adapted subroutine for multilayer soilc and litter (fast decaying)
    !     J. Melton       carbon pools
    !
    !     19  Jan 2016  - Implemented new LUC litter and soil C pools
    !     J. Melton
    !
    !     14  Jan 2016  - Converted to F90 and made it so it can handle > 3 soil layers
    !     J. Melton
    !
    !     30  Jul 2015  - Based on work by Yuanqiao Wu, respiration was found to
    !     J. Melton       behave incorrectly if the soil froze as it thought the water
    !                     was leaving the soil. This is now fixed.
    !
    !     17  Jan 2014  - Moved parameters to global file (classicParams.f90)
    !     J. Melton
    !
    !     23  Jul 2013  - add in module for parameters
    !     J. Melton
    !     J. Melton and V.Arora - changed tanhq10 parameters, they were switched
    !               25 Sep 2012
    !     J. Melton 31 Aug 2012 - remove isnow, it is not used.
    !     J. Melton 23 Aug 2012 - bring in isand, converting sand to
    !                             int was missing some gridcells assigned
    !                             to bedrock in classb

    ! COMBAK PERLAY
    use classicParams,        only : icc, ignd, zero, tanhq10, a_hetr, &
                                     bsratelt_g, bsratesc_g
    implicit none
    integer, intent(in) :: ilg   !<
    integer, intent(in) :: il1   !< il1=1
    integer, intent(in) :: il2   !< il2=ilg
    integer :: i, j, k
    integer, intent(in) :: isand(ilg,ignd) !<
    real, intent(in) :: litrmass(ilg,1)    !< litter mass for bare in \f$kg c/m^2\f$
    real, intent(in) :: soilcmas(ilg,1)    !< soil carbon mass for bare in \f$kg c/m^2\f$
    real, intent(in) :: tbar(ilg,ignd)     !< soil temperature, k
    real, intent(in) :: thliq(ilg,ignd)    !< liquid soil moisture content in soil layers
    real, intent(in) :: zbotw(ilg,ignd)    !< bottom of soil layers
    real, intent(out) :: litres(ilg)        !< litter respiration over the given unvegetated sub-area in umol co2/m2.s
    real, intent(out) :: socres(ilg)        !< soil c respiration over the given unvegetated sub-area in umol co2/m2.s
    real, intent(in) :: frac(ilg)          !< fraction of ground
    real, intent(in) :: delzw(ilg,ignd)  !<
    real, intent(in) :: thice(ilg,ignd) !<
    real :: zcarb_g          !<
    real :: litrq10          !<
    real :: soilcq10         !<
    real :: litrtemp(ilg)    !< litter temperature
    real :: solctemp(ilg)    !< soil carbon pool temperature
    real :: q10func          !<
    real, intent(in) :: psisat(ilg,ignd) !< saturation matric potential
    real :: grksat(ilg,ignd) !< saturation hyd. conductivity
    real, intent(in) :: b(ilg,ignd)      !< parameter b of clapp and hornberger
    real, intent(in) :: thpor(ilg,ignd)  !< porosity
    real :: beta             !<
    real :: fracarb(ilg,ignd)!< fraction of carbon in each soil layer
    real :: zcarbon          !<
    real :: tempq10l(ilg)    !<
    real :: socmoscl(ilg)    !< soil moisture scalar for soil carbon decomposition
    real :: scmotrm(ilg,ignd)!<
    real :: ltrmoscl(ilg)    !< soil moisture scalar for litter decomposition
    real :: psi(ilg,ignd)    !<
    real :: tempq10s(ilg)    !<
    real :: fcoeff           !<
    !     initialize required arrays to zero
    fracarb(:,:) = 0.0  ! fraction of carbon in each soil layer
    solctemp(:) = 0.0     ! soil carbon pool temperature
    socmoscl(:) = 0.0     ! soil moisture scalar for soil carbon decomposition
    ltrmoscl(:) = 0.0     ! soil moisture scalar for litter decomposition
    litres(:) = 0.0       ! litter resp. rate
    tempq10l(:) = 0.0
    socres(:) = 0.0       ! soil c resp. rate
    tempq10s(:) = 0.0
    !     initialization ends
    !     ------------------------------------------------------------------
    !>     Estimate temperature of the litter and soil carbon pools.
    !!     Over the bare fraction there is no live root. So we make the
    !!     simplest assumption that litter temperature is same as temperature of the top soil layer.
    do i = il1,il2 ! loop 210
      litrtemp(i) = tbar(i,1)
    end do ! loop 210
    !>     We estimate the temperature of the soil c pool assuming that soil carbon over the bare fraction is distributed exponentially. note
    !!     that bare fraction may contain dead roots from different pfts all of which may be distributed differently. For simplicity we do not
    !!     track each pft's dead root biomass and assume that distribution of soil carbon over the bare fraction can be described by a single
    !!     parameter.
    do i = il1,il2 ! loop 245
      zcarbon = 3.0 / a_hetr
      zcarb_g = 0.0
      do j = 1,ignd
        zcarb_g = zcarb_g + delzw(i,j)
      end do
      zcarbon = min(zcarbon,zcarb_g)
      fcoeff = exp( - a_hetr * zcarbon)
      fracarb(i,1) = 1.0 - (exp( - a_hetr * zbotw(i,1)) - fcoeff)/(1.0 - fcoeff)
      do j = 2,ignd
        if (zcarbon <= zbotw(i,j) - delzw(i,j) + 0.0001) then
          fracarb(i,j) = 0.0
        else if (zcarbon <= zbotw(i,j)) then
          fracarb(i,j) = (exp( - a_hetr * zbotw(i,j) - delzw(i,j)) - exp( - a_hetr * zcarbon))/(1. - exp( - a_hetr * zcarbon))
        else
          fracarb(i,j) = (exp( - a_hetr * zbotw(i,j) - delzw(i,j)) - exp( - a_hetr * zbotw(i,j)))/(1. - exp( - a_hetr * zcarbon))
        end if
      end do ! loop 245
      ! -------------
      solctemp(i) = sum(tbar(i,:) * fracarb(i,:))/ sum(fracarb(i,:))
    end do ! loop 240
    !>     find moisture scalar for soil c decomposition
    !!     this is modelled as function of logarithm of matric potential.
    !!     we find values for all soil layers, and then find an average value
    !!     based on fraction of carbon present in each layer.
    do j = 1,ignd ! loop 260
      do i = il1,il2 ! loop 270
        if (isand(i,j) == - 3 .or. isand(i,j) == - 4) then
          scmotrm (i,j) = 0.2
          psi (i,j) = 10000.0 ! set to large number so that
          ! ltrmoscl becomes 0.2
        else ! i.e., sand/=-3 or -4

          ! We don't place a lower limit on psi as it is only used here and the
          ! value of psi <= psisat is just used to select a scomotrm value.
          if (thice(i,j) <= thpor(i,j)) then ! flag new limits
            psi(i,j) = psisat(i,j) * (thliq(i,j)/(thpor(i,j) - thice(i,j))) ** ( - b(i,j))
          else
            ! if the whole pore space is ice then suction is assumed to be very high.
            psi(i,j) = 10000.0
          end if
          if (psi(i,j) >= 10000.0) then
            scmotrm(i,j) = 0.2
          else if (psi(i,j) < 10000.0 .and.  psi(i,j) > 6.0) then
            scmotrm(i,j) = 1.0 - 0.8 * ( (log10(psi(i,j)) - log10(6.0))/(log10(10000.0) - log10(6.0)) )
          else if (psi(i,j) <= 6.0 .and.  psi(i,j) >= 4.0) then
            scmotrm(i,j) = 1.0
          else if (psi(i,j) < 4.0 .and. psi(i,j) > psisat(i,j) ) then
            scmotrm(i,j) = 1.0 - 0.5 * ( (log10(4.0) - log10(psi(i,j))) /(log10(4.0) - log10(psisat(i,j))) )
          else if (psi(i,j) <= psisat(i,j) ) then
            scmotrm(i,j) = 0.5
          end if
        end if ! if sand==-3 or -4
        scmotrm(i,j) = max(0.0,min(scmotrm(i,j),1.0))
      end do ! loop 270
    end do ! loop 260
    do i = il1,il2 ! loop 290
      socmoscl(i) = sum(scmotrm(i,:) * fracarb(i,:)) / sum(fracarb(i,:))
      socmoscl(i) = max(0.2,min(socmoscl(i),1.0))
    end do ! loop 290
    !>     find moisture scalar for litter decomposition
    !!     the difference between moisture scalar for litter and soil c
    !!     is that the litter decomposition is not constrained by high
    !!     soil moisture (assuming that litter is always exposed to air).
    !!     in addition, we use moisture content of the top soil layer
    !!     as a surrogate for litter moisture content. so we use only
    !!     psi(i,1) calculated in loops 260 and 270 above.
    do i = il1,il2 ! loop 300
      if (psi(i,1) > 10000.0) then
        ltrmoscl(i) = 0.2
      else if (psi(i,1) <= 10000.0 .and.  psi(i,1) > 6.0) then
        ltrmoscl(i) = 1.0 - 0.8 * ( (log10(psi(i,1)) - log10(6.0))/(log10(10000.0) - log10(6.0)) )
      else if (psi(i,1) <= 6.0) then
        ltrmoscl(i) = 1.0
      end if
      ltrmoscl(i) = max(0.2,min(ltrmoscl(i),1.0))
    end do ! loop 300
    !!     use temperature of the litter and soil c pools, and their soil
    !!     moisture scalars to find respiration rates from these pools
    do i = il1,il2 ! loop 330
      if (frac(i) > zero) then
        !       first find the q10 response function to scale base respiration
        !       rate from 15 c to current temperature, we do litter first
        tempq10l(i) = litrtemp(i) - 273.16
        litrq10 = tanhq10(1) + tanhq10(2) * ( tanh( tanhq10(3) * (tanhq10(4) - tempq10l(i))  ) )
        q10func = litrq10 ** (0.1 * (litrtemp(i) - 273.16 - 15.0))
        litres(i) = ltrmoscl(i) * litrmass(i,1) * bsratelt_g * 2.64 * q10func ! 2.64 converts bsratelt_g from kg c/kg c.year
        ! to u-mol co2/kg c.s
        !       respiration from soil c pool
        tempq10s(i) = solctemp(i) - 273.16
        soilcq10 = tanhq10(1) + tanhq10(2) * ( tanh( tanhq10(3) * (tanhq10(4) - tempq10s(i))  ) )
        q10func = soilcq10 ** (0.1 * (solctemp(i) - 273.16 - 15.0))
        socres(i) = socmoscl(i) * soilcmas(i,1) * bsratesc_g * 2.64 * q10func ! 2.64 converts bsratesc_g from kg c/kg c.year
        ! to u-mol co2/kg c.s
      end if
    end do ! loop 330

    !   use classicParams,         only : icc,ignd,zero,tanhq10,a_hetr, &
    !                                     bsratelt_g,bsratesc_g,r_depthredu, &
    !                                     tcrit,frozered
    !
    !   implicit none
    !
    !   integer, intent(in) :: ilg   !<
    !   integer, intent(in) :: il1   !< il1=1
    !   integer, intent(in) :: il2   !< il2=ilg
    !   integer, intent(in) :: isand(:,:) !<
    !   integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !                                 !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !                                 !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !                                 !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !                                 !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    !
    !   real, intent(in) :: litrmass(:,:)   !< litter mass for the CTEM pfts + bare in \f$kg c/m^2\f$
    !   real, intent(in) :: soilcmas(:,:)   !< soil carbon mass for the CTEM pfts + bare in \f$kg c/m^2\f$
    !   real, intent(in) :: tbar(:,:)    !< soil temperature, k
    !   real, intent(in) :: thliq(:,:)   !< liquid soil moisture content in soil layers
    !   real, intent(in) :: frac(:)         !< fraction of bare ground (fg)
    !   real, intent(in) :: delzw(:,:)   !<
    !   real, intent(in) :: thice(:,:)   !<
    !   real, intent(in) :: psisat(:,:)  !< saturation matric potential
    !   real, intent(in) :: b(:,:)       !< parameter b of clapp and hornberger
    !   real, intent(in) :: thpor(:,:)   !< porosity
    !   real, intent(in) :: tracerLitrMass(:,:)   !< Tracer litter mass for the CTEM pfts + bare in \f$kg c/m^2\f$
    !   real, intent(in) :: tracerSoilCMass(:,:)   !< Tracer soil carbon mass for the CTEM pfts + bare in \f$kg c/m^2\f$
    !
    !   real, intent(out) :: litres(ilg,ignd)      !< litter respiration over the given unvegetated sub-area in umol co2/m2.s
    !   real, intent(out) :: socres(ilg,ignd)      !< soil c respiration over the given unvegetated sub-area in umol co2/m2.s
    !   real, intent(out) :: ltResTracer(ilg,ignd)      !< Tracer litter respiration over the given unvegetated sub-area in umol co2/m2.s
    !   real, intent(out) :: sCResTracer(ilg,ignd)      !< Tracer soil C respiration over the given unvegetated sub-area in umol co2/m2.s
    !
    !   ! Local
    !   real :: litrq10          !<
    !   real :: soilcq10         !<
    !   real :: q10func          !<
    !   real :: tempq10l         !<
    !   real :: socmoscl(ilg,ignd)    !< soil moisture scalar for soil carbon decomposition
    !   real :: ltrmoscl(ilg,ignd)    !< soil moisture scalar for litter decomposition
    !   real :: psi(ilg,ignd)    !<
    !   real :: tempq10s         !<
    !   real :: reduceatdepth
    !   integer :: i, j, k
    !
    !   ! -------------------------------------------------------------------------
    !
    !   ! initialize required arrays to zero
    !
    !   socmoscl(:,:)=0.0   ! soil moisture scalar for soil carbon decomposition
    !   ltrmoscl(:,:)=0.0   ! soil moisture scalar for litter decomposition
    !   litres(:,:)=0.0     ! litter resp. rate
    !   socres(:,:)=0.0     ! soil c resp. rate
    !   ltResTracer = 0.
    !   sCResTracer = 0.
    !
    !   !     initialization ends
    !
    !   !> Find moisture scalar for soil c decomposition
    !   !! this is modelled as function of logarithm of matric potential.
    !   !! we find values for all soil layers, and then find an average value
    !   !! based on fraction of carbon present in each layer.
    !
    !   do j = 1,ignd
    !     do i = il1,il2
    !
    !       if (isand(i,j) == -3 .or. isand(i,j) == -4) then
    !         socmoscl(i,j)=0.2
    !         psi (i,j) = 10000.0 ! set to large number so that ltrmoscl becomes 0.2
    !       else ! i.e., sand/=-3 or -4
    !
    !         ! We don't place a lower limit on psi as it is only used here and the
    !         ! value of psi <= psisat is just used to select a scomotrm value.
    !         if (thice(i,j) <= thpor(i,j)) then ! flag new limits
    !           psi(i,j) = psisat(i,j) * (thliq(i,j) / (thpor(i,j) - thice(i,j)))**(-b(i,j))
    !         else
    !           ! if the whole pore space is ice then suction is assumed to be very high.
    !           psi(i,j) = 10000.0
    !         end if
    !
    !         if (psi(i,j) >= 10000.0) then
    !           socmoscl(i,j) = 0.2
    !         else if (psi(i,j) < 10000.0 .and.  psi(i,j) > 6.0) then
    !           socmoscl(i,j) = 1.0 - 0.8 * ( (log10(psi(i,j)) - log10(6.0)) &
    !                                       / (log10(10000.0) - log10(6.0)) )
    !         else if (psi(i,j) <= 6.0 .and. psi(i,j) >= 4.0) then
    !           socmoscl(i,j) = 1.0
    !         else if (psi(i,j) < 4.0 .and. psi(i,j) > psisat(i,j)) then
    !           socmoscl(i,j) = 1.0 - 0.5 * ((log10(4.0) - log10(psi(i,j))) &
    !                                      / (log10(4.0) - log10(psisat(i,j))) )
    !         else if ( psi(i,j) <= psisat(i,j)) then
    !           socmoscl(i,j) = 0.5
    !         end if
    !       end if
    !
    !       socmoscl(i,j) = max(0.2,min(socmoscl(i,j),1.0))
    !
    ! 270   continue
    ! 260 continue
    !
    !   !> Find moisture scalar for litter decomposition for top soil layer
    !   !> The difference between moisture scalar for litter and soil c
    !   !! is that the litter decomposition is not constrained by high
    !   !! soil moisture (assuming that litter is always exposed to air).
    !   !! In addition, we use moisture content of the top soil layer
    !   !! as a surrogate for litter moisture content. So we use only
    !   !! psi(i,1) calculated in loops 260 and 270 above.
    !
    !   do i = il1,il2
    !     if (psi(i,1) > 10000.0) then
    !       ltrmoscl(i,1) = 0.2
    !     else if (psi(i,1) <= 10000.0 .and.  psi(i,1) > 6.0) then
    !       ltrmoscl(i,1) = 1.0 - 0.8 * ((log10(psi(i,1)) - log10(6.0)) &
    !                                   /(log10(10000.0) - log10(6.0)) )
    !     else if (psi(i,1) <= 6.0) then
    !       ltrmoscl(i,1) = 1.0
    !     end if
    !     ltrmoscl(i,1) = max(0.2,min(ltrmoscl(i,1),1.0))
    ! 300   continue
    !
    !   !> Treat the lower litter layers like the soil C ones:
    !   ltrmoscl(:,2:ignd) = socmoscl(:,2:ignd)
    !
    ! !> Use temperature of the litter and soil c pools,and their soil
    ! !! moisture scalars to find respiration rates from these pools
    !
    !   do i = il1,il2
    !     do j = 1,ignd
    !       if (frac(i)>zero) then
    !
    !         !! First find the q10 response function to scale base respiration
    !         !! rate from 15 c to current temperature, we do litter first
    !         tempq10l = tbar(i,j) - 273.16
    !         litrq10 = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - tempq10l)))
    !
    !         !! Apply a step reduction in q10func when the soil layer freezes. This is to reflect the
    !         !! lost mobility of the population due to less liquid water
    !         if (tbar(i,j) - 273.16 > tcrit) then
    !           q10func = litrq10**(0.1* (tbar(i,j) - 273.16 - 15.0))
    !         else
    !           q10func = litrq10**(0.1 * (tbar(i,j) - 273.16 - 15.0)) * frozered
    !         end if
    !
    !         !> Reduce the respiration at depth due to unresolved depth dependent processes including
    !         !! soil microbial population dynamics, pore-scale oxygen availability, mineral sorption,
    !         !! priming effects, and other unrepresented processes. This is following Lawrence et al.
    !         !! Enviro Res Lett 2015 \cite Lawrence2015-tj. We apply this for all soils.
    !         reduceatdepth = exp(-delzw(i,j) / r_depthredu)
    !         ! 2.64 converts bsratelt_g from kg c/kg c.year to u-mol co2/kg c.s
    !         litres(i,j) = ltrmoscl(i,j) * litrmass(i,j) * bsratelt_g * 2.64 &
    !                       * q10func * reduceatdepth
    !
    !         if (useTracer > 0) ltResTracer(i,j) = ltrmoscl(i,j) * tracerLitrMass(i,j) &
    !                                                         * bsratelt_g * 2.64 * q10func * reduceatdepth
    !
    !         ! Respiration from soil c pool
    !         tempq10s = tbar(i,j) - 273.16
    !         soilcq10 = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - tempq10s)))
    !
    !         if (tbar(i,j) - 273.16 > tcrit) then
    !             q10func = soilcq10**(0.1 * (tbar(i,j) - 273.16 - 15.0))
    !         else
    !             q10func = soilcq10**(0.1 * (tbar(i,j) - 273.16 - 15.0)) * frozered
    !         end if
    !         ! 2.64 converts bsratelt_g from kg c/kg c.year to u-mol co2/kg c.s
    !         socres(i,j) = socmoscl(i,j) * soilcmas(i,j) * bsratesc_g * 2.64 * q10func * reduceatdepth
    !
    !         if (useTracer > 0) sCResTracer(i,j) = socmoscl(i,j) * tracerSoilCMass(i,j) &
    !                                                         * bsratesc_g * 2.64 * q10func * reduceatdepth
    !
    !       end if
    !
    ! 340  continue
    ! 330 continue
    ! COMBAK PERLAY  !

    return

  end subroutine hetresg
  !! @}
  ! ------------------------------------------------------------------------------------

  !> \ingroup heterotrophicrespirationmod_hetresv
  !! @{
  !> Heterotrophic Respiration Subroutine For Vegetated Fraction
  !> @author Vivek Arora, Joe Melton, Yuanqiao Wu
  ! COMBAK PERLAY
  subroutine hetresv (fcan, fct, litrmass, soilcmas, &! In
                      delzw, thpor, il1, il2, &! In
                      ilg, tbar, psisat, thliq, &! In
                      roottemp, zbotw, sort, b, &! In
                      isand, thice, ipeatland, &! In
                      ltresveg, scresveg) ! Out
    ! subroutine hetresv ( fcan, fct, litrmass, soilcmas,  & ! In
    !                     delzw, thpor, il1, il2,  & ! In
    !                     ilg, tbar, psisat, thliq,  & ! In
    !                     sort, b, isand, thice, ipeatland, & ! In
    !                     useTracer, tracerLitrMass, tracerSoilCMass, & ! In
    !                     ltresveg, scresveg, ltResTracer, sCResTracer) ! Out
    ! COMBAK PERLAY
    !               Canadian Terrestrial Ecosystem Model (CTEM)
    !           Heterotrophic Respiration Subtoutine For Vegetated Fraction

    !     16  oct. 2001 - this subroutine calculates heterotrophic respiration
    !     V. Arora        for a given sub-area, from litter and soil carbon
    !                     pools.

    !     change history:

    !
    !      7  Jun 2016  - Bring in step function for reduction in resp when soil freezes.
    !     J. Melton
    !
    !      6  Jun 2016  - Het resp reduces at depth in soil column.
    !     J. Melton
    !
    !      8  Feb 2016  - Adapted subroutine for multilayer soilc and litter (fast decaying)
    !     J. Melton       carbon pools
    !
    !     14  Jan 2016  - Converted to F90 and made it so it can handle > 3 soil layers
    !     J. Melton

    !     10  April 2015 -Bring in peatland scheme
    !     Y. Wu
    !
    !     30  Jul 2015  - Based on work by Yuanqiao Wu, respiration was found to
    !     J. Melton       behave incorrectly if the soil froze as it thought the water
    !                     was leaving the soil. This is now fixed.
    !     17  Jan 2014  - Moved parameters to global file (classicParams.f90)
    !     J. Melton

    !     22  Jul 2013  - Add in module for parameters
    !     J. Melton

    !     J. Melton and V.Arora - changed tanhq10 parameters, they were switched
    !               25 sep 2012
    !     J. Melton 23 aug 2012 - bring in isand, converting sand to
    !                             int was missing some gridcells assigned
    !                             to bedrock in classb
    !     ------

    ! COMBAK PERLAY
    use classicParams,  only : icc, ignd, kk, zero, bsratelt, &
                               bsratesc, abar, tanhq10, &
                               alpha_hetres
    implicit none
    integer, intent(in) :: ilg       !<
    integer, intent(in) :: il1       !< il1=1
    integer, intent(in) :: il2       !< il2=ilg
    integer :: i, j, k
    integer, intent(in) :: sort(icc) !< index for correspondence between 9 pfts and 12 values in the parameters vectors
    integer, intent(in) :: isand(ilg,ignd) !<
    integer, intent(in)  :: ipeatland(ilg) !<
    real, intent(in) :: fcan(ilg,icc)      !< fractional coverage of ctem's 9 pfts
    real, intent(in) :: fct(ilg)           !< sum of all fcan, fcan & fct are not used at this time but could be used at some later stage
    real, intent(in) :: litrmass(ilg,icc)!< litter mass for the 9 pfts in \f$kg c/m^2\f$
    real, intent(in) :: tbar(ilg,ignd)     !< soil temperature, k
    real, intent(in) :: soilcmas(ilg,icc)!< soil carbon mass for the 9 pfts e in \f$kg c/m^2\f$
    real, intent(in) :: thliq(ilg,ignd)    !< liquid soil moisture content in 3 soil layers
    real, intent(in) :: roottemp(ilg,icc)  !< root temperature as estimated in mainres subroutine
    real, intent(in) :: zbotw(ilg,ignd)    !< bottom of soil layers
    real, intent(out) :: ltresveg(ilg,icc)  !< litter respiration for the given sub-area in umol co2/m2.s, for ctem's 9 pfts
    real, intent(out) :: scresveg(ilg,icc)  !< soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's 9 pfts
    real, intent(in) :: thice(ilg,ignd)   !< frozen soil moisture content in the soil layers
    real, intent(in) :: delzw(ilg,ignd)    !<
    real :: zcarb_g            !<
    real :: litrq10            !<
    real :: soilcq10           !<
    real :: litrtemp(ilg,icc)  !< litter temperature
    real :: solctemp(ilg,icc)  !< soil carbon pool temperature
    real :: q10func            !<
    real, intent(in) :: psisat(ilg,ignd)   !< saturation matric potential
    real :: grksat(ilg,ignd)   !< saturation hyd. conductivity
    real, intent(in) :: b(ilg,ignd)        !< parameter b of clapp and hornberger
    real, intent(in) :: thpor(ilg,ignd)    !< porosity
    real :: fracarb(ilg,icc,ignd) !< fraction of carbon in each soil layer for each vegetation
    real :: zcarbon            !<
    real :: tempq10l(ilg,icc)  !<
    real :: socmoscl(ilg,icc)  !< soil moisture scalar for soil carbon decomposition
    real :: scmotrm(ilg,ignd)  !< soil carbon moisture term
    real :: ltrmoscl(ilg)      !< soil moisture scalar for litter decomposition
    real :: psi(ilg,ignd)      !<
    real :: tempq10s(ilg,icc)  !<
    real :: fcoeff             !<
    !     ------------------------------------------------------------------
    !!     Constants and parameters are located in classicParams.f90
    !     ---------------------------------------------------------------
    !     initialize required arrays to zero
    do j = 1,icc ! loop 100
      do i = il1,il2 ! loop 110
        litrtemp(i,j) = 0.0       ! litter temperature
        tempq10l(i,j) = 0.0
        solctemp(i,j) = 0.0       ! soil carbon pool temperature
        tempq10s(i,j) = 0.0
        socmoscl(i,j) = 0.0       ! soil moisture scalar for soil carbon decomposition
        ltresveg(i,j) = 0.0       ! litter resp. rate for each pft
        scresveg(i,j) = 0.0       ! soil c resp. rate for each pft
        do k = 1,ignd ! loop 120
          scmotrm(i,k) = 0.0        ! soil carbon moisture term
          fracarb(i,j,k) = 0.0    ! fraction of carbon in each soil layer for each vegetation
        end do ! loop 120
      end do ! loop 110
    end do ! loop 100
    do i = il1,il2
      ltrmoscl(i) = 0.0           ! soil moisture scalar for litter decomposition
    end do ! loop 130
    !     initialization ends
    !     ------------------------------------------------------------------
    !> estimate temperature of the litter and soil carbon pools. litter
    !! temperature is weighted average of temperatue of top soil layer
    !! (where the stem and leaf litter sits) and root temperature, because
    !! litter pool is made of leaf, stem, and root litter.
    do j = 1,icc ! loop 200
      do i = il1,il2 ! loop 210
        if (fcan(i,j) > 0.) then
          litrtemp(i,j) = alpha_hetres * tbar(i,1) + roottemp(i,j) * (1.0 - alpha_hetres)
        end if
      end do ! loop 210
    end do ! loop 200
    !> estimation of soil carbon pool temperature is not straight forward.
    !! ideally soil c pool temperature should be set same as root temperature,
    !! since soil c profiles are similar to root profiles. but in the event
    !! when the roots die then we may run into trouble. so we find the
    !! temperature of the soil c pool assuming that soil carbon is
    !! exponentially distributed, just like roots. but rather than using
    !! the parameter of this exponential profile from our variable root
    !! distribution we use fixed vegetation-dependent parameters.
    do j = 1,icc ! loop 230
      do i = il1,il2 ! loop 240
        if (fcan(i,j) > 0.) then
          zcarbon = 3.0 / abar(sort(j))
          zcarb_g = 0.0
          do k = 1,ignd
            zcarb_g = zcarb_g + delzw(i,k)
          end do
          zcarbon = min(zcarbon,zcarb_g)
          fcoeff = exp( - abar(sort(j)) * zcarbon)
          fracarb(i,j,1) = 1.0 - (exp( - abar(sort(j)) * zbotw(i,1)) - fcoeff)/(1.0 - fcoeff)
          do k = 2,ignd
            if (zcarbon <= zbotw(i,k) - delzw(i,k) + 0.0001) then
              fracarb(i,j,k) = 0.0
            else if (zcarbon <= zbotw(i,k)) then
              fracarb(i,j,k) = (exp( - abar(sort(j)) * zbotw(i,k) - delzw(i,k)) - &
                               exp( - abar(sort(j)) * zcarbon))/(1. - exp( - abar(sort(j)) * zcarbon))
            else
              fracarb(i,j,k) = (exp( - abar(sort(j)) * zbotw(i,k) - delzw(i,k)) - &
                               exp( - abar(sort(j)) * zbotw(i,k)))/(1. - exp( - abar(sort(j)) * zcarbon))
            end if
          end do ! loop 245
          solctemp(i,j) = sum(tbar(i,:) * fracarb(i,j,:))/ sum(fracarb(i,j,:))
        end if
      end do ! loop 240
    end do ! loop 230
    !> find moisture scalar for soil c decomposition
    !! this is modelled as function of logarithm of matric potential.
    !! we find values for all soil layers, and then find an average value
    !! based on fraction of carbon present in each layer. this makes
    !! moisture scalar a function of vegetation type.
    do j = 1,ignd ! loop 260
      do i = il1,il2 ! loop 270
        if (isand(i,j) == - 3 .or. isand(i,j) == - 4) then
          scmotrm (i,j) = 0.2
          psi (i,j) = 10000.0 ! set to large number so that
          ! ltrmoscl becomes 0.2
        else ! i.e., sand/=-3 or -4
          !           FLAG- check on this as I had to change a fair amount what YW had, JM. Sep 21 2016.
          !           Also not sure if it is needed?
          !           JM - Turn off for now, we'll see how testing looks. Nov 2016.
          !           EC - Re-implemented as peatland testing shows that in some situations, can get an invalid operation
          !                if thpor+0.005-thice < 0. Note: same approach as in hetresPeat.  Feb 06 2017.
          if (ipeatland(i) > 0) then
            if (thliq(i,j) + thice(i,j) + 0.01 < thpor(i,j) .and. tbar(i,j) < 273.16) then
              psi(i,j) = 0.001
            else if (thice(i,j) > thpor(i,j) ) then
              psi(i,j) = 0.001   ! set to saturation
            else
              psi(i,j) = psisat(i,j) * (thliq(i,j)/(thpor(i,j) - thice(i,j))) ** ( - b(i,j))
            end if
          else
            if (thice(i,j) <= thpor(i,j)) then ! flag new limits
              psi(i,j)   = psisat(i,j) * (thliq(i,j)/(thpor(i,j) - thice(i,j))) ** ( - b(i,j))
            else
              ! if the whole pore space is ice then suction is assumed to be very high.
              psi(i,j) = 10000.0
            end if
          end if
          if (psi(i,j) >= 10000.0) then
            scmotrm(i,j) = 0.2
          else if (psi(i,j) < 10000.0 .and.  psi(i,j) > 6.0) then
            scmotrm(i,j) = 1.0 - 0.8 * ( (log10(psi(i,j)) - log10(6.0))/(log10(10000.0) - log10(6.0)) )
          else if (psi(i,j) <= 6.0 .and. psi(i,j) >= 4.0) then
            scmotrm(i,j) = 1.0
          else if (psi(i,j) < 4.0 .and. psi(i,j) > psisat(i,j) ) then
            scmotrm(i,j) = 1.0 - 0.5 * ( (log10(4.0) - log10(psi(i,j))) /(log10(4.0) - log10(psisat(i,j))) )
          else if (psi(i,j) <= psisat(i,j) ) then
            scmotrm(i,j) = 0.5
          end if
          scmotrm(i,j) = max(0.2,min(1.0,scmotrm(i,j)))
        end if ! sand==-3 or -4
      end do ! loop 270
    end do ! loop 260
    do j = 1,icc ! loop 280
      do i = il1,il2
        if (fcan(i,j) > 0.) then
          socmoscl(i,j) = sum(scmotrm(i,:) * fracarb(i,j,:)) / sum(fracarb(i,j,:))
          socmoscl(i,j) = max(0.2,min(1.0,socmoscl(i,j)))
        end if
      end do ! loop 290
    end do ! loop 280
    !> find moisture scalar for litter decomposition
    !! the difference between moisture scalar for litter and soil c
    !! is that the litter decomposition is not constrained by high
    !! soil moisture (assuming that litter is always exposed to air).
    !! in addition, we use moisture content of the top soil layer
    !! as a surrogate for litter moisture content. so we use only
    !! psi(i,1) calculated in loops 260 and 270 above.
    do i = il1,il2 ! loop 300
      if (ipeatland(i) == 0) then   ! not peatland
        if (psi(i,1) > 10000.0) then
          ltrmoscl(i) = 0.2
        else if (psi(i,1) <= 10000.0 .and. psi(i,1) > 6.0) then
          ltrmoscl(i) = 1.0 -  0.8 * ( (log10(psi(i,1)) - log10(6.0))/(log10(10000.0) - log10(6.0)) )
        else if (psi(i,1) <= 6.0) then
          ltrmoscl(i) = 1.0
        end if
        ltrmoscl(i) = max(0.2,min(1.0,ltrmoscl(i)))
      else  ! is peatland
        !
        !    test psi optimal at psisat
        !    peatland microbals performs better towards wet environment,
        !    for b = 2.3,thpor = 0.98 as soil layer 1,
        !    thliq = 0.01  0.1   0.2    0.3    0.4    0.5   0.6   0.7    0.8     0.9
        !    psi   =  391  1.0   0.38  0.15   0.08   0.05   0.03  0.022  0.016  0.012
        !
        !    set the upper boundary at 500,optimal psi between 0.05 and 0.03
        !    (Mayono et al. 2013)
        !
        !    limit of ltrmoscalms at saturation
        if (psi(i,1) >= 10000.0) then
          ltrmoscl(i) = 0.2
        else if (psi(i,1) <= 10000.0 .and. psi(i,1) > 6.0) then
          ltrmoscl(i) = 1.0 - 0.8 * ((log10(psi(i,1)) - log10(6.0)) &
                        /(log10(10000.0) - log10(6.0))) ** 1.
        else if (psi(i,1) <= 6.0 .and. psi(i,1) > 4.0) then
          ltrmoscl(i) = 1.0
        else if (psi(i,1) <= 4.0 .and. psi(i,1) > psisat(i,1)) then
          ltrmoscl(i) = 1.0 - 0.99 * ((log10(4.0) - log10(psi(i,1)))/ &
                        (log10(4.0) - log10(psisat(i,1))))
        else if (psi(i,1) <= psisat(i,1)) then
          ltrmoscl(i) = 0.01
        end if
        ltrmoscl(i) = max(0.0,min(ltrmoscl(i),1.0))
      end if  ! peatland
    end do ! loop 300
    !!    use temperature of the litter and soil c pools, and their soil
    !!     moisture scalars to find respiration rates from these pools
    do j = 1,icc ! loop 320
      do i = il1,il2
        if (fcan(i,j) > 0.) then
          !         first find the q10 response function to scale base respiration
          !         rate from 15 c to current temperature,we do litter first
          tempq10l(i,j) = litrtemp(i,j) - 273.16
          litrq10 = tanhq10(1) + tanhq10(2) * ( tanh( tanhq10(3) * (tanhq10(4) - tempq10l(i,j))  ) )
          q10func = litrq10 ** (0.1 * (litrtemp(i,j) - 273.16 - 15.0))
          ltresveg(i,j) = ltrmoscl(i) * litrmass(i,j) * bsratelt(sort(j)) * 2.64 * q10func ! 2.64 converts bsratelt from kg c/kg c.year
          ! to u-mol co2/kg c.s
          !         respiration from soil c pool
          tempq10s(i,j) = solctemp(i,j) - 273.16
          soilcq10 = tanhq10(1) + tanhq10(2) * ( tanh( tanhq10(3) * (tanhq10(4) - tempq10s(i,j))  ) )
          q10func = soilcq10 ** (0.1 * (solctemp(i,j) - 273.16 - 15.0))
          scresveg(i,j) = socmoscl(i,j) * soilcmas(i,j) * bsratesc(sort(j)) * 2.64 * q10func  ! 2.64 converts bsratesc from kg c/kg c.year
          ! to u-mol co2/kg c.s
        end if
      end do ! loop 330
    end do ! loop 320

    !   use classicParams,  only : icc,ignd,kk,zero,bsratelt, &
    !                              bsratesc,abar,tanhq10, &
    !                              alpha_hetres,r_depthredu,tcrit,frozered
    !
    !   implicit none
    !
    !   ! Arguments:
    !   integer, intent(in) :: ilg                              !< il1=1
    !   integer, intent(in) :: il1                              !< il1=1
    !   integer, intent(in) :: il2                              !< il2=ilg
    !   integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !                                 !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !                                 !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !                                 !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !                                 !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    !   real, dimension(ilg,icc,ignd),intent(in) :: litrmass   !< litter mass for the 9 pfts + bare [ \f$kg C/m^2\f$ ]
    !   real, dimension(ilg,icc,ignd),intent(in) :: soilcmas   !< soil carbon mass for the 9 pfts + bare [ \f$kg C/m^2\f$ ]
    !   real, dimension(ilg,ignd),intent(in) :: thpor          !< Soil total porosity [ \f$(cm^3 cm^{-3})\f$ ] - daily average
    !   real, dimension(ilg,ignd),intent(in) :: tbar           !< Soil temperature [ K ]
    !   real, dimension(ilg,ignd),intent(in) :: psisat         !< Saturated soil matric potential [ m ]
    !   real, dimension(ilg,ignd),intent(in) :: b              !< Clapp and Hornberger empirical “b” parameter [ ]
    !   real, dimension(ilg,ignd),intent(in) :: thliq          !< liquid soil moisture content in soil layers [ \f$(cm^3 cm^{-3})\f$ ]
    !   real, dimension(ilg,ignd),intent(in) :: thice         !< frozen soil moisture content in 3 soil layers in canopy covered subarea [ \f$(cm^3 cm^{-3})\f$ ]
    !   real, dimension(ilg),intent(in) :: fct                 !< Sum of all fcans [ ]
    !   real, dimension(ilg,icc),intent(in) :: fcan            !< fractional coverage of ctem's 9 pfts [ ]
    !   integer, dimension(ilg,ignd),intent(in) :: isand       !< flag for soil/bedrock/ice/glacier
    !   integer, dimension(icc),intent(in) :: sort             !< index for correspondence between CTEM pfts and all values in the parameters vectors
    !   real, dimension(ilg,ignd),intent(in) :: delzw          !< thickness of the permeable soil layers [ m ]
    !   integer, intent(in) ::  ipeatland(ilg)                  !<
    !   real, intent(in) :: tracerLitrMass(:,:,:)   !< Tracer litter mass for the CTEM pfts + bare in \f$kg c/m^2\f$
    !   real, intent(in) :: tracerSoilCMass(:,:,:)   !< Tracer soil carbon mass for the CTEM pfts + bare in \f$kg c/m^2\f$
    !
    !   real, dimension(ilg,icc,ignd),intent(out) :: ltresveg  !< litter respiration for the given vegetated sub-area [ \f$u-mol co2/m2.sec\f$ ]
    !   real, dimension(ilg,icc,ignd),intent(out) :: scresveg  !< soil carbon respiration for the given vegetated sub-area [ \f$u-mol co2/m2.sec\f$ ]
    !   real, intent(out) :: ltResTracer(ilg,icc,ignd)      !< Tracer litter respiration over the given unvegetated sub-area in umol co2/m2.s
    !   real, intent(out) :: sCResTracer(ilg,icc,ignd)      !< Tracer soil C respiration over the given unvegetated sub-area in umol co2/m2.s
    !
    !
    !   ! Local vars:
    !   integer :: i,j,k
    !   real :: litrq10
    !   real :: soilcq10
    !   real :: q10func
    !   real :: tempq10l
    !   real :: tempq10s
    !   real :: reduceatdepth
    !   real, dimension(ilg,ignd) :: socmoscl
    !   real, dimension(ilg,ignd) :: psi
    !   real, dimension(ilg,ignd) :: ltrmoscl
    !
    !   ! Initialize required arrays to zero
    !   socmoscl(:,:)=0.0      ! soil carbon moisture term
    !   ltrmoscl(:,:)=0.0     ! soil moisture scalar for litter decomposition
    !   ltresveg(:,:,:)=0.0   ! litter resp. rate for each pft
    !   scresveg(:,:,:)=0.0   ! soil c resp. rate for each pft
    !   ltResTracer = 0.
    !   sCResTracer = 0.
    !
    !   !! Find moisture scalar for soil c decomposition
    !   !! This is modelled as function of logarithm of matric potential.
    !   !! we find values for all soil layers, and then find an average value
    !   !! based on fraction of carbon present in each layer. this makes
    !   !! moisture scalar a function of vegetation type.
    !
    !   do j = 1,ignd
    !     do i = il1,il2
    !       if (isand(i,j)==-3.or.isand(i,j)==-4) then
    !         socmoscl(i,j)=0.2
    !         psi (i,j) = 10000.0 ! set to large number so that ltrmoscl becomes 0.2
    !       else ! i.e., sand/=-3 or -4
    !
    !         if (ipeatland(i) > 0) then
    !           if (thliq(i,j) + thice(i,j) + 0.01 < thpor(i,j) .and. tbar(i,j) < 273.16) then
    !             psi(i,j) = 0.001
    !           else if (thice(i,j) > thpor(i,j)) then
    !             psi(i,j) = 0.001   ! set to saturation
    !           else
    !             psi(i,j) = psisat(i,j) * (thliq(i,j) / (thpor(i,j) - thice(i,j)))**(-b(i,j))
    !           end if
    !         else
    !           if (thice(i,j) <= thpor(i,j)) then ! flag new limits
    !             psi(i,j) = psisat(i,j) * (thliq(i,j) / (thpor(i,j) - thice(i,j)))**(-b(i,j))
    !           else
    !             ! if the whole pore space is ice then suction is assumed to be very high.
    !             psi(i,j) = 10000.0
    !           end if
    !         end if
    !
    !         if (psi(i,j) >= 10000.0) then
    !           socmoscl(i,j) = 0.2
    !         else if (psi(i,j) < 10000.0 .and.  psi(i,j) > 6.0) then
    !           socmoscl(i,j) = 1.0 - 0.8 * ((log10(psi(i,j)) - log10(6.0)) &
    !                                      / (log10(10000.0) - log10(6.0)))
    !         else if (psi(i,j) <= 6.0 .and. psi(i,j) >= 4.0) then
    !           socmoscl(i,j) = 1.0
    !         else if (psi(i,j) < 4.0 .and. psi(i,j) > psisat(i,j)) then
    !           socmoscl(i,j) = 1.0 -0.5 * ((log10(4.0) - log10(psi(i,j))) &
    !                                    / (log10(4.0) - log10(psisat(i,j))))
    !         else if (psi(i,j) <= psisat(i,j)) then
    !           socmoscl(i,j) = 0.5
    !         end if
    !       end if ! sand==-3 or -4
    !
    !       socmoscl(i,j) = max(0.2,min(1.0,socmoscl(i,j)))
    !
    ! 270     continue
    ! 260   continue
    !
    !   !! Find moisture scalar for litter decomposition
    !   !! the difference between moisture scalar for litter and soil c
    !   !! is that the litter decomposition is not constrained by high
    !   !! soil moisture (assuming that litter is always exposed to air).
    !   !! in addition, we use moisture content of the top soil layer
    !   !! as a surrogate for litter moisture content. so we use only
    !   !! psi(i,1) calculated in loops 260 and 270 above.
    !
    !   !  FLAG right now assume that the first soil layer litter behaves like this.
    !   !  the layers below are more impeded by soil moisture (same as soil C). JM Feb 8 2016.
    !   do i = il1,il2
    !     if (ipeatland(i) == 0) then   ! not peatland
    !       if (psi(i,1) > 10000.0) then
    !         ltrmoscl(i,1) = 0.2
    !       else if (psi(i,1) <= 10000.0 .and. psi(i,1) > 6.0) then
    !         ltrmoscl(i,1) = 1.0 - 0.8 *((log10(psi(i,1)) - log10(6.0)) &
    !                                   / (log10(10000.0) - log10(6.0)))
    !       else if (psi(i,1) <= 6.0) then
    !         ltrmoscl(i,1) = 1.0
    !       end if
    !       ltrmoscl(i,1) = max(0.2,min(1.0,ltrmoscl(i,1)))
    !     else  ! is peatland
    !       !  FLAG per layer implementation untested !! JM Mar 14 2019.
    !       !    test psi optimal at psisat
    !       !    peatland microbals performs better towards wet environment,
    !       !    for b = 2.3, thpor = 0.98 as soil layer 1,
    !       !    thliq = 0.01  0.1   0.2    0.3    0.4    0.5   0.6   0.7    0.8     0.9
    !       !    psi   =  391  1.0   0.38  0.15   0.08   0.05   0.03  0.022  0.016  0.012
    !       !
    !       !    set the upper boundary at 500,optimal psi between 0.05 and 0.03
    !       !    (Mayono et al. 2013)
    !       !
    !       !    limit of ltrmoscalms at saturation
    !       if (psi(i,1) >= 10000.0) then
    !         ltrmoscl(i,1) = 0.2
    !       else if (psi(i,1) <= 10000.0 .and. psi(i,1) > 6.0) then
    !         ltrmoscl(i,1) =1.0 - 0.8 * ((log10(psi(i,1)) - log10(6.0)) &
    !                                   / (log10(10000.0) - log10(6.0)))**1.
    !       else if (psi(i,1) <= 6.0 .and. psi(i,1) > 4.0) then
    !         ltrmoscl(i,1) = 1.0
    !       else if (psi(i,1) <= 4.0 .and. psi(i,1) > psisat(i,1)) then
    !         ltrmoscl(i,1) = 1.0 - 0.99*((log10(4.0) - log10(psi(i,1))) &
    !                                   / (log10(4.0) - log10(psisat(i,1))))
    !       else if (psi(i,1) <= psisat(i,1)) then
    !         ltrmoscl(i,1) = 0.01
    !       end if
    !       ltrmoscl(i,1) = max(0.0,min(ltrmoscl(i,1),1.0))
    !     end if  ! peatland
    ! 300   continue
    !
    !   !> Set the lower levels to have the same moisture sensitivity to soil C.
    !   ltrmoscl(:,2:ignd)=socmoscl(:,2:ignd)
    !
    !   !< Use temperature of the litter and soil c pools, and their soil
    !   !! moisture scalars to find respiration rates from these pools
    !
    !   do j = 1,icc
    !     do i = il1,il2
    !       do k= 1,ignd
    !         if (fcan(i,j) > 0.) then
    !
    !           ! First find the q10 response function to scale base respiration
    !           ! rate from 15 c to current temperature, we do litter first
    !
    !           tempq10l = tbar(i,k) - 273.16
    !           litrq10 = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - tempq10l)))
    !
    !           !! Apply a step reduction in q10func when the soil layer freezes. This is to reflect the
    !           !! lost mobility of the population due to less liquid water
    !           if (tbar(i,k) - 273.16 > tcrit) then
    !               q10func = litrq10**(0.1 * (tbar(i,k) - 273.16 - 15.0))
    !           else
    !               q10func = litrq10**(0.1 * (tbar(i,k) - 273.16 - 15.0)) * frozered
    !           end if
    !
    !           !! Reduce the respiration at depth due to unresolved depth dependent processes including
    !           !! soil microbial population dynamics, pore-scale oxygen availability, mineral sorption,
    !           !! priming effects, and other unrepresented processes. This is following Lawrence et al.
    !           !! Enviro Res Lett 2015. We apply this for all soils.
    !           reduceatdepth = exp(-delzw(i,k) / r_depthredu)
    !
    !           ! 2.64 converts bsratelt from kg c/kg c.year to u-mol co2/kg c.s
    !           ltresveg(i,j,k) = ltrmoscl(i,k) * litrmass(i,j,k) * bsratelt(sort(j)) * 2.64 &
    !                                           * q10func * reduceatdepth
    !
    !           if (useTracer > 0) ltResTracer(i,j,k) = ltrmoscl(i,k) * tracerLitrMass(i,j,k) &
    !                                                           * bsratelt(sort(j)) * 2.64 * q10func &
    !                                                           * reduceatdepth
    !
    !           ! Respiration from soil c pool
    !           tempq10s = tbar(i,k) - 273.16
    !           soilcq10 = tanhq10(1) + tanhq10(2) *(tanh(tanhq10(3) * (tanhq10(4) - tempq10s)))
    !
    !           if (tbar(i,k) - 273.16 > tcrit) then
    !             q10func = soilcq10**(0.1 * (tbar(i,k) - 273.16 - 15.0))
    !           else
    !             q10func = soilcq10**(0.1 * (tbar(i,k) - 273.16 - 15.0)) * frozered
    !           end if
    !           ! 2.64 converts bsratesc from kg c/kg c.year to u-mol co2/kg c.s
    !           scresveg(i,j,k) = socmoscl(i,k) * soilcmas(i,j,k) * bsratesc(sort(j)) * 2.64 &
    !                                           * q10func * reduceatdepth
    !
    !           if (useTracer > 0) sCResTracer(i,j,k) = socmoscl(i,k) * tracerSoilCMass(i,j,k) &
    !                                                           * bsratesc(sort(j)) * 2.64 * q10func &
    !                                                           * reduceatdepth
    !
    !         end if
    ! 340   continue
    ! 330 continue
    ! 320 continue
    ! COMBAK PERLAY

    return

  end subroutine hetresv
  !! @}

  ! ---------------------------------------------------------------

  !> \ingroup heterotrophicrespirationmod_updatePoolsHetResp
  !! @{ Find vegetation and tile averaged litter and soil C respiration rates
  !! using values from canopy over ground and canopy over snow subareas.
  !! Also adds the moss and peat soil respiration to the tile level quantities.
  !! Next the litter and soil C pools are updated based on litter and soil C respiration rates.
  !! The humidified litter is then transferred to the soil C pool.
  !! Soil respiration is estimated as the sum of heterotrophic respiration and root maintenance respiration.
  !! For peatlands, we additionally add moss values to the grid (litter respiration
  !! and moss root respiration).
  !> @author Vivek Arora, Joe Melton
  subroutine updatePoolsHetResp (il1, il2, ilg, fcancmx, ltresveg, scresveg, & ! In
                                 ipeatland, fg, litresmoss, socres_peat, & ! In
                                 sort, spinfast, rmrveg, rmr, leapnow, & ! In
                                 useTracer, ltResTracer, sCResTracer, litResMossTracer, soCResPeatTracer, & ! In
                                 litrmass, soilcmas, Cmossmas, & ! In / Out
                                 tracerLitrMass, tracerSoilCMass, tracerMossCMass, & ! In / Out
                                 hetrsveg, litres, socres, hetrores, humtrsvg, soilresp, & ! Out
                                 humiftrs, hutrstep_g, litrfallmoss, ltrestepmoss, &
                                 humstepmoss, socrestep) ! Out

    use classicParams, only : icc, iccp1, iccp2, ignd, zero, humicfac, humicfac_bg, &
                              deltat, humicfacmoss, rmortmoss

    implicit none

    ! Inputs
    integer, intent(in) :: il1             !< il1=1
    integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg
    integer, intent(in) ::  spinfast        !< spinup factor for soil carbon whose default value is 1. as this factor increases the
    !< soil c pool will come into equilibrium faster. reasonable value for spinfast is
    !< between 5 and 10. when spinfast/=1 then the balcar subroutine is not run.
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    integer, intent(in) :: ipeatland(ilg)       !< Peatland flag: 0 = not a peatland, 1 = bog, 2 = fen
    logical, intent(in) :: leapnow            !< true if this year is a leap year. Only used if the switch 'leap' is true.
    real, intent(in) :: fcancmx(ilg,icc)      !< max. fractional coverage of CTEM's pfts, but this can be
    !< modified by land-use change, and competition between pfts
    ! COMBAK PERLAY
    real, intent(in) :: ltresveg(ilg,iccp2)     !< fluxes for each pft: litter respiration for each pft + bare fraction
    real, intent(in) :: scresveg(ilg,iccp2)     !< soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
    ! real, intent(in) :: ltresveg(ilg,iccp2,ignd)     !< fluxes for each pft: litter respiration for each pft + bare fraction
    ! real, intent(in) :: scresveg(ilg,iccp2,ignd)     !< soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
    ! COMBAK PERLAY
    real, intent(in) :: fg(ilg)             !< Fraction of grid cell that is bare ground.
    integer, intent(in) :: sort(icc)
    real, intent(in) :: litresmoss(ilg)    !< moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: socres_peat(ilg)   !< heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: rmrveg(ilg,icc)    !< Maintenance respiration for root for the CTEM pfts in u mol co2/m2. sec
    real, intent(in) :: rmr(ilg)           !< Root maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: ltResTracer(ilg,iccp2,ignd)  !< Tracer fluxes for each pft: litter respiration for each pft + bare fraction
    real, intent(in) :: sCResTracer(ilg,iccp2,ignd)  !< Tracer soil carbon respiration for the given sub-area in umol co2/m2.s,for ctem's pfts
    real, intent(in) :: litResMossTracer(ilg)    !< Tracer moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$) TODO: Units
    real, intent(in) :: soCResPeatTracer(ilg)    !< Tracer heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    ! Updates
    ! COMBAK PERLAY
    real, intent(inout) :: litrmass (ilg,iccp2)  !< Litter mass for each of the CTEM pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, intent(inout) :: soilcmas(ilg,iccp2)   !< Soil carbon mass for each of the CTEM pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    ! real, intent(inout) :: litrmass (ilg,iccp2,ignd)  !< Litter mass for each of the CTEM pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    ! real, intent(inout) :: soilcmas(ilg,iccp2,ignd)   !< Soil carbon mass for each of the CTEM pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    ! COMBAK PERLAY
    real, intent(inout) :: Cmossmas(ilg)         !< Moss biomass C pool (kgC/m2)

    real, intent(inout) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(inout) :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(inout) :: tracerMossCMass(:)      !< Tracer mass in moss biomass, \f$kg C/m^2\f$

    ! Outputs
    real, intent(out) :: hutrstep_g(ilg)    !< Tile sum of humification from vascualr litter (kgC/m2/timestep)
    real, intent(out) :: litrfallmoss(ilg)  !< moss litter fall (kgC/m2/timestep)
    real, intent(out) :: ltrestepmoss(ilg)  !< litter respiration from moss (kgC/m2/timestep)
    real, intent(out) :: humstepmoss(ilg)   !< moss humification (kgC/m2/timestep)
    real, intent(out) :: socrestep(ilg)     !< heterotrophic respiration from soil (kgC/m2/timestep)
    real, intent(out) :: hetrsveg(ilg,iccp1) !< Vegetation averaged litter and soil C respiration rates (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: litres(ilg)      !< Litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: socres(ilg)      !< Soil carbon respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: hetrores(ilg)    !< Heterotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    ! COMBAK PERLAY
    real, intent(out) :: humtrsvg(ilg,iccp2)     !< transfer of humidified litter from litter to soil c pool per PFT. (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    ! real, intent(out) :: humtrsvg(ilg,iccp2,ignd)     !< transfer of humidified litter from litter to soil c pool per PFT. (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    ! COMBAK PERLAY
    real, intent(out) :: soilresp(ilg)    !< Soil respiration. This includes root respiration and respiration from
    !! litter and soil carbon pools. Note that soilresp is different from
    !! socres,which is respiration from the soil C pool.(\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: humiftrs(ilg)    !< Transfer of humidified litter from litter to soil C pool (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    ! Local
    integer :: k, j, i
    real :: ltrestep !<
    real :: screstep !<
    ! COMBAK PERLAY
    real :: hutrstep(ilg,iccp2) !<
    ! real :: hutrstep(ilg,iccp2,ignd) !<
    ! COMBAK PERLAY
    real :: tracerhutrstep !<
    real :: soilrsvg(ilg,iccp2) !<

    ! -----------

    !> Find vegetation averaged litter and soil C respiration rates
    !! using values from canopy over ground and canopy over snow subareas
    hetrsveg(:,:) = 0.0
    do j = 1,icc ! loop 340
      do i = il1,il2
        if (fcancmx(i,j) > zero) then
          ! COMBAK PERLAY
          hetrsveg(i,j) =  hetrsveg(i,j) + ltresveg(i,j) + scresveg(i,j)
          ! do k = 1,ignd
          !   ! hetrsveg is kept per PFT and tile (not per layer) at the moment.
          !   hetrsveg(i,j) =  hetrsveg(i,j) + ltresveg(i,j,k) + scresveg(i,j,k)
          ! end do
          ! COMBAK PERLAY
        end if
      end do ! loop 350
    end do ! loop 340

    !> Find litter and soil c respiration rates averaged over the bare
    !! fraction of the grid cell using values from ground and snow over ground sub-areas.
    do i = il1,il2 ! loop 355
      if (fg(i) > zero) then
        ! COMBAK PERLAY
        hetrsveg(i,iccp1) = hetrsveg(i,iccp1) + ltresveg(i,iccp1) + scresveg(i,iccp1)
        ! do k = 1,ignd
        !   ! hetrsveg is kept per PFT and tile (not per layer) at the moment.
        !   hetrsveg(i,iccp1) = hetrsveg(i,iccp1) + ltresveg(i,iccp1,k) + scresveg(i,iccp1,k)
        ! end do
        ! COMBAK PERLAY
        ! nepveg(i,iccp1)=0.-hetrsveg(i,iccp1)
      end if
    end do ! loop 355

    !> Find grid averaged litter and soil c respiration rates
    litres(:) = 0.
    socres(:) = 0.
    do j = 1,icc ! loop 360
      do i = il1,il2
        ! COMBAK PERLAY
        litres(i) = litres(i) + fcancmx(i,j) * ltresveg(i,j)
        socres(i) = socres(i) + fcancmx(i,j) * scresveg(i,j)
        ! do k = 1,ignd
        !   litres(i) = litres(i) + fcancmx(i,j) * ltresveg(i,j,k)
        !   socres(i) = socres(i) + fcancmx(i,j) * scresveg(i,j,k)
        ! end do
        ! COMBAK PERLAY
      end do ! loop 370
    end do ! loop 360

    !> Add moss and peat soil respiration to the grid
    !! In loop 360/370 we have aggregated across all pfts for litres and socres, In loop 380 we add
    !! bareground (iccp1) and LUC pool (iccp2) values to the grid sum if it's not peatland.
    !! If it is a peatland, we add litresmoss to
    !! the grid sum but no bareground values as we assume peatlands have no bareground.
    do i = il1,il2 ! loop 380
      if (ipeatland(i) == 0) then
        ! COMBAK PERLAY
        litres(i) = litres(i) + fg(i) * ltresveg(i,iccp1)
        socres(i) = socres(i) + fg(i) * scresveg(i,iccp1)
        ! do k = 1,ignd
        !   litres(i) = litres(i) + fg(i) * ltresveg(i,iccp1,k)
        !   socres(i) = socres(i) + fg(i) * scresveg(i,iccp1,k)
        ! end do
        ! COMBAK PERLAY
      else  ! peatlands
        litres(i) = litres(i) + litresmoss(i) ! add the moss litter, which is assumed to cover whole tile.
        !
        ! CAUTION says Vivek
        ! Note that the following line overwrites socres(i) calculated above in loop 370, although
        ! socres(i) is based on scresveg(i,j) which isn't modified. Similarly, the implementation
        ! of peatlands also means hetrores(i), calculated below, is now inconsistent with hetrsveg(i,j).
        ! The implementation of peatlands is based on overwriting variables without making things
        ! consistent. Of course, overwriting variables is not advised because it makes things confusing.
        !
        socres(i) = socres_peat(i) ! since this is only peatland on this tile, use just the peat value.

      end if

      hetrores(i) = litres(i) + socres(i)

    end do ! loop 380

    !> Update the litter and soil c pools based on litter and soil c respiration rates
    !! found above. also transfer humidified litter to the soil c pool.
    do j = 1,iccp2 ! loop 420
      do i = il1,il2
        ! COMBAK PERLAY
        !> Convert u mol co2/m2.sec -> \f$(kg C/m^2)\f$ respired over the model time step
        ltrestep = ltresveg(i,j) * deltat / 963.62
        screstep = scresveg(i,j) * deltat / 963.62

        !> Update litter and soil c pools
        if (j < iccp1) then
          litrmass(i,j) = litrmass(i,j) - (ltrestep * (1.0 + humicfac(sort(j))))
          hutrstep(i,j) = humicfac(sort(j)) * ltrestep
        else
          !> Next we add bareground and LUC pool litter mass and humification for non-peatlands.
          if (ipeatland(i) == 0) then
            litrmass(i,j) = litrmass(i,j) - (ltrestep * (1.0 + humicfac_bg))
            hutrstep(i,j) = humicfac_bg * ltrestep
            ! else for peatlands:
            ! In peatlands there is no bareground litter mass since it is the moss layer.
          end if
        end if

        humtrsvg(i,j) = hutrstep(i,j) * (963.62 / deltat) ! u-mol co2/m2.sec

        soilcmas(i,j) = soilcmas(i,j) + real(spinfast) * (hutrstep(i,j) -  screstep)

        if (litrmass(i,j) < zero) litrmass(i,j) = 0.0
        if (soilcmas(i,j) < zero) soilcmas(i,j) = 0.0

        !       do k = 1,ignd
        !
        !         !> Convert u mol co2/m2.sec -> \f$(kg C/m^2)\f$ respired over the model time step
        !         ltrestep = ltresveg(i,j,k) * deltat / 963.62
        !         screstep = scresveg(i,j,k) * deltat / 963.62
        !
        !         !> Update litter and soil c pools
        !         if (j < iccp1) then
        !           litrmass(i,j,k) = litrmass(i,j,k) - (ltrestep * (1.0 + humicfac(sort(j))))
        !           hutrstep(i,j,k) = humicfac(sort(j)) * ltrestep
        !         else
        !           !> Next we add bareground and LUC pool litter mass and humification for non-peatlands.
        !           if (ipeatland(i) == 0) then
        !             litrmass(i,j,k) = litrmass(i,j,k) - (ltrestep * (1.0 + humicfac_bg))
        !             hutrstep(i,j,k) = humicfac_bg * ltrestep
        !           ! else for peatlands:
        !             ! In peatlands there is no bareground litter mass since it is the moss layer.
        !           end if
        !         end if
        !
        !         humtrsvg(i,j,k) = hutrstep(i,j,k) * (963.62 / deltat) ! u-mol co2/m2.sec
        !
        !         soilcmas(i,j,k) = soilcmas(i,j,k) + real(spinfast) * (hutrstep(i,j,k) -  screstep)
        !
        !         if (litrmass(i,j,k) < zero) litrmass(i,j,k) = 0.0
        !         if (soilcmas(i,j,k) < zero) soilcmas(i,j,k) = 0.0
        !
        !         ! Now perform the same calculations to the tracer pools
        !         if (useTracer > 0) then
        !           if (j < iccp1) then ! vegetated
        !             tracerLitrMass(i,j,k) = tracerLitrMass(i,j,k) - (ltResTracer(i,j,k) * deltat / 963.62 &
        !                                                             * (1.0 + humicfac(sort(j))))
        !             tracerhutrstep = humicfac(sort(j)) * ltResTracer(i,j,k) * deltat / 963.62
        !           else ! bare ground
        !             if (ipeatland(i) == 0) then ! non-peatlands
        !               tracerLitrMass(i,j,k) = tracerLitrMass(i,j,k) - (ltResTracer(i,j,k) * deltat / 963.62 &
        !                                                              * (1.0 + humicfac_bg))
        !               tracerhutrstep = humicfac_bg * ltResTracer(i,j,k) * deltat / 963.62
        !             end if
        !           end if
        !
        !           tracerSoilCMass(i,j,k) = tracerSoilCMass(i,j,k) + real(spinfast) &
        !                                           * (tracerhutrstep -  sCResTracer(i,j,k) * deltat / 963.62)
        !
        !           if (litrmass(i,j,k) < zero) tracerLitrMass(i,j,k) = 0.0
        !           if (soilcmas(i,j,k) < zero) tracerSoilCMass(i,j,k) = 0.0
        !         end if
        !
        ! 435   continue
        ! COMBAK PERLAY
      end do ! loop 430
    end do ! loop 420

    !> Estimate soil respiration. this is sum of heterotrophic respiration and root maintenance respiration.
    soilrsvg(:,:) = 0.
    do j = 1,icc ! loop 440
      do i = il1,il2
        ! COMBAK PERLAY
        soilrsvg(i,j) = soilrsvg(i,j) + ltresveg(i,j) + scresveg(i,j)
        !       do k = 1,ignd
        !         ! soilrsvg kept as per pft/per tile for now (not per layer)
        !         soilrsvg(i,j) = soilrsvg(i,j) + ltresveg(i,j,k) + scresveg(i,j,k)
        ! 450   continue
        ! COMBAK PERLAY
        soilrsvg(i,j) = soilrsvg(i,j) + rmrveg(i,j)
      end do ! loop 445
    end do ! loop 440

    !> But over the bare fraction and LUC product pool there is no live root.

    do i = il1,il2 ! loop 460
      ! COMBAK PERLAY
      soilrsvg(i,iccp1) = soilrsvg(i,iccp1) + ltresveg(i,iccp1) + scresveg(i,iccp1)
      !     do k = 1,ignd
      !       soilrsvg(i,iccp1) = soilrsvg(i,iccp1) + ltresveg(i,iccp1,k) + scresveg(i,iccp1,k)
      ! 465 continue
      ! COMBAK PERLAY
    end do ! loop 460

    !> Find grid averaged humification and soil respiration rates
    soilresp(:) = 0.0
    humiftrs(:) = 0.0
    hutrstep_g(:) = 0.0
    do i = il1,il2 ! loop 470
      do j = 1,icc
        soilresp(i) = soilresp(i) + fcancmx(i,j) * soilrsvg(i,j)
        ! COMBAK PERLAY
        hutrstep_g(i) = hutrstep_g(i) + fcancmx(i,j) * hutrstep(i,j)
        humiftrs(i) = humiftrs(i) + fcancmx(i,j) * humtrsvg(i,j)
        ! do k = 1,ignd
        !   hutrstep_g(i) = hutrstep_g(i) + fcancmx(i,j) * hutrstep(i,j,k)
        !   humiftrs(i) = humiftrs(i) + fcancmx(i,j) * humtrsvg(i,j,k)
        ! end do
        ! COMBAK PERLAY
      end do ! loop 480

      !> After aggregation of humification and soil respiration rates for non-peatlands
      !! to the grid/tile level, the same must be done for bareground (iccp1).
      !! For peatlands, we additionally add moss values to the grid (litter respiration
      !! and moss root respiration). Note in loop 430 iccp1 is passed for peatlands

      if (ipeatland(i) == 0) then ! non peatland

        soilresp(i) = soilresp(i) + fg(i) * soilrsvg(i,iccp1)
        ! COMBAK PERLAY
        humiftrs(i) = humiftrs(i) + fg(i) * humtrsvg(i,iccp1)
        ! do k = 1,ignd
        !   humiftrs(i) = humiftrs(i) + fg(i) * humtrsvg(i,iccp1,k)
        ! end do
        ! COMBAK PERLAY

        ! Set all peatland vars to 0.
        litrfallmoss(i) = 0.
        ltrestepmoss(i) = 0.
        humstepmoss(i) = 0.
        socrestep(i) = 0.

      else ! peatland

        ! CAUTION says Vivek
        ! Here again soilresp(i) is overwritten with socres(i)=socres_peat(i) as calculated
        ! above in loop 380. This makes soilresp(i) inconsistent with soilrsvg(i,j) for
        ! peatland gridcells/tile.

        soilresp(i) = socres(i) + litres(i) + rmr(i) ! moss root and litter respiration. No bareground !

        ! Calculate moss time step C fluxes, '/365*deltat' convert year-1
        ! to time step-1, 'deltat/963.62' convert umol CO2/m2/s to kgC/m2/deltat.
        ! note that hutrstep_g aggregation for icc was done in loop 480
        if (leapnow) then
          litrfallmoss(i) = Cmossmas(i) * rmortmoss / 366. * deltat ! kgC/m2/day(dt)
          ! if (useTracer > 0) FLAG,not connected up !
        else
          litrfallmoss(i) = Cmossmas(i) * rmortmoss / 365. * deltat ! kgC/m2/day(dt)
        end if
        ltrestepmoss(i) = litresmoss(i) * (1.0 / 963.62) * deltat   ! kgC/m2/dt
        socrestep(i) = socres(i) * (1.0 / 963.62) * deltat     ! kgC/m2/dt
        soilresp(i) = soilresp(i) * (1.0 / 963.62) * deltat   ! kgC/m2/dt
        humstepmoss(i) = humicfacmoss * ltrestepmoss(i)        ! kgC/m2/dt
        hutrstep_g(i) = hutrstep_g(i) + humstepmoss(i)     ! kgC/m2/dt
        humiftrs(i)  = humiftrs(i) + humstepmoss(i) * (963.62 / deltat)! umol/m2/s

      end if
    end do ! loop 470

  end subroutine updatePoolsHetResp
  !! @}
  ! ------------------------------------------------------------------------------
  !> \namespace heterotrophicrespirationmod
  !! Heterotrophic Respiration Module (Vegetated and Bare Ground)
  !!
  !! Central module for all heterotrophic respiration-related operations
  !!
  !! Heterotrophic respiration, \f$R_\mathrm{h}\f$ (\f$mol\, CO_2\, m^{-2}\, s^{-1}\f$), in CTEM is
  !! based on respiration from the litter (which includes contributions from the stem, leaf
  !! and root components), \f$R_{h, D}\f$, and soil carbon, \f$R_{h, H}\f$, pools,
  !!
  !! \f[ \label{hetres_all} R_\mathrm{h}=R_{h, D}+R_{h, H}. \hspace{10pt}[Eqn 1] \f]
  !!
  !! Heterotrophic respiration is regulated by soil temperature and moisture and is
  !! calculated on a daily time step. The original heterotrophic respiration scheme is
  !! described in Arora (2003) \cite Arora2003-3b7 while the modified parametrization used in CTEM v. 2.0
  !! is detailed in Melton and Arora (2014) \cite Melton2014-xy.
  !! Respiration from the litter and soil carbon pools takes the following basic form
  !!
  !! \f[ R_{\mathrm{h}, i} = 2.64 \times 10^{-6}\, \varsigma_i C_i f_{15}(Q_{10}) f(\Psi)_i f(z),
  !! \nonumber \\ i = \mathrm{D}, \mathrm{H}. \hspace{10pt}[Eqn 2] \f]
  !!
  !! The soil carbon and litter respiration depends on the amount of carbon in these components
  !! (\f$C_\mathrm{H}\f$ and \f$C_\mathrm{D}\f$; \f$kg\, C\, m^{-2}\f$) and on a PFT-dependent
  !! respiration rate specified at \f$15\, {C}\f$ (\f$\varsigma_\mathrm{H}\f$ and
  !! \f$\varsigma_\mathrm{D}\f$; \f$kg\, C\, (kg\, C)^{-1}\, yr^{-1}\f$; see also
  !! classicParams.f90). The constant \f$2.64 \times 10^{-6}\f$ converts units from
  !! \f$kg\, C\, m^{-2}\, yr^{-1}\f$ to \f$mol\, CO_2\, m^{-2}\, s^{-1}\f$.
  !!
  !! The effect of soil moisture is accounted for via dependence on soil matric
  !! potential (\f$f(\Psi)\f$), described later. The temperature dependency of
  !! microbial soil respiration rates has been estimated by several different
  !! formulations, ranging from simple \f$Q_{10}\f$ (exponential) to Arrhenius-type
  !! formulations (see review by Lloyd and Taylor (1994) \cite Lloyd1994-ct). In CTEM, soil temperature
  !! influences heterotrophic respiration through a temperature-dependent
  !! \f$Q_{10}\f$ function (\f$f_{15}(Q_{10})\f$). The value of \f$Q_{10}\f$
  !! itself is assumed to be a function of temperature following a hyperbolic
  !! tan function:
  !!
  !! \f[ Q_{10} = 1.44 + 0.56\, \tanh[0.075 (46.0 - T_i)], \nonumber\\ i
  !! = \mathrm{D}, \mathrm{H}, \hspace{10pt}[Eqn 3]\f]
  !!
  !! where \f$T_{\{D, H\}}\f$ is the temperature of either the litter or soil
  !! carbon pool (\f$C\f$), respectively. The parametrization is a compromise
  !! between the temperature-independent \f$Q_{10}\f$ commonly found in many
  !! terrestrial ecosystem models (Cox, 2011) \cite Cox2001-am and the temperature-dependent
  !! \f$Q_{10}\f$ of Kirschbaum (1995) \cite Kirschbaum1995-db. While a constant \f$Q_{10}\f$ yields
  !! an indefinitely increasing respiration rate with increasing temperature, the
  !! formulation of Kirschbaum (1995) \cite Kirschbaum1995-db gives a continuously increasing
  !! \f$Q_{10}\f$ under decreasing temperature, which leads to unreasonably high
  !! soil and litter carbon pools at high latitudes in CTEM. The CTEM
  !! parametrization avoids these issues with a \f$Q_{10}\f$ value of about 2.0
  !! for temperatures less than \f$20\, C\f$, while a decreasing value of
  !! \f$Q_{10}\f$ at temperatures above \f$20\, C\f$ ensures that the
  !! respiration rate does not increase indefinitely. As the soil temperature decreases below
  !! \f$ T_{crit} \f$ (typically \f$1\, C\f$) a step function is applied to the \f$f_{15}(Q_{10})\f$
  !! function to reflect lost mobility of the microbial populations due to less liquid water as:
  !! \f[ f_{15}(Q_{10, (T_i < T_{crit})}) = f_{15}(Q_{10}) * 0.1 \f]
  !!
  !!   \image html "Q10_response_sm.png" "Q10 response"
  !!
  !! The soil detrital pools are explictly tracked per soil layer.
  !!
  !! The response of heterotrophic respiration to soil moisture is formulated through
  !! soil matric potential (\f$\Psi\f$; \f$MPa\f$). While soil matric potential values
  !! are usually negative, the formulation uses absolute values to allow its logarithm
  !! to be taken. Absolute values of soil matric potential are high when soil is dry
  !! and low when it is wet. The primary premise of soil moisture control on heterotrophic
  !! respiration is that heterotrophic respiration is constrained both when the soils
  !! are dry (due to reduced microbial activity) and when they are wet (due to impeded
  !! oxygen supply to microbes) with optimum conditions in-between. The exception is the
  !! respiration from the litter component of the first soil layer, which is assumed to be continually exposed
  !! to air, and thus never oxygen deprived, even when soil moisture content is high
  !! (\f$0.04 > \vert \Psi \vert \geq \vert \Psi_{sat} \vert\f$, where \f$\Psi_{sat}\f$
  !! is the soil matric potential at saturation). The soil moisture dependence for each
  !! soil layer thus varies between 0 and 1 with matric potential as follows:
  !!
  !! for \f$0.04 > \vert\Psi_i\vert \geq \vert\Psi_{sat, i}\vert\f$
  !!
  !! \f[ f(\Psi_i)_\mathrm{H, (D, i>1)} = 1 - 0.5  \frac{\log(0.04)-\log\vert\Psi_i\vert}
  !! {\log(0.04)-\log\vert\Psi_{sat, i}\vert} \hspace{10pt}[Eqn 4]\f]
  !!
  !! \f[f(\Psi_i)_{D, i=1} = 1\nonumber; \hspace{10pt}[Eqn 5]\f]
  !!
  !! for \f$0.06 \geq \vert\Psi_i\vert \geq 0.04\f$
  !! \f[ f(\Psi_i)_{D, H} = 1; \hspace{10pt}[Eqn 6]\f]
  !!
  !! for \f$100.0 \geq \vert\Psi_i\vert > 0.06\f$
  !! \f[ f(\Psi_i)_{D, H} = 1 - 0.8\frac{\log\vert\Psi_i\vert-\log(0.06)}{\log(100)-\log(0.06)}; \hspace{10pt}[Eqn 7]\f]
  !!
  !! for \f$\vert\Psi_i\vert > 100.0\f$
  !! \f[ f(\Psi_i)_{D, H}=0.2. \hspace{10pt}[Eqn 8]\f]
  !!
  !! Respiration also is reduced at depth in soil (\f$f(z)\f$) following Lawrence et al. (2015)
  !! \cite Lawrence2015-tj. This term is meant to represent unresolved depth dependent soil
  !! processes (such as oxygen availability, microbial community changes, etc.). The reduction
  !! in respiration with depth per soil layer is dependent upon the layer depth and
  !! a term, \f$z_t\f$, which is given a value of 10.0 (see ctem_params.f90) as,
  !!
  !! \f[ f(z_i) =\exp (-z_i / z_t) \hspace{10pt}[Eqn 9]\f]
  !!
  !!   \image html "decr_resp_wit_depth_sm.png" "Decrease in respiration with depth"
  !!
  !! Heterotrophic respiration for bare ground is treated separately in CTEM. The carbon
  !! contributions to the bare ground litter and soil carbon pools come via processes
  !! such as creation of bare ground due to fire, competition between PFTs and land use
  !! change. The heterotrophic respiration is sensitive to temperature and moisture in
  !! the same manner as vegetated areas using Eqs. (2)--(8). The
  !! base respiration rates of \f$\varsigma_{D, bare}\f$ and \f$\varsigma_{H, bare}\f$ are
  !! set to \f$0.5605\f$ and \f$0.02258\, kg\, C\, (kg\, C)^{-1}\, yr^{-1}\f$, respectively.
  !!
  !! The amount of humidified litter, which is transferred from the litter to the soil
  !! carbon pool (\f$C_{\mathrm{D} \rightarrow \mathrm{H}}\f$) is modelled as a fraction
  !! of litter respiration (\f$R_{h, D}\f$) as
  !!
  !! \f[ \label{cdtoh} C_{\mathrm{D} \rightarrow \mathrm{H}} = \chi\, R_{h, D} \hspace{10pt}[Eqn 10] \f]
  !!
  !! where \f$\chi\f$ (see also ctem_params.f90) is the PFT-dependent humification factor
  !! and varies between 0.4 and 0.5. For crops, \f$\chi\f$ is set to 0.1 to account for
  !! reduced transfer of humidified litter to the soil carbon pool which leads to loss in
  !! soil carbon when natural vegetation is converted to croplands. Over the bare ground
  !! fraction \f$\chi\f$ is set to 0.45.
  !!
  !! With heterotrophic respiration known, net ecosystem productivity (\f$NEP\f$) is
  !! calculated as
  !! \f[ NEP = G_{canopy} - R_\mathrm{m} - R_\mathrm{g} - R_\mathrm{h}. \hspace{10pt}[Eqn 11] \f]
  !!

end module heterotrophicRespirationMod
