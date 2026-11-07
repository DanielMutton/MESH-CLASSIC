!> \file
!> Peatland specific parameterizations (photosynthesis and heterotrophic respiration)
module peatlandsMod

  ! J. Melton. Sep 26, 2016

  implicit none

  ! Subroutines contained in this module:
  public  :: mossPht
  public  :: hetresPeat
  public  :: peatDayEnd
  public  :: peatDepth

contains

  ! ------------------------------------------------------------------

  !> \ingroup peatlandsmod_mossPht
  !! @{
  !> Moss photosynthesis subroutine (equations are in module-level description)
  !> @author Yuanqiao Wu
  subroutine mossPht (il1, il2, iday, qswnv, thliq, co2conc, tsurfk, zsnow, &
                      pres, Cmossmas, dmoss, anmoss, rmlmoss, cevapmoss, &
                      ievapmoss, ipeatland, daylength, pdd)

    ! History

    ! J. Melton Sep 26 2016
    !   - Bring into rest of model and model formatting, convert to doxygen compatible code,
    !     trim out extraneous.
    ! ----------

    use classicParams,  only : rmlmoss25, tau25m, ektau, gasc, kc25, ko25, ec, ej, eo, evc, sj, &
                               hj, alpha_moss, thpmoss, thmmoss, ilg, ignd, TFREZ, RHOW

    implicit none

    ! arguments:

    integer, intent(in) ::  il1                     !<
    integer, intent(in) ::  il2                     !<
    integer, intent(in) ::  iday                    !<
    integer, intent(in) ::  ipeatland(ilg)          !<
    real, dimension(ilg), intent(in) :: qswnv       !< visible short wave radiation = qswnv in energBalNoVegSolve
    !! and qswnvg in energBalVegSolve (W/m2)
    real, dimension(ilg,ignd), intent(in) :: thliq  !<
    real, dimension(ilg), intent(in) :: zsnow       !<
    real, dimension(ilg), intent(in) :: co2conc     !<
    real, dimension(ilg), intent(in) :: daylength   !<
    real, dimension(ilg), intent(in) :: tsurfk      !< grid average ground surface temprature in K
    real, dimension(ilg), intent(in) :: pres        !<
    real, dimension(ilg), intent(in) :: Cmossmas    !< unit kg moss C updated in ctem
    real, dimension(ilg), intent(in) :: dmoss       !< unit m, depth of living moss. assume = 2 cm
    !! can be related to mmoss as a variable

    real, dimension(ilg), intent(inout) :: pdd(ilg) !<

    integer, dimension(ilg), intent(out) ::  ievapmoss   !< Value is 0 is no evaporation from moss, 1 otherwise
    real, dimension(ilg), intent(out) :: anmoss     !< net photosynthesis (umol CO2/m2/s)
    real, dimension(ilg), intent(out) :: rmlmoss    !< moss autotrophic respiration (umol CO2/m2/s)
    real, dimension(ilg), intent(out) :: cevapmoss  !< evaporation coefficent for moss surface

    ! Local variables

    integer :: i               !>
    integer :: j               !>
    real :: mmoss(ilg)      !< unit kg moss mass in kg, update later

    integer :: pheno(ilg)   !< phenology flag of mosses, 1 = photosynthesis,
    !! 0 = does not photosynthesis
    real :: parm(ilg)        !< par at the ground (moss layer) surface umol/m2/s
    real :: tsurf(ilg)       !< grid average ground surface temperature in C

    real :: wmoss(ilg)       !< water content extraporated from the surface
    !! humidity qg and thliq of the first soil layer
    !! unit kg water/ kg dw
    real :: wmosmin(ilg)     !< residual water content kg water /kg moss
    real :: wmosmax(ilg)     !< maximum water content kg water /kg moss
    real :: fwmoss(ilg)      !< relative water content of mosses in g fw /g dw
    real :: dsmoss(ilg)      !< degree of moss saturation = relative water
    !! content/maximum relative water content

    real :: g_moss(ilg)      !< moss conductance umol CO2/m2/s (based on
    !! Williams and Flanagan, 1998 for Sphagnum)
    real :: mwce (ilg)       !< moisture function of dark respiration of moss
    real :: tmoss(ilg)       !< moss temperature extraporated from the tbar 1
    !! and grid averaged ground surface temperature
    !! tsurf
    real :: tmossk(ilg)      !< moss temperature in K
    real :: q10rmlmos(ilg)   !< temperature function of the moss dark respiration
    real :: gamma(ilg)       !< compensation point for gross photosynthesis (Pa)
    real :: o2(ilg)          !< partial presure of oxygen (Pa)
    real :: co2a(ilg)        !< partical pressure of co2 (pa) same as in PHTSYN
    real :: tau(ilg)         !< arrhenius funciton of temperature
    real :: kc(ilg)          !< kinetic coeffficient of CO2 for photosynthesis
    real :: ko(ilg)          !< kinetic coeffficient of O2 for photosynthesis
    real :: bc(ilg)          !< coefficient used for wc
    real :: vcmax25(ilg)     !< seasonal varied maximum carboylation at 25
    !! sphagnum (fig. 6, Williams and Flanagan, 1998)
    real :: vcmax(ilg)       !< max carboxylation rate (umol/m2/s)
    real :: jmax25(ilg)      !< maximum electorn transport rate at 25 degrees (umol/m2/s)
    real :: jmax(ilg)        !< maximum electorn transport rate (umol/m2/s)
    real :: wj(ilg)          !< net co2 assimilation rate limited by
    !! electron transport (umol/m2/s)=jE in PHTSYN
    real :: wc(ilg)          !< net co2 assimilation rate limited by
    real :: ws(ilg)          !< net co2 assimilation rate limited by
    !! sucrose availability (umol/m2/s)=JE IN PHTSYN
    real :: photon(ilg)      !< electron transport rate (umol/m2/s)
    real :: term1(ilg)       !< temporary terms for photosynthesis calculations
    real :: term2(ilg)       !< temporary terms for photosynthesis calculations
    real :: term3(ilg)       !< temporary terms for photosynthesis calculations
    real :: psna(ilg)        !< coefficients for quadratric solution of net photosynthesis
    real :: psnb(ilg)        !< coefficients for quadratric solution of net photosynthesis
    real :: psne(ilg)        !< coefficients for quadratric solution of net photosynthesis
    real :: mI(ilg)          !< coefficients of the solutions for net psn
    real :: mII(ilg)         !< coefficients of the solutions for net psn

    ! Local parameters:
    real, parameter :: tref = 298.16    !< unit K

    ! ...........................
    ! Begin calculations:

    !     PHOTOSYNTHESIS COUPLING OR CURVATURE COEFFICIENTS
    ! real, parameter :: BETA1 = 0.950
    ! real, parameter :: BETA2 = 0.990

    !>    Find the light level (parm) at the ground surface for moss photosynthesis,
    !!    and a scaling factor degree of saturation of the moss layer phenology
    !!    parm is in umol/m2/s and converted from qswnv in W/m2

    do i = il1,il2 ! loop 100
      parm(i) = qswnv(i) * 4.6
      o2(i)  = 20.9/100.0 * pres(i)
      co2a(i) = co2conc(i)/1000000.0 * pres(i)
      wj(i)     = 0.0
      ws(i)     = 0.0
      wc(i)     = 0.0
      anmoss(i) = 0.0
      rmlmoss(i) = 0.0
      tsurf(i) = tsurfk(i) - tfrez
      tmossk(i) = tsurfk(i)
      tmoss(i) = tsurf(i)
    end do ! loop 100

    do i =  il1,il2 ! loop 200

      if (ipeatland(i) > 0) then ! only do for peatland regions.

        !>
        !!    phenology  water factor on mosses, grow when temperature > -4.0
        !!    and snowpack < 0.15 m

        if (zsnow(i) > 0.05 .or. tsurf (i) < 0.5) then
          pheno(i) = 0
        else
          pheno(i) = 1
        end if

        !>    ** water content used for the living moss (depth dmoss)
        !!    dmoss is an input and site specific. Preferably make dmoss a function
        !!    of Cmoss and ipeatland (different species in fens and bogs)
        !!    observed range of wmoss: 5 to 40 in Robrek (2007, 2009), 5 to 25
        !!    (Flanagen and Williams 1998)
        !!    dmoss is between 2.5 to 5cm based on the species (Lamberty et al. 2006)

        mmoss(i) = Cmossmas(i)/0.46
        wmoss(i) = thliq(i,1) * rhow/(mmoss(i)/dmoss(i))
        wmosmax(i) = min(45.0,thpmoss * dmoss(i) * rhow/mmoss(i))  ! FLAG dmoss was (1) not (i)!!
        wmosmin(i) = max(5.0,thmmoss * dmoss(i) * rhow/mmoss(i))  ! FLAG dmoss was (1) not (i)!!
        wmoss(i) = min(wmosmax(i),max(wmosmin(i),wmoss(i)))
        fwmoss(i) = wmoss(i) + 1.     ! g fresh weight /g dry weight

        !>    ** moss conductance g_moss umol CO2/m2/s
        !!   (Williams and Flanagan, 1998 for Sphagnum). follow MWM, fwmoss is
        !!    the mosswat_fd in MWM. Empirical equation is only valid up to
        !!   fwmoss=13, above 13 apply a linear extension to the equation.

        if (fwmoss(i) <= 13.0) then
          g_moss(i) = - 0.195 + 0.134 * fwmoss(i) - 0.0256 * (fwmoss(i)) &
                      ** 2 + 0.00228 * (fwmoss(i)) ** 3 - 0.0000984 * &
                      (fwmoss(i)) ** 4 + 0.00000168 * (fwmoss(i)) ** 5
        else
          g_moss (i) = - 0.000447 * fwmoss(i) + 0.0489
        end if

        g_moss (i) = g_moss(i) * 1000000.0
        g_moss(i) = max(0.0,g_moss(i))

        !>    ** moss surface evaporation coefficient
        !!    controled by the degree of saturation in moss, pass to energBalVegSolve and energBalNoVegSolve
        !!    apply a similar equation of soil surface cevap in energyBudgetPrep
        !!   CEVAP = 0.25*[1 – cos(THLIQ*pi/THFC)]^2

        dsmoss(i) = (wmoss(i) - wmosmin(i))/(wmosmax(i) - wmosmin(i))
        if (dsmoss(i) <   0.001) then
          ievapmoss(i) = 0
          cevapmoss(i) = 0.
        else if (dsmoss(i) >= 1.0) then
          ievapmoss(i) = 1
          cevapmoss(i) = 1.0
        else
          ievapmoss(i) = 1
          cevapmoss(i) = 0.25 * (1.0 - cos(3.14159 * dsmoss(i))) ** 2
        end if


        !>    ** moss water content effect on dark respiration
        !!   in MWM and PDM an optimal wmoss is at 5.8 gw/gdw(fig. 2e, Frolking et al., 1996)
        !!   Recent studies show weak but significant increases of sphagnum dark respiration
        !!   with moss water content above 5.8 gw/gdw (Adkinson and Humphreys, 2011 and ref.)
        !!   this change has improved the ER simulation greatly
        if (wmoss(i) < 0.4) then
          mwce (i) = 0.0
        else if (wmoss(i) < 5.8 .and. wmoss(i) > 0.4) then
          mwce(i) = 0.35 * wmoss(i) ** (2.0/3.0) - 0.14
        else
          !              mwce(i) =-0.04*wmoss(i)+1.232 ! in MWM and PDM
          mwce(i) = 0.01 * wmoss(i) + 0.942
        end if

        !>    ** moss dark respiration
        !!    observed range of rmlmoss 0.60 to 1.60 umol/m2/s (e.g. Adkinson 2006)

        q10rmlmos(i) = (3.22 - (0.046 * tmoss(i))) ** ((tmoss(i) - 25.0)/10.0)
        rmlmoss(i) = rmlmoss25 * mwce(i) * q10rmlmos(i)

        !>   ** moss photosynthesis
        !!   calculate bc (coefficient used for Wc, limited by Rubisco)

        tau(i) = tau25m * exp((tmossk(i) - tref) * ektau/(tref * gasc * tmossk(i)))
        gamma(i) = 0.5 * o2(i)/ tau(i)
        kc(i) = kc25 * exp((tmossk(i) - tref) * ec/(tref * gasc * tmossk(i)))
        ko(i) = ko25 * exp((tmossk(i) - tref) * eo/(tref * gasc * tmossk(i)))
        bc(i)  = kc(i) * (1.0 + (o2(i)/ ko(i)))

        !    seasonal change of Vcmax, sphagnum (fig. 6, Williams and Flanagan, 1998)
        !    from May 1st to september 1st Vmax is maximum
        !         if (iday < 121)        then
        !              vcmax25(i) = 6.0
        !         else if (iday > 245)   then
        !              vcmax25(i) = 7.0
        !         else
        !              vcmax25(i) = 13.5
        !         end if

        !!     use a function that is also valid for the southern hemisphere
        if (daylength(i) > 14.0 .and. pdd(i) > 200. .and. pdd(i) < 2000. ) then
          vcmax25(i) = 14.0
        else
          vcmax25(i) = 6.5
        end if

        vcmax(i) = vcmax25(i) * exp((tmossk(i) - tref) * evc/(tref * gasc * tmossk(i)))

        !>     calculate ws (phototysnthesis rate limited by transport capacity)
        !!    = js in photosynCanopyConduct

        ! if (coszs(i)>0.0)     THEN
        if (qswnv(i) > 0.) then ! FLAG JM version.
          ws(i) = 0.5 * vcmax(i)
        end if


        !>    calculate the maximum electron transport rate Jmax (umol/m2/s)
        !!   1.67 = vcmax25m/jmax25m ratio

        jmax25(i) = 1.67 * vcmax25(i)
        term1(i) = exp(((tmossk(i)/tref) - 1.) * ej/(gasc * tmossk(i)))
        term2(i) = 1. + exp(((tref * sj) - hj)/(tref * gasc))  ! FLAG: constant, can be pre-calculated !
        term3(i) = 1. + exp(((sj * tmossk(i)) - hj)/(gasc * tmossk(i)))
        jmax(i) = jmax25(i) * term1(i) * term2(i) * term3(i)
        if (jmax(i) > 0.0) then
          photon(i) = alpha_moss * parm(i)/sqrt(1.0 + (alpha_moss ** 2 * parm(i) ** 2/(jmax(i) ** 2)))
        else
          photon(i) = 0.0     ! electron trasport rate in mosses
        end if

        !>    calculate Wj, Wc (Farquhar and Caemmerer 1982)
        !!   wj = light limited, = je in PHTSYN

        wj(i) = photon(i) * (co2a(i) - gamma(i))/(4. * co2a(i) + (8. * gamma(i)))

        !>    carboxylase(rubisco) limitation = jc in PHTSYN
        wc(i) = vcmax(i) * (co2a(i) - gamma(i))/(co2a(i) + bc(i))

        !>    Choose the minimum of Wj and Wj both having the form:
        !!   W = (a Ci - ad) / (e Ci + b)
        !!   Then set a, b, d and e for the quadratic solution for net photosynthesis.

        if (wj(i) < wc(i)) then
          psnb(i) = 8. * gamma(i)
          psna(i) = photon(i)
          psne(i) = 4.0
        else if (wc(i) < wj(i)) then
          psnb(i) = bc(i)
          psna(i) = vcmax(i)
          psne(i) = 1.0
        end if

        !>    Calculate net and gross photosynthesis by solve the quadratic equation
        !!   first root of solution is net photosynthesis An= min(Wj,Wc) - Rd
        !!   gross photosynthesis GPP = min(Wc,Wj) = An + Rd

        if (psna(i) > 0.0) then
          mI(i) = rmlmoss(i) - (psnb(i) * g_moss(i)/pres(i)/psne(i)) &
                  - (co2a(i) * g_moss(i)/pres(i)) - (psna(i)/psne(i))

          mII(i) = (psna(i) * co2a(i) * g_moss(i)/pres(i)/psne(i)) - &
                   (rmlmoss(i) * co2a(i) * g_moss(i)/pres(i)) - (rmlmoss(i) * &
                   psnb(i) * g_moss(i)/pres(i)/psne(i)) - (psna(i) * gamma(i) * &
                   g_moss(i)/pres(i)/psne(i))
        else
          mI(i) = 0.0
          mII(i) = 0.0
        end if
        anmoss(i) = ( - mI(i) - (mI(i) * mI(i) - 4 * mII(i)) ** 0.5)/2
        anmoss(i) = pheno(i) * min(anmoss(i),ws(i))

      end if ! ipeatland.

    end do ! loop 200

    return
  end subroutine mossPht
  !! @}

  ! ---------------------------------------------------------------------------------------------------

  !> \ingroup peatlandsmod_hetresPeat
  !! @{
  !> Grid average peat soil heterotrophic respiration subroutine (equations are in module-level description)
  !> @author Yuanqiao Wu
  subroutine  hetresPeat (il1, il2, ilg, ipeatland, & ! In
                          isand, litrmsmoss, peatdep, wtable, & ! In
                          tbar, thliq, thice, thpor, & ! In
                          bi, zbotw, delzw, psisat, &! In
                          useTracer, tracerMossLitrMass, & ! In
                          litresms, socresp, resoxic, resanoxic, & ! Out
                          litResMossTracer, soCResPeatTracer) ! Out

    !   History:

    ! J. Melton Sep 26 2016
    !   - Bring into rest of model and model formatting, convert to doxygen compatible code
    !   ----------------------------------

    use classicParams, only :icc, ignd, tanhq10, dctmin, dcbaset, bsrateltms, TFREZ

    implicit none

    integer, intent(in) :: il1, il2, ilg
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    integer, dimension(ilg,ignd), intent(in) :: isand   !<
    integer, dimension(ilg), intent(in) :: ipeatland    !< peatland flag, 0 = not peatland, 1 = bog, 2 = fen
    real, dimension(ilg), intent(in) :: peatdep         !<
    real, dimension(ilg), intent(in) :: litrmsmoss      !<
    real, dimension(ilg,ignd), intent(in) :: thliq      !<
    real, dimension(ilg,ignd), intent(in) :: thice      !<
    real, dimension(ilg,ignd), intent(in) :: thpor      !<
    real, dimension(ilg,ignd), intent(in) :: bi         !<
    real, dimension(ilg,ignd), intent(in) :: zbotw      !<
    real, dimension(ilg,ignd), intent(in) :: delzw      !<
    real, dimension(ilg,ignd), intent(in) :: psisat     !< saturated matrix potential in soil (m)
    real, dimension(ilg), intent(in) :: wtable          !< water table (m)
    real, dimension(ilg,ignd), intent(in) :: tbar       !< soil temperature (K)
    real, intent(in) :: tracerMossLitrMass(:)   !< Tracer mass in moss litter, \f$kg C/m^2\f$

    real, dimension(ilg), intent(out) :: litresms        !< moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: socresp         !< soil C respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: resoxic         !< respiration rate of the oxic compartment (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: resanoxic       !< respiration rate of the oxic compartment (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: litResMossTracer(ilg)    !< Tracer moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: soCResPeatTracer(ilg)    !< Tracer heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    !     internal variables
    integer :: i,j
    real :: Cso (ilg)            !< carbon mass in the oxic compartment (\f$kg C m^{-2}\f$)
    real :: Csa (ilg)            !< carbon mass in the anxic compartment (\f$kg C m^{-2}\f$)
    ! real :: tracerCso (ilg)      !< Tracer carbon mass in the oxic compartment (\f$kg C m^{-2}\f$)
    ! real :: tracerCsa (ilg)      !< Tracer carbon mass in the anxic compartment (\f$kg C m^{-2}\f$)
    real :: fto(ilg)             !< temperature factor of the oxic soil respiration
    real :: fta(ilg)             !< temperature factor of the anoxic soil respiration
    real :: tsoilo(ilg)          !< average temperature in the oxic compartment(C)
    real :: tsoila(ilg)          !< average temperature of the anoxic compartment (C)
    real :: ewtable(ilg)         !< effective water table depth (m)
    real :: ratescpo(ilg)        !< oxic respiration rate constant(\f$\mu mol CO_2 [kg C]^{-2} s^{-1}\f$)
    real :: ratescpa(ilg)        !< anoxic respiration rate constant(\f$\mu mol CO_2 [kg C]^{-2} s^{-1}\f$)
    real :: psi(ilg,ignd)        !< matrix potential of soil layers (Pa)
    real :: q10funcms(ilg)       !< q10 fuction for moss litter respiration
    real :: ltrmosclms(ilg)      !< moisture scale for litter respiration
    real :: litrtempms(ilg)      !< temperature of the moss litter
    real :: litpsims(ilg)        !< matrix potential of litter
    real :: litrq10ms(ilg)       !< q10 coefficient as a functon of T (as in CTEM)
    real :: soilq10o(ilg)        !< q10 coefficient of oxic soil
    real :: soilq10a(ilg)        !< q10 coefficient of anoxic soil
    integer :: lewtable(ilg)     !< layer index of the water table layer

    !    ------------------------------------------------------------------
    !
    !    initialization
    litresms(:) = 0.0
    socresp(:) = 0.0
    resoxic(:) = 0.0
    resanoxic(:) = 0.0
    tsoilo(:) = 0.0
    tsoila(:) = 0.0
    litResMossTracer = 0.
    soCResPeatTracer = 0.

    !> Calculate soil respiration from peat

    !> Find the effective water table depth and the layer index to divide
    !! the peat soil into two compartments
    do i = il1,il2

      if (ipeatland(i) > 0) then ! only operate on peat tiles

        ewtable(i) = wtable(i)
        if (ewtable(i) <= 0.0) then ! flooded
          lewtable(i) = 0
        else if (ewtable(i) <= 0.1) then ! below ground surface but in the moss layer
          lewtable(i) = 1
        else ! deeper in the soil column.
          do j = 1,ignd
            if (ewtable(i) > zbotw(i,j)) then
              lewtable(i) = j + 1
            end if
          end do
        end if

        !> Find the temperature in litter, oxic soil and anoxic soil in kelvin
        !! lewtable is the layer index of the water table layer, lewtable = 0
        !! indicates WTD is above the ground surface.
        !! Set the oxic layer temperature to dctmin (minimum soil respiration
        !! temperature) when the entire soil is in the anoxic zone.

        if (lewtable(i) == 0) then  ! WT is at or above the surface

          do j = 1,ignd ! so find the temp of the total soil column.
            ! FLAG JM - this might not be appropriate for runs with the newer deeper soils and 20 layers ! Sep 30 2016.
            tsoila(i) = tsoila(i) + tbar(i,j) * delzw(i,j)
          end do
          tsoila(i) = tsoila(i) / zbotw(i,ignd)
          tsoilo(i) = dctmin

        else if (lewtable(i) == 1) then     ! WT is at the first layer

          tsoilo(i) = tbar(i,1)
          do j = lewtable(i) + 1,ignd
            tsoila(i) = tsoila(i) + tbar(i,j) * delzw(i,j)
          end do
          tsoila(i) = (tsoila(i) + tbar(i,1) * (zbotw(i,lewtable(i)) - ewtable(i))) &
                      / (zbotw(i,ignd) - ewtable(i))

        else                                ! WT is below layer 1

          do j = 1,lewtable(i) - 1
            tsoilo(i) = tsoilo(i) + tbar(i,j) * delzw(i,j)
          end do
          tsoilo(i) = (tsoilo(i) + tbar(i,lewtable(i)) * (ewtable(i) - zbotw(i,lewtable(i) - 1))) &
                      / ewtable(i)
          do j = ignd,lewtable(i) + 1, - 1
            tsoila(i) = tsoila(i) + tbar(i,j) * delzw(i,j)
          end do
          tsoila(i) = (tsoila(i) + tbar(i,lewtable(i)) &
                      * (zbotw(i,lewtable(i)) - ewtable(i)) ) &
                      / (zbotw(i,ignd) - ewtable(i))

        end if

        !> Calculate the temperature multiplier (ftsocres) for oxic and anoxic
        !! soil compartments

        soilq10o(i) = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - (tsoilo(i) - tfrez))))
        soilq10a(i) = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - (tsoila(i) - tfrez))))

        fto(i) = soilq10o(i) ** (0.1 * (tsoilo(i) - tfrez - 15.0))
        fta(i) = soilq10a(i) ** (0.1 * (tsoila(i) - tfrez - 15.0))

        !> Find the heterotrophic respiration rate constant in the oxic and
        !! anoxic (unit in yr-1), based on Fig.2b in Frolking 2001

        if (ipeatland(i) == 1) then ! bogs
          if (ewtable(i) < 0.0) then ! flooded

            ratescpo(i) = 0.0
            ratescpa(i) = - 0.183 * exp( - 18.0 * peatdep(i)) + 0.03 * peatdep(i) + 0.0134

          else if (ewtable(i) < 0.30 .and. ewtable(i) >= 0.0) then ! within the first 30 cm of surface

            ratescpo(i) = 0.009 * (1 - exp( - 20. * ewtable(i))) + 0.015 * ewtable(i)
            ratescpa(i) = 0.009 * exp( - 20. * ewtable(i)) - 0.183 * exp( - 18. * peatdep(i)) - 0.015 * ewtable(i) + 0.0044

          else if (ewtable(i) >= 0.30) then ! deeper in the soil column

            ratescpo(i) = 0.0134 - 0.183 * exp( - 18. * ewtable(i)) + 0.003 * ewtable(i)
            ratescpa(i) = - 0.183 * exp( - 18. * peatdep(i)) + 0.003 * (peatdep(i) - wtable(i)) &
                          + 0.183 * exp( - 18. * ewtable(i))
            ! ratescpa(i)=-0.183*exp(-18*peatdep(i))+0.003*(peatdep(i)-wtable(i))+0.183*exp(-18*ewtable(i))-0.004504   ! for continuity
          end if
        else if (ipeatland(i) == 2) then ! fens
          if (ewtable(i) < 0.0) then ! flooded

            ratescpo(i) = 0.0
            ratescpa(i) = 0.01512 - 1.12 * exp( - 25. * peatdep(i))

          else if (ewtable(i) < 0.30 .and. ewtable(i) >= 0.0) then ! within the first 30 cm of surface

            ratescpo(i) = - 0.01 * exp( - 40. * ewtable(i)) + 0.015 * ewtable(i) + 0.01
            ratescpa(i) = abs( - 0.01 * exp( - 40. * ewtable(i)) - 1.12 * exp( &
                          - 25. * peatdep(i)) + 0.015 * ewtable(i) + 0.005119)

          else if (ewtable(i) >= 0.30) then ! deeper in the soil column

            ratescpo(i) = 0.01512 - 1.12 * exp( - 25 * ewtable(i))
            ratescpa(i) = - 1.12 * (exp( - 25. * peatdep(i)) - exp( - 25. * ewtable(i)))
          end if
        end if

        !>  Convert respiration rates from kg c/kg c.year to u-mol co2/kgC/s
        ratescpo(i) = 2.64 * ratescpo(i)
        ratescpa(i) = 2.64 * ratescpa(i)

        !> Find the carbon storage in oxic and anoxic compartments (Cso. Csa)
        !! The water table depth delineates the oxic and anoxic compartments.
        !! functions (R**2 = 0.9999) determines the carbon content of each
        !! compartment from a peat bulk density profile based on unpulished
        !! data from P.J.H. Richard (described in fig. 1, Frokling et al.(2001)
        !! conversion of peat into carbon with 48.7% (Mer Bleue unpublished data,
        !! Moore)

        Cso(i) = peatStorage(ewtable(i))
        Csa(i) = peatStorage(peatdep(i)) - Cso(i)

        !> Find the soil respiration rate in Cso and Csa umol/m2/s.
        !! Moisture multiplier (0.025) indicates rate reduction in decomposition due
        !! to anoxia (Frolking et al. 2001), only applied to anoxic layer

        resoxic(i)   = ratescpo(i) * Cso(i) * fto(i)
        resanoxic(i) = ratescpa(i) * Csa(i) * fta(i) * 0.025
        socresp(i)   = resoxic(i) + resanoxic(i)


        !> Calcualte litter respiration of moss

        !> First find the matrix potential of the soil layers
        do j = 1,ignd
          if (isand(i,j) == - 3 .or. isand(i,j) == - 4) then  ! ice or rock
            psi(i,j) = 10000.0 ! a large number so that ltrmoscl = 0.2
          else ! soils
            if (thliq(i,j) + thice(i,j) + 0.01 < thpor(i,j) .and. tbar(i,j) < 273.16) then
              psi(i,j) = 0.001
            else if (thice(i,j) > thpor(i,j)) then
              psi(i,j) = 0.001   ! set to saturation
            else
              psi(i,j) = psisat(i,j) * (thliq(i,j) / (thpor(i,j) - thice(i,j))) ** ( - bi(i,j))
            end if
          end if
        end do ! loop 60

        !> Litter in peatlands can be saturated so we limit the rate by high
        !! moisuture level similar to soil in CTEM, but less effectively (the
        !! min moisture factor is at 0.5 for moss litter but at 0.2 for soil).

        litpsims(i) = psi(i,1)

        ! limit of ltrmoscalms at saturation
        if (litpsims(i) > 10000.0) then
          ltrmosclms(i) = 0.2
        else if (litpsims(i) <= 10000.0 .and. litpsims(i) >  6.0) then
          ltrmosclms(i) = 1.0 - 0.8 * ((log10(litpsims(i)) - log10(6.0)) / (log10(10000.0) - log10(6.0))) ** 1
        else if (litpsims(i) <= 6.0 .and. litpsims(i) > 4.0) then
          ltrmosclms(i) = 1.0
        else if (litpsims(i) <= 4.0 .and. litpsims(i) > psisat(i,1)) then
          ltrmosclms(i) = 1.0 - 0.99 * ((log10(4.0) - log10(litpsims(i))) / (log10(4.0) - log10(psisat(i,1))))
        else if (litpsims(i) <= psisat(i,1)) then
          ltrmosclms(i) = 0.01
        end if

        ltrmosclms(i) = max(0.0,min(ltrmosclms(i),1.0))

        !> Find the temperature factor for moss litter respiration
        litrtempms(i) = tbar(i,1) - tfrez
        litrq10ms(i) = tanhq10(1) + tanhq10(2) * (tanh(tanhq10(3) * (tanhq10(4) - litrtempms(i))))
        q10funcms(i) = litrq10ms(i) ** (0.1 * (litrtempms(i) - 15.0))

        !> Calculate the litter respiration rate in mosses and converts it
        !! from kg c/kg c.year to u-mol co2/kg c.s using 2.64
        litresms(i) = ltrmosclms(i) * litrmsmoss(i) * bsrateltms * 2.64 * q10funcms(i)

      end if ! ipeatland
    end do ! i loop
    return

  end subroutine hetresPeat

  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup peatlandsmod_peatDayEnd
  !! @{
  !> At the end of the day update the degree days for moss photosynthesis and the peat bottom layer depth
  !> @author Yuanqiao Wu
  subroutine peatDayEnd (nml)

    use ctemStateVars, only : ctem_tile, vgat
    use classStateVars, only : class_gat
    use classicParams, only : ignd, TFREZ

    implicit none

    integer, intent(in) :: nml

    integer :: i

    real, pointer, dimension(:)    :: taaccgat_t     !< Daily mean air temperature [K]
    real, pointer, dimension(:)    :: pddgat         !< peatland degree days above 0 deg C.
    integer, pointer, dimension(:) :: ipeatlandgat   !< peatland flag, 0 = not peatland, 1 = bog, 2 = fen
    real, pointer, dimension(:,:)  :: dlzwgat        !< Permeable thickness of soil layer [m]
    real, pointer, dimension(:)    :: peatdepgat     !< Depth of peat column [m]
    real, pointer, dimension(:)    :: sdepgat        !< Depth to bedrock in the soil profile [m]

    taaccgat_t       => ctem_tile%taaccgat_t
    pddgat           => vgat%pdd
    ipeatlandgat     => vgat%ipeatland
    dlzwgat          => class_gat%dlzwgat
    peatdepgat       => vgat%peatdep
    sdepgat          => vgat%sdepgat

    !> Calculate degree days for mossPht Vmax seasonality (only once per day)
    do   i = 1,nml
      if (taaccgat_t(i) > tfrez) then
        pddgat(i) = pddgat(i) + taaccgat_t(i) - tfrez
      end if

      !> Update peatland bottom layer depth
      if (ipeatlandgat(i) > 0) then
        dlzwgat(i,ignd) = peatdepgat(i) - 0.90
        sdepgat(i) = peatdepgat(i)
      end if
    end do

  end subroutine peatDayEnd

  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup peatlandsmod_peatDepth
  !! @{
  !> Calculate the peat depth based on equation 18 in Wu, Verseghy, Melton 2016 GMD.
  !! @author Y. Wu, J. Melton
  !!
  real function peatDepth (gavgscms)

    implicit none

    real, intent(in) :: gavgscms  !< Grid averaged soil c mass, \f$kg c/m^2\f$

    ! Calculate the peat depth based on equation 18 in Wu, Verseghy, Melton 2016 GMD.
    peatDepth = ( - 72067.0 + sqrt((72067.0 ** 2.0) - (4.0 * 4056.6 &
                * ( - gavgscms * 1000 / 0.487)))) / (2 * 4056.6)
  end function peatDepth
  !! @}

  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup peatlandsmod_peatStorage
  !! @{ Finds the carbon storage in the peat based on depth (or oxic and anoxic compartments)
  !! The water table depth delineates the oxic and anoxic compartments.
  !! functions (R**2 = 0.9999) determines the carbon content of each
  !! compartment from a peat bulk density profile based on unpulished
  !! data from P.J.H. Richard (described in fig. 1, Frokling et al.(2001)
  !! conversion of peat into carbon with 48.7% (Mer Bleue unpublished data,
  !! Moore)
  !>
  !! @author Y. Wu, J. Melton
  !!
  real function peatStorage (depth)

    implicit none

    real, intent(in) :: depth  !< Peat compartment depth (m). Either total column or oxic/anoxic.

    peatStorage = (4056.6 * depth ** 2 + 72067.0 * depth) * 0.487 / 1000.0

  end function peatStorage
  !! @}
  ! ---------------------------------------------------------------------------------------------------

  !> \namespace peatlandsmod
  !! Peatland specific processes 
  !! @author Y. Wu, J. Melton
  !!
  !! The peatland module is published in Geoscientific Model Development (Wu et al. 2016) \cite Wu2016-zt.
  !!
  !! To account for the eco-hydrological and biogeochemical interactions among
  !! vegetation, atmosphere and soil in peatlands, the following modifications
  !! were made to the coupled CLASS3.6--CTEM2.0 modelling framework:
  !!
  !! 1. The top soil layer was characterized as a moss layer with a higher heat
  !! and hydraulic capacity than a mineral soil layer. The moss layer buffers the
  !! exchange of energy and water at the soil surface and regulates the soil
  !! temperature and moisture (Turetsky et al. 2012) \cite Turetsky2012-qh.
  !!
  !! 2. Three peatland vascular PFTs (evergreen shrubs, deciduous shrubs and
  !! sedges) as well as mosses were added to the existing nine CTEM PFTs. These
  !! peatland-specific PFTs are adapted to cold climate and inundated soil with
  !! optimized plant structure (shoot/root ratio, rooting depth), growth strategy
  !! and metabolic acclimations to light, water and temperature.
  !!
  !! 3. We considered the soil inundation stress on microbial respiration in the
  !! litter C pool. The original CTEM assumed that litter respiration was not
  !! affected by oxygen deficit as a result of flooding, since litter was always
  !! assumed to have access to air. This assumption does not hold for peatlands
  !! where high water table positions occur routinely.
  !!
  !! 4. To provide the framework for future runs coupled to the global earth system model, we separated the soil C balance and heterotrophic respiration
  !! (HR) calculations for peatland and non-peatland fractions for each grid cell
  !! in the global model. Over the non-peatland fraction, we use the original CTEM
  !! approach that aggregates the HR from each PFT weighted by the fractional
  !! cover. Over the peatland fraction the soil C pool and decomposition are
  !! controlled by the water table position, following the two-compartment
  !! approach used in the MWM (St. Hilaire et al. 2010) \cite St-Hilaire2010-5e9.
  !!
  !! The standard configuration of soil layers in CLASS consists of three layers with
  !! thickness of 0.10, 0.25, and 3.75, m. Organic soil in CLASS was parameterized
  !! by Letts et al. (2000) \cite Letts2000-pg as fibric, hemic and sapric peat in the three soil
  !! layers respectively, representing fresh, moderately decomposed and highly
  !! decomposed organic matter. Tests of CLASS on peatlands revealed improved
  !! performance in the energy simulations for fens and bogs with this organic
  !! soil parameterization. However, the model overestimated energy and water
  !! fluxes at bog surfaces during dry periods due to the neglect of the moss
  !! cover (Comer et al., 2000) \cite Comer2000-mz.
  !!
  !! To take into account the interaction amongst the moss and the soil layers and
  !! the overlying atmosphere for energy and water transfer, we added a new soil
  !! layer 0.10 m thick above the fibric organic soil to represent living and dead
  !! peatland bryophytes, such as Sphagnum mosses and true mosses
  !! (Bryopsida). The physical characteristics of mosses differ from those of
  !! either the shoots or the roots of vascular plants (Rice et al., 2008). In
  !! particular, mosses can hold more than 30 g of water per gram of biomass
  !! (Robroek et al., 2009). More than 90 \% of the moss leaf volume is
  !! occupied by the water-holding hyaline cells (Rice et al., 2008), which retain
  !! water even when the water table depth declines to 1--10 m below the surface
  !! (Hayward and Clymo, 1982).
  !!
  !! The parameter values of the moss layer for water and energy properties were
  !! derived from a number of recent experiments measuring the hydraulic
  !! properties of mosses (Price et al., 2008 \cite Price2008-fr; Price and Whittington, 2010
  !! McCarter and Price, 2012). Living mosses range from 2--3 to over
  !! 5 cm in height (Rice et al., 2008) and have lower values of dry bulk density
  !! and field capacity than fibric peat (Price et al., 2008) \cite Price2008-fr. Compared to fibric
  !! peat, the saturated hydraulic conductivity of living moss is higher by orders
  !! of magnitude (Price et al., 2008) \cite Price2008-fr and the thermal conductivity is more
  !! affected by the water content (O'Donnell et al., 2009). To fully account for
  !! the effect of mosses, we set the depth of the living moss (\f$z_{m}\f$)
  !! within the top soil (i.e. moss) layer to 3 cm for fens and 4 cm for bogs,
  !! and interpolated its water content \f$w_m\f$ (kg water /(kg dry mass)) from the water content of the overall layer
  !! \f$\theta_{l, 1}\f$ (m\f$^3\f$ water /(m soil)\f$^3\f$) and the depth of the
  !! living moss:
  !!
  !! \f$w_m =\frac{z_m\theta_{l, 1}\rho_w}{B_m} \f$
  !!
  !! where the dry moss biomass (\f$B_m) \f$ is converted from moss C
  !! (C\f$_m \f$) using the standard conversion factor of 0.46 kg C per kg dry
  !! biomass, \f$\theta _{l, 1} \f$ (m\f$^3\f$\, m\f$^{-3}\f$) is the liquid water
  !! content of the top soil layer, and \f$\rho_w \f$ is the density of water
  !! (1000 kg m\f$^{-3} \f$). The maximum and minimum moss water contents were
  !! estimated from a number of observed moss water contents (e.g. Williams and Flanagan, 1998; Robroek et al., 2009). In CLASS, evaporation at the soil
  !! surface is controlled by a soil evaporation efficiency coefficient \f$\beta \f$
  !! (Verseghy, 2012). This parameter is calculated from the liquid water content
  !! and the field capacity of the first soil layer following Lee and
  !! Pielke (1992). For peatlands, \f$\beta \f$ was assumed to be regulated by the
  !! relative moisture of the living moss rather than the ratio of relative liquid
  !! water content of the first soil layer:
  !!
  !! \f$ \beta = 0.25 [ 1- \cos \left( \frac{w_m
  !! -w_{m, min}}{w_m-w_{m, max}} \right)]^{2} \f$
  !!
  !! where \f$w_m \f$, \f$w_{m, max} \f$, and \f$w_{m, min} \f$ are the water content
  !! and the maximum and minimum water contents of the living moss in kg water / (kg dry mass).
  !!
  !!
  !! Moss photosynthesis subroutine
  !!
  !! Mosses are an important contributor to the primary production and the C
  !! sequestration in peatlands, owing to the low decomposability of the moss
  !! tissue. Sphagnum in peatlands grows at
  !! 20--1600 g biomass m\f$^{-2} \f$/yr and accounts for about 50 \% of
  !! the total peat volume (Turetsky, 2003). We have modified CTEM to include a
  !! moss C pool and moss litter pool along with the related C fluxes, i.e.
  !! photosynthesis, autotrophic respiration, heterotrophic respiration, and
  !! humification. The net photosynthesis of moss (\f$G_m)\f$ is calculated
  !! from the gross photosynthesis (\f$G_{0, m})\f$ and dark respiration
  !! (\f$R_{d, m}\f$):
  !!
  !! \f$ G_m=G_{0, m}-R_{d, m} \f$.
  !!
  !! The moss photosynthesis and dark respiration are calculated using the
  !! Farquhar~(1989) biochemical approach following the MWM (St-Hilaire et al.,
  !! 2010) \cite St-Hilaire2010-5e9 and CTEM (Melton and Arora, 2016) \cite Melton2016-zx, with modifications for integration
  !! with CLASS--CTEM and moss phenology. The leaf-level gross photosynthesis rate
  !! \f$G_{0, m} \f$ (\f$\mu \f$ mol CO\f$_2 \f$ /m \f$^2 \f$/s) is obtained as
  !! the minimum of the transportation limited photosynthesis rates (\f$J_s)\f$ and
  !! the first root of the quadratic solution of the light-limited rate (\f$J_e)\f$
  !! and the Rubisco limited rate (\f$J_c)\f$. A logistic factor (\f$\varsigma\f$) is
  !! added with values 0 or 1 to introduce a seasonal control of moss
  !! photosynthesis. In the MWM, spring photosynthesis starts when the snow depth
  !! is below 0.05 m and the soil temperature at 5 cm depth goes above
  !! 0.5 C (Moore et al., 2006). Since in our case CLASS sets the
  !! minimum depth for melting, discontinuous snow to 0.10 m, this limits the
  !! spring photosynthesis to starting only once the snow is completely melted.
  !!
  !! \f$ G_{0, m} = \quad \varsigma \min \left(J_{s}, \frac{(J_{c}+J_{e})\pm
  !! \sqrt{{(J_c +J_ e)}^2 - 4(J_c +J_e)} }{2}\right) \f$
  !!
  !! The dark respiration in mosses (\f$R_{d, m})\f$ is calculated as a function
  !! of the base dark respiration rate (\f$R_{d, m, 0})\f$, which has a value of
  !! 1.1 (\f$\mu \f$ mol CO\f$_2 \f$ /m \f$^2 \f$/s) (Adkinson and Humphreys, 2011) scaled
  !! by the moss moisture (\f$f_{m, rd}\f$) and soil temperature functions
  !! (\f$f_{T, rd})\f$. The moss moisture function is based on the volumetric
  !! water content of the moss, \f$\theta_m\f$ (kg water /(kg dry
  !! mass)). The MWM models the relation between water content in mosses and dark
  !! respiration with optimal water content at 5.8 g water per g dry weight,
  !! following the approach in Frolking et al. (1996). We modified the relation
  !! for water content above the optimal water content, based on a recent
  !! discovery of a weak linear positive relation between the dark respiration
  !! rate and the water content above the optimal water content during the late
  !! summer and fall (Adkinson and Humphreys, 2011):
  !!
  !! \f$R_{d, m} = R_{d, m, 0} f_{m, rd} f_{T, {rd}}\f$
  !!
  !! \f$ f_{T, {rd}} = (3.22-(0.046 \cdot
  !! T_{moss})^{(T_{moss}-25/10)} \f$
  !!
  !! \f$ f_{m, rd} = 0\f$ for\f$ \theta_{m} <0.4 \f$
  !!
  !! \f$ f_{m, rd} = 0.35 \theta_{m}^{2/3}-0.14\f$ for\f$ 0.4\le \theta_m<5.8 \f$
  !!
  !! \f$ f_{m, rd} = 0.01 \theta_{m} +0.942\f$ for \f$5.8<\theta _{m} \f$
  !!
  !! Photosynthetic photon flux density (PPFD) is measured by the
  !! photosynthetically active radiation (PAR), which is defined as the solar
  !! radiation between 0.4 to 0.7 \f$\mu \f$ mol that can be used by plants via
  !! photosynthesis. In the coupled CLASS--CTEM system, the PAR received by the
  !! moss (PAR\f$_m\f$, unit \f$\mu \f$mol photons m\f$^-2\f$/s is
  !! converted from the visible short-wave radiation reaching the ground (\f$K_{\ast g}\f$,
  !! unit W/m\f$^2\f$) in CLASS by a factor of
  !! 4.6 \f$\mu\f$ mol /m\f$^2\f$/s per W/m\f$^2\f$ (McCree, 1972).
  !! \f$K_{\ast g}\f$ is a function of the incoming short-wave radiation
  !! (\f$K\downarrow \f$, unit: W/m\f$^2\f$), the surface albedo (\f$\alpha_g\f$),
  !! and the canopy transmissivity (\f$\tau _c\f$):
  !!
  !! \f$K_{\ast g} = K \downarrow \tau_{c}\left( 1-\alpha_{g}\right)\f$
  !!
  !! The energy uptake by the moss layer is thus a function of the total incoming
  !! short-wave radiation, the aggregated LAI of the PFTs present, the snow depth,
  !! the fractional vegetation cover, and the soil water content (Verseghy, 2012).
  !! In peatland C models that do not consider vegetation dynamics, the
  !! transmissivity of the vegetation canopy is usually assumed to be constant
  !! (e.g. St-Hilaire et al., 2010) \cite St-Hilaire2010-5e9. Compared with such models, CLASS enables a
  !! more detailed representation of light incident on the moss surface since it
  !! includes partitioning of direct/diffuse and visible/near-IR radiation,
  !! PFT-specific transmissivities, and time-varying LAI and fractional PFT
  !! coverages (Verseghy, 2012) \cite Verseghy2012-c0e.
  !!
  !!
  !! Peatland soil heterotrophic respiration
  !!
  !! Over the non-peatland fraction, HR is calculated
  !! as the sum of the respiration from litter and soil carbon pools as in the
  !! original version of CTEM (Arora, 2003). The soil C pool over the
  !! non-peatland areas is assumed to be exponentially distributed with depth
  !! (Arora, 2003). In peatlands a large amount of humic soil is generally
  !! located in the permanently saturated zone and the bulk density increases
  !! with soil depth (Loisel and Garneau, 2010). Thus, the assumption of
  !! exponentially decreasing distribution of C content with increasing soil
  !! depth is not valid in peatlands. We used a quadratic equation to calculate
  !! the distribution of soil C content over depth based on an empirically
  !! determined bulk density profile (Frolking et al., 2001).
  !!
  !! HR over the peatland fraction of a grid cell is modelled using a two-pool
  !! approach with a flexible boundary between the pools that depends on the
  !! depth of the water table:
  !!
  !!
  !! \f$ R_{o}=C_{SOM, o}k_{o}f_{T, {o}} \f$
  !!
  !! \f$ R_{a}=C_{SOM, a}k_{a}f_{T, {a}}f_{anoxic} \f$
  !!
  !! where o and a denote the oxic and anoxic portions of the soil C pool respectively. The respiration rate \f$R\f$ (unit:
  !! \f$\mu mol C m^{-2}\f$/s) is obtained from the respiration rate
  !! coefficient \f$k\f$ (\f$\mu\f$ mol C / kg C / s), the temperature
  !! functions \f$f_T\f$, the soil C mass \f$C_{SOM}\f$ (kg), and a scaling factor
  !! \f$f_{anoxic}\f$ after Frolking et al.~(2010, 2001), which represents the
  !! inhibition of microbial respiration under anoxic conditions. The value of
  !! this parameter is uncertain, varying in those two papers between 0.001, 0.025
  !! and 0.1. Based on calibration runs using two of the data sets described below
  !! (MB-Bog and AB-Fen), we adopted a value of 0.025. \f$Q_{10}\f$ is calculated using
  !! a hyperbolic tan function of the soil temperatures (\f$T_s)\f$ of the oxic
  !! and anoxic zones (Melton and Arora, 2016), which are in turn functions of
  !! water table depth. The \f$Q_{10}\f$ values of the anoxic and the oxic
  !! zones of the soil are indicated as \f$Q_{10, a}\f$ and \f$Q_{10, o}\f$.
  !! The values of \f$k\f$, \f$f_T\f$, and \f$C_{SOM}\f$ are updated along with the
  !! water table depth (\f$z_{wt}\f$, unit: m, positive downward) and the peat
  !! depth (\f$z_p\f$, unit: m) at each CTEM time step. The equations for \f$k\f$
  !! and \f$C_{SOM}\f$ are derived from Fig. 2 in Frolking et al. (2001), and
  !! parameterized differently for fens and bogs :
  !!
  !!  \f$ f_{T, o} =Q_{10, o}^{\left(\int\limits_{0}^{z_{wt}} T_{j} -15\right)/10} \f$
  !!
  !!  \f$ f_{T, {a}} =Q_{10, {a}}^{\left(\int\limits_{z_{wt}}^{z_\mathrm{p}} T_{j} -15\right)/10} \f$
  !!
  !! \f$ Q_{10}=1.44+0.56 {tanh}[0.075\left( 46.0-T_{s} \right)] \f$
  !!
  !! \f$ T_{s, o}=\int\limits_{0}^{z_{wt}} T_{j} \, /(z_{wt}) \f$
  !!
  !! \f$ T_{s, a}=\int\limits_{z_{wt}}^{z_{p}} T_{j} /(z_{p}-z_{wt}) \f$
  !!
  !! \f$ k_{o}= 0 \f$ for \f$ z_{wt}<0 \f$
  !!
  !!  \f$ k_{o}= k_{1}\left( 1-e^{k_{2}z_{wt}} \right)+k_{3}z_{wt} \f$ for \f$0.3 > z_t \ge 0\f$
  !!
  !!  \f$ k_{4}e^{k_{5}z_{wt}}+k_{6}z_{wt}+k_{7} \f$ for \f$ z_{wt}\ge 0.3 \f$
  !!
  !! \f$ k_{a}= k_{4}e^{k_{5}z_{p}}+{10k}_{6}z_{p}+k_{7} \f$ for \f$ z_{wt}<0 \f$
  !!
  !! \f$ k_{a} = \f$ | \f$ k_1 e^{ k_2 z_{wt} } - k_4 e^{k_5 z_p} -k_3 z_{wt}+k_8 \f$ for \f$0.3>z_{wt} \ge 0 \f$
  !!
  !! \f$ k_{a}=  k_{4}{(e}^{k_{5}z_{P-}}e^{k_{5}z_{wt}})+k_{6}\left( z_p-z_{wt} \right)\f$ for \f$ z_{wt}\ge 0.3 \f$
  !!
  !! \f$ C_{SOM, o}= 0.487\ast (k_{9}z_{wt}^{2}+k_{10}z_{wt}) \f$
  !!
  !!  \f$ C_{SOM, a}= C_{SOM}-C_{SOM, o} \f$
  !!
  !! where 0.487 is a parameter that converts from soil mass to soil C content.
  !! The variation of \f$k_{o}\f$ and \f$k_{a}\f$ with water table depth for
  !! bogs and fens is shown in Fig. 2 of (Wu et al. 2016) \cite Wu2016-zt. It will be noted that there is a sharp
  !! transition in decomposition rate at a depth of 0.3 m, reflecting the work of
  !! Frolking et al. (2001). As noted above, this value is widely accepted
  !! as a representative estimate of the depth dividing the acrotelm and catotelm.
  !! In reality, of course, this depth will vary among peatlands.
  !!
  !! As only organic soil is considered in peatlands, the peat soil C is updated
  !! from the humification (C\f$_{hum}\f$, kg C/ m\f$^2\f$ / day) and soil
  !! respiration from the oxic (\f$R_o\f$ in kg C/ m\f$^2\f$ / day) and
  !! anoxic (\f$R_a\f$ in kg C/ m\f$^2\f$ / day) components during the
  !! time step:
  !!
  !! \f$ \frac{{dC}_{SOM}}{dt} = C_{hum}-R_{o}-R_{a} \f$
  !!
  !! \f$C_{hum}\f$ is calculated as a PFT-dependent fraction of the decomposition
  !! rate. At the end of each time step, the peat depth (i.e. the depth
  !! of the organic soil) \f$z_p\f$ is updated from the updated peat C mass
  !! (C\f$_{SOM}\f$ in kg) by solving the quadratic equation
  !!
  !! \f$ z_{p}=
  !! \frac{{-k}_{10}+\sqrt{k_{10}+\frac{4k_{9}{C}_{SOM}}{0.487}}.
  !! }{2k_{9}} \f$
  !!
  !! The water table depth \f$z_{wt}\f$ is deduced by searching for a soil layer
  !! below, which the soil is saturated and above which the soil moisture is at or
  !! below the retention capacity with respect to gravitational drainage. Within
  !! this soil layer \f$j\f$, \f$z_{wt}\f$ is calculated as
  !!
  !! \f$ z_{wt}=z_{{b}, j}-\Delta z\left[
  !! \frac{\theta_{{l}, j}+\theta_{i, j}-\theta
  !! _{{ret}, j}}{\theta_{{p}, j}-\theta_{{ret}, j}} \right] \f$
  !!
  !! where \f$\Delta z \f$ is the thickness of soil layer (unit: m),
  !! \f$\theta_l \f$ and \f$\theta_i\f$ are the liquid and frozen water contents
  !! (unit, m\f$^3\f$ / m\f$^3\f$), \f$\theta_{ret}\f$ and \f$\theta_p\f$ are the
  !! water retention capacity and the porosity, and \f$z_b\f$ (unit: m) is the
  !! bottom depth of the soil layer.
  !!

  !> \file
end module peatlandsMod
