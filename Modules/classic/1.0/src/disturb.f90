!> \file
!> Calculates disturbance as both natural and human-influenced fires.
module disturbance_scheme

  implicit none

  ! Subroutines
  public  :: disturb
  public  :: burntobare

contains

  ! ------------------------------------------------------------------
  !> \ingroup disturbance_scheme_disturb
  !! @{
  !> Calculates whether fire occurs, burned area, amount of C emitted and litter generated
  !> @author Vivek Arora and Joe Melton

  subroutine disturb (thliq, THLW, THFC, uwind, useTracer, & ! In
                      vwind, lightng, fcancmx, isand, & ! In
                      rmatctem, ilg, il1, il2, sort, & ! In
                      grclarea, thice, popdin, lucemcom, & ! In
                      dofire, currlat, iday, fsnow, & ! In
                      stemmass, rootmass, gleafmas, bleafmas, litrmass, & ! In/out
                      tracerStemMass, tracerRootMass, tracerGLeafMass, tracerBLeafMass, tracerLitrMass, & ! In/Ou
                      stemltdt, rootltdt, glfltrdt, blfltrdt, & ! Out (Primary)
                      glcaemls, rtcaemls, stcaemls, & ! Out (Primary)
                      blcaemls, ltrcemls, burnfrac, & ! Out (Primary)
                      pstemmass, pgleafmass, emit_co2, emit_ch4, & ! Out (Primary)
                      emit_co, emit_nmhc, emit_h2, emit_nox, & ! Out (Secondary)
                      emit_n2o, emit_pm25, emit_tpm, emit_tc, & ! Out (Secondary)
                      emit_oc, emit_bc, burnvegf, bterm_veg, & ! Out (Secondary)
                      mterm_veg, lterm, smfunc_veg) ! Out (Secondary)

    !
    !     20  Nov 2018  - Remove extnprob from going in. It's redundant now because
    !     Vivek Arora     it's calculated inside as a f(POPDIN). Clean up arguments
    !                     list so that secondary output can be cleanly removed for AGCM.
    !                     Pass in ILG back again as an argument.
    !
    !     12  Jan 2016  - Change to allow more than 3 soil layers
    !     J. Melton
    !
    !     11  Jun 2015  - Clean up to make all calculations PFT dependent and add in a
    !     V. Arora        capability to do a fire count method for comparison with
    !                     our usual parameterization
    !
    !     26  Mar 2014  - Split subroutine into two and create module. Move all fcancmx
    !     J. Melton       adjustments into subroutine burntobare and call from competition
    !
    !     20  Feb 2014  - Adapt to deal with competition on. Bring in code that makes
    !     J. Melton       bare fractions from competition module. Moved parameters to
    !                     classicParams.f90
    !
    !     4   Jan 2014  - Convert to f90 and include a saturation effect for lightning
    !     J. Melton       strikes

    !     25  Jul 2013  - Add in module for common params, cleaned up code
    !     J. Melton

    !     24  Sep 2012  - add in checks to prevent calculation of non-present
    !     J. Melton       pfts

    !     09  May 2012  - addition of emission factors and revising of the
    !     J. Melton       fire scheme

    !     15  May 2003  - this subroutine calculates the litter generated
    !     V. Arora        and c emissions from leaves, stem, and root
    !                     components due to fire. c emissions from burned
    !                     litter are also estimated. at present no other
    !                     form of disturbance is modelled.


    use classicParams, only : ignd, icc, ican, zero, kk, pi, c2dom, crop, &
                              iccp2, standreplace, tolrance, bmasthrs_fire, &
                              lwrlthrs, hgrlthrs, parmlght, parblght, reparea, popdthrshld, &
                              f0, maxsprd, frco2glf, frco2blf, &
                              frltrglf, frltrblf, frco2stm, frltrstm, frco2rt, frltrrt, &
                              frltrbrn, emif_co2, emif_co, emif_ch4, emif_nmhc, emif_h2, &
                              emif_nox, emif_n2o, emif_pm25, emif_tpm, emif_tc, emif_oc, emif_bc, &
                              grass, extnmois_veg, extnmois_duff, iccp1, nol2pfts

    implicit none

    integer, intent(in) :: ilg !<
    integer, intent(in) :: il1 !< il1=1
    integer, intent(in) :: il2 !< il2=ilg
    integer, intent(in) :: isand(ilg,ignd) !<
    integer, intent(in) :: sort(icc) !< index for correspondence between 9 pfts and size 12 of parameters vectors
    integer, intent(in) :: iday
    logical, intent(in) :: dofire !< boolean, if true allow fire, if false no fire
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.

    real, intent(in) :: thliq(ilg,ignd)    !< liquid soil moisture content
    real, intent(in) :: THLW(ilg,ignd)   !< wilting point soil moisture content
    real, intent(in) :: THFC(ilg,ignd)  !< field capacity soil moisture content
    real, intent(in) :: uwind(ilg)         !< wind speed, \f$m/s\f$
    real, intent(in) :: vwind(ilg)         !< wind speed, \f$m/s\f$
    real, intent(in) :: fcancmx(ilg,icc)   !< fractional coverages of ctem's 9 pfts
    real, intent(in) :: lightng(ilg)       !< total \f$lightning, flashes/(km^2 . year)\f$ it is assumed that cloud
    !< to ground lightning is some fixed fraction of total lightning.
    real, intent(in) :: rmatctem(ilg,icc,ignd) !< fraction of roots in each soil layer for each pft
    real, intent(in) :: thice(ilg,ignd)   !< frozen soil moisture content over canopy fraction
    real, intent(in) :: popdin(ilg)       !< population density \f$(people / km^2)\f$
    real, intent(in) :: lucemcom(ilg)     !< land use change (luc) related combustion emission losses, \f$u-mol co2/m2.sec\f$
    real, intent(in) :: currlat(ilg)      !<
    real, intent(in) :: fsnow(ilg)        !< fraction of snow simulated by class

    real, intent(inout) :: stemmass(ilg,icc)  !< stem mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, intent(inout) :: rootmass(ilg,icc)  !< root mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, intent(inout) :: gleafmas(ilg,icc)  !< green leaf mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, intent(inout) :: bleafmas(ilg,icc)  !< brown leaf mass

    ! COMBAK PERLAY
    real, intent(inout) :: litrmass(ilg,iccp2)!< litter mass for each of the CTEM pfts
    ! real, intent(inout) :: litrmass(:,:,:)!< litter mass for each of the CTEM pfts
    ! COMBAK PERLAY

    real, intent(inout) :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$tracer C units/m^2\f$

    real, intent(out) :: stemltdt(ilg,icc) !< stem litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(out) :: rootltdt(ilg,icc) !< root litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(out) :: glfltrdt(ilg,icc) !< green leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(out) :: pstemmass(ilg,icc)
    real, intent(out) :: pgleafmass(ilg,icc)
    real, intent(out) :: burnvegf(ilg,icc)    !< per PFT fraction burned of that PFTs area
    real, intent(out) :: bterm_veg(ilg,icc)    !< biomass fire probability term, Vivek
    real, intent(out) :: mterm_veg(ilg,icc)    !< moisture fire probability term, Vivek
    real, intent(out) :: lterm(ilg)           !< lightning fire probability term
    real, intent(out) :: smfunc_veg(ilg,icc)   !< soil moisture dependence on fire spread rate, Vivek

    !     note the following c burned will be converted to a trace gas
    !     emission or aerosol on the basis of emission factors.
    real, intent(out) :: glcaemls(ilg,icc) !< green leaf carbon emission losses, \f$kg c/m^2\f$
    real, intent(out) :: rtcaemls(ilg,icc) !< root carbon emission losses, \f$kg c/m^2\f$
    real, intent(out) :: stcaemls(ilg,icc) !< stem carbon emission losses, \f$kg c/m^2\f$
    real, intent(out) :: ltrcemls(ilg,icc) !< litter carbon emission losses, \f$kg c/m^2\f$
    real, intent(out) :: blfltrdt(ilg,icc) !< brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(out) :: blcaemls(ilg,icc) !< brown leaf carbon emission losses, \f$kg c/m^2\f$
    real, intent(out) :: burnfrac(ilg)     !< total areal :: fraction burned, (%)

    !     emitted compounds from biomass burning (kg {species} / m2 / s)
    real, intent(out) :: emit_co2(ilg,icc) !< carbon dioxide
    real, intent(out) :: emit_co(ilg,icc)  !< carbon monoxide
    real, intent(out) :: emit_ch4(ilg,icc) !< methane
    real, intent(out) :: emit_nmhc(ilg,icc)!< non-methane hydrocarbons
    real, intent(out) :: emit_h2(ilg,icc)  !< hydrogen gas
    real, intent(out) :: emit_nox(ilg,icc) !< nitrogen oxides
    real, intent(out) :: emit_n2o(ilg,icc) !< nitrous oxide
    real, intent(out) :: emit_pm25(ilg,icc)!< particulate matter less than 2.5 um in diameter
    real, intent(out) :: emit_tpm(ilg,icc) !< total particulate matter
    real, intent(out) :: emit_tc(ilg,icc)  !< total carbon
    real, intent(out) :: emit_oc(ilg,icc)  !< organic carbon
    real, intent(out) :: emit_bc(ilg,icc)  !< black carbon

    ! Local
    integer :: i, j, k, m, k1, k2, n
    real :: prbfrhuc(ilg)      !< probability of fire due to human causes
    real :: extnprob(ilg)      !< fire extinguishing probability
    real :: burnarea(ilg)     !< total area burned, \f$km^2\f$
    real :: biomass(ilg,icc)  !< total biomass for fire purposes
    real :: drgtstrs(ilg,icc) !< soil dryness factor for pfts
    real :: betadrgt(ilg,ignd)!< dryness term for soil layers
    real :: avgdryns(ilg)     !< avg. dryness over the vegetated fraction
    real :: fcsum(ilg)        !< total vegetated fraction
    real :: avgbmass(ilg)     !< avg. veg. biomass over the veg. fraction of grid cell
    real :: c2glgtng(ilg)     !< cloud-to-ground lightning
    real :: betalght(ilg)     !< 0-1 lightning term
    real :: y(ilg)            !< logistic dist. for fire prob. due to lightning
    real :: temp(ilg)         !<
    real :: betmsprd(ilg)     !< beta moisture for calculating fire spread rate
    real :: smfunc(ilg)       !< soil moisture function used for fire spread rate
    real :: wind(ilg)         !< wind speed in km/hr
    real :: wndfunc(ilg)      !< wind function for fire spread rate
    real :: sprdrate(ilg)     !< fire spread rate
    real :: lbratio(ilg)      !< length to breadth ratio of fire
    real :: arbn1day(ilg)     !< area burned in 1 day, \f$km^2\f$
    real :: areamult(ilg)     !< multiplier to find area burned
    real :: vegarea(ilg)      !< total vegetated area in a grid cell
    real :: grclarea(ilg)     !< gcm grid cell area, \f$km^2\f$
    real :: tot_emit          !< sum of all pools to be converted to emissions/aerosols \f$(g c/m^2)\f$
    real :: tot_emit_dom      !< tot_emit converted to \f$kg dom / m^2\f$
    real :: natural_ignitions(ilg) !<
    real :: anthro_ignitions(ilg)  !<
    real :: fire_supp(ilg)         !<
    real :: num_ignitions(ilg,icc)!<
    real :: num_fires(ilg,icc)    !<
    real :: hb_interm    !< interm calculation
    real :: hbratio(ilg) !< head to back ratio of ellipse
    real :: surface_duff_f(ilg)  !< fraction of biomass that is in the surface duff (grass brown leaves + litter)
    real :: pftareab(ilg,icc)    !< areas of different pfts in a grid cell, before fire, \f$km^2\f$
    !< pft area before fire \f$(km^2)\f$
    real :: ymin, ymax, slope
    real :: soilterm  !< temporary variable
    real :: duffterm  !< temporary variable
    real :: extn_par1 !< parameter used in calculation of fire extinguishing probability
    real :: duff_frac_veg(ilg,icc)!< duff fraction for each PFT, Vivek
    real :: probfire_veg(ilg,icc) !< PFT fire probability term, Vivek
    real :: sprdrate_veg(ilg,icc) !< per PFT fire spread rate
    real :: arbn1day_veg(ilg,icc) !< per PFT area burned in 1 day
    real :: burnarea_veg(ilg,icc) !< per PFT area burned over fire duration
    logical :: fire_veg(ilg,icc)  !< fire occuring logical, Vivek
    real :: soilterm_veg, duffterm_veg, betmsprd_veg, betmsprd_duff      ! temporary variables
    real :: tracerLitTemp, tracerEmitTemp ! temp tracer variables
    real :: fractionPFTburned  !< temp variable.

    !> initialize required arrays to zero, or assign value

    do j = 1,icc
      do i = il1,il2
        stemltdt(i,j) = 0.0
        rootltdt(i,j) = 0.0
        glfltrdt(i,j) = 0.0
        blfltrdt(i,j) = 0.0
        biomass(i,j) = 0.0
        drgtstrs(i,j) = 0.0
        glcaemls(i,j) = 0.0
        blcaemls(i,j) = 0.0
        stcaemls(i,j) = 0.0
        rtcaemls(i,j) = 0.0
        ltrcemls(i,j) = 0.0

        emit_co2(i,j) = 0.0
        emit_co(i,j) = 0.0
        emit_ch4(i,j) = 0.0
        emit_nmhc(i,j) = 0.0
        emit_h2(i,j) = 0.0
        emit_nox(i,j) = 0.0
        emit_n2o(i,j) = 0.0
        emit_pm25(i,j) = 0.0
        emit_tpm(i,j) = 0.0
        emit_tc(i,j) = 0.0
        emit_oc(i,j) = 0.0
        emit_bc(i,j) = 0.0
        burnvegf(i,j) = 0.0

        bterm_veg(i,j) = 0.0
        mterm_veg(i,j) = 0.0
        probfire_veg(i,j) = 0.0
        fire_veg(i,j) = .false.
        duff_frac_veg(i,j) = 0.0
        smfunc_veg(i,j) = 0.0
        sprdrate_veg(i,j) = 0.0
        arbn1day_veg(i,j) = 0.0
        burnarea_veg(i,j) = 0.0
      end do ! loop 150
    end do ! loop 140

    do k = 1,ignd
      do i = il1,il2
        betadrgt(i,k) = 1.0
      end do ! loop 170
    end do ! loop 160

    do i = il1,il2

      avgbmass(i) = 0.0
      avgdryns(i) = 0.0
      fcsum(i) = 0.0
      c2glgtng(i) = 0.0
      betalght(i) = 0.0
      y(i) = 0.0
      lterm(i) = 0.0
      burnarea(i) = 0.0
      burnfrac(i) = 0.0
      betmsprd(i) = 0.0
      smfunc(i) = 0.0
      wind(i) = 0.0
      wndfunc(i) = 0.0
      sprdrate(i) = 0.0
      lbratio(i) = 0.0
      arbn1day(i) = 0.0
      areamult(i) = 0.0
      vegarea(i) = 0.0
      surface_duff_f(i) = 0.0

    end do ! loop 180

    if (dofire) then

    !> initialization ends

    !> Find pft areas before
    do j = 1,icc
      do i = il1,il2
        pftareab(i,j) = fcancmx(i,j) * grclarea(i)  !> area in \f$km^2\f$
      end do ! loop 83
    end do ! loop 82

    !     ------------------------------------------------------------------

    !> Find the probability of fire as a product of three functions
    !> with dependence on total biomass, soil moisture, and lightning

    !> ### 1 Dependence on total biomass

    do j = 1,icc
      do i = il1,il2
        if (.not. crop(j)) then !> don't allow it to bring in crops since they are not allowed to burn.

          !> Root biomass is not used to initiate fire. For example if
          !> the last fire burned all grass leaves, and some of the roots
          !> were left, its unlikely these roots could catch fire.
          !> Here we ignore the LUC litrmass on iccp2 and the litter on the bare ground
          !> (iccp1). We only consider the litrmass on layer 1 as the rest are buried.
          ! COMBAK PERLAY
          biomass(i,j) = gleafmas(i,j) + bleafmas(i,j) + stemmass(i,j) &
                         + litrmass(i,j)
          ! biomass(i,j) = gleafmas(i,j) + bleafmas(i,j) + stemmass(i,j) &
          !              + litrmass(i,j,1)
          ! COMBAK PERLAY

          !> Find average biomass over the vegetated fraction
          avgbmass(i) = avgbmass(i) + biomass(i,j) * fcancmx(i,j)

          !> Sum up the vegetated area
          fcsum(i) = fcsum(i) + fcancmx(i,j)

        end if
      end do ! loop 210
    end do ! loop 200

    !> calculate bterm for individual PFTs as well

    do j = 1,icc
      do i = il1,il2
        bterm_veg(i,j) = min(1.0,max(0.0,(biomass(i,j) - bmasthrs_fire(1))/(bmasthrs_fire(2) - bmasthrs_fire(1))))
      end do ! loop 252
    end do ! loop 251

    !> ### 2 Dependence on soil moisture

    !> This is calculated in a way such that the more dry the root zone
    !> of a pft type is, and more fractional area is covered with that
    !> pft, the more likely it is that fire will get started. that is
    !> the dryness factor is weighted by fraction of roots in soil
    !> layers, as well as according to the fractional coverage of
    !> different pfts. the assumption here is that if there is less
    !> moisture in root zone, then it is more likely the vegetation
    !> will be dry and thus the likeliness of fire is more.
    !!
    !! First find the dryness factor for each soil layer.
    !!
    !! If there is snow on the ground, similarly if not soil,
    !! do not allow fire so set betadrgt to
    !! 0 for all soil layers otherwise calculate as per normal.
    do i = il1,il2
      if (fsnow(i) == 0.) then
        do j = 1,ignd
          if (isand(i,j) > - 2) then
            betadrgt(i,j) = min(1.0,max(0.0,(thliq(i,j) + thice(i,j) - THLW(i,j))/(THFC(i,j) - THLW(i,j))))
          else
            betadrgt(i,j) = 1.0
          end if
        end do
      end if
    end do

    !> Now find weighted value of this dryness factor averaged over the rooting depth, for each pft

    do j = 1,icc
      do i = il1,il2
        if (.not. crop(j)) then

          drgtstrs(i,j) =  sum((betadrgt(i,:)) * rmatctem(i,j,:))
          drgtstrs(i,j) = min(1.0,max(0.0,drgtstrs(i,j)/sum(rmatctem(i,j,:))))

          !> Next find this dryness factor averaged over the vegetated fraction
          !! \f$avgdryns(i) = avgdryns(i) + drgtstrs(i,j)*fcancmx(i,j)\f$
          !!
          !! The litter and brown leaves are not affected by the soil water potential
          !! therefore they will react only to the moisture conditions (our proxy here
          !! is the upper soil moisture). If it is dry they increase the probability of
          !! fire corresponding to the proportion of total C they contribute. Only allow
          !! if there is no snow.
          !           if (biomass(i,j) > 0. .and. fsnow(i) == 0.) then
          !             ! The surface duff calculation ignores the litter on the bare fraction.
          !             surface_duff_f(i) = surface_duff_f(i) + (bleafmas(i,j)+litrmass(i,j)) &
          !                                                      /biomass(i,j) * fcancmx(i,j)
          !
          !           end if

        end if
      end do ! loop 330
    end do ! loop 320

    !> Use average root zone vegetation dryness to find likelihood of
    !! fire due to moisture.
    !!
    !! calculate mterm for each PFT

    do j = 1,icc
      do i = il1,il2
        !> duff fraction for each PFT, Vivek
        !> We only consider the first layer litter
        if (biomass(i,j) > 0. .and. fsnow(i) == 0.) then
          ! COMBAK PERLAY
          duff_frac_veg(i,j) = (bleafmas(i,j) + litrmass(i,j)) / biomass(i,j)
          ! duff_frac_veg(i,j) = (bleafmas(i,j)+litrmass(i,j,1)) / biomass(i,j)
          ! COMBAK PERLAY
        end if

        !> \f$drgtstrs(i,j)\f$ is \f$\phi_{root}\f$ in Melton and Arora GMDD (2015) paper
        soilterm_veg = 1.0 - tanh((1.75 * drgtstrs(i,j)/extnmois_veg) ** 2)
        duffterm_veg = 1.0 - tanh((1.75 * betadrgt(i,1)/extnmois_duff) ** 2)

        if (fcancmx(i,j) > zero) then
          mterm_veg(i,j) = soilterm_veg * (1. - duff_frac_veg(i,j)) + duffterm_veg * duff_frac_veg(i,j)
        else
          mterm_veg(i,j) = 0.0   !> no fire likelihood due to moisture if no vegetation
        end if

      end do ! loop 382
    end do ! loop 381

    !> ### 3 dependence on lightning
    !!
    !! Dependence on lightning is modelled in a simple way which implies that
    !! a large no. of lightning flashes are more likely to cause fire than
    !! few lightning flashes.

    do i = il1,il2

      !> New approximation of Price and Rind equation. It was developed from a more complete dataset than Prentice and Mackerras.
      !! Lightning comes in in units of \f$flashes/km^2/yr\f$ so divide by 12 to make per month.
      ! c2glgtng(i)=lightng(i) ! FLAG FireMIP lightning is already C2G ! If using 'normal' lightning climatology use line below.
      c2glgtng(i) = 0.219913 * exp(0.0058899 * abs(currlat(i))) * lightng(i)

      betalght(i) = min(1.0,max(0.0,(c2glgtng(i) - lwrlthrs)/(hgrlthrs - lwrlthrs)))
      y(i) = 1.0/( 1.0 + exp((parmlght - betalght(i))/parblght) )

      ! No need to calculate each time, once settled on parameters, precalc and moved into a parameter. JM. Feb 19 2014.
      ymin = 1.0/( 1.0 + exp((parmlght - 0.0)/parblght) )
      ymax = 1.0/( 1.0 + exp((parmlght - 1.0)/parblght) )
      slope = abs(0.0 - ymin) + abs(1.0 - ymax)

      temp(i) = y(i) + (0.0 - ymin) + betalght(i) * slope

      !> Determine the probability of fire due to human causes
      !! this is based upon the population density from the .popd read-in file

      prbfrhuc(i) = min(1.0,(popdin(i)/popdthrshld) ** 0.43) ! From Kloster et al. (2010)

      lterm(i) = max(0.0,min(1.0,temp(i) + (1.0 - temp(i)) * prbfrhuc(i) ))

      !> ----------------------- Number of fire calculations ----------------------\\
      !>
      !> This is not used in CTEM in general but we keep the code here for testing purposes
      !!
      !! calculate natural and anthorpogenic ignitions/km2.day
      !! the constant 0.25 assumes not all c2g lightning hits causes ignitions, only 0.25 of them do
      !! the constant (1/30.4) converts c2g lightning from flashes/km2.month to flashes/km2.day
      !! MAKE SURE LIGHTNING IS IN UNITS OF FLASHES/KM2.MONTH
      !!
      !! Eqs. (4) and (5) of Li et al. 2012 doi:10.5194/bg-9-2761-2012 + also see corrigendum

      ! natural_ignitions(i)=c2glgtng(i) * 0.25 * (1/30.4)
      ! anthro_ignitions(i)=9.72E-4 * (1/30.4) * popdin(i) * (6.8 * (popdin(i)**(-0.6)) )

      !> calculate fire suppression also as a function of population density.
      !! Li et al. (2012) formulation is quite similar to what we already use based
      !! on Kloster et al. (2010, I think) but we use Kloster's formulation together
      !! with our fire extingishing probability. Here, the fire suppression formulation
      !! is just by itself

      ! fire_supp(i)=0.99 - ( 0.98*exp(-0.025*popdin(i)*(i)) )

      !> ----------------------- Number of fire calculations ----------------------//
    end do ! loop 400

    !> calculate fire probability for each PFT. Recall that lightning term is still grid averaged

    do j = 1,icc
      do i = il1,il2
        probfire_veg(i,j) = bterm_veg(i,j) * mterm_veg(i,j) * lterm(i)
        if (probfire_veg(i,j) > zero) fire_veg(i,j) = .true.

        !> ----------------------- Number of fire calculations ----------------------\\
        !!
        !! This is not used in CTEM in general but we keep the code here for testing purposes
        !!
        !! calculate total number of ignitions based natural and anthorpogenic ignitions
        !! for the whole grid cell

        ! num_ignitions(i,j) = ( natural_ignitions(i) + anthro_ignitions(i) ) * fcancmx(i,j)*grclarea(i)

        !> finally calculate number of fire,noting that not all ignitions turn into fire
        !! because moisture and biomass may now allow that to happen,and some of those
        !! will be suppressed due to fire fighting efforts

        ! num_fires(i,j) = num_ignitions(i,j)*(1-fire_supp(i))*bterm_veg(i,j)*mterm_veg(i,j)
        !>
        !> ----------------------- Number of fire calculations ----------------------//

      end do ! loop 422
    end do ! loop 421

    !> Calculate area burned for each PFT, make sure it's not greater than the
    !! area available, then find area and fraction burned for the whole gridcell

    do j = 1,icc
      do i = il1,il2
        if (fire_veg(i,j) ) then

          !> soil moisture dependence on fire spread rate
          betmsprd_veg = (1. - min(1., (drgtstrs(i,j)/extnmois_veg) )) ** 2
          betmsprd_duff = (1. - min(1., (betadrgt(i,1)/extnmois_duff) )) ** 2
          smfunc_veg(i,j) = betmsprd_veg * (1 - duff_frac_veg(i,j)) + betmsprd_duff * duff_frac_veg(i,j)

          !> wind speed, which is gridcell specific
          wind(i) = sqrt(uwind(i) ** 2.0 + vwind(i) ** 2.0)
          wind(i) = wind(i) * 3.60     ! change m/s to km/hr

          !> Length to breadth ratio from Li et al. (2012)
          lbratio(i) = 1.0 + 10.0 * (1.0 - exp( - 0.06 * wind(i)))

          !> head to back ratio from Li et al. (2012).
          hb_interm = (lbratio(i) ** 2 - 1.0) ** 0.5
          hbratio(i) = (lbratio(i) + hb_interm) / (lbratio(i) - hb_interm)

          !> dependence of spread rate on wind
          wndfunc(i) = (2.0 * lbratio(i)) / (1.0 + 1.0 / hbratio(i)) * f0

          !> fire spread rate per PFT
          n = sort(j)
          sprdrate_veg(i,j) = maxsprd(n) * smfunc_veg(i,j) * wndfunc(i)

          !> area burned in 1 day for that PFT
          arbn1day_veg(i,j) = (pi * 24.0 * 24.0 * sprdrate_veg(i,j) ** 2)/(4.0 * lbratio(i)) * (1.0 + 1.0 / hbratio(i)) ** 2

          !> fire extinguishing probability as a function of grid-cell averaged population density

          !> account for low suppression in Savanna regions, see above for
          !> increase in ignition due to cultural practices
          extn_par1 = - 0.015 ! Value based on Vivek's testing. Jul 14 2016. old = -0.025 FLAG MOVE TO PARAM !

          ! change 0.9 to 1.0 in eqn. A78 of Melton and Arora, 2016, GMD competition paper
          extnprob(i) = max(0.0,1.0 - exp(extn_par1 * popdin(i)))
          extnprob(i) = 0.5 + extnprob(i)/2.0

          !> area multipler to calculate area burned over the duration of the fire
          areamult(i) = ((1.0 - extnprob(i)) * (2.0 - extnprob(i)))/ extnprob(i) ** 2

          !> per PFT area burned, \f$km^2\f$
          burnarea_veg(i,j) = arbn1day_veg(i,j) * areamult(i) * (grclarea(i) * fcancmx(i,j) * probfire_veg(i,j))/reparea


          !          ------- Area burned based on number of fire calculations ----------------------\\

          !       This is not used in CTEM in general but we keep the code here for testing purposes

          !         the constant 4 is suppose to address the fact that Li et al. (2012) suggested to
          !         double the fire spread rates. However, if we do that than our usual calculations
          !         based on CTEM's original parameterization will be affected. Rather than do that
          !         we just use a multiplier of 4, since doubling fire spread rates means 4 times the
          !         area burned
          !
          !          burnarea_veg(i,j)=arbn1day_veg(i,j)*num_fires(i,j)*2.0  ! flag test was 4 !
          !
          !          ------- Area burned based on number of fire calculations ----------------------//

          !         if area burned greater than area of PFT, set it to area of PFT

          if (burnarea_veg(i,j) > grclarea(i) * fcancmx(i,j) ) then
            burnarea_veg(i,j) = grclarea(i) * fcancmx(i,j)
          end if
        end if
      end do ! loop 501
    end do ! loop 500

    !> Calculate gridcell area burned and fraction

    do j = 1,icc
      do i = il1,il2
        burnarea(i) = burnarea(i) + burnarea_veg(i,j)
      end do ! loop 511
    end do ! loop 510

    do i = il1,il2
      burnfrac(i) = burnarea(i)/grclarea(i)
    end do ! loop 512

    !> Finally estimate amount of litter generated from each pft, and
    !> each vegetation component (leaves, stem, and root) based on their
    !> resistance to combustion. Update the veg pools due to combustion.

    do j = 1,icc
      n = sort(j)
      do i = il1,il2
        if (pftareab(i,j) > zero) then

          !> Set aside these pre-disturbance stem and root masses for use
          !> in burntobare subroutine.
          pstemmass(i,j) = stemmass(i,j)
          pgleafmass(i,j) = gleafmas(i,j)

          fractionPFTburned = burnarea_veg(i,j) / pftareab(i,j)

          glfltrdt(i,j) = frltrglf(n) * gleafmas(i,j) * fractionPFTburned
          blfltrdt(i,j) = frltrblf(n) * bleafmas(i,j) * fractionPFTburned
          stemltdt(i,j) = frltrstm(n) * stemmass(i,j) * fractionPFTburned
          rootltdt(i,j) = frltrrt(n)  * rootmass(i,j) * fractionPFTburned
          glcaemls(i,j) = frco2glf(n) * gleafmas(i,j) * fractionPFTburned
          blcaemls(i,j) = frco2blf(n) * bleafmas(i,j) * fractionPFTburned
          stcaemls(i,j) = frco2stm(n) * stemmass(i,j) * fractionPFTburned
          rtcaemls(i,j) = frco2rt(n)  * rootmass(i,j) * fractionPFTburned

          ! The burned litter comes from the first layer
          ! COMBAK PERLAY
          ltrcemls(i,j) = frltrbrn(n) * litrmass(i,j) * fractionPFTburned
          ! ltrcemls(i,j)= frltrbrn(n) *litrmass(i,j,1) * fractionPFTburned
          ! COMBAK PERLAY

          !> Update the pools:
          gleafmas(i,j) = gleafmas(i,j) - glfltrdt(i,j) - glcaemls(i,j)
          bleafmas(i,j) = bleafmas(i,j) - blfltrdt(i,j) - blcaemls(i,j)
          stemmass(i,j) = stemmass(i,j) - stemltdt(i,j) - stcaemls(i,j)
          rootmass(i,j) = rootmass(i,j) - rootltdt(i,j) - rtcaemls(i,j)

          ! The burned litter is placed on the top litter layer except for the root litter which
          ! goes into litter layers according to the root distribution.
          ! COMBAK PERLAY
          litrmass(i,j) = litrmass(i,j) + glfltrdt(i,j) + blfltrdt(i,j) + stemltdt(i,j) &
                          + rootltdt(i,j) - ltrcemls(i,j)
          !
          ! litrmass(i,j,1) = litrmass(i,j,1) + glfltrdt(i,j) + blfltrdt(i,j) + stemltdt(i,j) &
          !                 + rootltdt(i,j) * rmatctem(i,j,1) - ltrcemls(i,j)
          !
          ! do k = 2,ignd
          !   litrmass(i,j,k) = litrmass(i,j,k) + rootltdt(i,j) * rmatctem(i,j,k)
          ! end do
          ! COMBAK PERLAY
          ! tmpzzz=0.
          if (useTracer > 0) then

            ! The burned litter comes out of the first litter layer. Update it first here
            ! as later we change the tracerLitrMass value so doing it later would change
            ! the result.
            tracerLitTemp = frltrbrn(n) * tracerLitrMass(i,j,1) * fractionPFTburned
            tracerLitrMass(i,j,1) = tracerLitrMass(i,j,1) - tracerLitTemp

            ! Update the tracer pools for fire losses. Need to calculte the litter and emissions
            ! terms since they are dependent upon pool size.
            tracerLitTemp = frltrglf(n) * tracerGLeafMass(i,j) * fractionPFTburned
            tracerEmitTemp = frco2glf(n) * tracerGLeafMass(i,j) * fractionPFTburned
            tracerGLeafMass(i,j) = tracerGLeafMass(i,j) - tracerLitTemp - tracerEmitTemp

            ! Add the green leaf litter to top litter layer
            tracerLitrMass(i,j,1) = tracerLitrMass(i,j,1) + tracerLitTemp

            tracerLitTemp = frltrblf(n) * tracerBLeafMass(i,j) * fractionPFTburned
            tracerEmitTemp = frco2blf(n) * tracerBLeafMass(i,j) * fractionPFTburned
            tracerBLeafMass(i,j) = tracerBLeafMass(i,j) - tracerLitTemp - tracerEmitTemp

            ! Add the brown leaf litter to top litter layer
            tracerLitrMass(i,j,1) = tracerLitrMass(i,j,1) + tracerLitTemp

            tracerLitTemp = frltrstm(n) * tracerStemMass(i,j) * fractionPFTburned
            tracerEmitTemp = frco2stm(n) * tracerStemMass(i,j) * fractionPFTburned
            tracerStemMass(i,j) = tracerStemMass(i,j) - tracerLitTemp - tracerEmitTemp

            ! Add the stem litter to top litter layer
            tracerLitrMass(i,j,1) = tracerLitrMass(i,j,1) + tracerLitTemp

            tracerLitTemp = frltrrt(n) * tracerRootMass(i,j) * fractionPFTburned
            tracerEmitTemp = frco2rt(n) * tracerRootMass(i,j) * fractionPFTburned
            tracerRootMass(i,j) = tracerRootMass(i,j) - tracerLitTemp - tracerEmitTemp

            ! Add the root litter,which goes into litter layers according to the root distribution.
            do k = 1,ignd
              tracerLitrMass(i,j,k) = tracerLitrMass(i,j,k) + tracerLitTemp * rmatctem(i,j,k)
            end do

          end if

          !> Output the burned area per PFT (the units here are burned fraction of each PFTs area. So
          !> if a PFT has 50% gridcell cover and 50% of that burns it will have a burnvegf of 0.5 (which
          !> then translates into a gridcell fraction of 0.25). This units is for consistency
          !! outside of this subroutine.
          burnvegf(i,j) = burnarea_veg(i,j) /  pftareab(i,j)

        end if
      end do ! loop 530
    end do ! loop 520
  end if ! dofire check.
  
    !> Even if no fire we still enter here and perform the calculations for the emissions
    !> since we might have some from luc.
    !>
    !> We also estimate \f$CO_2\f$ emissions from each
    !> of these components. Note that the litter which is generated due
    !! to disturbance is uniformly distributed over the entire area of
    !! a given pft, and this essentially thins the vegetation biomass.
    !! If PFTCompetition is not on, this does not change the vegetation fractions,
    !! if competition is on a fraction will become bare. That is handled in
    !! burntobare subroutine called from competition subroutine.
    !!
    do j = 1,icc ! loop 620
      n = sort(j)
      do i = il1,il2 ! loop 630

        !> Calculate the emissions of trace gases and aerosols based upon how much plant matter was burnt

        !> Sum all pools that will be converted to emissions/aerosols \f$(g c/m^2 /day)\f$
        tot_emit = (glcaemls(i,j) + blcaemls(i,j) + rtcaemls(i,j) + stcaemls(i,j) + ltrcemls(i,j)) * 1000.0

        !> Add in the emissions due to luc fires (deforestation)
        !! the luc emissions are converted from \f$umol co_2 m-2 s-1 to
        !! g c m-2\f$ (day-1) before adding to tot_emit. NOTE: This tot_emit is
        !! not what is used to calculate NBP so the LUC emissions are accounted for
        !! in ctem.f90, not here. This here is for model outputting.
        tot_emit = tot_emit + (lucemcom(i) / 963.62 * 1000.0)

        !> Convert burnt plant matter from carbon to dry organic matter using
        !> a conversion factor, assume all parts of the plant has the same
        !> ratio of carbon to dry organic matter. units: \f$kg dom / m^2 / day\f$
        tot_emit_dom = tot_emit / c2dom

        !> Convert the dom to emissions/aerosols using emissions factors units: \f$g species / m^2 /d\f$
        !!
        !! Also convert units to \f$kg species / m^2 / s^{-1}\f$
        ! g {species} / m2/d * d/86400s * 1kg/1000g = kg {species} / m2 / s

        emit_co2(i,j)  = emif_co2(n) * tot_emit_dom * 1.1574e-8
        emit_co(i,j)   = emif_co(n)  * tot_emit_dom * 1.1574e-8
        emit_ch4(i,j)  = emif_ch4(n) * tot_emit_dom * 1.1574e-8
        emit_nmhc(i,j) = emif_nmhc(n) * tot_emit_dom * 1.1574e-8
        emit_h2(i,j)   = emif_h2(n) * tot_emit_dom * 1.1574e-8
        emit_nox(i,j)  = emif_nox(n) * tot_emit_dom * 1.1574e-8
        emit_n2o(i,j)  = emif_n2o(n) * tot_emit_dom * 1.1574e-8
        emit_pm25(i,j) = emif_pm25(n) * tot_emit_dom * 1.1574e-8
        emit_tpm(i,j)  = emif_tpm(n) * tot_emit_dom * 1.1574e-8
        emit_tc(i,j)   = emif_tc(n) * tot_emit_dom * 1.1574e-8
        emit_oc(i,j)   = emif_oc(n) * tot_emit_dom * 1.1574e-8
        emit_bc(i,j)   = emif_bc(n) * tot_emit_dom * 1.1574e-8

      end do ! loop 630
    end do ! loop 620

    ! FLAG for the optimization of popd effect on fire I am taking the lterm out as the 'temp' var. So I have made tmp be dimension
    ! ilg and I overwrite the lterm (which is written to an output file) here in its place.

    lterm = temp

  end subroutine disturb
  !> @}
  ! ------------------------------------------------------------------------------------

  !> \ingroup disturbance_scheme_burntobare
  !! @{
  !> Update fractional coverages of pfts to take into account the area
  !! burnt by fire. Adjust all pools with new densities in their new
  !! areas and increase bare fraction. And while we are doing this
  !! also run a small check to make sure grid averaged quantities do not get messed up.
  !> @author Joe Melton

  subroutine burntobare (il1, il2, nilg, sort, pvgbioms, pgavltms, pgavscms, & ! In
                         burnvegf, pstemmass, pgleafmass, useTracer, & ! In
                         fcancmx, stemmass, rootmass, gleafmas, bleafmas, & ! In/Out
                         litrmass, soilcmas, nppveg, tracerLitrMass, tracerSoilCMass, & ! In/Out
                         tracerGLeafMass, tracerBLeafMass, tracerStemMass, tracerRootMass) ! In/Out



    use classicParams, only : crop, icc, seed, standreplace, grass, zero, &
                              iccp1, tolrance, numcrops, iccp2, ignd

    implicit none

    integer, intent(in) :: il1
    integer, intent(in) :: il2
    integer, intent(in) :: nilg       !< no. of grid cells in latitude circle (this is passed in as either ilg or nlat depending on mos/comp)
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    integer, dimension(icc), intent(in) :: sort             !< index for correspondence between 9 ctem pfts and
    !< size 12 of parameter vectors
    real, dimension(nilg), intent(in) :: pvgbioms           !< initial veg biomass
    real, dimension(nilg), intent(in) :: pgavltms           !< initial litter mass
    real, dimension(nilg), intent(in) :: pgavscms           !< initial soil c mass
    real, dimension(nilg,icc), intent(in) :: burnvegf       !< per PFT fraction burned of that PFTs area
    real, dimension(nilg,icc), intent(in)    :: pstemmass   !< grid averaged stemmass prior to disturbance, \f$kg c/m^2\f$
    real, dimension(nilg,icc), intent(in)    :: pgleafmass  !< grid averaged rootmass prior to disturbance, \f$kg c/m^2\f$

    real, dimension(nilg,icc), intent(inout) :: fcancmx     !< initial fractions of the ctem pfts
    real, dimension(nilg,icc), intent(inout) :: gleafmas    !< green leaf carbon mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, dimension(nilg,icc), intent(inout) :: bleafmas    !< brown leaf carbon mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, dimension(nilg,icc), intent(inout) :: stemmass    !< stem carbon mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, dimension(nilg,icc), intent(inout) :: rootmass    !< roots carbon mass for each of the 9 ctem pfts, \f$kg c/m^2\f$
    real, dimension(nilg,icc), intent(inout) :: nppveg      !< npp for individual pfts,  \f$u-mol co_2/m^2.sec\f$
    ! COMBAK PERLAY
    real, dimension(nilg,iccp2), intent(inout) :: soilcmas  !< soil carbon mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    real, dimension(nilg,iccp2), intent(inout) :: litrmass  !< litter carbon mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    ! real, dimension(nilg,iccp2,ignd), intent(inout) :: soilcmas  !< soil carbon mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    ! real, dimension(nilg,iccp2,ignd), intent(inout) :: litrmass  !< litter carbon mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    ! COMBAK PERLAY
    real, intent(inout) :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

    logical, dimension(nilg) :: shifts_occur      !< true if any fractions changed
    integer :: i, j, n, k
    real :: pftfraca_old
    real :: term                                  !< temp variable for change in fraction due to fire
    real, dimension(nilg) :: pbarefra             !< bare fraction prior to fire
    real, dimension(nilg) :: barefrac             !< bare fraction of grid cell
    ! COMBAK PERLAY
    real, dimension(nilg) :: litr_lost            !< litter that is transferred to bare
    real, dimension(nilg) :: soilc_lost           !< soilc that is transferred to bare
    ! real, dimension(nilg,ignd) :: litr_lost            !< litter that is transferred to bare
    ! real, dimension(nilg,ignd) :: soilc_lost           !< soilc that is transferred to bare
    ! COMBAK PERLAY
    real :: tracerLitrLost(nilg,ignd)       !< tracer litter that is transferred to bare
    real :: tracerSoilCLost(nilg,ignd)      !< tracer soilc that is transferred to bare
    real, dimension(nilg) :: vgbiomas_temp        !< grid averaged vegetation biomass for internal checks, \f$kg c/m^2\f$
    real, dimension(nilg) :: gavgltms_temp        !< grid averaged litter mass for internal checks, \f$kg c/m^2\f$
    real, dimension(nilg) :: gavgscms_temp        !< grid averaged soil c mass for internal checks, \f$kg c/m^2\f$
    real, dimension(nilg,icc) :: pftfracb         !< pft fractions before accounting for creation of bare ground
    real, dimension(nilg,icc) :: pftfraca         !< pft fractions after accounting for creation of bare ground
    real :: frac_chang                            !< pftfracb - pftfraca

    ! -----------------------------------------

    !> Do some initializations
    do i = il1,il2
      shifts_occur(i) = .false.
      pbarefra(i) = 1.0
      barefrac(i) = 1.0
      vgbiomas_temp(i) = 0.0
      gavgltms_temp(i) = 0.0
      gavgscms_temp(i) = 0.0
      ! COMBAK PERLAY
      litr_lost(i) = 0.0
      soilc_lost(i) = 0.0
      ! litr_lost(i,:)=0.0
      ! soilc_lost(i,:)=0.0
      ! COMBAK PERLAY
      tracerLitrLost(i,:) = 0.0
      tracerSoilCLost(i,:) = 0.0
    end do ! loop 10
    !>
    !>  Account for disturbance creation of bare ground. This occurs with relatively low
    !>  frequency and is PFT dependent. We must adjust the amount of bare ground created
    !>  to ensure that we do not increase the density of the remaining vegetation.

    do i = il1,il2
      do j = 1,icc
        if (.not. crop(j)) then

          n = sort(j)

          pbarefra(i) = pbarefra(i) - fcancmx(i,j)
          pftfracb(i,j) = fcancmx(i,j)

          pftfraca(i,j) = max(seed,fcancmx(i,j) - burnvegf(i,j) * standreplace(n))

          fcancmx(i,j) = pftfraca(i,j)

          barefrac(i) = barefrac(i) - fcancmx(i,j)

        else  ! crops

          pbarefra(i) = pbarefra(i) - fcancmx(i,j)
          barefrac(i) = barefrac(i) - fcancmx(i,j)

        end if

      end do ! loop 25
    end do ! loop 20

    do i = il1,il2 ! loop 40
      do j = 1,icc ! loop 50
        if (.not. crop(j)) then

          !> Test the pftfraca to ensure it does not cause densification of the exisiting biomass
          !> Trees compare the stemmass while grass compares the root mass.
          if (pftfraca(i,j) /= pftfracb(i,j)) then

            shifts_occur(i) = .true.

            term = pftfracb(i,j)/pftfraca(i,j)

            if (.not. grass(j)) then
              if (stemmass(i,j) * term > pstemmass(i,j) .and. pstemmass(i,j) > 0.) then
                !> the pstemmass is from before the fire occurred,i.e. no thinning !
                pftfraca_old = pftfraca(i,j)
                pftfraca(i,j) = max(seed,stemmass(i,j) * pftfracb(i,j) / pstemmass(i,j))
                fcancmx(i,j) = pftfraca(i,j)

                !> adjust the bare frac to accomodate for the changes
                barefrac(i) = barefrac(i) + pftfraca_old - pftfraca(i,j)

              end if
            else ! grasses

              if (gleafmas(i,j) * term > pgleafmass(i,j) .and. pgleafmass(i,j) > 0.) then
                !> the pgleafmass is from before the fire occurred,i.e. no thinning !
                pftfraca_old = pftfraca(i,j)
                pftfraca(i,j) = max(seed,gleafmas(i,j) * pftfracb(i,j) / pgleafmass(i,j))
                fcancmx(i,j) = pftfraca(i,j)

                !> adjust the bare frac to accomodate for the changes
                barefrac(i) = barefrac(i) + pftfraca_old - pftfraca(i,j)

              end if
            end if

            term = pftfracb(i,j)/pftfraca(i,j)
            gleafmas(i,j) = gleafmas(i,j) * term
            bleafmas(i,j) = bleafmas(i,j) * term
            stemmass(i,j) = stemmass(i,j) * term
            rootmass(i,j) = rootmass(i,j) * term
            nppveg(i,j)  = nppveg(i,j) * term

            if (useTracer > 0) then ! Now same operation for tracer
              tracerGLeafMass(i,j) = tracerGLeafMass(i,j) * term
              tracerBLeafMass(i,j) = tracerBLeafMass(i,j) * term
              tracerStemMass(i,j) = tracerStemMass(i,j) * term
              tracerRootMass(i,j) = tracerRootMass(i,j) * term
            end if

            !>
            !> Soil and litter carbon are treated such that we actually transfer the carbon
            !> to the bare fraction since it would remain in place as a location was devegetated
            !> In doing so we do not adjust the litter or soilc density on the
            !> remaining vegetated fraction. But we do adjust it on the bare fraction to ensure
            !> our carbon balance works out.
            frac_chang = pftfracb(i,j) - pftfraca(i,j)
            ! COMBAK PERLAY
            litr_lost(i) = litr_lost(i) + litrmass(i,j) * frac_chang
            soilc_lost(i) = soilc_lost(i) + soilcmas(i,j) * frac_chang
            do k = 1,ignd
              !  litr_lost(i,k)= litr_lost(i,k) + litrmass(i,j,k) * frac_chang
              !  soilc_lost(i,k)= soilc_lost(i,k) + soilcmas(i,j,k) * frac_chang
              if (useTracer > 0) then ! Now same operation for tracer
                tracerLitrLost(i,k) = tracerLitrLost(i,k) + tracerLitrMass(i,j,k) * frac_chang
                tracerSoilCLost(i,k) = tracerSoilCLost(i,k) + tracerSoilCMass(i,j,k) * frac_chang
              end if
            end do
            ! COMBAK PERLAY
            ! else
            ! no changes
          end if
        end if  ! crop
      end do ! loop 50
    end do ! loop 40

    do i = il1,il2 ! loop 100
      if (shifts_occur(i)) then !> only do checks if we actually shifted fractions here.

        if (barefrac(i) >= zero .and. barefrac(i) > pbarefra(i)) then
          ! COMBAK PERLAY
          litrmass(i,iccp1) = (litrmass(i,iccp1) * pbarefra(i) + litr_lost(i)) / barefrac(i)
          soilcmas(i,iccp1) = (soilcmas(i,iccp1) * pbarefra(i) + soilc_lost(i)) / barefrac(i)
          do k = 1,ignd
            !  litrmass(i,iccp1,k) = (litrmass(i,iccp1,k)*pbarefra(i) + litr_lost(i,k)) / barefrac(i)
            !  soilcmas(i,iccp1,k) = (soilcmas(i,iccp1,k)*pbarefra(i) + soilc_lost(i,k)) / barefrac(i)
            if (useTracer > 0) then ! Now same operation for tracer
              tracerLitrMass(i,iccp1,k) = (tracerLitrMass(i,iccp1,k) * pbarefra(i) &
                                          + tracerLitrLost(i,k)) / barefrac(i)
              tracerSoilCMass(i,iccp1,k) = (tracerSoilCMass(i,iccp1,k) * pbarefra(i) &
                                           + tracerSoilCLost(i,k)) / barefrac(i)
            end if
          end do
          ! COMBAK PERLAY
        else if (barefrac(i) < 0.) then
          write(6, * )' In burntobare you have negative bare area,which should be impossible...'
          write(6, * )' bare is',barefrac(i),' original was',pbarefra(i)
          call errorHandler('disturb - burntobare', - 6)
        end if
      end if
    end do ! loop 100

    !> check if total grid average biomass density is same before and after adjusting fractions
    do i = il1,il2 ! loop 200

      if (shifts_occur(i)) then !> only do checks if we actually shifted fractions here.

        do j = 1,icc ! loop 250

          vgbiomas_temp(i) = vgbiomas_temp(i) + fcancmx(i,j) * (gleafmas(i,j) + &
                             bleafmas(i,j) + stemmass(i,j) + rootmass(i,j))
          ! COMBAK PERLAY
          gavgltms_temp(i) = gavgltms_temp(i) + fcancmx(i,j) * litrmass(i,j)
          gavgscms_temp(i) = gavgscms_temp(i) + fcancmx(i,j) * soilcmas(i,j)
          ! do k = 1,ignd
          !   gavgltms_temp(i)=gavgltms_temp(i)+fcancmx(i,j)*litrmass(i,j,k)
          !   gavgscms_temp(i)=gavgscms_temp(i)+fcancmx(i,j)*soilcmas(i,j,k)
          ! end do
          ! COMBAK PERLAY

        end do ! loop 250

        ! then add the bare ground in.
        ! COMBAK PERLAY
        gavgltms_temp(i) = gavgltms_temp(i) + barefrac(i) * litrmass(i,iccp1)
        gavgscms_temp(i) = gavgscms_temp(i) + barefrac(i) * soilcmas(i,iccp1)
        ! do k = 1,ignd  ! FLAG I think we can keep this as per grid like this. JM Feb 8 2016.
        !   gavgltms_temp(i)=gavgltms_temp(i)+ barefrac(i)*litrmass(i,iccp1,k)
        !   gavgscms_temp(i)=gavgscms_temp(i)+ barefrac(i)*soilcmas(i,iccp1,k)
        ! end do
        ! COMBAK PERLAY

        if (abs(vgbiomas_temp(i) - pvgbioms(i)) > tolrance) then
          write(6, * )'grid averaged biomass densities do not balance'
          write(6, * )'after fractional coverages are changed to take'
          write(6, * )'into account burn area'
          write(6, * )'vgbiomas_temp(',i,') = ',vgbiomas_temp(i)
          write(6, * )'pvgbioms(',i,') = ',pvgbioms(i)
          call errorHandler('disturb', - 7)
        end if
        if (abs(gavgltms_temp(i) - pgavltms(i)) > tolrance) then
          write(6, * )'grid averaged litter densities do not balance'
          write(6, * )'after fractional coverages are changed to take'
          write(6, * )'into account burn area'
          write(6, * )'gavgltms_temp(',i,') = ',gavgltms_temp(i)
          write(6, * )'pgavltms(',i,') = ',pgavltms(i)
          call errorHandler('disturb', - 8)
        end if
        if (abs(gavgscms_temp(i) - pgavscms(i)) > tolrance) then
          write(6, * )'grid averaged soilc densities do not balance'
          write(6, * )'after fractional coverages are changed to take'
          write(6, * )'into account burn area'
          write(6, * )'gavgscms_temp(',i,') = ',gavgscms_temp(i)
          write(6, * )'pgavscms(',i,') = ',pgavscms(i)
          call errorHandler('disturb', - 9)
        end if
      end if
    end do ! loop 200

    return

  end subroutine burntobare
  !! @}
  !> \namespace disturbance_scheme
  !! Calculates disturbance as both natural and human-influenced fires.
  !!
  !! CTEM represents disturbance as both natural and human-influenced fires.
  !! The original fire parametrization corresponding to CTEM v. 1.0 is described
  !! in Arora and Boer (2005) \cite Arora20052ac. The parametrization has since been adapted and used
  !! in several other DGVMs (Kloster et al. 2010 \cite Kloster2010-633, Kloster et al. 2012 \cite Kloster2012-c79,
  !! Migliavacca et al. 2013 \cite Migliavacca2013-eh, and Li et al. 2012 \cite Li20121c2).
  !!  CTEM v. 2.0 incorporates changes suggested in these studies as well as several new improvements.
  !!
  !! Fire in CTEM is simulated using a process-based scheme of intermediate complexity
  !! that accounts for all elements of the fire triangle: fuel load, combustibility
  !! of fuel, and an ignition source. CTEM represents the probability of a fire
  !! occurrence (\f$P_\mathrm{f}\f$), for a representative area of \f$500\, km^2\f$
  !! (\f$a_{rep}\f$), as
  !!
  !! \f[ \label{fieya} P_\mathrm{f} = P_\mathrm{b}P_\mathrm{i}P_\mathrm{m}, \qquad (Eqn 1)\f]
  !!
  !! where the right hand side terms represent the fire probabilities that are
  !! conditioned on (i) the availability of biomass as a fuel source
  !! (\f$P_\mathrm{b}\f$), (ii) the combustibility of the fuel based on its moisture
  !! content (\f$P_\mathrm{m}\f$), and (iii) the presence of an ignition source
  !! (\f$P_\mathrm{i}\f$). The probability of fire and the subsequent calculations
  !! are performed for each PFT present in a grid cell (but the PFT index
  !! \f$\alpha\f$ is omitted for clarity in Eqn 1). Since the CTEM
  !! parametrization is based on one fire per day per representative area, the
  !! representative area has to be sufficiently small that the requirement of only
  !! one fire per day is reasonable, yet sufficiently large such that it is not
  !! possible to burn the entire representative area in 1 day. Based on MODIS observed
  !! fire counts in Fig. 1 of Li et al. (2012) \cite Li20121c2, \f$500\, km^2\f$ is an appropriate size
  !! to not have more than one fire per day and still be a large enough area to be
  !! assumed representative of the grid cell as a whole.
  !!
  !! The \f$P_\mathrm{b}\f$ term depends on the aboveground biomass (\f$B_{ag}\f$)
  !! available for sustaining a fire (which includes the green and brown leaf mass,
  !! stem mass and litter mass, \f$B_{ag} = C_\mathrm{L} + C_\mathrm{S} + C_\mathrm{D}\f$).
  !! Below a lower threshold of aboveground biomass (\f$B_{low}\f$; \f$0.2\, kg\, C\, m^{-2}\f$
  !! similar to Moorcroft et al. (2001) \cite Moorcroft2001-co, and Kucharik et al. (2000) \cite Kucharik2000-xk), fire is not
  !! sustained and thus has a probability of 0. Above a biomass of
  !! \f$1.0\, kg\, C\, m^{-2}\f$ (\f$B_{high}\f$), \f$P_\mathrm{b}\f$ is set to 1 as the
  !! amount of fuel available is assumed sufficient for fire. \f$P_\mathrm{b}\f$ is
  !! then calculated using the aboveground biomass, \f$B_{ag}\f$ (\f$kg\, C\, m^{-2}\f$)
  !! with a linear variation between the upper and lower thresholds as
  !!
  !! \f[ \label{eqn:Pb} P_\mathrm{b}=\max\left[0, \min\left(1, \frac{B_{ag}-B_{low}}
  !! {B_{high} - B_{low}}\right)\right]. \qquad (Eqn 2)\f]
  !!
  !! The linear decrease of \f$P_\mathrm{b}\f$ from \f$B_{high}\f$ to \f$B_{low}\f$
  !! reflects the fragmentation of fuel that occurs as biomass decreases. Fuel
  !! fragmentation impacts upon area burned as it impedes the fire spread rate (Guyette et al. 2002)
  !! \cite Guyette2002-rc.
  !!
  !! The probability of fire based on the presence of ignition sources (\f$P_\mathrm{i}\f$)
  !! is influenced by both natural (lightning) and anthropogenic agents (either
  !! intentional or accidental). An initial lightning scalar, \f$\vartheta_F\f$,
  !! that varies between 0 and 1 is found as
  !!
  !! \f[ \vartheta_F = \max\left[0, \min \left(1, \frac{F_c2g - F_{low}}{F_{high}
  !! - F_{low}} \right)\right], \qquad (Eqn 3)\f]
  !!
  !! where \f$F_{low}\f$ and \f$F_{high}\f$ represent lower and upper thresholds of
  !! cloud-to-ground lightning strikes (\f$F_c2g\f$, \f$flashes\, km^{-2}\, month^{-1}\f$)
  !! , respectively. Similar to Eqn 2, below the lower threshold
  !! (\f$F_{low}\f$; \f$0.25\, flashes\, km^{-2}\, month^{-1}\f$), \f$\vartheta_F\f$
  !! is 0 implying lightning strikes are not sufficient to cause fire ignition,
  !! above the upper threshold (\f$F_{high}\f$; \f$10.0\, flashes\, km^{-2}\, month^{-1}\f$)
  !! \f$\vartheta_F\f$ is 1, as ignition sources now do not pose a constraint on fire.
  !! The amount of cloud-to-ground lightning, \f$F_c2g\f$, is a fraction of the total
  !! lightning based on the relationship derived by Price and Rind (1993) \cite Price1993-fm (approximation
  !! of their Eqs. 1 and 2) as
  !!
  !! \f[ F_c2g = 0.22 \exp (0.0059 \times \vert {\Phi}\vert) F_{tot}, \qquad (Eqn 4)\f]
  !!
  !! where \f$\Phi\f$ is the grid cell latitude in degrees and \f$F_{tot}\f$ is the
  !! total number of lightning \f$flashes\, km^{-2}\, month^{-1}\f$ (both cloud-to-cloud
  !! and cloud-to-ground). The probability of fire due to natural ignition,
  !! \f$P_i, n\f$, depends on the lightning scalar, \f$\vartheta_F\f$, as
  !!
  !! \f[ P_i, n = y(\vartheta_F) - y(0)(1 -  \vartheta_F) + \vartheta_F[1-y(1)]
  !! \nonumber\\ y(\vartheta_F) = \frac{1}{1 + \exp\left(\frac{0.8 - \vartheta_F}{0.1}\right)}. \qquad (Eqn 5)\f]
  !!
  !! Fire probability due to ignition caused by humans, \f$P_i, h\f$, is parametrized
  !! following Kloster et al. (2010) \cite Kloster2010-633 with a dependence on population density,
  !! \f$p_\mathrm{d}\f$ (\f$number of people\, km^{-2}\f$)
  !!
  !! \f[ \label{eqn:Ph} P_i, h = \min\left[1, \left(\frac{p_\mathrm{d}}{p_{thres}}\right)^{0.43}\right], \qquad (Eqn 6)\f]
  !!
  !! where \f$p_{thres}\f$ is a population threshold (\f$300\, people\, km^{-2}\f$) above which \f$P_{i, h}\f$
  !! is 1. The probability of fire conditioned on ignition, \f$P_\mathrm{i}\f$, is then
  !! the total contribution from both natural and human ignition sources
  !!
  !! \f[ \label{eqn:Pi} P_\mathrm{i} = \max[0, \min\{1, P_{i, n} + (1 - P_{i, n})P_{i, h}\}].\qquad (Eqn 7) \f]
  !!
  !! The population data used to calculate probability of fire ignition caused by humans
  !! and anthropogenic fire suppression (discussed further down in this section) is
  !! read in as a model geophysical field.
  !!
  !! The probability of fire due to the combustibility of the fuel, \f$P_\mathrm{m}\f$,
  !! is dependent on the soil moisture in vegetation's root zone and in the litter
  !! layer. The root-zone soil wetness (\f$\phi_{root}\f$, allocate.f90 Eqn. 1)
  !! is used as a surrogate for the vegetation moisture content and the soil wetness
  !! of the top soil layer as a surrogate for the litter moisture content. If a grid
  !! cell is covered by snow, \f$P_\mathrm{m}\f$ is set to zero. The probability of
  !! fire conditioned on soil wetness in vegetation's rooting zone, \f$P_{m, V}\f$,
  !! is then
  !!
  !! \f[ P_{m, V} = 1-\tanh
  !! \left[\left( \frac{1.75\ \phi_{root}} {E_\mathrm{V}}\right)^2\right], \qquad (Eqn 8)\f]
  !!
  !! where \f$E_\mathrm{V}\f$ is the extinction soil wetness above which \f$P_{f, V}\f$
  !! is reduced to near zero and is set to 0.30.
  !!
  !! The probability of fire based on the moisture content in the \f$\textit{duff}\f$
  !! layer, \f$P_{m, D}\f$, which includes the brown leaf mass (grasses only) and
  !! litter mass (\f$B_{duff} = C_{L, b} + C_\mathrm{D}\f$; \f$kg\, C\, m^{-2}\f$), is
  !! calculated in a similar way but uses the soil wetness of the first soil layer,
  !! (\f$\phi_1\f$, photosynCanopyConduct.f90 Eqn. 7), as a surrogate for the moisture in the duff
  !! layer itself as
  !!
  !! \f[ P_{m, D} = 1 -\tanh\left[\left(\frac{1.75 \phi_1}{E_{\mathrm{D}}}\right)^2\right], \qquad (Eqn 9)\f]
  !!
  !! where the extinction soil wetness for the litter layer, \f$E_{\mathrm{D}}\f$,
  !! is set to 0.50, which yields a higher probability of fire for the litter layer
  !! than for the vegetation for the same soil wetness. \f$P_\mathrm{m}\f$ is then
  !! the weighted average of \f$P_{m, V}\f$ and \f$P_{m, D}\f$ given by
  !!
  !! \f[ \label{eqn:Pf} P_\mathrm{m} = P_{m, V} (1-f_{duff}) + P_{m, D} f_{duff}
  !! \nonumber \\ f_{duff}=\frac{B_{duff}}{B_{ag}}\qquad (Eqn 10)\f]
  !!
  !! where \f$f_{duff}\f$ is the duff fraction of aboveground combustible biomass.
  !!
  !! The area burned (\f$a\f$) is assumed to be elliptical in shape for fires based
  !! upon the wind speed and properties of an ellipse
  !!
  !! \f[ a(t)=\pi \frac{l}{2}\frac{w}{2}= \frac{\pi}{2} (v_\mathrm{d}+v_\mathrm{u})v_\mathrm{p}t^2, \qquad (Eqn 11)\f]
  !!
  !! where \f$l\f$ (\f$m\f$) and \f$w\f$ (\f$m\f$) are the lengths of major and minor
  !! axes of the elliptical area burnt; \f$v_\mathrm{d}\f$ (\f$km\, h^{-1}\f$) and
  !! \f$v_\mathrm{u}\f$ (\f$km\, h^{-1}\f$) are the fire spread rates in the downwind
  !! and upwind directions, respectively; \f$v_\mathrm{p}\f$ (\f$km\, h^{-1}\f$) is
  !! the fire spread rate perpendicular to the wind direction and \f$t\f$ is the
  !! time (\f$h\f$).
  !!
  !! The fire spread rate in the downwind direction (\f$v_\mathrm{d}\f$) is represented
  !! as
  !!
  !! \f[ \label{firespreadrate} v_\mathrm{d} = v_{d, max}\, g(u)\, h(\phi_{r, d})\qquad (Eqn 12)\f]
  !!
  !! where \f$v_{d, max}\f$ (\f$km\, h^{-1}\f$) is the PFT-specific maximum fire spread
  !! rate from Li et al. (2012)\cite Li20121c2, which is set to zero for crop PFTs (see also
  !! classicParams.f90). The functions \f$g(u)\f$ accounts for the effect of wind
  !! speed and \f$ h(\phi_{r, d})\f$ accounts for the effect of rooting zone and
  !!  duff soil wetness on the fire spread rate, as discussed below.
  !!
  !! The wind speed (\f$u\f$; \f$km\, h^{-1}\f$) is used to determine the length
  !! (\f$l\f$) to breadth (\f$w\f$) ratio, \f$L_\mathrm{b}\f$, of the elliptical area
  !! burned by fire
  !!
  !! \f[ \label{lb} L_\mathrm{b}= \frac{l}{w} = \frac{v_\mathrm{d} + v_\mathrm{u}}{2v_\mathrm{p}}
  !! = 1 + 10 [1 -\exp(-0.06 u)] \qquad (Eqn 13)\f]
  !!
  !! and its head to back ratio, \f$H_\mathrm{b}\f$, following Li et al. (2012) \cite Li20121c2, as
  !!
  !! \f[ \label{hb} H_\mathrm{b} = \frac{v_\mathrm{d}}{v_\mathrm{u}} =
  !! \frac{L_\mathrm{b} + (L_\mathrm{b}^2 - 1)^{0.5}}{L_\mathrm{b} -
  !! (L_\mathrm{b}^2 - 1)^{0.5}},
  !! \qquad (Eqn 14)\f]
  !!
  !! which help determine the fire spread rate in the direction perpendicular to
  !! the wind speed and in the downward direction. Equations 13 and
  !! 14 are combined to estimate the wind scalar \f$g(u)\f$ as
  !!
  !! \f[ g(u)= g(0) \frac{2.0 L_\mathrm{b}}{(1 + 1/H_\mathrm{b})} \nonumber\\
  !! \frac{g(u)}{g(0)}=\frac{v_\mathrm{d}}{v_\mathrm{p}} = \frac{2.0 L_\mathrm{b}}
  !! {(1 + 1/H_\mathrm{b})}, \qquad (Eqn 15)
  !! \f]
  !!
  !! which varies between 0.05 and 1. The lower limit is imposed by the \f$g(0)\f$ term,
  !! which has a value of 0.05 and represents the fire spread rate in the absence of
  !! wind (\f$u = 0\f$); the upper limit is assigned a maximum value of 1. The fire
  !! spread rate in the absence of wind is essentially the spread rate in the
  !! direction perpendicular to the wind speed (\f$v_\mathrm{p}\f$). The value of the
  !! \f$g(0)\f$ term is derived by considering the case where the wind speed becomes
  !! very large. As \f$u\f$ \f$\rightarrow \infty\f$ then \f$L_\mathrm{b}
  !! \rightarrow 11\f$ and \f$H_\mathrm{b} \rightarrow 482\f$, while
  !! \f$g(\infty)=1\f$ due to its definition, which yields \f$g(0) = 0.0455
  !! \approx 0.05\f$.
  !!
  !! The dependence of fire spread rate on the rooting zone and duff soil wetness,
  !! \f$h(\phi_{r, d})\f$ is represented as
  !!
  !! \f[ h(\phi_{r, d})= h(\phi_{root})(1-f_{duff}) + h(\phi_{1})f_{duff}\nonumber
  !! \\ h(\phi_{root})= \left(1-min \left(1, \frac{\phi_{root}}{E_\mathrm{V}} \right)
  !! \right)^2\nonumber \\ h(\phi_{1})= \left(1-min \left(1, \frac{\phi_{1}}{E_\mathrm{D}} \right) \right)^2.
  !! \qquad (Eqn 16)\f]
  !!
  !! Both \f$h(\phi_{root})\f$ and \f$h(\phi_{1})\f$ gradually decrease from 1
  !! (when soil wetness is 0 and soil moisture does not constrain fire spread rate)
  !! to 0 when soil wetness exceeds the respective extinction wetness thresholds,
  !! \f$E_\mathrm{V}\f$ and \f$E_\mathrm{D}\f$.
  !!
  !! With fire spread rate determined, and the geometry of the burned area defined,
  !! the area burned in 1 day, \f$a_{1{\mathrm{d}}}\f$ (\f$km^2\, day^{-1}\f$),
  !! following Li et al. (2012) \cite Li20121c2, is calculated as
  !!
  !! \f[ a_{1{\mathrm{d}}} = \frac{\pi v_\mathrm{d}^2 t^2}{4L_\mathrm{b}}\left(1 +
  !! \frac{1}{H_\mathrm{b}}\right)^2 \nonumber\\ = \frac{\pi v_\mathrm{d}^2 (24^2)}
  !! {4L_\mathrm{b}}\left(1 + \frac{1}{H_\mathrm{b}}\right)^2\label{aburned}
  !! \qquad (Eqn 17)\f]
  !!
  !! by setting \f$t\f$ equal to \f$24\, h\f$.
  !!
  !! The fire extinguishing probability, \f$q\f$, is used to calculate the duration
  !! (\f$\tau\f$, \f$days\f$) of the fire, which in turn is used to calculated the
  !! area burned over the duration of the fire, \f${a_{\tau d}}\f$. \f$q\f$ is
  !! represented following Kloster et al. (2010) \cite Kloster2010-633 as
  !!
  !! \f[ q = 0.5 + \frac{\max\left[0, 0.9 - \exp(-0.025\, p_\mathrm{d})\right]}{2}, \qquad (Eqn 18)
  !! \f]
  !!
  !! which yields a value of \f$q\f$ that varies from 0.5 to 0.95 as population density,
  !! \f$p_\mathrm{d}\f$ (\f$\frac{\text{number of people}}{km^{2}}\f$), increases from zero to
  !! infinity. Higher population density thus implies a higher probability of fire
  !! being extinguished. \f$q\f$ represents the probability that a fire will be
  !! extinguished on the same day it initiated and the probability that it will
  !! continue to the next day is (\f$1-q\f$). Assuming individual days are independent,
  !! the probability that the fire will still be burning on day \f$\tau\f$ is
  !! \f$(1-q)^\tau\f$. The probability that a fire will last exactly \f$\tau\f$ days,
  !! \f$P(\tau)\f$, is the product of the probability that the fire still exists at day
  !! \f$\tau\f$ and the probability it will be extinguished on that day hence
  !! \f$P(\tau) = q(1-q)^\tau\f$. This yields an exponential distribution of fire
  !! duration whose expected value is
  !!
  !! \f[ \overline{\tau} = E(\tau) = \sum_{\tau=0}^\infty\, \tau\, q(1-q)^{\tau}
  !! = \frac{1-q}{q}.\qquad (Eqn 19)
  !! \f]
  !!
  !! Based on this fire duration and the area burned in 1 day (Eqn. 17),
  !! the area burned over the duration of the fire (\f$a_{\tau \mathrm{d}}\f$)
  !! (but still implemented in 1 day since the model does not track individual fires
  !! over their duration, \f$km^2\, day^{-1}\f$) is calculated as
  !!
  !! \f[ a_{\tau \mathrm{d}} =E(a_{1{\mathrm{d}}} \tau^2)=\sum_{\tau=0}^
  !! \infty\, a_{1{\mathrm{d}}}\, \tau^2  q(1-q)^{\tau} \\ = a_{1{\mathrm{d}}}\,
  !! \frac{(1-q) (2-q)}{q^2}.\nonumber\qquad (Eqn 20)
  !! \f]
  !!
  !! Finally, and reintroducing the PFT index \f$\alpha\f$, the area burned
  !! is extrapolated for a PFT \f$\alpha\f$ (\f$A_{\mathrm{b}, \alpha}\f$, \f$km^2\,
  !! day^{-1}\f$) to the whole grid cell as
  !!
  !! \f[A_{\mathrm{b}, \alpha}=P_{f, \alpha}\, a_{\tau \mathrm{d}, \alpha}
  !! \frac{A_\mathrm{g}f_\alpha}{a_{rep}}, \qquad (Eqn 21) \f]
  !!
  !! where \f$A_\mathrm{g}\f$ is area of a grid cell (\f$km^2\f$), \f$f_\alpha\f$
  !! the fractional coverage of PFT \f$\alpha\f$ and \f$a_{rep}\f$ the representative
  !! area of \f$500\, km^2\f$, as mentioned earlier. Area burned over the whole grid
  !! cell (\f$A_\mathrm{b}\f$, \f$km^2\, day^{-1}\f$) is then calculated as the sum
  !! of area burned for individual PFTs,
  !!
  !! \f[ A_\mathrm{b}=\sum_{\alpha=1}^{N}A_{\mathrm{b}, \alpha}.\qquad (Eqn 22)\f]
  !!
  !! Fire emits \f$CO_2\f$, other trace gases, and aerosols as biomass is burned while
  !! plant mortality and damage due to fire contribute to the litter pool. The
  !! emissions of a trace gas/aerosol species \f$j\f$ from PFT \f$\alpha\f$,
  !! \f$E_{\alpha, j}\f$ (\f$g species (m^{-2} grid cell area) day^{-1}\f$) are
  !! obtained from a vector of carbon densities \f$\vec{C}_{\alpha} = (C_\mathrm{L},
  !! C_\mathrm{S}, C_\mathrm{R}, C_\mathrm{D})_\alpha\f$ (\f$kg\, C\, m^{-2}\f$) for
  !! its leaf, stem, root and litter components, multiplied by a vector of combustion
  !! factors \f$mho_{\alpha} = (mho_\mathrm{L}, mho_\mathrm{S}, mho_\mathrm{R},
  !! mho_\mathrm{D})_\alpha\f$, which determines what fraction of leaf, stem, root
  !! and litter components gets burned, multiplied by a vector of emissions factors
  !! \f$\Upsilon_{j} = (\Upsilon_\mathrm{L}, \Upsilon_\mathrm{S}, \Upsilon_\mathrm{R},
  !! \Upsilon_\mathrm{D})_j\f$ (\f$g species (kg\, C\, dry organic matter)^{-1}\f$), and
  !! by the area burned \f$A_{\mathrm{b}, \alpha}\f$ for that PFT.
  !!
  !! The dot product of \f$\vec{C}_{\alpha}\f$, \f$\Upsilon_{j}\f$ and \f$mho_{\alpha}
  !! \f$ thus yields emissions per unit grid cell area of species \f$j\f$ from PFT
  !! \f$\alpha\f$,
  !!
  !! \f[ \label{emiss_combust_factor} {E_{\alpha, j}}= ((\vec{C}_\alpha\cdot mho_{\alpha}
  !! )\cdot \Upsilon_{j}) \frac{A_{\mathrm{b}, \alpha}}{A_\mathrm{g}}\frac{1000}{450}, \qquad (Eqn 23)
  !! \f]
  !!
  !! where the constant 1000 converts \f$\vec{C}_\alpha\f$ from \f$kg\, C\, m^{-2}\f$
  !! to \f$g\, C\, m^{-2}\f$ and the constant 450 (\f$g\, C\, (kg dry organic matter)^{-1}
  !! \f$) converts biomass from carbon units to dry organic matter (Li et al. 2012)\cite Li20121c2.
  !! The corresponding loss of carbon (\f$kg\, C\, m^{-2}\, day^{-1}\f$) from the three
  !! live vegetation components (L, S, R) and the litter pool (D) of PFT \f$\alpha\f$
  !! is given by
  !!
  !! \f[ \label{emiss_combust_loss} H_{\alpha, i}= C_{\alpha, i}
  !! mho_i\left(\frac{A_{\mathrm{b}, \alpha}}{A_\mathrm{g}}\right)\quad i={L, S, R, D}.\qquad (Eqn 24)
  !! \f]
  !!
  !! The PFT-specific combustion factors for leaf (\f$mho_\mathrm{L}\f$), stem
  !! (\f$mho_{\mathrm{S}}\f$), root (\f$mho_{\mathrm{R}}\f$) and litter
  !! (\f$mho_{\mathrm{D}}\f$) components are summarized in classicParams.f90.
  !! Emission factors for all species of trace gases and aerosols
  !! (\f$CO_2\f$, \f$CO\f$, \f$CH_4\f$, \f$H_2\f$, \f$NHMC\f$, \f$NO_x\f$,
  !! \f$N_2O\f$, total particulate matter, particulate matter less than \f$2.5\,
  !! \mu m\f$ in diameter, and black and organic carbon) are read in from the CLASSIC parameters namelist file.
  !!
  !! Litter generated by fire is based on similar mortality factors, which reflect
  !! a PFT's susceptibility to damage due to fire \f$\vec{\Theta}_{\alpha} =
  !! (\Theta_\mathrm{L}, \Theta_\mathrm{S}, \Theta_\mathrm{R})_\alpha\f$ (fraction).
  !! The contribution to litter pool of each PFT due to plant mortality associated
  !! with fire (\f$kg\, C\, m^{-2}\, day^{-1}\f$) is calculated as
  !!
  !! \f[ \label{eqn_using_mort_factors} {M_{\alpha}}= (\vec{C}_\alpha \cdot
  !! \Theta_{\alpha} ) \frac{A_{\mathrm{b}, \alpha}}{A_\mathrm{g}}, \qquad (Eqn 25)
  !! \f]
  !!
  !! which is the sum of contribution from individual live vegetation pools
  !!
  !! \f[ \label{eqn_using_mort_factors_individual} M_{\alpha, i}= C_{\alpha, i}
  !! \Theta_{\alpha, i} \left(\frac{A_{\mathrm{b}, \alpha}}{A_\mathrm{g}} \right)
  !! \quad i={L, S, R}.\qquad (Eqn 26)
  !! \f]
  !!
  !! The carbon loss terms associated with combustion of vegetation components and
  !! litter (\f$H_{\alpha, i}, i={L, S, R, D}\f$) and mortality of vegetation components
  !! (\f$M_{\alpha, i}, i={L, S, R}\f$) due to fire are used in \ref CTEMRateChgEqns Eqns 1 and 3, which
  !! describe the rate of change of carbon in model's five pools (however, listed
  !! there without the PFT subscript \f$\alpha\f$). The PFT-specific mortality factors
  !! for leaf (\f$\Theta_\mathrm{L}\f$), stem (\f$\Theta_{\mathrm{S}}\f$) and root
  !! (\f$\Theta_\mathrm{R}\f$) components are listed in classicParams.f90.
  !!
  !! When CTEM is run with prescribed PFT fractional cover, the area of PFTs does not
  !! change and the fire-related emissions of \f$CO_2\f$, other trace gases and aerosols
  !! , and generation of litter act to thin the remaining biomass. When competition
  !! between PFTs for space is allowed, fire both thins the remaining biomass and
  !! through plant mortality creates bare ground, which is subsequently available for
  !! colonization. The creation of bare ground depends on the susceptibility of each
  !! PFT to stand replacing fire (\f$\zeta_\mathrm{r}\f$, fraction) (see also
  !! classicParams.f90) and the PFT area burned. The fire-related mortality rate,
  !! \f$m_{dist}\f$ (\f$day^{-1}\f$), used in mortality.f90 Eqn. 1, is then
  !!
  !! \f[ \label{m_dist} m_{dist, \alpha} = \zeta_{\mathrm{r}, \alpha}
  !! \frac{A_{\mathrm{b}, \alpha}}{f_\alpha A_\mathrm{g}}.\qquad (Eqn 27)
  !! \f]
  !!
  !! After bare ground generation associated with fire, the thinned biomass
  !! is spread uniformly over the remaining fraction of a PFT. However,
  !! it is ensured that the carbon density of the remaining biomass does not
  !! increase to a value above what it was before the fire occurred.
  !!

  !> \file
end module disturbance_scheme
