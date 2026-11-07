!> \file
!> Methane subroutines for wetland areal determination and methane emission and upland (oxic) soils methane uptake
module methaneProcesses

  implicit none

  ! Subroutines
  public :: wetland_methane
  public :: soil_ch4uptake

contains

  ! ------------------------------------------------------------------
  !> \ingroup methaneprocesses_soil_ch4uptake
  !! @{
  !> Calculates the uptake of atmospheric methane by oxic soils.
  !> @author J. Melton
  subroutine soil_ch4uptake (il1, il2, ilmos, ilg, tbar, & ! In
                             bi, thlq, thic, psis, & ! In
                             fcan, wetfdyn, wetfrac, & ! In
                             isand, atm_CH4, thp, & ! In
                             CH4_soills) ! Out

    ! History:

    !  06  Dec 2018  - Pass ilg back in as an argument + minor reorganization
    !  V. Arora        of arguments

    !  J. Melton. Dec 22 2015 - Coded up based on C. Curry (2007) Modelling the
    !  the soil consumption of atmospheric methane at the global scale. Global
    !  Biogeo. Cycl. v. 21 GB4012 doi: 10.1029/2006GB002818.

    use classicParams, only : ignd, ican, nlat, wtCH4, D_air, g_0, betaCH4, k_o, &
                              GRAV, RHOW, RHOICE, icp1

    implicit none

    ! Arguments:
    real, dimension(ilg,ignd), intent(in) :: tbar     !< Temperature of soil layers (K) - daily average
    real, dimension(ilg,ignd), intent(in) :: THP      !< Total porosity \f$(cm^3 cm^{-3})\f$ - daily average
    real, dimension(ilg,ignd), intent(in) :: BI       !< Clapp and Hornberger b-term (-)
    real, dimension(ilg,ignd), intent(in) :: THLQ     !< Fractional water content (-) - daily average
    real, dimension(ilg,ignd), intent(in) :: THIC     !< Fractional ice content (-) - daily average
    real, dimension(ilg,ignd), intent(in) :: PSIS     !< Soil moisture suction at saturation (m)
    real, dimension(ilg,icp1), intent(in) :: FCAN     !< Fractional coverage of vegetation (-)
    real, dimension(nlat), intent(in) :: wetfrac      !< Prescribed fraction of wetlands in a grid cell
    real, dimension(ilg), intent(in) :: wetfdyn       !< Dynamic gridcell wetland fraction determined using slope and soil moisture
    real, dimension(ilg), intent(in) :: atm_CH4       !< Atmospheric \f$CH_4\f$ concentration at the soil surface (ppmv)
    integer, dimension(ilg,ignd), intent(in) :: isand !< flag for soil/bedrock/ice/glacier
    integer, intent(in), dimension(:) :: ilmos        !< Index of gridcell corresponding to current element of gathered vector of land surface variables [ ]
    integer, intent(in) :: il1
    integer, intent(in) :: il2
    integer, intent(in) :: ilg
    real, dimension(ilg), intent(out) :: CH4_soills   !< Methane uptake into the soil column \f$(\mu mol CH4 m^{-2} s^{-1})\f$

    ! Local variables:
    real :: Tsoil                           !< Temperature of soil layers \f$(\circ C)\f$
    real :: D_soil                          !< Diffusivity of CH4 in soil \f$(cm^2 s^{-1})\f$
    real :: G_T                             !< Temperature factor used in determining D_soil (-)
    real :: G_soil                          !< Soil moisture factor used in determining D_soil (-)
    real :: THP_air                         !< Air-filled porosity \f$(cm^3 cm^{-3})\f$
    real :: k_oxidr                         !< First-order oxidation rate constant \f$(s^-1)\f$
    real :: r_T                             !< Temperature factor used in determination of rate constant (-)
    real :: r_SM                            !< Soil moisture factor used in determination of rate constant (-)
    integer :: i, j, layer                  !< Counters
    real :: psi                             !< Soil moisture suction / matric potential (m)
    real :: r_C                             !< Factor to account for croplands
    real :: r_W                             !< Factor to account for wetlands
    real :: THP_tot                         !< temp variable for total porosity \f$(cm^3 cm^{-3})\f$

    !> ---------------------------------------------------------------------
    !> Begin
    !!
    !! The soil oxidation methane sink is assumed to only operate in the first model
    !! soil layer, thus we only consider that layer here.
    layer = 1

    do i = IL1,IL2

      if (isand(i,layer) <= - 1) cycle !> not soil so move on.

      !> Convert tbar to Tsoil (from K to deg C)
      Tsoil = tbar(i,layer) - 273.16

      !> Find the diffusion coefficient in soil (D_soil)

      !> First the temperature factor, G_T:
      G_T = 1.0 + 0.0055 * Tsoil

      !> Find the air filled porosity,THP_air:
      THP_air = THP(i,layer) - (THLQ(i,layer) + THIC(i,layer) * RHOICE/RHOW)
      THP_tot = THP(i,layer)

      !> Note: THP_air can fall to < 0 after snow melt
      if (THP_air  < 0.) then
        THP_air = 0.0
        THP_tot = (THLQ(i,layer) + THIC(i,layer) * RHOICE/RHOW)
      end if

      !> The BI  (Clapp and Hornberger b-term) is already calculated by CLASS as:
      !> BI = 15.9 * f_clay + 2.91, thus we use that value.

      !> G_soil is the influence of the soil texture, moisture, and porosity:
      G_soil = THP_tot ** (4./3.) * (THP_air / THP_tot) ** (1.5 + 3. / BI(i,layer))

      !> The diffusion coefficient of CH4 in soil is then:
      D_soil = D_air * G_T * G_soil

      !> Determine the first-order oxidation rate constant (k_oxidr)

      !> First find the temperature term, r_T (FLAG note that Charles' original code does not have the high temp limit !)

      if (Tsoil < 0.0 .and. Tsoil >= - 10.0) then
        r_T = (0.1 * Tsoil + 1.0) ** 2
      else if (Tsoil >= 0.0 .and. Tsoil < 43.3) then
        r_T = exp(0.0693 * Tsoil - 8.56E-7 * Tsoil ** 4)
      else !> all other temps (<-10 and >=43.3)
        r_T = 0.
      end if

      !> Next find the term based on soil moisture (suction)

      !> Find the soil water potential for the uppermost layer
      !> need the absolute value.
      psi = abs(PSIS(i,layer) * (THLQ(i,layer)/THP_tot) ** ( - BI(i,layer)))

      !> Convert units from m to kPa
      psi = psi * GRAV

      if (psi < 200.) then !> 0.2 MPa in paper (NOTE: In Charles's code this is \f$\leq\f$, but is < in paper)
        r_SM = 1.0
      else if (psi >= 200. .and. psi <= 1.E5) then !> 0.2 and 100 Mpa in paper
        r_SM = (1. - (log10(psi) - log10(200.)) / (log10(1.E5) - log10(200.))) ** betaCH4
      else !> psi > 100 MPa.
        r_SM = 0.
      end if

      k_oxidr = k_o * r_T * r_SM

      !> Find the flux correction for croplands

      r_C = 1.0 - (0.75 * FCAN(i,3))

      !> Find the flux correction due to wetlands
      ! Presently we assume the wetland fraction is split across 

      if (wetfrac(ilmos(i)) > 0.) then  !> Use the prescribed wetland fractions
        r_W = 1.0 - wetfrac(ilmos(i))
      else !> use the dynamically determined wetland area
        r_W = 1.0 - wetfdyn(ilmos(i))
      end if

      !> Find the surface flux (CH4_soills) for each tile, then for each gridcell

      CH4_soills(i) =  atm_CH4(i) * r_C * r_W * g_0 * sqrt(D_soil * k_oxidr)

      !> Convert from mg CH4 m^-2 s^-1 to umol CH4 m^-2 s^-1

      CH4_soills(i) = CH4_soills(i) * 1.E3 / wtCH4


    end do ! loop 10

  end subroutine soil_ch4uptake
  !> @}
  ! ------------------------------------------------------------------

  !> \ingroup methaneprocesses_wetland_methane
  !! @{
  !> Prognostic determination of wetland fractional coverage and emission of wetland methane subroutine
  !> @author B. Amiro, J. Melton, V. Arora
  subroutine wetland_methane (hetrores, il1, il2, ilg, wetfrac, &       ! inputs
                              thliqg, currlat, sand, slopefrac, ta, &   ! inputs
                              ch4WetSpec, wetfdyn, ch4WetDyn)         ! outputs

    !     20  Feb   2019 - Update code to allow namelist of params. Remove no CH4 flux if soils are
    !     J. Melton        frozen condition.
    !     31  Aug   2016 - Change how we find wetlands from discrete limits to
    !     V. Arora         smooth function
    !
    !      4  July  2014 - Convert to f90 and bring in classicParams
    !     J. Melton
    !
    !      9  June. 2010 - this subroutine calculates methane flux
    !     Brian Amiro      from wetlands of a grid cell (i.e. land only)

    use classicParams,    only : wtdryres, ratioch4, lat_thrshld1, &
                                 lat_thrshld2, soilw_thrshN, soilw_thrshE, &
                                 soilw_thrshS, ignd

    implicit none

    integer, intent(in) :: ilg                      !<
    integer, intent(in) :: il1                      !< il1=1
    integer, intent(in) :: il2                      !< il2=ilg
    real, dimension(ilg), intent(in) :: hetrores    !< heterotrophic respiration from main ctem program calculated as sum of litres + socres
    real, dimension(ilg), intent(in) :: ta          !< air temperature, k
    real, dimension(ilg), intent(in) :: wetfrac     !< prescribed fraction of wetlands in a grid cell
    real, dimension(ilg), intent(in) :: currlat     !< centre latitude of grid cells in degrees
    real, dimension(ilg,8), intent(in) :: slopefrac !< Fraction of gridcell flatter than slope thresholds (0.025, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3 and 0.35 percent slope thresholds)
    real, dimension(ilg,ignd), intent(in) :: thliqg !< liquid soil moisture content (fraction)
    real, dimension(ilg,ignd), intent(in) :: sand   !< percentage sand in soil layers
    real, dimension(ilg), intent(out) :: ch4WetSpec    !< methane flux from wetlands calculated using hetrores in umol ch4/m2.s
    real, dimension(ilg), intent(out) :: wetfdyn    !< dynamic gridcell wetland fraction determined using  slope and soil moisture
    real, dimension(ilg), intent(out) :: ch4WetDyn    !< methane flux from wetlands calculated using hetrores and wetfdyn, in umol ch4/m2.s

    ! local variables
    real, dimension(ilg) :: wetresp !< heterotrophic wetland respiration
    real :: upvslow    ! ratio of wetland to upland respiration for this location
    real :: porosity
    real :: soil_wetness
    integer :: i

    real, parameter :: alpha = 0.45  ! Determines shape of curve.
    real :: low_mois_lim
    real :: mid_mois_lim
    real :: upp_mois_lim
    real :: x1
    real :: x2
    real :: y1
    real :: y2
    real :: slope
    real :: intercept

    !>
    !> ---------------------------------------------------------------
    !> Constants and parameters are located in classicParams.f90
    !> -----------------------------------------------------------------
    !>
    !> initialize required arrays to zero
    !>
    wetresp(:) = 0.0
    wetfdyn(:) = 0.0
    ch4WetSpec(:) = 0.0
    ch4WetDyn(:) = 0.0
    !>
    !> --------------------------------------------------
    !>
    !> Estimate the methane flux from wetlands for each grid cell scaling by the wetland fraction in a grid cell

    !> Set up the latitude bounds based on the parameters read in from the namelist file.
    !! If soil wetness meets a latitude specific threshold then the slope based wetland
    !! fraction is wet and is an actual wetland, else not. As well the ratio of upland to
    !! wetland respiration is assumed to be latitidionally varying. This is intended to
    !! reflect the differnt wetland types in the tropics (floodplain types) vs. those in
    !! higher latitudes (peatlands). This is a very coarse approximation and should be
    !! replaced once we have CH4 in our peatland module. The choice of wtdryres instead
    !! of ratioCH4 to try and mimic this is arbitrary, either parameter could be used since
    !! they are simply multiplicative.

    do i = il1,il2

      if (currlat(i) >= lat_thrshld1) then ! Northern high lats, all area north of lat_thrshld1
        low_mois_lim = soilw_thrshN(1)
        mid_mois_lim = soilw_thrshN(2)
        ! upp_mois_lim = soilw_thrshN(3)
        upvslow = wtdryres(1)
      else if (currlat(i) < lat_thrshld1 .and. currlat(i) >= lat_thrshld2) then ! Tropics
        low_mois_lim = soilw_thrshE(1)
        mid_mois_lim = soilw_thrshE(2)
        ! upp_mois_lim = soilw_thrshE(3)
        upvslow = wtdryres(2)
      else ! S. Hemi, everything else below lat_thrshld2
        low_mois_lim = soilw_thrshS(1)
        mid_mois_lim = soilw_thrshS(2)
        ! upp_mois_lim = soilw_thrshS(3)
        upvslow = wtdryres(3)
      end if

      !> First calculate for the specified wetland fractions read in from OBSWETFFile

      wetresp(i) = hetrores(i) * upvslow * wetfrac(i)
      ch4WetSpec(i) = ratioch4 * wetresp(i)

      !>
      !> Next dynamically find the wetland locations and determine their methane emissions
      !>

      porosity = ( - 0.126 * sand(i,1) + 48.9) / 100.0 ! top soil layer porosity
      soil_wetness = (thliqg(i,1) / porosity)
      soil_wetness = max(0.0,min(soil_wetness,1.0))


      !    implement Vivek's new way of modelling WETFDYN

      x1 = low_mois_lim * (1. - alpha) + mid_mois_lim * alpha
      x2 = 1.
      y1 = 0.
      y2 = slopefrac(i,5)
      slope = (y2 - y1) / (x2 - x1)
      intercept = slope * x1 * ( - 1)

      wetfdyn(i) = min(1.0,max(0.0,slope * soil_wetness + intercept))

      !    new dynamic calculation
      !    same as ch4WetSpec & 2, but wetfrac replaced by wetfdyn

      wetresp(i) = hetrores(i) * upvslow * wetfdyn(i)
      ch4WetDyn(i) = ratioch4 * wetresp(i)

    end do

    return
  end subroutine wetland_methane
  !! @}
  !> \namespace methaneprocesses
  !! Methane subroutines for wetland areal determination and methane emission and upland (oxic) soils methane uptake

  !> \file
end module methaneProcesses
