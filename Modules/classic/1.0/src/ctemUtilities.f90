!> \file
!> Central module for CTEM (biogeochemical)-related utilities
module ctemUtilities

  implicit none

  public :: genSortIndex
  public :: dayEndCTEMPreparation
  public :: accumulateForCTEM
  public :: ctemInit

contains

  ! --------------------------------------------------------------------------------------------------------------------
  !> \ingroup ctemutilities_genSortIndex
  !! @{
  !> Generate the sort index for correspondence between the CTEM pfts and the
  !! array of values in the parameter vectors (e.g. for 9 CTEM, the array is of
  !! size 12)
  !! @author V.Arora, J. Melton
  function genSortIndex ()

    use classicParams, only : ican, l2max, icc, nol2pfts

    implicit none

    integer :: genSortIndex(icc)

    integer :: icount, j, m, n

    icount = 0
    do j = 1,ican
      do m = 1,nol2pfts(j)
        n = (j - 1) * l2max + m
        icount = icount + 1
        genSortIndex(icount) = n
      end do ! loop 96
    end do ! loop 95

  end function genSortIndex
  !! @}

  ! --------------------------------------------------------------------------------------------------------------------
  !> \ingroup ctemutilities_dayEndCTEMPreparation
  !! @{
  !> Prepare the CTEM input (physics) variables at the end of the day.
  !! @author V.Arora, J. Melton
  subroutine dayEndCTEMPreparation (nml, nday, ILMOS)

    use classicParams, only : icc, ignd
    use ctemStateVars, only : vgat, ctem_tile

    implicit none

    integer, intent(in) :: nml     !< Counter representing number of mosaic tiles on modelled domain that are land
    integer, intent(in) :: nday    !< Number of short (physics) timesteps in one day. e.g., if physics timestep is 15 min this is 48.
    integer, intent(in)  :: ILMOS(:) !< Index of grid cell corresponding to current element
    !< of gathered vector of land surface variables [ ]

    integer, pointer :: altotcount_ctm(:) ! nlat  !< Counter used for calculating total albedo
    real, pointer    :: fsinacc_gat(:) !(ilg)    !<
    real, pointer    :: flutacc_gat(:) !(ilg)    !<
    real, pointer    :: flinacc_gat(:) !(ilg)    !<
    real, pointer    :: altotacc_gat(:) !(ilg)   !<
    real, pointer    :: netrad_gat(:) !(ilg)     !<
    real, pointer :: anmossac_t(:)
    real, pointer :: rmlmossac_t(:)
    real, pointer :: gppmossac_t(:)
    real, pointer :: fsnowacc_t(:)       !<
    real, pointer :: taaccgat_t(:)       !<
    real, pointer :: uvaccgat_t(:)       !<
    real, pointer :: vvaccgat_t(:)       !<
    real, pointer :: tbaraccgat_t(:,:) !<
    real, pointer :: thliqacc_t(:,:)  !<
    real, pointer :: thiceacc_t(:,:)  !<
    real, pointer :: ancgvgac_t(:,:)
    real, pointer :: rmlcgvga_t(:,:)
    integer, pointer :: ipeatlandgat(:)

    integer :: i, j
    real :: fsstar_gat
    real :: flstar_gat

    anmossac_t        => ctem_tile%anmossac_t
    rmlmossac_t       => ctem_tile%rmlmossac_t
    gppmossac_t       => ctem_tile%gppmossac_t
    altotcount_ctm    => vgat%altotcount_ctm
    fsinacc_gat       => vgat%fsinacc_gat
    flutacc_gat       => vgat%flutacc_gat
    flinacc_gat       => vgat%flinacc_gat
    altotacc_gat      => vgat%altotacc_gat
    netrad_gat        => vgat%netrad_gat
    ipeatlandgat      => vgat%ipeatland
    tbaraccgat_t      => ctem_tile%tbaraccgat_t
    thliqacc_t        => ctem_tile%thliqacc_t
    thiceacc_t        => ctem_tile%thiceacc_t
    ancgvgac_t        => ctem_tile%ancgvgac_t
    rmlcgvga_t        => ctem_tile%rmlcgvga_t
    fsnowacc_t        => ctem_tile%fsnowacc_t
    taaccgat_t        => ctem_tile%taaccgat_t
    uvaccgat_t        => ctem_tile%uvaccgat_t
    vvaccgat_t        => ctem_tile%vvaccgat_t


    do i = 1,nml

      ! net radiation and precipitation estimates for ctem's bioclim

      uvaccgat_t(i) = uvaccgat_t(i)/real(nday)
      vvaccgat_t(i) = vvaccgat_t(i)/real(nday)
      fsinacc_gat(i) = fsinacc_gat(i)/real(nday)
      flinacc_gat(i) = flinacc_gat(i)/real(nday)
      flutacc_gat(i) = flutacc_gat(i)/real(nday)

      if (altotcount_ctm(ilmos(i)) > 0) then
        altotacc_gat(i) = altotacc_gat(i)/real(altotcount_ctm(ilmos(i)))
      else
        altotacc_gat(i) = 0.
      end if

      fsstar_gat = fsinacc_gat(i) * (1. - altotacc_gat(i))
      flstar_gat = flinacc_gat(i) - flutacc_gat(i)
      netrad_gat(i) = fsstar_gat + flstar_gat

      fsnowacc_t(i) = fsnowacc_t(i)/real(nday)
      taaccgat_t(i) = taaccgat_t(i)/real(nday)

      do j = 1,ignd
        tbaraccgat_t(i,j) = tbaraccgat_t(i,j)/real(nday)
        thliqacc_t(i,j) = thliqacc_t(i,j)/real(nday)
        thiceacc_t(i,j) = thiceacc_t(i,j)/real(nday)
      end do ! loop 831

      do j = 1,icc
        ancgvgac_t(i,j) = ancgvgac_t(i,j)/real(nday)
        rmlcgvga_t(i,j) = rmlcgvga_t(i,j)/real(nday)
      end do ! loop 832

      !     -daily average moss C fluxes for ctem.f-------------------\
      !     Capitulum biomass = 0.22 kg/m2 in hummock, 0.1 kg/m2 in lawn
      !     stem biomass = 1.65 kg/m2 in hummock, 0.77 kg/m2 in lawn (Bragazza et al.2004)
      !     the ratio between stem and capitulum = 7.5 and 7.7
      if (ipeatlandgat(i) > 0) then
        anmossac_t(i) = anmossac_t(i)/real(nday)
        rmlmossac_t(i) = rmlmossac_t(i)/real(nday)
        gppmossac_t(i) = gppmossac_t(i)/real(nday)
      end if

    end do ! nml loop

  end subroutine dayEndCTEMPreparation
  !! @}
  ! --------------------------------------------------------------------------------------------------------------------
  !> \ingroup ctemutilities_accumulateForCTEM
  !! @{
  !> Accumulate the CTEM input (physics) variables at the end of each physics timestep
  !! @author V.Arora, J. Melton
  subroutine accumulateForCTEM (nml, ILMOS)

    use classicParams,  only : icc, ignd, DELT, SBC
    use classStateVars, only : class_gat, class_rot
    use ctemStateVars,  only : vgat, ctem_tile

    implicit none

    integer, intent(in) :: nml     !< Counter representing number of mosaic tiles on modelled domain that are land
    integer, intent(in)  :: ILMOS(:) !< Index of grid cell corresponding to current element
    !< of gathered vector of land surface variables [ ]

    integer :: i, j

    real, pointer :: FSIHGAT(:) !< Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, pointer :: FSVHGAT(:) !< Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, pointer :: ancsmoss(:)
    real, pointer :: angsmoss(:)
    real, pointer :: ancmoss(:)
    real, pointer :: angmoss(:)
    real, pointer :: rmlcsmoss(:)
    real, pointer :: rmlgsmoss(:)
    real, pointer :: rmlcmoss(:)
    real, pointer :: rmlgmoss(:)
    real, pointer :: ALIRGAT(:) !< Diagnosed total near-infrared albedo of land surface [ ]
    real, pointer :: ALVSGAT(:) !< Diagnosed total visible albedo of land surface [ ]
    real, pointer :: FSSROW(:)  !< Shortwave radiation \f$[W m^{-2} ]\f$
    real, pointer :: GTGAT(:)   !< Diagnosed effective surface black-body temperature [K]
    real, pointer :: FDLGAT(:)  !< Downwelling longwave radiation at bottom of atmosphere (i.e. incident on modelled land surface elements \f$[W m^{-2} ]\f$
    real, pointer :: PREGAT(:)  !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer :: FSNOGAT(:) !< Diagnosed fractional snow coverage [ ]
    real, pointer :: TAGAT(:)   !< Air temperature at reference height [K]
    real, pointer :: VLGAT(:)   !< Meridional component of wind velocity \f$[m s^{-1} ]\f$
    real, pointer :: ULGAT(:)   !< Zonal component of wind velocity \f$[m s^{-1} ]\f$
    real, pointer :: FSGGGAT(:) !< Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, pointer :: TBARGAT(:,:) !< Temperature of soil layers [K]
    real, pointer :: FSGSGAT(:) !< Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, pointer :: FSGVGAT(:) !< Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer :: THICGAT(:,:) !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer :: THLQGAT(:,:) !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer :: ancsveggat(:,:)
    real, pointer :: ancgveggat(:,:)
    real, pointer :: rmlcsveggat(:,:)
    real, pointer :: rmlcgveggat(:,:)
    real, pointer :: gppmossgat(:)
    real, pointer :: anmossgat(:)
    real, pointer :: rmlmossgat(:)
    real, pointer :: FC(:)      !<
    real, pointer :: FG(:)      !<
    real, pointer :: FCS(:)     !<
    real, pointer :: FGS(:)     !<
    integer, pointer :: altotcount_ctm(:) ! nlat  !< Counter used for calculating total albedo
    real, pointer :: fsinacc_gat(:) !(ilg)    !<
    real, pointer :: flutacc_gat(:) !(ilg)    !<
    real, pointer :: flinacc_gat(:) !(ilg)    !<
    real, pointer :: altotacc_gat(:) !(ilg)   !<
    real, pointer :: netrad_gat(:) !(ilg)     !<
    real, pointer :: preacc_gat(:) !(ilg)     !<
    real, pointer :: anmossac_t(:)
    real, pointer :: rmlmossac_t(:)
    real, pointer :: gppmossac_t(:)
    real, pointer :: fsnowacc_t(:)       !<
    real, pointer :: taaccgat_t(:)       !<
    real, pointer :: uvaccgat_t(:)       !<
    real, pointer :: vvaccgat_t(:)       !<
    real, pointer :: tbaraccgat_t(:,:) !<
    real, pointer :: thliqacc_t(:,:)  !<
    real, pointer :: thiceacc_t(:,:)  !<
    real, pointer :: ancgvgac_t(:,:)
    real, pointer :: rmlcgvga_t(:,:)
    integer, pointer :: ipeatlandgat(:)

    anmossac_t        => ctem_tile%anmossac_t
    rmlmossac_t       => ctem_tile%rmlmossac_t
    gppmossac_t       => ctem_tile%gppmossac_t
    altotcount_ctm    => vgat%altotcount_ctm
    fsinacc_gat       => vgat%fsinacc_gat
    flutacc_gat       => vgat%flutacc_gat
    flinacc_gat       => vgat%flinacc_gat
    altotacc_gat      => vgat%altotacc_gat
    netrad_gat        => vgat%netrad_gat
    preacc_gat        => vgat%preacc_gat
    ipeatlandgat      => vgat%ipeatland
    tbaraccgat_t      => ctem_tile%tbaraccgat_t
    thliqacc_t        => ctem_tile%thliqacc_t
    thiceacc_t        => ctem_tile%thiceacc_t
    ancgvgac_t        => ctem_tile%ancgvgac_t
    rmlcgvga_t        => ctem_tile%rmlcgvga_t
    fsnowacc_t        => ctem_tile%fsnowacc_t
    taaccgat_t        => ctem_tile%taaccgat_t
    uvaccgat_t        => ctem_tile%uvaccgat_t
    vvaccgat_t        => ctem_tile%vvaccgat_t
    ancsmoss         => vgat%ancsmoss
    angsmoss         => vgat%angsmoss
    ancmoss          => vgat%ancmoss
    angmoss          => vgat%angmoss
    rmlcsmoss        => vgat%rmlcsmoss
    rmlgsmoss        => vgat%rmlgsmoss
    rmlcmoss         => vgat%rmlcmoss
    rmlgmoss         => vgat%rmlgmoss
    FC => class_gat%FC
    FG => class_gat%FG
    FCS => class_gat%FCS
    FGS => class_gat%FGS
    FSIHGAT => class_gat%FSIHGAT
    FSVHGAT => class_gat%FSVHGAT
    ALIRGAT => class_gat%ALIRGAT
    ALVSGAT => class_gat%ALVSGAT
    FSSROW => class_rot%FSSROW
    GTGAT => class_gat%GTGAT
    FDLGAT => class_gat%FDLGAT
    PREGAT => class_gat%PREGAT
    FSNOGAT => class_gat%FSNOGAT
    TAGAT => class_gat%TAGAT
    VLGAT => class_gat%VLGAT
    ULGAT => class_gat%ULGAT
    FSGGGAT => class_gat%FSGGGAT
    TBARGAT => class_gat%TBARGAT
    FSGSGAT => class_gat%FSGSGAT
    FSGVGAT => class_gat%FSGVGAT
    THICGAT => class_gat%THICGAT
    THLQGAT => class_gat%THLQGAT
    ancsveggat        => vgat%ancsveg
    ancgveggat        => vgat%ancgveg
    rmlcsveggat       => vgat%rmlcsveg
    rmlcgveggat       => vgat%rmlcgveg
    anmossgat        => vgat%anmoss
    rmlmossgat       => vgat%rmlmoss
    gppmossgat       => vgat%gppmoss

    do i = 1,nml

      fsinacc_gat(i) = fsinacc_gat(i) + FSSROW(ilmos(I)) 
      flinacc_gat(i) = flinacc_gat(i) + fdlgat(i)
      flutacc_gat(i) = flutacc_gat(i) + sbc * gtgat(i) ** 4
      preacc_gat(i) = preacc_gat(i) + pregat(i) * delt
      fsnowacc_t(i) = fsnowacc_t(i) + fsnogat(i)
      taaccgat_t(i) = taaccgat_t(i) + tagat(i)
      vvaccgat_t(i) = vvaccgat_t(i) + vlgat(i)
      uvaccgat_t(i) = uvaccgat_t(i) + ulgat(i)
      if (FSSROW(ilmos(I)) > 0.) then
        altotacc_gat(i) = altotacc_gat(i) + (FSSROW(ilmos(I)) - &
                          (FSGVGAT(I) + FSGSGAT(I) + FSGGGAT(I))) &
                          /FSSROW(ilmos(I))
        altotcount_ctm = altotcount_ctm + 1
      end if

      do j = 1,ignd
        tbaraccgat_t(i,j) = tbaraccgat_t(i,j) + tbargat(i,j)
        thliqacc_t(i,j) = thliqacc_t(i,j) + THLQGAT(i,j)
        thiceacc_t(i,j) = thiceacc_t(i,j) + THICGAT(i,j)
      end do ! loop 710

      do j = 1,icc
        ancgvgac_t(i,j) = ancgvgac_t(i,j) + (1. - fsnogat(i)) * ancgveggat(i,j) + fsnogat(i) * ancsveggat(i,j)
        rmlcgvga_t(i,j) = rmlcgvga_t(i,j) + (1. - fsnogat(i)) * rmlcgveggat(i,j) + fsnogat(i) * rmlcsveggat(i,j)
      end do ! loop 713

      !    -accumulate moss C fluxes to tile level then daily----
      if (ipeatlandgat(i) > 0) then
        anmossgat(i) = fcs(i) * ancsmoss(i) + fgs(i) * angsmoss(i) + fc(i) * ancmoss(i) + fg(i) * angmoss(i)
        rmlmossgat(i) = fcs(i) * rmlcsmoss(i) + fgs(i) * rmlgsmoss(i) + fc(i) * rmlcmoss(i) + fg(i) * rmlgmoss(i)
        gppmossgat(i) = anmossgat(i) + rmlmossgat(i)

        anmossac_t(i) = anmossac_t(i)   + anmossgat(i)
        rmlmossac_t(i) = rmlmossac_t(i)  + rmlmossgat(i)
        gppmossac_t(i) = gppmossac_t(i)  + gppmossgat(i)

      end if
    end do

  end subroutine accumulateForCTEM
  !! @}
  !
  ! --------------------------------------------------------------------------------------------------------------------
  !> \ingroup ctemutilities_ctemInit
  !! @{
  !> Find mosaic tile (grid) average vegetation biomass, litter mass, and soil c mass.
  !! Also initialize additional variables which are used by CTEM (biogeochemical processes).
  !! @author V.Arora, J. Melton
  subroutine ctemInit (nltest, nmtest)

    use classicParams,  only : icc, ilg, ignd, iccp1
    use ctemStateVars,  only : vrot, ctem_tile, vgat
    use classStateVars, only : class_rot
    use generalUtils,   only : findDaylength
    use peatlandsMod,   only : peatStorage

    implicit none

    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest

    integer :: i, m, j, k

    real, pointer, dimension(:,:,:) :: co2i1cgrow
    real, pointer, dimension(:,:,:) :: co2i1csrow
    real, pointer, dimension(:,:,:) :: co2i2cgrow
    real, pointer, dimension(:,:,:) :: co2i2csrow
    real, pointer, dimension(:,:,:) :: slairow
    real, pointer, dimension(:) :: fsnowacc_t
    real, pointer, dimension(:) :: taaccgat_t
    real, pointer, dimension(:,:) :: ancgvgac_t
    real, pointer, dimension(:,:) :: rmlcgvga_t
    real, pointer, dimension(:,:)  :: todfrac
    real, pointer, dimension(:,:) :: thliqacc_t
    real, pointer, dimension(:,:) :: thiceacc_t
    real, pointer, dimension(:,:) :: vgbiomasrow
    real, pointer, dimension(:,:) :: gavgltmsrow
    real, pointer, dimension(:,:) :: gavgscmsrow
    real, pointer, dimension(:,:) :: gavglairow
    real, pointer, dimension(:,:) :: lucemcomrow
    real, pointer, dimension(:,:) :: lucltrinrow
    real, pointer, dimension(:,:) :: lucsocinrow
    integer, pointer, dimension(:,:,:) :: colddaysrow
    real, pointer, dimension(:,:,:) :: flhrlossrow
    real, pointer, dimension(:,:,:) :: stmhrlosrow
    real, pointer, dimension(:,:,:) :: rothrlosrow
    real, pointer, dimension(:,:,:) :: tymaxlairow
    real, pointer, dimension(:,:,:) :: grwtheffrow
    real, pointer, dimension(:,:,:) :: lystmmasrow
    real, pointer, dimension(:,:,:) :: lyrotmasrow
    integer, pointer, dimension(:,:) :: ipeatlandrow ! This is first set in read_from_ctm.
    real, pointer, dimension(:,:,:) :: fcanrot !<
    real, pointer, dimension(:,:,:) :: gleafmasrow        !
    real, pointer, dimension(:,:,:) :: bleafmasrow        !
    real, pointer, dimension(:,:,:) :: stemmassrow        !
    real, pointer, dimension(:,:,:) :: rootmassrow        !
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: litrmassrow
    real, pointer, dimension(:,:,:) :: soilcmasrow
    ! real, pointer, dimension(:,:,:,:) :: litrmassrow
    ! real, pointer, dimension(:,:,:,:) :: soilcmasrow
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:) :: peatdeprow
    real, pointer, dimension(:) :: anmossac_t
    real, pointer, dimension(:) :: rmlmossac_t
    real, pointer, dimension(:) :: gppmossac_t
    real, pointer, dimension(:,:) :: litrmsmossrow
    real, pointer, dimension(:,:) :: Cmossmasrow
    real, pointer, dimension(:) :: dayl_maxrow
    real, pointer, dimension(:,:) :: sdeprot    !< Depth to bedrock in the soil profile
    real, pointer, dimension(:) :: RADJROW      !< Latitude of grid cell (positive north of equator) [rad]

    ipeatlandrow      => vrot%ipeatland
    co2i1cgrow        => vrot%co2i1cg
    co2i1csrow        => vrot%co2i1cs
    co2i2cgrow        => vrot%co2i2cg
    co2i2csrow        => vrot%co2i2cs
    slairow           => vrot%slai
    fsnowacc_t        => ctem_tile%fsnowacc_t
    taaccgat_t        => ctem_tile%taaccgat_t
    ancgvgac_t        => ctem_tile%ancgvgac_t
    rmlcgvga_t        => ctem_tile%rmlcgvga_t
    todfrac           => vgat%todfrac
    thliqacc_t        => ctem_tile%thliqacc_t
    thiceacc_t        => ctem_tile%thiceacc_t
    vgbiomasrow       => vrot%vgbiomas
    gavglairow        => vrot%gavglai
    gavgltmsrow       => vrot%gavgltms
    gavgscmsrow       => vrot%gavgscms
    lucemcomrow       => vrot%lucemcom
    lucltrinrow       => vrot%lucltrin
    lucsocinrow       => vrot%lucsocin
    colddaysrow       => vrot%colddays
    flhrlossrow       => vrot%flhrloss
    stmhrlosrow       => vrot%stmhrlos
    rothrlosrow       => vrot%rothrlos
    tymaxlairow       => vrot%tymaxlai
    grwtheffrow       => vrot%grwtheff
    lystmmasrow       => vrot%lystmmas
    lyrotmasrow       => vrot%lyrotmas
    fcanrot           => class_rot%fcanrot
    gleafmasrow       => vrot%gleafmas
    bleafmasrow       => vrot%bleafmas
    stemmassrow       => vrot%stemmass
    rootmassrow       => vrot%rootmass
    litrmassrow       => vrot%litrmass
    soilcmasrow       => vrot%soilcmas
    fcancmxrow        => vrot%fcancmx
    peatdeprow        => vrot%peatdep
    anmossac_t        => ctem_tile%anmossac_t
    rmlmossac_t       => ctem_tile%rmlmossac_t
    gppmossac_t       => ctem_tile%gppmossac_t
    litrmsmossrow     => vrot%litrmsmoss
    Cmossmasrow       => vrot%Cmossmas
    dayl_maxrow       => vrot%dayl_max
    sdeprot           => class_rot%sdeprot
    RADJROW           => class_rot%RADJROW

    ! ------

    ! Initialize to zero:
    co2i1csrow(:,:,:) = 0.0     ! intercellular co2 concentrations
    co2i1cgrow(:,:,:) = 0.0
    co2i2csrow(:,:,:) = 0.0
    co2i2cgrow(:,:,:) = 0.0
    slairow(:,:,:) = 0.0        ! if bio2str is not called we need to initialize this to zero
    fsnowacc_t(:) = 0.0         ! daily accu. fraction of snow
    taaccgat_t(:) = 0.0
    ancgvgac_t(:,:) = 0.0    ! daily accu. net photosyn.
    rmlcgvga_t(:,:) = 0.0    ! daily accu. leaf respiration
    todfrac(:,:) = 0.0
    thliqacc_t(:,:) = 0.0
    thiceacc_t(:,:) = 0.0
    vgbiomasrow(:,:) = 0.0
    gavglairow(:,:) = 0.0
    gavgltmsrow(:,:) = 0.0
    gavgscmsrow(:,:) = 0.0
    lucemcomrow(:,:) = 0.0      ! land use change combustion emission losses
    lucltrinrow(:,:) = 0.0      ! land use change inputs to litter pool
    lucsocinrow(:,:) = 0.0      ! land use change inputs to soil c pool
    colddaysrow(:,:,:) = 0      ! cold days counter for ndl dcd and crops
    flhrlossrow(:,:,:) = 0.0     ! fall/harvest loss
    stmhrlosrow(:,:,:) = 0.0     ! stem harvest loss for crops
    rothrlosrow(:,:,:) = 0.0     ! root death for crops
    tymaxlairow(:,:,:) = 0.0

    do i = 1,nltest ! loop 115
      do m = 1,nmtest
        do j = 1,icc
          vgbiomasrow(i,m) = vgbiomasrow(i,m) + fcancmxrow(i,m,j) * &
                             (gleafmasrow(i,m,j) + stemmassrow(i,m,j) + &
                             rootmassrow(i,m,j) + bleafmasrow(i,m,j))
          ! COMBAK PERLAY
          gavgltmsrow(i,m) = gavgltmsrow(i,m) + fcancmxrow(i,m,j) * &
                             litrmassrow(i,m,j)
          gavgscmsrow(i,m) = gavgscmsrow(i,m) + fcancmxrow(i,m,j) * &
                             soilcmasrow(i,m,j)
          ! do k = 1,ignd
          !   gavgltmsrow(i,m)=gavgltmsrow(i,m)+fcancmxrow(i,m,j)* &
          !       &                       litrmassrow(i,m,j,k)
          !   gavgscmsrow(i,m)=gavgscmsrow(i,m)+fcancmxrow(i,m,j)* &
          !       &         soilcmasrow(i,m,j,k)
          ! end do ! ignd
          ! COMBAK PERLAY
          !grwtheffrow(i,m,j) = 100.0   ! set growth efficiency to some large number
          ! so that no growth related mortality occurs in
          ! first year
          lystmmasrow(i,m,j) = stemmassrow(i,m,j)
          lyrotmasrow(i,m,j) = rootmassrow(i,m,j)

        end do ! loop 116
      end do
    end do ! loop 115

    do i = 1,nltest ! loop 117
      do m = 1,nmtest
        if (ipeatlandrow(i,m) == 0) then ! NON-peatland tile
          ! COMBAK PERLAY
          gavgltmsrow(i,m) = gavgltmsrow(i,m) + (1.0 - sum(fcanrot(i,m,:))) * litrmassrow(i,m,iccp1)
          gavgscmsrow(i,m) = gavgscmsrow(i,m) + (1.0 - sum(fcanrot(i,m,:))) * soilcmasrow(i,m,iccp1)
          ! do k = 1,ignd
          !   gavgltmsrow(i,m)=gavgltmsrow(i,m)+ (1.0-sum(fcanrot(i,m,:)))*litrmassrow(i,m,iccp1,k)
          !   gavgscmsrow(i,m)=gavgscmsrow(i,m)+ (1.0-sum(fcanrot(i,m,:)))*soilcmasrow(i,m,iccp1,k)
          ! end do
          ! COMBAK PERLAY
        else ! peatland tile
          gavgltmsrow(i,m) = gavgltmsrow(i,m) + litrmsmossrow(i,m)
          peatdeprow(i,m) = sdeprot(i,m) ! the peatdepth is set to the soil depth

          ! The soil carbon on the peatland tiles is assigned based on depth. This
          ! is the same relation as found in hetresPeat subroutine.
          gavgscmsrow(i,m) = peatStorage(peatdeprow(i,m))
          ! gavgscmsrow(i,m) = 0.487*(4056.6*peatdeprow(i,m)**2+ &
          !             72067.0*peatdeprow(i,m))/1000

          vgbiomasrow(i,m) = vgbiomasrow(i,m) + Cmossmasrow(i,m)
        end if
      end do
    end do ! loop 117

    !    Also initialize the accumulators for moss daily C fluxes.

    anmossac_t  = 0.0
    rmlmossac_t = 0.0
    gppmossac_t = 0.0

    !    ----------------------------YW March 25, 2015 --------------------/

    ! Lastly, find the maximum daylength at this location for day 172 = June 21st - summer solstice.
    do i = 1,nltest
      if (radjrow(i) > 0.) then
        dayl_maxrow(i) = findDaylength(172.0,radjrow(i)) ! following rest of code, radjrow is always given index of 1 offline.
      else ! S. Hemi so do N.Hemi winter solstice Dec 21
        dayl_maxrow(i) = findDaylength(355.0,radjrow(i)) ! following rest of code, radjrow is always given index of 1 offline.
      end if
    end do

  end subroutine ctemInit
  !! @}
  !> \namespace ctemutilities
  !! Central module for CTEM (biogeochemical)-related utilities
end module ctemUtilities
