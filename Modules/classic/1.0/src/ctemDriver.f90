!> \file
!> Principle driver for CTEM
!! @author V. Arora, J. Melton, Y. Peng, R. Shrestha
!!
module ctemDriver

  implicit none

  ! subroutines contained in this module:

  public  :: ctem
  public  :: calcNPP
  public  :: calcNEP
  public  :: calcNBP

contains
  !> \ingroup ctemdriver_ctem
  !! @{
  subroutine ctem (fsnow, sand, ilmos, & ! In
                   ilg, il1, il2, iday, radj, & ! In
                   ta, delzw, ancgveg, rmlcgveg, & ! In
                   zbotw, doMethane, & ! In
                   uwind, vwind, lightng, tbar, & ! In
                   soildpth, spinfast, todfrac, & ! In
                   netrad, precip, psisat, & ! In
                   grclarea, popdin, isand, & ! In
                   wetfrac, slopefrac, bi, & ! In
                   thpor, currlat, ch4conc, & ! In
                   THFC, THLW, thliq, thice, & ! In
                   ipeatland, anmoss, rmlmoss, gppmoss, & ! In
                   wtable, maxAnnualActLyr, & ! In
                   PFTCompetition, dofire, lnduseon, inibioclim, & ! In
                   leapnow, useTracer, tracerCO2, & ! In
                   pfcancmx, nfcancmx, & ! In/Out
                   stemmass, rootmass, litrmass, gleafmas, & ! In/ Out
                   bleafmas, soilcmas, ailcg, ailc, & ! In/ Out
                   zolnc, rmatctem, rmatc, ailcb, & ! In/ Out
                   flhrloss, pandays, lfstatus, grwtheff, & ! In/ Out
                   lystmmas, lyrotmas, tymaxlai, vgbiomas, & ! In/ Out
                   gavgltms, gavgscms, stmhrlos, slai, & ! In/ Out
                   bmasveg, cmasvegc, colddays, rothrlos, & ! In/ Out
                   fcanmx, alvisc, alnirc, gavglai, & ! In/ Out
                   Cmossmas, litrmsmoss, peatdep, fcancmx, &! In/ Out
                   geremort, intrmort, pstemmass, pgleafmass, &! In/ Out
                   tcurm, srpcuryr, dftcuryr, lambda, &! In/ Out
                   tmonth, anpcpcur, anpecur, gdd5cur, &! In/ Out
                   surmncur, defmncur, srplscur, defctcur, &! In/ Out
                   aridity, srplsmon, defctmon, anndefct, &! In/ Out
                   annsrpls, annpcp, dry_season_length, &! In/ Out
                   pftexist, twarmm, tcoldm, gdd5, nppveg,&! In/ Out
                   tracerStemMass, tracerRootMass, tracerGLeafMass, tracerBLeafMass, & ! In/Out
                   tracerSoilCMass, tracerLitrMass, tracerMossCMass, tracerMossLitrMass, & ! In/Out
                   npp, nep, hetrores, autores, & ! Out (Primary)
                   soilresp, rm, rg, nbp, & ! Out (Primary)
                   litres, socres, gpp, dstcemls1, & ! Out (Primary)
                   litrfall, humiftrs, veghght, rootdpth, & ! Out (Primary)
                   rml, rms, rmr, tltrleaf, & ! Out (Primary)
                   tltrstem, tltrroot, leaflitr, roottemp, & ! Out (Primary)
                   burnfrac, lucemcom, lucltrin, & ! Out (Primary)
                   lucsocin, dstcemls3, & ! Out (Primary)
                   ch4WetSpec, ch4WetDyn, wetfdyn, ch4soills, & ! Out (Primary)
                   paicgat, slaicgat, & ! Out (Primary)
                   emit_co2, emit_ch4, reprocost, blfltrdt, glfltrdt, &  ! Out (Primary)
                   glcaemls, blcaemls, rtcaemls, stcaemls, ltrcemls, &  ! Out (Primary)
                   ntchlveg, ntchsveg, ntchrveg, &  ! Out (Primary)
                   emit_co, emit_nmhc, smfunc_veg, & ! Out (Secondary)
                   emit_h2, emit_nox, emit_n2o, emit_pm25, & ! Out (Secondary)
                   emit_tpm, emit_tc, emit_oc, emit_bc, & ! Out (Secondary)
                   bterm_veg, lterm, mterm_veg, burnvegf, & ! Out (Secondary)
                   litrfallveg, humtrsvg, ltstatus, & ! Out (Secondary)
                   afrleaf, afrstem, afrroot, wtstatus, & ! Out (Secondary)
                   rmlveg, rmsveg, rmrveg, rgveg, & ! Out (Secondary)
                   vgbiomas_veg, gppveg, nepveg, nbpveg, & ! Out (Secondary)
                   hetrsveg, autoresveg, ltresveg, scresveg, & ! Out (Secondary)
                   nppmoss, armoss, & ! Out (Secondary)
                   colrate, mortrate) ! Out (Secondary)
    !
    !             Canadian Terrestrial Ecosystem Model (CTEM)
    !             Main Ctem Subroutine Compatible With CLASS

    !     29 Mar 2019 - Rewrote for clarity, removed sub-area TBAR, THLIQ, THICE, ANVEG, rmlveg,
    !     J. Melton     was needlessly confusing and didn't have impact. Added subroutines to make
    !                   this a driver, not a main calculation subroutine.
    !     28  Nov 2018  - Clean up argument list before implementation in AGCM
    !     V. Arora
    !
    !     14  Mar 2016  - Remove grid cell area calculation to the driver. This will
    !     J. Melton       harmonize this subroutine with the coupled code.
    !
    !      8  Feb 2016  - Adapted subroutine for multilayer soilc and litter (fast decaying)
    !     J. Melton       carbon pools
    !
    !     3   Feb 2016  - Bring in onetile_perPFT switch so now in mosaic mode the
    !     J. Melton       tiles in a grid cell are treated independently by competition
    !                     and LUC.
    !     19  Jan 2016  - Implemented new LUC litter and soil C pools
    !     J. Melton
    !
    !     3   Jul 2014  - Bring in wetland and wetland methane code
    !     R. Shrestha
    !
    !     12  Jun 2014  - Bring in a constant reproductive cost, remove expnbaln,
    !     J. Melton       add a smoothing function for lambda calculation for competition,
    !                     made it so NEP and NBP work with competition on.
    !
    !     17  Jan 2014  - Moved parameters to global file (classicParams.f90)
    !     J. Melton
    !
    !     Dec 6   2012   Make it so competition and luc can function in both
    !     J. Melton      composite and mosaic modes.
    !
    !     sep 25  2012   Add competition_map and competition_unmap
    !     Y. Peng
    !
    !     Sep 12  2012
    !     J. Melton     Add in bottom limit to bleafmas to ensure it does no
    !                   slowly decay to infintesimaly small number.
    !
    !     Aug 23  2012
    !     J. Melton     Pass in isand to ensure soil levels are properly
    !                   labelled as bedrock if assigned so in classb
    !
    !     Jan 10  2012
    !     Yiran         Re-test bioclim and existence for competition
    !
    !     19  Sep. 2001 - This is the main terrestrial carbon model subrouti
    !     V. Arora
    !                     all primary ctem subroutines are called from here,
    !                     except phtsyn which is called from tsolvc
    !
    !     08 June  2001 - Add calls to three new subroutines (bioclim,
    !     V. Arora        existence, and competition) to add dynamic
    !                     competition between pfts. changes are also made
    !                     to land use change (luc) subroutine and the manner
    !                     in which distrubance is handles. fire now creates
    !                     bare ground which is subsequently available for
    !                     colonization. other changes are also made to keep
    !                     everything consistent with changing vegetation
    !                     fractions.
    !    -----------------------------------------------------------------

    use classicParams,             only : kk, pi, zero, icp1, &
                                          iccp1, ican, nlat, &
                                          ignd, icc, nmos, l2max, grescoef, &
                                          humicfac, laimin, laimax, lambdamax, &
                                          crop, repro_fraction, &
                                          rmortmoss, humicfacmoss, GRAV, RHOW, RHOICE, &
                                          classpfts, ctempfts, iccp2, humicfac_bg, &
                                          nol2pfts, deltat, tolrance
    use landuseChange,             only : luc
    use competitionScheme,         only : bioclim, existence, competition, expansion
    use disturbance_scheme,        only : disturb
    use heterotrophicRespirationMod, only : heterotrophicRespiration, updatePoolsHetResp
    use peatlandsMod,              only : peatDayEnd, peatDepth
    use ctemUtilitiesGAT,             only : genSortIndex
    use autotrophicRespiration,    only : mainres
    use balanceCarbon,             only : balcar, prepBalanceC
    use mortality,                 only : mortalty, updatePoolsMortality
    use turnover,                  only : turnoverStemRoot, updatePoolsTurnover
    use applyAllometry,            only : allometry
    use tracerModule,              only : prepTracer, doTracerBalance, checkTracerBalance
    use methaneProcesses,          only : soil_ch4uptake, wetland_methane

    ! COMBAK PERLAY
    ! use soilCProcesses, only : turbation
    ! COMBAK PERLAY
    use allocateCarbon,            only : allocate, updatePoolsAllocateRepro

    implicit none

    ! inputs

    logical, intent(in) :: lnduseon                         !< logical switch to run the land use change subroutine or not.
    logical, intent(in) :: PFTCompetition                   !< logical boolean telling if competition between pfts is on or not
    logical, intent(in) :: dofire                           !< boolean, if true allow fire, if false no fire.
    logical, intent(in) :: leapnow                          !< true if this year is a leap year. Only used if the switch 'leap' is true.
    logical, intent(in) :: doMethane                        !< true if you wish to do the methane calculations.
    integer, intent(in) :: iday                             !< day of year
    integer, intent(in) ::  spinfast                        !< spinup factor for soil carbon whose default value is 1. as this factor increases the
    !< soil c pool will come into equilibrium faster. reasonable value for spinfast is
    !< between 5 and 10. when spinfast/=1 then the balcar subroutine is not run.
    integer, intent(in) :: ilg                              !< ilg=no. of grid cells in latitude circle
    integer, intent(in) :: il1                              !< il1=1
    integer, intent(in) :: il2                              !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in), dimension(:) :: ilmos              !< Index of gridcell corresponding to current element of gathered vector of land surface variables [ ]
    integer, dimension(ilg,ignd), intent(in) :: isand       !<
    integer, dimension(ilg), intent(in) :: ipeatland        !< Peatland flag: 0 = not a peatland, 1 = bog, 2 = fen
    real, dimension(ilg), intent(in) :: fsnow               !< fraction of snow simulated by class
    real, dimension(ilg,ignd), intent(in) :: sand           !< percentage sand
    real, dimension(ilg), intent(in) :: radj                !< latitude in radians
    real, dimension(ilg,ignd), intent(in) ::  tbar          !< Soil temperature, K
    real, dimension(ilg,ignd), intent(in) :: psisat         !< Saturated soil matric potential (m)
    real, dimension(ilg,ignd), intent(in) :: bi             !< Brooks and Corey/Clapp and Hornberger b term
    real, dimension(ilg,ignd), intent(in) :: thpor          !< Soil total porosity \f$(cm^3 cm^{-3})\f$ - daily average
    real, dimension(ilg), intent(in) :: ta                  !< air temp, K
    real, dimension(ilg,ignd), intent(in) :: delzw          !< thicknesses of the soil layers
    real, dimension(ilg,ignd), intent(in) :: zbotw          !< bottom of soil layers
    real, dimension(ilg), intent(in) :: soildpth            !< soil depth (m)
    real, dimension(ilg,ignd), intent(in) :: thliq          !< liquid mois. content of soil layers
    real, dimension(ilg,ignd), intent(in) :: thice          !< Frozen soil moisture content
    real, dimension(ilg), intent(in) ::  grclarea           !< area of the grid cell, \f$km^2\f$
    real, dimension(ilg), intent(in) ::  currlat            !< centre latitude of grid cells in degrees
    real, dimension(ilg), intent(in) :: uwind               !< u wind speed, m/s
    real, dimension(ilg), intent(in) :: vwind               !< v wind speed, m/s
    real, dimension(ilg), intent(in) ::  precip             !< daily precipitation (mm/day)
    real, dimension(ilg), intent(in) ::  netrad             !< daily net radiation (w/m2)
    real, dimension(ilg), intent(in) :: lightng             !< total lightning frequency, flashes/km2.year
    real, dimension(ilg,icc), intent(in) :: todfrac         !< max. fractional coverage of ctem's 9 pfts by the end of the day, for use by land use subroutine
    real, dimension(ilg), intent(in) :: ch4conc             !< Atmospheric \f$CH_4\f$ concentration at the soil surface (ppmv)
    real, dimension(ilg), intent(in) :: wetfrac             !< Prescribed fraction of wetlands in a grid cell
    real, dimension(ilg,8), intent(in) :: slopefrac         !<
    real, dimension(ilg), intent(in) :: anmoss              !< moss net photoysnthesis -daily averaged C fluxes rates (umol/m2/s)
    real, dimension(ilg), intent(in) :: rmlmoss             !< moss maintainance respiration -daily averaged C fluxes rates (umol/m2/s)
    real, dimension(ilg), intent(in) :: gppmoss             !< moss GPP -daily averaged C fluxes rates (umol/m2/s)
    real, dimension(ilg), intent(in) :: wtable              !< water table (m)
    real, dimension(ilg,icc), intent(in) :: ancgveg         !< net photosynthetic rate for CTEM's pfts
    real, dimension(ilg,icc), intent(in) :: rmlcgveg        !< leaf respiration rate for CTEM's pfts
    real, dimension(ilg,ignd), intent(in) :: THFC           !<
    real, dimension(ilg,ignd), intent(in) :: THLW           !<
    real, dimension(ilg), intent(in) :: maxAnnualActLyr     !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    real, intent(in) :: tracerCO2(:)       !< Tracer CO2 value read in from tracerCO2File, units vary (simple: ppm, 14C \f$\Delta ^{14}C\f$)
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.

    !
    !     updates
    !
    logical, intent(inout) :: pftexist(ilg,icc)             !< True if PFT is present in tile.
    logical, intent(inout) :: inibioclim                    !< switch telling if bioclimatic parameters are being initialized from scratch (false)
    !< or being initialized from some spun up values(true).
    integer, dimension(ilg,icc), intent(inout) :: pandays   !< days with positive net photosynthesis (an) for use in the phenology subroutine
    integer, dimension(ilg,2), intent(inout) :: colddays    !< cold days counter for tracking days below a certain temperature threshold for ndl dcd and crop pfts.
    integer, dimension(ilg,icc), intent(inout) :: lfstatus  !< leaf phenology status
    real, dimension(ilg,icc), intent(inout) :: fcancmx      !< max. fractional coverage of CTEM's pfts, but this can be
    !< modified by land-use change,and competition between pfts
    real, dimension(ilg,icc), intent(inout) :: pfcancmx        !< previous year's fractional coverages of pfts
    real, dimension(ilg,icc), intent(inout) :: nfcancmx        !< next year's fractional coverages of pfts
    real, dimension(ilg,ican,ignd), intent(inout) :: rmatc  !< fraction of roots for each of class' 4 pfts in each soil layer
    real, dimension(ilg), intent(inout) :: surmncur         !< number of months with surplus water for current year
    real, dimension(ilg), intent(inout) :: defmncur         !< number of months with water deficit for current year
    real, dimension(ilg), intent(inout) :: tcurm            !< temperature of the current month (c)
    real, dimension(ilg), intent(inout) :: annpcp           !< annual precipitation (mm)
    real, dimension(ilg), intent(inout) :: dry_season_length !< length of the dry season (months)
    real, dimension(ilg), intent(inout) :: twarmm           !< temperature of the warmest month (c)
    real, dimension(ilg), intent(inout) :: tcoldm           !< temperature of the coldest month (c)
    real, dimension(ilg), intent(inout) :: gdd5             !< growing degree days above 5 c
    real, dimension(ilg), intent(inout) :: aridity          !< aridity index, ratio of potential evaporation to precipitation
    real, dimension(ilg), intent(inout) :: srplsmon         !< number of months in a year with surplus water i.e. precipitation more than potential evaporation
    real, dimension(ilg), intent(inout) :: defctmon         !< number of months in a year with water deficit i.e. precipitation less than potential evaporation
    real, dimension(12,ilg), intent(inout) :: tmonth        !< monthly temperatures
    real, dimension(ilg), intent(inout) :: anpcpcur         !< annual precipitation for current year (mm)
    real, dimension(ilg), intent(inout) :: anpecur          !< annual potential evaporation for current year (mm)
    real, dimension(ilg), intent(inout) :: gdd5cur          !< growing degree days above 5 c for current year
    real, dimension(ilg), intent(inout) :: srplscur         !< water surplus for the current month
    real, dimension(ilg), intent(inout) :: defctcur         !< water deficit for the current month
    real, dimension(ilg), intent(inout) :: srpcuryr         !< water surplus for the current year
    real, dimension(ilg), intent(inout) :: dftcuryr         !< water deficit for the current year
    real, dimension(ilg), intent(inout) :: anndefct         !< annual water deficit (mm)
    real, dimension(ilg), intent(inout) :: annsrpls         !< annual water surplus (mm)
    real, dimension(ilg,icc), intent(inout) :: stemmass     !< stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: rootmass     !< root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    ! COMBAK PERLAY
    ! real, dimension(ilg,iccp2,ignd), intent(inout) :: litrmass   !< litter mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, dimension(ilg,iccp2), intent(inout) :: litrmass   !< litter mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$ ! COMBAK PERLAY
    ! COMBAK PERLAY
    real, dimension(ilg,icc), intent(inout) :: gleafmas     !< green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: bleafmas     !< brown leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    ! COMBAK PERLAY
    ! real, dimension(ilg,iccp2,ignd), intent(inout) :: soilcmas   !< soil carbon mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, dimension(ilg,iccp2), intent(inout) :: soilcmas   !< soil carbon mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$ ! COMBAK PERLAY
    ! COMBAK PERLAY
    real, dimension(ilg,icc), intent(inout) :: ailcg        !< Green LAI for ctem's pfts \f$(m^2 leaf/m^2 ground)\f$
    real, dimension(ilg,ican), intent(inout) :: ailc        !< lumped lai for class' 4 pfts
    real, dimension(ilg,icc,ignd), intent(inout) :: rmatctem !< fraction of roots for each of ctem's 9 pfts in each soil layer
    real, dimension(ilg,ican), intent(inout) :: zolnc       !< lumped log of roughness length for class' 4 pfts
    real, dimension(ilg,icc), intent(inout) :: ailcb        !< brown lai for ctem's 9 pfts. for now we assume only grasses can have brown lai
    real, dimension(ilg), intent(inout) :: vgbiomas         !< grid averaged vegetation biomass, \f$(kg C/m^2)\f$
    real, dimension(ilg), intent(inout) :: gavgltms         !< grid averaged litter mass, \f$(kg C/m^2)\f$
    real, dimension(ilg), intent(inout) :: gavgscms         !< grid averaged soil c mass, \f$(kg C/m^2)\f$
    real, dimension(ilg), intent(inout) :: gavglai          !< grid averaged green leaf area index
    real, dimension(ilg,icc), intent(inout) :: bmasveg      !< total (gleaf + stem + root) biomass for each ctem pft, \f$(kg C/m^2)\f$
    real, dimension(ilg,ican), intent(inout) :: cmasvegc    !< total canopy mass for each of the 4 class pfts. recall that class requires canopy
    !< mass as an input,and this is now provided by ctem. \f$kg/m^2\f$.
    real, dimension(ilg,icp1), intent(inout) :: fcanmx      !< fractional coverage of class' 4 pfts
    real, dimension(ilg,ican), intent(inout) :: alvisc      !< visible albedo for class' 4 pfts
    real, dimension(ilg,ican), intent(inout) :: alnirc      !< near ir albedo for class' 4 pfts
    real, dimension(ilg,icc), intent(inout) :: pstemmass    !< stem mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(ilg,icc), intent(inout) :: pgleafmass   !< root mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(ilg,icc), intent(inout) :: flhrloss     !< fall or harvest loss for deciduous trees and crops, respectively, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: stmhrlos     !< stem harvest loss for crops, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: rothrlos     !< root death as crops are harvested, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: grwtheff     !< growth efficiency. change in biomass per year per unit max. lai (\f$(kg C/m^2)\f$)/(m2/m2),
    !< for use in mortality subroutine
    real, dimension(ilg,icc), intent(inout) :: lystmmas     !< stem mass at the end of last year
    real, dimension(ilg,icc), intent(inout) :: lyrotmas     !< root mass at the end of last year
    real, dimension(ilg,icc), intent(inout) :: tymaxlai     !< this year's maximum lai
    real, dimension(ilg,icc), intent(inout) ::  geremort    !<
    real, dimension(ilg,icc), intent(inout) :: intrmort     !<
    real, dimension(ilg,icc), intent(inout) ::  burnvegf    !< per PFT fraction burned of that PFT's area
    real, dimension(ilg), intent(inout) :: popdin           !< population density \f$(people / km^2)\f$
    real, dimension(ilg), intent(inout) :: Cmossmas         !< moss biomass C pool (kgC/m2)
    real, dimension(ilg), intent(inout) :: litrmsmoss       !< moss litter C pool (kgC/m2)
    real, dimension(ilg), intent(inout) :: peatdep          !< peat depth (m)

    real, intent(inout) :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$tracer C units/m^2\f$
    real, intent(inout) :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, intent(inout) :: tracerMossCMass(:)      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, intent(inout) :: tracerMossLitrMass(:)   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, dimension(ilg,icc), intent(inout) :: slai           !< storage/imaginary lai for phenology purposes
    real, dimension(ilg,icc), intent(inout) :: nppveg         !< NPP for individual pfts, (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(inout) :: lambda         !< Fraction of npp that is to be used for
    !! horizontal expansion (lambda) during the next
    !! day (i.e. this will be determining
    !! the colonization rate in competition)
    !
    !   outputs
    !
    real, dimension(ilg), intent(out) :: ch4soills          !< Methane uptake into the soil column \f$(mg CH_4 m^{-2} s^{-1})\f$
    real, dimension(ilg), intent(out) :: rml                !< Tile level leaf maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: gpp                !< Tile level gross primary productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real, dimension(ilg,icc), intent(out) :: veghght        !< vegetation height (meters)
    real, dimension(ilg,icc), intent(out) :: rootdpth       !< 99% soil rooting depth (meters) both veghght & rootdpth can be used as diagnostics
    !< to see how vegetation grows above and below ground, respectively
    real, dimension(ilg), intent(out) :: npp                !< Tile-level net primary productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: nep                !< Tile-level net ecosystem productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: hetrores           !< Tile-level heterotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: autores             !< Tile level autotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: soilresp           !< Soil respiration. This includes root respiration and respiration from
    !! litter and soil carbon pools. Note that soilresp is different from
    !! socres,which is respiration from the soil C pool.(\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: rm                 !< Tile level maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: rg                 !< Tile level growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: nbp                !< Tile level net biome productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: dstcemls1          !< carbon emission losses due to disturbance (fire at present) from vegetation
    real, dimension(ilg), intent(out) :: litrfall           !< total litter fall (from leaves, stem, and root) due to all causes (mortality, turnover, and disturbance)
    real, dimension(ilg), intent(out) :: humiftrs           !< Transfer of humidified litter from litter to soil C pool (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: lucemcom           !< land use change (luc) related combustion emission losses, u-mol co2/m2.sec
    real, dimension(ilg), intent(out) :: lucltrin           !< luc related inputs to litter pool, u-mol co2/m2.sec
    real, dimension(ilg), intent(out) :: lucsocin           !< luc related inputs to soil c pool, u-mol co2/m2.sec
    real, dimension(ilg), intent(out) :: dstcemls3          !< carbon emission losses due to disturbance (fire at present) from litter pool
    real, dimension(ilg), intent(out) :: rms                !< Tile level stem maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: rmr                !< Tile level root maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: litres             !< Litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: socres             !< Soil carbon respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: rmsveg         !< Maintenance respiration for stem for the CTEM pfts in u mol co2/m2. sec
    real, dimension(ilg,icc), intent(out) :: rmrveg         !< Maintenance respiration for root for the CTEM pfts in u mol co2/m2. sec
    real, dimension(ilg,icc), intent(out) :: rmlveg         !< Leaf maintenance respiration per PFT (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: gppveg         !< Gross primary productivity per PFT (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: rgveg          !< PFT level growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,iccp1), intent(out) :: nepveg       !<
    real, dimension(ilg,iccp1), intent(out) :: nbpveg       !<
    real, dimension(ilg,iccp2,ignd), intent(out) :: ltresveg     !< fluxes for each pft: litter respiration for each pft + bare fraction
    real, dimension(ilg,iccp2,ignd), intent(out) :: scresveg     !< soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
    real, dimension(ilg,iccp1), intent(out) :: hetrsveg     !< Vegetation averaged litter and soil C respiration rates (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,iccp2,ignd), intent(out) :: humtrsvg     !< transfer of humidified litter from litter to soil c pool per PFT. (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: autoresveg     !< PFT level autotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: litrfallveg    !<
    real, dimension(ilg,icc), intent(out) :: roottemp       !< root temperature, k
    real, dimension(ilg,icc), intent(out) :: emit_co2       !< carbon dioxide emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_co        !< carbon monoxide emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_ch4       !< methane emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_nmhc      !< non-methane hydrocarbons emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_h2        !< hydrogen gas emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_nox       !< nitrogen oxides emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_n2o       !< nitrous oxide emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_pm25      !< particulate matter less than 2.5 um in diameter emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_tpm       !< total particulate matter emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_tc        !< total carbon emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_oc        !< organic carbon emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_bc        !< black carbon emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: bterm_veg      !< biomass term for fire probabilty calc
    real, dimension(ilg), intent(out) :: lterm              !< lightning term for fire probabilty calc
    real, dimension(ilg,icc), intent(out) :: mterm_veg      !< moisture term for fire probabilty calc
    real, dimension(ilg), intent(out) :: ch4WetSpec            !<
    real, dimension(ilg), intent(out) :: wetfdyn            !<
    real, dimension(ilg), intent(out) :: ch4WetDyn            !<
    real, dimension(ilg,icc), intent(out) :: colrate        !< colonization rate (1/day)
    real, dimension(ilg,icc), intent(out) :: mortrate       !< mortality rate
    real, dimension(ilg,icc), intent(out) :: afrleaf        !< allocation fraction for leaves
    real, dimension(ilg,icc), intent(out) :: afrstem        !< allocation fraction for stem
    real, dimension(ilg,icc), intent(out) :: afrroot        !< allocation fraction for root
    real, dimension(ilg,icc), intent(out) :: wtstatus       !< soil water status used for calculating allocation fractions
    real, dimension(ilg,icc), intent(out) :: ltstatus       !< light status used for calculating allocation fractions
    real, dimension(ilg), intent(out) :: burnfrac           !< areal :: fraction burned due to fire for every grid cell (%)
    real, dimension(ilg,icc), intent(out) :: leaflitr       !< leaf litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$). this leaf litter does not
    !< include litter generated due to mortality/fire
    real, dimension(ilg,icc), intent(out) :: smfunc_veg     !< soil moisture dependence on fire spread rate
    real, dimension(ilg,icc), intent(out) :: tltrleaf       !< total leaf litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: tltrstem       !< total stem litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: tltrroot       !< total root litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,ican), intent(out) :: paicgat       !<
    real, dimension(ilg,ican), intent(out) :: slaicgat      !<
    real, dimension(ilg,icc), intent(out) :: vgbiomas_veg   !<
    real, dimension(ilg), intent(out) :: armoss             !< autotrophic respiration of moss (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: nppmoss            !< net primary production of moss (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: reprocost(ilg,icc) !< Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: blfltrdt(ilg,icc)  !< brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(out) :: glfltrdt(ilg,icc)  !< green leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(out) :: glcaemls(ilg,icc)  !< green leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, intent(out) :: blcaemls(ilg,icc)  !< brown leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, intent(out) :: rtcaemls(ilg,icc)  !< root carbon emission disturbance losses, \f$kg c/m^2\f$
    real, intent(out) :: stcaemls(ilg,icc)  !< stem carbon emission disturbance losses, \f$kg c/m^2\f$
    real, intent(out) :: ltrcemls(ilg,icc)  !< litter carbon emission disturbance losses, \f$kg c/m^2\f$
    real, intent(out) :: ntchlveg(ilg,icc)  !< fluxes for each pft: Net change in leaf biomass, u-mol CO2/m2.sec
    real, intent(out) :: ntchsveg(ilg,icc)  !< fluxes for each pft: Net change in stem biomass, u-mol CO2/m2.sec
    real, intent(out) :: ntchrveg(ilg,icc)  !< fluxes for each pft: Net change in root biomass,
    !! the net change is the difference between allocation and
    !! autotrophic respiratory fluxes, u-mol CO2/m2.sec

    ! ---------------------------------------------
    ! Local variables:

    integer :: i, j, k, n, m
    integer :: sort(icc)
    real :: yesfrac_comp(ilg,icc) !<
    real :: fc(ilg)  !< Fraction of grid cell that is covered by canopy
    real :: fg(ilg)  !< Fraction of grid cell that is bare ground.
    real :: pglfmass(ilg,icc)  !< Prior timestep green leaf mass for each of the CTEM pfts, \f$(kg C/m^2)\f$
    real :: pblfmass(ilg,icc)  !< Prior timestep brown leaf mass for each of the CTEM pfts, \f$(kg C/m^2)\f$
    real :: pstemass(ilg,icc)  !< Prior timestep stem mass for each of the CTEM pfts, \f$(kg C/m^2)\f$
    real :: protmass(ilg,icc)  !< Prior timestep root mass for each of the CTEM pfts, \f$(kg C/m^2)\f$
    real :: plitmass(ilg,iccp2)!< Prior timestep litter mass for each of the CTEM pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real :: psocmass(ilg,iccp2)!< Prior timestep soil carbon mass for each of the CTEM pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real :: pvgbioms(ilg)      !< Prior timestep
    real :: pgavltms(ilg)      !< Prior timestep
    real :: pgavscms(ilg)      !< Prior timestep
    real :: pheanveg(ilg,icc) !<
    real :: rootlitr(ilg,icc) !< root litter \f$(kg C/m^2)\f$
    real :: stemlitr(ilg,icc) !< stem litter \f$(kg C/m^2)\f$
    real :: stemltrm(ilg,icc) !<
    real :: rootltrm(ilg,icc) !<
    real :: glealtrm(ilg,icc) !<
    real :: stemltdt(ilg,icc)  !<
    real :: rootltdt(ilg,icc)  !<
    real :: dscemlv1(ilg,icc)  !<
    real :: dscemlv2(ilg,icc)  !<
    real :: add2allo(ilg,icc)  !<
    real :: repro_cost_g(ilg)  !< Tile-level cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: litresmoss(ilg)    !< moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: socres_peat(ilg)   !< heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: resoxic(ilg)       !< oxic respiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: resanoxic(ilg)     !< anoxic respiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: litrfallmoss(ilg)  !< moss litter fall (kgC/m2/timestep)
    real :: ltrestepmoss(ilg)  !< litter respiration from moss (kgC/m2/timestep)
    real :: humstepmoss(ilg)   !< moss humification (kgC/m2/timestep)
    real :: pCmossmas(ilg)     !< moss biomass C at the previous time step (kgC/m2)
    real :: plitrmsmoss(ilg)   !< moss litter C at the previous time step (kgC/m2)
    real :: nppmosstep(ilg)    !< moss npp (kgC/m2/timestep)
    real :: socrestep(ilg)     !< heterotrophic respiration from soil (kgC/m2/timestep)
    real :: hutrstep_g(ilg)    !< grid sum of humification from vascualr litter (kgC/m2/timestep)

    real :: rmsTracer(ilg,icc)   !< Tracer maintenance respiration for stem for the CTEM pfts (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: rmrTracer(ilg,icc)   !< Tracer maintenance respiration for root for the CTEM pfts both (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: ltResTracer(ilg,iccp2,ignd)  !< Tracer fluxes for each pft: litter respiration for each pft + bare fraction
    real :: sCResTracer(ilg,iccp2,ignd)  !< Tracer soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
    real :: litResMossTracer(ilg)    !< Tracer moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: soCResPeatTracer(ilg)    !< Tracer heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: tracerNPP(ilg,icc) !< tracer NPP for individual pfts, (\f$tracer C units m^{-2} s^{-1}\f$)
    real :: tracerLeafLitr(ilg,icc)  !< Tracer leaf litter generated by normal turnover, cold
    !! and drought stress,and leaf fall/harvest, \f$tracer C units/m^2\f$
    real :: tracerStemLitr(ilg,icc) !< Tracer stem litter \f$tracer C units/m^2\f$
    real :: tracerRootLitr(ilg,icc) !< Tracer root litter \f$tracer C units/m^2\f$
    real :: tracerReproCost(ilg,icc)  !< Tracer cost of making reproductive tissues, only non-zero when
    !! NPP is positive (\f$tracer C units m^{-2} s^{-1}\f$)
    real :: tracerStemMort(ilg,icc) !< Tracer stem litter from mortality \f$tracer C units/m^2\f$
    real :: tracerRootMort(ilg,icc) !< Tracer stem litter from mortality \f$tracer C units/m^2\f$
    real :: tracerGLeafMort(ilg,icc) !< Tracer stem litter from mortality \f$tracer C units/m^2\f$
    real :: tracerValue(ilg)       !< Tracer CO2 value updated from the prepTracer subroutine, units vary (simple: ppm, 14C \f$\Delta ^{14}C\f$)
    real :: tracerRML(ilg,icc)     !< Tracer leaf maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: tracerGPP(ilg,icc) !< Tracer GPP (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    ! ---------------------------------------------------------------------------
    ! Begin calculations

    !> Generate the sort index for correspondence between CTEM pfts and the
    !>  values in the parameter vectors
    sort = genSortIndex()

    !> Set up the tracer for today
    if (useTracer > 0) call prepTracer(il1,il2,ilg,tracerCO2, & ! In
        &                              tracerValue) ! Out

    if (PFTCompetition) then

      !> Calculate bioclimatic parameters for estimating pfts existence
      call  bioclim(iday, ta, precip, netrad, &
                    1, il2, ilg, leapnow, &
                    tcurm, srpcuryr, dftcuryr, inibioclim, &
                    tmonth, anpcpcur, anpecur, gdd5cur, &
                    surmncur, defmncur, srplscur, defctcur, &
                    twarmm, tcoldm, gdd5, aridity, &
                    srplsmon, defctmon, anndefct, annsrpls, &
                    annpcp, dry_season_length)

      if (inibioclim) then

        !> If first day of year then based on updated bioclimatic parameters
        !! find if pfts can exist or not.
        !! If .not. inibioclim then it is the first year of a run that you do not have the
        !! climatological means already in the CTM file. After one
        !! year inibioclim is set to true and the climatological means
        !! are used from the first year.
        !!
        call existence(iday, 1, il2, ilg, &
                       sort, twarmm, tcoldm, &
                       gdd5, aridity, srplsmon, defctmon, &
                       anndefct, annsrpls, annpcp, pftexist, &
                       dry_season_length)

        !> Call competition subroutine which on the basis of previous day's
        !! NPP estimates changes in fractional coverage of pfts
        call competition(iday, 1, il2, ilg, nppveg, dofire, leapnow, useTracer, & ! In
                         pftexist, geremort, intrmort, pgleafmass, rmatctem, & ! In
                         grclarea, lambda, burnvegf, sort, pstemmass, & ! In
                         gleafmas, bleafmas, stemmass, rootmass, & ! In/Out
                         litrmass, soilcmas, fcancmx, fcanmx, & ! In/Out
                         tracerGLeafMass, tracerBLeafMass, tracerStemMass, tracerRootMass, & ! In/Out
                         tracerLitrMass, tracerSoilCMass, & ! In/Out
                         vgbiomas, gavgltms, gavgscms, bmasveg, & ! In/Out
                         add2allo, colrate, mortrate) ! Out

      end if ! inibioclim
    end if  ! if (PFTCompetition)
    !>
    !! If landuse is on, then implelement luc, change fractional coverages,
    !! move biomasses around, and estimate luc related combustion emission losses.
    !!
    if (lnduseon) then

      do j = 1,icc
        do i = il1,il2
          yesfrac_comp(i,j) = fcancmx(i,j)
        end do
      end do

      call luc(il1, il2, ilg, PFTCompetition, leapnow, useTracer, & ! In
               grclarea, iday, todfrac, yesfrac_comp, .true., & ! In
               pfcancmx, nfcancmx, &! In/Out
               gleafmas, bleafmas, stemmass, rootmass, &! In/Out
               litrmass, soilcmas, vgbiomas, gavgltms, &! In/Out
               gavgscms, fcancmx, fcanmx, tracerLitrMass, tracerSoilCMass, & ! In/Out
               tracerGLeafMass, tracerBLeafMass, tracerStemMass, tracerRootMass, & ! In / Out
               lucemcom, lucltrin, lucsocin) ! Out
    else
      lucemcom = 0.
      lucltrin = 0.
      lucsocin = 0.
    end if ! lnduseon

    !> Store green and brown leaf, stem, and root biomass, and litter and
    !! soil c pool mass in arrays. knowing initial sizes of all pools and
    !! final sizes at the end of this subroutine, we check for conservation of mass.
    plitmass = 0.0
    psocmass = 0.0
    do j = 1,icc
      do i = il1,il2
        pglfmass(i,j) = gleafmas(i,j)    !< green leaf mass from last time step
        pblfmass(i,j) = bleafmas(i,j)    !< brown leaf mass from last time step
        pstemass(i,j) = stemmass(i,j)    !< stem mass from last time step
        protmass(i,j) = rootmass(i,j)    !< root mass from last time step
        ! COMBAK PERLAY
        plitmass(i,j) = plitmass(i,j) + litrmass(i,j)    ! litter mass from last time step
        psocmass(i,j) = psocmass(i,j) + soilcmas(i,j)    ! soil c mass from last time step
        ! do k = 1,ignd ! FLAG at this stage keep as per pft and per tile. JM Feb8 2016.
        !   plitmass(i,j)=plitmass(i,j) + litrmass(i,j,k)    ! litter mass from last time step
        !   psocmass(i,j)=psocmass(i,j) + soilcmas(i,j,k)    ! soil c mass from last time step
        !   plitmasspl(i,j,k)= litrmass(i,j,k)
        !   psocmasspl(i,j,k)= soilcmas(i,j,k)
        ! end do
        ! COMBAK PERLAY
      end do ! loop 140
    end do ! loop 130
    !
    do i = il1,il2
      pvgbioms(i) = vgbiomas(i)          !< vegetation biomass from last time step
      pgavltms(i) = gavgltms(i)          !< litter mass from last time step
      pgavscms(i) = gavgscms(i)          !< soil c mass from last time step
      do j = iccp1,iccp2 ! do over the bare fraction and the LUC pool
        ! COMBAK PERLAY
        plitmass(i,j) = plitmass(i,j) + litrmass(i,j)  ! litter mass over bare fraction
        psocmass(i,j) = psocmass(i,j) + soilcmas(i,j)  ! soil c mass over bare fraction
        ! do k = 1,ignd ! FLAG at this stage keep as per pft and per tile. JM Feb8 2016.
        !   plitmass(i,j)=plitmass(i,j) + litrmass(i,j,k)  ! litter mass over bare fraction
        !   psocmass(i,j)=psocmass(i,j) + soilcmas(i,j,k)  ! soil c mass over bare fraction
        ! end do
        ! COMBAK PERLAY
      end do

      pCmossmas(i)  = Cmossmas(i)
      plitrmsmoss(i) = litrmsmoss(i)
    end do ! loop 145

    ! Find the canopy covered fraction and the bare fraction of the tiles:
    do i = il1,il2
      fc(i) = sum(fcancmx(i,:))
      fg(i) = 1.0 - fc(i)
    end do

    !     ------------------------------------------------------------------
    !> Initialization ends

    !> Autotrophic respiration
    !!
    !! Leaf respiration is calculated in phtsyn subroutine, while stem
    !! and root maintenance respiration are calculated here. We use air
    !! temperature as a surrogate for  stem temperature
    !!
    !! Find stem and root maintenance respiration in umol co2/m2/sec
    !! If the tracer is being calculated then we also determine the
    !! tracer flux.
    call mainres(fcancmx, fc, stemmass, rootmass, & ! In
                 il1, il2, ilg, leapnow, & ! In
                 ta, tbar, rmatctem, sort, isand, & ! In
                 useTracer, tracerStemMass, tracerRootMass, & ! In
                 rmsveg, rmrveg, roottemp, & ! Out
                 rmsTracer, rmrTracer) ! Out

    !> Calculate NPP (net primary productivity), difference between GPP and
    !! autotrophic respriation, for each pft and tile. Also sets the net
    !! photosynthesis to be used by phenology and determines the total
    !! autotrophic respiration fluxes at the PFT and tile-levels. Calculates
    !! values for both upland and peatland sites.
    call calcNPP(il1, il2, ilg, ancgveg, lfstatus, rmlcgveg, & ! In
                 slai, ailcg, sort, rmsveg, rmrveg, ipeatland, & ! In
                 anmoss, fcancmx, rmlmoss, gppmoss, useTracer, & ! In
                 gleafmas, tracerGLeafMass, tracerValue, rmsTracer, rmrTracer, & ! In
                 pheanveg, rmlveg, rml, rms, rmr, rm, rg, npp, gpp, & ! Out
                 autores, autoresveg, rgveg, nppmoss, armoss, & ! Out
                 gppveg, nppmosstep, nppveg, tracerNPP, tracerRML, tracerGPP) ! Out

    !! Find heterotrophic respiration rates (umol co2/m2/sec)
    ! If tracer is being used then calculate the tracer fluxes.

    ! CAUTION says Vivek
    ! Note that ipeatland is passed to hetresv to calculate psi (matric potential)
    ! in a different way if ipeatland(i)=1. psi is used to calculated soil moisture
    ! dependence scalars for litter and soil cabron respiration in hetresv. Yet,
    ! soil carbon respiration for gridcells/tiles which are peatlands is overwritten
    ! below in loop 380 by socres_peat(i) making passing of ipeatland to hetresv and
    ! calculation of psi in a different way useless.

    call heterotrophicRespiration(il1, il2, ilg, ipeatland, fcancmx, fc, & ! In
                                  litrmass, soilcmas, delzw, thpor, tbar, & ! In
                                  psisat, thliq, sort, bi, isand, thice, & ! In
                                  fg, litrmsmoss, peatdep, wtable, zbotw, roottemp, & ! In
                                  useTracer, tracerLitrMass, tracerSoilCMass, tracerMossLitrMass, & ! In
                                  ltresveg, scresveg, litresmoss, socres_peat, & ! Out
                                  resoxic, resanoxic , & ! Out
                                  ltResTracer, sCResTracer, litResMossTracer, soCResPeatTracer) ! Out

    !> Find vegetation and tile averaged litter and soil C respiration rates
    !! using values from canopy over ground and canopy over snow subareas.
    !! Also adds the moss and peat soil respiration to the tile level quantities.
    !! Next the litter and soil C pools are updated based on litter and soil C respiration rates.
    !! The humidified litter is then transferred to the soil C pool.
    !! Soil respiration is estimated as the sum of heterotrophic respiration and root maintenance respiration.
    !! For peatlands, we additionally add moss values to the grid (litter respiration
    !! and moss root respiration).
    call updatePoolsHetResp(il1, il2, ilg, fcancmx, ltresveg, scresveg, & ! In
                            ipeatland, fg, litresmoss, socres_peat, & ! In
                            sort, spinfast, rmrveg, rmr, leapnow, & ! In
                            useTracer, ltResTracer, sCResTracer, litResMossTracer, soCResPeatTracer, & ! In
                            litrmass, soilcmas, Cmossmas, & ! In / Out
                            tracerLitrMass, tracerSoilCMass, tracerMossCMass, & ! In / Out
                            hetrsveg, litres, socres, hetrores, humtrsvg, soilresp, & ! Out
                            humiftrs, hutrstep_g, litrfallmoss, ltrestepmoss, &
                            humstepmoss, socrestep) ! Out

    !> Calculate NEP (net ecosystem productivity), difference between NPP and
    !! heterotrophic respriation, for each pft and tile
    call calcNEP(il1, il2, ilg, nppveg, hetrsveg, fg, npp, hetrores, & ! In
                 nep, nepveg) ! Out


    if (doMethane) then 
      
      !> Find CH4 wetland area (if not prescribed) and emissions:
      call  wetland_methane(hetrores, il1, il2, ilg, & ! In
                          wetfrac, thliq, currlat, sand, &   ! In
                          slopefrac, ta, & ! In
                          ch4WetSpec, wetfdyn, ch4WetDyn) ! Out

      !> Calculate the methane that is oxidized by the soil sink
      call soil_ch4uptake(il1, il2, ilmos, ilg, tbar, & ! In
                        bi, thliq, thice, psisat, & ! In
                        fcanmx, wetfdyn, wetfrac, & ! In
                        isand, ch4conc, thpor, & ! In
                        ch4soills) ! Out
    end if 
    
    !> Estimate allocation fractions for leaf, stem, and root components.
    call allocate(lfstatus, thliq, ailcg, & ! In
                  il1, il2, ilg, & ! In
                  rmatctem, gleafmas, stemmass, & ! In
                  rootmass, sort, fcancmx, & ! In
                  isand, THFC, THLW, & ! In
                  afrleaf, afrstem, afrroot, & ! Out
                  wtstatus, ltstatus) ! Out

    call updatePoolsAllocateRepro(il1, il2, ilg, sort, ailcg, lfstatus, nppveg, PFTCompetition, & ! In
                                  pftexist, gppveg, rmsveg, rmrveg, rmlveg, fcancmx, & ! In
                                  useTracer, tracerNPP, rmsTracer, rmrTracer, tracerRML, tracerGPP, & ! In
                                  lambda, afrleaf, afrstem, afrroot, gleafmas, stemmass, & ! In /Out
                                  rootmass, bleafmas, tracerGLeafMass, tracerStemMass, & ! In/Out
                                  tracerRootMass, tracerBLeafMass, & ! In/Out
                                  reprocost, ntchlveg, ntchsveg, ntchrveg, & ! Out
                                  repro_cost_g, tracerReproCost) ! Out

    !! The phenology subroutine determines leaf status for each pft and calculates leaf litter.
    !! the phenology subroutine uses soil temperature (tbar) and root temperature. however,
    !! since CTEM doesn't make the distinction between canopy over ground, and canopy over
    !! snow sub-areas for phenology purposes (for example, leaf onset is not assumed to occur
    !! at different times over these sub-areas) we use average soil and root temperature in
    !! the phenology subroutine.
    call phenolgy(il1, il2, ilg, leapnow, tbar, thice, &! In
                  thliq, THLW, THFC, ta, &! In
                  pheanveg, iday, radj, roottemp, &! In
                  rmatctem, stemmass, rootmass, sort, &! In
                  fcancmx, isand, useTracer, &! In
                  lfstatus, pandays, colddays, gleafmas, &
                  bleafmas, tracerGLeafMass, tracerBLeafMass, & ! In/Out
                  flhrloss, leaflitr, tracerLeafLitr) ! Out

    !> While leaf litter is calculated in the phenology subroutine, stem
    !! and root turnover is calculated in the turnoverStemRoot subroutine.
    call turnoverStemRoot(stemmass, rootmass, lfstatus, ailcg, & ! In
                          il1, il2, ilg, leapnow, useTracer, &! In
                          sort, fcancmx, tracerStemMass, tracerRootMass, &! In
                          stmhrlos, rothrlos, & ! In/Out
                          stemlitr, rootlitr, tracerStemLitr, tracerRootLitr) ! Out

    !> Update green leaf biomass for trees and crops, brown leaf biomass for grasses,
    !! stem and root biomass for litter deductions, and update litter pool with leaf
    !! litter calculated in the phenology subroutine and stem and root litter
    !! calculated in the turnoverStemRoot subroutine. Also add the reproduction
    !!  carbon directly to the litter pool. We only add to non-perennially frozen soil
    !! layers so first check which layers are unfrozen and then do the allotment
    !! appropriately. For defining which layers are frozen, we use the active layer depth.
    call updatePoolsTurnover(il1, il2, ilg, reprocost, rmatctem, useTracer, tracerReproCost, & ! In
                             stemmass, rootmass, litrmass, rootlitr, & ! In/Out
                             gleafmas, bleafmas, leaflitr, stemlitr, & ! In/Out
                             tracerGLeafMass, tracerBLeafMass, tracerLitrMass, & ! In/Out
                             tracerRootMass, tracerStemMass, & ! In/Out
                             tracerLeafLitr, tracerStemLitr, tracerRootLitr) ! In/Out

    !> Call the mortality subroutine which calculates mortality due to reduced growth and aging.
    !! Exogenous mortality due to fire and other disturbances and the subsequent litter
    !! that is generated is calculated in the disturb subroutine.
    !!
    !! Set maxage >0 in classicParams.f90 to switch on mortality due to age and
    !! reduced growth. Mortality is linked to the competition parameterization and generates bare fraction.
    call mortalty(stemmass, rootmass, ailcg, gleafmas, & ! In
                  bleafmas, il1, il2, ilg, & ! In
                  leapnow, iday, sort, fcancmx, & ! In
                  useTracer, tracerStemMass, tracerRootMass, tracerGLeafMass, & ! In
                  lystmmas, lyrotmas, tymaxlai, grwtheff, & ! In/Out
                  stemltrm, rootltrm, glealtrm, geremort, & ! Out
                  intrmort, tracerStemMort, tracerRootMort, tracerGLeafMort) ! Out

    !> Update leaf, stem, and root biomass pools to take into loss due to mortality, and put the
    !! litter into the litter pool. the mortality for green grasses doesn't generate litter, instead they turn brown.
    call updatePoolsMortality(il1, il2, ilg, stemltrm, rootltrm, useTracer, & ! In
                              rmatctem, tracerStemMort, tracerRootMort, tracerGLeafMort, & ! In
                              stemmass, rootmass, litrmass, & ! In/Out
                              glealtrm, gleafmas, bleafmas, tracerLitrMass, & ! In/Out
                              tracerStemMass, tracerRootMass, tracerGLeafMass, tracerBLeafMass) ! In/Out

    !> Call the disturbance subroutine which calculates mortality due to fire and other disturbances.
    !> the primary output from from disturbance subroutine is litter generated, c emissions due to fire
    !! and area burned, which may be used to estimate change in fractional coverages.
    !!
    !! Disturbance is spatial and requires area of gcm grid cell and areas of different pfts present in
    !! a given grid cell. however, when ctem is operated at a point scale then it is assumed that the
    !! spatial scale is 1 hectare = 10, 000 m2. the disturbance subroutine may be stopped from simulating
    !! any fire by specifying fire extingushing probability equal to 1.
    call disturb(thliq, THLW, THFC, uwind, useTracer, & ! In
                 vwind, lightng, fcancmx, isand, & ! In
                 rmatctem, ilg, il1, il2, sort, & ! In
                 grclarea, thice, popdin, lucemcom, & ! In
                 dofire, currlat, iday, fsnow, & ! In
                 stemmass, rootmass, gleafmas, bleafmas, litrmass, & ! In/Out
                 tracerStemMass, tracerRootMass, tracerGLeafMass, tracerBLeafMass, tracerLitrMass, & ! In/Ou
                 stemltdt, rootltdt, glfltrdt, blfltrdt, & ! Out (Primary)
                 glcaemls, rtcaemls, stcaemls, & ! Out (Primary)
                 blcaemls, ltrcemls, burnfrac, & ! Out (Primary)
                 pstemmass, pgleafmass, emit_co2, emit_ch4, & ! Out (Primary)
                 emit_co, emit_nmhc, emit_h2, emit_nox, & ! Out (Secondary)
                 emit_n2o, emit_pm25, emit_tpm, emit_tc, & ! Out (Secondary)
                 emit_oc, emit_bc, burnvegf, bterm_veg, & ! Out (Secondary)
                 mterm_veg, lterm, smfunc_veg) ! Out (Secondary)

    !> Calculate NBP (net biome production) for each pft by taking into account
    !! C emission losses. The disturbance routine produces emissions due to fire
    !! and while the land use change subroutine calculates emissions due to LUC.
    !! The LUC related combustion flux is assumed to be spread uniformly over the
    !! tile as it is no longer associated with any one PFT. To calculate the NBP
    !! we do not subtract LUC emissions from the PFT-level NBP but we do subtract
    !! it from the per tile NBP.
    call calcNBP(il1, il2, ilg, deltat, nepveg, fcancmx, & ! In
                 lucemcom, ltresveg, scresveg, nep, & ! In
                 glcaemls, blcaemls, stcaemls, rtcaemls, ltrcemls, & ! In/Out
                 nbpveg, dstcemls1, dstcemls3, nbp) ! Out

    !> Allow cryoturbation and bioturbation to move the soil C between
    !! layers. Since this is neither consuming nor adding C, this does not
    !! affect our C balance in balcar. There is also an internal C balance check.
    ! COMBAK PERLAY
    ! call turbation(il1, il2, delzw, zbotw, isand, maxAnnualActLyr, spinfast, &! In
    !               litrmass, soilcmas)! In/Out
    !
    ! if (useTracer > 0) call turbation(il1, il2, delzw, zbotw, isand, maxAnnualActLyr, spinfast, &! In
    !                             tracerLitrMass, tracerSoilCMass) ! In/Out
    ! COMBAK PERLAY

    !> Prepare for the carbon balance check. Calculate total litter fall from each
    !! component (leaves, stem, and root) from all causes (normal turnover, drought
    !! and cold stress for leaves, mortality, and disturbance), calculate grid-average
    !! vegetation biomass, litter mass, and soil carbon mass, and litter fall rate.
    !! Also add the bare ground values to the grid-average. If a peatland, we assume no bareground and
    !! add the moss values instead. Note: peatland soil C is not aggregated from plants but updated
    !! by humification and respiration from the previous stored value
    call prepBalanceC(il1, il2, ilg, fcancmx, glealtrm, glfltrdt, &  ! In
                      blfltrdt, stemltrm, stemltdt, rootltrm, rootltdt, & ! In
                      ipeatland, nppmosstep, pgavscms, humstepmoss, & ! In
                      ltrestepmoss, stemlitr, rootlitr, rootmass, & ! In
                      litrmass, soilCmas, hutrstep_g, stemmass, bleafmas, & ! In
                      gleafmas, socrestep, fg, litrfallmoss, & ! In
                      leaflitr, Cmossmas, litrmsmoss, & ! In/Out
                      tltrleaf, tltrstem, tltrroot, & ! Out
                      vgbiomas, litrfall, gavgltms, litrfallveg, &! Out
                      gavgscms, vgbiomas_veg) ! Out

    !> At this stage we have all required fluxes in u-mol co2/m2.sec and initial (loop 140 and 145)
    !! and updated sizes of all pools (in \f$(kg C/m^2)\f$). Now we call the balcar subroutine and make sure
    !! that C in leaves, stem, root, litter and soil C pool balances within a certain tolerance.
    if (spinfast == 1) then
      call  balcar(gleafmas, stemmass, rootmass, bleafmas, &
                   litrmass, soilcmas, ntchlveg, ntchsveg, &
                   ntchrveg, tltrleaf, tltrstem, tltrroot, &
                   glcaemls, blcaemls, stcaemls, rtcaemls, &
                   ltrcemls, ltresveg, scresveg, humtrsvg, &
                   pglfmass, pblfmass, pstemass, protmass, &
                   plitmass, psocmass, vgbiomas, reprocost, &
                   pvgbioms, gavgltms, pgavltms, gavgscms, &
                   pgavscms, dstcemls3, repro_cost_g, &
                   autores, hetrores, gpp, &
                   litres, socres, dstcemls1, &
                   litrfall, humiftrs, &
                   il1, il2, ilg, &
                   ipeatland, Cmossmas, pCmossmas, &
                   nppmosstep, litrfallmoss, litrmsmoss, &
                   plitrmsmoss, ltrestepmoss, humstepmoss)

      !> Check for mass balance for the tracer if the tracer is being used and
      !! doTracerBalance is true.
      if (useTracer > 0 .and. doTracerBalance) call checkTracerBalance(il1, il2)

    end if

    !> Finally find vegetation structural attributes which can be passed
    !! to the land surface scheme using leaf, stem, and root biomass.
    !>
    call allometry(gleafmas, bleafmas, stemmass, rootmass, & ! In
                   il1, il2, ilg, zbotw, & ! In
                   soildpth, fcancmx, & ! In
                   ipeatland, maxAnnualActLyr, & ! In
                   ailcg, ailcb, ailc, zolnc, & ! Out
                   rmatc, rmatctem, slai, bmasveg, & ! Out
                   cmasvegc, veghght, rootdpth, alvisc, & ! Out
                   alnirc, paicgat, slaicgat) ! Out

    !> Calculation of gavglai is moved from loop 1100 to here since ailcg is updated by allometry
    gavglai (:) = 0.0
    do j = 1,icc
      do i = il1,il2
        gavglai(i) = gavglai(i) + fcancmx(i,j) * ailcg(i,j)
      end do
    end do

    !> At the end of the day, find the depth of the peat, update the degree days for moss photosynthesis and the peat bottom layer depth
    do i = il1,il2
      if (ipeatland(i) > 0) peatdep(i) = peatDepth(gavgscms(i))
    end do

    call peatDayEnd(il2)

    return

  end subroutine ctem
  !! @}
  ! ----------------------------------------------------------------------------
  !> \ingroup ctemdriver_calcnpp
  !! @{
  !> Calculate NPP (net primary productivity), difference between GPP and
  !! autotrophic respriation, for each pft and tile. Also sets the net
  !! photosynthesis to be used by phenology and determines the total
  !! autotrophic respiration fluxes at the PFT and tile-levels. Calculates
  !! values for both upland and peatland sites.
  !! @author V. Arora, J. Melton
  subroutine calcNPP (il1, il2, ilg, ancgveg, lfstatus, rmlcgveg, & ! In
                      slai, ailcg, sort, rmsveg, rmrveg, ipeatland, & ! In
                      anmoss, fcancmx, rmlmoss, gppmoss, useTracer, & ! In
                      gleafmas, tracerGLeafMass, tracerValue, rmsTracer, rmrTracer, & ! In
                      pheanveg, rmlveg, rml, rms, rmr, rm, rg, npp, gpp, & ! Out
                      autores, autoresveg, rgveg, nppmoss, armoss, & ! Out
                      gppveg, nppmosstep, nppveg, tracerNPP, tracerRML, tracerGPP) ! Out

    use classicParams, only : icc, iccp1, kn, zero, grescoefmoss, deltat, &
                              grescoef
    use autotrophicRespiration, only : growthRespiration

    implicit none

    ! arguments
    integer, intent(in) :: il1             !< il1=1
    integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg
    integer, intent(in) :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    real, intent(in) :: ancgveg(:,:)      !< net photosynthetic rate for CTEM's pfts
    integer, intent(in) :: lfstatus(:,:)  !< leaf phenology status
    real, intent(in) :: rmlcgveg(:,:)     !< leaf respiration rate for CTEM's pfts
    real, intent(in) :: slai(:,:)         !< storage/imaginary lai for phenology purposes
    real, intent(in) :: ailcg(:,:)        !< green lai for ctem's 9 pfts
    integer, intent(in) :: sort(:)        !< index for correspondence between biogeochem pfts and the number of values in parameters vectors in run params file.
    real, intent(in) :: rmsveg(:,:)       !< Maintenance respiration for stem for the CTEM pfts in u mol co2/m2. sec
    real, intent(in) :: rmrveg(:,:)       !< Maintenance respiration for root for the CTEM pfts in u mol co2/m2. sec
    integer, intent(in) :: ipeatland(:)   !< Peatland flag: 0 = not a peatland, 1 = bog, 2 = fen
    real, intent(in) :: anmoss(:)         !< moss net photoysnthesis -daily averaged C fluxes rates (umol/m2/s)
    real, intent(in) :: fcancmx(:,:)      !< max. fractional coverage of CTEM's pfts, but this can be
    !! modified by land-use change,and competition between pfts
    real, intent(in) :: rmlmoss(:)        !< moss maintainance respiration -daily averaged C fluxes rates (umol/m2/s)
    real, intent(in) :: gppmoss(:)        !< moss GPP -daily averaged C fluxes rates (umol/m2/s)
    real, intent(in) :: gleafmas(:,:)     !< green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: tracerGLeafMass(:,:) !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, intent(in) :: tracerValue(:)       !< Tracer CO2 value read in from tracerCO2File, units vary (simple: ppm, 14C \f$\Delta ^{14}C\f$)
    real, intent(in) :: rmsTracer(:,:)   !< Tracer maintenance respiration for stem for the CTEM pfts (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: rmrTracer(:,:)   !< Tracer maintenance respiration for root for the CTEM pfts both (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real, intent(out) :: pheanveg(ilg,icc) !<
    real, intent(out) :: rmlveg(ilg,icc)   !< Leaf maintenance respiration per PFT (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: gppveg(ilg,icc)   !< Gross primary productivity per PFT (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rml(ilg)          !< Tile level leaf maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rms(ilg)         !< Tile level stem maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rmr(ilg)         !< Tile level Root maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rm(ilg)          !< Tile level maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rg(ilg)          !< Tile level growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: npp(ilg)         !< Tile-level net primary productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: gpp(ilg)         !< Tile level gross primary productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: autores(ilg)     !< Tile level autotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: autoresveg(ilg,icc)  !< PFT level autotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rgveg(ilg,icc)   !< PFT level growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: nppmoss(ilg)     !< Net primary production of moss (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: armoss(ilg)      !< autotrophic respiration of moss (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: nppmosstep(ilg)  !< moss NPP (kgC/m2/timestep)
    real, intent(out) :: nppveg(ilg,icc)  !< NPP for individual pfts, (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: tracerNPP(ilg,icc) !< tracer NPP for individual pfts, (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: tracerRML(ilg,icc) !< Tracer leaf maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: tracerGPP(ilg,icc) !< Tracer GPP (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    ! Local
    integer :: j,i
    real :: anveg(ilg,icc)    !<
    real :: term                 !<
    real :: rgmoss(ilg)          !< moss growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: rmveg(ilg,icc)       !<
    real :: frac                 !< temp var.
    real :: tracerRG(ilg,icc)    !< Tracer growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: tracerRM             !<  Tracer total maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    !--------

    ! NOTE: This next bit is a little tricky. Remember ancgveg is the net photosynthesis
    ! so it is ancgveg = gpp - rmlcgveg. Here we assign the daily mean net photosynthesis
    ! (ancgveg) and mean rml (rmlveg) to the anveg and rmlveg variables
    ! and their sum to gpp veg since gpp = anveg + rml if we both have some coverage of the
    ! PFT (fcancmx > 1) and the leaves are not imaginary (lfstatus /= 4).
    ! If no real leaves are in existence we leave rml, anveg and gpp set to 0.
    ! However, if we have real leaves, but they are quite small we are going to
    ! reduce rml, but we need to use the original rml to find the gppveg otherwise
    ! the gpp will not be corrected properly. We store the original ancgveg for
    ! phenology to test if the leaves should be coming out.

    pheanveg = 0.
    rmlveg = 0.
    anveg = 0.
    gppveg = 0.
    do j = 1,icc
      do i = il1,il2

        if (fcancmx(i,j) > zero) then

          pheanveg(i,j) = ancgveg(i,j) ! to be used for phenology purposes

          if (lfstatus(i,j) /= 4) then ! real leaves so use values

            anveg(i,j) = ancgveg(i,j)
            rmlveg(i,j) = rmlcgveg(i,j)
            gppveg(i,j) = anveg(i,j) + rmlveg(i,j)

            if (slai(i,j) > ailcg(i,j)) then
              term = ((1.0 / kn(sort(j))) * (1.0 - exp( - kn(sort(j)) * ailcg(i,j))) &
                     / (1.0 / kn(sort(j))) * (1.0 - exp( - kn(sort(j)) * slai(i,j))))
              rmlveg(i,j) = rmlveg(i,j) * term
            end if
            ! else
            ! the leaves were imaginary so leave variables set to the initialized 0 value.
          end if
        end if
      end do ! loop 190
    end do ! loop 180

    !! Find total maintenance respiration values and net primary productivity.
    rml(:) = 0.
    rms(:) = 0.
    rmr(:) = 0.
    rm(:) = 0.
    rg(:) = 0.
    npp(:) = 0.
    gpp(:) = 0.
    tracerNPP(:,:) = 0.

    do j = 1,icc
      do i = il1,il2
        rmveg(i,j)  = rmlveg(i,j) + rmrveg(i,j) + rmsveg(i,j)
        nppveg(i,j) = gppveg(i,j) - rmveg(i,j)

        !> Now that we know maintenance respiration from leaf, stem, and root
        !! and gpp, we can find growth respiration for each vegetation type
        call growthRespiration(il1,il2,ilg,sort,useTracer,nppveg,tracerNPP,&
                              rgveg,tracerRG)

        nppveg(i,j) = nppveg(i,j) - rgveg(i,j)
        
        if (useTracer > 0) then
          !> Determine the NPP for the tracer. To find the value of rml for the tracer
          !! we need to scale the 'normal' rml value by the proportion of tracer to
          !! 12C. Because the calculation of rml does not include leaf mass (only fpar
          !! so by extentions LAI), there is no simple means to calculate tracerRML.
          if (gleafmas(i,j) > 0.) then
            frac = tracerGLeafMass(i,j) / gleafmas(i,j)
          else
            frac = 0.
          end if
          tracerRML(i,j) = rmlveg(i,j) * frac

          tracerRM = tracerRML(i,j) + rmrTracer(i,j) + rmsTracer(i,j)

          !> Now we can find the tracerGPP by multiplying the gppveg value by the
          !! tracerValue.
          tracerGPP(i,j) = gppveg(i,j) * tracerValue(i)

          !> And from that calculate the tracerNPP,which is what we are needing.
          tracerNPP(i,j) = tracerGPP(i,j) - tracerRM

          ! tracerRG = 0.
          ! if (tracerNPP(i,j) > 0.) tracerRG = grescoef(sort(j)) * tracerNPP(i,j)

          tracerNPP(i,j) = tracerNPP(i,j) - tracerRG(i,j)

        end if

        !> Calculate grid/tile-averaged rates of rm, rg, npp, and gpp
        rml(i) = rml(i) + fcancmx(i,j) * rmlveg(i,j)
        rms(i) = rms(i) + fcancmx(i,j) * rmsveg(i,j)
        rmr(i) = rmr(i) + fcancmx(i,j) * rmrveg(i,j)
        rm(i) = rm(i) + fcancmx(i,j) * rmveg(i,j)
        rg(i) = rg(i) + fcancmx(i,j) * rgveg(i,j)
        npp(i) = npp(i) + fcancmx(i,j) * nppveg(i,j)
        gpp(i) = gpp(i) + fcancmx(i,j) * gppveg(i,j)
        autores(i) = rg(i) + rm(i)
        autoresveg(i,j) = rmveg(i,j) + rgveg(i,j)

      end do ! loop 280
    end do ! loop 270
    !
    !>    Add moss GPP and rml to the grid/tile average C fluxes
    !>    for grid cells which have peatlands
    !
    do i = il1,il2
      if (ipeatland(i) > 0) then
        rgmoss(i) = anmoss(i) * grescoefmoss
        rml(i) = rml(i) + rmlmoss(i)
        rm(i) = rm(i)  + rmlmoss(i)
        rg(i) = rg(i)  + rgmoss (i)
        armoss(i) = rmlmoss(i) + rgmoss(i)
        nppmoss(i) = anmoss(i) - rgmoss(i)
        npp(i) = npp(i) + nppmoss (i)
        gpp(i) = gpp(i) + gppmoss(i)
        autores(i) = autores(i) + armoss(i)

        nppmosstep(i) = nppmoss(i) * (1.0/963.62) * deltat    ! kgC/m2/dt

      else
        nppmosstep(i) = 0.
      end if
    end do ! loop 335

  end subroutine calcNPP
  !! @}
  ! ----------------------------------------------------------------------------

  !> \ingroup ctemdriver_calcnep
  !! @{
  !> Calculate NEP (net ecosystem productivity), difference between NPP and
  !! heterotrophic respriation, for each pft and tile
  !!
  !! @author V. Arora, J. Melton
  subroutine calcNEP (il1, il2, ilg, nppveg, hetrsveg, fg, npp, hetrores, & ! In
                      nep, nepveg) ! Out

    use classicParams,   only : icc, iccp1, zero

    implicit none

    ! arguments
    integer, intent(in) :: il1             !< il1=1
    integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg
    real, intent(in) :: nppveg(:,:)    !< NPP for individual pfts,  (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: hetrsveg(:,:) !< Vegetation averaged litter and soil C respiration rates (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: fg(:)              !< Fraction of tile that is bare ground.
    real, intent(in) :: npp(:)            !< Tile-level net primary productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: hetrores(:)       !< Tile-level heterotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real, intent(out) :: nep(ilg) !< Tile-level net ecosystem productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: nepveg(ilg,iccp1) !< PFT-level net ecosystem productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    integer :: i, j

    nep(:) = 0.
    do i = il1,il2
      do j = 1,icc
        nepveg(i,j) = nppveg(i,j) - hetrsveg(i,j)
      end do
      if (fg(i) > zero) then
        nepveg(i,iccp1) = 0. - hetrsveg(i,iccp1)
      end if
      nep(i) = npp(i) - hetrores(i)
    end do

  end subroutine calcNEP
  !! @}
  ! ----------------------------------------------------------------------------
  !> \ingroup ctemdriver_calcnbp
  !! @{
  !> Calculate NBP (net biome production) for each pft by taking into account
  !! C emission losses. The disturbance routine produces emissions due to fire
  !! and while the land use change subroutine calculates emissions due to LUC.
  !! The LUC related combustion flux is assumed to be spread uniformly over the
  !! tile as it is no longer associated with any one PFT. To calculate the NBP
  !! we do not subtract LUC emissions from the PFT-level NBP but we do subtract
  !! it from the per tile NBP.
  !! @author V. Arora, J. Melton
  subroutine calcNBP (il1, il2, ilg, deltat, nepveg, fcancmx, & ! In
                      lucemcom, ltresveg, scresveg, nep, & ! In
                      glcaemls, blcaemls, stcaemls, rtcaemls, ltrcemls, & ! In/Out
                      nbpveg, dstcemls1, dstcemls3, nbp) ! Out

    use classicParams, only : icc, iccp1, iccp2

    implicit none

    ! arguments
    integer, intent(in) :: il1             !< il1=1
    integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg
    real, intent(in)    :: fcancmx(:,:)    !< max. fractional coverage of ctem's 9 pfts, but this can be
    !! modified by land-use change,and competition between pfts
    real, intent(in) :: deltat             !< CTEM (biogeochemical) time step (days)
    real, intent(in) :: nepveg(:,:)        !< Net ecosystem productivity,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(in) :: lucemcom(:)        !< Land use change (LUC) related combustion emission losses,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(in) :: ltresveg(:,:,:)    !< Litter respiration for each pft, bare fraction, and LUC product pool, \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(in) :: scresveg(:,:,:)    !< Soil carbon respiration for each pft, bare fraction, and LUC product pool, \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(in) :: nep(:)             !< Net ecosystem productivity, tile average,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$

    real, intent(inout) ::  glcaemls(:,:)  !< Green leaf carbon emission losses, \f$(kg C/m^2)\f$
    real, intent(inout) ::  blcaemls(:,:)  !< Brown leaf carbon emission losses, \f$(kg C/m^2)\f$
    real, intent(inout) ::  rtcaemls(:,:)  !< Root carbon emission losses, \f$(kg C/m^2)\f$
    real, intent(inout) ::  stcaemls(:,:)  !< Stem carbon emission losses, \f$(kg C/m^2)\f$
    real, intent(inout) ::  ltrcemls(:,:)  !< Litter carbon emission losses, \f$(kg C/m^2)\f$

    real, intent(out) :: nbpveg(ilg,iccp1) !< Net biome productivity,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(out) :: dstcemls1(ilg)    !< grid ave. carbon emission losses due to disturbance, vegetation,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(out) :: dstcemls3(ilg)    !< grid ave. carbon emission losses due to disturbance, litter,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(out) :: nbp(ilg)          !< Net biome productivity, tile average, \f$\mu mol CO_2 m^{-2} s^{-1}\f$

    ! Local vars
    integer :: i, j
    real :: dscemlv1(ilg,icc)  !< Disturbance emission losses from plants, \f$(kg C/m^2)\f$
    real :: dscemlv2(ilg,icc)  !< Disturbance emission losses from plants and litter, \f$(kg C/m^2)\f$
    real :: dstcemls2(ilg)     !< grid ave. carbon emission losses due to disturbance, total,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$

    !--

    do i = il1,il2
      do j = 1,icc
        dscemlv1(i,j) = glcaemls(i,j) + blcaemls(i,j) + stcaemls(i,j) + rtcaemls(i,j)
        dscemlv2(i,j) = dscemlv1(i,j) + ltrcemls(i,j)

        ! Convert \f$(kg C/m^2)\f$ emitted in one day into u mol co2/m2.sec before
        ! subtracting emission losses from nep.
        nbpveg(i,j) = nepveg(i,j) - dscemlv2(i,j) * (963.62/deltat)

      end do ! loop 101

      ! For accounting purposes, we also need to account for the bare fraction
      ! NBP. Since there is no fire on the bare, we use 0.
      nbpveg(i,iccp1) = nepveg(i,iccp1)   - 0.

    end do ! loop 100
    !>
    !! Calculate grid. averaged rate of carbon emissions due to fire in u-mol co2/m2.sec.
    !! Convert all emission losses from \f$(kg C/m^2)\f$ emitted in 1 day to u-mol co2/m2.sec.
    !! Calculate grid averaged carbon emission losses from litter.
    !!
    dstcemls1 = 0.0
    dstcemls2 = 0.0
    do j = 1,icc
      do i = il1,il2
        dstcemls1(i) = dstcemls1(i) + fcancmx(i,j) * dscemlv1(i,j) * (963.62 / deltat)
        dstcemls2(i) = dstcemls2(i) + fcancmx(i,j) * dscemlv2(i,j) * (963.62 / deltat)
        glcaemls(i,j) = glcaemls(i,j) * (963.62/deltat)
        blcaemls(i,j) = blcaemls(i,j) * (963.62/deltat)
        stcaemls(i,j) = stcaemls(i,j) * (963.62/deltat)
        rtcaemls(i,j) = rtcaemls(i,j) * (963.62/deltat)
        ltrcemls(i,j) = ltrcemls(i,j) * (963.62/deltat)
      end do ! loop 104
    end do ! loop 103

    ! For the tile-level NBP, we include the disturbance emissions as well as
    ! respiration from the paper (litter) and furniture (soil carbon) pools (LUC
    ! product pools). Also include here the instantaneous emissions due to LUC.
    nbp(:) = 0.
    dstcemls3 = 0.0
    do i = il1,il2
      nbp(i) = nep(i) - dstcemls2(i) - (ltresveg(i,iccp2,1) + scresveg(i,iccp2,1)) - lucemcom(i)
      dstcemls3(i) = dstcemls2(i) - dstcemls1(i)  ! litter is total - vegetation.
    end do ! loop 105

  end subroutine calcNBP
  !! @}

  !> \namespace ctemdriver
  !! Central module that contains the ctem driver and associated subroutines.
  !!
  !! The basic model structure of CTEM includes three live vegetation components
  !! (leaf (L), stem (S) and root (R)) and two dead carbon pools (litter or
  !! detritus (D) and soil carbon (H)). The amount of carbon in these pools
  !! (\f$C_\mathrm{L}\f$, \f$C_\mathrm{S}\f$, \f$C_\mathrm{R}\f$, \f$C_\mathrm{D}\f$,
  !! \f$C_\mathrm{H}\f$, \f$kgC m^{-2}\f$) is tracked prognostically through the
  !! fluxes in and out of them. The rate change equations for carbon in these
  !! pools are summarized in Sect. \ref{rate_change_eqns} after the processes
  !! leading to the calculation of fluxes in and out of these pools are introduced
  !! in the following sections.
  !!
end module ctemDriver
