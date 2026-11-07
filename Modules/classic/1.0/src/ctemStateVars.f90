!> \file
!> Contains the biogeochemistry-related variable type structures.
!! @author J. Melton
!! Variable types herein:
!!
!! 1. c_switch (ctem_switches) - Switches for running the model, read from the joboptions file
!! 2. vrot (veg_rot) - CTEM's 'rot' vars
!! 3. vgat (veg_gat) - CTEM's 'gat' vars
!! 4. ctem_tile (ctem_tile_level) - CTEM's variables per tile
!! 5. ctem_mo (ctem_monthly) - CTEM's variables monthly averaged (per pft)
!! 6. ctem_grd_mo (ctem_gridavg_monthly)  - CTEM's grid average monthly values
!! 7. ctem_tile_mo (ctem_tileavg_monthly) - CTEM's variables per tile monthly values
!! 8. ctem_yr (ctem_annual) - CTEM's average annual values (per PFT)
!! 9. ctem_grd_yr (ctem_gridavg_annual) - CTEM's grid average annual values
!! 10. ctem_tile_yr (ctem_tileavg_annual) - CTEM's variables per tile annual values

module ctemStateVars

  ! J. Melton Apr 2015

  use classicParams,  only : nlat, nmos, ilg, ican, ignd, icp1, icc, iccp2, iccp1, monthend, &
                             mmday, modelpft, l2max, deltat, abszero, monthdays, seed, crop, NBS
  use model_variables, only : vs

  implicit none

  public :: allocCtemVars ! Allocate the biogeochemistry (CTEM) variables in preparation for a simulation
  public :: initRowVars     ! Initializes 'row' variables
  public :: resetMonthEnd   ! Resets monthly variables at month end in preparation for next month
  public :: resetYearEnd    ! Resets annual variables in preparation for next year
  public :: resetMosaicAccum ! Resets physics accumulator variables (used as input to CTEM) after CTEM has been called


  !=================================================================================
  !> Switches for running the model, read from the joboptions file
  type ctem_switches

    logical :: projectedGrid    !< True if you have a projected lon lat grid, false if not. Projected grids can only have
    !! regions referenced by the indexes, not coordinates, when running a sub-region
    logical :: ctem_on          !< True if this run includes the biogeochemistry parameterizations (CTEM)
    integer :: metLoop          !< no. of times the meteorological data is to be looped over. this
    !< option is useful to equilibrate CTEM's C pools
    logical :: leap             !< set to true if the meteorological file includes leap years

    integer :: spinfast         !< set this to a number >1 (up to ~10) to spin up soil carbon pool faster, set to 1
    !< for final production simulations
    integer :: useTracer        !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    character(350) :: tracerCO2file !< Tracer CO2 file, this file needs to be correctly chosen for the useTracer option. It uses the transientCO2
    !! and fixedYearCO2 switches to determine how it is read in.
    integer :: readMetStartYear !< First year of meteorological forcing to read in from the met file
    integer :: readMetEndYear   !< Last year of meteorological forcing to read in from the met file

    logical :: transientCO2     !< use \f$CO_2\f$ time series, if false, fixedYearCO2 is used
    character(350) :: CO2File   !< Location of the netcdf file containing atmospheric \f$CO_2\f$ values
    integer :: fixedYearCO2     !< set the year to use for atmospheric \f$CO_2\f$ if transientCO2 is false. (ppmv)

    logical :: doMethane        !< Setting to true will enable wetland CH4 and soil uptake of CH4 subroutines and also require the 
                                !! [CH4] input file.
    logical :: transientCH4     !< use \f$CH_4\f$ time series, if false, fixedYearCH4 is used
    character(350) :: CH4File   !< Location of the netcdf file containing atmospheric \f$CH_4\f$ values
    integer :: fixedYearCH4     !< set the year to use for atmospheric \f$CH_4\f$ if transientCH4 is false. (ppmv)

    logical :: dofire           !< boolean, if true allow fire disturbance, if false no fire occurs.
    logical :: transientPOPD    !< if set true use time series of population density data to calculate
    !< fire extinguishing probability and probability of fire due to human causes, false use
    !< value from single year (fixedYearPOPD)
    character(350) :: POPDFile  !< Location of the netcdf file containing population density values
    integer :: fixedYearPOPD    !< set the year to use for population density values if transientPOPD is false. (\f$people/km^2\f$)

    logical :: transientLGHT    !< use lightning strike time series, otherwise use fixedYearLGHT
    character(350) :: LGHTFile  !< Location of the netcdf file containing lightning strike values
    integer :: fixedYearLGHT    !< set the year to use for lightning strikes if transientLGHT is false.

    logical :: PFTCompetition   !< logical boolean telling if competition between pfts is on or not
    logical :: start_bare       !< set this to true if competition is on, and if you wish to start from bare ground.
    !< if this is set to false, the init file info will be used to set up the run.
    !< NOTE: This still keeps the crop fractions (while setting all pools to zero)
    logical :: inibioclim       !< switch telling if bioclimatic parameters are being initialized
    !< from scratch (false) or being initialized from some spun up
    !< values(true).

    logical :: lnduseon         !< If true then the land cover is read in from LUCFile and changes annually
    character(350) :: LUCFile   !< Location of the netcdf file containing land use change information
    integer :: fixedYearLUC     !< Set the year to use for land cover if lnduseon is false. If set to -9999,
    !< we use the PFT distribution found in the initialization file. Any other year
    !< we search for that year in the LUCFile

    logical :: transientOBSWETF    !< use observed wetland fraction time series, otherwise use fixedYearOBSWETF
    character(350) :: OBSWETFFile  !< Location of the netcdf file containing observed wetland fraction
    integer :: fixedYearOBSWETF    !< set the year to use for observed wetland fraction if transientOBSWETF is false.

    character(350) :: metFileFss        !< location of the incoming shortwave radiation meteorology file
    character(350) :: metFileFdl        !< location of the incoming longwave radiation meteorology file
    character(350) :: metFilePre        !< location of the precipitation meteorology file
    character(350) :: metFileTa         !< location of the air temperature meteorology file
    character(350) :: metFileQa         !< location of the specific humidity meteorology file
    character(350) :: metFileUv         !< location of the wind speed meteorology file
    character(350) :: metFilePres       !< location of the atmospheric pressure meteorology file
    character(350) :: init_file         !< location of the netcdf initialization file
    character(350) :: rs_file_to_overwrite !< location of the netcdf file that will be written for the restart file
    character(350) :: runparams_file    !< location of the namelist file containing the model parameters
    character(350) :: Comment           !< Comment about the run that will be written to the output netcdfs

    character(350) :: output_directory  !< Directory where the output netcdfs will be placed
    character(350) :: xmlFile           !< location of the xml file that outlines the possible netcdf output files

    logical :: doperpftoutput           !< Switch for making extra output files that are at the per PFT level
    logical :: dopertileoutput          !< Switch for making extra output files that are at the per tile level
    logical :: doChecksums              !< Switch for doing checksum calculations for the run

    logical :: doAnnualOutput           !< Switch for making annual output files
    logical :: doMonthOutput            !< Switch for making monthly output files
    integer :: jmosty                   !< Year to start writing out the monthly output files. If you want to write monthly outputs right
    !< from the start then put in a negative number (like -9999)

    logical :: doDayOutput              !< Switch for making daily output files
    integer :: jdstd                    !< day of the year to start writing the daily output
    integer :: jdendd                   !< day of the year to stop writing the daily output
    integer :: jdsty                    !< simulation year (iyear) to start writing the daily output
    integer :: jdendy                   !< simulation year (iyear) to stop writing the daily output

    logical :: doHhOutput               !< Switch for making half hourly output files
    integer :: jhhstd                   !< day of the year to start writing the half-hourly output
    integer :: jhhendd                  !< day of the year to stop writing the half-hourly output
    integer :: jhhsty                   !< simulation year (iyear) to start writing the half-hourly output
    integer :: jhhendy                  !< simulation year (iyear) to stop writing the half-hourly output

    ! Physics switches:

    integer :: idisp    !< if idisp=0, vegetation displacement heights are ignored,
    !< because the atmospheric model considers these to be part
    !< of the "terrain".
    !< if idisp=1, vegetation displacement heights are calculated.

    integer :: izref    !< if izref=1, the bottom of the atmospheric model is taken
    !< to lie at the ground surface.
    !< if izref=2, the bottom of the atmospheric model is taken
    !< to lie at the local roughness height.

    integer :: islfd    !< if islfd=0, drcoef is called for surface stability corrections
    !< and the original gcm set of screen-level diagnostic calculations
    !< is done.
    !< if islfd=1, drcoef is called for surface stability corrections
    !< and sldiag is called for screen-level diagnostic calculations.
    !< if islfd=2, flxsurfz is called for surface stability corrections
    !< and diasurf is called for screen-level diagnostic calculations.

    integer :: ipcp     !< if ipcp=1, the rainfall-snowfall cutoff is taken to lie at 0 c.
    !< if ipcp=2, a linear partitioning of precipitation betweeen
    !< rainfall and snowfall is done between 0 c and 2 c.
    !< if ipcp=3, rainfall and snowfall are partitioned according to
    !< a polynomial curve between 0 c and 6 c.

    integer :: ITC !< itc, itcg and itg are switches to choose the iteration scheme to
    !< be used in calculating the canopy or ground surface temperature
    !< respectively.  if the switch is set to 1, a bisection method is
    !< used; if to 2, the newton-raphson method is used.
    integer :: ITCG !< itc, itcg and itg are switches to choose the iteration scheme to
    !< be used in calculating the canopy or ground surface temperature
    !< respectively.  if the switch is set to 1, a bisection method is
    !< used; if to 2, the newton-raphson method is used.
    integer :: ITG !< itc, itcg and itg are switches to choose the iteration scheme to
    !< be used in calculating the canopy or ground surface temperature
    !< respectively.  if the switch is set to 1, a bisection method is
    !< used; if to 2, the newton-raphson method is used.

    integer :: IPAI !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !< plant area index, vegetation height, canopy albedo, snow albedo
    !< and soil albedo respectively calculated by class are used.
    !< if any of these switches is set to 1, the value of the
    !< corresponding parameter calculated by class is overridden by
    !< a user-supplied input value.
    integer :: IHGT !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !< plant area index, vegetation height, canopy albedo, snow albedo
    !< and soil albedo respectively calculated by class are used.
    !< if any of these switches is set to 1, the value of the
    !< corresponding parameter calculated by class is overridden by
    !< a user-supplied input value.
    integer :: IALC !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !< plant area index, vegetation height, canopy albedo, snow albedo
    !< and soil albedo respectively calculated by class are used.
    !< if any of these switches is set to 1, the value of the
    !< corresponding parameter calculated by class is overridden by
    !< a user-supplied input value.
    integer :: IALS !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !< plant area index, vegetation height, canopy albedo, snow albedo
    !< and soil albedo respectively calculated by class are used.
    !< if any of these switches is set to 1, the value of the
    !< corresponding parameter calculated by class is overridden by
    !< a user-supplied input value.
    integer :: IALG !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !< plant area index, vegetation height, canopy albedo, snow albedo
    !< and soil albedo respectively calculated by class are used.
    !< if any of these switches is set to 1, the value of the
    !< corresponding parameter calculated by class is overridden by
    !< a user-supplied input value.
    integer :: isnoalb !< if isnoalb is set to 0, the original two-band snow albedo algorithms are used.
    !< if it is set to 1, the new four-band routines are used.

  end type ctem_switches

  type (ctem_switches), save, target :: c_switch

  !=================================================================================
  !> CTEM's 'rot' vars
  type veg_rot

    logical, allocatable, dimension(:,:,:) :: pftexist  !< logical array indicating pfts exist (t) or not (f)
    integer, allocatable, dimension(:,:,:) :: lfstatus  !< leaf phenology status
    integer, allocatable, dimension(:,:,:) :: pandays   !< days with positive net photosynthesis (an) for use in
    !< the phenology subroutine
    real, allocatable, dimension(:,:,:) :: gleafmas     !< green leaf mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: bleafmas     !< brown leaf mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: stemmass     !< stem mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: rootmass     !< root mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: pstemmass    !< stem mass from previous timestep, is value before fire. used by burntobare subroutine
    real, allocatable, dimension(:,:,:) :: pgleafmass   !< root mass from previous timestep, is value before fire. used by burntobare subroutine
    real, allocatable, dimension(:,:,:) :: fcancmx      !< max. fractional coverage of ctem's pfts, but this can be
    !< modified by land-use change,and competition between pfts
    real, allocatable, dimension(:,:,:) :: ailcg        !< Green LAI for CTEM's pfts
    real, allocatable, dimension(:,:,:) :: ailcgs       !< Green LAI for canopy over snow sub-area
    real, allocatable, dimension(:,:,:) :: fcancs       !< Fraction of canopy over snow for ctem's pfts
    real, allocatable, dimension(:,:,:) :: fcanc        !< Fractional coverage of carbon pfts, canopy over snow
    real, allocatable, dimension(:,:,:) :: co2i1cg      !< Intercellular CO2 conc for pfts for canopy over ground subarea (Pa) - for single/sunlit leaf
    real, allocatable, dimension(:,:,:) :: co2i1cs      !< Same as above but for shaded leaf (above being co2i1cg)
    real, allocatable, dimension(:,:,:) :: co2i2cg      !< Intercellular CO2 conc for pfts for canopy over snowsubarea (pa) - for single/sunlit leaf
    real, allocatable, dimension(:,:,:) :: co2i2cs      !< Same as above but for shaded leaf (above being co2i2cg)
    real, allocatable, dimension(:,:,:) :: ancsveg      !< Net photosynthetic rate for CTEM's pfts for canopy over snow subarea
    real, allocatable, dimension(:,:,:) :: ancgveg      !< Net photosynthetic rate for CTEM's pfts for canopy over ground subarea
    real, allocatable, dimension(:,:,:) :: rmlcsveg     !< Leaf respiration rate for CTEM' pfts forcanopy over snow subarea
    real, allocatable, dimension(:,:,:) :: rmlcgveg     !< Leaf respiration rate for CTEM' pfts forcanopy over ground subarea
    real, allocatable, dimension(:,:,:) :: slai         !< storage/imaginary lai for phenology purposes
    real, allocatable, dimension(:,:,:) :: ailcb        !< brown lai for ctem's 9 pfts. for now we assume only grasses can have brown lai
    real, allocatable, dimension(:,:,:) :: flhrloss     !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, allocatable, dimension(:,:,:) :: grwtheff     !< growth efficiency. change in biomass per year per unit max.
    !< lai (\f$kg c/m^2\f$)/(m2/m2),for use in mortality subroutine
    real, allocatable, dimension(:,:,:) :: lystmmas     !< stem mass at the end of last year
    real, allocatable, dimension(:,:,:) :: lyrotmas     !< root mass at the end of last year
    real, allocatable, dimension(:,:,:) :: tymaxlai     !< this year's maximum lai
    real, allocatable, dimension(:,:,:) :: stmhrlos     !< stem harvest loss for crops, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: vgbiomas_veg !< vegetation biomass for each pft
    real, allocatable, dimension(:,:,:) :: emit_co2     !< carbon dioxide
    real, allocatable, dimension(:,:,:) :: emit_co      !< carbon monoxide
    real, allocatable, dimension(:,:,:) :: emit_ch4     !< methane
    real, allocatable, dimension(:,:,:) :: emit_nmhc    !< non-methane hydrocarbons
    real, allocatable, dimension(:,:,:) :: emit_h2      !< hydrogen gas
    real, allocatable, dimension(:,:,:) :: emit_nox     !< nitrogen oxides
    real, allocatable, dimension(:,:,:) :: emit_n2o     !< nitrous oxide
    real, allocatable, dimension(:,:,:) :: emit_pm25    !< particulate matter less than 2.5 um in diameter
    real, allocatable, dimension(:,:,:) :: emit_tpm     !< total particulate matter
    real, allocatable, dimension(:,:,:) :: emit_tc      !< total carbon
    real, allocatable, dimension(:,:,:) :: emit_oc      !< organic carbon
    real, allocatable, dimension(:,:,:) :: emit_bc      !< black carbon
    real, allocatable, dimension(:,:,:) :: burnvegf     !< per PFT fraction burned of that PFT's area
    real, allocatable, dimension(:,:,:) :: smfuncveg    !<
    real, allocatable, dimension(:,:,:) :: bterm        !< biomass term for fire probabilty calc
    real, allocatable, dimension(:,:,:) :: mterm        !< moisture term for fire probabilty calc
    real, allocatable, dimension(:,:,:) :: bmasveg      !< total (gleaf + stem + root) biomass for each ctem pft, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: veghght      !< vegetation height (meters)
    real, allocatable, dimension(:,:,:) :: rootdpth     !< 99% soil rooting depth (meters)
    !< both veghght & rootdpth can be used as diagnostics to see
    !< how vegetation grows above and below ground, respectively
    real, allocatable, dimension(:,:,:) :: tltrleaf     !< total leaf litter fall rate (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:,:) :: tltrstem     !< total stem litter fall rate (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:,:) :: tltrroot     !< total root litter fall rate (u-mol co2/m2.sec)
    real, allocatable, dimension(:,:,:) :: leaflitr     !< leaf litter fall rate (u-mol co2/m2.sec). this leaf litter
    !< does not include litter generated due to mortality/fire
    real, allocatable, dimension(:,:,:) :: roottemp     !< root temperature, k
    real, allocatable, dimension(:,:,:) :: afrleaf      !< allocation fraction for leaves
    real, allocatable, dimension(:,:,:) :: afrstem      !< allocation fraction for stem
    real, allocatable, dimension(:,:,:) :: afrroot      !< allocation fraction for root
    real, allocatable, dimension(:,:,:) :: wtstatus     !< soil water status used for calculating allocation fractions
    real, allocatable, dimension(:,:,:) :: ltstatus     !< light status used for calculating allocation fractions
    real, allocatable, dimension(:,:,:) :: gppveg       !< ! gross primary productity for each pft
    real, allocatable, dimension(:,:,:) :: nppveg       !< npp for individual pfts, u-mol co2/m2.sec
    real, allocatable, dimension(:,:,:) :: autoresveg   !<
    real, allocatable, dimension(:,:,:) :: rmlvegacc    !<
    real, allocatable, dimension(:,:,:) :: rmsveg       !< stem maintenance resp. rate for each pft
    real, allocatable, dimension(:,:,:) :: rmrveg       !< root maintenance resp. rate for each pft
    real, allocatable, dimension(:,:,:) :: rgveg        !< growth resp. rate for each pft
    real, allocatable, dimension(:,:,:) :: litrfallveg  !< litter fall in \f$kg c/m^2\f$ for each pft
    real, allocatable, dimension(:,:,:) :: rothrlos     !< root death as crops are harvested, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: pfcancmx     !< previous year's fractional coverages of pfts
    real, allocatable, dimension(:,:,:) :: nfcancmx     !< next year's fractional coverages of pfts
    real, allocatable, dimension(:,:,:) :: anveg        !< net photosynthesis rate for each pft
    real, allocatable, dimension(:,:,:) :: rmlveg       !< leaf maintenance resp. rate for each pft


    ! allocated with nlat,nmos:
    real, allocatable, dimension(:,:) :: gavglai               !< grid averaged green leaf area index
    real, allocatable, dimension(:,:) :: co2conc               !< ATMOS. CO2 CONC. IN PPM
    real, allocatable, dimension(:,:) :: ch4conc               !<
    real, allocatable, dimension(:,:) :: canres                !<
    real, allocatable, dimension(:,:) :: vgbiomas              !< grid averaged vegetation biomass, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: gavgltms              !< grid averaged litter mass, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: gavgscms              !< grid averaged soil c mass, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: burnfrac              !< areal :: fraction burned due to fire for every grid cell (%)
    real, allocatable, dimension(:,:) :: popdin                !< population density \f$(people / km^2)\f$
    real, allocatable, dimension(:,:) :: lterm                 !< lightning term for fire probabilty calc
    real, allocatable, dimension(:,:) :: extnprob              !< fire extingusinging probability
    real, allocatable, dimension(:,:) :: prbfrhuc              !< probability of fire due to human causes
    real, allocatable, dimension(:,:) :: rml                   !< leaf maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: rms                   !< stem maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: rmr                   !< root maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: ch4WetSpec               !< methane flux from wetlands calculated using hetrores in umol ch4/m2.s
    real, allocatable, dimension(:,:) :: wetfdyn               !< dynamic wetland fraction
    real, allocatable, dimension(:,:) :: wetfrac_pres          !< Prescribed wetland fraction read in from OBSWETFFile
    real, allocatable, dimension(:,:) :: ch4WetDyn                !< methane flux from wetlands calculated using hetrores and wetfdyn, in umol ch4/m2.s
    real, allocatable, dimension(:,:) :: ch4_soills            !< Methane uptake into the soil column (\f$mg CH_4 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: lucemcom              !< land use change (luc) related combustion emission losses (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: lucltrin              !< luc related inputs to litter pool (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: lucsocin              !< luc related inputs to soil c pool (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: npp                   !< net primary productivity
    real, allocatable, dimension(:,:) :: nep                   !< net ecosystem productivity
    real, allocatable, dimension(:,:) :: nbp                   !< net biome productivity
    real, allocatable, dimension(:,:) :: gpp                   !< gross primary productivity
    real, allocatable, dimension(:,:) :: hetrores              !< heterotrophic respiration
    real, allocatable, dimension(:,:) :: autores               !< autotrophic respiration
    real, allocatable, dimension(:,:) :: soilcresp             !<
    real, allocatable, dimension(:,:) :: rm                    !< maintenance respiration
    real, allocatable, dimension(:,:) :: rg                    !< growth respiration
    real, allocatable, dimension(:,:) :: litres                !< litter respiration
    real, allocatable, dimension(:,:) :: socres                !< soil carbon respiration
    real, allocatable, dimension(:,:) :: dstcemls              !< carbon emission losses due to disturbance, mainly fire
    real, allocatable, dimension(:,:) :: litrfall              !< total litter fall (from leaves, stem, and root) due to
    !< all causes (mortality,turnover,and disturbance)
    real, allocatable, dimension(:,:) :: humiftrs              !< transfer of humidified litter from litter to soil c pool
    real, allocatable, dimension(:,:) :: cfluxcg               !<
    real, allocatable, dimension(:,:) :: cfluxcs               !<
    real, allocatable, dimension(:,:) :: dstcemls3             !< carbon emission losses due to disturbance (fire at present) from litter pool
    real, allocatable, dimension(:,:) :: tcanrs                !<
    real, allocatable, dimension(:,:) :: tsnors                !<
    real, allocatable, dimension(:,:) :: tpndrs                !<
    real, allocatable, dimension(:,:) :: uvaccrow_m            !<
    real, allocatable, dimension(:,:) :: vvaccrow_m            !<
    real, allocatable, dimension(:,:) :: qevpacc_m_save        !<
    real, allocatable, dimension(:,:) :: twarmm                !< temperature of the warmest month (c)
    real, allocatable, dimension(:,:) :: tcoldm                !< temperature of the coldest month (c)
    real, allocatable, dimension(:,:) :: gdd5                  !< growing degree days above 5 c
    real, allocatable, dimension(:,:) :: aridity               !< aridity index, ratio of potential evaporation to precipitation
    real, allocatable, dimension(:,:) :: srplsmon              !< number of months in a year with surplus water i.e. precipitation more than potential evaporation
    real, allocatable, dimension(:,:) :: defctmon              !< number of months in a year with water deficit i.e. precipitation less than potential evaporation
    real, allocatable, dimension(:,:) :: anndefct              !< annual water deficit (mm)
    real, allocatable, dimension(:,:) :: annsrpls              !< annual water surplus (mm)
    real, allocatable, dimension(:,:) :: annpcp                !< annual precipitation (mm)
    real, allocatable, dimension(:,:) :: dry_season_length     !< length of dry season (months)

    integer, allocatable, dimension(:,:) :: ipeatland !< Peatland flag: 0 = not a peatland, 1 = bog, 2 = fen
    real, allocatable, dimension(:,:) :: litrmsmoss
    real, allocatable, dimension(:,:) :: Cmossmas
    real, allocatable, dimension(:,:) :: dmoss
    real, allocatable, dimension(:,:) :: nppmoss
    real, allocatable, dimension(:,:) :: rmlmoss
    real, allocatable, dimension(:,:) :: gppmoss
    real, allocatable, dimension(:,:) :: anmoss
    real, allocatable, dimension(:,:) :: armoss
    real, allocatable, dimension(:,:) :: peatdep
    real, allocatable, dimension(:,:) :: pdd

    ! allocated with nlat,nmos,ican:
    real, allocatable, dimension(:,:,:) :: zolnc            !< lumped log of roughness length for class' 4 pfts
    real, allocatable, dimension(:,:,:) :: ailc             !< lumped lai for class' 4 pfts
    real, allocatable, dimension(:,:,:) :: cmasvegc         !< total canopy mass for each of the 4 class pfts. recall that
    !< class requires canopy mass as an input,and this is now provided by ctem. \f$kg/m^2\f$.
    real, allocatable, dimension(:,:,:) :: alvsctm          !<
    real, allocatable, dimension(:,:,:) :: paic             !< plant area index for class' 4 pfts. this is the sum of leaf
    !< area index and stem area index.
    real, allocatable, dimension(:,:,:) :: slaic            !< storage lai. this will be used as min. lai that class sees
    !< so that it doesn't blow up in its stomatal conductance calculations.
    real, allocatable, dimension(:,:,:) :: alirctm          !<
    real, allocatable, dimension(:,:,:) :: csum             !<

    ! allocated with nlat,nmos,ican,ignd:
    real, allocatable, dimension(:,:,:,:) :: rmatc       !< fraction of roots for each of class' 4 pfts in each soil layer

    ! allocated with nlat,nmos,icc,ignd:
    real, allocatable, dimension(:,:,:,:) :: rmatctem     !< fraction of roots for each of ctem's 9 pfts in each soil layer

    ! allocated with nlat,nmos,iccp1:
    real, allocatable, dimension(:,:,:) :: nepveg      !< net ecosystem productity for each pft
    real, allocatable, dimension(:,:,:) :: nbpveg      !< net biome productity for bare fraction OR net biome productity for each pft
    real, allocatable, dimension(:,:,:) :: hetroresveg !<

    ! allocated with nlat,nmos,iccp2,ignd:
    ! COMBAK PERLAY
    real, allocatable, dimension(:,:,:) :: litrmass    !< litter mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: soilcmas    !< soil carbon mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: litresveg   !<
    real, allocatable, dimension(:,:,:) :: soilcresveg !<
    real, allocatable, dimension(:,:,:) :: humiftrsveg !<
    ! real, allocatable, dimension(:,:,:,:) :: litrmass    !< litter mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    ! real, allocatable, dimension(:,:,:,:) :: soilcmas    !< soil carbon mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    ! real, allocatable, dimension(:,:,:,:) :: litresveg   !<
    ! real, allocatable, dimension(:,:,:,:) :: soilcresveg !<
    ! real, allocatable, dimension(:,:,:,:) :: humiftrsveg !<
    ! COMBAK PERLAY

    ! allocated with nlat,nmos,{some number}:
    integer, allocatable, dimension(:,:,:) :: colddays          !< cold days counter for tracking days below a certain
    !< temperature threshold for ndl dcd and crop pfts.
    real, allocatable, dimension(:,:,:)  :: slopefrac          !< prescribed fraction of wetlands based on slope
    !< only(0.025,0.05,0.1,0.15,0.20,0.25,0.3 and 0.35 percent slope thresholds)

    ! allocated with nlat:
    real, allocatable, dimension(:)    :: dayl_max        !< maximum daylength for that location (hours)
    real, allocatable, dimension(:)    :: dayl            !< daylength for that location (hours)

  end type veg_rot

  type (veg_rot),save,target :: vrot

  !=================================================================================
  !> CTEM's 'gat' vars
  type veg_gat

    ! This is the basic data structure that contains the state variables
    ! for the Plant Functional Type (PFT). The dimensions are ilg,{icc,iccp1,iccp2}

    real, allocatable, dimension(:,:) :: gleafmas   !< Green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: bleafmas   !< Brown leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: stemmass   !< Stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: rootmass   !< Root mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: pstemmass  !< Stem mass from previous timestep, is value before fire. used by burntobare subroutine
    real, allocatable, dimension(:,:) :: pgleafmass !< Root mass from previous timestep, is value before fire. used by burntobare subroutine
    real, allocatable, dimension(:,:) :: fcancmx    !< max. fractional coverage of ctem's 9 pfts, but this can be
    !< modified by land-use change,and competition between pfts
    real, allocatable, dimension(:) :: gavglai        !< grid averaged green leaf area index

    real, allocatable, dimension(:) :: lightng        !< total lightning frequency, flashes/km2.year

    real, allocatable, dimension(:,:) :: zolnc     !< lumped log of roughness length for class' 4 pfts
    real, allocatable, dimension(:,:) :: ailc      !< lumped lai for class' 4 pfts

    real, allocatable, dimension(:,:) :: ailcg      !< green lai for ctem's 9 pfts
    real, allocatable, dimension(:,:) :: ailcgs     !< GREEN LAI FOR CANOPY OVER SNOW SUB-AREA
    real, allocatable, dimension(:,:) :: fcancs     !< FRACTION OF CANOPY OVER SNOW FOR CTEM's 9 PFTs
    real, allocatable, dimension(:,:) :: fcanc      !< FRACTIONAL COVERAGE OF 8 CARBON PFTs, CANOPY OVER GROUND

    real, allocatable, dimension(:)   :: co2conc    !< ATMOS. CO2 CONC. IN PPM
    real, allocatable, dimension(:)   :: ch4conc    !<

    real, allocatable, dimension(:,:) :: co2i1cg    !< INTERCELLULAR CO2 CONC FOR 8 PFTs FOR CANOPY OVER GROUND SUBAREA (Pa) - FOR SINGLE/SUNLIT LEAF
    real, allocatable, dimension(:,:) :: co2i1cs    !< SAME AS ABOVE BUT FOR SHADED LEAF (above being co2i1cg)
    real, allocatable, dimension(:,:) :: co2i2cg    !< INTERCELLULAR CO2 CONC FOR 8 PFTs FOR CANOPY OVER SNOWSUBAREA (Pa) - FOR SINGLE/SUNLIT LEAF
    real, allocatable, dimension(:,:) :: co2i2cs    !< SAME AS ABOVE BUT FOR SHADED LEAF (above being co2i2cg)
    real, allocatable, dimension(:,:) :: ancsveg    !< net photosynthetic rate for ctems 9 pfts for canopy over snow subarea
    real, allocatable, dimension(:,:) :: ancgveg    !< net photosynthetic rate for ctems 9 pfts for canopy over ground subarea
    real, allocatable, dimension(:,:) :: rmlcsveg   !< leaf respiration rate for ctems 9 pfts forcanopy over snow subarea
    real, allocatable, dimension(:,:) :: rmlcgveg   !< leaf respiration rate for ctems 9 pfts forcanopy over ground subarea
    real, allocatable, dimension(:,:) :: slai       !< storage/imaginary lai for phenology purposes
    real, allocatable, dimension(:,:) :: ailcb      !< brown lai for ctem's 9 pfts. for now we assume only grasses can have brown lai
    real, allocatable, dimension(:)   :: canres     !<
    real, allocatable, dimension(:,:) :: flhrloss   !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1

    real, allocatable, dimension(:,:) :: grwtheff   !< growth efficiency. change in biomass per year per unit max.
    !< lai (\f$kg c/m^2\f$)/(m2/m2),for use in mortality subroutine
    real, allocatable, dimension(:,:) :: lystmmas   !< stem mass at the end of last year
    real, allocatable, dimension(:,:) :: lyrotmas   !< root mass at the end of last year
    real, allocatable, dimension(:,:) :: tymaxlai   !< this year's maximum lai
    real, allocatable, dimension(:)   :: vgbiomas   !< grid averaged vegetation biomass, \f$kg c/m^2\f$
    real, allocatable, dimension(:)   :: gavgltms   !< grid averaged litter mass, \f$kg c/m^2\f$
    real, allocatable, dimension(:)   :: gavgscms   !< grid averaged soil c mass, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: stmhrlos   !< stem harvest loss for crops, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: rmatc    !< fraction of roots for each of class' 4 pfts in each soil layer
    real, allocatable, dimension(:,:,:) :: rmatctem !< fraction of roots for each of ctem's 9 pfts in each soil layer
    ! COMBAK PERLAY
    real, allocatable, dimension(:,:) :: litrmass   !< Litter mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: soilcmas   !< Soil carbon mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    ! real, allocatable, dimension(:,:,:) :: litrmass   !< Litter mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    ! real, allocatable, dimension(:,:,:) :: soilcmas   !< Soil carbon mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    ! COMBAK PERLAY
    real, allocatable, dimension(:,:) :: vgbiomas_veg !< vegetation biomass for each pft

    real, allocatable, dimension(:,:) :: emit_co2   !< carbon dioxide (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:,:) :: emit_co    !< carbon monoxide (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:,:) :: emit_ch4   !< methane (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:,:) :: emit_nmhc  !< non-methane hydrocarbons (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:,:) :: emit_h2    !< hydrogen gas (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:,:) :: emit_nox   !< nitrogen oxides (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:,:) :: emit_n2o   !< nitrous oxide (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:,:) :: emit_pm25  !< particulate matter less than 2.5 um in diameter (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:,:) :: emit_tpm   !< total particulate matter (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:,:) :: emit_tc    !< total carbon (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:,:) :: emit_oc    !< organic carbon (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:,:) :: emit_bc    !< black carbon (kg <species> $m^{-2}$$s^{-1}$)
    real, allocatable, dimension(:)   :: burnfrac   !< areal :: fraction burned due to fire for every grid cell (%)
    real, allocatable, dimension(:,:) :: burnvegf   !< per PFT fraction burned of that PFT's area
    real, allocatable, dimension(:,:) :: smfuncveg  !<
    real, allocatable, dimension(:)   :: popdin     !< population density (people / \f$km^2\f$)
    real, allocatable, dimension(:,:) :: bterm      !< biomass term for fire probabilty calc
    real, allocatable, dimension(:)   :: lterm      !< lightning term for fire probabilty calc
    real, allocatable, dimension(:,:) :: mterm      !< moisture term for fire probabilty calc
    real, allocatable, dimension(:,:) :: glcaemls  !< green leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: blcaemls  !< brown leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: rtcaemls  !< root carbon emission disturbance losses, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: stcaemls  !< stem carbon emission disturbance losses, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: ltrcemls  !< litter carbon emission disturbance losses, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: ntchlveg  !< fluxes for each pft: Net change in leaf biomass, u-mol CO2/m2.sec
    real, allocatable, dimension(:,:) :: ntchsveg  !< fluxes for each pft: Net change in stem biomass, u-mol CO2/m2.sec
    real, allocatable, dimension(:,:) :: ntchrveg  !< fluxes for each pft: Net change in root biomass,
    !! the net change is the difference between allocation and
    !! autotrophic respiratory fluxes, u-mol CO2/m2.sec

    real, allocatable, dimension(:)   :: extnprob   !< fire extingusinging probability
    real, allocatable, dimension(:)   :: prbfrhuc   !< probability of fire due to human causes
    real, allocatable, dimension(:)   :: dayl_max   !< maximum daylength for that location (hours)
    real, allocatable, dimension(:)   :: dayl       !< daylength for that location (hours)

    real, allocatable, dimension(:,:) :: bmasveg    !< total (gleaf + stem + root) biomass for each ctem pft, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: cmasvegc   !< total canopy mass for each of the 4 class pfts. recall that
    !< class requires canopy mass as an input,and this is now provided by ctem. \f$kg/m^2\f$.
    real, allocatable, dimension(:,:) :: veghght    !< vegetation height (meters)
    real, allocatable, dimension(:,:) :: rootdpth   !< 99% soil rooting depth (meters)
    !< both veghght & rootdpth can be used as diagnostics to see
    !< how vegetation grows above and below ground, respectively
    real, allocatable, dimension(:)   :: rml        !< leaf maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:)   :: rms        !< stem maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: tltrleaf   !< total leaf litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: blfltrdt !< brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, allocatable, dimension(:,:) :: glfltrdt !< brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, allocatable, dimension(:,:) :: tltrstem   !< total stem litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: tltrroot   !< total root litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: leaflitr   !< leaf litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$). this leaf litter
    !< does not include litter generated due to mortality/fire
    real, allocatable, dimension(:,:) :: roottemp   !< root temperature, k
    real, allocatable, dimension(:,:) :: afrleaf    !< allocation fraction for leaves
    real, allocatable, dimension(:,:) :: afrstem    !< allocation fraction for stem
    real, allocatable, dimension(:,:) :: afrroot    !< allocation fraction for root
    real, allocatable, dimension(:,:) :: wtstatus   !< soil water status used for calculating allocation fractions
    real, allocatable, dimension(:,:) :: ltstatus   !< light status used for calculating allocation fractions
    real, allocatable, dimension(:)   :: rmr        !< root maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)

    real, allocatable, dimension(:,:) :: slopefrac      !< prescribed fraction of wetlands based on slope
    !< only(0.025,0.05,0.1,0.15,0.20,0.25,0.3 and 0.35 percent slope thresholds)
    real, allocatable, dimension(:)   :: wetfrac_pres  !< Prescribed fraction of wetlands in a grid cell
    real, allocatable, dimension(:)   :: ch4WetSpec       !< methane flux from wetlands calculated using hetrores (\f$\mu mol CH_4 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:)   :: wetfdyn       !< dynamic wetland fraction
    real, allocatable, dimension(:)   :: ch4WetDyn       !< methane flux from wetlands calculated using hetrores
    !< and wetfdyn, (\f$\mu mol CH_4 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:)   :: ch4_soills    !< Methane uptake into the soil column (\f$mg CH_4 m^{-2} s^{-1}\f$)

    real, allocatable, dimension(:)   :: lucemcom   !< land use change (luc) related combustion emission losses, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:)   :: lucltrin   !< luc related inputs to litter pool, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:)   :: lucsocin   !< luc related inputs to soil c pool, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)

    real, allocatable, dimension(:)   :: npp        !< net primary productivity
    real, allocatable, dimension(:)   :: nep        !< net ecosystem productivity
    real, allocatable, dimension(:)   :: nbp        !< net biome productivity
    real, allocatable, dimension(:)   :: gpp        !< gross primary productivity
    real, allocatable, dimension(:)   :: hetrores   !< heterotrophic respiration
    real, allocatable, dimension(:)   :: autores    !< autotrophic respiration
    real, allocatable, dimension(:)   :: soilcresp  !<
    real, allocatable, dimension(:)   :: rm         !< maintenance respiration
    real, allocatable, dimension(:)   :: rg         !< growth respiration
    real, allocatable, dimension(:)   :: litres     !< litter respiration
    real, allocatable, dimension(:)   :: socres     !< soil carbon respiration
    real, allocatable, dimension(:)   :: dstcemls   !< carbon emission losses due to disturbance, mainly fire
    real, allocatable, dimension(:)   :: litrfall   !< total litter fall (from leaves, stem, and root) due to
    !< all causes (mortality,turnover,and disturbance)
    real, allocatable, dimension(:)   :: humiftrs   !< transfer of humidified litter from litter to soil c pool

    real, allocatable, dimension(:,:) :: gppveg     !< gross primary productity for each pft
    real, allocatable, dimension(:,:) :: nepveg     !< net ecosystem productity for bare fraction expnbaln(i)=0.0 amount
    !< of c related to spatial expansion Not used JM Jun 2014
    !< OR net ecosystem productity for each pft

    integer, allocatable, dimension(:) :: ipeatland !< Peatland flag: 0 = not a peatland, 1 = bog, 2 = fen
    real, allocatable, dimension(:) :: peatdep      !< Depth of peat column (m)
    real, allocatable, dimension(:) :: anmoss     !< net photosynthetic rate of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: rmlmoss    !< maintenance respiration rate of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: gppmoss    !< gross primaray production of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: nppmoss    !< net primary production of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: armoss     !< autotrophic respiration of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: litrmsmoss !< moss litter mass, \f$kg C/m^2\f$
    real, allocatable, dimension(:) :: Cmossmas   !< C in moss biomass, \f$kg C/m^2\f$
    real, allocatable, dimension(:) :: dmoss      !< depth of living moss (m)
    real, allocatable, dimension(:) :: pdd        !< peatland degree days above 0 deg C.
    real, allocatable, dimension(:) :: ancsmoss   !< moss net photosynthesis in canopy snow subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: angsmoss   !< moss net photosynthesis in snow ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: ancmoss    !< moss net photosynthesis in canopy ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: angmoss    !< moss net photosynthesis in bare ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: rmlcsmoss  !< moss maintenance respiration in canopy snow subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: rmlgsmoss  !< moss maintenance respiration in ground snow subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: rmlcmoss   !< moss maintenance respiration in canopy ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:) :: rmlgmoss   !< moss maintenance respiration in bare ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)

    real, allocatable, dimension(:,:) :: nbpveg     !< net biome productity for bare fraction OR net biome productity for each pft
    real, allocatable, dimension(:,:) :: nppveg     !< npp for individual pfts, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, allocatable, dimension(:,:) :: hetroresveg !<
    real, allocatable, dimension(:,:) :: autoresveg !<
    ! COMBAK PERLAY
    real, allocatable, dimension(:,:) :: litresveg  !<
    real, allocatable, dimension(:,:) :: soilcresveg !<
    real, allocatable, dimension(:,:) :: humiftrsveg !<
    ! real, allocatable, dimension(:,:,:) :: litresveg  !<
    ! real, allocatable, dimension(:,:,:) :: soilcresveg !<
    ! real, allocatable, dimension(:,:,:) :: humiftrsveg !<
    ! COMBAK PERLAY
    real, allocatable, dimension(:,:) :: rmlvegacc  !<
    real, allocatable, dimension(:,:) :: rmsveg     !< stem maintenance resp. rate for each pft
    real, allocatable, dimension(:,:) :: rmrveg     !< root maintenance resp. rate for each pft
    real, allocatable, dimension(:,:) :: rgveg      !< growth resp. rate for each pft
    real, allocatable, dimension(:,:) :: litrfallveg !< litter fall in \f$kg c/m^2\f$ for each pft
    real, allocatable, dimension(:,:) :: reprocost   !< Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real, allocatable, dimension(:,:) :: rothrlos !< root death as crops are harvested, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: pfcancmx !< previous year's fractional coverages of pfts
    real, allocatable, dimension(:,:) :: nfcancmx !< next year's fractional coverages of pfts
    real, allocatable, dimension(:,:) :: alvsctm  !<
    real, allocatable, dimension(:,:) :: paic     !< plant area index for class' 4 pfts. this is the sum of leaf
    !< area index and stem area index.
    real, allocatable, dimension(:,:) :: slaic    !< storage lai. this will be used as min. lai that class sees
    !< so that it doesn't blow up in its stomatal conductance calculations.
    real, allocatable, dimension(:,:) :: alirctm  !<
    real, allocatable, dimension(:)   :: cfluxcg  !<
    real, allocatable, dimension(:)   :: cfluxcs  !<
    real, allocatable, dimension(:)   :: dstcemls3 !< carbon emission losses due to disturbance (fire at present) from litter pool
    real, allocatable, dimension(:,:) :: anveg    !< net photosynthesis rate for each pft
    real, allocatable, dimension(:,:) :: rmlveg   !< leaf maintenance resp. rate for each pft

    real, allocatable, dimension(:) :: twarmm            !< temperature of the warmest month (c)
    real, allocatable, dimension(:) :: tcoldm            !< temperature of the coldest month (c)
    real, allocatable, dimension(:) :: gdd5              !< growing degree days above 5 c
    real, allocatable, dimension(:) :: aridity           !< aridity index, ratio of potential evaporation to precipitation
    real, allocatable, dimension(:) :: srplsmon          !< number of months in a year with surplus water i.e. precipitation more than potential evaporation
    real, allocatable, dimension(:) :: defctmon          !< number of months in a year with water deficit i.e. precipitation less than potential evaporation
    real, allocatable, dimension(:) :: anndefct          !< annual water deficit (mm)
    real, allocatable, dimension(:) :: annsrpls          !< annual water surplus (mm)
    real, allocatable, dimension(:) :: annpcp            !< annual precipitation (mm)
    real, allocatable, dimension(:) :: dry_season_length !< length of dry season (months)

    ! These go into CTEM and are used to keep track of the bioclim limits.
    real, allocatable, dimension(:) :: tcurm     !< temperature of the current month (c)
    real, allocatable, dimension(:) :: srpcuryr  !< water surplus for the current year
    real, allocatable, dimension(:) :: dftcuryr  !< water deficit for the current year
    real, allocatable, dimension(:,:) :: tmonth  !< monthly temperatures
    real, allocatable, dimension(:) :: anpcpcur  !< annual precipitation for current year (mm)
    real, allocatable, dimension(:) :: anpecur   !< annual potential evaporation for current year (mm)
    real, allocatable, dimension(:) :: gdd5cur   !< growing degree days above 5 c for current year
    real, allocatable, dimension(:) :: surmncur  !< number of months with surplus water for current year
    real, allocatable, dimension(:) :: defmncur  !< number of months with water deficit for current year
    real, allocatable, dimension(:) :: srplscur  !< water surplus for the current month
    real, allocatable, dimension(:) :: defctcur  !< water deficit for the current month

    real, allocatable, dimension(:,:) :: geremort !< growth efficiency related mortality (1/day)
    real, allocatable, dimension(:,:) :: intrmort !< intrinsic (age related) mortality (1/day)
    real, allocatable, dimension(:,:) :: lambda   !< Used to determine the colonization rate
    real, allocatable, dimension(:,:) :: cc       !< colonization rate & mortality rate
    real, allocatable, dimension(:,:) :: mm       !< colonization rate & mortality rate

    logical, allocatable, dimension(:,:) :: pftexist !< logical array indicating pfts exist (t) or not (f)
    integer, allocatable, dimension(:,:) :: colddays !< cold days counter for tracking days below a certain
    !< temperature threshold for ndl dcd and crop pfts.
    integer, allocatable, dimension(:,:) :: lfstatus !< leaf phenology status
    integer, allocatable, dimension(:,:) :: pandays  !< days with positive net photosynthesis (an) for use in
    !< the phenology subroutine
    real, allocatable, dimension(:) :: grclarea      !< area of the grid cell, \f$km^2\f$

    integer, allocatable, dimension(:) :: altotcount_ctm ! nlat  !< Counter used for calculating total albedo
    real, allocatable, dimension(:,:)  :: todfrac  !(ilg,icc)   !< Max. fractional coverage of ctem's 9 pfts by the end of the day, for use by land use subroutine
    real, allocatable, dimension(:)    :: fsinacc_gat !(ilg)    !<
    real, allocatable, dimension(:)    :: flutacc_gat !(ilg)    !<
    real, allocatable, dimension(:)    :: flinacc_gat !(ilg)    !<
    ! real, allocatable, dimension(:)    :: pregacc_gat !(ilg)    !<
    real, allocatable, dimension(:)    :: altotacc_gat !(ilg)   !<
    real, allocatable, dimension(:)    :: netrad_gat !(ilg)     !<
    real, allocatable, dimension(:)    :: preacc_gat !(ilg)     !<
    real, allocatable, dimension(:)    :: sdepgat !(ilg)        !<
    !     real, allocatable, dimension(:,:)  :: rgmgat !(ilg,ignd)    !<
    real, allocatable, dimension(:,:)  :: sandgat !(ilg,ignd)   !<
    real, allocatable, dimension(:,:)  :: claygat !(ilg,ignd)   !<
    real, allocatable, dimension(:,:)  :: orgmgat !(ilg,ignd)   !<
    real, allocatable, dimension(:)    :: socigat !(ilg)
    real, allocatable, dimension(:)    :: xdiffusgat !(ilg)
    real, allocatable, dimension(:)    :: faregat !(ilg)

  end type veg_gat

  type (veg_gat), save, target :: vgat

  !=================================================================================
  type tracersType
    !   Simple tracer variables. Only written to if useTracer > 0.

    ! NOTE: Units may vary depending on the tracer used, see convertTracerUnits in tracer.f90

    ! Pools:
    ! allocated with nlat, nmos, ...:
    real, allocatable, dimension(:,:) :: mossCMassrot      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, allocatable, dimension(:,:) :: mossLitrMassrot   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, allocatable, dimension(:,:) :: tracerCO2rot     !< Atmopspheric tracer CO2 concentration (units vary)

    ! allocated with nlat, nmos, icc:
    real, allocatable, dimension(:,:,:) :: gLeafMassrot      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg C/m^2\f$
    real, allocatable, dimension(:,:,:) :: bLeafMassrot      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg C/m^2\f$
    real, allocatable, dimension(:,:,:) :: stemMassrot       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: rootMassrot       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    ! allocated with nlat, nmos, iccp2, ignd:
    real, allocatable, dimension(:,:,:,:) :: litrMassrot       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:,:) :: soilCMassrot      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

    ! allocated with ilg, ...:
    real, allocatable, dimension(:) :: mossCMassgat      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, allocatable, dimension(:) :: mossLitrMassgat   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, allocatable, dimension(:) :: tracerCO2gat     !< Atmopspheric tracer CO2 concentration (units vary)
    ! allocated with nlat, nmos, icc:
    real, allocatable, dimension(:,:) :: gLeafMassgat      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: bLeafMassgat      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: stemMassgat       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:) :: rootMassgat       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    ! allocated with nlat, nmos, iccp2, ignd:
    real, allocatable, dimension(:,:,:) :: litrMassgat       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, allocatable, dimension(:,:,:) :: soilCMassgat      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

  end type tracersType

  type (tracersType), save, target :: tracer

  !=================================================================================
  !> CTEM's variables per tile
  type ctem_tile_level

    !   Tile-level variables (denoted by an ending of "_t")

    real, allocatable, dimension(:) :: fsnowacc_t       !<
    real, allocatable, dimension(:) :: tcansacc_t       !<
    real, allocatable, dimension(:) :: tcanoaccgat_t    !<
    real, allocatable, dimension(:) :: taaccgat_t       !<
    real, allocatable, dimension(:) :: uvaccgat_t       !<
    real, allocatable, dimension(:) :: vvaccgat_t       !<
    real, allocatable, dimension(:) :: anmossac_t       !< daily averaged moss net photosynthesis accumulated (\f$\mu mol /m^2 /s\f$)
    real, allocatable, dimension(:) :: rmlmossac_t      !< daily averaged moss maintainence respiration (\f$\mu mol /m^2 /s\f$)
    real, allocatable, dimension(:) :: gppmossac_t      !< daily averaged gross primary production (\f$\mu mol /m^2 /s\f$)

    ! allocated with ilg, ignd:
    real, allocatable, dimension(:,:) :: tbaraccgat_t !<
    real, allocatable, dimension(:,:) :: thliqacc_t  !<
    real, allocatable, dimension(:,:) :: thiceacc_t  !< Added in place of YW's thicaccgat_m. EC Dec 23 2016.

    ! allocated with ilg, icc:
    real, allocatable, dimension(:,:) :: ancgvgac_t  !<
    real, allocatable, dimension(:,:) :: rmlcgvga_t  !<

  end type ctem_tile_level

  type (ctem_tile_level), save, target :: ctem_tile

  !=================================================================================
  !> CTEM's variables monthly averaged (per pft)
  type ctem_monthly

    !     Tile-level monthly variables (denoted by name ending in "_mo_t")

    ! allocated with nlat, nmos, icc/iccp1/iccp2:

    real, allocatable, dimension(:,:,:) :: laimaxg_mo    !<
    real, allocatable, dimension(:,:,:) :: stemmass_mo   !<
    real, allocatable, dimension(:,:,:) :: leafmass_mo   !<
    real, allocatable, dimension(:,:,:) :: rootmass_mo   !<
    real, allocatable, dimension(:,:,:) :: litrfallveg_mo !<
    real, allocatable, dimension(:,:,:) :: humiftrsveg_mo !<
    real, allocatable, dimension(:,:,:) :: npp_mo        !<
    real, allocatable, dimension(:,:,:) :: gpp_mo        !<
    real, allocatable, dimension(:,:,:) :: vgbiomas_mo   !<
    real, allocatable, dimension(:,:,:) :: autores_mo    !<
    real, allocatable, dimension(:,:,:) :: soilres_mo    !<
    real, allocatable, dimension(:,:,:) :: totcmass_mo   !<
    real, allocatable, dimension(:,:,:) :: nep_mo        !<
    real, allocatable, dimension(:,:,:) :: hetrores_mo   !<
    real, allocatable, dimension(:,:,:) :: nbp_mo        !<
    real, allocatable, dimension(:,:,:) :: emit_co2_mo  !<
    real, allocatable, dimension(:,:,:) :: emit_co_mo   !<
    real, allocatable, dimension(:,:,:) :: emit_ch4_mo  !<
    real, allocatable, dimension(:,:,:) :: emit_nmhc_mo !<
    real, allocatable, dimension(:,:,:) :: emit_h2_mo   !<
    real, allocatable, dimension(:,:,:) :: emit_nox_mo  !<
    real, allocatable, dimension(:,:,:) :: emit_n2o_mo  !<
    real, allocatable, dimension(:,:,:) :: emit_pm25_mo !<
    real, allocatable, dimension(:,:,:) :: emit_tpm_mo  !<
    real, allocatable, dimension(:,:,:) :: emit_tc_mo   !<
    real, allocatable, dimension(:,:,:) :: emit_oc_mo   !<
    real, allocatable, dimension(:,:,:) :: emit_bc_mo   !<
    real, allocatable, dimension(:,:,:) :: burnfrac_mo  !<
    real, allocatable, dimension(:,:,:) :: bterm_mo     !<
    real, allocatable, dimension(:,:,:) :: mterm_mo     !<
    real, allocatable, dimension(:,:,:) :: smfuncveg_mo !<
    ! allocated with nlat, nmos, icc/iccp1/iccp2, ignd:
    ! COMBAK PERLAY
    real, allocatable, dimension(:,:,:) :: litrmass_mo   !<
    real, allocatable, dimension(:,:,:) :: soilcmas_mo   !<
    real, allocatable, dimension(:,:,:) :: litres_mo     !<
    real, allocatable, dimension(:,:,:) :: soilcres_mo   !<
    ! real, allocatable, dimension(:,:,:,:) :: litrmass_mo   !<
    ! real, allocatable, dimension(:,:,:,:) :: soilcmas_mo   !<
    ! real, allocatable, dimension(:,:,:,:) :: litres_mo     !<
    ! real, allocatable, dimension(:,:,:,:) :: soilcres_mo   !<
    ! COMBAK PERLAY

  end type ctem_monthly

  type (ctem_monthly), save, target :: ctem_mo

  !=================================================================================
  !> CTEM's grid average monthly values
  type ctem_gridavg_monthly

    !  Grid averaged monthly variables (denoted by name ending in "_mo_g")

    ! allocated with nlat:
    real, allocatable, dimension(:) :: laimaxg_mo_g  !<
    real, allocatable, dimension(:) :: stemmass_mo_g !<
    real, allocatable, dimension(:) :: leafmass_mo_g !<
    real, allocatable, dimension(:) :: rootmass_mo_g !<
    real, allocatable, dimension(:) :: litrfall_mo_g !<
    real, allocatable, dimension(:) :: humiftrs_mo_g !<
    real, allocatable, dimension(:) :: npp_mo_g      !<
    real, allocatable, dimension(:) :: gpp_mo_g      !<
    real, allocatable, dimension(:) :: nep_mo_g      !<
    real, allocatable, dimension(:) :: nbp_mo_g      !<
    real, allocatable, dimension(:) :: hetrores_mo_g !<
    real, allocatable, dimension(:) :: autores_mo_g  !<
    real, allocatable, dimension(:) :: soilres_mo_g  !<
    real, allocatable, dimension(:) :: vgbiomas_mo_g !<
    real, allocatable, dimension(:) :: totcmass_mo_g !<
    real, allocatable, dimension(:) :: emit_co2_mo_g !<
    real, allocatable, dimension(:) :: emit_co_mo_g  !<
    real, allocatable, dimension(:) :: emit_ch4_mo_g !<
    real, allocatable, dimension(:) :: emit_nmhc_mo_g !<
    real, allocatable, dimension(:) :: emit_h2_mo_g  !<
    real, allocatable, dimension(:) :: emit_nox_mo_g !<
    real, allocatable, dimension(:) :: emit_n2o_mo_g !<
    real, allocatable, dimension(:) :: emit_pm25_mo_g !<
    real, allocatable, dimension(:) :: emit_tpm_mo_g !<
    real, allocatable, dimension(:) :: emit_tc_mo_g  !<
    real, allocatable, dimension(:) :: emit_oc_mo_g  !<
    real, allocatable, dimension(:) :: emit_bc_mo_g  !<
    real, allocatable, dimension(:) :: smfuncveg_mo_g !<
    real, allocatable, dimension(:) :: luc_emc_mo_g  !<
    real, allocatable, dimension(:) :: lucltrin_mo_g !<
    real, allocatable, dimension(:) :: lucsocin_mo_g !<
    real, allocatable, dimension(:) :: burnfrac_mo_g !<
    real, allocatable, dimension(:) :: bterm_mo_g    !<
    real, allocatable, dimension(:) :: lterm_mo_g    !<
    real, allocatable, dimension(:) :: mterm_mo_g    !<
    real, allocatable, dimension(:) :: ch4WetSpec_mo_g  !<
    real, allocatable, dimension(:) :: wetfdyn_mo_g  !<
    real, allocatable, dimension(:) :: wetfpres_mo_g  !<
    real, allocatable, dimension(:) :: ch4WetDyn_mo_g  !<
    real, allocatable, dimension(:) :: ch4soills_mo_g !<
    real, allocatable, dimension(:) :: cProduct_mo_g          !< Carbon in the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2}]\f$
    real, allocatable, dimension(:) :: fProductDecomp_mo_g    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$

    ! allocated with nlat, ignd:
    ! COMBAK PERLAY
    real, allocatable, dimension(:) :: litrmass_mo_g !<
    real, allocatable, dimension(:) :: soilcmas_mo_g !<
    real, allocatable, dimension(:) :: litres_mo_g   !<
    real, allocatable, dimension(:) :: soilcres_mo_g !<
    ! real, allocatable, dimension(:,:) :: litrmass_mo_g !<
    ! real, allocatable, dimension(:,:) :: soilcmas_mo_g !<
    ! real, allocatable, dimension(:,:) :: litres_mo_g   !<
    ! real, allocatable, dimension(:,:) :: soilcres_mo_g !<
    ! COMBAK PERLAY

  end type ctem_gridavg_monthly

  type (ctem_gridavg_monthly), save, target :: ctem_grd_mo

  !=================================================================================
  !> CTEM's variables per tile monthly values
  type ctem_tileavg_monthly

    !     Tile-level monthly variables (denoted by name ending in "_mo_t")

    ! allocated with nlat, nmos:

    real, allocatable, dimension(:,:) :: laimaxg_mo_t  !<
    real, allocatable, dimension(:,:) :: stemmass_mo_t !<
    real, allocatable, dimension(:,:) :: leafmass_mo_t !<
    real, allocatable, dimension(:,:) :: rootmass_mo_t !<
    real, allocatable, dimension(:,:) :: litrfall_mo_t !<
    real, allocatable, dimension(:,:) :: humiftrs_mo_t !<
    real, allocatable, dimension(:,:) :: npp_mo_t      !<
    real, allocatable, dimension(:,:) :: gpp_mo_t      !<
    real, allocatable, dimension(:,:) :: vgbiomas_mo_t !<
    real, allocatable, dimension(:,:) :: autores_mo_t  !<
    real, allocatable, dimension(:,:) :: soilres_mo_t  !<
    real, allocatable, dimension(:,:) :: totcmass_mo_t !<
    real, allocatable, dimension(:,:) :: nep_mo_t      !<
    real, allocatable, dimension(:,:) :: hetrores_mo_t !<
    real, allocatable, dimension(:,:) :: nbp_mo_t      !<
    real, allocatable, dimension(:,:) :: emit_co2_mo_t !<
    real, allocatable, dimension(:,:) :: emit_co_mo_t  !<
    real, allocatable, dimension(:,:) :: emit_ch4_mo_t !<
    real, allocatable, dimension(:,:) :: emit_nmhc_mo_t !<
    real, allocatable, dimension(:,:) :: emit_h2_mo_t  !<
    real, allocatable, dimension(:,:) :: emit_nox_mo_t !<
    real, allocatable, dimension(:,:) :: emit_n2o_mo_t !<
    real, allocatable, dimension(:,:) :: emit_pm25_mo_t !<
    real, allocatable, dimension(:,:) :: emit_tpm_mo_t !<
    real, allocatable, dimension(:,:) :: emit_tc_mo_t  !<
    real, allocatable, dimension(:,:) :: emit_oc_mo_t  !<
    real, allocatable, dimension(:,:) :: emit_bc_mo_t  !<
    real, allocatable, dimension(:,:) :: burnfrac_mo_t !<
    real, allocatable, dimension(:,:) :: smfuncveg_mo_t !<
    real, allocatable, dimension(:,:) :: bterm_mo_t    !<
    real, allocatable, dimension(:,:) :: luc_emc_mo_t  !<
    real, allocatable, dimension(:,:) :: lterm_mo_t    !<
    real, allocatable, dimension(:,:) :: lucsocin_mo_t !<
    real, allocatable, dimension(:,:) :: mterm_mo_t    !<
    real, allocatable, dimension(:,:) :: lucltrin_mo_t !<
    real, allocatable, dimension(:,:) :: ch4WetSpec_mo_t  !<
    real, allocatable, dimension(:,:) :: wetfdyn_mo_t  !<
    real, allocatable, dimension(:,:) :: wetfpres_mo_t  !<
    real, allocatable, dimension(:,:) :: ch4WetDyn_mo_t  !<
    real, allocatable, dimension(:,:) :: ch4soills_mo_t !<
    real, allocatable, dimension(:,:) :: wind_mo_t     !<
    real, allocatable, dimension(:,:) :: fProductDecomp_mo_t    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$

    ! allocated with nlat, nmos, ignd:
    ! COMBAK PERLAY
    real, allocatable, dimension(:,:) :: litrmass_mo_t !<
    real, allocatable, dimension(:,:) :: soilcmas_mo_t !<
    real, allocatable, dimension(:,:) :: litres_mo_t   !<
    real, allocatable, dimension(:,:) :: soilcres_mo_t !<
    ! real, allocatable, dimension(:,:,:) :: litrmass_mo_t !<
    ! real, allocatable, dimension(:,:,:) :: soilcmas_mo_t !<
    ! real, allocatable, dimension(:,:,:) :: litres_mo_t   !<
    ! real, allocatable, dimension(:,:,:) :: soilcres_mo_t !<
    ! COMBAK PERLAY

  end type ctem_tileavg_monthly

  type (ctem_tileavg_monthly), save, target :: ctem_tile_mo


  !=================================================================================
  !> CTEM's average annual values (per PFT)
  type ctem_annual

    ! c      Annual output for CTEM mosaic variables:
    ! c      (denoted by name ending in "_yr_m")
    !
    ! allocated with nlat, nmos, icc/iccp1/iccp2:

    real, allocatable, dimension(:,:,:) :: laimaxg_yr   !<
    real, allocatable, dimension(:,:,:) :: stemmass_yr  !<
    real, allocatable, dimension(:,:,:) :: rootmass_yr  !<
    real, allocatable, dimension(:,:,:) :: npp_yr       !<
    real, allocatable, dimension(:,:,:) :: gpp_yr       !<
    real, allocatable, dimension(:,:,:) :: vgbiomas_yr  !<
    real, allocatable, dimension(:,:,:) :: autores_yr   !<
    real, allocatable, dimension(:,:,:) :: totcmass_yr !<
    real, allocatable, dimension(:,:,:) :: nep_yr     !<
    real, allocatable, dimension(:,:,:) :: hetrores_yr !<
    real, allocatable, dimension(:,:,:) :: nbp_yr     !<
    real, allocatable, dimension(:,:,:) :: emit_co2_yr  !<
    real, allocatable, dimension(:,:,:) :: emit_co_yr   !<
    real, allocatable, dimension(:,:,:) :: emit_ch4_yr  !<
    real, allocatable, dimension(:,:,:) :: emit_nmhc_yr !<
    real, allocatable, dimension(:,:,:) :: emit_h2_yr   !<
    real, allocatable, dimension(:,:,:) :: emit_nox_yr  !<
    real, allocatable, dimension(:,:,:) :: emit_n2o_yr  !<
    real, allocatable, dimension(:,:,:) :: emit_pm25_yr !<
    real, allocatable, dimension(:,:,:) :: emit_tpm_yr  !<
    real, allocatable, dimension(:,:,:) :: emit_tc_yr   !<
    real, allocatable, dimension(:,:,:) :: emit_oc_yr   !<
    real, allocatable, dimension(:,:,:) :: emit_bc_yr   !<
    real, allocatable, dimension(:,:,:) :: bterm_yr     !<
    real, allocatable, dimension(:,:,:) :: mterm_yr     !<
    real, allocatable, dimension(:,:,:) :: burnfrac_yr  !<
    real, allocatable, dimension(:,:,:) :: smfuncveg_yr !<
    real, allocatable, dimension(:,:,:) :: veghght_yr   !<

    ! allocated with nlat, nmos, iccp2, ignd:
    ! COMBAK PERLAY
    real, allocatable, dimension(:,:,:) :: litrmass_yr !<
    real, allocatable, dimension(:,:,:) :: soilcmas_yr !<
    real, allocatable, dimension(:,:,:) :: litres_yr  !<
    real, allocatable, dimension(:,:,:) :: soilcres_yr !<
    ! real, allocatable, dimension(:,:,:,:) :: litrmass_yr !<
    ! real, allocatable, dimension(:,:,:,:) :: soilcmas_yr !<
    ! real, allocatable, dimension(:,:,:,:) :: litres_yr  !<
    ! real, allocatable, dimension(:,:,:,:) :: soilcres_yr !<
    ! COMBAK PERLAY

  end type ctem_annual

  type (ctem_annual), save, target :: ctem_yr

  !=================================================================================
  !> CTEM's grid average annual values
  type ctem_gridavg_annual

    ! Annual output for CTEM grid-averaged variables:
    ! (denoted by name ending in "_yr_g")

    ! allocated with nlat:
    real, allocatable, dimension(:) :: laimaxg_yr_g  !<
    real, allocatable, dimension(:) :: stemmass_yr_g !<
    real, allocatable, dimension(:) :: rootmass_yr_g !<
    real, allocatable, dimension(:) :: npp_yr_g      !<
    real, allocatable, dimension(:) :: gpp_yr_g      !<
    real, allocatable, dimension(:) :: nep_yr_g      !<
    real, allocatable, dimension(:) :: nbp_yr_g      !<
    real, allocatable, dimension(:) :: hetrores_yr_g !<
    real, allocatable, dimension(:) :: autores_yr_g  !<
    real, allocatable, dimension(:) :: vgbiomas_yr_g !<
    real, allocatable, dimension(:) :: totcmass_yr_g !<
    real, allocatable, dimension(:) :: emit_co2_yr_g !<
    real, allocatable, dimension(:) :: emit_co_yr_g  !<
    real, allocatable, dimension(:) :: emit_ch4_yr_g !<
    real, allocatable, dimension(:) :: emit_nmhc_yr_g !<
    real, allocatable, dimension(:) :: emit_h2_yr_g  !<
    real, allocatable, dimension(:) :: emit_nox_yr_g !<
    real, allocatable, dimension(:) :: emit_n2o_yr_g !<
    real, allocatable, dimension(:) :: emit_pm25_yr_g !<
    real, allocatable, dimension(:) :: emit_tpm_yr_g !<
    real, allocatable, dimension(:) :: emit_tc_yr_g  !<
    real, allocatable, dimension(:) :: emit_oc_yr_g  !<
    real, allocatable, dimension(:) :: emit_bc_yr_g  !<
    real, allocatable, dimension(:) :: smfuncveg_yr_g !<
    real, allocatable, dimension(:) :: luc_emc_yr_g  !<
    real, allocatable, dimension(:) :: lucltrin_yr_g !<
    real, allocatable, dimension(:) :: lucsocin_yr_g !<
    real, allocatable, dimension(:) :: burnfrac_yr_g !<
    real, allocatable, dimension(:) :: bterm_yr_g    !<
    real, allocatable, dimension(:) :: lterm_yr_g    !<
    real, allocatable, dimension(:) :: mterm_yr_g    !<
    real, allocatable, dimension(:) :: ch4WetSpec_yr_g  !<
    real, allocatable, dimension(:) :: wetfdyn_yr_g  !<
    real, allocatable, dimension(:) :: ch4WetDyn_yr_g  !<
    real, allocatable, dimension(:) :: ch4soills_yr_g !<
    real, allocatable, dimension(:) :: veghght_yr_g  !<
    real, allocatable, dimension(:) :: peatdep_yr_g  !<
    real, allocatable, dimension(:) :: cProduct_yr_g          !< Carbon in the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2}]\f$
    real, allocatable, dimension(:) :: fProductDecomp_yr_g    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$

    ! allocated with nlat, ignd:
    ! COMBAK PERLAY
    real, allocatable, dimension(:) :: litrmass_yr_g !<
    real, allocatable, dimension(:) :: soilcmas_yr_g !<
    real, allocatable, dimension(:) :: litres_yr_g   !<
    real, allocatable, dimension(:) :: soilcres_yr_g !<
    ! real, allocatable, dimension(:,:) :: litrmass_yr_g !<
    ! real, allocatable, dimension(:,:) :: soilcmas_yr_g !<
    ! real, allocatable, dimension(:,:) :: litres_yr_g   !<
    ! real, allocatable, dimension(:,:) :: soilcres_yr_g !<
    ! COMBAK PERLAY

  end type ctem_gridavg_annual

  type (ctem_gridavg_annual), save, target :: ctem_grd_yr

  !=================================================================================
  !> CTEM's variables per tile annual values
  type ctem_tileavg_annual

    ! c      Annual output for CTEM mosaic variables:
    ! c      (denoted by name ending in "_yr_m")
    !
    ! allocated with nlat, nmos:
    real, allocatable, dimension(:,:) :: laimaxg_yr_t  !<
    real, allocatable, dimension(:,:) :: stemmass_yr_t !<
    real, allocatable, dimension(:,:) :: rootmass_yr_t !<
    real, allocatable, dimension(:,:) :: npp_yr_t      !<
    real, allocatable, dimension(:,:) :: gpp_yr_t      !<
    real, allocatable, dimension(:,:) :: vgbiomas_yr_t !<
    real, allocatable, dimension(:,:) :: autores_yr_t  !<
    real, allocatable, dimension(:,:) :: totcmass_yr_t !<
    real, allocatable, dimension(:,:) :: nep_yr_t      !<
    real, allocatable, dimension(:,:) :: hetrores_yr_t !<
    real, allocatable, dimension(:,:) :: nbp_yr_t      !<
    real, allocatable, dimension(:,:) :: emit_co2_yr_t !<
    real, allocatable, dimension(:,:) :: emit_co_yr_t  !<
    real, allocatable, dimension(:,:) :: emit_ch4_yr_t !<
    real, allocatable, dimension(:,:) :: emit_nmhc_yr_t !<
    real, allocatable, dimension(:,:) :: emit_h2_yr_t  !<
    real, allocatable, dimension(:,:) :: emit_nox_yr_t !<
    real, allocatable, dimension(:,:) :: emit_n2o_yr_t !<
    real, allocatable, dimension(:,:) :: emit_pm25_yr_t !<
    real, allocatable, dimension(:,:) :: emit_tpm_yr_t !<
    real, allocatable, dimension(:,:) :: emit_tc_yr_t  !<
    real, allocatable, dimension(:,:) :: emit_oc_yr_t  !<
    real, allocatable, dimension(:,:) :: emit_bc_yr_t  !<
    real, allocatable, dimension(:,:) :: burnfrac_yr_t !<
    real, allocatable, dimension(:,:) :: smfuncveg_yr_t !<
    real, allocatable, dimension(:,:) :: bterm_yr_t    !<
    real, allocatable, dimension(:,:) :: luc_emc_yr_t  !<
    real, allocatable, dimension(:,:) :: lterm_yr_t    !<
    real, allocatable, dimension(:,:) :: lucsocin_yr_t !<
    real, allocatable, dimension(:,:) :: mterm_yr_t    !<
    real, allocatable, dimension(:,:) :: lucltrin_yr_t !<
    real, allocatable, dimension(:,:) :: ch4WetSpec_yr_t  !<
    real, allocatable, dimension(:,:) :: wetfdyn_yr_t  !<
    real, allocatable, dimension(:,:) :: ch4WetDyn_yr_t  !<
    real, allocatable, dimension(:,:) :: ch4soills_yr_t !<
    real, allocatable, dimension(:,:) :: veghght_yr_t  !<
    real, allocatable, dimension(:,:) :: peatdep_yr_t  !<
    real, allocatable, dimension(:,:) :: fProductDecomp_yr_t    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$

    ! allocated with nlat, nmos, ignd:
    ! COMBAK PERLAY
    real, allocatable, dimension(:,:) :: litrmass_yr_t !<
    real, allocatable, dimension(:,:) :: soilcmas_yr_t !<
    real, allocatable, dimension(:,:) :: litres_yr_t   !<
    real, allocatable, dimension(:,:) :: soilcres_yr_t !<
    ! real, allocatable, dimension(:,:,:) :: litrmass_yr_t !<
    ! real, allocatable, dimension(:,:,:) :: soilcmas_yr_t !<
    ! real, allocatable, dimension(:,:,:) :: litres_yr_t   !<
    ! real, allocatable, dimension(:,:,:) :: soilcres_yr_t !<
    ! COMBAK PERLAY

  end type ctem_tileavg_annual

  type (ctem_tileavg_annual), save, target :: ctem_tile_yr

contains

  !=================================================================================
  !> \ingroup ctemstatevars_allocCtemVars
  !! @{
  !> Allocate the biogeochemistry (CTEM) variables in preparation for a simulation
  subroutine allocCtemVars

    use classicParams,    only : ican, icc, iccp2, iccp1, ilg, nlat, nmos, ignd

    implicit none

    !-----------

    ! allocated with nlat, nmos, icc:
	!found in model_variables
	!nlat = vs%active_cell_count
!	nlat = vs%cell_count
!	nmos = vs%gru_count
!	print *, 'GRU_count ', vs%gru_count
!	print *, 'cell_count ', vs%cell_count
!	print *, 'landtile count ', vs%landtile_count
!	print *, 'active cell count ', vs%active_cell_count
	!ilg = vs%landtile_count
!	ilg = nlat*nmos

    allocate(vrot%pftexist(nlat,nmos,icc), &
             vrot%lfstatus(nlat,nmos,icc), &
             vrot%pandays (nlat,nmos,icc), &
             vrot%gleafmas(nlat,nmos,icc), &
             vrot%bleafmas(nlat,nmos,icc), &
             vrot%stemmass(nlat,nmos,icc), &
             vrot%rootmass(nlat,nmos,icc), &
             vrot%pstemmass (nlat,nmos,icc), &
             vrot%pgleafmass (nlat,nmos,icc), &
             vrot%fcancmx (nlat,nmos,icc), &
             vrot%ailcg   (nlat,nmos,icc), &
             vrot%ailcgs  (nlat,nmos,icc), &
             vrot%fcancs  (nlat,nmos,icc), &
             vrot%fcanc   (nlat,nmos,icc), &
             vrot%co2i1cg (nlat,nmos,icc), &
             vrot%co2i1cs (nlat,nmos,icc), &
             vrot%co2i2cg (nlat,nmos,icc), &
             vrot%co2i2cs (nlat,nmos,icc), &
             vrot%ancsveg (nlat,nmos,icc), &
             vrot%ancgveg (nlat,nmos,icc), &
             vrot%rmlcsveg(nlat,nmos,icc), &
             vrot%rmlcgveg(nlat,nmos,icc), &
             vrot%slai    (nlat,nmos,icc), &
             vrot%ailcb   (nlat,nmos,icc), &
             vrot%flhrloss(nlat,nmos,icc), &
             vrot%grwtheff(nlat,nmos,icc), &
             vrot%lystmmas(nlat,nmos,icc), &
             vrot%lyrotmas(nlat,nmos,icc), &
             vrot%tymaxlai(nlat,nmos,icc), &
             vrot%stmhrlos(nlat,nmos,icc), &
             vrot%vgbiomas_veg(nlat,nmos,icc), &
             vrot%emit_co2(nlat,nmos,icc), &
             vrot%emit_co (nlat,nmos,icc), &
             vrot%emit_ch4(nlat,nmos,icc), &
             vrot%emit_nmhc(nlat,nmos,icc), &
             vrot%emit_h2 (nlat,nmos,icc), &
             vrot%emit_nox(nlat,nmos,icc), &
             vrot%emit_n2o(nlat,nmos,icc), &
             vrot%emit_pm25(nlat,nmos,icc), &
             vrot%emit_tpm(nlat,nmos,icc), &
             vrot%emit_tc (nlat,nmos,icc), &
             vrot%emit_oc (nlat,nmos,icc), &
             vrot%emit_bc (nlat,nmos,icc), &
             vrot%burnvegf(nlat,nmos,icc), &
             vrot%smfuncveg(nlat,nmos,icc), &
             vrot%bterm   (nlat,nmos,icc), &
             vrot%mterm   (nlat,nmos,icc), &
             vrot%bmasveg (nlat,nmos,icc), &
             vrot%veghght (nlat,nmos,icc), &
             vrot%rootdpth(nlat,nmos,icc), &
             vrot%tltrleaf(nlat,nmos,icc), &
             vrot%tltrstem(nlat,nmos,icc), &
             vrot%tltrroot(nlat,nmos,icc), &
             vrot%leaflitr(nlat,nmos,icc), &
             vrot%roottemp(nlat,nmos,icc), &
             vrot%afrleaf (nlat,nmos,icc), &
             vrot%afrstem (nlat,nmos,icc), &
             vrot%afrroot (nlat,nmos,icc), &
             vrot%wtstatus(nlat,nmos,icc), &
             vrot%ltstatus(nlat,nmos,icc), &
             vrot%gppveg  (nlat,nmos,icc), &
             vrot%nppveg  (nlat,nmos,icc), &
             vrot%autoresveg(nlat,nmos,icc), &
             vrot%rmlvegacc (nlat,nmos,icc), &
             vrot%rmsveg  (nlat,nmos,icc), &
             vrot%rmrveg  (nlat,nmos,icc), &
             vrot%rgveg   (nlat,nmos,icc), &
             vrot%litrfallveg(nlat,nmos,icc), &
             vrot%rothrlos(nlat,nmos,icc), &
             vrot%pfcancmx(nlat,nmos,icc), &
             vrot%nfcancmx(nlat,nmos,icc), &
             vrot%anveg   (nlat,nmos,icc), &
             vrot%rmlveg  (nlat,nmos,icc), &

             tracer%gLeafMassrot(nlat,nmos,icc), &
             tracer%bLeafMassrot(nlat,nmos,icc), &
             tracer%stemMassrot(nlat,nmos,icc), &
             tracer%rootMassrot(nlat,nmos,icc), &

             ! allocated with nlat,nmos:
             vrot%gavglai (nlat,nmos), &
             vrot%co2conc (nlat,nmos), &
             vrot%ch4conc (nlat,nmos), &
             vrot%canres  (nlat,nmos), &
             vrot%vgbiomas(nlat,nmos), &
             vrot%gavgltms(nlat,nmos), &
             vrot%gavgscms(nlat,nmos), &
             vrot%burnfrac(nlat,nmos), &
             vrot%popdin  (nlat,nmos), &
             vrot%lterm   (nlat,nmos), &
             vrot%extnprob(nlat,nmos), &
             vrot%prbfrhuc(nlat,nmos), &
             vrot%rml     (nlat,nmos), &
             vrot%rms     (nlat,nmos), &
             vrot%rmr     (nlat,nmos), &
             vrot%ch4WetSpec (nlat,nmos), &
             vrot%wetfdyn (nlat,nmos), &
             vrot%wetfrac_pres(nlat,nmos), &
             vrot%ch4WetDyn (nlat,nmos), &
             vrot%ch4_soills(nlat,nmos), &
             vrot%lucemcom(nlat,nmos), &
             vrot%lucltrin(nlat,nmos), &
             vrot%lucsocin(nlat,nmos), &
             vrot%npp     (nlat,nmos), &
             vrot%nep     (nlat,nmos), &
             vrot%nbp     (nlat,nmos), &
             vrot%gpp     (nlat,nmos), &
             vrot%hetrores(nlat,nmos), &
             vrot%autores (nlat,nmos), &
             vrot%soilcresp(nlat,nmos), &
             vrot%rm      (nlat,nmos), &
             vrot%rg      (nlat,nmos), &
             vrot%litres  (nlat,nmos), &
             vrot%socres  (nlat,nmos), &
             vrot%dstcemls(nlat,nmos), &
             vrot%litrfall(nlat,nmos), &
             vrot%humiftrs(nlat,nmos), &
             vrot%cfluxcg (nlat,nmos), &
             vrot%cfluxcs (nlat,nmos), &
             vrot%dstcemls3 (nlat,nmos), &
             vrot%tcanrs  (nlat,nmos), &
             vrot%tsnors  (nlat,nmos), &
             vrot%tpndrs  (nlat,nmos), &
             vrot%uvaccrow_m(nlat,nmos), &
             vrot%vvaccrow_m (nlat,nmos), &
             vrot%qevpacc_m_save(nlat,nmos), &
             vrot%twarmm  (nlat,nmos), &
             vrot%tcoldm  (nlat,nmos), &
             vrot%gdd5    (nlat,nmos), &
             vrot%aridity (nlat,nmos), &
             vrot%srplsmon(nlat,nmos), &
             vrot%defctmon(nlat,nmos), &
             vrot%anndefct(nlat,nmos), &
             vrot%annsrpls(nlat,nmos), &
             vrot%annpcp  (nlat,nmos), &
             vrot%dry_season_length(nlat,nmos), &
             vrot%ipeatland(nlat,nmos), &
             vrot%litrmsmoss(nlat,nmos), &
             vrot%Cmossmas(nlat,nmos), &
             vrot%dmoss(nlat,nmos), &
             vrot%nppmoss(nlat,nmos), &
             vrot%rmlmoss(nlat,nmos), &
             vrot%gppmoss(nlat,nmos), &
             vrot%anmoss(nlat,nmos), &
             vrot%armoss(nlat,nmos), &
             vrot%peatdep(nlat,nmos), &
             vrot%pdd(nlat,nmos), &

             tracer%tracerCO2rot(nlat,nmos), &
             tracer%mossCMassrot(nlat,nmos), &
             tracer%mossLitrMassrot(nlat,nmos), &

    ! allocated with nlat, nmos, ican:
             vrot%zolnc(nlat,nmos,ican), &
             vrot%ailc(nlat,nmos,ican), &
             vrot%cmasvegc(nlat,nmos,ican), &
             vrot%alvsctm(nlat,nmos,ican), &
             vrot%paic(nlat,nmos,ican), &
             vrot%slaic(nlat,nmos,ican), &
             vrot%alirctm(nlat,nmos,ican), &
             vrot%csum(nlat,nmos,ican), &

    ! allocated with nlat, nmos, ican, ignd:
             vrot%rmatc(nlat,nmos,ican,ignd), &

    ! allocated with nlat, nmos, icc, ignd:
             vrot%rmatctem(nlat,nmos,icc,ignd), &

    ! allocated with nlat, nmos, iccp1:
             vrot%nepveg(nlat,nmos,iccp1), &
             vrot%nbpveg(nlat,nmos,iccp1), &
             vrot%hetroresveg(nlat,nmos,iccp1), &

    ! allocated with nlat, nmos, iccp2, ignd:
    ! COMBAK PERLAY
             vrot%litrmass(nlat,nmos,iccp2), &
             vrot%soilcmas(nlat,nmos,iccp2), &
             vrot%litresveg(nlat,nmos,iccp2), &
             vrot%soilcresveg(nlat,nmos,iccp2), &
             vrot%humiftrsveg(nlat,nmos,iccp2), &
             ! vrot%litrmass(nlat,nmos,iccp2,ignd), &
             ! vrot%soilcmas(nlat,nmos,iccp2,ignd), &
             ! vrot%litresveg(nlat,nmos,iccp2,ignd), &
             ! vrot%soilcresveg(nlat,nmos,iccp2,ignd), &
             ! vrot%humiftrsveg(nlat,nmos,iccp2,ignd), &
             ! COMBAK PERLAY
             tracer%litrMassrot(nlat,nmos,iccp2,ignd), &
             tracer%soilCMassrot(nlat,nmos,iccp2,ignd), &

             ! allocated with nlat, nmos, {some number}:
             vrot%colddays(nlat,nmos,2), &
             vrot%slopefrac(nlat,nmos,8), &

             ! allocated with nlat:
             vrot%dayl_max(nlat), &
             vrot%dayl(nlat), &
             vgat%altotcount_ctm(nlat))

    ! Now on to the veg_gat vars

    !ilg = nlat * nmos

    ! allocated with ilg

    allocate(vgat%grclarea(ilg), &
             vgat%gavglai (ilg), &
             vgat%lightng (ilg), &
             vgat%co2conc (ilg), &
             vgat%ch4conc (ilg), &
             vgat%canres (ilg), &
             vgat%vgbiomas (ilg), &
             vgat%gavgltms (ilg), &
             vgat%gavgscms (ilg), &
             vgat%lterm (ilg), &
             vgat%ch4WetSpec (ilg), &
             vgat%wetfdyn (ilg), &
             vgat%ch4WetDyn (ilg), &

             vgat%ch4_soills (ilg), &
             vgat%lucemcom (ilg), &
             vgat%lucltrin (ilg), &
             vgat%lucsocin (ilg), &
             vgat%npp (ilg), &
             vgat%nep (ilg), &
             vgat%nbp (ilg), &
             vgat%gpp (ilg), &
             vgat%hetrores (ilg), &
             vgat%autores (ilg), &
             vgat%soilcresp (ilg), &
             vgat%rm (ilg), &
             vgat%rg (ilg), &
             vgat%litres (ilg), &
             vgat%socres (ilg), &
             vgat%dstcemls (ilg), &
             vgat%litrfall (ilg), &
             vgat%humiftrs (ilg), &
             vgat%rml (ilg), &
             vgat%rms (ilg), &
             vgat%rmr (ilg), &
             vgat%burnfrac (ilg), &
             vgat%popdin (ilg), &
             vgat%extnprob (ilg), &
             vgat%prbfrhuc (ilg), &
             vgat%dayl_max (ilg), &
             vgat%dayl (ilg), &
             vgat%wetfrac_pres (ilg), & 
             vgat%twarmm (ilg), &
             vgat%tcoldm (ilg), &
             vgat%gdd5 (ilg), &
             vgat%aridity (ilg), &
             vgat%srplsmon (ilg), &
             vgat%defctmon (ilg), &
             vgat%anndefct (ilg), &
             vgat%annsrpls (ilg), &
             vgat%annpcp (ilg), &
             vgat%dry_season_length (ilg), &
             vgat%cfluxcg (ilg), &
             vgat%cfluxcs (ilg), &
             vgat%dstcemls3 (ilg), &
             vgat%tcurm (ilg), &
             vgat%srpcuryr (ilg), &
             vgat%dftcuryr (ilg), &
             vgat%anpcpcur (ilg), &
             vgat%anpecur (ilg), &
             vgat%gdd5cur (ilg), &
             vgat%surmncur (ilg), &
             vgat%defmncur (ilg), &
             vgat%srplscur (ilg), &
             vgat%defctcur (ilg), &
             vgat%ipeatland (ilg), &
             vgat%anmoss (ilg), &
             vgat%rmlmoss (ilg), &
             vgat%gppmoss (ilg), &
             vgat%nppmoss (ilg), &
             vgat%armoss  (ilg), &
             vgat%litrmsmoss (ilg), &
             vgat%Cmossmas (ilg), &
             vgat%dmoss (ilg), &
             vgat%pdd   (ilg), &
             vgat%peatdep(ilg), &
             vgat%ancsmoss (ilg), &
             vgat%angsmoss (ilg), &
             vgat%ancmoss (ilg), &
             vgat%angmoss (ilg), &
             vgat%rmlcsmoss (ilg), &
             vgat%rmlgsmoss (ilg), &
             vgat%rmlcmoss (ilg), &
             vgat%rmlgmoss (ilg), &
             vgat%fsinacc_gat(ilg), &
             vgat%flutacc_gat(ilg), &
             vgat%flinacc_gat(ilg), &
             vgat%altotacc_gat(ilg), &
             vgat%netrad_gat(ilg), &
             vgat%preacc_gat(ilg), &
             vgat%sdepgat(ilg), &
             vgat%xdiffusgat(ilg), & ! the corresponding ROW is CLASS's XDIFFUS
             vgat%faregat(ilg), &    ! the ROT is FAREROT
             tracer%mossCMassgat(ilg), &
             tracer%mossLitrMassgat(ilg), &
             tracer%tracerCO2gat(ilg), &

             ! allocated with ilg, icc
             vgat%gleafmas (ilg,icc), &
             vgat%bleafmas (ilg,icc), &
             vgat%stemmass (ilg,icc), &
             vgat%rootmass (ilg,icc), &
             vgat%pstemmass (ilg,icc), &
             vgat%pgleafmass (ilg,icc), &
             vgat%fcancmx (ilg,icc), &
             vgat%ailcg (ilg,icc), &
             vgat%ailcgs (ilg,icc), &
             vgat%fcancs (ilg,icc), &
             vgat%fcanc (ilg,icc), &
             vgat%co2i1cg (ilg,icc), &
             vgat%co2i1cs (ilg,icc), &
             vgat%co2i2cg (ilg,icc), &
             vgat%co2i2cs (ilg,icc), &
             vgat%ancsveg (ilg,icc), &
             vgat%ancgveg (ilg,icc), &
             vgat%rmlcsveg (ilg,icc), &
             vgat%rmlcgveg (ilg,icc), &
             vgat%slai (ilg,icc), &
             vgat%ailcb (ilg,icc), &
             vgat%flhrloss (ilg,icc), &
             vgat%grwtheff (ilg,icc), &
             vgat%lystmmas (ilg,icc), &
             vgat%lyrotmas (ilg,icc), &
             vgat%tymaxlai (ilg,icc), &
             vgat%stmhrlos (ilg,icc), &
             vgat%vgbiomas_veg (ilg,icc), &
             vgat%emit_co2 (ilg,icc), &
             vgat%emit_co (ilg,icc), &
             vgat%emit_ch4 (ilg,icc), &
             vgat%emit_nmhc (ilg,icc), &
             vgat%emit_h2 (ilg,icc), &
             vgat%emit_nox (ilg,icc), &
             vgat%emit_n2o (ilg,icc), &
             vgat%emit_pm25 (ilg,icc), &
             vgat%emit_tpm (ilg,icc), &
             vgat%emit_tc (ilg,icc), &
             vgat%emit_oc (ilg,icc), &
             vgat%emit_bc (ilg,icc), &
             vgat%burnvegf (ilg,icc), &
             vgat%smfuncveg (ilg,icc), &
             vgat%bterm (ilg,icc), &
             vgat%mterm (ilg,icc), &
             vgat%bmasveg (ilg,icc), &
             vgat%veghght (ilg,icc), &
             vgat%rootdpth (ilg,icc), &
             vgat%tltrleaf (ilg,icc), &
             vgat%blfltrdt (ilg,icc), &
             vgat%glfltrdt (ilg,icc), &
             vgat%glcaemls(ilg,icc), &
             vgat%blcaemls(ilg,icc), &
             vgat%rtcaemls(ilg,icc), &
             vgat%stcaemls(ilg,icc), &
             vgat%ltrcemls(ilg,icc), &
             vgat%ntchlveg(ilg,icc), &
             vgat%ntchsveg(ilg,icc), &
             vgat%ntchrveg(ilg,icc), &
             vgat%tltrstem (ilg,icc), &
             vgat%tltrroot (ilg,icc), &
             vgat%leaflitr (ilg,icc), &
             vgat%roottemp (ilg,icc), &
             vgat%afrleaf (ilg,icc), &
             vgat%afrstem (ilg,icc), &
             vgat%afrroot (ilg,icc), &
             vgat%wtstatus (ilg,icc), &
             vgat%ltstatus (ilg,icc), &
             vgat%gppveg (ilg,icc), &
             vgat%nppveg (ilg,icc), &
             vgat%autoresveg (ilg,icc), &
             vgat%rmlvegacc (ilg,icc), &
             vgat%rmsveg (ilg,icc), &
             vgat%rmrveg (ilg,icc), &
             vgat%rgveg (ilg,icc), &
             vgat%litrfallveg (ilg,icc), &
             vgat%reprocost(ilg,icc), &
             vgat%anveg (ilg,icc), &
             vgat%rmlveg (ilg,icc), &
             vgat%geremort (ilg,icc), &
             vgat%intrmort (ilg,icc), &
             vgat%lambda (ilg,icc), &
             vgat%cc (ilg,icc), &
             vgat%mm (ilg,icc), &
             vgat%pftexist (ilg,icc), &
             vgat%rothrlos (ilg,icc), &
             vgat%pfcancmx (ilg,icc), &
             vgat%nfcancmx (ilg,icc), &
             vgat%lfstatus (ilg,icc), &
             vgat%pandays (ilg,icc), &
             vgat%todfrac(ilg,icc), &
             tracer%gLeafMassgat(ilg,icc), &
             tracer%bLeafMassgat(ilg,icc), &
             tracer%stemMassgat(ilg,icc), &
             tracer%rootMassgat(ilg,icc), &

             ! allocated with ilg, ican:
             vgat%zolnc (ilg,ican), &
             vgat%ailc (ilg,ican), &
             vgat%cmasvegc (ilg,ican), &
             vgat%alvsctm (ilg,ican), &
             vgat%paic (ilg,ican), &
             vgat%slaic (ilg,ican), &
             vgat%alirctm (ilg,ican), &

             ! allocated with ilg, iccp1:
             vgat%nepveg (ilg,iccp1), &
             vgat%nbpveg (ilg,iccp1), &
             vgat%hetroresveg (ilg,iccp1), &

             ! allocated with ilg, iccp2, ignd:
             ! COMBAK PERLAY
             vgat%litrmass (ilg,iccp2), &
             vgat%soilcmas (ilg,iccp2), &
             vgat%litresveg (ilg,iccp2), &
             vgat%soilcresveg (ilg,iccp2), &
             vgat%humiftrsveg (ilg,iccp2), &
             ! vgat%litrmass (ilg,iccp2,ignd), &
             ! vgat%soilcmas (ilg,iccp2,ignd), &
             ! vgat%litresveg (ilg,iccp2,ignd), &
             ! vgat%soilcresveg (ilg,iccp2,ignd), &
             ! vgat%humiftrsveg (ilg,iccp2,ignd), &
             ! COMBAK PERLAY
             tracer%litrMassgat(ilg,iccp2,ignd), &
             tracer%soilCMassgat(ilg,iccp2,ignd), &

             ! allocated with ilg, ignd:
             !          vgat%rgmgat(ilg,ignd), &
             vgat%sandgat(ilg,ignd), &
             vgat%claygat(ilg,ignd), &
             vgat%orgmgat(ilg,ignd), &
             vgat%socigat(ilg), &

             ! allocated with ilg, ican, ignd:
             vgat%rmatc (ilg,ican,ignd), &

             ! allocated with ilg, icc, ignd:
             vgat%rmatctem (ilg,icc,ignd), &

             ! allocated with ilg, icc, {some number}:
             vgat%colddays (ilg,2), &
             vgat%slopefrac (ilg,8), &
             vgat%tmonth (12,ilg), &

             ! allocated with ilg
             ctem_tile%fsnowacc_t (ilg), &
             ctem_tile%tcansacc_t (ilg), &
             ctem_tile%tcanoaccgat_t (ilg), &
             ctem_tile%taaccgat_t (ilg), &
             ctem_tile%uvaccgat_t (ilg), &
             ctem_tile%vvaccgat_t (ilg), &
             ctem_tile%anmossac_t (ilg), &
             ctem_tile%rmlmossac_t (ilg), &
             ctem_tile%gppmossac_t (ilg), &
             ctem_tile%tbaraccgat_t (ilg,ignd), &
             ctem_tile%thliqacc_t (ilg,ignd), &
             ctem_tile%thiceacc_t (ilg,ignd), &
             ctem_tile%ancgvgac_t (ilg,icc), &
             ctem_tile%rmlcgvga_t (ilg,icc), &

             ctem_mo%laimaxg_mo (nlat,nmos,icc), &
             ctem_mo%stemmass_mo (nlat,nmos,icc), &
             ctem_mo%leafmass_mo (nlat,nmos,icc), &
             ctem_mo%rootmass_mo (nlat,nmos,icc), &
             ctem_mo%litrfallveg_mo (nlat,nmos,icc), &
             ctem_mo%npp_mo (nlat,nmos,icc), &
             ctem_mo%gpp_mo (nlat,nmos,icc), &
             ctem_mo%vgbiomas_mo (nlat,nmos,icc), &
             ctem_mo%autores_mo (nlat,nmos,icc), &
             ctem_mo%soilres_mo (nlat,nmos,icc), &
             ctem_mo%humiftrsveg_mo (nlat,nmos,iccp2), &
             ctem_mo%totcmass_mo (nlat,nmos,iccp1), &
             ! COMBAK PERLAY
             ctem_mo%litrmass_mo (nlat,nmos,iccp2), &
             ctem_mo%soilcmas_mo (nlat,nmos,iccp2), &
             ctem_mo%litres_mo (nlat,nmos,iccp2), &
             ctem_mo%soilcres_mo (nlat,nmos,iccp2), &
             ! ctem_mo%litrmass_mo (nlat,nmos,iccp2,ignd), &
             ! ctem_mo%soilcmas_mo (nlat,nmos,iccp2,ignd), &
             ! ctem_mo%litres_mo (nlat,nmos,iccp2,ignd), &
             ! ctem_mo%soilcres_mo (nlat,nmos,iccp2,ignd), &
             ! COMBAK PERLAY
             ctem_mo%nep_mo (nlat,nmos,iccp1), &
             ctem_mo%hetrores_mo (nlat,nmos,iccp1), &
             ctem_mo%nbp_mo (nlat,nmos,iccp1), &
             ctem_mo%emit_co2_mo (nlat,nmos,icc), &
             ctem_mo%emit_co_mo (nlat,nmos,icc), &
             ctem_mo%emit_ch4_mo (nlat,nmos,icc), &
             ctem_mo%emit_nmhc_mo (nlat,nmos,icc), &
             ctem_mo%emit_h2_mo (nlat,nmos,icc), &
             ctem_mo%emit_nox_mo (nlat,nmos,icc), &
             ctem_mo%emit_n2o_mo (nlat,nmos,icc), &
             ctem_mo%emit_pm25_mo (nlat,nmos,icc), &
             ctem_mo%emit_tpm_mo (nlat,nmos,icc), &
             ctem_mo%emit_tc_mo (nlat,nmos,icc), &
             ctem_mo%emit_oc_mo (nlat,nmos,icc), &
             ctem_mo%emit_bc_mo (nlat,nmos,icc), &
             ctem_mo%burnfrac_mo (nlat,nmos,icc), &
             ctem_mo%bterm_mo (nlat,nmos,icc), &
             ctem_mo%mterm_mo (nlat,nmos,icc), &
             ctem_mo%smfuncveg_mo (nlat,nmos,icc), &

             ctem_grd_mo%laimaxg_mo_g (nlat), &
             ctem_grd_mo%stemmass_mo_g (nlat), &
             ctem_grd_mo%leafmass_mo_g (nlat), &
             ctem_grd_mo%rootmass_mo_g (nlat), &
             ! COMBAK PERLAY
             ctem_grd_mo%litrmass_mo_g (nlat), &
             ctem_grd_mo%soilcmas_mo_g (nlat), &
             ctem_grd_mo%litres_mo_g (nlat), &
             ctem_grd_mo%soilcres_mo_g (nlat), &
             ! ctem_grd_mo%litrmass_mo_g (nlat,ignd), &
             ! ctem_grd_mo%soilcmas_mo_g (nlat,ignd), &
             ! ctem_grd_mo%litres_mo_g (nlat,ignd), &
             ! ctem_grd_mo%soilcres_mo_g (nlat,ignd), &
             ! COMBAK PERLAY
             ctem_grd_mo%litrfall_mo_g (nlat), &
             ctem_grd_mo%humiftrs_mo_g (nlat), &
             ctem_grd_mo%npp_mo_g (nlat), &
             ctem_grd_mo%gpp_mo_g (nlat), &
             ctem_grd_mo%nep_mo_g (nlat), &
             ctem_grd_mo%nbp_mo_g (nlat), &
             ctem_grd_mo%hetrores_mo_g (nlat), &
             ctem_grd_mo%autores_mo_g (nlat), &
             ctem_grd_mo%soilres_mo_g (nlat), &
             ctem_grd_mo%vgbiomas_mo_g (nlat), &
             ctem_grd_mo%totcmass_mo_g (nlat), &
             ctem_grd_mo%emit_co2_mo_g (nlat), &
             ctem_grd_mo%emit_co_mo_g (nlat), &
             ctem_grd_mo%emit_ch4_mo_g (nlat), &
             ctem_grd_mo%emit_nmhc_mo_g (nlat), &
             ctem_grd_mo%emit_h2_mo_g (nlat), &
             ctem_grd_mo%emit_nox_mo_g (nlat), &
             ctem_grd_mo%emit_n2o_mo_g (nlat), &
             ctem_grd_mo%emit_pm25_mo_g (nlat), &
             ctem_grd_mo%emit_tpm_mo_g (nlat), &
             ctem_grd_mo%emit_tc_mo_g (nlat), &
             ctem_grd_mo%emit_oc_mo_g (nlat), &
             ctem_grd_mo%emit_bc_mo_g (nlat), &
             ctem_grd_mo%smfuncveg_mo_g (nlat), &
             ctem_grd_mo%luc_emc_mo_g (nlat), &
             ctem_grd_mo%lucltrin_mo_g (nlat), &
             ctem_grd_mo%lucsocin_mo_g (nlat), &
             ctem_grd_mo%burnfrac_mo_g (nlat), &
             ctem_grd_mo%bterm_mo_g (nlat), &
             ctem_grd_mo%lterm_mo_g(nlat), &
             ctem_grd_mo%mterm_mo_g (nlat), &
             ctem_grd_mo%ch4WetSpec_mo_g (nlat), &
             ctem_grd_mo%wetfdyn_mo_g (nlat), &
             ctem_grd_mo%wetfpres_mo_g (nlat), &
             ctem_grd_mo%ch4WetDyn_mo_g (nlat), &
             ctem_grd_mo%ch4soills_mo_g (nlat), &
             ctem_grd_mo%cProduct_mo_g (nlat), &
             ctem_grd_mo%fProductDecomp_mo_g (nlat), &

             ctem_tile_mo%laimaxg_mo_t (nlat,nmos), &
             ctem_tile_mo%stemmass_mo_t (nlat,nmos), &
             ctem_tile_mo%leafmass_mo_t (nlat,nmos), &
             ctem_tile_mo%rootmass_mo_t (nlat,nmos), &
             ctem_tile_mo%litrfall_mo_t (nlat,nmos), &
             ctem_tile_mo%humiftrs_mo_t (nlat,nmos), &
             ctem_tile_mo%npp_mo_t (nlat,nmos), &
             ctem_tile_mo%gpp_mo_t (nlat,nmos), &
             ctem_tile_mo%vgbiomas_mo_t (nlat,nmos), &
             ctem_tile_mo%autores_mo_t (nlat,nmos), &
             ctem_tile_mo%soilres_mo_t (nlat,nmos), &
             ctem_tile_mo%totcmass_mo_t (nlat,nmos), &
             ! COMBAK PERLAY
             ctem_tile_mo%litrmass_mo_t (nlat,nmos), &
             ctem_tile_mo%soilcmas_mo_t (nlat,nmos), &
             ctem_tile_mo%litres_mo_t (nlat,nmos), &
             ctem_tile_mo%soilcres_mo_t (nlat,nmos), &
             ! ctem_tile_mo%litrmass_mo_t (nlat,nmos,ignd), &
             ! ctem_tile_mo%soilcmas_mo_t (nlat,nmos,ignd), &
             ! ctem_tile_mo%litres_mo_t (nlat,nmos,ignd), &
             ! ctem_tile_mo%soilcres_mo_t (nlat,nmos,ignd), &
             ! COMBAK PERLAY
             ctem_tile_mo%nep_mo_t (nlat,nmos), &
             ctem_tile_mo%hetrores_mo_t (nlat,nmos), &
             ctem_tile_mo%nbp_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_co2_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_co_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_ch4_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_nmhc_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_h2_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_nox_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_n2o_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_pm25_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_tpm_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_tc_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_oc_mo_t (nlat,nmos), &
             ctem_tile_mo%emit_bc_mo_t (nlat,nmos), &
             ctem_tile_mo%burnfrac_mo_t (nlat,nmos), &
             ctem_tile_mo%smfuncveg_mo_t (nlat,nmos), &
             ctem_tile_mo%bterm_mo_t (nlat,nmos), &
             ctem_tile_mo%luc_emc_mo_t (nlat,nmos), &
             ctem_tile_mo%lterm_mo_t (nlat,nmos), &
             ctem_tile_mo%lucsocin_mo_t (nlat,nmos), &
             ctem_tile_mo%mterm_mo_t (nlat,nmos), &
             ctem_tile_mo%lucltrin_mo_t (nlat,nmos), &
             ctem_tile_mo%ch4WetSpec_mo_t (nlat,nmos), &
             ctem_tile_mo%wetfdyn_mo_t (nlat,nmos), &
             ctem_tile_mo%wetfpres_mo_t (nlat,nmos), &
             ctem_tile_mo%ch4WetDyn_mo_t (nlat,nmos), &
             ctem_tile_mo%ch4soills_mo_t (nlat,nmos), &
             ctem_tile_mo%wind_mo_t (nlat,nmos), &
             ctem_tile_mo%fProductDecomp_mo_t (nlat,nmos), &

             ctem_yr%laimaxg_yr (nlat,nmos,icc), &
             ctem_yr%stemmass_yr (nlat,nmos,icc), &
             ctem_yr%rootmass_yr (nlat,nmos,icc), &
             ctem_yr%npp_yr (nlat,nmos,icc), &
             ctem_yr%gpp_yr (nlat,nmos,icc), &
             ctem_yr%vgbiomas_yr (nlat,nmos,icc), &
             ctem_yr%autores_yr (nlat,nmos,icc), &
             ctem_yr%totcmass_yr (nlat,nmos,iccp1), &
             ! COMBAK PERLAY
             ctem_yr%litrmass_yr (nlat,nmos,iccp2), &
             ctem_yr%soilcmas_yr (nlat,nmos,iccp2), &
             ctem_yr%litres_yr (nlat,nmos,iccp2), &
             ctem_yr%soilcres_yr (nlat,nmos,iccp2), &
             ! ctem_yr%litrmass_yr (nlat,nmos,iccp2,ignd), &
             ! ctem_yr%soilcmas_yr (nlat,nmos,iccp2,ignd), &
             ! ctem_yr%litres_yr (nlat,nmos,iccp2,ignd), &
             ! ctem_yr%soilcres_yr (nlat,nmos,iccp2,ignd), &
             ! COMBAK PERLAY
             ctem_yr%nep_yr (nlat,nmos,iccp1), &
             ctem_yr%hetrores_yr (nlat,nmos,iccp1), &
             ctem_yr%nbp_yr (nlat,nmos,iccp1), &
             ctem_yr%emit_co2_yr (nlat,nmos,icc), &
             ctem_yr%emit_co_yr (nlat,nmos,icc), &
             ctem_yr%emit_ch4_yr (nlat,nmos,icc), &
             ctem_yr%emit_nmhc_yr (nlat,nmos,icc), &
             ctem_yr%emit_h2_yr (nlat,nmos,icc), &
             ctem_yr%emit_nox_yr (nlat,nmos,icc), &
             ctem_yr%emit_n2o_yr (nlat,nmos,icc), &
             ctem_yr%emit_pm25_yr (nlat,nmos,icc), &
             ctem_yr%emit_tpm_yr (nlat,nmos,icc), &
             ctem_yr%emit_tc_yr (nlat,nmos,icc), &
             ctem_yr%emit_oc_yr (nlat,nmos,icc), &
             ctem_yr%emit_bc_yr (nlat,nmos,icc), &
             ctem_yr%bterm_yr (nlat,nmos,icc), &
             ctem_yr%mterm_yr (nlat,nmos,icc), &
             ctem_yr%burnfrac_yr (nlat,nmos,icc), &
             ctem_yr%smfuncveg_yr (nlat,nmos,icc), &
             ctem_yr%veghght_yr (nlat,nmos,icc), &

             ctem_grd_yr%laimaxg_yr_g (nlat), &
             ctem_grd_yr%stemmass_yr_g (nlat), &
             ctem_grd_yr%rootmass_yr_g (nlat), &
             ! COMBAK PERLAY
             ctem_grd_yr%litrmass_yr_g (nlat), &
             ctem_grd_yr%soilcmas_yr_g (nlat), &
             ctem_grd_yr%litres_yr_g (nlat), &
             ctem_grd_yr%soilcres_yr_g (nlat), &
             ! ctem_grd_yr%litrmass_yr_g (nlat,ignd), &
             ! ctem_grd_yr%soilcmas_yr_g (nlat,ignd), &
             ! ctem_grd_yr%litres_yr_g (nlat,ignd), &
             ! ctem_grd_yr%soilcres_yr_g (nlat,ignd), &
             ! COMBAK PERLAY
             ctem_grd_yr%npp_yr_g (nlat), &
             ctem_grd_yr%gpp_yr_g (nlat), &
             ctem_grd_yr%nep_yr_g (nlat), &
             ctem_grd_yr%nbp_yr_g (nlat), &
             ctem_grd_yr%hetrores_yr_g (nlat), &
             ctem_grd_yr%autores_yr_g (nlat), &
             ctem_grd_yr%vgbiomas_yr_g (nlat), &
             ctem_grd_yr%totcmass_yr_g (nlat), &
             ctem_grd_yr%emit_co2_yr_g (nlat), &
             ctem_grd_yr%emit_co_yr_g (nlat), &
             ctem_grd_yr%emit_ch4_yr_g (nlat), &
             ctem_grd_yr%emit_nmhc_yr_g (nlat), &
             ctem_grd_yr%emit_h2_yr_g (nlat), &
             ctem_grd_yr%emit_nox_yr_g (nlat), &
             ctem_grd_yr%emit_n2o_yr_g (nlat), &
             ctem_grd_yr%emit_pm25_yr_g (nlat), &
             ctem_grd_yr%emit_tpm_yr_g (nlat), &
             ctem_grd_yr%emit_tc_yr_g (nlat), &
             ctem_grd_yr%emit_oc_yr_g (nlat), &
             ctem_grd_yr%emit_bc_yr_g (nlat), &
             ctem_grd_yr%smfuncveg_yr_g (nlat), &
             ctem_grd_yr%luc_emc_yr_g (nlat), &
             ctem_grd_yr%lucltrin_yr_g (nlat), &
             ctem_grd_yr%lucsocin_yr_g (nlat), &
             ctem_grd_yr%burnfrac_yr_g (nlat), &
             ctem_grd_yr%bterm_yr_g (nlat), &
             ctem_grd_yr%lterm_yr_g (nlat), &
             ctem_grd_yr%mterm_yr_g (nlat), &
             ctem_grd_yr%ch4WetSpec_yr_g (nlat), &
             ctem_grd_yr%wetfdyn_yr_g (nlat), &
             ctem_grd_yr%ch4WetDyn_yr_g (nlat), &
             ctem_grd_yr%ch4soills_yr_g (nlat), &
             ctem_grd_yr%veghght_yr_g (nlat), &
             ctem_grd_yr%peatdep_yr_g (nlat), &
             ctem_grd_yr%cProduct_yr_g (nlat), &
             ctem_grd_yr%fProductDecomp_yr_g (nlat), &

             ctem_tile_yr%laimaxg_yr_t (nlat,nmos), &
             ctem_tile_yr%stemmass_yr_t (nlat,nmos), &
             ctem_tile_yr%rootmass_yr_t (nlat,nmos), &
             ctem_tile_yr%npp_yr_t (nlat,nmos), &
             ctem_tile_yr%gpp_yr_t (nlat,nmos), &
             ctem_tile_yr%vgbiomas_yr_t (nlat,nmos), &
             ctem_tile_yr%autores_yr_t (nlat,nmos), &
             ctem_tile_yr%totcmass_yr_t (nlat,nmos), &
             ! COMBAK PERLAY
             ctem_tile_yr%litrmass_yr_t (nlat,nmos), &
             ctem_tile_yr%soilcmas_yr_t (nlat,nmos), &
             ctem_tile_yr%litres_yr_t (nlat,nmos), &
             ctem_tile_yr%soilcres_yr_t (nlat,nmos), &
             ! ctem_tile_yr%litrmass_yr_t (nlat,nmos,ignd), &
             ! ctem_tile_yr%soilcmas_yr_t (nlat,nmos,ignd), &
             ! ctem_tile_yr%litres_yr_t (nlat,nmos,ignd), &
             ! ctem_tile_yr%soilcres_yr_t (nlat,nmos,ignd), &
             ! COMBAK PERLAY
             ctem_tile_yr%nep_yr_t (nlat,nmos), &
             ctem_tile_yr%hetrores_yr_t (nlat,nmos), &
             ctem_tile_yr%nbp_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_co2_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_co_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_ch4_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_nmhc_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_h2_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_nox_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_n2o_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_pm25_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_tpm_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_tc_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_oc_yr_t (nlat,nmos), &
             ctem_tile_yr%emit_bc_yr_t (nlat,nmos), &
             ctem_tile_yr%burnfrac_yr_t (nlat,nmos), &
             ctem_tile_yr%smfuncveg_yr_t (nlat,nmos), &
             ctem_tile_yr%bterm_yr_t (nlat,nmos), &
             ctem_tile_yr%luc_emc_yr_t (nlat,nmos), &
             ctem_tile_yr%lterm_yr_t (nlat,nmos), &
             ctem_tile_yr%lucsocin_yr_t (nlat,nmos), &
             ctem_tile_yr%mterm_yr_t (nlat,nmos), &
             ctem_tile_yr%lucltrin_yr_t (nlat,nmos), &
             ctem_tile_yr%ch4WetSpec_yr_t (nlat,nmos), &
             ctem_tile_yr%wetfdyn_yr_t (nlat,nmos), &
             ctem_tile_yr%ch4WetDyn_yr_t (nlat,nmos), &
             ctem_tile_yr%ch4soills_yr_t (nlat,nmos), &
             ctem_tile_yr%veghght_yr_t (nlat,nmos), &
             ctem_tile_yr%peatdep_yr_t (nlat,nmos), &
             ctem_tile_yr%fProductDecomp_yr_t (nlat,nmos))
		
		
  end subroutine allocCtemVars
  !! @}

  ! -----------------------------------------------------

  !> \ingroup ctemstatevars_initRowVars
  !! @{
  !> Initializes 'row' variables
  subroutine initRowVars

    implicit none

    ! nlat, nmos
    vrot%co2conc    = 0.0
    vrot%npp        = 0.0
    vrot%nep        = 0.0
    vrot%hetrores   = 0.0
    vrot%autores    = 0.0
    vrot%soilcresp  = 0.0
    vrot%rm         = 0.0
    vrot%rg         = 0.0
    vrot%nbp        = 0.0
    vrot%litres     = 0.0
    vrot%socres     = 0.0
    vrot%gpp        = 0.0
    vrot%dstcemls   = 0.0
    vrot%dstcemls3  = 0.0
    vrot%litrfall   = 0.0
    vrot%humiftrs   = 0.0
    vrot%canres     = 0.0
    vrot%rml        = 0.0
    vrot%rms        = 0.0
    vrot%rmr        = 0.0
    vrot%lucemcom   = 0.0
    vrot%lucltrin   = 0.0
    vrot%lucsocin   = 0.0
    vrot%burnfrac   = 0.0
    vrot%lterm      = 0.0
    vrot%cfluxcg    = 0.0
    vrot%cfluxcs    = 0.0
    vrot%ch4WetSpec = 0.0
    vrot%wetfdyn    = 0.0
    vrot%ch4WetDyn  = 0.0
    vrot%ch4_soills = 0.0
    vrot%nppmoss    = 0.0
    vrot%rmlmoss    = 0.0
    vrot%gppmoss    = 0.0
    vrot%anmoss     = 0.0
    vrot%armoss     = 0.0
    vrot%peatdep    = 0.0
    vrot%pdd        = 0.0

    vrot%ZOLNC      = 0.0
    vrot%AILC       = 0.0
    vrot%CMASVEGC   = 0.0
    vrot%ALVSCTM    = 0.0
    vrot%ALIRCTM    = 0.0
    vrot%CSUM       = 0.0
    vrot%PAIC       = 0.0
    vrot%SLAIC      = 0.0
    vrot%RMATC      = 0.0

    vrot%smfuncveg  = 0.0
    vrot%gleafmas   = 0.
    vrot%bleafmas   = 0.
    vrot%stemmass   = 0.
    vrot%rootmass   = 0.
    vrot%pstemmass  = 0.
    vrot%pgleafmass = 0.
    vrot%litrfallveg = 0.
    vrot%bterm       = 0.0
    vrot%mterm       = 0.0
    vrot%ailcg       = 0.0
    vrot%ailcgs      = 0.0
    vrot%fcancs      = 0.0
    vrot%fcanc       = 0.0
    vrot%fcancmx     = 0.0
    vrot%co2i1cg     = 0.0
    vrot%co2i1cs     = 0.0
    vrot%co2i2cg     = 0.0
    vrot%co2i2cs     = 0.0
    vrot%ancsveg     = 0.0
    vrot%ancgveg     = 0.0
    vrot%rmlcsveg    = 0.0
    vrot%rmlcgveg    = 0.0
    vrot%stemmass    = 0.0
    vrot%rootmass    = 0.0
    vrot%ailcb       = 0.0
    vrot%grwtheff    = 0.0
    vrot%bmasveg     = 0.0
    vrot%tltrleaf    = 0.0
    vrot%tltrstem    = 0.0
    vrot%tltrroot    = 0.0
    vrot%leaflitr    = 0.0
    vrot%roottemp    = 0.0
    vrot%afrleaf     = 0.0
    vrot%afrstem     = 0.0
    vrot%afrroot     = 0.0
    vrot%wtstatus    = 0.0
    vrot%ltstatus    = 0.0
    vrot%pfcancmx    = 0.0
    vrot%nfcancmx    = 0.0
    vrot%nppveg      = 0.0
    vrot%veghght     = 0.0
    vrot%rootdpth    = 0.0
    vrot%gleafmas    = 0.0
    vrot%bleafmas    = 0.0
    vrot%anveg       = 0.0
    vrot%rmlveg      = 0.0
    vrot%rmlvegacc   = 0.0
    vrot%rmsveg      = 0.0
    vrot%rmrveg      = 0.0
    vrot%rgveg       = 0.0
    vrot%vgbiomas_veg = 0.0
    vrot%gppveg      = 0.0
    vrot%autoresveg  = 0.0
    vrot%emit_co2    = 0.0
    vrot%emit_co     = 0.0
    vrot%emit_ch4    = 0.0
    vrot%emit_nmhc   = 0.0
    vrot%emit_h2     = 0.0
    vrot%emit_nox    = 0.0
    vrot%emit_n2o    = 0.0
    vrot%emit_pm25   = 0.0
    vrot%emit_tpm    = 0.0
    vrot%emit_tc     = 0.0
    vrot%emit_oc     = 0.0
    vrot%emit_bc     = 0.0
    vrot%burnvegf    = 0.0

    vrot%rmatctem    = 0.0
    vrot%hetroresveg = 0.0
    vrot%nepveg      = 0.0
    vrot%nbpveg      = 0.0
    vrot%litrmass    = 0.0
    vrot%soilcmas    = 0.0
    vrot%litresveg   = 0.0
    vrot%soilcresveg = 0.0
    vrot%humiftrsveg = 0.0

  end subroutine initRowVars
  !! @}

  !==================================================

  !> \ingroup ctemstatevars_resetMonthEnd
  !! @{
  !> Resets monthly variables at month end in preparation for next month
  subroutine resetMonthEnd (nltest, nmtest)

    use classicParams,    only : iccp2, icc, iccp1, ignd

    implicit none

    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest

    integer :: i, m, j

    ! These are assigned to mid-month, but are not accumulated so can be
    ! zeroed out at the same time as the other month-end vars.
    do i = 1,nltest
      ctem_grd_mo%stemmass_mo_g(i) = 0.0
      ctem_grd_mo%leafmass_mo_g(i) = 0.0
      ctem_grd_mo%rootmass_mo_g(i) = 0.0
      ! COMBAK PERLAY
      ctem_grd_mo%litrmass_mo_g(i) = 0.0
      ctem_grd_mo%soilcmas_mo_g(i) = 0.0
      ! ctem_grd_mo%litrmass_mo_g(i,1:ignd)=0.0
      ! ctem_grd_mo%soilcmas_mo_g(i,1:ignd)=0.0
      ! COMBAK PERLAY
      ctem_grd_mo%vgbiomas_mo_g(i) = 0.0
      ctem_grd_mo%totcmass_mo_g(i) = 0.0
      do m = 1,nmtest
        ctem_tile_mo%stemmass_mo_t(i,m) = 0.0
        ctem_tile_mo%leafmass_mo_t(i,m) = 0.0
        ctem_tile_mo%rootmass_mo_t(i,m) = 0.0
        ! COMBAK PERLAY
        ctem_tile_mo%litrmass_mo_t(i,m) = 0.0
        ctem_tile_mo%soilcmas_mo_t(i,m) = 0.0
        ! ctem_tile_mo%litrmass_mo_t(i,m,1:ignd)=0.0
        ! ctem_tile_mo%soilcmas_mo_t(i,m,1:ignd)=0.0
        ! COMBAK PERLAY
        ctem_tile_mo%vgbiomas_mo_t(i,m) = 0.0
        ctem_tile_mo%totcmass_mo_t(i,m) = 0.0
        do j = 1,icc
          ctem_mo%stemmass_mo(i,m,j) = 0.0
          ctem_mo%leafmass_mo(i,m,j) = 0.0
          ctem_mo%rootmass_mo(i,m,j) = 0.0
          ! COMBAK PERLAY
          ctem_mo%litrmass_mo(i,m,j) = 0.0
          ctem_mo%soilcmas_mo(i,m,j) = 0.0
          ! ctem_mo%litrmass_mo(i,m,j,1:ignd)=0.0
          ! ctem_mo%soilcmas_mo(i,m,j,1:ignd)=0.0
          ! COMBAK PERLAY
          ctem_mo%vgbiomas_mo(i,m,j) = 0.0
          ctem_mo%totcmass_mo(i,m,j) = 0.0
        end do
        ctem_mo%totcmass_mo(i,m,iccp1) = 0.0
        ! COMBAK PERLAY
        ctem_mo%litrmass_mo(i,m,iccp1) = 0.0
        ctem_mo%soilcmas_mo(i,m,iccp1) = 0.0
        ctem_mo%litrmass_mo(i,m,iccp2) = 0.0
        ctem_mo%soilcmas_mo(i,m,iccp2) = 0.0
        ! ctem_mo%litrmass_mo(i,m,iccp1,1:ignd)=0.0
        ! ctem_mo%soilcmas_mo(i,m,iccp1,1:ignd)=0.0
        ! ctem_mo%litrmass_mo(i,m,iccp2,1:ignd)=0.0
        ! ctem_mo%soilcmas_mo(i,m,iccp2,1:ignd)=0.0
        ! COMBAK PERLAY
      end do
    end do

    ! Now zero out the month end vars.
    do i = 1,nltest
      ! Grid avg
      ctem_grd_mo%laimaxg_mo_g(i) = 0.0
      ctem_grd_mo%npp_mo_g(i) = 0.0
      ctem_grd_mo%gpp_mo_g(i) = 0.0
      ctem_grd_mo%nep_mo_g(i) = 0.0
      ctem_grd_mo%nbp_mo_g(i) = 0.0
      ctem_grd_mo%hetrores_mo_g(i) = 0.0
      ctem_grd_mo%autores_mo_g(i) = 0.0
      ctem_grd_mo%soilres_mo_g(i) = 0.0
      ! COMBAK PERLAY
      ctem_grd_mo%litres_mo_g(i) = 0.0
      ctem_grd_mo%soilcres_mo_g(i) = 0.0
      ! ctem_grd_mo%litres_mo_g(i,1:ignd)=0.0
      ! ctem_grd_mo%soilcres_mo_g(i,1:ignd)=0.0
      ! COMBAK PERLAY
      ctem_grd_mo%litrfall_mo_g(i) = 0.0
      ctem_grd_mo%humiftrs_mo_g(i) = 0.0
      ctem_grd_mo%emit_co2_mo_g(i) = 0.0
      ctem_grd_mo%emit_co_mo_g(i) = 0.0
      ctem_grd_mo%emit_ch4_mo_g(i) = 0.0
      ctem_grd_mo%emit_nmhc_mo_g(i) = 0.0
      ctem_grd_mo%emit_h2_mo_g(i) = 0.0
      ctem_grd_mo%emit_nox_mo_g(i) = 0.0
      ctem_grd_mo%emit_n2o_mo_g(i) = 0.0
      ctem_grd_mo%emit_pm25_mo_g(i) = 0.0
      ctem_grd_mo%emit_tpm_mo_g(i) = 0.0
      ctem_grd_mo%emit_tc_mo_g(i) = 0.0
      ctem_grd_mo%emit_oc_mo_g(i) = 0.0
      ctem_grd_mo%emit_bc_mo_g(i) = 0.0
      ctem_grd_mo%smfuncveg_mo_g(i) = 0.0
      ctem_grd_mo%luc_emc_mo_g(i) = 0.0
      ctem_grd_mo%lucsocin_mo_g(i) = 0.0
      ctem_grd_mo%lucltrin_mo_g(i) = 0.0
      ctem_grd_mo%burnfrac_mo_g(i) = 0.0
      ctem_grd_mo%bterm_mo_g(i)    = 0.0
      ctem_grd_mo%lterm_mo_g(i)    = 0.0
      ctem_grd_mo%mterm_mo_g(i)    = 0.0
      ctem_grd_mo%ch4WetSpec_mo_g(i)  = 0.0
      ctem_grd_mo%wetfdyn_mo_g(i)  = 0.0
      ctem_grd_mo%wetfpres_mo_g(i)  = 0.0
      ctem_grd_mo%ch4WetDyn_mo_g(i)  = 0.0
      ctem_grd_mo%ch4soills_mo_g(i)  = 0.0
      ctem_grd_mo%cProduct_mo_g(i)  = 0.0
      ctem_grd_mo%fProductDecomp_mo_g(i)  = 0.0

      do m = 1,nmtest
        ! Tile avg
        ctem_tile_mo%laimaxg_mo_t(i,m) = 0.0
        ctem_tile_mo%npp_mo_t(i,m) = 0.0
        ctem_tile_mo%gpp_mo_t(i,m) = 0.0
        ctem_tile_mo%nep_mo_t(i,m) = 0.0
        ctem_tile_mo%nbp_mo_t(i,m) = 0.0
        ctem_tile_mo%hetrores_mo_t(i,m) = 0.0
        ctem_tile_mo%autores_mo_t(i,m) = 0.0
        ctem_tile_mo%soilres_mo_t(i,m) = 0.0
        ! COMBAK PERLAY
        ctem_tile_mo%litres_mo_t(i,m) = 0.0
        ctem_tile_mo%soilcres_mo_t(i,m) = 0.0
        ! ctem_tile_mo%litres_mo_t(i,m,1:ignd)=0.0
        ! ctem_tile_mo%soilcres_mo_t(i,m,1:ignd)=0.0
        ! COMBAK PERLAY
        ctem_tile_mo%litrfall_mo_t(i,m) = 0.0
        ctem_tile_mo%humiftrs_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_co2_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_co_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_ch4_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_nmhc_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_h2_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_nox_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_n2o_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_pm25_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_tpm_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_tc_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_oc_mo_t(i,m) = 0.0
        ctem_tile_mo%emit_bc_mo_t(i,m) = 0.0
        ctem_tile_mo%smfuncveg_mo_t(i,m) = 0.0
        ctem_tile_mo%luc_emc_mo_t(i,m) = 0.0
        ctem_tile_mo%lucsocin_mo_t(i,m) = 0.0
        ctem_tile_mo%lucltrin_mo_t(i,m) = 0.0
        ctem_tile_mo%burnfrac_mo_t(i,m) = 0.0
        ctem_tile_mo%bterm_mo_t(i,m)    = 0.0
        ctem_tile_mo%lterm_mo_t(i,m)    = 0.0
        ctem_tile_mo%mterm_mo_t(i,m)    = 0.0
        ctem_tile_mo%ch4WetSpec_mo_t(i,m)  = 0.0
        ctem_tile_mo%wetfdyn_mo_t(i,m)  = 0.0
        ctem_tile_mo%wetfpres_mo_t(i,m)  = 0.0
        ctem_tile_mo%ch4WetDyn_mo_t(i,m)  = 0.0
        ctem_tile_mo%ch4soills_mo_t(i,m)  = 0.0
        ctem_tile_mo%wind_mo_t(i,m) = 0.0
        ctem_tile_mo%fProductDecomp_mo_t(i,m) = 0.0

        do j = 1,icc
          ! per pft
          ctem_mo%laimaxg_mo(i,m,j) = 0.0
          ctem_mo%npp_mo(i,m,j) = 0.0
          ctem_mo%gpp_mo(i,m,j) = 0.0
          ctem_mo%nep_mo(i,m,j) = 0.0
          ctem_mo%nbp_mo(i,m,j) = 0.0
          ctem_mo%hetrores_mo(i,m,j) = 0.0
          ctem_mo%autores_mo(i,m,j) = 0.0
          ctem_mo%soilres_mo(i,m,j) = 0.0
          ! COMBAK PERLAY
          ctem_mo%litres_mo(i,m,j) = 0.0
          ctem_mo%soilcres_mo(i,m,j) = 0.0
          ! ctem_mo%litres_mo(i,m,j,1:ignd)=0.0
          ! ctem_mo%soilcres_mo(i,m,j,1:ignd)=0.0
          ! COMBAK PERLAY
          ctem_mo%litrfallveg_mo(i,m,j) = 0.0
          ctem_mo%humiftrsveg_mo(i,m,j) = 0.0
          ctem_mo%emit_co2_mo(i,m,j) = 0.0
          ctem_mo%emit_co_mo(i,m,j) = 0.0
          ctem_mo%emit_ch4_mo(i,m,j) = 0.0
          ctem_mo%emit_nmhc_mo(i,m,j) = 0.0
          ctem_mo%emit_h2_mo(i,m,j) = 0.0
          ctem_mo%emit_nox_mo(i,m,j) = 0.0
          ctem_mo%emit_n2o_mo(i,m,j) = 0.0
          ctem_mo%emit_pm25_mo(i,m,j) = 0.0
          ctem_mo%emit_tpm_mo(i,m,j) = 0.0
          ctem_mo%emit_tc_mo(i,m,j) = 0.0
          ctem_mo%emit_oc_mo(i,m,j) = 0.0
          ctem_mo%emit_bc_mo(i,m,j) = 0.0
          ctem_mo%burnfrac_mo(i,m,j) = 0.0
          ctem_mo%bterm_mo(i,m,j) = 0.0
          ctem_mo%mterm_mo(i,m,j) = 0.0
          ctem_mo%smfuncveg_mo(i,m,j) = 0.0
        end do

        ctem_mo%nep_mo(i,m,iccp1) = 0.0
        ctem_mo%nbp_mo(i,m,iccp1) = 0.0
        ctem_mo%hetrores_mo(i,m,iccp1) = 0.0
        ctem_mo%humiftrsveg_mo(i,m,iccp1) = 0.0
        ctem_mo%humiftrsveg_mo(i,m,iccp2) = 0.0

        ! COMBAK PERLAY
        ctem_mo%litres_mo(i,m,iccp1) = 0.0
        ctem_mo%soilcres_mo(i,m,iccp1) = 0.0

        ctem_mo%litres_mo(i,m,iccp2) = 0.0
        ctem_mo%soilcres_mo(i,m,iccp2) = 0.0
        ! ctem_mo%litres_mo(i,m,iccp1,1:ignd)=0.0
        ! ctem_mo%soilcres_mo(i,m,iccp1,1:ignd)=0.0
        ! ctem_mo%humiftrsveg_mo(i,m,iccp1)=0.0
        !
        ! ctem_mo%litres_mo(i,m,iccp2,1:ignd)=0.0
        ! ctem_mo%soilcres_mo(i,m,iccp2,1:ignd)=0.0
        ! ctem_mo%humiftrsveg_mo(i,m,iccp2)=0.0
        ! COMBAK PERLAY

      end do ! nmtest
    end do ! nltest

  end subroutine resetMonthEnd
  !! @}
  !==================================================

  !> \ingroup ctemstatevars_resetYearEnd
  !! @{
  !> Resets annual variables in preparation for next year
  subroutine resetYearEnd (nltest, nmtest)

    use classicParams,   only : iccp2, icc, iccp1, ignd

    implicit none

    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest

    integer :: i, m, j

    do i = 1,nltest
      ! Grid avg
      ctem_grd_yr%laimaxg_yr_g(i) = 0.0
      ctem_grd_yr%stemmass_yr_g(i) = 0.0
      ctem_grd_yr%rootmass_yr_g(i) = 0.0
      ! COMBAK PERLAY
      ctem_grd_yr%litrmass_yr_g(i) = 0.0
      ctem_grd_yr%soilcmas_yr_g(i) = 0.0
      ctem_grd_yr%litres_yr_g(i) = 0.0
      ctem_grd_yr%soilcres_yr_g(i) = 0.0
      ! ctem_grd_yr%litrmass_yr_g(i,1:ignd)=0.0
      ! ctem_grd_yr%soilcmas_yr_g(i,1:ignd)=0.0
      ! ctem_grd_yr%litres_yr_g(i,1:ignd)=0.0
      ! ctem_grd_yr%soilcres_yr_g(i,1:ignd)=0.0
      ! COMBAK PERLAY
      ctem_grd_yr%vgbiomas_yr_g(i) = 0.0
      ctem_grd_yr%totcmass_yr_g(i) = 0.0
      ctem_grd_yr%veghght_yr_g(i) = 0.0
      ctem_grd_yr%npp_yr_g(i) = 0.0
      ctem_grd_yr%gpp_yr_g(i) = 0.0
      ctem_grd_yr%nep_yr_g(i) = 0.0
      ctem_grd_yr%nbp_yr_g(i) = 0.0
      ctem_grd_yr%hetrores_yr_g(i) = 0.0
      ctem_grd_yr%autores_yr_g(i) = 0.0
      ctem_grd_yr%emit_co2_yr_g(i) = 0.0
      ctem_grd_yr%emit_co_yr_g(i) = 0.0
      ctem_grd_yr%emit_ch4_yr_g(i) = 0.0
      ctem_grd_yr%emit_nmhc_yr_g(i) = 0.0
      ctem_grd_yr%emit_h2_yr_g(i) = 0.0
      ctem_grd_yr%emit_nox_yr_g(i) = 0.0
      ctem_grd_yr%emit_n2o_yr_g(i) = 0.0
      ctem_grd_yr%emit_pm25_yr_g(i) = 0.0
      ctem_grd_yr%emit_tpm_yr_g(i) = 0.0
      ctem_grd_yr%emit_tc_yr_g(i) = 0.0
      ctem_grd_yr%emit_oc_yr_g(i) = 0.0
      ctem_grd_yr%emit_bc_yr_g(i) = 0.0
      ctem_grd_yr%smfuncveg_yr_g(i) = 0.0
      ctem_grd_yr%luc_emc_yr_g(i) = 0.0
      ctem_grd_yr%lucsocin_yr_g(i) = 0.0
      ctem_grd_yr%lucltrin_yr_g(i) = 0.0
      ctem_grd_yr%burnfrac_yr_g(i) = 0.0
      ctem_grd_yr%bterm_yr_g(i) = 0.0
      ctem_grd_yr%lterm_yr_g(i) = 0.0
      ctem_grd_yr%mterm_yr_g(i) = 0.0
      ctem_grd_yr%ch4WetSpec_yr_g(i)  = 0.0
      ctem_grd_yr%wetfdyn_yr_g(i)  = 0.0
      ctem_grd_yr%ch4WetDyn_yr_g(i)  = 0.0
      ctem_grd_yr%ch4soills_yr_g(i)  = 0.0
      ctem_grd_yr%peatdep_yr_g(i)  = 0.0
      ctem_grd_yr%cProduct_yr_g(i)  = 0.0
      ctem_grd_yr%fProductDecomp_yr_g(i)  = 0.0

      do m = 1,nmtest
        ! Tile avg
        ctem_tile_yr%laimaxg_yr_t(i,m) = 0.0
        ctem_tile_yr%stemmass_yr_t(i,m) = 0.0
        ctem_tile_yr%rootmass_yr_t(i,m) = 0.0
        ! COMBAK PERLAY
        ctem_tile_yr%litrmass_yr_t(i,m) = 0.0
        ctem_tile_yr%soilcmas_yr_t(i,m) = 0.0
        ctem_tile_yr%litres_yr_t(i,m) = 0.0
        ctem_tile_yr%soilcres_yr_t(i,m) = 0.0
        ! ctem_tile_yr%litrmass_yr_t(i,m,1:ignd)=0.0
        ! ctem_tile_yr%soilcmas_yr_t(i,m,1:ignd)=0.0
        ! ctem_tile_yr%litres_yr_t(i,m,1:ignd)=0.0
        ! ctem_tile_yr%soilcres_yr_t(i,m,1:ignd)=0.0
        ! COMBAK PERLAY
        ctem_tile_yr%vgbiomas_yr_t(i,m) = 0.0
        ctem_tile_yr%totcmass_yr_t(i,m) = 0.0
        ctem_tile_yr%veghght_yr_t(i,m) = 0.0
        ctem_tile_yr%npp_yr_t(i,m) = 0.0
        ctem_tile_yr%gpp_yr_t(i,m) = 0.0
        ctem_tile_yr%nep_yr_t(i,m) = 0.0
        ctem_tile_yr%nbp_yr_t(i,m) = 0.0
        ctem_tile_yr%hetrores_yr_t(i,m) = 0.0
        ctem_tile_yr%autores_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_co2_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_co_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_ch4_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_nmhc_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_h2_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_nox_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_n2o_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_pm25_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_tpm_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_tc_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_oc_yr_t(i,m) = 0.0
        ctem_tile_yr%emit_bc_yr_t(i,m) = 0.0
        ctem_tile_yr%smfuncveg_yr_t(i,m) = 0.0
        ctem_tile_yr%luc_emc_yr_t(i,m) = 0.0
        ctem_tile_yr%lucsocin_yr_t(i,m) = 0.0
        ctem_tile_yr%lucltrin_yr_t(i,m) = 0.0
        ctem_tile_yr%burnfrac_yr_t(i,m) = 0.0
        ctem_tile_yr%bterm_yr_t(i,m) = 0.0
        ctem_tile_yr%lterm_yr_t(i,m) = 0.0
        ctem_tile_yr%mterm_yr_t(i,m) = 0.0
        ctem_tile_yr%ch4WetSpec_yr_t(i,m)  = 0.0
        ctem_tile_yr%wetfdyn_yr_t(i,m)  = 0.0
        ctem_tile_yr%ch4WetDyn_yr_t(i,m)  = 0.0
        ctem_tile_yr%ch4soills_yr_t(i,m)  = 0.0
        ctem_tile_yr%peatdep_yr_t(i,m)  = 0.0
        ctem_tile_yr%fProductDecomp_yr_t(i,m) = 0.0

        do j = 1,icc
          ! per pft
          ctem_yr%laimaxg_yr(i,m,j) = 0.0
          ctem_yr%stemmass_yr(i,m,j) = 0.0
          ctem_yr%rootmass_yr(i,m,j) = 0.0
          ! COMBAK PERLAY
          ctem_yr%litrmass_yr(i,m,j) = 0.0
          ctem_yr%soilcmas_yr(i,m,j) = 0.0
          ctem_yr%litres_yr(i,m,j) = 0.0
          ctem_yr%soilcres_yr(i,m,j) = 0.0
          ! ctem_yr%litrmass_yr(i,m,j,1:ignd)=0.0
          ! ctem_yr%soilcmas_yr(i,m,j,1:ignd)=0.0
          ! ctem_yr%litres_yr(i,m,j,1:ignd)=0.0
          ! ctem_yr%soilcres_yr(i,m,j,1:ignd)=0.0
          ! COMBAK PERLAY
          ctem_yr%vgbiomas_yr(i,m,j) = 0.0
          ctem_yr%totcmass_yr(i,m,j) = 0.0
          ctem_yr%veghght_yr(i,m,j) = 0.0
          ctem_yr%npp_yr(i,m,j) = 0.0
          ctem_yr%gpp_yr(i,m,j) = 0.0
          ctem_yr%nep_yr(i,m,j) = 0.0
          ctem_yr%nbp_yr(i,m,j) = 0.0
          ctem_yr%hetrores_yr(i,m,j) = 0.0
          ctem_yr%autores_yr(i,m,j) = 0.0
          ctem_yr%emit_co2_yr(i,m,j) = 0.0
          ctem_yr%emit_co_yr(i,m,j) = 0.0
          ctem_yr%emit_ch4_yr(i,m,j) = 0.0
          ctem_yr%emit_nmhc_yr(i,m,j) = 0.0
          ctem_yr%emit_h2_yr(i,m,j) = 0.0
          ctem_yr%emit_nox_yr(i,m,j) = 0.0
          ctem_yr%emit_n2o_yr(i,m,j) = 0.0
          ctem_yr%emit_pm25_yr(i,m,j) = 0.0
          ctem_yr%emit_tpm_yr(i,m,j) = 0.0
          ctem_yr%emit_tc_yr(i,m,j) = 0.0
          ctem_yr%emit_oc_yr(i,m,j) = 0.0
          ctem_yr%emit_bc_yr(i,m,j) = 0.0
          ctem_yr%bterm_yr(i,m,j) = 0.0
          ctem_yr%mterm_yr(i,m,j) = 0.0
          ctem_yr%burnfrac_yr(i,m,j) = 0.0
          ctem_yr%smfuncveg_yr(i,m,j) = 0.0
        end do

        ctem_yr%hetrores_yr(i,m,iccp1) = 0.0
        ctem_yr%nep_yr(i,m,iccp1) = 0.0
        ctem_yr%nbp_yr(i,m,iccp1) = 0.0
        ctem_yr%totcmass_yr(i,m,iccp1) = 0.0

        ! COMBAK PERLAY
        ctem_yr%litres_yr(i,m,iccp1) = 0.0
        ctem_yr%soilcres_yr(i,m,iccp1) = 0.0
        ctem_yr%litrmass_yr(i,m,iccp1) = 0.0
        ctem_yr%soilcmas_yr(i,m,iccp1) = 0.0

        ctem_yr%litrmass_yr(i,m,iccp2) = 0.0
        ctem_yr%soilcmas_yr(i,m,iccp2) = 0.0
        ctem_yr%litres_yr(i,m,iccp2) = 0.0
        ctem_yr%soilcres_yr(i,m,iccp2) = 0.0
        ! ctem_yr%litres_yr(i,m,iccp1,1:ignd)=0.0
        ! ctem_yr%soilcres_yr(i,m,iccp1,1:ignd)=0.0
        ! ctem_yr%litrmass_yr(i,m,iccp1,1:ignd)=0.0
        ! ctem_yr%soilcmas_yr(i,m,iccp1,1:ignd)=0.0
        !
        ! ctem_yr%litrmass_yr(i,m,iccp2,1:ignd)=0.0
        ! ctem_yr%soilcmas_yr(i,m,iccp2,1:ignd)=0.0
        ! ctem_yr%litres_yr(i,m,iccp2,1:ignd)=0.0
        ! ctem_yr%soilcres_yr(i,m,iccp2,1:ignd)=0.0
        ! COMBAK PERLAY

      end do ! nmtest
    end do ! nltest

  end subroutine resetYearEnd
  !! @}
  !==================================================
  !> \ingroup ctemstatevars_resetMosaicAccum
  !! @{
  !> Resets physics accumulator variables (used as input to CTEM) after CTEM has been called
  subroutine resetMosaicAccum

    implicit none

    vgat%fsinacc_gat(:) = 0.
    vgat%flinacc_gat(:) = 0.
    vgat%flutacc_gat(:) = 0.
    vgat%preacc_gat(:) = 0.
    ctem_tile%fsnowacc_t(:) = 0.0
    ctem_tile%tcanoaccgat_t(:) = 0.0
    ctem_tile%tcansacc_t(:) = 0.0
    ctem_tile%taaccgat_t(:) = 0.0
    vgat%altotacc_gat(:) = 0.0
    vgat%altotcount_ctm(:) = 0
    ctem_tile%thliqacc_t(:,:) = 0.0
    ctem_tile%thiceacc_t(:,:) = 0.0
    ctem_tile%ancgvgac_t(:,:) = 0.0
    ctem_tile%rmlcgvga_t(:,:) = 0.0

    !-reset peatland accumulators-------------------------------
    ctem_tile%anmossac_t(:)  = 0.0
    ctem_tile%rmlmossac_t(:) = 0.0
    ctem_tile%gppmossac_t(:) = 0.0

  end subroutine resetMosaicAccum
  !! @}
  !=================================================================================


  !> \namespace ctemstatevars
  !> Contains the biogeochemistry-related variable type structures.
  !! 1. c_switch - switches for running CTEM, read from the joboptions file
  !! 2. vrot - CTEM's 'rot' vars
  !! 3. vgat - CTEM's 'gat' vars
  !! 4. ctem_grd - CTEM's grid average variables
  !! 5. ctem_tile - CTEM's variables per tile
  !! 6. ctem_mo - CTEM's variables monthly averaged (per pft)
  !! 7. ctem_grd_mo - CTEM's grid average monthly values
  !! 8. ctem_tile_mo - CTEM's variables per tile monthly values
  !! 9. ctem_yr - CTEM's average annual values (per PFT)
  !! 10. ctem_grd_yr - CTEM's grid average annual values
  !! 11. ctem_tile_yr - CTEM's variables per tile annual values
  !
end module ctemStateVars
