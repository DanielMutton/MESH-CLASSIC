module CLASSICv1_0_mesh

    use mesh_io

    ! Local variables
    integer :: lastDOY             !< Initialized to 365 days, can be overwritten later is leap = true and it is a leap year.
    integer :: JLAT     !< Integer index corresponding to latitude of grid cell
    integer :: NLANDCS  !< Number of modelled areas that contain subareas of canopy over snow
    integer :: NLANDGS  !< Number of modelled areas that contain subareas of snow over bare ground
    integer :: NLANDC   !< Number of modelled areas that contain subareas of canopy over bare ground
    integer :: NLANDG   !< Number of modelled areas that contain subareas of bare ground
    integer :: NLANDI   !< Number of modelled areas that are ice sheets
    character(len = 2), dimension(5) :: CANOPY_LEVEL_MAP = (/'NL', 'BL', 'CR', 'GR', 'UR'/)

    !> Constants.
    !* ICTEM: Number of CTEM vegetation categories (1: not using CTEM).
    ! RETRIEVED FROM RUNCLASS CONSTANTS - DANIEL MUTTON
!    integer, parameter :: ICTEM = 1

!   real, dimension(:), allocatable :: FRZCGAT                              !ADDED BY DANIEL MUTTON FOR RUNCLASS36_CONFIG

    !> MESH variables (frozen soil).
    real, dimension(:), allocatable :: FRZCGAT                              !< Coefficient for frozen soils algorithm.

    !> MESH variables (CLASS variables).
    character(len = 4) :: TITLE1, TITLE2, TITLE3, TITLE4, TITLE5, TITLE6
    character(len = 4) :: NAME1, NAME2, NAME3, NAME4, NAME5, NAME6
    character(len = 4) :: PLACE1, PLACE2, PLACE3, PLACE4, PLACE5, PLACE6
    real, dimension(:), allocatable :: GGEOGRD                              !< Variable to initialize GGEOGAT from file.
    integer, dimension(:), allocatable :: MIDGAT                            !< Vectorized Mosaic ID (MID) flag (to permit masking from MESH inputs).
    integer, dimension(:), allocatable :: IWFGAT                            !< Vectorized IWF flag (to permit different IWF for different tiles or GRUs).

    !> MESH variables (WATROF, LATFLW).
    real, dimension(:), allocatable :: DDGAT                                !< Drainage density (WATROF, LATFLW).
    real, dimension(:), allocatable :: XDGAT                                !< Extinction coefficient (WATROF, LATFLW).
    real, dimension(:), allocatable :: MANNGAT                              !< Manning's 'n' for overland routing (WATROF, LATFLW).
    real, dimension(:), allocatable :: KSGAT                                !< Hydraulic conductivity (WATROF, LATFLW).

    !> MESH variables (PDMROF, LATFLW).
    real, dimension(:), allocatable :: CMINPDM
    real, dimension(:), allocatable :: CMAXPDM
    real, dimension(:), allocatable :: BPDM
    real, dimension(:), allocatable :: K1PDM
    real, dimension(:), allocatable :: K2PDM
    real, dimension(:), allocatable :: ZPNDPRECS, ZPONDPREC, ZPONDPREG, ZPNDPREGS
    real, dimension(:), allocatable :: UM1CS, UM1C, UM1G, UM1GS
    real, dimension(:), allocatable :: QM1CS, QM1C, QM1G, QM1GS
    real, dimension(:), allocatable :: QM2CS, QM2C, QM2G, QM2GS
    real, dimension(:), allocatable :: UMQ
    real, dimension(:), allocatable :: FSTRCS, FSTRC, FSTRG, FSTRGS

    !> MESH variables (ICEBAL).
    real, dimension(:), allocatable :: FREZTH                               !< Freeze threshold (ICEBAL).
    real, dimension(:), allocatable :: SWELIM                               !< SWE limit when snow is convereted to part of the glacier (ICEBAL).
    real, dimension(:), allocatable :: SNDENLIM                             !< Snow density limit when snow is converted to part of the glacier (ICEBAL).
    real, dimension(:), allocatable :: ICE
    real, dimension(:), allocatable :: TICE

    private

    public :: CLASSICv1_0_init, CLASSICv1_0_tile, CLASSICv1_0_finalize

    save

    contains

    subroutine CLASSICv1_0_init(shd, fls)

        use classicParams
        use classStateVars, only: class_gat, allocClassVars
        use ctemStateVars, only: c_switch, vgat, allocCtemVars
        use ctemUtilitiesGAT, only: ctemInitGAT
        use generalUtilsGAT, only: initRandomSeed
        use landuseChange, only: initializeLandCoverGAT
        use applyAllometry, only: allometry

        use PBSM_module
        use mpi_module
        use model_files_variables
        use sa_mesh_common
        use model_dates
        use FLAGS
        use mesh_io_constants
        use input_fields
        use parse_utilities

        implicit none

        !> Arguments.
        type(ShedGridParams) :: shd
        type(fl_ids) :: fls
        integer :: z, l, k, ik, jk, m, j, i, iun, ierr, NMELT
        character(len = DEFAULT_LINE_LENGTH) CLASSOUTTITLE, CLASSOUTNAME, CLASSOUTPLACE, line
        real :: FRAC, factor_a, factor_m, valid_max, valid_min
        logical, allocatable :: value_mask(:)
        integer, dimension(:), allocatable :: INFILTYPE                                                 !previously described in RUNCLASS36 variables.
        real, dimension(:), allocatable :: SI, TSI, SNOWMELTD, SNOWMELTD_LAST, &                        !However, causes errors when added to CLASS variables so I defined within subroutine - Daniel Mutton
            SNOWINFIL, CUMSNOWINFILCS, MELTRUNOFF, CUMSNOWINFILGS, t0_ACC
        real, dimension(:, :), allocatable :: FRZCROW, CMAXROW, CMINROW, BROW, K1ROW, K2ROW             !Type Hydrological Variables from RunCLASS36 Variables - Daniel Mutton
        character(len = DEFAULT_FIELD_LENGTH) code

        !> Return if the process is not marked active.
        if (.false.) then
            return
        else
            call print_new_section("CLASSIC (1.0) is active.")
            call increase_tab()
        end if




		
		
        !> Reset constants, flags and options.
        ican = NO_DATA_INT
        icc = NO_DATA_INT
        l2max = NO_DATA_INT
        c_switch%ctem_on = .false.
        c_switch%spinfast = NO_DATA_INT
        c_switch%useTracer = 0
        c_switch%tracerCO2file = NO_DATA_CHAR
        c_switch%transientCO2 = .false.
        c_switch%CO2File = NO_DATA_CHAR
        c_switch%fixedYearCO2 = NO_DATA_INT
        c_switch%doMethane = .false.
        c_switch%transientCH4 = .false.
        c_switch%CH4File = NO_DATA_CHAR
        c_switch%fixedYearCH4 = NO_DATA_INT
        c_switch%transientOBSWETF = .false.
        c_switch%OBSWETFFile = NO_DATA_CHAR
        c_switch%fixedYearOBSWETF = NO_DATA_INT
        c_switch%dofire = .false.
        c_switch%transientPOPD = .false.
        c_switch%POPDFile = NO_DATA_CHAR
        c_switch%fixedYearPOPD = NO_DATA_INT
        c_switch%transientLGHT= .false.
        c_switch%LGHTFile = NO_DATA_CHAR
        c_switch%fixedYearLGHT = NO_DATA_INT
        c_switch%PFTCompetition = .false.
        c_switch%inibioclim = .false.
        c_switch%start_bare = .false.
        c_switch%lnduseon = .false.
        c_switch%LUCFile = NO_DATA_CHAR
        c_switch%fixedYearLUC = NO_DATA_INT
        c_switch%doChecksums = .false.

        !> Assign constants, flags and options from field list.
        do i = 1, size(input_field_list)
            associate(field => input_field_list(i))
                select case (field%field_name)

                    !> CLASSIC PFT basics (classicPFTbasic).
                    case ('ICAN')
                        call get_field_value(field, ican)
                    case ('ICC')
                        call get_field_value(field, icc)
                    case ('L2MAX')
                        call get_field_value(field, l2max)
                end select
            end associate
        end do

        !> Check to see if CLASSIC run parameters are defined.
        ierr = 0
        call check_field_assigned(ican, 'ICAN', 4, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(icc, 'ICC', 9, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(l2max, 'L2MAX', 3, error_status = z); if (z /= 0) ierr = z
        c_switch%ctem_on = .true.
        if (c_switch%ctem_on) then
            call check_field_assigned(c_switch%spinfast, 'spinfast', 1, error_status = z); if (z /= 0) ierr = z
            if (c_switch%useTracer /= 0) then
                write(code, *) c_switch%useTracer
                call print_remark("The option 'useTracer' is enabled with option '" // trim(adjustl(code)) // "'.")
                call print_warning("The flag 'useTracer' is configured but the feature is not avaiable.")
                call check_field_assigned(c_switch%tracerCO2file, 'tracerCO2file', is_error = .true., error_status = z); if (z /= 0) ierr = z
            end if
            if (c_switch%transientCO2) then
                call print_warning("The flag 'transientCO2' is enabled but the feature is not avaiable.")
                call check_field_assigned(c_switch%CO2File, 'CO2File', is_error = .true., error_status = z); if (z /= 0) ierr = z
            else
                call print_remark("The flag 'transientCO2' is not enabled.")
                call check_field_assigned(c_switch%fixedYearCO2, 'fixedYearCO2', 1850, error_status = z); if (z /= 0) ierr = z
            end if
            if (c_switch%doMethane) then
                call print_remark("The flag 'doMethane' is enabled.")
                if (c_switch%transientCH4) then
                    call print_warning("The flag 'transientCH4' is enabled but the feature is not avaiable.")
                    call check_field_assigned(c_switch%CH4File, 'CH4File', is_error = .true., error_status = z); if (z /= 0) ierr = z
                else
                    call print_remark("The flag 'transientCH4' is not enabled.")
                    call check_field_assigned(c_switch%fixedYearCH4, 'fixedYearCH4', 1850, error_status = z); if (z /= 0) ierr = z
                end if
            end if
            if (c_switch%transientOBSWETF) then
                call print_warning("The flag 'transientOBSWETF' is enabled but the feature is not avaiable.")
                call check_field_assigned(c_switch%OBSWETFFile, 'OBSWETFFile', is_error = .true., error_status = z); if (z /= 0) ierr = z
            else
                call print_remark("The flag 'transientOBSWETF' is not enabled.")
                call check_field_assigned(c_switch%fixedYearOBSWETF, 'fixedYearOBSWETF', -9999, error_status = z); if (z /= 0) ierr = z
            end if
            if (c_switch%dofire) call print_remark("The flag 'dofire' is enabled.")
            if (c_switch%transientPOPD) then
                call print_warning("The flag 'transientPOPD' is enabled but the feature is not avaiable.")
                call check_field_assigned(c_switch%POPDFile, 'POPDFile', is_error = .true., error_status = z); if (z /= 0) ierr = z
            else
                call print_remark("The flag 'transientPOPD' is not enabled.")
                call check_field_assigned(c_switch%fixedYearPOPD, 'fixedYearPOPD', 1850, error_status = z); if (z /= 0) ierr = z
            end if
            if (c_switch%transientLGHT) then
                call print_warning("The flag 'transientLGHT' is enabled but the feature is not avaiable.")
                call check_field_assigned(c_switch%LGHTFile, 'LGHTFile', is_error = .true., error_status = z); if (z /= 0) ierr = z
            else
                call print_remark("The flag 'transientLGHT' is not enabled.")
                call check_field_assigned(c_switch%fixedYearLGHT, 'fixedYearLGHT', 1850, error_status = z); if (z /= 0) ierr = z
            end if
            if (c_switch%PFTCompetition) then
                call print_remark("The flag 'PFTCompetition' is enabled.")
                if (c_switch%inibioclim) call print_remark("The flag 'inibioclim' is enabled.")
                if (c_switch%start_bare) call print_remark("The flag 'start_bare' is enabled.")
            end if
            if (c_switch%lnduseon) then
                call print_warning("The flag 'lnduseon' is enabled but the feature is not avaiable.")
                call check_field_assigned(c_switch%LUCFile, 'LUCFile', is_error = .true., error_status = z); if (z /= 0) ierr = z
            else
                call print_remark("The flag 'lnduseon' is not enabled.")
                call check_field_assigned(c_switch%fixedYearLUC, 'fixedYearLUC', -9999, error_status = z); if (z /= 0) ierr = z
            end if
        end if
        if (c_switch%doChecksums) call print_remark("The flag 'doChecksums' is enabled.")

        !> Check for errors.
        if (ierr /= 0) then
            call reset_tab()
            call print_error("The basic numbers of plant functional types are not assigned or were missing during initialization.")
            call program_abort()
        end if

        !> Set CLASSIC dimensions.
        call readInPFTNums

        !> Allocate the arrays that store the CTEM parameter values
        call allocateParamsCLASSIC

        !> Transfer dimensions from MESH.
        nlat = vs%cell_count
        nmos = vs%gru_count
        ilg = vs%landtile_count
        ignd = vs%sol_count
        delt = ic%dts

        !> Allocate CLASSIC variables.
        call allocClassVars
        call allocCtemVars

        !> Transfer maps.
        class_gat%ILMOS = shd%lc%ILMOS
        class_gat%JLMOS = shd%lc%JLMOS
!        class_gat%IWMOS = shd%wc%ILMOS
!        class_gat%JWMOS = shd%wc%JLMOS

        !> Reset variables.
        VKC = NO_DATA_REAL
        CT = NO_DATA_REAL
        VMIN = NO_DATA_REAL
        TCW = NO_DATA_REAL
        TCICE = NO_DATA_REAL
        TCSAND = NO_DATA_REAL
        TCCLAY = NO_DATA_REAL
        TCOM = NO_DATA_REAL
        TCDRYS = NO_DATA_REAL
        RHOSOL = NO_DATA_REAL
        RHOOM = NO_DATA_REAL
        HCPW = NO_DATA_REAL
        HCPICE = NO_DATA_REAL
        HCPSOL = NO_DATA_REAL
        HCPOM = NO_DATA_REAL
        HCPSND = NO_DATA_REAL
        HCPCLY = NO_DATA_REAL
        SPHW = NO_DATA_REAL
        SPHICE = NO_DATA_REAL
        SPHVEG = NO_DATA_REAL
        RHOW = NO_DATA_REAL
        RHOICE = NO_DATA_REAL
        TCGLAC = NO_DATA_REAL
        CLHMLT = NO_DATA_REAL
        CLHVAP = NO_DATA_REAL
        ZOLNG = NO_DATA_REAL
        ZOLNS = NO_DATA_REAL
        ZOLNI = NO_DATA_REAL
        ZORATG = NO_DATA_REAL
        ALVSI = NO_DATA_REAL
        ALIRI = NO_DATA_REAL
        ALVSO = NO_DATA_REAL
        ALIRO = NO_DATA_REAL
        ALBRCK = NO_DATA_REAL
        GROWYR = NO_DATA_REAL
        ZORAT = NO_DATA_REAL
        CANEXT = NO_DATA_REAL
        XLEAF = NO_DATA_REAL
        RSMN = NO_DATA_REAL
        QA50 = NO_DATA_REAL
        PSGA = NO_DATA_REAL
        PSGB = NO_DATA_REAL
        VPDA = NO_DATA_REAL
        VPDB = NO_DATA_REAL
        THPORG = NO_DATA_REAL
        THRORG = NO_DATA_REAL
        THMORG = NO_DATA_REAL
        BORG = NO_DATA_REAL
        PSISORG = NO_DATA_REAL
        GRKSORG = NO_DATA_REAL
        DELT = NO_DATA_REAL
        ALVSWC = NO_DATA_REAL
        ALIRWC = NO_DATA_REAL
        CXTLRG = NO_DATA_REAL
        ALWV = NO_DATA_REAL
        ALWN = NO_DATA_REAL
        ALDV = NO_DATA_REAL
        ALDN = NO_DATA_REAL
        WSNCAP = NO_DATA_REAL
        REFF0_LAND = NO_DATA_REAL
        ZSNMIN = NO_DATA_REAL
        ZSNMAX2 = NO_DATA_REAL
        modelpft = NO_DATA_INT
        classpfts = NO_DATA_CHAR
        ctempfts = NO_DATA_CHAR
        kn = NO_DATA_REAL
        omega = NO_DATA_REAL
        omega_compete = NO_DATA_REAL
        epsilonl = NO_DATA_REAL
        epsilonl_compete = NO_DATA_REAL
        epsilons = NO_DATA_REAL
        epsilons_compete = NO_DATA_REAL
        epsilonr = NO_DATA_REAL
        epsilonr_compete = NO_DATA_REAL
        rtsrmin = NO_DATA_REAL
        aldrlfon = NO_DATA_REAL
        caleaf = NO_DATA_REAL
        castem = NO_DATA_REAL
        caroot = NO_DATA_REAL
        abar = NO_DATA_REAL
        avertmas = NO_DATA_REAL
        alpha = NO_DATA_REAL
        prcnslai = NO_DATA_REAL
        minslai = NO_DATA_REAL
        mxrtdpth = NO_DATA_REAL
        albvis = NO_DATA_REAL
        albnir = NO_DATA_REAL
        tcoldmin = NO_DATA_REAL
        tcoldmax = NO_DATA_REAL
        twarmmax = NO_DATA_REAL
        gdd5lmt = NO_DATA_REAL
        aridlmt = NO_DATA_REAL
        dryseasonlmt = NO_DATA_REAL
        bio2sap = NO_DATA_REAL
        bioclimrt = NO_DATA_REAL
        grescoef = NO_DATA_REAL
        humicfac = NO_DATA_REAL
        humicfac_bg = NO_DATA_REAL
        laimin = NO_DATA_REAL
        laimax = NO_DATA_REAL
        lambdamax = NO_DATA_REAL
        repro_fraction = NO_DATA_REAL
        bmasthrs_fire = NO_DATA_REAL
        extnmois_veg = NO_DATA_REAL
        extnmois_duff = NO_DATA_REAL
        lwrlthrs = NO_DATA_REAL
        hgrlthrs = NO_DATA_REAL
        parmlght = NO_DATA_REAL
        parblght = NO_DATA_REAL
        reparea = NO_DATA_REAL
        popdthrshld = NO_DATA_REAL
        f0 = NO_DATA_REAL
        maxsprd = NO_DATA_REAL
        frco2glf = NO_DATA_REAL
        frco2blf = NO_DATA_REAL
        frltrglf = NO_DATA_REAL
        frltrblf = NO_DATA_REAL
        frco2stm = NO_DATA_REAL
        frltrstm = NO_DATA_REAL
        frco2rt = NO_DATA_REAL
        frltrrt = NO_DATA_REAL
        frltrbrn = NO_DATA_REAL
        standreplace = NO_DATA_REAL
        emif_co2 = NO_DATA_REAL
        emif_co = NO_DATA_REAL
        emif_ch4 = NO_DATA_REAL
        emif_nmhc = NO_DATA_REAL
        emif_h2 = NO_DATA_REAL
        emif_nox = NO_DATA_REAL
        emif_n2o = NO_DATA_REAL
        emif_pm25 = NO_DATA_REAL
        emif_tpm = NO_DATA_REAL
        emif_tc = NO_DATA_REAL
        emif_oc = NO_DATA_REAL
        emif_bc = NO_DATA_REAL
        bsratelt = NO_DATA_REAL
        bsratesc = NO_DATA_REAL
        tanhq10 = NO_DATA_REAL
        alpha_hetres = NO_DATA_REAL
        bsratelt_g = NO_DATA_REAL
        bsratesc_g = NO_DATA_REAL
        a_hetr = NO_DATA_REAL
        r_depthredu = NO_DATA_REAL
        tcrit = NO_DATA_REAL
        frozered = NO_DATA_REAL
        combust = NO_DATA_REAL
        paper = NO_DATA_REAL
        furniture = NO_DATA_REAL
        bmasthrs = NO_DATA_REAL
        bsrtstem = NO_DATA_REAL
        bsrtstem_compete = NO_DATA_REAL
        bsrtroot = NO_DATA_REAL
        bsrtroot_compete = NO_DATA_REAL
        minlvfr = NO_DATA_REAL
        rmlmoss25 = NO_DATA_REAL
        tau25m = NO_DATA_REAL
        ektau = NO_DATA_REAL
        kc25 = NO_DATA_REAL
        ko25 = NO_DATA_REAL
        ec = NO_DATA_REAL
        ej = NO_DATA_REAL
        eo = NO_DATA_REAL
        evc = NO_DATA_REAL
        sj = NO_DATA_REAL
        hj = NO_DATA_REAL
        alpha_moss = NO_DATA_REAL
        dctmin = NO_DATA_REAL
        dcbaset = NO_DATA_REAL
        bsrateltms = NO_DATA_REAL
        zolnmoss = NO_DATA_REAL
        thpmoss = NO_DATA_REAL
        thrmoss = NO_DATA_REAL
        thmmoss = NO_DATA_REAL
        bmoss = NO_DATA_REAL
        psismoss = NO_DATA_REAL
        grksmoss = NO_DATA_REAL
        hcpmoss = NO_DATA_REAL
        grescoefmoss = NO_DATA_REAL
        rmortmoss = NO_DATA_REAL
        humicfacmoss = NO_DATA_REAL
        dayschk = NO_DATA_INT
        drgta = NO_DATA_REAL
        eta = NO_DATA_REAL
        kappa = NO_DATA_REAL
        flhrspan = NO_DATA_REAL
        fracbofg = NO_DATA_REAL
        harvthrs = NO_DATA_REAL
        specsla = NO_DATA_REAL
        thrprcnt = NO_DATA_REAL
        lwrthrsh = NO_DATA_REAL
        cdlsrtmx = NO_DATA_REAL
        kmort1 = NO_DATA_REAL
        mxmortge = NO_DATA_REAL
        mxmortge_compete = NO_DATA_REAL
        maxage = NO_DATA_REAL
        maxage_compete = NO_DATA_REAL
        lfespany = NO_DATA_REAL
        drlsrtmx = NO_DATA_REAL
        drlsrtmx_compete = NO_DATA_REAL
        colda = NO_DATA_REAL
        coldlmt = NO_DATA_INT
        coldthrs = NO_DATA_REAL
        roothrsh = NO_DATA_REAL
        cryodiffus = NO_DATA_REAL
        biodiffus = NO_DATA_REAL
        kterm = NO_DATA_REAL
        eftime = NO_DATA_REAL
        D_air = NO_DATA_REAL
        g_0 = NO_DATA_REAL
        betaCH4 = NO_DATA_REAL
        k_o = NO_DATA_REAL
        stemlife = NO_DATA_REAL
        rootlife = NO_DATA_REAL
        stmhrspn = NO_DATA_REAL
        ratioch4 = NO_DATA_REAL
        wtdryres = NO_DATA_REAL
        lat_thrshld1 = NO_DATA_REAL
        lat_thrshld2 = NO_DATA_REAL
        soilw_thrshN = NO_DATA_REAL
        soilw_thrshE = NO_DATA_REAL
        soilw_thrshS = NO_DATA_REAL
        kn = NO_DATA_REAL
        tup = NO_DATA_REAL
        tlow = NO_DATA_REAL
        alpha_phtsyn = NO_DATA_REAL
        omega_phtsyn = NO_DATA_REAL
!?        isc4 = NO_DATA_INT
        mm = NO_DATA_REAL
        bb = NO_DATA_REAL
        vpd0 = NO_DATA_REAL
        sn = NO_DATA_INT
        smscale = NO_DATA_REAL
        vmax = NO_DATA_REAL
        reqiter = NO_DATA_INT
        co2imax = NO_DATA_REAL
        beta1 = NO_DATA_REAL
        beta2 = NO_DATA_REAL
        inico2i = NO_DATA_REAL
        chi = NO_DATA_REAL
        rmlcoeff = NO_DATA_REAL
        gamma_w = NO_DATA_REAL
        gamma_m = NO_DATA_REAL

        !> Assign variables from field list.
        do i = 1, size(input_field_list)
            associate(field => input_field_list(i))
                select case (field%field_name)
                end select
            end associate
        end do

        !> Check to see if CLASSIC run parameters are defined.
        ierr = 0
        call check_field_assigned(VKC, 'VKC', 0.4, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(CT   , 'CT', 1.15E-3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(VMIN , 'VMIN', 0.1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(TCW  , 'TCW', 0.57, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(TCICE, 'TCICE', 2.24, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(TCSAND, 'TCSAND', 2.5, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(TCCLAY, 'TCCLAY', 2.5, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(TCOM  , 'TCOM', 0.25, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(TCDRYS, 'TCDRYS', 0.275, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(RHOSOL, 'RHOSOL', 2.65E3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(RHOOM , 'RHOOM', 1.3E3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(HCPW  , 'HCPW', 4.187E6, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(HCPICE, 'HCPICE', 1.9257E6, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(HCPSOL, 'HCPSOL', 2.25E6, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(HCPOM , 'HCPOM', 2.5E6, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(HCPSND, 'HCPSND', 2.13E6, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(HCPCLY, 'HCPCLY', 2.38E6, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(SPHW  , 'SPHW', 4.186E3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(SPHICE, 'SPHICE', 2.1E3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(SPHVEG, 'SPHVEG', 2.7E3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(RHOW  , 'RHOW', 1.0E3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(RHOICE, 'RHOICE', 0.917E3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(TCGLAC, 'TCGLAC', 2.24, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(CLHMLT, 'CLHMLT', 0.334E6, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(CLHVAP, 'CLHVAP', 2.501E6, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ZOLNG  , 'ZOLNG', -4.605, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ZOLNS  , 'ZOLNS', -6.908, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ZOLNI  , 'ZOLNI', -6.215, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ZORATG , 'ZORATG', 3.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ALVSI  , 'ALVSI', 0.95, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ALIRI  , 'ALIRI', 0.73, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ALVSO  , 'ALVSO', 0.05, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ALIRO  , 'ALIRO', 0.3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ALBRCK , 'ALBRCK', 0.16, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(GROWYR , 'GROWYR', reshape((/ &
            213.0, 213.0, 213.0, 213.0, 213.0,  213.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, 75.0, 106.0, 136.0, 167.0, 167.0, 167.0, &
            273.0, 273.0, 273.0, 273.0, 273.0, 273.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, 135.0, 166.0,196.0, 196.0, 196.0, 196.0, &
            121.0, 121.0, 121.0, 121.0, 121.0, 121.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, 275.0, 244.0,214.0, 214.0, 214.0, 214.0, &
            151.0, 151.0, 151.0, 151.0, 151.0, 151.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, 305.0, 274.0,244.0, 244.0, 244.0, 244.0, &
            213.0, 213.0, 213.0, 213.0, 213.0, 213.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 75.0, 106.0, 136.0, 167.0, 167.0, 167.0, 167.0, &
            273.0, 273.0, 273.0, 273.0, 273.0, 273.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 135.0, 166.0, 196.0, 196.0, 196.0, 196.0, 196.0, &
            121.0, 121.0, 121.0, 121.0, 121.0, 121.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 275.0, 244.0, 214.0, 214.0, 214.0, 214.0, 214.0, &
            151.0, 151.0, 151.0, 151.0, 151.0, 151.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 305.0, 274.0, 244.0, 244.0, 244.0, 244.0, 244.0/), (/18, 4, 2/)), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ZORAT  , 'ZORAT', (/1.0, 1.0, 1.0, 1.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(CANEXT , 'CANEXT', (/-0.5, -1.5, -0.8, -0.8/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(XLEAF  , 'XLEAF', (/0.0247, 0.0204, 0.0456, 0.0456/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(RSMN, 'RSMN', (/200.0, 125.0, 85.0, 100.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(QA50, 'QA50', (/30.0,  40.0, 30.0, 30.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(PSGA, 'PSGA', (/100.0, 100.0, 100.0, 100.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(PSGB, 'PSGB', (/5.0, 5.0, 5.0, 5.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(VPDA, 'VPDA', (/0.65, 0.5, 0.5, 0.5/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(VPDB, 'VPDB', (/1.05, 0.6, 1.0, 1.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(THPORG , 'THPORG', (/0.93, 0.88, 0.83/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(THRORG , 'THRORG', (/0.275, 0.620, 0.705/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(THMORG , 'THMORG', (/0.04, 0.15, 0.22/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(BORG , 'BORG', (/2.7, 6.1, 12.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(PSISORG, 'PSISORG', (/0.0103, 0.0102, 0.0101/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(GRKSORG, 'GRKSORG', (/2.8E-4, 2.0E-6, 1.0E-7/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(DELT, 'DELT', 1800.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ALVSWC, 'ALVSWC', 0.27, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ALIRWC, 'ALIRWC', 0.38, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(CXTLRG, 'CXTLRG', 1.0E20, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ALWV, 'ALWV', (/ &
            0.25, 0.23, 0.21, 0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, &
            0.13, 0.12, 0.11, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ALWN, 'ALWN', (/ &
            0.5, 0.46, 0.42, 0.4, 0.38, 0.36, 0.34, 0.32, 0.3, 0.28, &
            0.26, 0.24, 0.22, 0.2, 0.18, 0.16, 0.14, 0.12, 0.1, 0.08/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ALDV, 'ALDV', (/ &
            0.36, 0.34, 0.32, 0.31, 0.3, 0.29, 0.28, 0.27, 0.26, 0.25, &
            0.24, 0.23, 0.22, 0.2, 0.18, 0.16, 0.14, 0.12, 0.1, 0.08/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ALDN, 'ALDN', (/ &
            0.61, 0.57, 0.53, 0.51, 0.49, 0.48, 0.45, 0.43, 0.41, 0.39, &
            0.37, 0.35, 0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.16/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(WSNCAP, 'WSNCAP', 0.04, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(REFF0_LAND, 'REFF0_LAND', 5.452563E-5, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ZSNMIN, 'ZSNMIN', 5.E-2, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ZSNMAX2, 'ZSNMAX2', 25.E-2, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(modelpft, 'MODELPFT', (/1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(classpfts, 'CLASSPFTS', (/'NdlTr', 'BdlTr', 'Crops', 'Grass'/), error_status = z); if (z /= 0) ierr = z
        !Moved classic v1.5 PFTs to the end, as they're considered inactive for these runs
		call check_field_assigned(ctempfts, 'CTEMPFTS', (/ &
            'NdlEvgTr', 'NdlDcdTr', 'BdlEvgTr', 'BdlDCoTr', 'BdlDDrTr', &
            'CropC3  ', 'CropC4  ', &
            'GrassC3 ', 'GrassC4 ', &
			'BdlDCoSh', 'BdlEvgSh','Sedge   '/), &
            error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(kn, 'KN', (/0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.40, 0.48, 0.0, 0.46, 0.44, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(omega, 'OMEGA', (/0.8, 0.5, 0.0, 0.8, 0.8, 0.8, 0.05, 0.05, 0.0, 1.0, 1.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(omega_compete, 'OMEGA_COMPETE', (/0.8, 0.5, 0.0, 0.8, 0.45, 0.8, 0.05, 0.05, 0.0, 1.0, 1.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(epsilonl, 'EPSILONL', (/0.2, 0.06, 0.0, 0.35, 0.35, 0.25, 0.8, 0.8, 0.0, 0.01, 0.01, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(epsilonl_compete, 'EPSILONL_COMPETE', (/0.19, 0.45, 0.0, 0.39, 0.5, 0.3, 0.8, 0.8, 0.0, 0.1, 0.1, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(epsilons, 'EPSILONS', (/0.15, 0.05, 0.0, 0.05, 0.1, 0.1, 0.15, 0.15, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(epsilons_compete, 'EPSILONS_COMPETE', (/0.4, 0.34, 0.0, 0.21, 0.35, 0.1, 0.15, 0.15, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(epsilonr, 'EPSILONR', (/0.65, 0.89, 0.0, 0.6, 0.55, 0.65, 0.05, 0.05, 0.0, 0.99, 0.99, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(epsilonr_compete, 'EPSILONR_COMPETE', (/0.41, 0.21, 0.0, 0.4, 0.15, 0.6, 0.05, 0.05, 0.0, 0.9, 0.9, 0.0/), error_status = z); if (z /= 0) ierr = z
        if (ierr == 0) then
            i = 0
            if (c_switch%PFTCompetition) then

                !> Check competition allocation parameters.
                i = count(abs(epsilons_compete + epsilonr_compete + epsilonl_compete - 1.0) > abszero)
                if (i > 0) then
                    where (abs(epsilons_compete + epsilonr_compete + epsilonl_compete - 1.0) > abszero)
                        epsilonl_compete = 1.0 - epsilons_compete - epsilonr_compete
                    end where
                end if
            else

                !> Check non-competition allocation parameters.
                i = count(abs(epsilonl + epsilons + epsilonr - 1.0) > abszero)
                if (i > 0) then
                    where (abs(epsilonl + epsilons + epsilonr - 1.0) > abszero)
                        epsilonl = 1.0 - epsilons - epsilonr
                    end where
                end if
            end if
            if (i > 0) then
                write(code, *) i
                call print_remark( &
                    "Adjusting where the sum of epsilon for stem, root and leaf do not add up to 1.0 (" // trim(adjustl(code)) // &
                    " PFTs).")
            end if
        end if
        call check_field_assigned(rtsrmin, 'RTSRMIN', (/0.16, 0.16, 0.0, 0.16, 0.16, 0.32, 0.16, 0.16, 0.0, 0.5, 0.5, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(aldrlfon, 'ALDRLFON', (/1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(caleaf, 'CALEAF', (/0.275, 0.3, 0.0, 0.2, 0.25, 0.25, 0.4, 0.4, 0.0, 0.45, 0.45, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(castem, 'CASTEM', (/0.475, 0.45, 0.0, 0.37, 0.4, 0.4, 0.15, 0.15, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(caroot, 'CAROOT', (/0.25, 0.25, 0.0, 0.43, 0.35, 0.35, 0.45, 0.45, 0.0, 0.55, 0.55, 0.0/), error_status = z); if (z /= 0) ierr = z
		call check_field_assigned(abar, 'ABAR', (/4.70, 5.86, 0.0, 3.87, 3.46, 3.97, 3.97, 3.97, 0.0, 5.86, 4.92, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(avertmas, 'AVERTMAS', (/1.85, 1.45, 0.0, 2.45, 2.1, 2.1, 0.1, 0.1, 0.0, 0.7, 0.7, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(alpha, 'ALPHA', (/0.8, 0.8, 0.0, 0.8, 0.8, 0.8, 0.8, 0.8, 0.0, 0.8, 0.8, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(prcnslai, 'PRCNSLAI', (/7.5, 7.5, 0.0, 7.5, 7.5, 7.5, 7.5, 7.5, 0.0, 2.5, 2.5, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(minslai, 'MINSLAI', (/0.3, 0.3, 0.0, 0.3, 0.3, 0.3, 0.2, 0.2, 0.0, 0.2, 0.2, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(mxrtdpth, 'MXRTDPTH', (/3.0, 3.0, 0.0, 5.0, 5.0, 3.0, 2.0, 2.0, 0.0, 1.0, 1.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(albvis, 'ALBVIS', (/3.0, 3.0, 0.0, 3.0, 5.0, 5.0, 5.5, 5.5, 0.0, 5.0, 6.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(albnir, 'ALBNIR', (/19.0, 19.0, 0.0, 23.0, 29.0, 29.0, 34.0, 34.0, 0.0, 30.0, 34.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(tcoldmin, 'TCOLDMIN', (/-999.9, -999.9, 0.0, 2.5, -35.0, 4.0, -999.9, -999.9, 0.0, -999.9, -999.9, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(tcoldmax, 'TCOLDMAX', (/18.0, -28.0, 0.0, 999.9, 16.0, 900.0, 999.9, 999.9, 0.0, 999.9, 999.9, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(twarmmax, 'TWARMMAX', (/99.9, 25.0, 0.0, 99.9, 99.9, 99.9, 99.9, 99.9, 0.0, 99.9, 99.9, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(gdd5lmt, 'GDD5LMT', (/375.0, 600.0, 0.0, 1200.0, 300.0, 9.9, 9.9, 9.9, 0.0, 9.9, 9.9, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(aridlmt, 'ARIDLMT', (/9.9, 9.9, 0.0, 9.9, 9.9, 0.9, 9.9, 9.9, 0.0, 9.9, 9.9, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(dryseasonlmt, 'DRYSEASONLMT', (/9.0, 99.9, 0.0, 99.9, 99.9, 5.5, 99.9, 99.9, 0.0, 99.9, 99.9, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bio2sap, 'BIO2SAP', (/0.32, 0.2, 0.0, 0.08, 0.14, 0.13, 0.0, 0.0, 0.0, 0.2, 0.2, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bioclimrt, 'BIOCLIMRT', 0.25, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(grescoef, 'GRESCOEF', (/0.15, 0.15, 0.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.0, 0.15, 0.15, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(humicfac, 'HUMICFAC', (/0.42, 0.42, 0.0, 0.53, 0.48, 0.48, 0.1, 0.1, 0.0, 0.42, 0.42, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(humicfac_bg, 'HUMICFAC_BG', 0.45, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(laimin, 'LAIMIN', (/1.0, 1.0, 0.0, 1.5, 1.0, 1.0, 1.0, 1.0,  0.0, 0.01, 0.01, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(laimax, 'LAIMAX', (/4.0, 3.0, 0.0, 6.0, 5.0, 5.0, 8.0, 8.0, 0.0, 4.0, 4.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(lambdamax, 'LAMBDAMAX', 0.1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(repro_fraction, 'REPRO_FRACTION', 0.1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bmasthrs_fire, 'BMASTHRS_FIRE', (/0.4, 1.2/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(extnmois_veg, 'EXTNMOIS_VEG', 0.3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(extnmois_duff, 'EXTNMOIS_DUFF', 0.5, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(lwrlthrs, 'LWRLTHRS', 0.25, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(hgrlthrs, 'HGRLTHRS', 10.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(parmlght, 'PARMLGHT', 0.8, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(parblght, 'PARBLGHT', 0.1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(reparea, 'REPAREA', 500.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(popdthrshld, 'POPDTHRSHLD', 300.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(f0, 'F0', 0.05, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(maxsprd, 'MAXSPRD', (/0.38, 0.38, 0.0, 0.28, 0.28, 0.28, 0.0, 0.0, 0.0, 0.51, 0.75, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(frco2glf, 'FRCO2GLF', (/0.42, 0.42, 0.0, 0.42, 0.42, 0.42, 0.0, 0.0, 0.0, 0.48, 0.48, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(frco2blf, 'FRCO2BLF', (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.54, 0.54, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(frltrglf, 'FRLTRGLF', (/0.2, 0.2, 0.0, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.10, 0.10, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(frltrblf, 'FRLTRBLF', (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.06, 0.06, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(frco2stm, 'FRCO2STM', (/0.12, 0.12, 0.0, 0.12, 0.06, 0.06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(frltrstm, 'FRLTRSTM', (/0.6, 0.6, 0.0, 0.6, 0.4, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(frco2rt, 'FRCO2RT', (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(frltrrt, 'FRLTRRT', (/0.1, 0.1, 0.0, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.25, 0.25, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(frltrbrn, 'FRLTRBRN', (/0.3, 0.3, 0.0, 0.36, 0.36, 0.36, 0.0, 0.0, 0.0, 0.42, 0.42, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(standreplace, 'STANDREPLACE', (/0.2, 0.2, 0.0, 0.5, 0.2, 0.15, 0.0, 0.0, 0.0, 0.25, 0.25, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_co2, 'EMIF_CO2', (/1576.0, 1576.0, 0.0, 1604.0, 1576.0, 1654.0, 1576.0, 1654.0, 0.0, 1576.0, 1654.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_co, 'EMIF_CO', (/106.0, 106.0, 0.0, 103.0, 106.0, 64.0, 106.0,  64.0, 0.0, 106.0, 64.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_ch4, 'EMIF_CH4', (/4.8, 4.8, 0.0, 5.8, 4.8, 2.4, 4.8, 2.4, 0.0, 4.8, 2.4, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_nmhc, 'EMIF_NMHC', (/5.7, 5.7, 0.0, 6.4, 5.7, 3.7, 5.7, 3.7, 0.0, 5.7, 3.7, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_h2, 'EMIF_H2', (/1.8, 1.8, 0.0, 2.54, 1.8, 0.98, 1.8, 0.98, 0.0, 1.8, 0.98, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_nox, 'EMIF_NOX', (/3.24, 3.24, 0.0, 2.90, 3.24, 2.49, 3.24, 2.49, 0.0, 3.24, 2.49, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_n2o, 'EMIF_N2O', (/0.26, 0.26, 0.0, 0.23, 0.26, 0.2, 0.26, 0.2, 0.0, 0.26, 0.2, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_pm25, 'EMIF_PM25', (/12.7, 12.7, 0.0, 10.5, 12.7, 5.2, 12.7, 5.2, 0.0, 12.7, 5.2, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_tpm, 'EMIF_TPM', (/17.6, 17.6, 0.0, 14.7, 17.6, 8.5, 17.6, 8.5, 0.0, 17.6, 8.5, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_tc, 'EMIF_TC', (/8.3, 8.3, 0.0, 7.2, 8.3, 3.4, 8.3, 3.4, 0.0, 8.3, 3.4, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_oc, 'EMIF_OC', (/9.1, 9.1, 0.0, 6.7, 9.1, 3.2, 9.1, 3.2, 0.0, 9.1, 3.2, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(emif_bc, 'EMIF_BC', (/0.56, 0.56, 0.0, 0.56, 0.56, 0.47, 0.56, 0.47, 0.0, 0.56, 0.47, 0.0/), error_status = z); if (z /= 0) ierr = z
		call check_field_assigned(bsratelt, 'BSRATELT', (/0.4, 0.4, 0.0, 0.6, 0.7, 0.7, 0.4, 0.3, 0.0, 0.526, 0.526, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bsratesc, 'BSRATESC', (/0.026, 0.026, 0.0, 0.0208, 0.0208, 0.0208, 0.035, 0.035, 0.0, 0.0125, 0.0125, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(tanhq10, 'TANHQ10', (/2.0, 0.0, 0.075, 46.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(alpha_hetres, 'ALPHA_HETRES', 0.7, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bsratelt_g, 'BSRATELT_G', 0.5605, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bsratesc_g, 'BSRATESC_G', 0.02258, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(a_hetr, 'A_HETR', 4.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(r_depthredu, 'R_DEPTHREDU', 10.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(tcrit, 'TCRIT', -1.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(frozered, 'FROZERED', 0.1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(combust, 'COMBUST', (/0.15, 0.3, 0.45/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(paper, 'PAPER', (/0.7, 0.7, 0.55/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(furniture, 'FURNITURE', (/0.15, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bmasthrs, 'BMASTHRS', (/4.0, 1.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bsrtstem, 'BSRTSTEM', (/0.09, 0.055, 0.0, 0.06, 0.0335, 0.03, 0.0365, 0.0365, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bsrtstem_compete, 'BSRTSTEM_COMPETE', (/0.07, 0.055, 0.0, 0.05, 0.0335, 0.035, 0.0365, 0.0365, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bsrtroot, 'BSRTROOT', (/0.5, 0.285, 0.0, 0.65, 0.225, 0.055, 0.16, 0.16, 0.0, 0.1, 0.1, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bsrtroot_compete, 'BSRTROOT_COMPETE', (/0.5, 0.285, 0.0, 0.4, 0.225, 0.15, 0.16, 0.16, 0.0, 0.1, 0.1, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(minlvfr, 'MINLVFR', 0.05, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(rmlmoss25, 'RMLMOSS25', 1.1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(tau25m, 'TAU25M', 2321.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ektau, 'EKTAU', -29000.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(kc25, 'KC25', 30.2, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ko25, 'KO25', 24800.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ec, 'EC', 63500.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ej, 'EJ', 37000.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(eo, 'EO', 35000.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(evc, 'EVC', 53000.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(sj, 'SJ', 719.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(hj, 'HJ', 220300.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(alpha_moss, 'ALPHA_MOSS', 0.21, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(dctmin, 'DCTMIN', 269.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(dcbaset, 'DCBASET', 283.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bsrateltms, 'BSRATELTMS', 0.2, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(zolnmoss, 'ZOLNMOSS', -6.57, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(thpmoss, 'THPMOSS', 0.98, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(thrmoss, 'THRMOSS', 0.2, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(thmmoss, 'THMMOSS', 0.01, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bmoss, 'BMOSS', 2.3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(psismoss, 'PSISMOSS', 0.0103, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(grksmoss, 'GRKSMOSS', 1.83E-3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(hcpmoss, 'HCPMOSS', 2.5E6, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(grescoefmoss, 'GRESCOEFMOSS', 0.15, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(rmortmoss, 'RMORTMOSS', 0.7, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(humicfacmoss, 'HUMICFACMOSS', 2.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(dayschk, 'DAYSCHK', (/7, 7, 0, 7, 7, 7, 7, 7, 0, 7, 7, 0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(drgta, 'DRGTA', (/3.0, 3.0, 0.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.0, 3.0, 3.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        !DM changes from 7.0 to 3
		call check_field_assigned(eta, 'ETA', (/10.0, 30.8, 0.0, 31.0, 50.0, 30.0, 3.0, 3.0, 0.0, 3.0, 3.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(kappa, 'KAPPA', (/1.6, 1.6, 0.0, 1.6, 1.6, 1.6, 1.6, 1.6, 0.0, 1.2, 1.2, 0.0/), error_status = z); if (z /= 0) ierr = z
        !DM changed from /17.0, 45.0/
		call check_field_assigned(flhrspan, 'FLHRSPAN', (/15.0, 30.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(fracbofg, 'FRACBOFG', 0.55, error_status = z); if (z /= 0) ierr = z
		!DM raised to 5.5 from 4.5
		call check_field_assigned(harvthrs, 'HARVTHRS', (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.5, 3.5, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        !DM raised from 0 to 30
		call check_field_assigned(specsla, 'SPECSLA', (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(thrprcnt, 'THRPRCNT', (/40.0, 40.0, 0.0, 40.0, 50.0, 50.0, 50.0, 50.0, 0.0, 40.0, 40.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(lwrthrsh, 'LWRTHRSH', (/-50.0, -5.0, 0.0, 5.0, 8.0, 5.0, 5.0, 5.0, 0.0, 0.1, 5.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(cdlsrtmx, 'CDLSRTMX', (/0.1, 0.3, 0.0, 0.3, 0.4, 0.15, 0.15, 0.15, 0.0, 0.15, 0.15, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(kmort1, 'KMORT1', 0.3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(mxmortge, 'MXMORTGE', (/0.005, 0.005, 0.0, 0.005, 0.005, 0.005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(mxmortge_compete, 'MXMORTGE_COMPETE', (/0.005, 0.005, 0.0, 0.005, 0.005, 0.005, 0.0,  0.0,  0.0, 0.05,  0.1, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(maxage, 'MAXAGE', (/250.0, 400.0, 0.0, 600.0, 250.0, 500.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(maxage_compete, 'MAXAGE_COMPETE', (/800.0, 500.0, 0.0, 700.0, 450.0, 500.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(lfespany, 'LFESPANY', (/5.0, 1.0, 0.0, 1.5, 1.0, 1.0, 1.75, 1.75, 0.0, 1.0, 1.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(drlsrtmx, 'DRLSRTMX', (/0.0025, 0.005, 0.0, 0.005, 0.005, 0.025, 0.005, 0.005, 0.0, 0.05, 0.05, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(drlsrtmx_compete, 'DRLSRTMX_COMPETE', (/0.006, 0.005, 0.0, 0.01, 0.025, 0.03, 0.005, 0.005, 0.0, 0.02, 0.02, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(colda, 'COLDA', (/3.0, 3.0, 0.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.0, 3.0, 3.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(coldlmt, 'COLDLMT', (/7, 5/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(coldthrs, 'COLDTHRS', (/-5.0 , 3.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(roothrsh, 'ROOTHRSH', 8.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(cryodiffus, 'CRYODIFFUS', 1.36986E-6, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(biodiffus, 'BIODIFFUS', 2.73972E-7, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(kterm, 'KTERM', 3.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(eftime, 'EFTIME', 10.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(D_air, 'D_AIR', 0.196, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(g_0, 'G_0', 6.79051E-3, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(betaCH4, 'BETACH4', 0.8, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(k_o, 'K_O', 5.03E-5, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(stemlife, 'STEMLIFE', (/86.3, 86.3, 0.0, 80.5, 80.5, 75.8, 20.0, 20.0, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(rootlife, 'ROOTLIFE', (/13.8,13.2, 0.0, 12.7,10.9, 9.8, 3.0, 3.0, 0.0, 3.0, 3.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(stmhrspn, 'STMHRSPN', 17.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ratioch4, 'RATIOCH4', 0.125, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(wtdryres, 'WTDRYRES', (/0.23, 0.37, 0.23/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(lat_thrshld1, 'LAT_THRSHLD1', 40.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(lat_thrshld2, 'LAT_THRSHLD2', -35.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(soilw_thrshN, 'SOILW_THRSHN', (/0.45, 0.65/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(soilw_thrshE, 'SOILW_THRSHE', (/0.55, 0.85/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(soilw_thrshS, 'SOILW_THRSHS', (/0.70, 0.85/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(tup, 'TUP', (/307.1, 307.1, 0.0, 318.1, 310.1, 310.1, 315.1, 315.1, 0.0, 313.1, 323.1, 0.0/), error_status = z); if (z /= 0) ierr = z
        !DM raised Tlow(1) and (7) by 5
		call check_field_assigned(tlow, 'TLOW', (/283.1, 268.1, 0.0, 273.1, 273.1, 273.1, 283.1, 278.1, 0.0, 272.1, 283.1, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(alpha_phtsyn, 'ALPHA_PHTSYN', (/0.08, 0.08, 0.0, 0.08, 0.08, 0.08, 0.08, 0.04, 0.0, 0.08, 0.04, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(omega_phtsyn, 'OMEGA_PHTSYN', (/0.15, 0.15, 0.0, 0.15, 0.15, 0.15, 0.15, 0.17, 0.0, 0.15, 0.17, 0.0/), error_status = z); if (z /= 0) ierr = z
!?        call check_field_assigned(isc4, 'ISC4', (/'off', 'off', 'off', 'off', 'off', 'off', 'off', 'on ', 'off', 'off', 'on ', 'off'/), error_status = z); if (z /= 0) ierr = z
        isc4 = (/.false., .false., .false., .false., .false., .false., .false., .true., .false., .false., .true., .false./)
        call check_field_assigned(mm, 'MM', (/9.0, 9.0, 0.0, 12.0,12.0,12.0, 12.0, 6.0, 0.0, 12.0, 6.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(bb, 'BB', (/0.01, 0.01, 0.0, 0.01, 0.01, 0.01, 0.01, 0.04, 0.0, 0.01, 0.04, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(vpd0, 'VPD0', (/2000.0, 2000.0, 0.0, 2000.0, 2000.0, 2000.0, 1500.0, 1500.0, 0.0, 1500.0, 1500.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        !DM changed from 2 to 6
		call check_field_assigned(sn, 'SN', (/2, 2, 0, 4, 2, 2, 6, 2, 0, 2, 2, 0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(smscale, 'SMSCALE', (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(vmax, 'VMAX', (/75.0E-6, 47.0E-6, 0.0, 35.0E-6, 75.0E-6, 40.0E-6, 75.0E-6, 40.0E-6, 0.0, 55.0E-6, 40.0E-6, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(reqiter, 'REQITER', 10, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(co2imax, 'CO2IMAX', 2000.00, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(beta1, 'BETA1', 0.95, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(beta2, 'BETA2', 0.99, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(inico2i, 'INICO2I', (/0.65, 0.65, 0.0, 0.65, 0.65, 0.65, 0.65, 0.37, 0.0, 0.65, 0.37, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(chi, 'CHI', (/0.01,  0.01, 0.0, 0.17,  0.17, 0.17, -0.30, -0.30, 0.0, -0.30, -0.30, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(rmlcoeff, 'RMLCOEFF', (/0.015, 0.021, 0.0, 0.025, 0.015, 0.015, 0.015, 0.025, 0.0, 0.013, 0.025, 0.0/), error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(gamma_w, 'GAMMA_W', 0.35, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(gamma_m, 'GAMMA_M', 0.95, error_status = z); if (z /= 0) ierr = z

		print *, "printing Vmax: ", vmax
        !> Check for errors.
        if (ierr /= 0) then
            call reset_tab()
            call print_error("Certain CLASSIC constants or parameters are not assigned or were missing during initialization.")
            call program_abort()
        end if

        !> Initialize the CTEM parameters.
        call readInParams

        !> Allocate MESH variables.
        allocate( &
            FRZCGAT(ilg), GGEOGRD(1), MIDGAT(ilg), IWFGAT(ilg), DDGAT(ilg), XDGAT(ilg), MANNGAT(ilg), KSGAT(ilg), &
            FREZTH(ilg), SWELIM(ilg), SNDENLIM(ilg), ICE(ilg), TICE(ilg))
        FRZCGAT = 0.0
        GGEOGRD = 0.0
        MIDGAT = 0
        IWFGAT = 1
        DDGAT = 0.0
        XDGAT = 0.0
        MANNGAT = 0.0
        KSGAT = 0.0
        FREZTH = 0.0
        SWELIM = 0.0
        SNDENLIM = 0.0

    JLAT = 1
    c_switch%leap = .true.
    lastDOY = 365

    vgat%wetfrac_pres = -9999. !< If transientOBSWETF or fixedYearOBSWETF != -9999 this variable will be overwritten with

    !RUNCLASS36_config init subroutine translated - Daniel Mutton

    ierr = 0
        if (.not. allocated(vs%tile%fsin)) then
            call print_error("The driving variable '" // VN_FSIN // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (.not. allocated(vs%tile%flin)) then
            call print_error("The driving variable '" // VN_FLIN // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (.not. allocated(vs%tile%ta)) then
            call print_error("The driving variable '" // VN_TA // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (.not. allocated(vs%tile%qa)) then
            call print_error("The driving variable '" // VN_QA // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (.not. allocated(vs%tile%pres)) then
            call print_error("The driving variable '" // VN_PRES // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (.not. allocated(vs%tile%uv)) then
            call print_error("The driving variable '" // VN_UV // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (c_switch%IPCP == 4) then
            if (.not. allocated(vs%tile%prern) .or. .not. allocated(vs%tile%presno)) then
                call print_error( &
                    "'IPCP 4' is active but the driving variables '" // VN_PRERN // "' and '" // VN_PRESNO // &
                    "' are not active or not associated with an input file.")
                ierr = 1
            else if (allocated(vs%tile%pre)) then
                call print_info( &
                    "'IPCP 4' is active with the '" // VN_PRERN // "' and '" // VN_PRESNO // "' variables. The '" // VN_PRE // &
                    "' variable is also active but inputs on the field are not being used.")
            end if
        else if (.not. allocated(vs%tile%pre)) then
            call print_error("The driving variable '" // VN_PRE // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (ierr /= 0) then
            call reset_tab()
            call print_error( &
                "The variables required to drive the model are not active or have not been associated with an input file.")
            call program_abort()
        end if

!       NA = shd%NA
!        NTYPE = shd%lc%NTYPE
!        NSL = shd%lc%IGND
!        NML = shd%lc%NML
        DELT = ic%dts
        if (DELT > 1800) then

        !> Print a warning if the time-step exceeds recommendations from the CLASS manual.
        write(line, *) DELT*60
        call print_warning( &
            "CLASS is designed to run at a time-step of 30 minutes or less. The current time-step of " // &
            trim(adjustl(line)) // " minutes may lead to the appearance of numerical instabilities in simulated states.")
        end if

    !> Call 'CLASSD' to initialize constants.
        !call CLASSD                                                                                                        !uncomment when fully understand if CLASSD provides unique variables - Daniel Mutton

            !> Read an initial value for geothermal flux from file.
        if (GGEOFLAG == 1) then
            iun = fls%fl(mfk%f18)%iun
            open(iun, file = trim(adjustl(fls%fl(mfk%f18)%fn)), status = 'old', action = 'read', iostat = ierr)
            read(iun, *) GGEOGRD(1)
            close(iun)
        else
            GGEOGRD(1) = 0.0
        end if

        !> Initialize PBSM or allocate and initialize variables used in CLASS even if PBSM is not enabled.
        if (pbsm%PROCESS_ACTIVE) then
            call PBSM_init(fls, shd)
        else

            !> Variables used in CLASSZ.
            allocate(pbsm%vs%Drift(ilg), pbsm%vs%Subl(ilg))
            pbsm%vs%Drift = 0.0; pbsm%vs%Subl = 0.0
        end if

        !> Input options and variables.
        TITLE1 = ''
        TITLE2 = ''
        TITLE3 = ''
        TITLE4 = ''
        TITLE5 = ''
        TITLE6 = ''
        NAME1 = ''
        NAME2 = ''
        NAME3 = ''
        NAME4 = ''
        NAME5 = ''
        NAME6 = ''
        PLACE1 = ''
        PLACE2 = ''
        PLACE3 = ''
        PLACE4 = ''
        PLACE5 = ''
        PLACE6 = ''
        c_switch%IDISP = NO_DATA_INT
        c_switch%IZREF = NO_DATA_INT
        c_switch%ISLFD = NO_DATA_INT
        c_switch%IPCP = NO_DATA_INT
        IWFGAT = NO_DATA_INT
        c_switch%ISNOALB = NO_DATA_INT
        c_switch%ITC = NO_DATA_INT
        c_switch%ITCG = NO_DATA_INT
        c_switch%ITG = NO_DATA_INT
        c_switch%IPAI = NO_DATA_INT
        c_switch%IHGT = NO_DATA_INT
        c_switch%IALC = NO_DATA_INT
        c_switch%IALS = NO_DATA_INT
        c_switch%IALG = NO_DATA_INT
!        ICTEMMOD = NO_DATA_INT
        SOILINIFLAG = NO_DATA_INT
        NUDGESDEPFLAG = NO_DATA_INT
        FROZENSOILINFILFLAG = NO_DATA_INT
        class_gat%ZRFMGAT = NO_DATA_REAL
        class_gat%ZRFHGAT = NO_DATA_REAL
        class_gat%ZBLDGAT = NO_DATA_REAL
        class_gat%GCGAT = NO_DATA_REAL
        class_gat%FCANGAT = NO_DATA_REAL
        class_gat%LNZ0GAT = NO_DATA_REAL
        class_gat%ALVCGAT = NO_DATA_REAL
        class_gat%ALICGAT = NO_DATA_REAL
        class_gat%PAMXGAT = NO_DATA_REAL
        class_gat%PAMNGAT = NO_DATA_REAL
        class_gat%CMASGAT = NO_DATA_REAL
        class_gat%ROOTGAT = NO_DATA_REAL
        class_gat%RSMNGAT = NO_DATA_REAL
        class_gat%QA50GAT = NO_DATA_REAL
        class_gat%VPDAGAT = NO_DATA_REAL
        class_gat%VPDBGAT = NO_DATA_REAL
        class_gat%PSGAGAT = NO_DATA_REAL
        class_gat%PSGBGAT = NO_DATA_REAL
        class_gat%DRNGAT = NO_DATA_REAL
        vgat%sdepgat = NO_DATA_REAL
        vgat%faregat = NO_DATA_REAL
        DDGAT = NO_DATA_REAL
        class_gat%XSLPGAT = NO_DATA_REAL
        XDGAT = NO_DATA_REAL
        MANNGAT = NO_DATA_REAL
        KSGAT = NO_DATA_REAL
        MIDGAT = NO_DATA_INT
        vgat%sandgat = NO_DATA_REAL
        vgat%claygat = NO_DATA_REAL
        vgat%orgmgat = NO_DATA_REAL
        class_gat%TBARGAT = NO_DATA_REAL
        class_gat%TBASGAT = NO_DATA_REAL
        class_gat%QACGAT = NO_DATA_REAL
        class_gat%TACGAT = NO_DATA_REAL
        class_gat%TCANGAT = NO_DATA_REAL
        class_gat%TSNOGAT = NO_DATA_REAL
        class_gat%TPNDGAT = NO_DATA_REAL
        class_gat%THLQGAT = NO_DATA_REAL
        class_gat%THICGAT = NO_DATA_REAL
        class_gat%TSFSGAT = NO_DATA_REAL
        class_gat%ZPNDGAT = NO_DATA_REAL
        class_gat%RCANGAT = NO_DATA_REAL
        class_gat%SCANGAT = NO_DATA_REAL
        class_gat%CMAIGAT = NO_DATA_REAL
        class_gat%SNOGAT = NO_DATA_REAL
        class_gat%ALBSGAT = NO_DATA_REAL
        class_gat%RHOSGAT = NO_DATA_REAL
        class_gat%WSNOGAT = NO_DATA_REAL
        class_gat%GROGAT = NO_DATA_REAL
        class_gat%ZPLGGAT = NO_DATA_REAL
        class_gat%ZPLSGAT = NO_DATA_REAL
        class_gat%ZSNLGAT = NO_DATA_REAL
        FREZTH = NO_DATA_REAL
        SWELIM = NO_DATA_REAL
        SNDENLIM = NO_DATA_REAL

        !> Allocate the filter for assigning distributed fields when MID (mosaic ID) is used.
        allocate(value_mask(ilg))
        value_mask = .true.

         !> Assign input variables from field list.
        do I = 1, size (input_field_list)

            !> Reset adjustment factors.
            factor_m = 1.0
            factor_a = 0.0
            valid_max = huge(valid_max)
            valid_min = tiny(valid_min)

            !> Assign recognized fields.
            associate(field => input_field_list(i))
                select case (field%field_name)
                    case ('TITLE', 'CLASSOUTTITLE')
                        call get_field_value(field, CLASSOUTTITLE)
                        TITLE1 = CLASSOUTTITLE(1:4)
                        TITLE2 = CLASSOUTTITLE(5:8)
                        TITLE3 = CLASSOUTTITLE(9:12)
                        TITLE4 = CLASSOUTTITLE(13:16)
                        TITLE5 = CLASSOUTTITLE(17:20)
                        TITLE6 = CLASSOUTTITLE(21:24)
                    case ('NAME', 'CLASSOUTNAME')
                        call get_field_value(field, CLASSOUTNAME)
                        NAME1 = CLASSOUTNAME(1:4)
                        NAME2 = CLASSOUTNAME(5:8)
                        NAME3 = CLASSOUTNAME(9:12)
                        NAME4 = CLASSOUTNAME(13:16)
                        NAME5 = CLASSOUTNAME(17:20)
                        NAME6 = CLASSOUTNAME(21:24)
                    case ('PLACE', 'CLASSOUTPLACE')
                        call get_field_value(field, CLASSOUTPLACE)
                        PLACE1 = CLASSOUTPLACE(1:4)
                        PLACE2 = CLASSOUTPLACE(5:8)
                        PLACE3 = CLASSOUTPLACE(9:12)
                        PLACE4 = CLASSOUTPLACE(13:16)
                        PLACE5 = CLASSOUTPLACE(17:20)
                        PLACE6 = CLASSOUTPLACE(21:24)
                    case ('IDISP')
                        call get_field_value(field, c_switch%IDISP)
                    case ('IZREF')
                        call get_field_value(field, c_switch%IZREF)
                    case ('ISLFD')
                        call get_field_value(field, c_switch%ISLFD)
                    case ('IPCP')
                        call get_field_value(field, c_switch%IPCP)
                    case ('ITC')
                        call get_field_value(field, c_switch%ITC)
                    case ('ITCG')
                        call get_field_value(field, c_switch%ITCG)
                    case ('ITG')
                        call get_field_value(field, c_switch%ITG)
                    case ('IWF')
                        call get_field_value(field, IWFGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('IPAI')
                        call get_field_value(field, c_switch%IPAI)
                    case ('IHGT')
                        call get_field_value(field, c_switch%IHGT)
                    case ('IALC')
                        call get_field_value(field, c_switch%IALC)
                    case ('IALS')
                        call get_field_value(field, c_switch%IALS)
                    case ('IALG')
                        call get_field_value(field, c_switch%IALG)
                    case ('ISNOALB')
                        call get_field_value(field, c_switch%ISNOALB)
!                    case ('ICTEMMOD')
!                        call get_field_value(field, ICTEMMOD)
                    case ('SOILINIFLAG')
                        call get_field_value(field, SOILINIFLAG)
                    case ('NUDGESDEPFLAG')
                        call get_field_value(field, NUDGESDEPFLAG)
                    case ('FROZENSOILINFILFLAG')
                        call get_field_value(field, FROZENSOILINFILFLAG)
                    case ('ZRFM', 'ZREFMOMENTUM')
                        call get_field_value(field, class_gat%ZRFMGAT, DIM_NAME_G, valid_min = 0.0, value_mask = value_mask)
                    case ('ZRFH', 'ZREFTHERMODYN')
                        call get_field_value(field, class_gat%ZRFHGAT, DIM_NAME_G, valid_min = 0.0, value_mask = value_mask)
                    case ('ZBLD', 'ZREFBLEND')
                        call get_field_value(field, class_gat%ZBLDGAT, DIM_NAME_G, valid_min = 0.0, value_mask = value_mask)
                    case ('GRIDCOVER')
                        call get_field_value( &
                            field, class_gat%GCGAT, DIM_NAME_G, valid_min = -1.0, valid_max = -1.0, value_mask = value_mask)
                    case ('FCAN')
                        call get_field_value( &
                            field, class_gat%FCANGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, &
                            value_mask = value_mask)
                    case ('LNZ0')
                        call get_field_value( &
                            field, class_gat%LNZ0GAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, &
                            value_mask = value_mask)
                    case ('ALVC')
                        call get_field_value( &
                            field, class_gat%ALVCGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, &
                            value_mask = value_mask)
                    case ('ALIC')
                        call get_field_value( &
                            field, class_gat%ALICGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, &
                            value_mask = value_mask)
                    case ('PAMX', 'LAMX')
                        call get_field_value( &
                            field, class_gat%PAMXGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP(1:4), &
                            value_mask = value_mask)
                    case ('PAMN', 'LAMN')
                        call get_field_value( &
                            field, class_gat%PAMNGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP(1:4), &
                            value_mask = value_mask)
                    case ('CMAS')
                        call get_field_value( &
                            field, class_gat%CMASGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP(1:4), &
                            value_mask = value_mask)
                    case ('ROOT')
                        call get_field_value( &
                            field, class_gat%ROOTGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP(1:4), &
                            value_mask = value_mask)
                    case ('RSMN')
                        call get_field_value( &
                            field, class_gat%RSMNGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP(1:4), &
                            value_mask = value_mask)
                    case ('QA50')
                        call get_field_value( &
                            field, class_gat%QA50GAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP(1:4), &
                            value_mask = value_mask)
                    case ('VPDA')
                        call get_field_value( &
                            field, class_gat%VPDAGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP(1:4), &
                            value_mask = value_mask)
                    case ('VPDB')
                        call get_field_value( &
                            field, class_gat%VPDBGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP(1:4), &
                            value_mask = value_mask)
                    case ('PSGA')
                        call get_field_value( &
                            field, class_gat%PSGAGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP(1:4), &
                            value_mask = value_mask)
                    case ('PSGB')
                        call get_field_value( &
                            field, class_gat%PSGBGAT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP(1:4), &
                            value_mask = value_mask)
                    case ('DRN')
                        call get_field_value( &
                            field, class_gat%DRNGAT, DIM_NAME_G, valid_min = 0.0, valid_max = 1.0, value_mask = value_mask)
                    case ('SDEP', 'SDEPTH')
                        call get_field_value(field, vgat%sdepgat, DIM_NAME_G, valid_min = 0.0, value_mask = value_mask)
                    case ('FARE')
                        call get_field_value_real1d( &
                            field, vgat%faregat, DIM_NAME_G, valid_min = 1.0, valid_max = 1.0, value_mask = value_mask)
                    case ('DD', 'DDEN')
                        if (index(to_lowercase(field%meta%units), 'km') > 0) then
                            !> Unit conversion to m m**-2 if units are km km**-2 (WATROF/WATDRN expects m m**-2).
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'm m**-2'.")
                            factor_m = 1.0/1000.0
                        end if
                    call get_field_value(field, DDGAT, DIM_NAME_G, const_mul = factor_m, value_mask = value_mask)
                    case ('XSLP')
                        call get_field_value(field, class_gat%XSLPGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('GRKEF', 'GRKF', 'XDRAINH')
                        call get_field_value(field, XDGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('MANN')
                        call get_field_value(field, MANNGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('KSAT')
                        call get_field_value(field, KSGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('MOSAICID', 'MID')
                        call get_field_value(field, MIDGAT, DIM_NAME_G, valid_min = 1)

                        !> Check for mosaic ID to control assigning distributed parameters.
                        if (any(MIDGAT >= 100 .and. MIDGAT < 1000)) then
                            call print_remark( &
                                "There are tile values of variable '" // trim(field%label) // "' > 100 and < 1000. " // &
                                "Those tiles will not be reassigned new values if the parameter inputs include " // &
                                "spatially distributed fields.")
                            where (MIDGAT >= 100 .and. MIDGAT < 1000) value_mask = .false.
                        end if
                    case ('SOLSAND', 'SAND')
                        call get_field_value(field, vgat%sandgat, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('SOLCLAY', 'CLAY')
                        call get_field_value(field, vgat%claygat, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('SOLORGCONTENT', 'ORGM')
                        call get_field_value(field, vgat%orgmgat, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('TSOL', 'TBAR', 'TILE_TBAR')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value( &
                            field, class_gat%TBARGAT, (/DIM_NAME_G, DIM_NAME_L/), const_add = factor_a, value_mask = value_mask)
                    case ('TBAS', 'TILE_TBAS')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value(field, class_gat%TBASGAT, DIM_NAME_G, const_add = factor_a, value_mask = value_mask)
                    case ('QAC', 'TILE_QAC')
                        call get_field_value(field, class_gat%QACGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('TAC', 'TILE_TAC')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value(field, class_gat%TACGAT, DIM_NAME_G, const_add = factor_a, value_mask = value_mask)
                    case ('TCAN', 'TILE_TCAN')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value(field, class_gat%TCANGAT, DIM_NAME_G, const_add = factor_a, value_mask = value_mask)
                    case ('TSNO', 'TILE_TSNO')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value(field, class_gat%TSNOGAT, DIM_NAME_G, const_add = factor_a, value_mask = value_mask)
                    case ('TPND', 'TILE_TPND')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value(field, class_gat%TPNDGAT, DIM_NAME_G, const_add = factor_a, value_mask = value_mask)
                    case ('THLQSOL', 'THLQ', 'TILE_THLQ')
                        call get_field_value(field, class_gat%THLQGAT, (/DIM_NAME_G, DIM_NAME_L/), valid_min = 0.04, value_mask = value_mask)
                    case ('THICSOL', 'THIC', 'TILE_THIC')
                        call get_field_value(field, class_gat%THICGAT, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('TSFS', 'TILE_TSFS')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value( &
                            field, class_gat%TSFSGAT, (/DIM_NAME_G, DIM_NAME_L/), const_add = factor_a, value_mask = value_mask)
                    case ('ZPND', 'TILE_ZPND')
                        call get_field_value(field, class_gat%ZPNDGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('LQWSCAN', 'RCAN', 'TILE_RCAN')
                        call get_field_value(field, class_gat%RCANGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('FZWSCAN', 'SCAN', 'TILE_SNCAN')
                        call get_field_value(field, class_gat%SCANGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('CMAI', 'TILE_CMAS')
                        call get_field_value(field, class_gat%CMAIGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('SNO', 'TILE_SNO')
                        call get_field_value(field, class_gat%SNOGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('ALBSNO', 'ALBS', 'TILE_ALBS')
                        call get_field_value(field, class_gat%ALBSGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('RHOSNO', 'RHOS', 'TILE_RHOS')
                        call get_field_value(field, class_gat%RHOSGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('LQWSSNO', 'WSNO', 'TILE_WSNO')
                        call get_field_value(field, class_gat%WSNOGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('GRO', 'TILE_GRO')
                        call get_field_value(field, class_gat%GROGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('ZPLG')
                        call get_field_value(field, class_gat%ZPLGGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('ZPLS')
                        call get_field_value(field, class_gat%ZPLSGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('ZSNL')
                        call get_field_value(field, class_gat%ZSNLGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('FREZTH', 'ICEBAL_FREEZE_THRESHOLD')
                        call get_field_value(field, FREZTH, DIM_NAME_G)
                    case ('SWELIM', 'ICEBAL_SWE_LIMIT', 'SNDEPLIM', 'ICEBAL_SNOW_DEPTH_LIMIT')
                        call get_field_value(field, SWELIM, DIM_NAME_G)
                    case ('SNDENLIM', 'ICEBAL_SNOW_DENSITY_LIMIT')
                        call get_field_value(field, SNDENLIM, DIM_NAME_G)
                    case ('SOCI')
                        call get_field_value(field, vgat%socigat, DIM_NAME_G, value_mask = value_mask)

                    !> CTEM parameters.
                    ! ic(ic) = 1, 2, 3, 4 ; "physics (CLASS) PFTs"
                    ! iccp2(iccp2) = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 ; "biogeochemical (CTEM) PFTs plus bareground and land use change product pools"
                    ! icctem(icctem) = 1, 2, 3, 4, 5, 6, 7, 8, 9 ; "biogeochemical (CTEM) PFTs"
                    ! icp1(icp1) = 1, 2, 3, 4, 5 ; "physics (CLASS) PFTs plus bareground"
                    case ('IPEATLAND') !< "Peatland flag: 0 = not a peatland, 1= bog, 2 = fen"
                        call get_field_value(field, vgat%ipeatland, DIM_NAME_G, value_mask = value_mask)
                    case ('BLEAFMAS') !< "Brown leaf mass"
                        call get_field_value(field, vgat%bleafmas, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('DMOSS') !< "Depth of living moss"
                        call get_field_value(field, vgat%dmoss, DIM_NAME_G, value_mask = value_mask)
                    case ('FCANCMX') !< "PFT fractional coverage per grid cell"
                        call get_field_value(field, vgat%fcancmx, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('GLEAFMAS') !< "Green leaf mass"
                        call get_field_value(field, vgat%gleafmas, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('LFSTATUS') !< "Leaf status, see Phenology"
                        call get_field_value(field, vgat%lfstatus, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('LITRMASS') !< "Litter mass per soil layer"
                        call get_field_value(field, vgat%litrmass, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('LITRMSMOSS') !< "Moss litter mass"
                        call get_field_value(field, vgat%litrmsmoss, DIM_NAME_G, value_mask = value_mask)
                    case ('PANDAYS') !< "Days with +ve new photosynthesis, see Phenology"
                        call get_field_value(field, vgat%pandays, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('ROOTMASS') !< "Root mass"
                        call get_field_value(field, vgat%rootmass, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    ! slope(slope) = 1, 2, 3, 4, 5, 6, 7, 8 ; "wetland slope fractions for 0.025, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3 and 0.35 percent slope threshold"
                    case ('SLOPEFRAC') !< "Slope-based fraction for dynamic wetlands"
                        call get_field_value(field, vgat%slopefrac, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('SOILCMAS') !< "Soil C mass per soil layer"
                        call get_field_value(field, vgat%soilcmas, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('STEMMASS') !< "Stem mass"
                        call get_field_value(field, vgat%stemmass, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('MAXANNUALACTLYR') !< "Active layer thickness maximum over the e-folding period specified by parameter eftime"
                        call get_field_value(field, class_gat%maxannualactlyrgat, DIM_NAME_G, value_mask = value_mask)
                    case ('GRWTHEFF') !< "Growth efficiency. Change in biomass per year per unit max. LAI,for use in mortality subroutine"
                        call get_field_value(field, vgat%grwtheff, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('CO2', 'CO2CONC')
                        call get_field_value(field, vgat%co2conc, DIM_NAME_G, value_mask = value_mask)
                end select
            end associate
        end do
		

        !> Assign defaults or derived fields where possible.
        ierr = 0
        z = 0
		class_gat%actLyrThisYrGAT = 0

        !> Options.
        call check_field_assigned(c_switch%IDISP, 'IDISP', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%IZREF, 'IZREF', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%ISLFD, 'ISLFD', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%IPCP, 'IPCP', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(IWFGAT, 'IWF', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%ISNOALB, 'ISNOALB', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%ITC, 'ITC', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%ITCG, 'ITCG', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%ITG, 'ITG', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%IPAI, 'IPAI', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%IHGT, 'IHGT', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%IALC, 'IALC', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%IALS, 'IALS', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(c_switch%IALG, 'IALG', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(SOILINIFLAG, 'SOILINIFLAG', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(NUDGESDEPFLAG, 'NUDGESDEPFLAG', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(FROZENSOILINFILFLAG, 'FROZENSOILINFILFLAG', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%ZRFHGAT, 'ZRFH', 28.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%ZRFMGAT, 'ZRFM', 28.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%ZBLDGAT, 'ZBLD', 50.0, error_status = z); if (z /= 0) ierr = z
        if (z /= 0) then
            ierr = z
        else if (NUDGESDEPFLAG == 1) then

            !> Nudge SDEP to the nearest interface between soil layers in the profile.
            call print_remark( &
                "'NUDGESDEPFLAG' is enabled. The permeable depth 'SDEP' will be nudged at each location to the nearest " // &
                "interface between soil layers in the soil profile.")
            where (vgat%sdepgat <= shd%lc%sl%zbot(1))

                !> Within the first layer, nudge to the first boundary and not to zero so it remains active.
                vgat%sdepgat = shd%lc%sl%zbot(1)
            elsewhere (vgat%sdepgat > shd%lc%sl%zbot(ignd))

                !> Limit depths deeper than the soil profile to the last boundary for consistency.
                vgat%sdepgat = shd%lc%sl%zbot(ignd)
            end where
            do j = 2, ignd

                !> Nudge in-between depths to the closest boundary.
                where (vgat%sdepgat > shd%lc%sl%zbot(j - 1) .and. vgat%sdepgat <= (0.5*shd%lc%sl%delz(j) + shd%lc%sl%zbot(j - 1)))
                    vgat%sdepgat = shd%lc%sl%zbot(j - 1)
                elsewhere (pm%tile%sdep > (0.5*shd%lc%sl%delz(j) + shd%lc%sl%zbot(j - 1)) .and. vgat%sdepgat <= shd%lc%sl%zbot(j))
                    vgat%sdepgat = shd%lc%sl%zbot(j)
                end where
            end do
        end if
        call check_field_assigned(vgat%faregat, 'FARE', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(DDGAT, 'DD', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%XSLPGAT, 'XSLP', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(XDGAT, 'GRKEF', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(MANNGAT, 'MANN', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(KSGAT, 'KSAT', is_error = .true., error_status = z); if (z /= 0) ierr = z
        if (any(class_gat%DRNGAT < 1.0 .and. vgat%sdepgat <= 4.1)) then
            call print_remark( &
                "The configuration contains drainage index 'DRN' values less than 1.0. " // &
                "In most cases, 'DRN' can be set to 1.0 except in cases of deep soils where it is desired to " // &
                "suppress drainage from the bottom of the soil profile or in soils with a high water table.")
        end if
        call check_field_assigned(MIDGAT, 'MOSAICID', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(vgat%sandgat, 'SAND', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else
            if (any(vgat%sdepgat == 0.0 .and. vgat%sandgat(:, 1) > -2.5)) then

                !> Check the first layer for impermable soils not identified as such.
                write(line, *) count(vgat%sdepgat == 0.0 .and. vgat%sandgat(:, 1) > -2.5)
                call print_remark( &
                    "The configuration contains " // trim(adjustl(line)) // " tiles where the permeable depth 'SDEP' is zero. " // &
                    "The soil properties in those areas will be set to those of 'exposed rock' ('SAND' -3).")
                do j = 1, ignd
                    where (vgat%sdepgat == 0.0 .and. vgat%sandgat(:, 1) > -2.5) vgat%sandgat(:, j) = -3.0
                end do
            end if

            !> Check for unidentified impermeable soils at depth.
            do j = 2, ignd
                where (vgat%sdepgat < (shd%lc%sl%ZBOT(j - 1) + 0.001) .and. vgat%sandgat(:, j) > -2.5) vgat%sandgat(:, j) = -3.0
            end do
        end if
        call check_field_assigned(vgat%claygat, 'CLAY', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else if (any(vgat%sandgat <= -2)) then
            call print_remark( &
                "The configuration contains areas of non-mineral soils ('SAND' -2, -3 or -4). " // &
                "The 'CLAY' fraction will be set to zero in those areas.")
            where (vgat%sandgat <= -2) vgat%claygat = 0.0
        end if
        call check_field_assigned(vgat%orgmgat, 'ORGM', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else
            if (any(vgat%sandgat <= -2)) then
                call print_remark( &
                    "The configuration contains areas of non-mineral soils ('SAND' -2, -3 or -4). " // &
                    "The 'ORGM' fraction will be set to zero in those areas.")
                where (vgat%sandgat <= -2) vgat%orgmgat = 0.0
            end if
            if (any((vgat%orgmgat /= 1.0 .or. vgat%orgmgat /= 2.0 .or. vgat%orgmgat /= 3.0) .and. vgat%sandgat == -2)) then
                call print_error("The configuration contains areas of organic soils ('SAND' -2) but improper 'ORGM' class types.")
                ierr = 1
            end if
        end if
        call check_field_assigned(class_gat%ZPLGGAT, 'ZPLG', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%ZPLSGAT, 'ZPLS', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%ZSNLGAT, 'ZSNL', is_error = .true., error_status = z); if (z /= 0) ierr = z
        if (any(vgat%sandgat == -4)) then
            call check_field_assigned(FREZTH, 'FREZTH', -2.0, error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(SWELIM, 'SWELIM', 100.0, error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(SNDENLIM, 'SNDENLIM', 900.0, error_status = z); if (z /= 0) ierr = z
        else
            if (any(FREZTH /= NO_DATA_REAL)) then
                call print_remark( &
                    "The variable '" // 'FREZTH' // "' is assigned but not active (no glaciers are active).")
            end if
            if (any(SWELIM /= NO_DATA_REAL)) then
                call print_remark( &
                    "The variable '" // 'SWELIM' // "' is assigned but not active (no glaciers are active).")
            end if
            if (any(SNDENLIM /= NO_DATA_REAL)) then
                call print_remark( &
                    "The variable '" // 'SNDENLIM' // "' is assigned but not active (no glaciers are active).")
            end if
        end if

        !> States.
        call check_field_assigned(class_gat%TCANGAT, 'TCAN', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else if (any(class_gat%TCANGAT < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TCAN' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (class_gat%TCANGAT < (TFREZ - 100.0)) class_gat%TCANGAT = class_gat%TCANGAT + TFREZ
        end if
        call check_field_assigned(class_gat%TACGAT, 'TAC', error_status = z)
        if (z /= 0 .and. all(class_gat%TCANGAT /= NO_DATA_REAL)) then
!+ Set unassigned TAC to air temperature

            !> Set the field to 'TCAN'.
            call print_remark("Setting the unassigned values of 'TAC' from 'TCAN'.")
            where (class_gat%TACGAT == NO_DATA_REAL) class_gat%TACGAT = class_gat%TCANGAT
        else if (any(class_gat%TACGAT < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TAC' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (class_gat%TACGAT < (TFREZ - 100.0)) class_gat%TACGAT = class_gat%TACGAT + TFREZ
        end if
        call check_field_assigned(class_gat%QACGAT, 'QAC', 0.5e-2, error_status = z); if (z /= 0) ierr = z
!+ Set unassigned QAC to air specific humidity
        call check_field_assigned(class_gat%RCANGAT, 'RCAN', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%SCANGAT, 'SCAN', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%CMAIGAT, 'RCAN', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%CMAIGAT, 'CMAI', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%GROGAT, 'GRO', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%ZPNDGAT, 'ZPND', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(class_gat%TPNDGAT, 'TPND', 0.0, error_status = z)
        if (any(class_gat%TPNDGAT < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TPND' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (class_gat%TPNDGAT < (TFREZ - 100.0)) class_gat%TPNDGAT = class_gat%TPNDGAT + TFREZ
        end if
        call check_field_assigned(class_gat%SNOGAT, 'SNO', 0.0, error_status = z)
        if (any(class_gat%SNOGAT > 0.0)) then
            call print_remark("Starting the simulation in snow-covered conditions could lead to a bias in simulated states.")
            call check_field_assigned(class_gat%ALBSGAT, 'ALBS', is_error = .true., error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(class_gat%RHOSGAT, 'RHOS', is_error = .true., error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(class_gat%TSNOGAT, 'TSNO', is_error = .true., error_status = z)
            if (z /= 0) then
                ierr = z
            else if (any(class_gat%TSNOGAT < (TFREZ - 100.0))) then

                !> Implicit conversion.
                call print_remark( &
                    "There may be tile values of the variable '" // 'TSNO' // &
                    "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
                where (class_gat%TSNOGAT < (TFREZ - 100.0)) class_gat%TSNOGAT = class_gat%TSNOGAT + TFREZ
            end if
            call check_field_assigned(class_gat%WSNOGAT, 'WSNO', 0.0, error_status = z); if (z /= 0) ierr = z
        else
            class_gat%ALBSGAT = 0.0
            class_gat%RHOSGAT = 0.0
            class_gat%TSNOGAT = TFREZ
            class_gat%WSNOGAT = 0.0
        end if
        call check_field_assigned(class_gat%TBARGAT, 'TBAR', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else if (any(class_gat%TBARGAT < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TBAR' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (class_gat%TBARGAT < (TFREZ - 100.0)) class_gat%TBARGAT = class_gat%TBARGAT + TFREZ
        end if
        call check_field_assigned(class_gat%THICGAT, 'THIC', 0.0, error_status = z)
        if (z == 0) then
            if (any(vgat%sandgat <= -3)) then
                call print_remark( &
                    "Setting the volumetric frozen water content 'THIC' of non-soil layers ('SAND' -3 or -4) to zero.")
                where (vgat%sandgat <= -3) class_gat%THICGAT = 0.0
            end if
        end if
        call check_field_assigned(class_gat%THLQGAT, 'THLQ', 0.04, error_status = z)
        if (z == 0) then
            if (any(vgat%sandgat <= -3)) then
                call print_remark( &
                    "Setting the volumetric liquid water content 'THLQ' of non-soil layers ('SAND' -3 or -4) to zero.")
                where (vgat%sandgat <= -3) class_gat%THLQGAT = 0.0
            end if
            if (any(vgat%sandgat >= -2 .and. class_gat%THLQGAT < 0.04)) then
                write(line, *) count(vgat%sandgat >= -2 .and. class_gat%THLQGAT < 0.04)
                call print_warning( &
                    "The configuration contains " // trim(adjustl(line)) // " tiles where the initial value for the " // &
                    "volumetric liquid water content 'THLQ' is lower than the permissible limit of 0.04.")
            end if
        end if
        call check_field_assigned(class_gat%TBASGAT, 'TBAS', error_status = z)
        if (z /= 0 .and. all(class_gat%TBARGAT /= NO_DATA_REAL)) then
            call print_remark("Assigning the unassigned values of 'TBAS' from the bottom layer value of 'TBAR'.")
            where (class_gat%TBASGAT == NO_DATA_REAL) class_gat%TBASGAT = class_gat%TBARGAT(:, ignd)
        else if (any(class_gat%TBASGAT < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TBAS' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (class_gat%TBASGAT < (TFREZ - 100.0)) class_gat%TBASGAT = class_gat%TBASGAT + TFREZ
        end if
        call check_field_assigned(class_gat%TSFSGAT, 'TSFS', error_status = z)
        if (z /= 0 .and. all(class_gat%TBARGAT /= NO_DATA_REAL)) then
            call print_remark("Assigning the unassigned values of 'TSFS' to a background field and from 'TBAR'.")
            do j = 1, 2

                !> Assign snow-covers the freezing point of water.
                where (class_gat%TSFSGAT(:, j) == NO_DATA_REAL) class_gat%TSFSGAT(:, j) = TFREZ
            end do
            do j = 3, 4

                !> Assign ground covers to the first layer 'TBAR'.
                where (class_gat%TSFSGAT(:, j) == NO_DATA_REAL) class_gat%TSFSGAT(:, j) = class_gat%TBARGAT(:, 1)
            end do
        else if (any(class_gat%TSFSGAT < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TSFS' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (class_gat%TSFSGAT < (TFREZ - 100.0)) class_gat%TSFSGAT = class_gat%TSFSGAT + TFREZ
        end if

        !> Old variables. Find a better way to connect - Daniel Mutton
        !vs%tile%cmas(il1:il2) = class_gat%CMAIGAT(il1:il2)
        !vs%tile%lqwssno(il1:il2) = class_gat%WSNOGAT(il1:il2)
        !vs%tile%qacan(il1:il2) = class_gat%QACGAT(il1:il2)
        !vs%tile%tcan(il1:il2) = class_gat%TCANGAT(il1:il2)
        !vs%tile%tacan(il1:il2) = class_gat%TACGAT(il1:il2)
        !vs%tile%tsno(il1:il2) = class_gat%TSNOGAT(il1:il2)
        !vs%tile%tpnd(il1:il2) = class_gat%TPNDGAT(il1:il2)
        !vs%tile%zpnd(il1:il2) = class_gat%ZPNDGAT(il1:il2)
        !vs%tile%lqwscan(il1:il2) = class_gat%RCANGAT(il1:il2)
        !vs%tile%fzwscan(il1:il2) = class_gat%SCANGAT(il1:il2)
        !vs%tile%sno(il1:il2) = class_gat%SNOGAT(il1:il2)
        !vs%tile%albsno(il1:il2) = class_gat%ALBSGAT(il1:il2)
        !vs%tile%rhosno(il1:il2) = class_gat%RHOSGAT(il1:il2)
        !vs%tile%gro(il1:il2) = class_gat%GROGAT(il1:il2)
        !vs%tile%tsfs(il1:il2, :) = class_gat%TSFSGAT(il1:il2, :)
        !vs%tile%tsol(il1:il2, :) = class_gat%TBARGAT(il1:il2, :)
        !vs%tile%thlqsol(il1:il2, :) = class_gat%THLQGAT(il1:il2, :)
        !vs%tile%thicsol(il1:il2, :) = class_gat%THICGAT(il1:il2, :)
        !vs%tile%tbas(il1:il2) = class_gat%TBASGAT(il1:il2)

        !> Derived or other fields.
        class_gat%ITCTGAT = 0

        !> FROZENSOILINFILFLAG 1.
!        allocate(FRZCGAT(NML))
        allocate(INFILTYPE(ilg), SI(ilg), TSI(ilg), &
                 SNOWMELTD(ilg), SNOWMELTD_LAST(ilg), SNOWINFIL(ilg), &
                 CUMSNOWINFILCS(ilg), MELTRUNOFF(ilg), CUMSNOWINFILGS(ilg))

        !todo: move to read_parameters
        if (FROZENSOILINFILFLAG == 1) then
            NMELT = 1
            INFILTYPE = 2 !> INITIALIZED WITH UNLIMITED INFILTRATION
            SNOWMELTD = 0.0
            SNOWINFIL = 0.0
            CUMSNOWINFILCS = 0.0
            CUMSNOWINFILGS = 0.0
            MELTRUNOFF = 0.0
            SI = 0.20
            TSI = -0.10
            do k = il1, il2
                FRZCGAT(k) = FRZCROW(class_gat%ILMOS(k), class_gat%JLMOS(k))
            end do
        end if

        !> *********************************************************************
        !> Call CLASSBG to set more CLASS variables
        !> *********************************************************************

        !> CLASS output files.
        !if (WF_NUM_POINTS > 0) call CLASSOUT_open_files(shd)                                                                               !commented out - Daniel Mutton

        !> Check if coordinates should be assigned from parameter inputs (e.g., point mode).
        if (all(vs%grid%lat == NO_DATA_REAL) .or. all(vs%grid%lon == NO_DATA_REAL)) then
            call print_remark("Reading coordination location information from parameter input files.")
            call get_field_value(input_field_list, DIM_NAMES_OF_Y, vs%grid%lat)
            call get_field_value(input_field_list, DIM_NAMES_OF_X, vs%grid%lon)
        end if
        call check_field_assigned(vs%grid%lat, 'DEGLAT', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else

            !> Convert latitudes to radians.
            do i = 1, vs%tile%dim_length
                class_gat%DLATGAT(i) = vs%grid%lat(vs%tile%from_cell(i))
                class_gat%RADJGAT(i) = vs%grid%lat(vs%tile%from_cell(i))*PI/180.0
            end do
        end if
        call check_field_assigned(vs%grid%lon, 'DEGLON', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else

            !> Convert longitudes from (-180:180) to (0:360).
!+            where (vs%grid%lon < 0.0) vs%grid%lon = vs%grid%lon + 360.0
            where (vs%grid%lon > 180.0) vs%grid%lon = vs%grid%lon - 360.0
            do i = 1, vs%tile%dim_length
                class_gat%DLONGAT(i) = vs%grid%lon(vs%tile%from_cell(i))
            end do
        end if

        class_gat%Z0ORGAT = 0.0
        class_gat%GGEOGAT(:) = GGEOGRD(1)
        class_gat%ZDMGAT = 10.0
        class_gat%ZDHGAT = 2.0

        !> Check for errors.
        if (ierr /= 0) then
            call reset_tab()
            call print_error("Variables required by the module are not assigned or were missing during initialization.")
            call program_abort()
        end if

        !> IWF 2 (PDMROF) and IWF 3 (LATFLOW).
        allocate(&
            CMINPDM(ilg), CMAXPDM(ilg), BPDM(ilg), K1PDM(ilg), K2PDM(ilg), &
            ZPNDPRECS(ilg), ZPONDPREC(ilg), ZPONDPREG(ilg), ZPNDPREGS(ilg), &
            UM1CS(ilg), UM1C(ilg), UM1G(ilg), UM1GS(ilg), &
            QM1CS(ilg), QM1C(ilg), QM1G(ilg), QM1GS(ilg), &
            QM2CS(ilg), QM2C(ilg), QM2G(ilg), QM2GS(ilg), &
            UMQ(ilg), &
            FSTRCS(ilg), FSTRC(ilg), FSTRG(ilg), FSTRGS(ilg))
        ZPNDPRECS = 0.0; ZPONDPREC = 0.0; ZPONDPREG = 0.0; ZPNDPREGS = 0.0
        UM1CS = 0.0; UM1C = 0.0; UM1G = 0.0; UM1GS = 0.0
        QM1CS = 0.0; QM1C = 0.0; QM1G = 0.0; QM1GS = 0.0
        QM2CS = 0.0; QM2C = 0.0; QM2G = 0.0; QM2GS = 0.0
        UMQ = 0.0
        FSTRCS = 0.0; FSTRC = 0.0; FSTRG = 0.0; FSTRGS = 0.0
        if (any(IWFGAT == 2) .or. any(IWFGAT == 3)) then
            call print_new_section("Prairie pothole ponding runoff (PDMROF) is active.")
            CMINPDM = NO_DATA_REAL
            CMAXPDM = NO_DATA_REAL
            BPDM = NO_DATA_REAL
            K1PDM = NO_DATA_REAL
            K2PDM = NO_DATA_REAL

            !> Assign input variables from field list.
            do i = 1, size(input_field_list)
                associate(field => input_field_list(i))
                    select case (field%field_name)
                        case ('CMIN')
                            call get_field_value(field, CMINPDM, DIM_NAME_G)
                        case ('CMAX')
                            call get_field_value(field, CMAXPDM, DIM_NAME_G)
                        case ('B')
                            call get_field_value(field, BPDM, DIM_NAME_G)
                        case ('K1')
                            call get_field_value(field, K1PDM, DIM_NAME_G)
                        case ('K2')
                            call get_field_value(field, K2PDM, DIM_NAME_G)
                    end select
                end associate
            end do

            !> Assign defaults or derived fields where possible.
            ierr = 0
            z = 0
            call check_field_assigned(CMINPDM, 'CMIN', 0.0, error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(CMAXPDM, 'CMAX', is_error = .true., error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(BPDM, 'B', is_error = .true., error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(K1PDM, 'K1', 0.0, error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(K2PDM, 'K2', 0.0, error_status = z); if (z /= 0) ierr = z

            !> Check for errors.
            if (ierr /= 0) then
                call reset_tab()
                call print_error("Variables required by the module are not assigned or were missing during initialization.")
                call program_abort()
            end if
        end if

    !END OF RUNCLASS CONFIG, RETURNING TO CLASSIC

    if (c_switch%ctem_on) then
      !> Read in the inputs for a run with biogeochemical component turned on
      !> Regardless of whether lnduseon or not, we need to check the land cover that was read in
      !! and assign the CLASS PFTs as they are not read in when ctem_on.
      call initializeLandCoverGAT(ilg, il1, il2)
    end if

    !> In preparation for the use of the random number generator by disaggMet,
    !! we need to provide a seed to allow repeatable results.
    call initRandomSeed()
        !> As the last step in the initialization sequence, the subroutine soilProperties is
        !! called, to assign soil thermal and hydraulic properties on the basis of the
        !! textural information read in for each of the soil layers.

        class_gat%DELZ = shd%lc%sl%delz
        class_gat%ZBOT = shd%lc%sl%zbot
        call soilPropertiesGAT(class_gat%THPGAT, class_gat%THRGAT, class_gat%THMGAT, class_gat%BIGAT, class_gat%PSISGAT, class_gat%GRKSGAT, & ! Formerly CLASSB
                               class_gat%THRAGAT, class_gat%HCPSGAT, class_gat%TCSGAT, class_gat%THFCGAT, class_gat%THLWGAT, class_gat%PSIWGAT, &
                               class_gat%DLZWGAT, class_gat%ZBTWGAT, &
                               class_gat%ALGWVGAT, class_gat%ALGWNGAT, class_gat%ALGDVGAT, class_gat%ALGDNGAT, &
                               vgat%sandgat, vgat%claygat, vgat%orgmgat, vgat%SOCIGAT, class_gat%DELZ, class_gat%ZBOT, &
                               vgat%sdepgat, class_gat%ISNDGAT, class_gat%IGDRGAT, &
                               vs%landtile_count, il1, il2, IGND, vgat%ipeatland)

        pm%tile%thpor(il1:il2, :) = class_gat%THPGAT(il1:il2, :)
        pm%tile%thlret(il1:il2, :) = class_gat%THRGAT(il1:il2, :)
        pm%tile%thlmin(il1:il2, :) = class_gat%THMGAT(il1:il2, :)
        pm%tile%bi(il1:il2, :) = class_gat%BIGAT(il1:il2, :)
        pm%tile%psisat(il1:il2, :) = class_gat%PSISGAT(il1:il2, :)
        pm%tile%grksat(il1:il2, :) = class_gat%GRKSGAT(il1:il2, :)
        pm%tile%thlrat(il1:il2, :) = class_gat%THRAGAT(il1:il2, :)
        pm%tile%hcps(il1:il2, :) = class_gat%HCPSGAT(il1:il2, :)
        pm%tile%tcs(il1:il2, :) = class_gat%TCSGAT(il1:il2, :)
        pm%tile%thfc(il1:il2, :) = class_gat%THFCGAT(il1:il2, :)
        pm%tile%psiwlt(il1:il2, :) = class_gat%PSIWGAT(il1:il2, :)
        vs%tile%dzsolhyd(il1:il2, :) = class_gat%DLZWGAT(il1:il2, :)
        vs%tile%zsolhyd(il1:il2) = sum(class_gat%ZBTWGAT(il1:il2, :), 2)

        ! ctem initializations.
        if (c_switch%ctem_on) then
            call ctemInitGAT (il1, il2)
            call allometry( &
                vgat%gleafmas, vgat%bleafmas, vgat%stemmass, vgat%rootmass, & ! In
                1, ilg, ilg, class_gat%ZBTWGAT, & ! In
                vgat%sdepgat, vgat%fcancmx, & ! In
                vgat%ipeatland, class_gat%maxAnnualActLyrGAT, & ! In
                vgat%ailcg, vgat%ailcb, vgat%ailc, vgat%zolnc, & ! Out
                vgat%rmatc, vgat%rmatctem, vgat%slai, vgat%bmasveg, & ! Out
                vgat%cmasvegc, vgat%veghght, vgat%rootdpth, vgat%alvsctm, & ! Out
                vgat%alirctm, vgat%paic, vgat%slaic) ! Out
        end if   ! if (ctem_on)
		
        call add_available_output_field((/'npp'/), vgat%npp, (/DIM_NAME_G/))

    end subroutine

    subroutine CLASSICv1_0_tile(shd, fls)

        !> Variables and related functions.
        use classicParams
        use classStateVars, only: class_gat, initDiagnosticVars
        use ctemStateVars, only: c_switch, vgat, vrot, ctem_tile, tracer, resetMosaicAccum

        !> Utilities.
        use ctemUtilitiesGAT, only: accumulateForCTEM, dayEndCTEMPreparation
        use generalUtilsGAT, only: findLeapYears, findCloudinessGAT, findPermafrostVarsGAT, findDaylength

        !> Process modules.
        use ctemDriver, only: ctem
        use tracerModule, only: decay14C
        use fourBandAlbedo, only: fourBandDriver

        use PBSM_module
        use mpi_module
        use model_dates
        use FLAGS
        use model_files_variables      !For 'fls%GENDIR_OUT'.\
        use model_dates                !For current dates and counter
        use sa_mesh_common
        use model_files_variables
        use shd_variables

        implicit none

        double precision  :: time
        integer           :: ierr, o, cell, blocks, remainder
        real SDOWN1, FSDOWN2, FSDOWN3, RDAY, DECL, HOUR, COSZ, &
        CUMSNOWINFILCS(ilg), CUMSNOWINFILGS(ilg), INFILTYPE(ilg), MELTRUNOFF(ilg), &
        SNOWINFIL(ilg), longitude, latitude

        !> Local variables.
        integer k, ik, j, i, N, INGD1, NMELT

        ! leap year flag (if the switch 'leap' is true, this will be used, otherwise it remains false)
        logical :: leapnow

        type(ShedGridParams) :: shd
        type(fl_ids) :: fls
        type(GridParams) :: cs
        integer, pointer :: readMetStartYear    !< First year of meteorological forcing to read in from the met file

		!DM
		do i = il1,il2
		  vgat%dayl(i) = findDaylength(real(ic%now%jday), class_gat%radjgat(i))
		end do
		!print *,  "Day Length", vgat%dayl
		!print *,  "Day Max", vgat%dayl_max

        !> N is only used for debugging purposes.
        N = ic%ts_count

            !> Transfer driving variables.
            class_gat%FSVHGAT(il1:il2) = vs%tile%fsin(il1:il2)/2.0
            class_gat%FSIHGAT(il1:il2) = vs%tile%fsin(il1:il2)/2.0
            class_gat%FSSGAT(il1:il2) = vs%tile%fsin(il1:il2)
            class_gat%FDLGAT(il1:il2) = vs%tile%flin(il1:il2)
            if (c_switch%IPCP == 4) then
                class_gat%RPCPGAT(il1:il2) = vs%tile%prern(il1:il2)
                class_gat%SPCPGAT(il1:il2) = vs%tile%presno(il1:il2)

                !> Required to activate 'CLASSI'.
                class_gat%PREGAT(il1:il2) = class_gat%RPCPGAT(il1:il2) + class_gat%SPCPGAT(il1:il2)
            else
                class_gat%PREGAT(il1:il2) = vs%tile%pre(il1:il2)
            end if
            class_gat%TAGAT(il1:il2) = vs%tile%ta(il1:il2)
            class_gat%ULGAT(il1:il2) = vs%tile%uv(il1:il2)
            class_gat%VLGAT(il1:il2) = 0.0 !DOUBLE CHECK THIS IS CORRECT, DANIEL MUTTON
            class_gat%PRESGAT(il1:il2) = vs%tile%pres(il1:il2)
            class_gat%QAGAT(il1:il2) = max(vs%tile%qa(il1:il2), 1.0E-6)
            class_gat%VMODGAT(il1:il2) = max(VMIN, class_gat%ULGAT(il1:il2))
            class_gat%FSSBGAT(il1:il2, 1) = class_gat%FSVHGAT(il1:il2)
            class_gat%FSSBGAT(il1:il2, 2) = class_gat%FSIHGAT(il1:il2)


        !> Check if we are on the first timestep of the day.
        if (ic%ts_daily == 1) then
            !> Check if this is the first day of the year.
            if (ic%now%jday == 1) then

                !> Check if this year is a leap year, and if so adjust the monthdays, monthend and mmday values.
                if (c_switch%leap) call findLeapYears(ic%now%year, leapnow, lastDOY)

                !> If needed, update values that were read in from the accessory input files (popd, wetlands, lightning...)
                if (c_switch%ctem_on) then

                    if (c_switch%lnduseon) then
!?              call updateInput('LUC',runyr)
                    else ! If landuse change is not on, then set the next years landcover to be
                        ! the same as this years.
                        vgat%nfcancmx = vgat%pfcancmx
                    end if

                    ! Reset the peatland degree days counter (JM:FLAG is this ok for S. Hemi to be Jan 1?)
                    vgat%pdd = 0
                end if
            end if ! first day
        end if   ! first timestep

      !> The cosine of the solar zenith angle COSZ is calculated from the day of
      !> the year, the hour, the minute and the latitude using basic radiation geometry,
      !> and (avoiding vanishingly small numbers) is assigned to CSZROW.  The fractional
      !> cloud cover FCLOROW is commonly not available so a rough estimate is
      !> obtained by setting it to 1 when precipitation is occurring, and to the fraction
      !> of incoming diffuse radiation XDIFFUS otherwise (assumed to be 1 when the sun
      !> is at the horizon, and 0.10 when it is at the zenith). These calculations are
      !> done in findCloudiness

        call findCloudinessGAT( &
            class_gat%RADJGAT, class_gat%CSZGAT, class_gat%PREGAT, vgat%xdiffusgat, class_gat%FCLOGAT, &
            ilg, il1, il2, ic%now%mins, ic%now%hour, ic%now%day, lastDOY)

      !> atmosphericVarsCalc evaluates a series of derived atmospheric variables

        call atmosphericVarsCalc( &
            class_gat%VPDGAT, class_gat%TADPGAT, class_gat%PADRGAT, class_gat%RHOAGAT, class_gat%RHSIGAT, & ! Formerly CLASSI
            class_gat%RPCPGAT, class_gat%TRPCGAT, class_gat%SPCPGAT, class_gat%TSPCGAT, class_gat%TAGAT, class_gat%QAGAT, &
            class_gat%PREGAT, class_gat%RPCPGAT, class_gat%SPCPGAT, class_gat%PRESGAT, &
            c_switch%IPCP, ilg, il1, il2)
		
      ! FLAG New for Jason's Albedo
      ! Calculate snow fall flux in kg/m2/s
        if (c_switch%ISNOALB == 1) then
            do k = il1, il2
                class_gat%PCSNGAT(k) = class_gat%PCSNGAT(k) + class_gat%SPCPGAT(k)* class_gat%RHSIGAT(k) * DELT  ! snowfall rate * density of fresh snow
            end do
        end if
      ! End flag

      !    * INITIALIZATION OF DIAGNOSTIC VARIABLES SPLIT OUT OF classGather
      !    * FOR CONSISTENCY WITH GCM APPLICATIONS.
        call initDiagnosticVars(ilg, ilg)

      !========================================================================

        if (ipid /= 0 .or. izero == 0) then

            !> Initialize diagnostic variables for PBSM.
            if (pbsm%PROCESS_ACTIVE) then
                pbsm%vs%Drift(il1:il2) = 0.0
                pbsm%vs%Subl(il1:il2) = 0.0
                pbsm%vs%Deposition(il1:il2) = 0.0
            end if
            ICE = 0.0
            TICE = 0.0

      !> energyWaterBalanceCheck does the initial calculations for the energy and water balance checks
      call energyWaterBalanceCheck(0, class_gat%CTVSTP, class_gat%CTSSTP, class_gat%CT1STP, class_gat%CT2STP, class_gat%CT3STP, & ! Formerly CLASSZ
                                   class_gat%WTVSTP, class_gat%WTSSTP, class_gat%WTGSTP, &
                                   class_gat%FSGVGAT, class_gat%FLGVGAT, class_gat%HFSCGAT, class_gat%HEVCGAT, class_gat%HMFCGAT, class_gat%HTCCGAT, &
                                   class_gat%FSGSGAT, class_gat%FLGSGAT, class_gat%HFSSGAT, class_gat%HEVSGAT, class_gat%HMFNGAT, class_gat%HTCSGAT, &
                                   class_gat%FSGGGAT, class_gat%FLGGGAT, class_gat%HFSGGAT, class_gat%HEVGGAT, class_gat%HMFGGAT, class_gat%HTCGAT, &
                                   class_gat%PCFCGAT, class_gat%PCLCGAT, class_gat%QFCFGAT, class_gat%QFCLGAT, class_gat%ROFCGAT, class_gat%WTRCGAT, &
                                   class_gat%PCPNGAT, class_gat%QFNGAT, class_gat%ROFNGAT, class_gat%WTRSGAT, class_gat%PCPGGAT, class_gat%QFGGAT, &
                                   class_gat%QFCGAT, class_gat%ROFGAT, class_gat%WTRGGAT, class_gat%CMAIGAT, class_gat%RCANGAT, class_gat%SCANGAT, &
                                   class_gat%TCANGAT, class_gat%SNOGAT, class_gat%WSNOGAT, class_gat%TSNOGAT, class_gat%THLQGAT, class_gat%THICGAT, &
                                   class_gat%HCPSGAT, class_gat%THPGAT, class_gat%DLZWGAT, class_gat%TBARGAT, class_gat%ZPNDGAT, class_gat%TPNDGAT, &
                                   class_gat%DELZ, class_gat%FCS, class_gat%FGS, class_gat%FC, class_gat%FG, &
                                   il1, il2, ILG, IGND, N)

      !> radiationDriver manages the calculation of albedos and other surface parameters
    call radiationDriver(class_gat%FC, class_gat%FG, class_gat%FCS, class_gat%FGS, class_gat%ALVSCN, class_gat%ALIRCN, & ! Formerly CLASSA
                           class_gat%ALVSG, class_gat%ALIRG, class_gat%ALVSCS, class_gat%ALIRCS, class_gat%ALVSSN, class_gat%ALIRSN, &
                           class_gat%ALVSGC, class_gat%ALIRGC, class_gat%ALVSSC, class_gat%ALIRSC, class_gat%TRVSCN, class_gat%TRIRCN, &
                           class_gat%TRVSCS, class_gat%TRIRCS, class_gat%FSVF, class_gat%FSVFS, &
                           class_gat%RAICAN, class_gat%RAICNS, class_gat%SNOCAN, class_gat%SNOCNS, class_gat%FRAINC, class_gat%FSNOWC, &
                           class_gat%FRAICS, class_gat%FSNOCS, class_gat%DISP, class_gat%DISPS, class_gat%ZOMLNC, class_gat%ZOMLCS, &
                           class_gat%ZOELNC, class_gat%ZOELCS, class_gat%ZOMLNG, class_gat%ZOMLNS, class_gat%ZOELNG, class_gat%ZOELNS, &
                           class_gat%CHCAP, class_gat%CHCAPS, class_gat%CMASSC, class_gat%CMASCS, class_gat%CWLCAP, class_gat%CWFCAP, &
                           class_gat%CWLCPS, class_gat%CWFCPS, class_gat%RC, class_gat%RCS, class_gat%RBCOEF, class_gat%FROOT, &
                           class_gat%FROOTS, class_gat%ZPLIMC, class_gat%ZPLIMG, class_gat%ZPLMCS, class_gat%ZPLMGS, class_gat%ZSNOW, &
                           class_gat%WSNOGAT, class_gat%ALVSGAT, class_gat%ALIRGAT, class_gat%HTCCGAT, class_gat%HTCSGAT, class_gat%HTCGAT, &
                           class_gat%ALTG, class_gat%ALSNO, class_gat%TRSNOWC, class_gat%TRSNOWG, &
                           class_gat%WTRCGAT, class_gat%WTRSGAT, class_gat%WTRGGAT, class_gat%CMAIGAT, class_gat%FSNOGAT, &
                           class_gat%FCANGAT, class_gat%LNZ0GAT, class_gat%ALVCGAT, class_gat%ALICGAT, class_gat%PAMXGAT, class_gat%PAMNGAT, &
                           class_gat%CMASGAT, class_gat%ROOTGAT, class_gat%RSMNGAT, class_gat%QA50GAT, class_gat%VPDAGAT, class_gat%VPDBGAT, &
                           class_gat%PSGAGAT, class_gat%PSGBGAT, class_gat%PAIDGAT, class_gat%HGTDGAT, class_gat%ACVDGAT, class_gat%ACIDGAT, &
                           class_gat%ASVDGAT, class_gat%ASIDGAT, class_gat%AGVDGAT, class_gat%AGIDGAT, &
                           class_gat%ALGWVGAT, class_gat%ALGWNGAT, class_gat%ALGDVGAT, class_gat%ALGDNGAT, &
                           class_gat%THLQGAT, class_gat%THICGAT, class_gat%TBARGAT, class_gat%RCANGAT, class_gat%SCANGAT, class_gat%TCANGAT, &
                           class_gat%GROGAT, class_gat%SNOGAT, class_gat%TSNOGAT, class_gat%RHOSGAT, class_gat%ALBSGAT, class_gat%ZBLDGAT, &
                           class_gat%Z0ORGAT, class_gat%ZSNLGAT, class_gat%ZPLGGAT, class_gat%ZPLSGAT, &
                           class_gat%FCLOGAT, class_gat%TAGAT, class_gat%VPDGAT, class_gat%RHOAGAT, class_gat%CSZGAT, &
                           class_gat%FSDBGAT, class_gat%FSFBGAT, class_gat%REFGAT, class_gat%BCSNGAT, &
                           class_gat%FSVHGAT, class_gat%RADJGAT, class_gat%DLONGAT, class_gat%RHSIGAT, class_gat%DELZ, class_gat%DLZWGAT, &
                           class_gat%ZBTWGAT, class_gat%THPGAT, class_gat%THMGAT, class_gat%PSISGAT, class_gat%BIGAT, class_gat%PSIWGAT, &
                           class_gat%HCPSGAT, class_gat%ISNDGAT, &
                           vgat%fcancmx, ICC, c_switch%ctem_on, vgat%rmatc, vgat%zolnc, &
                           vgat%cmasvegc, vgat%ailc, vgat%paic, NOL2PFTS, &
                           vgat%slaic, vgat%ailcg, vgat%ailcgs, vgat%fcanc, vgat%fcancs, &
                           ic%now%jday, ILG, 1, ilg, NBS, &
                           JLAT, N, ICAN, ICAN + 1, IGND, c_switch%IDISP, c_switch%IZREF, &
                           IWFGAT, c_switch%IPAI, c_switch%IHGT, c_switch%IALC, c_switch%IALS, c_switch%IALG, &
                           c_switch%ISNOALB, vgat%alvsctm, vgat%alirctm, vgat%ipeatland, &
                           il1, il2)

            !> FROZENSOILINFILFLAG 1.
            if (FROZENSOILINFILFLAG == 1) then

                !> bruce davison - only increase NMELT if we don't start the run on January 1st, otherwise t0_ACC allocation is too large
                !> and the model crashes if the compiler is checking for array bounds when t0_ACC is passed into CLASSW with size NMELT
                if (ic%now%jday == 1 .and. ic%ts_daily == 48) then
                    if (ic%start%jday == 1 .and. ic%ts_count < 49) then
                        continue ! NMELT should stay = 1
                    else
                        NMELT = NMELT + 1
                    end if
                    CUMSNOWINFILCS = 0.0
                    CUMSNOWINFILGS = 0.0
                    INFILTYPE = 2
                end if
            end if

      !> energyBudgetDriver calls the subroutines associated with the surface energy balance calculations
    call energyBudgetDriver(class_gat%TBARC, class_gat%TBARG, class_gat%TBARCS, class_gat%TBARGS, class_gat%THLIQC, class_gat%THLIQG, & ! Formerly CLASST
                              class_gat%THICEC, class_gat%THICEG, class_gat%HCPC, class_gat%HCPG, class_gat%TCTOPC, class_gat%TCBOTC, class_gat%TCTOPG, class_gat%TCBOTG, &
                              class_gat%GZEROC, class_gat%GZEROG, class_gat%GZROCS, class_gat%GZROGS, class_gat%G12C, class_gat%G12G, class_gat%G12CS, class_gat%G12GS, &
                              class_gat%G23C, class_gat%G23G, class_gat%G23CS, class_gat%G23GS, class_gat%QFREZC, class_gat%QFREZG, class_gat%QMELTC, class_gat%QMELTG, &
                              class_gat%EVAPC, class_gat%EVAPCG, class_gat%EVAPG, class_gat%EVAPCS, class_gat%EVPCSG, class_gat%EVAPGS, class_gat%TCANO, class_gat%TCANS, &
                              class_gat%RAICAN, class_gat%SNOCAN, class_gat%RAICNS, class_gat%SNOCNS, class_gat%CHCAP, class_gat%CHCAPS, class_gat%TPONDC, class_gat%TPONDG, &
                              class_gat%TPNDCS, class_gat%TPNDGS, class_gat%TSNOCS, class_gat%TSNOGS, class_gat%WSNOCS, class_gat%WSNOGS, class_gat%RHOSCS, class_gat%RHOSGS, &
                              class_gat%ITCTGAT, class_gat%CDHGAT, class_gat%CDMGAT, class_gat%HFSGAT, class_gat%TFXGAT, class_gat%QEVPGAT, class_gat%QFSGAT, &
                              class_gat%PETGAT, class_gat%GAGAT, class_gat%EFGAT, class_gat%GTGAT, class_gat%QGGAT, &
                              class_gat%SFCTGAT, class_gat%SFCUGAT, class_gat%SFCVGAT, class_gat%SFCQGAT, class_gat%SFRHGAT, &
                              class_gat%GTBS, class_gat%SFCUBS, class_gat%SFCVBS, class_gat%USTARBS, &
                              class_gat%FSGVGAT, class_gat%FSGSGAT, class_gat%FSGGGAT, class_gat%FLGVGAT, class_gat%FLGSGAT, class_gat%FLGGGAT, &
                              class_gat%HFSCGAT, class_gat%HFSSGAT, class_gat%HFSGGAT, class_gat%HEVCGAT, class_gat%HEVSGAT, class_gat%HEVGGAT, class_gat%HMFCGAT, class_gat%HMFNGAT, &
                              class_gat%HTCCGAT, class_gat%HTCSGAT, class_gat%HTCGAT, class_gat%QFCFGAT, class_gat%QFCLGAT, class_gat%DRGAT, class_gat%wtableGAT, class_gat%ILMOGAT, &
                              class_gat%UEGAT, class_gat%HBLGAT, class_gat%TACGAT, class_gat%QACGAT, class_gat%ZRFMGAT, class_gat%ZRFHGAT, class_gat%ZDMGAT, class_gat%ZDHGAT, &
                              class_gat%VPDGAT, class_gat%TADPGAT, class_gat%RHOAGAT, class_gat%FSVHGAT, class_gat%FSIHGAT, class_gat%FDLGAT, class_gat%ULGAT, class_gat%VLGAT, &
                              class_gat%TAGAT, class_gat%QAGAT, class_gat%PADRGAT, class_gat%FC, class_gat%FG, class_gat%FCS, class_gat%FGS, class_gat%RBCOEF, &
                              class_gat%FSVF, class_gat%FSVFS, class_gat%PRESGAT, class_gat%VMODGAT, class_gat%ALVSCN, class_gat%ALIRCN, class_gat%ALVSG, class_gat%ALIRG, &
                              class_gat%ALVSCS, class_gat%ALIRCS, class_gat%ALVSSN, class_gat%ALIRSN, class_gat%ALVSGC, class_gat%ALIRGC, class_gat%ALVSSC, class_gat%ALIRSC, &
                              class_gat%TRVSCN, class_gat%TRIRCN, class_gat%TRVSCS, class_gat%TRIRCS, class_gat%RC, class_gat%RCS, class_gat%WTRGGAT, class_gat%groundHeatFlux, class_gat%QLWOGAT, &
                              class_gat%FRAINC, class_gat%FSNOWC, class_gat%FRAICS, class_gat%FSNOCS, class_gat%CMASSC, class_gat%CMASCS, class_gat%DISP, class_gat%DISPS, &
                              class_gat%ZOMLNC, class_gat%ZOELNC, class_gat%ZOMLNG, class_gat%ZOELNG, class_gat%ZOMLCS, class_gat%ZOELCS, class_gat%ZOMLNS, class_gat%ZOELNS, &
                              class_gat%TBARGAT, class_gat%THLQGAT, class_gat%THICGAT, class_gat%TPNDGAT, class_gat%ZPNDGAT, class_gat%TBASGAT, class_gat%TCANGAT, class_gat%TSNOGAT, &
                              class_gat%ZSNOW, class_gat%RHOSGAT, class_gat%WSNOGAT, class_gat%THPGAT, class_gat%THRGAT, class_gat%THMGAT, class_gat%THFCGAT, class_gat%THLWGAT, &
                              class_gat%TRSNOWC, class_gat%TRSNOWG, class_gat%ALSNO, class_gat%FSSBGAT, class_gat%FROOT, class_gat%FROOTS, &
                              class_gat%RADJGAT, class_gat%PREGAT, class_gat%HCPSGAT, class_gat%TCSGAT, class_gat%TSFSGAT, class_gat%DELZ, class_gat%DLZWGAT, class_gat%ZBTWGAT, &
                              class_gat%FTEMP, class_gat%FVAP, class_gat%RIB, class_gat%ISNDGAT, &
                              vgat%ailcg, vgat%ailcgs, vgat%fcanc, vgat%fcancs, vgat%co2conc, vgat%co2i1cg, &
                              vgat%co2i1cs, vgat%co2i2cg, vgat%co2i2cs, class_gat%CSZGAT, vgat%xdiffusgat, vgat%slai, ICC, &
                              c_switch%ctem_on, vgat%rmatctem, vgat%fcancmx, L2MAX, NOL2PFTS, vgat%cfluxcg, &
                              vgat%cfluxcs, vgat%ancsveg, vgat%ancgveg, vgat%rmlcsveg, vgat%rmlcgveg, &
                              class_gat%TCSNOW, class_gat%GSNOW, c_switch%ITC, c_switch%ITCG, c_switch%ITG, ILG, 1, ilg, JLAT, N, ICAN, &
                              IGND, c_switch%IZREF, c_switch%ISLFD, NLANDCS, NLANDGS, NLANDC, NLANDG, NLANDI, &
                              NBS, c_switch%ISNOALB, vgat%dayl, vgat%dayl_max, &
                              vgat%ipeatland, vgat%ancsmoss, vgat%angsmoss, vgat%ancmoss, vgat%angmoss, &
                              vgat%rmlcsmoss, vgat%rmlgsmoss, vgat%rmlcmoss, vgat%rmlgmoss, &
                              vgat%Cmossmas, vgat%dmoss, ic%now%jday, vgat%pdd, &
                              il1, il2)

            !> Single column blowing snow calculations (PBSM).
            call PBSM_within_tile( &
                class_gat%ZSNOW, class_gat%WSNOGAT, class_gat%SNOGAT, class_gat%RHOSGAT, class_gat%TSNOGAT, class_gat%HTCSGAT, &
                class_gat%TSNOCS, class_gat%TSNOGS, &
                class_gat%RHOSCS, class_gat%RHOSGS,&
                class_gat%WSNOCS, class_gat%WSNOGS, &
                class_gat%FC, class_gat%FG, class_gat%FCS, class_gat%FGS, &
                class_gat%SFCTGAT, class_gat%SFCUGAT, class_gat%SFCQGAT, &
                class_gat%ZRFMGAT, class_gat%ZOMLCS, class_gat%ZOMLNS, &
                ilg, &
                fls, shd)

      !> waterBudgetDriver calls the subroutines associated with the surface water balance calculations
    call waterBudgetDriver(class_gat%THLQGAT, class_gat%THICGAT, class_gat%TBARGAT, class_gat%TCANGAT, class_gat%RCANGAT, class_gat%SCANGAT, & ! Formerly CLASSW
                             class_gat%ROFGAT, class_gat%TROFGAT, class_gat%SNOGAT, class_gat%TSNOGAT, class_gat%RHOSGAT, class_gat%ALBSGAT, &
                             class_gat%WSNOGAT, class_gat%ZPNDGAT, class_gat%TPNDGAT, class_gat%GROGAT, class_gat%TBASGAT, class_gat%GFLXGAT, &
                             class_gat%PCFCGAT, class_gat%PCLCGAT, class_gat%PCPNGAT, class_gat%PCPGGAT, class_gat%QFCFGAT, class_gat%QFCLGAT, &
                             class_gat%QFNGAT, class_gat%QFGGAT, class_gat%QFCGAT, class_gat%HMFCGAT, class_gat%HMFGGAT, class_gat%HMFNGAT, &
                             class_gat%HTCCGAT, class_gat%HTCSGAT, class_gat%HTCGAT, class_gat%ROFCGAT, class_gat%ROFNGAT, class_gat%ROVGGAT, &
                             class_gat%WTRSGAT, class_gat%WTRGGAT, class_gat%ROFOGAT, class_gat%ROFSGAT, class_gat%ROFBGAT, &
                             class_gat%TROOGAT, class_gat%TROSGAT, class_gat%TROBGAT, class_gat%QFSGAT, class_gat%QFXGAT, class_gat%RHOAGAT, &
                             class_gat%TBARC, class_gat%TBARG, class_gat%TBARCS, class_gat%TBARGS, class_gat%THLIQC, class_gat%THLIQG, &
                             class_gat%THICEC, class_gat%THICEG, class_gat%HCPC, class_gat%HCPG, class_gat%RPCPGAT, class_gat%TRPCGAT, &
                             class_gat%SPCPGAT, class_gat%TSPCGAT, class_gat%PREGAT, class_gat%TAGAT, class_gat%RHSIGAT, class_gat%GGEOGAT, &
                             class_gat%FC, class_gat%FG, class_gat%FCS, class_gat%FGS, class_gat%TPONDC, class_gat%TPONDG, &
                             class_gat%TPNDCS, class_gat%TPNDGS, class_gat%EVAPC, class_gat%EVAPCG, class_gat%EVAPG, class_gat%EVAPCS, &
                             class_gat%EVPCSG, class_gat%EVAPGS, class_gat%QFREZC, class_gat%QFREZG, class_gat%QMELTC, class_gat%QMELTG, &
                             class_gat%RAICAN, class_gat%SNOCAN, class_gat%RAICNS, class_gat%SNOCNS, class_gat%FSVF, class_gat%FSVFS, &
                             class_gat%CWLCAP, class_gat%CWFCAP, class_gat%CWLCPS, class_gat%CWFCPS, class_gat%TCANO, &
                             class_gat%TCANS, class_gat%CHCAP, class_gat%CHCAPS, class_gat%CMASSC, class_gat%CMASCS, class_gat%ZSNOW, &
                             class_gat%GZEROC, class_gat%GZEROG, class_gat%GZROCS, class_gat%GZROGS, class_gat%G12C, class_gat%G12G, &
                             class_gat%G12CS, class_gat%G12GS, class_gat%G23C, class_gat%G23G, class_gat%G23CS, class_gat%G23GS, &
                             class_gat%TSNOCS, class_gat%TSNOGS, class_gat%WSNOCS, class_gat%WSNOGS, class_gat%RHOSCS, class_gat%RHOSGS, &
                             class_gat%ZPLIMC, class_gat%ZPLIMG, class_gat%ZPLMCS, class_gat%ZPLMGS, class_gat%TSFSGAT, &
                             class_gat%TCTOPC, class_gat%TCBOTC, class_gat%TCTOPG, class_gat%TCBOTG, class_gat%FROOT, class_gat%FROOTS, &
                             class_gat%THPGAT, class_gat%THRGAT, class_gat%THMGAT, class_gat%BIGAT, class_gat%PSISGAT, class_gat%GRKSGAT, &
                             class_gat%THRAGAT, class_gat%THFCGAT, class_gat%DRNGAT, class_gat%HCPSGAT, class_gat%DELZ, &
                             class_gat%DLZWGAT, class_gat%ZBTWGAT, class_gat%XSLPGAT, XDGAT, class_gat%WFSFGAT, KSGAT, &
                             MANNGAT, DDGAT, &
                             class_gat%ISNDGAT, class_gat%IGDRGAT, &
                             IWFGAT, ILG, 1, ilg, N, &
                             JLAT, ICAN, IGND, IGND + 1, IGND + 2, &
                             NLANDCS, NLANDGS, NLANDC, NLANDG, NLANDI, &
                             il1, il2)

            !> Distribution of blowing snow mass between tiles (PBSM).
            call PBSM_within_grid( &
                class_gat%TSNOGAT, class_gat%ZSNOW, &
                class_gat%RHOSGAT, class_gat%SNOGAT, class_gat%TSNOCS, class_gat%RHOSCS, class_gat%TSNOGS, &
                class_gat%RHOSGS, &
                class_gat%GCGAT, vgat%faregat, &
                class_gat%WSNOCS, class_gat%WSNOGS, class_gat%FCS, class_gat%FGS, class_gat%FC, class_gat%FG, &
                class_gat%TROOGAT, class_gat%ROFOGAT, class_gat%TROFGAT, class_gat%ROFGAT, class_gat%ROFNGAT, class_gat%PCPGGAT, class_gat%HTCSGAT, class_gat%WSNOGAT, &
                ilg, &
                fls, shd)

            class_gat%ROFGAT = class_gat%ROFGAT - UMQ

      !> energyWaterBalanceCheck completes the energy and water balance checks for the current time step
    call energyWaterBalanceCheck(1, class_gat%CTVSTP, class_gat%CTSSTP, class_gat%CT1STP, class_gat%CT2STP, class_gat%CT3STP, & ! Formerly CLASSZ
                                   class_gat%WTVSTP, class_gat%WTSSTP, class_gat%WTGSTP, &
                                   class_gat%FSGVGAT, class_gat%FLGVGAT, class_gat%HFSCGAT, class_gat%HEVCGAT, class_gat%HMFCGAT, class_gat%HTCCGAT, &
                                   class_gat%FSGSGAT, class_gat%FLGSGAT, class_gat%HFSSGAT, class_gat%HEVSGAT, class_gat%HMFNGAT, class_gat%HTCSGAT, &
                                   class_gat%FSGGGAT, class_gat%FLGGGAT, class_gat%HFSGGAT, class_gat%HEVGGAT, class_gat%HMFGGAT, class_gat%HTCGAT, &
                                   class_gat%PCFCGAT, class_gat%PCLCGAT, class_gat%QFCFGAT, class_gat%QFCLGAT, class_gat%ROFCGAT, class_gat%WTRCGAT, &
                                   class_gat%PCPNGAT, class_gat%QFNGAT, class_gat%ROFNGAT, class_gat%WTRSGAT, class_gat%PCPGGAT, class_gat%QFGGAT, &
                                   class_gat%QFCGAT, class_gat%ROFGAT, class_gat%WTRGGAT, class_gat%CMAIGAT, class_gat%RCANGAT, class_gat%SCANGAT, &
                                   class_gat%TCANGAT, class_gat%SNOGAT, class_gat%WSNOGAT, class_gat%TSNOGAT, class_gat%THLQGAT, class_gat%THICGAT, &
                                   class_gat%HCPSGAT, class_gat%THPGAT, class_gat%DLZWGAT, class_gat%TBARGAT, class_gat%ZPNDGAT, class_gat%TPNDGAT, &
                                   class_gat%DELZ, class_gat%FCS, class_gat%FGS, class_gat%FC, class_gat%FG, &
                                   il1, il2, ILG, IGND, N)

        end if !(ipid /= 0 .or. izero == 0) then

        !> CLASS output files.
        !if (WF_NUM_POINTS > 0) call CLASSOUT_update_files(shd)
        !> Copy internal variables back to MESH.
        if (allocated(vs%tile%prern)) vs%tile%prern(il1:il2) = class_gat%RPCPGAT(il1:il2)*RHOW !from [m s-1] to [kg ms-2 s-1].
        if (allocated(vs%tile%presno)) vs%tile%presno(il1:il2) = class_gat%SPCPGAT(il1:il2)*class_gat%RHSIGAT(il1:il2) !from [m s-1] to [kg ms-2 s-1].
        vs%tile%lqwscan(il1:il2) = class_gat%RCANGAT(il1:il2)
        vs%tile%fzwscan(il1:il2) = class_gat%SCANGAT(il1:il2)
        vs%tile%cmas(il1:il2) = class_gat%CMAIGAT(il1:il2)
        vs%tile%tacan(il1:il2) = class_gat%TACGAT(il1:il2)
        vs%tile%tcan(il1:il2) = class_gat%TCANGAT(il1:il2)
        vs%tile%qacan(il1:il2) = class_gat%QACGAT(il1:il2)
        vs%tile%gro(il1:il2) = class_gat%GROGAT(il1:il2)
        vs%tile%sno(il1:il2) = class_gat%SNOGAT(il1:il2)
        vs%tile%albsno(il1:il2) = class_gat%ALBSGAT(il1:il2)
        vs%tile%fsno(il1:il2) = class_gat%FSNOGAT(il1:il2)
        vs%tile%drainsno(il1:il2) = class_gat%ROFNGAT(il1:il2)
        vs%tile%rhosno(il1:il2) = class_gat%RHOSGAT(il1:il2)
        vs%tile%lqwssno(il1:il2) = class_gat%WSNOGAT(il1:il2)
        vs%tile%tsno(il1:il2) = class_gat%TSNOGAT(il1:il2)
        vs%tile%alvs(il1:il2) = class_gat%ALVSGAT(il1:il2)
        vs%tile%alir(il1:il2) = class_gat%ALIRGAT(il1:il2)
        vs%tile%gte(il1:il2) = class_gat%GTGAT(il1:il2)
        vs%tile%zpnd(il1:il2) = class_gat%ZPNDGAT(il1:il2)
        vs%tile%tpnd(il1:il2) = class_gat%TPNDGAT(il1:il2)
        vs%tile%pndcaf(il1:il2) = &
            FSTRCS(il1:il2)*class_gat%FCS(il1:il2) + FSTRC(il1:il2)*class_gat%FC(il1:il2) + &
            FSTRG(il1:il2)*class_gat%FG(il1:il2) + FSTRGS(il1:il2)*class_gat%FGS(il1:il2)
        vs%tile%potevp(il1:il2) = class_gat%PETGAT(il1:il2)
        vs%tile%et(il1:il2) = class_gat%QFSGAT(il1:il2)
        vs%tile%ovrflw(il1:il2) = class_gat%ROFOGAT(il1:il2)
        vs%tile%qevp(il1:il2) = class_gat%QEVPGAT(il1:il2)
        vs%tile%qsens(il1:il2) = class_gat%HFSGAT(il1:il2)
        vs%tile%gzero(il1:il2) = class_gat%GFLXGAT(il1:il2, 1)
        vs%tile%tsfs(il1:il2, :) = class_gat%TSFSGAT(il1:il2, :)
        vs%tile%tsurf(il1:il2) = &
            class_gat%TSFSGAT(il1:il2, 1)*class_gat%FCS(il1:il2) + class_gat%TSFSGAT(il1:il2, 2)*class_gat%FGS(il1:il2) + &
            class_gat%TSFSGAT(il1:il2, 3)*class_gat%FC(il1:il2) + class_gat%TSFSGAT(il1:il2, 4)*class_gat%FG(il1:il2)
        vs%tile%tbas(il1:il2) = class_gat%TBASGAT(il1:il2)
        vs%tile%latflw(il1:il2, :) = class_gat%ROFSGAT(il1:il2, :)
        vs%tile%thicsol(il1:il2, :) = class_gat%THICGAT(il1:il2, :)
        vs%tile%thlqsol(il1:il2, :) = class_gat%THLQGAT(il1:il2, :)
        vs%tile%tsol(il1:il2, :) = class_gat%TBARGAT(il1:il2, :)
        vs%tile%gflx(il1:il2, :) = class_gat%GFLXGAT(il1:il2, :)
        vs%tile%drainsol(il1:il2) = class_gat%ROFBGAT(il1:il2)

      if (c_switch%ctem_on) then

        !> Accumulate variables not already accumulated but which are required by CTEM.
        call accumulateForCTEM(ilg, class_gat%ILMOS)

        if (ic%now%day /= ic%next%day) then

          ! Find daily averages of accumulated variables for CTEM
          call dayEndCTEMPreparation(ilg, ic%ts_daily, class_gat%ILMOS)

          ! Call Canadian Terrestrial Ecosystem Model which operates at a daily time step,
          ! and uses daily accumulated values of variables simulated by CLASS.
          call ctem(ctem_tile%fsnowacc_t, vgat%sandgat, class_gat%ILMOS, & ! In
                    ilg, 1, ilg, ic%now%jday, class_gat%RADJGAT, &! In
                    ctem_tile%taaccgat_t, class_gat%DLZWGAT, ctem_tile%ancgvgac_t, ctem_tile%rmlcgvga_t, & ! In
                    class_gat%ZBTWGAT, c_switch%doMethane, & ! In
                    ctem_tile%uvaccgat_t, ctem_tile%vvaccgat_t, vgat%lightng, ctem_tile%tbaraccgat_t, &! In
                    vgat%sdepgat, c_switch%spinfast, vgat%todfrac, & ! In
                    vgat%netrad_gat, vgat%preacc_gat, class_gat%PSISGAT, &! In
                    vgat%grclarea, vgat%popdin, class_gat%ISNDGAT, &! In
                    vgat%wetfrac_pres, vgat%slopefrac, class_gat%BIGAT, &! In
                    class_gat%THPGAT, class_gat%DLATGAT, vgat%ch4conc, &! In
                    class_gat%THFCGAT, class_gat%THLWGAT, ctem_tile%thliqacc_t, ctem_tile%thiceacc_t, &! In
                    vgat%ipeatland, ctem_tile%anmossac_t, ctem_tile%rmlmossac_t, ctem_tile%gppmossac_t, &! In
                    class_gat%wtableGAT, class_gat%maxAnnualActLyrGAT, & ! In
                    c_switch%PFTCompetition, c_switch%dofire, c_switch%lnduseon, c_switch%inibioclim, & ! In
                    leapnow, c_switch%useTracer, tracer%tracerCO2gat, &! In
                    vgat%pfcancmx, vgat%nfcancmx, & ! In
                    vgat%stemmass, vgat%rootmass, vgat%litrmass, vgat%gleafmas, & ! In/Out
                    vgat%bleafmas, vgat%soilcmas, vgat%ailcg, vgat%ailc, & ! In/Out
                    vgat%zolnc, vgat%rmatctem, vgat%rmatc, vgat%ailcb, & ! In/Out
                    vgat%flhrloss, vgat%pandays, vgat%lfstatus, vgat%grwtheff, & ! In/Out
                    vgat%lystmmas, vgat%lyrotmas, vgat%tymaxlai, vgat%vgbiomas, & ! In/Out
                    vgat%gavgltms, vgat%gavgscms, vgat%stmhrlos, vgat%slai, & ! In/Out
                    vgat%bmasveg, vgat%cmasvegc, vgat%colddays, vgat%rothrlos, & ! In/Out
                    class_gat%FCANGAT, vgat%alvsctm, vgat%alirctm, vgat%gavglai, &! In/Out
                    vgat%Cmossmas, vgat%litrmsmoss, vgat%peatdep, vgat%fcancmx, &! In/Out
                    vgat%geremort, vgat%intrmort, vgat%pstemmass, vgat%pgleafmass, &! In/Out
                    vgat%tcurm, vgat%srpcuryr, vgat%dftcuryr, vgat%lambda, &! In/Out
                    vgat%tmonth, vgat%anpcpcur, vgat%anpecur, vgat%gdd5cur, &! In/Out
                    vgat%surmncur, vgat%defmncur, vgat%srplscur, vgat%defctcur, &! In/Out
                    vgat%aridity, vgat%srplsmon, vgat%defctmon, vgat%anndefct, &! In/Out
                    vgat%annsrpls, vgat%annpcp, vgat%dry_season_length, &! In/Out
                    vgat%pftexist, vgat%twarmm, vgat%tcoldm, vgat%gdd5, vgat%nppveg, &! In/Out
                    tracer%stemMassgat, tracer%rootMassgat, tracer%gLeafMassgat, tracer%bLeafMassgat, & ! In/Out
                    tracer%soilCMassgat, tracer%litrMassgat, tracer%mossCMassgat, tracer%mossLitrMassgat, & ! In/Out
                    vgat%npp, vgat%nep, vgat%hetrores, vgat%autores, &! Out (Primary)
                    vgat%soilcresp, vgat%rm, vgat%rg, vgat%nbp, &! Out (Primary)
                    vgat%litres, vgat%socres, vgat%gpp, vgat%dstcemls, &! Out (Primary)
                    vgat%litrfall, vgat%humiftrs, vgat%veghght, vgat%rootdpth, &! Out (Primary)
                    vgat%rml, vgat%rms, vgat%rmr, vgat%tltrleaf, &! Out (Primary)
                    vgat%tltrstem, vgat%tltrroot, vgat%leaflitr, vgat%roottemp, &! Out (Primary)
                    vgat%burnfrac, vgat%lucemcom, vgat%lucltrin, &! Out (Primary)
                    vgat%lucsocin, vgat%dstcemls3, &! Out (Primary)
                    vgat%ch4WetSpec, vgat%ch4WetDyn, vgat%wetfdyn, vgat%ch4_soills, &! Out (Primary)
                    vgat%paic, vgat%slaic, &! Out (Primary)
                    vgat%emit_co2, vgat%emit_ch4, vgat%reprocost, vgat%blfltrdt, vgat%glfltrdt, &! Out (Primary)
                    vgat%glcaemls, vgat%blcaemls, vgat%rtcaemls, vgat%stcaemls, vgat%ltrcemls, &  ! Out (Primary)
                    vgat%ntchlveg, vgat%ntchsveg, vgat%ntchrveg, &  ! Out (Primary)
                    vgat%emit_co, vgat%emit_nmhc, vgat%smfuncveg, &! Out (Secondary)
                    vgat%emit_h2, vgat%emit_nox, vgat%emit_n2o, vgat%emit_pm25, &! Out (Secondary)
                    vgat%emit_tpm, vgat%emit_tc, vgat%emit_oc, vgat%emit_bc, &! Out (Secondary)
                    vgat%bterm, vgat%lterm, vgat%mterm, vgat%burnvegf, &! Out (Secondary)
                    vgat%litrfallveg, vgat%humiftrsveg, vgat%ltstatus,  &! Out (Secondary)
                    vgat%afrleaf, vgat%afrstem, vgat%afrroot, vgat%wtstatus, &! Out (Secondary)
                    vgat%rmlvegacc, vgat%rmsveg, vgat%rmrveg, vgat%rgveg, &! Out (Secondary)
                    vgat%vgbiomas_veg, vgat%gppveg, vgat%nepveg, vgat%nbpveg, &! Out (Secondary)
                    vgat%hetroresveg, vgat%autoresveg, vgat%litresveg, vgat%soilcresveg, &! Out (Secondary)
                    vgat%nppmoss, vgat%armoss, &! Out (Secondary)
                    vgat%cc, vgat%mm)

                    !> MESH variables.
					!dm! *convertkgC
                    vs%tile%nep(il1:il2) = vgat%nep(il1:il2)*0.216
                    vs%tile%nbp(il1:il2) = vgat%nbp(il1:il2)*0.216
                    vs%tile%npp(il1:il2) = vgat%npp(il1:il2)*0.216
                    vs%tile%gpp(il1:il2) = vgat%gpp(il1:il2)*0.216
					vs%tile%Re(il1:il2) = vgat%autores(il1:il2)*0.216 + vgat%hetrores(il1:il2)*0.216
					
					!print *, "ailcg: ", vgat%ailcg
					!print *, "lfstatus: ", vgat%lfstatus
					!print *, "stemmass: ", vgat%stemmass
					!print *, "rootmass: ", vgat%rootmass
					
          !     reset mosaic accumulator arrays. These are scattered in ctems2 so we need
          !     to reset here,prior to ctems2.
          do i = 1,ilg
            ctem_tile%vvaccgat_t(i) = 0.0  !
            ctem_tile%uvaccgat_t(i) = 0.0  !
            do j = 1,ignd
              ctem_tile%tbaraccgat_t(i,j) = 0.0 !
            end do
          end do

          ! Once a year,calculate the 14C lost to decay if using the 14C tracer.
          if (c_switch%useTracer == 2 .and. &
              ic%now%jday /= ic%next%jday .and. ic%now%day /= ic%next%day) call decay14C(1,ilg)

        end if  ! if (ncount==nday)
      end if  ! if (ctem_on)

      if (c_switch%ISNOALB == 1) then
        call fourBandDriver(ilg, class_gat%SNOGAT, class_gat%GCGAT, class_gat%TSNOGAT, & ! in
                            class_gat%ZSNOW, class_gat%TCSNOW, class_gat%GSNOW, class_gat%PCSNGAT, class_gat%WSNOGAT, & ! in
                            class_gat%ALBSGAT, class_gat%RHOSGAT, & ! in/out
                            class_gat%TZSGAT, class_gat%REFGAT) ! out
      end if

        if (ic%now%day /= ic%next%day) then

            !> Reset mosaic accumulator arrays.
            if (c_switch%ctem_on) call resetMosaicAccum

            !> Determine the active layer depth and depth to the frozen water table.
            !> This only occurs once per day since they don't change rapidly.
            call findPermafrostVarsGAT( &
                class_gat%ftable, class_gat%actlyr, class_gat%maxAnnualActLyrGAT, class_gat%actLyrThisYrGAT, &
                class_gat%TBARGAT, class_gat%THLQGAT, class_gat%THICGAT, class_gat%ISNDGAT, class_gat%DLZWGAT, class_gat%DELZ, class_gat%THMGAT, class_gat%DLATGAT, &
                ilg, il1, il2, ignd, ic%now%jday)
      end if ! ncount eq nday

    end subroutine

    subroutine CLASSICv1_0_finalize(shd, fls)

        !> Modules.
        use ctemStateVars, only: c_switch
        use shd_variables
        use model_files_variables

        !> Variables.
        type(ShedGridParams) :: shd
        type(fl_ids) :: fls

    end subroutine

    end module