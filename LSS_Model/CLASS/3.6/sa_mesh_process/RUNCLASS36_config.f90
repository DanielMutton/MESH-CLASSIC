module RUNCLASS36_config

    use RUNCLASS36_constants
    use RUNCLASS36_variables

    implicit none

    !>  CONSTANTS AND TEMPORARY VARIABLES.
    real FSDOWN1, FSDOWN2, FSDOWN3, RDAY, &
        DECL, HOUR, COSZ

    integer NLANDCS, NLANDGS, NLANDC, NLANDG, NLANDI
    real, dimension(:), allocatable :: &
        FRZCGAT

    !> CANOPY AND SOIL INFORMATION (CLASS):
    !> THE LENGTH OF THESE ARRAYS IS DETERMINED BY THE NUMBER
    !> OF SOIL LAYERS (3) AND THE NUMBER OF BROAD VEGETATION
    !> CATEGORIES (4, OR 5 INCLUDING URBAN AREAS).
    !* ALL: DEFINITIONS IN CLASS DOCUMENTATION (CLASS.INI)
    real, dimension(:), allocatable :: XDGAT, &
        KSGAT

    !> ATMOSPHERIC AND GRID-CONSTANT INPUT VARIABLES:
    real, dimension(:), allocatable :: ZDMGRD, &
        ZDHGRD, RADJGRD, CSZGRD, &
        PADRGRD, VPDGRD, &
        TADPGRD, RHOAGRD, RPCPGRD, TRPCGRD, SPCPGRD, TSPCGRD, RHSIGRD, &
        FCLOGRD, DLONGRD, Z0ORGRD, GGEOGRD, UVGRD, XDIFFUS, &
        RPREGRD, SPREGRD, VMODGRD

    !> LAND SURFACE DIAGNOSTIC VARIABLES:
    real, dimension(:), allocatable :: &
        SFRHGAT, &
        QLWOGAT, FTEMP, &
        FVAP, RIB
    real, dimension(:), allocatable :: CDHGRD, CDMGRD, HFSGRD, &
        TFXGRD, QEVPGRD, QFSGRD, QFXGRD, PETGRD, GAGRD, EFGRD, GTGRD, &
        QGGRD, TSFGRD, ALVSGRD, ALIRGRD, FSNOGRD, SFCTGRD, SFCUGRD, &
        SFCVGRD, SFCQGRD, FSGVGRD, FSGSGRD, FSGGGRD, FLGVGRD, FLGSGRD, &
        FLGGGRD, HFSCGRD, HFSSGRD, HFSGGRD, HEVCGRD, HEVSGRD, HEVGGRD, &
        HMFCGRD, HMFNGRD, HTCCGRD, HTCSGRD, PCFCGRD, PCLCGRD, PCPNGRD, &
        PCPGGRD, QFGGRD, QFNGRD, QFCLGRD, QFCFGRD, ROFGRD, ROFOGRD, &
        ROFSGRD, ROFBGRD, ROFCGRD, ROFNGRD, ROVGGRD, WTRCGRD, WTRSGRD, &
        WTRGGRD, DRGRD, WTABGRD, ILMOGRD, UEGRD, HBLGRD

    real, dimension(:, :), allocatable :: HMFGGRD, HTCGRD, QFCGRD, GFLXGRD

    !> CROSS-CLASS VARIABLES (CLASS):
    !> ARRAYS DEFINED TO PASS INFORMATION BETWEEN THE THREE MAJOR
    !> SUBSECTIONS OF CLASS ("CLASSA", "CLASST" AND "CLASSW").
    real, dimension(:), allocatable :: RBCOEF, &
        ZSNOW, FSVF, FSVFS, ALVSCN, ALIRCN, ALVSG, &
        ALIRG, ALVSCS, ALIRCS, ALVSSN, ALIRSN, ALVSGC, ALIRGC, ALVSSC, &
        ALIRSC, TRVSCN, TRIRCN, TRVSCS, TRIRCS, RC, RCS, FRAINC, &
        FSNOWC, FRAICS, FSNOCS, CMASSC, CMASCS, DISP, DISPS, ZOMLNC, &
        ZOELNC, ZOMLNG, &
        ZOELNG, ZOMLCS, ZOELCS, ZOMLNS, ZOELNS, TRSNOW, CHCAP, CHCAPS, &
        GZEROC, GZEROG, GZROCS, GZROGS, G12C, G12G, G12CS, G12GS, G23C, &
        G23G, G23CS, G23GS, QFREZC, QFREZG, QMELTC, QMELTG, EVAPC, &
        EVAPCG, EVAPG, EVAPCS, EVPCSG, EVAPGS, TCANO, TCANS, RAICAN, &
        SNOCAN, RAICNS, SNOCNS, CWLCAP, CWFCAP, CWLCPS, CWFCPS, TSNOCS, &
        TSNOGS, RHOSCS, RHOSGS, WSNOCS, WSNOGS, TPONDC, TPONDG, TPNDCS, &
        TPNDGS, ZPLMCS, ZPLMGS, ZPLIMC, ZPLIMG
    real, dimension(:, :), allocatable :: TBARC, TBARG, TBARCS, &
        TBARGS, THLIQC, THLIQG, THICEC, THICEG, FROOT, HCPC, HCPG, &
        TCTOPC, TCBOTC, TCTOPG, TCBOTG

    !> BALANCE ERRORS (CLASS):
    !> DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
    !> BALANCES.
    real, dimension(:), allocatable :: CTVSTP, CTSSTP, CT1STP, &
        CT2STP, CT3STP, WTVSTP, WTSSTP, WTGSTP

    character(len = 2), dimension(5) :: CANOPY_LEVEL_MAP = (/'NL', 'BL', 'CR', 'GR', 'UR'/)

    contains

    subroutine RUNCLASS36_init(shd, fls)

        use mpi_module
        use model_files_variables
        use sa_mesh_common
        use model_dates
        use FLAGS
        use mesh_io_constants
        use input_fields
        use parse_utilities, only: to_lowercase

        !> For CLASS output.
        use RUNCLASS36_save_output

        use PBSM_module

        type(ShedGridParams) :: shd
        type(fl_ids) :: fls

        integer NA, NTYPE, NML, NSL, z, l, k, ik, jk, m, j, i, iun, ierr
        real FRAC, factor_a, factor_m, valid_max, valid_min
        character(len = DEFAULT_LINE_LENGTH) CLASSOUTTITLE, CLASSOUTNAME, CLASSOUTPLACE, line
        logical, allocatable :: value_mask(:)

        !> Return if the process is not active.
        if (.not. RUNCLASS36_flgs%PROCESS_ACTIVE) then
            return
        else
            call print_new_section("CLASS (RUNCLASS36) is active.")
            call increase_tab()
        end if

        !> Check for required variables.
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
        if (IPCP == 4) then
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

        !> Local variables.
        NA = shd%NA
        NTYPE = shd%lc%NTYPE
        NSL = shd%lc%IGND
        NML = shd%lc%NML
        DELT = ic%dts
        if (DELT > 1800) then

            !> Print a warning if the time-step exceeds recommendations from the CLASS manual.
            write(line, *) DELT*60
            call print_warning( &
                "CLASS is designed to run at a time-step of 30 minutes or less. The current time-step of " // &
                trim(adjustl(line)) // " minutes may lead to the appearance of numerical instabilities in simulated states.")
        end if

        !> ALLOCATE ALL VARIABLES

1114 format(/1x, 'Error allocating ', a, ' variables.', &
            /1x, 'Check that these bounds are within an acceptable range.', /)
1118 format(3x, a, ': ', i6)

        allocate(XDGAT(NML), KSGAT(NML), stat = ierr)

        if (ierr /= 0) then
            print 1114, 'canopy and soil info.'
            print 1118, 'Grid squares', NA
            print 1118, 'GRUs', NTYPE
            print 1118, 'Total tile elements', NML
            print 1118, 'Canopy types with urban areas', ICP1
            print 1118, 'Canopy types', ICAN
            print 1118, 'Soil layers', NSL
            stop
        end if

        !> WATROF FLAGS AND VARIABLES:
        allocate(DDGAT(NML), MANNGAT(NML), stat = ierr)
        if (ierr /= 0) then
            print 1114, 'WATROF'
            print 1118, 'Grid squares', NA
            print 1118, 'GRUs', NTYPE
            print 1118, 'Total tile elements', NML
            stop
        end if

        !> ATMOSPHERIC AND GRID-CONSTANT INPUT VARIABLES:
        allocate(ZDMGRD(NA), &
                 ZDHGRD(NA), RADJGRD(NA), &
                 CSZGRD(NA), &
                 PADRGRD(NA), VPDGRD(NA), &
                 TADPGRD(NA), RHOAGRD(NA), RPCPGRD(NA), TRPCGRD(NA), &
                 SPCPGRD(NA), TSPCGRD(NA), RHSIGRD(NA), &
                 FCLOGRD(NA), DLONGRD(NA), Z0ORGRD(NA), GGEOGRD(NA), UVGRD(NA), &
                 XDIFFUS(NA), &
                 RPREGRD(NA), SPREGRD(NA), VMODGRD(NA), &
                 stat = ierr)
        if (ierr /= 0) then
            print 1114, 'atmospheric and grid-cst.'
            print 1118, 'Grid squares', NA
            print 1118, 'GRUs', NTYPE
            print 1118, 'Total tile elements', NML
            stop
        end if

        !> LAND SURFACE DIAGNOSTIC VARIABLES:
        allocate(SFRHGAT(NML), &
                 QLWOGAT(NML), &
                 FTEMP(NML), FVAP(NML), RIB(NML), &
                 CDHGRD(NA), CDMGRD(NA), HFSGRD(NA), &
                 TFXGRD(NA), QEVPGRD(NA), QFSGRD(NA), QFXGRD(NA), PETGRD(NA), &
                 GAGRD(NA), EFGRD(NA), GTGRD(NA), &
                 QGGRD(NA), TSFGRD(NA), ALVSGRD(NA), ALIRGRD(NA), FSNOGRD(NA), &
                 SFCTGRD(NA), SFCUGRD(NA), &
                 SFCVGRD(NA), SFCQGRD(NA), FSGVGRD(NA), FSGSGRD(NA), &
                 FSGGGRD(NA), FLGVGRD(NA), FLGSGRD(NA), &
                 FLGGGRD(NA), HFSCGRD(NA), HFSSGRD(NA), HFSGGRD(NA), &
                 HEVCGRD(NA), HEVSGRD(NA), HEVGGRD(NA), &
                 HMFCGRD(NA), HMFNGRD(NA), HTCCGRD(NA), HTCSGRD(NA), &
                 PCFCGRD(NA), PCLCGRD(NA), PCPNGRD(NA), &
                 PCPGGRD(NA), QFGGRD(NA), QFNGRD(NA), QFCLGRD(NA), QFCFGRD(NA), &
                 ROFGRD(NA), ROFOGRD(NA), &
                 ROFSGRD(NA), ROFBGRD(NA), ROFCGRD(NA), ROFNGRD(NA), &
                 ROVGGRD(NA), WTRCGRD(NA), WTRSGRD(NA), &
                 WTRGGRD(NA), DRGRD(NA), WTABGRD(NA), ILMOGRD(NA), UEGRD(NA), &
                 HBLGRD(NA), &
                 HMFGGRD(NA, NSL), HTCGRD(NA, NSL), QFCGRD(NA, NSL), &
                 GFLXGRD(NA, NSL), stat = ierr)
        if (ierr /= 0) then
            print 1114, 'land surface diagnostic'
            print 1118, 'Grid squares', NA
            print 1118, 'GRUs', NTYPE
            print 1118, 'Total tile elements', NML
            print 1118, 'Soil layers', NSL
            stop
        end if

        !> CROSS-CLASS VARIABLES (CLASS):
        allocate(TBARC(NML, NSL), TBARG(NML, NSL), &
                 TBARCS(NML, NSL), &
                 TBARGS(NML, NSL), THLIQC(NML, NSL), &
                 THLIQG(NML, NSL), THICEC(NML, NSL), &
                 THICEG(NML, NSL), FROOT(NML, NSL), &
                 HCPC(NML, NSL), HCPG(NML, NSL), &
                 TCTOPC(NML, NSL), TCBOTC(NML, NSL), &
                 TCTOPG(NML, NSL), TCBOTG(NML, NSL), &
                 RBCOEF(NML), &
                 ZSNOW(NML), &
                 FSVF(NML), FSVFS(NML), ALVSCN(NML), &
                 ALIRCN(NML), ALVSG(NML), &
                 ALIRG(NML), ALVSCS(NML), ALIRCS(NML), &
                 ALVSSN(NML), ALIRSN(NML), ALVSGC(NML), &
                 ALIRGC(NML), ALVSSC(NML), &
                 ALIRSC(NML), TRVSCN(NML), TRIRCN(NML), &
                 TRVSCS(NML), TRIRCS(NML), RC(NML), &
                 RCS(NML), FRAINC(NML), &
                 FSNOWC(NML),FRAICS(NML),FSNOCS(NML), &
                 CMASSC(NML), CMASCS(NML), &
                 DISP(NML), DISPS(NML), ZOMLNC(NML), &
                 ZOELNC(NML), ZOMLNG(NML), &
                 ZOELNG(NML), ZOMLCS(NML), ZOELCS(NML), &
                 ZOMLNS(NML), ZOELNS(NML), TRSNOW(NML), &
                 CHCAP(NML), CHCAPS(NML), &
                 GZEROC(NML), GZEROG(NML), GZROCS(NML), &
                 GZROGS(NML), G12C(NML), G12G(NML), &
                 G12CS(NML), G12GS(NML), G23C(NML), &
                 G23G(NML), G23CS(NML), G23GS(NML), &
                 QFREZC(NML), QFREZG(NML), QMELTC(NML), &
                 QMELTG(NML), EVAPC(NML), &
                 EVAPCG(NML), EVAPG(NML), EVAPCS(NML), &
                 EVPCSG(NML), EVAPGS(NML), TCANO(NML), &
                 TCANS(NML), RAICAN(NML), &
                 SNOCAN(NML), RAICNS(NML), SNOCNS(NML), &
                 CWLCAP(NML), CWFCAP(NML), CWLCPS(NML), &
                 CWFCPS(NML), TSNOCS(NML), &
                 TSNOGS(NML), RHOSCS(NML), RHOSGS(NML), &
                 WSNOCS(NML), WSNOGS(NML), TPONDC(NML), &
                 TPONDG(NML), TPNDCS(NML), &
                 TPNDGS(NML), ZPLMCS(NML), ZPLMGS(NML), &
                 ZPLIMC(NML), ZPLIMG(NML), stat = ierr)
        if (ierr /= 0) then
            print 1114, 'cross-CLASS'
            print 1118, 'Grid squares', NA
            print 1118, 'GRUs', NTYPE
            print 1118, 'Total tile elements', NML
            print 1118, 'Soil layers', NSL
            stop
        end if

        !> BALANCE ERRORS (CLASS):
        allocate(CTVSTP(NML), CTSSTP(NML), &
                 CT1STP(NML), &
                 CT2STP(NML), CT3STP(NML), WTVSTP(NML), &
                 WTSSTP(NML), WTGSTP(NML), stat = ierr)
        if (ierr /= 0) then
            print 1114, 'balance error diagnostic'
            print 1118, 'Grid squares', NA
            print 1118, 'GRUs', NTYPE
            print 1118, 'Total tile elements', NML
            stop
        end if

        !> CTEM ERRORS (CLASS):
        allocate(CO2CONC(NML), COSZS(NML), XDIFFUSC(NML), CFLUXCG(NML), CFLUXCS(NML), &
                 AILCG(NML, ICTEM), AILCGS(NML, ICTEM), FCANC(NML, ICTEM), FCANCS(NML, ICTEM), &
                 CO2I1CG(NML, ICTEM), CO2I1CS(NML, ICTEM), CO2I2CG(NML, ICTEM), CO2I2CS(NML, ICTEM), &
                 SLAI(NML, ICTEM), FCANCMX(NML, ICTEM), ANCSVEG(NML, ICTEM), ANCGVEG(NML, ICTEM), &
                 RMLCSVEG(NML, ICTEM), RMLCGVEG(NML, ICTEM), &
                 AILC(NML, ICAN), PAIC(NML, ICAN), FIELDSM(NML, NSL), WILTSM(NML, NSL), &
                 RMATCTEM(NML, ICTEM, NSL), RMATC(NML, ICAN, NSL), NOL2PFTS(ICAN), stat = ierr)
        if (ierr /= 0) then
            print 1114, 'CTEM'
            print 1118, 'Grid squares', NA
            print 1118, 'GRUs', NTYPE
            print 1118, 'Total tile elements', NML
            print 1118, 'Canopy types', ICAN
            print 1118, 'Soil layers', NSL
            print 1118, 'CTEM flag', ICTEM
            stop
        end if

        !> Forcing input.
        allocate(cfi%FDL(NML), cfi%FSIH(NML), cfi%FSVH(NML), cfi%PRE(NML), cfi%PRES(NML), cfi%QA(NML), cfi%TA(NML), cfi%UL(NML), &
                 cfi%VL(NML), cfi%VMOD(NML))

        !> Prognostic variables.
        allocate(cpv%ALBS(NML), cpv%CMAI(NML), cpv%GRO(NML), cpv%QAC(NML), cpv%RCAN(NML), cpv%RHOS(NML), cpv%SNCAN(NML), &
                 cpv%SNO(NML), cpv%TAC(NML), cpv%TBAS(NML), cpv%TCAN(NML), cpv%TPND(NML), cpv%TSNO(NML), cpv%WSNO(NML), &
                 cpv%ZPND(NML))
        allocate(cpv%TBAR(NML, NSL), cpv%THIC(NML, NSL), cpv%THLQ(NML, NSL))
        allocate(cpv%TSFS(NML, 4))

        !> Land-surface variables.
        allocate(csfv%AGID(NML), csfv%AGVD(NML), csfv%ALGD(NML), csfv%ALGW(NML), csfv%ASID(NML), csfv%ASVD(NML), csfv%DRN(NML), &
                 csfv%FARE(NML), csfv%GRKF(NML), csfv%MID(NML), csfv%IWF(NML), csfv%SDEP(NML), csfv%WFCI(NML), csfv%WFSF(NML), &
                 csfv%XSLP(NML), csfv%ZPLG(NML), csfv%ZPLS(NML), csfv%ZSNL(NML))
        allocate(csfv%IGDR(NML))
        allocate(csfv%IORG(NML, NSL), csfv%ISND(NML, NSL))
        allocate(csfv%BI(NML, NSL), csfv%CLAY(NML, NSL), csfv%DELZW(NML, NSL), csfv%GRKS(NML, NSL), csfv%HCPS(NML, NSL), &
                 csfv%ORGM(NML, NSL), csfv%PSIS(NML, NSL), csfv%PSIW(NML, NSL), csfv%SAND(NML, NSL), csfv%TCS(NML, NSL), &
                 csfv%THFC(NML, NSL), csfv%THM(NML, NSL), csfv%THP(NML, NSL), csfv%THR(NML, NSL), csfv%THRA(NML, NSL), &
                 csfv%ZBTW(NML, NSL))
        allocate(csfv%ACID(NML, ICAN), csfv%ACVD(NML, ICAN), csfv%CMAS(NML, ICAN), csfv%HGTD(NML, ICAN), csfv%PAID(NML, ICAN), &
                 csfv%PAMN(NML, ICAN), csfv%PAMX(NML, ICAN), csfv%PSGA(NML, ICAN), csfv%PSGB(NML, ICAN), csfv%QA50(NML, ICAN), &
                 csfv%ROOT(NML, ICAN), csfv%RSMN(NML, ICAN), csfv%VPDA(NML, ICAN), csfv%VPDB(NML, ICAN))
        allocate(csfv%ALIC(NML, ICP1), csfv%ALVC(NML, ICP1), csfv%FCAN(NML, ICP1), csfv%LNZ0(NML, ICP1))

        !> Glacier variables.
        allocate(cglv%FREZTH(NML), cglv%SWELIM(NML), cglv%SNDENLIM(NML))

        !> Atmospheric variables.
        allocate(catv%CSZ(NML), catv%DLON(NML), catv%FCLO(NML), catv%GC(NML), catv%GGEO(NML), catv%PADR(NML), catv%RADJ(NML), &
                 catv%RHOA(NML), catv%RHSI(NML), catv%RPCP(NML), catv%RPRE(NML), catv%SPCP(NML), catv%SPRE(NML), catv%TADP(NML), &
                 catv%TRPC(NML), catv%TSPC(NML), catv%VPD(NML), catv%Z0OR(NML), catv%ZBLD(NML), catv%ZDH(NML), catv%ZDM(NML), &
                 catv%ZRFH(NML), catv%ZRFM(NML))

        !> Diagnostic variables.
        allocate(cdv%ITCT(NML, 6, 50))
        allocate(cdv%ALIR(NML), cdv%ALVS(NML), cdv%CDH(NML), cdv%CDM(NML), cdv%DR(NML), cdv%EF(NML), cdv%FCS(NML), cdv%FGS(NML), &
                 cdv%FC(NML), cdv%FG(NML), cdv%FLGG(NML), cdv%FLGS(NML), cdv%FLGV(NML), cdv%FSGG(NML), cdv%FSGS(NML), &
                 cdv%FSGV(NML), cdv%FSNO(NML), cdv%GA(NML), cdv%GTE(NML), cdv%HBL(NML), cdv%HEVC(NML), cdv%HEVG(NML), &
                 cdv%HEVS(NML), cdv%HFS(NML), cdv%HFSC(NML), cdv%HFSG(NML), cdv%HFSS(NML), cdv%HMFC(NML), cdv%HMFN(NML), &
                 cdv%HTCC(NML), cdv%HTCS(NML), cdv%ILMO(NML), cdv%PCFC(NML), cdv%PCLC(NML), cdv%PCPG(NML), cdv%PCPN(NML), &
                 cdv%PET(NML), cdv%QEVP(NML), cdv%QFCF(NML), cdv%QFCL(NML), cdv%QFG(NML), cdv%QFN(NML), cdv%QFS(NML), &
                 cdv%QFX(NML), cdv%QG(NML), cdv%ROF(NML), cdv%ROFB(NML), cdv%ROFC(NML), cdv%ROFN(NML), cdv%ROFO(NML), &
                 cdv%ROVG(NML), cdv%SFCQ(NML), cdv%SFCT(NML), cdv%SFCU(NML), cdv%SFCV(NML), cdv%TFX(NML), &
                 cdv%TROB(NML), cdv%TROF(NML), cdv%TROO(NML), cdv%TSF(NML), cdv%UE(NML), cdv%WTAB(NML), &
                 cdv%WTRC(NML), cdv%WTRG(NML), cdv%WTRS(NML), cdv%ICE(NML), cdv%TICE(NML))
        allocate(cdv%ROFS(NML, NSL), cdv%TROS(NML, NSL), &
                 cdv%GFLX(NML, NSL), cdv%HMFG(NML, NSL), cdv%HTC(NML, NSL), cdv%QFC(NML, NSL))

        !> Call 'CLASSD' to initialize constants.
        call CLASSD

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
            allocate(pbsm%vs%Drift(NML), pbsm%vs%Subl(NML))
            pbsm%vs%Drift = 0.0; pbsm%vs%Subl = 0.0

            !> Variables used in CLASSW.
            !> These are initialized in WPREP.
            allocate(ZSNOCS(NML), ZSNOGS(NML), ZSNOWC(NML), ZSNOWG(NML), &
                     HCPSCS(NML), HCPSGS(NML), HCPSC(NML), HCPSG(NML), &
                     TSNOWC(NML), TSNOWG(NML), &
                     RHOSC(NML), RHOSG(NML), &
                     XSNOWC(NML), XSNOWG(NML), XSNOCS(NML), XSNOGS(NML))
        end if

        !> Distribute variables.
!-        catv%ZRFM(il1:il2) = pm%tile%zrfm(il1:il2)
!-        catv%ZRFH(il1:il2) = pm%tile%zrfh(il1:il2)
!-        catv%ZBLD(il1:il2) = pm%tile%zbld(il1:il2)
!-        catv%GC(il1:il2) = pm%tile%gc(il1:il2)
!-        csfv%FARE(il1:il2) = pm%tile%fare(il1:il2)
!-        csfv%MID(il1:il2) = pm%tile%mid(il1:il2)
!-        csfv%IWF(il1:il2) = pm%tile%iwf(il1:il2)
!-        csfv%FCAN(il1:il2, :) = pm%tile%fcan(il1:il2, :)
!-        csfv%LNZ0(il1:il2, :) = pm%tile%lnz0(il1:il2, :)
!-        csfv%ALVC(il1:il2, :) = pm%tile%alvc(il1:il2, :)
!-        csfv%ALIC(il1:il2, :) = pm%tile%alic(il1:il2, :)
!-        csfv%PAMX(il1:il2, :) = pm%tile%lamx(il1:il2, :)
!-        csfv%PAMN(il1:il2, :) = pm%tile%lamn(il1:il2, :)
!-        csfv%CMAS(il1:il2, :) = pm%tile%cmas(il1:il2, :)
!-        csfv%ROOT(il1:il2, :) = pm%tile%root(il1:il2, :)
!-        csfv%RSMN(il1:il2, :) = pm%tile%rsmn(il1:il2, :)
!-        csfv%QA50(il1:il2, :) = pm%tile%qa50(il1:il2, :)
!-        csfv%VPDA(il1:il2, :) = pm%tile%vpda(il1:il2, :)
!-        csfv%VPDB(il1:il2, :) = pm%tile%vpdb(il1:il2, :)
!-        csfv%PSGA(il1:il2, :) = pm%tile%psga(il1:il2, :)
!-        csfv%PSGB(il1:il2, :) = pm%tile%psgb(il1:il2, :)
!-        csfv%DRN(il1:il2) = pm%tile%drn(il1:il2)
!-        csfv%SDEP(il1:il2) = pm%tile%sdep(il1:il2)
!-        csfv%XSLP(il1:il2) = pm%tile%xslp(il1:il2)
!-        DDGAT(il1:il2) = pm%tile%dd(il1:il2)
!-        MANNGAT(il1:il2) = pm%tile%mann(il1:il2)
!-        XDGAT(il1:il2) = pm%tile%grkf(il1:il2)
!-        KSGAT(il1:il2) = pm%tile%ks(il1:il2)
!-        csfv%ZSNL(il1:il2) = pm%tile%zsnl(il1:il2)
!-        csfv%ZPLG(il1:il2) = pm%tile%zplg(il1:il2)
!-        csfv%ZPLS(il1:il2) = pm%tile%zpls(il1:il2)
!-        if (allocated(RUNCLASS36_flgs%pm%tile%FREZTH)) then
!-            if (DIAGNOSEMODE) then
!-                call print_message('ICEBAL_FREEZE_THRESHOLD (FREZTH) override is ACTIVE.')
!-                if (allocated(RUNCLASS36_flgs%pm%constant%FREZTH)) then
!-                    write(line, FMT_GEN) 'Uniform value: ', RUNCLASS36_flgs%pm%constant%FREZTH
!-                    call print_message(line)
!-                end if
!-                if (allocated(RUNCLASS36_flgs%pm%gru%FREZTH)) then
!-                    write(line, FMT_GEN) 'GRU value: ', (RUNCLASS36_flgs%pm%gru%FREZTH(j), j = 1, NTYPE)
!-                    call print_message(line)
!-                end if
!-            end if
!-            cglv%FREZTH(il1:il2) = RUNCLASS36_flgs%pm%tile%FREZTH(il1:il2)
!-        else
!-            cglv%FREZTH(il1:il2) = -2.0
!-        end if
!-        if (allocated(RUNCLASS36_flgs%pm%tile%SWELIM)) then
!-            if (DIAGNOSEMODE) then
!-                call print_message('ICEBAL_SWE_LIMIT (SWELIM) override is ACTIVE.')
!-                if (allocated(RUNCLASS36_flgs%pm%constant%SWELIM)) then
!-                    write(line, FMT_GEN) 'Uniform value: ', RUNCLASS36_flgs%pm%constant%SWELIM
!-                    call print_message(line)
!-                end if
!-                if (allocated(RUNCLASS36_flgs%pm%gru%SWELIM)) then
!-                    write(line, FMT_GEN) 'GRU value: ', (RUNCLASS36_flgs%pm%gru%SWELIM(j), j = 1, NTYPE)
!-                    call print_message(line)
!-                end if
!-            end if
!-            cglv%SWELIM(il1:il2) = RUNCLASS36_flgs%pm%tile%SWELIM(il1:il2)
!-        else
!-            cglv%SWELIM(il1:il2) = 100.0
!-        end if
!-        if (allocated(RUNCLASS36_flgs%pm%tile%SNDENLIM)) then
!-            if (DIAGNOSEMODE) then
!-                call print_message('ICEBAL_SNOW_DENSITY_LIMIT (SNDENLIM) override is ACTIVE.')
!-                if (allocated(RUNCLASS36_flgs%pm%constant%SNDENLIM)) then
!-                    write(line, FMT_GEN) 'Uniform value: ', RUNCLASS36_flgs%pm%constant%SNDENLIM
!-                    call print_message(line)
!-                end if
!-                if (allocated(RUNCLASS36_flgs%pm%gru%SNDENLIM)) then
!-                    write(line, FMT_GEN) 'GRU value: ', (RUNCLASS36_flgs%pm%gru%SNDENLIM(j), j = 1, NTYPE)
!-                    call print_message(line)
!-                end if
!-            end if
!-            cglv%SNDENLIM(il1:il2) = RUNCLASS36_flgs%pm%tile%SNDENLIM(il1:il2)
!-        else
!-            cglv%SNDENLIM(il1:il2) = 900.0
!-        end if

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
        IDISP = NO_DATA_INT
        IZREF = NO_DATA_INT
        ISLFD = NO_DATA_INT
        IPCP = NO_DATA_INT
        ITC = NO_DATA_INT
        ITCG = NO_DATA_INT
        ITG = NO_DATA_INT
        csfv%IWF = NO_DATA_INT
        IPAI = NO_DATA_INT
        IHGT = NO_DATA_INT
        IALC = NO_DATA_INT
        IALS = NO_DATA_INT
        IALG = NO_DATA_INT
        ICTEMMOD = NO_DATA_INT
        SOILINIFLAG = NO_DATA_INT
        NUDGESDEPFLAG = NO_DATA_INT
        FROZENSOILINFILFLAG = NO_DATA_INT
        catv%ZRFM = NO_DATA_REAL
        catv%ZRFH = NO_DATA_REAL
        catv%ZBLD = NO_DATA_REAL
        catv%GC = NO_DATA_REAL
        csfv%FCAN = NO_DATA_REAL
        csfv%LNZ0 = NO_DATA_REAL
        csfv%ALVC = NO_DATA_REAL
        csfv%ALIC = NO_DATA_REAL
        csfv%PAMX = NO_DATA_REAL
        csfv%PAMN = NO_DATA_REAL
        csfv%CMAS = NO_DATA_REAL
        csfv%ROOT = NO_DATA_REAL
        csfv%RSMN = NO_DATA_REAL
        csfv%QA50 = NO_DATA_REAL
        csfv%VPDA = NO_DATA_REAL
        csfv%VPDB = NO_DATA_REAL
        csfv%PSGA = NO_DATA_REAL
        csfv%PSGB = NO_DATA_REAL
        csfv%DRN = NO_DATA_REAL
        csfv%SDEP = NO_DATA_REAL
        csfv%FARE = NO_DATA_REAL
        DDGAT = NO_DATA_REAL
        csfv%XSLP = NO_DATA_REAL
        XDGAT = NO_DATA_REAL
        MANNGAT = NO_DATA_REAL
        KSGAT = NO_DATA_REAL
        csfv%MID = NO_DATA_INT
        csfv%SAND = NO_DATA_REAL
        csfv%CLAY = NO_DATA_REAL
        csfv%ORGM = NO_DATA_REAL
        cpv%TBAR = NO_DATA_REAL
        cpv%TBAS = NO_DATA_REAL
        cpv%QAC = NO_DATA_REAL
        cpv%TAC = NO_DATA_REAL
        cpv%TCAN = NO_DATA_REAL
        cpv%TSNO = NO_DATA_REAL
        cpv%TPND = NO_DATA_REAL
        cpv%THLQ = NO_DATA_REAL
        cpv%THIC = NO_DATA_REAL
        cpv%TSFS = NO_DATA_REAL
        cpv%ZPND = NO_DATA_REAL
        cpv%RCAN = NO_DATA_REAL
        cpv%SNCAN = NO_DATA_REAL
        cpv%CMAI = NO_DATA_REAL
        cpv%SNO = NO_DATA_REAL
        cpv%ALBS = NO_DATA_REAL
        cpv%RHOS = NO_DATA_REAL
        cpv%WSNO = NO_DATA_REAL
        cpv%GRO = NO_DATA_REAL
        csfv%ZPLG = NO_DATA_REAL
        csfv%ZPLS = NO_DATA_REAL
        csfv%ZSNL = NO_DATA_REAL
        cglv%FREZTH = NO_DATA_REAL
        cglv%SWELIM = NO_DATA_REAL
        cglv%SNDENLIM = NO_DATA_REAL

        !> Allocate the filter for assigning distributed fields when MID (mosaic ID) is used.
        allocate(value_mask(NML))
        value_mask = .true.

        !> Assign input variables from field list.
        do i = 1, size(input_field_list)

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
                        call get_field_value(field, IDISP)
                    case ('IZREF')
                        call get_field_value(field, IZREF)
                    case ('ISLFD')
                        call get_field_value(field, ISLFD)
                    case ('IPCP')
                        call get_field_value(field, IPCP)
                    case ('ITC')
                        call get_field_value(field, ITC)
                    case ('ITCG')
                        call get_field_value(field, ITCG)
                    case ('ITG')
                        call get_field_value(field, ITG)
                    case ('IWF')
                        call get_field_value(field, csfv%IWF, DIM_NAME_G, value_mask = value_mask)
                    case ('IPAI')
                        call get_field_value(field, IPAI)
                    case ('IHGT')
                        call get_field_value(field, IHGT)
                    case ('IALC')
                        call get_field_value(field, IALC)
                    case ('IALS')
                        call get_field_value(field, IALS)
                    case ('IALG')
                        call get_field_value(field, IALG)
                    case ('ICTEMMOD')
                        call get_field_value(field, ICTEMMOD)
                    case ('SOILINIFLAG')
                        call get_field_value(field, SOILINIFLAG)
                    case ('NUDGESDEPFLAG')
                        call get_field_value(field, NUDGESDEPFLAG)
                    case ('FROZENSOILINFILFLAG')
                        call get_field_value(field, FROZENSOILINFILFLAG)
                    case ('ZREFMOMENTUM')
                        call get_field_value(field, catv%ZRFM, DIM_NAME_G, valid_min = 0.0, value_mask = value_mask)
                    case ('ZREFTHERMODYN')
                        call get_field_value(field, catv%ZRFH, DIM_NAME_G, valid_min = 0.0, value_mask = value_mask)
                    case ('ZREFBLEND')
                        call get_field_value(field, catv%ZBLD, DIM_NAME_G, valid_min = 0.0, value_mask = value_mask)
                    case ('GRIDCOVER')
                        call get_field_value( &
                            field, catv%GC, DIM_NAME_G, valid_min = -1.0, valid_max = -1.0, value_mask = value_mask)
                    case ('FCAN')
                        call get_field_value( &
                            field, csfv%FCAN, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('LNZ0')
                        call get_field_value( &
                            field, csfv%LNZ0, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('ALVC')
                        call get_field_value( &
                            field, csfv%ALVC, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('ALIC')
                        call get_field_value( &
                            field, csfv%ALIC, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('PAMX')
                        call get_field_value( &
                            field, csfv%PAMX, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('PAMN')
                        call get_field_value( &
                            field, csfv%PAMN, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('CMAS')
                        call get_field_value( &
                            field, csfv%CMAS, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('ROOT')
                        call get_field_value( &
                            field, csfv%ROOT, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('RSMN')
                        call get_field_value( &
                            field, csfv%RSMN, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('QA50')
                        call get_field_value( &
                            field, csfv%QA50, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('VPDA')
                        call get_field_value( &
                            field, csfv%VPDA, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('VPDB')
                        call get_field_value( &
                            field, csfv%VPDB, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('PSGA')
                        call get_field_value( &
                            field, csfv%PSGA, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('PSGB')
                        call get_field_value( &
                            field, csfv%PSGB, (/DIM_NAME_G, DIM_NAME_L/), level_map = CANOPY_LEVEL_MAP, value_mask = value_mask)
                    case ('DRN')
                        call get_field_value(field, csfv%DRN, DIM_NAME_G, valid_min = 0.0, valid_max = 1.0, value_mask = value_mask)
                    case ('SDEP', 'SDEPTH')
                        call get_field_value(field, csfv%SDEP, DIM_NAME_G, valid_min = 0.0, value_mask = value_mask)
                    case ('FARE')
                        call get_field_value( &
                            field, csfv%FARE, DIM_NAME_G, valid_min = 1.0, valid_max = 1.0, value_mask = value_mask)
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
                        call get_field_value(field, csfv%XSLP, DIM_NAME_G, value_mask = value_mask)
                    case ('GRKEF')
                        call get_field_value(field, XDGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('MANN')
                        call get_field_value(field, MANNGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('KSAT')
                        call get_field_value(field, KSGAT, DIM_NAME_G, value_mask = value_mask)
                    case ('MOSAICID', 'MID')
                        call get_field_value(field, csfv%MID, DIM_NAME_G, valid_min = 1)

                        !> Check for mosaic ID to control assigning distributed parameters.
                        if (any(csfv%MID >= 100 .and. csfv%MID < 1000)) then
                            call print_remark( &
                                "There are tile values of variable '" // trim(field%label) // "' > 100 and < 1000. " // &
                                "Those tiles will not be reassigned new values if the parameter inputs include " // &
                                "spatially distributed fields.")
                            where (csfv%MID >= 100 .and. csfv%MID < 1000) value_mask = .false.
                        end if
                    case ('SOLSAND', 'SAND')
                        call get_field_value(field, csfv%SAND, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('SOLCLAY', 'CLAY')
                        call get_field_value(field, csfv%CLAY, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('SOLORGCONTENT', 'ORGM')
                        call get_field_value(field, csfv%ORGM, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('TSOL', 'TBAR', 'TILE_TBAR')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value( &
                            field, cpv%TBAR, (/DIM_NAME_G, DIM_NAME_L/), const_add = factor_a, value_mask = value_mask)
                    case ('TBAS', 'TILE_TBAS')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value(field, cpv%TBAS, DIM_NAME_G, const_add = factor_a, value_mask = value_mask)
                    case ('QAC', 'TILE_QAC')
                        call get_field_value(field, cpv%QAC, DIM_NAME_G, value_mask = value_mask)
                    case ('TAC', 'TILE_TAC')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value(field, cpv%TAC, DIM_NAME_G, const_add = factor_a, value_mask = value_mask)
                    case ('TCAN', 'TILE_TCAN')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value(field, cpv%TCAN, DIM_NAME_G, const_add = factor_a, value_mask = value_mask)
                    case ('TSNO', 'TILE_TSNO')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value(field, cpv%TSNO, DIM_NAME_G, const_add = factor_a, value_mask = value_mask)
                    case ('TPND', 'TILE_TPND')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value(field, cpv%TPND, DIM_NAME_G, const_add = factor_a, value_mask = value_mask)
                    case ('THLQSOL', 'THLQ', 'TILE_THLQ')
                        call get_field_value(field, cpv%THLQ, (/DIM_NAME_G, DIM_NAME_L/), valid_min = 0.04, value_mask = value_mask)
                    case ('THICSOL', 'THIC', 'TILE_THIC')
                        call get_field_value(field, cpv%THIC, (/DIM_NAME_G, DIM_NAME_L/), value_mask = value_mask)
                    case ('TSFS', 'TILE_TSFS')
                        if (index(to_lowercase(field%meta%units), 'c') > 0) then

                            !> Unit conversion to Kelvin from degrees Celsius.
                            call print_remark( &
                                "Converting the units of the variable '" // trim(field%label) // &
                                "' from '" // trim(field%meta%units) // "' to 'Kelvin'.")
                            factor_a = TFREZ
                        end if
                        call get_field_value( &
                            field, cpv%TSFS, (/DIM_NAME_G, DIM_NAME_L/), const_add = factor_a, value_mask = value_mask)
                    case ('ZPND', 'TILE_ZPND')
                        call get_field_value(field, cpv%ZPND, DIM_NAME_G, value_mask = value_mask)
                    case ('LQWSCAN', 'RCAN', 'TILE_RCAN')
                        call get_field_value(field, cpv%RCAN, DIM_NAME_G, value_mask = value_mask)
                    case ('FZWSCAN', 'SCAN', 'TILE_SNCAN')
                        call get_field_value(field, cpv%SNCAN, DIM_NAME_G, value_mask = value_mask)
                    case ('CMAI', 'TILE_CMAS')
                        call get_field_value(field, cpv%CMAI, DIM_NAME_G, value_mask = value_mask)
                    case ('SNO', 'TILE_SNO')
                        call get_field_value(field, cpv%SNO, DIM_NAME_G, value_mask = value_mask)
                    case ('ALBSNO', 'ALBS', 'TILE_ALBS')
                        call get_field_value(field, cpv%ALBS, DIM_NAME_G, value_mask = value_mask)
                    case ('RHOSNO', 'RHOS', 'TILE_RHOS')
                        call get_field_value(field, cpv%RHOS, DIM_NAME_G, value_mask = value_mask)
                    case ('LQWSSNO', 'WSNO', 'TILE_WSNO')
                        call get_field_value(field, cpv%WSNO, DIM_NAME_G, value_mask = value_mask)
                    case ('GRO', 'TILE_GRO')
                        call get_field_value(field, cpv%GRO, DIM_NAME_G, value_mask = value_mask)
                    case ('ZPLG')
                        call get_field_value(field, csfv%ZPLG, DIM_NAME_G, value_mask = value_mask)
                    case ('ZPLS')
                        call get_field_value(field, csfv%ZPLS, DIM_NAME_G, value_mask = value_mask)
                    case ('ZSNL')
                        call get_field_value(field, csfv%ZSNL, DIM_NAME_G, value_mask = value_mask)
                    case ('FREZTH', 'ICEBAL_FREEZE_THRESHOLD')
                        call get_field_value(field, cglv%FREZTH, DIM_NAME_G)
                    case ('SWELIM', 'ICEBAL_SWE_LIMIT', 'SNDEPLIM', 'ICEBAL_SNOW_DEPTH_LIMIT')
                        call get_field_value(field, cglv%SWELIM, DIM_NAME_G)
                    case ('SNDENLIM', 'ICEBAL_SNOW_DENSITY_LIMIT')
                        call get_field_value(field, cglv%SNDENLIM, DIM_NAME_G)
                end select
            end associate
        end do

        !> Assign defaults or derived fields where possible.
        ierr = 0
        z = 0

        !> Options.
        call check_field_assigned(IDISP, 'IDISP', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(IZREF, 'IZREF', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ISLFD, 'ISLFD', 2, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(IPCP, 'IPCP', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ITC, 'ITC', 2, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ITCG, 'ITCG', 2, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ITG, 'ITG', 2, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%IWF, 'IWF', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(IPAI, 'IPAI', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(IHGT, 'IHGT', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(IALC, 'IALC', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(IALS, 'IALS', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(IALG, 'IALG', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(ICTEMMOD, 'ICTEMMOD', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(SOILINIFLAG, 'SOILINIFLAG', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(NUDGESDEPFLAG, 'NUDGESDEPFLAG', 0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(FROZENSOILINFILFLAG, 'FROZENSOILINFILFLAG', 0, error_status = z); if (z /= 0) ierr = z

        !> Constants.
        call check_field_assigned(catv%ZBLD, 'ZBLD', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(catv%ZRFH, 'ZRFH', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(catv%ZRFM, 'ZRFM', is_error = .true., error_status = z); if (z /= 0) ierr = z
        if (any(catv%ZRFH /= catv%ZRFM) .and. ISLFD /= 2) then
            write(line, *) ISLFD
            call print_warning( &
                "The momentum 'ZRFM' and thermodynamic 'ZRFH' reference levels are staggered while 'ISLFD' is set to " // &
                trim(adjustl(line)) // ", which enables the 'DRCOEF' routine that does not support this configuration. " // &
                "When running CLASS with staggered levels, 'ISLFD' should be set to 2 to use the 'FLXSURFZ' routine instead.")
        end if
!+        if (ZRFM < canopy height or ZRFH < canopy height): may cause model to crash
        if (any(catv%ZRFM < catv%ZBLD) .or. any(catv%ZRFH < catv%ZBLD)) then
            call print_warning( &
                "The momentum 'ZRFM' and thermodynamic 'ZRFH' reference levels are below the blending height 'ZBLD', " // &
                "which may lead to instabilities when averaging roughness lengths.")
        end if
        call check_field_assigned(catv%GC, 'GC', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%FCAN, 'FCAN', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%LNZ0, 'LNZ0', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%ALVC, 'ALVC', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%ALIC, 'ALIC', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%PAMX, 'PAMX', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%PAMN, 'PAMN', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%CMAS, 'CMAS', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%ROOT, 'ROOT', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%RSMN, 'RSMN', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%QA50, 'QA50', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%VPDA, 'VPDA', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%VPDB, 'VPDB', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%PSGA, 'PSGA', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%PSGB, 'PSGB', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%DRN, 'DRN', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%SDEP, 'SDEP', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else if (NUDGESDEPFLAG == 1) then

            !> Nudge SDEP to the nearest interface between soil layers in the profile.
            call print_remark( &
                "'NUDGESDEPFLAG' is enabled. The permeable depth 'SDEP' will be nudged at each location to the nearest " // &
                "interface between soil layers in the soil profile.")
            where (csfv%SDEP <= shd%lc%sl%zbot(1))

                !> Within the first layer, nudge to the first boundary and not to zero so it remains active.
                csfv%SDEP = shd%lc%sl%zbot(1)
            elsewhere (csfv%SDEP > shd%lc%sl%zbot(NSL))

                !> Limit depths deeper than the soil profile to the last boundary for consistency.
                csfv%SDEP = shd%lc%sl%zbot(NSL)
            end where
            do j = 2, NSL

                !> Nudge in-between depths to the closest boundary.
                where (csfv%SDEP > shd%lc%sl%zbot(j - 1) .and. csfv%SDEP <= (0.5*shd%lc%sl%delz(j) + shd%lc%sl%zbot(j - 1)))
                    csfv%SDEP = shd%lc%sl%zbot(j - 1)
                elsewhere (pm%tile%sdep > (0.5*shd%lc%sl%delz(j) + shd%lc%sl%zbot(j - 1)) .and. csfv%SDEP <= shd%lc%sl%zbot(j))
                    csfv%SDEP = shd%lc%sl%zbot(j)
                end where
            end do
        end if
        call check_field_assigned(csfv%FARE, 'FARE', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(DDGAT, 'DD', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%XSLP, 'XSLP', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(XDGAT, 'GRKEF', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(MANNGAT, 'MANN', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(KSGAT, 'KSAT', is_error = .true., error_status = z); if (z /= 0) ierr = z
        if (any(csfv%DRN < 1.0 .and. csfv%SDEP <= 4.1)) then
            call print_remark( &
                "The configuration contains drainage index 'DRN' values less than 1.0. " // &
                "In most cases, 'DRN' can be set to 1.0 except in cases of deep soils where it is desired to " // &
                "suppress drainage from the bottom of the soil profile or in soils with a high water table.")
        end if
        call check_field_assigned(csfv%MID, 'MOSAICID', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%SAND, 'SAND', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else
            if (any(csfv%SDEP == 0.0 .and. csfv%SAND(:, 1) > -2.5)) then

                !> Check the first layer for impermable soils not identified as such.
                write(line, *) count(csfv%SDEP == 0.0 .and. csfv%SAND(:, 1) > -2.5)
                call print_remark( &
                    "The configuration contains " // trim(adjustl(line)) // " tiles where the permeable depth 'SDEP' is zero. " // &
                    "The soil properties in those areas will be set to those of 'exposed rock' ('SAND' -3).")
                do j = 1, NSL
                    where (csfv%SDEP == 0.0 .and. csfv%SAND(:, 1) > -2.5) csfv%SAND(:, j) = -3.0
                end do
            end if

            !> Check for unidentified impermeable soils at depth.
            do j = 2, NSL
                where (csfv%SDEP < (shd%lc%sl%ZBOT(j - 1) + 0.001) .and. csfv%SAND(:, j) > -2.5) csfv%SAND(:, j) = -3.0
            end do
        end if
        call check_field_assigned(csfv%CLAY, 'CLAY', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else if (any(csfv%SAND <= -2)) then
            call print_remark( &
                "The configuration contains areas of non-mineral soils ('SAND' -2, -3 or -4). " // &
                "The 'CLAY' fraction will be set to zero in those areas.")
            where (csfv%SAND <= -2) csfv%CLAY = 0.0
        end if
        call check_field_assigned(csfv%ORGM, 'ORGM', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else
            if (any(csfv%SAND <= -2)) then
                call print_remark( &
                    "The configuration contains areas of non-mineral soils ('SAND' -2, -3 or -4). " // &
                    "The 'ORGM' fraction will be set to zero in those areas.")
                where (csfv%SAND <= -2) csfv%ORGM = 0.0
            end if
            if (any((csfv%ORGM /= 1.0 .or. csfv%ORGM /= 2.0 .or. csfv%ORGM /= 3.0) .and. csfv%SAND == -2)) then
                call print_error("The configuration contains areas of organic soils ('SAND' -2) but improper 'ORGM' class types.")
                ierr = 1
            end if
        end if
        call check_field_assigned(csfv%ZPLG, 'ZPLG', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%ZPLS, 'ZPLS', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(csfv%ZSNL, 'ZSNL', is_error = .true., error_status = z); if (z /= 0) ierr = z
        if (any(csfv%SAND == -4)) then
            call check_field_assigned(cglv%FREZTH, 'FREZTH', -2.0, error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(cglv%SWELIM, 'SWELIM', 100.0, error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(cglv%SNDENLIM, 'SNDENLIM', 900.0, error_status = z); if (z /= 0) ierr = z
        else
            if (any(cglv%FREZTH /= NO_DATA_REAL)) then
                call print_remark( &
                    "The variable '" // 'FREZTH' // "' is assigned but not active (no glaciers are active).")
            end if
            if (any(cglv%SWELIM /= NO_DATA_REAL)) then
                call print_remark( &
                    "The variable '" // 'SWELIM' // "' is assigned but not active (no glaciers are active).")
            end if
            if (any(cglv%SNDENLIM /= NO_DATA_REAL)) then
                call print_remark( &
                    "The variable '" // 'SNDENLIM' // "' is assigned but not active (no glaciers are active).")
            end if
        end if

        !> States.
        call check_field_assigned(cpv%TCAN, 'TCAN', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else if (any(cpv%TCAN < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TCAN' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (cpv%TCAN < (TFREZ - 100.0)) cpv%TCAN = cpv%TCAN + TFREZ
        end if
        call check_field_assigned(cpv%TAC, 'TAC', error_status = z)
        if (z /= 0 .and. all(cpv%TCAN /= NO_DATA_REAL)) then
!+ Set unassigned TAC to air temperature

            !> Set the field to 'TCAN'.
            call print_remark("Setting the unassigned values of 'TAC' from 'TCAN'.")
            where (cpv%TAC == NO_DATA_REAL) cpv%TAC = cpv%TCAN
        else if (any(cpv%TAC < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TAC' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (cpv%TAC < (TFREZ - 100.0)) cpv%TAC = cpv%TAC + TFREZ
        end if
        call check_field_assigned(cpv%QAC, 'QAC', 0.5e-2, error_status = z); if (z /= 0) ierr = z
!+ Set unassigned QAC to air specific humidity
        call check_field_assigned(cpv%RCAN, 'RCAN', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(cpv%SNCAN, 'SCAN', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(cpv%CMAI, 'RCAN', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(cpv%CMAI, 'CMAI', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(cpv%GRO, 'GRO', 0.0, error_status = z); if (z /= 0) ierr = z
!+        if (GRO > 0): GRO > 1 during growing season, otherwise zero
        call check_field_assigned(cpv%ZPND, 'ZPND', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(cpv%TPND, 'TPND', 0.0, error_status = z)
        if (any(cpv%TPND < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TPND' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (cpv%TPND < (TFREZ - 100.0)) cpv%TPND = cpv%TPND + TFREZ
        end if
        call check_field_assigned(cpv%SNO, 'SNO', 0.0, error_status = z)
        if (any(cpv%SNO > 0.0)) then
            call print_remark("Starting the simulation in snow-covered conditions could lead to a bias in simulated states.")
            call check_field_assigned(cpv%ALBS, 'ALBS', is_error = .true., error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(cpv%RHOS, 'RHOS', is_error = .true., error_status = z); if (z /= 0) ierr = z
            call check_field_assigned(cpv%TSNO, 'TSNO', is_error = .true., error_status = z)
            if (z /= 0) then
                ierr = z
            else if (any(cpv%TSNO < (TFREZ - 100.0))) then

                !> Implicit conversion.
                call print_remark( &
                    "There may be tile values of the variable '" // 'TSNO' // &
                    "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
                where (cpv%TSNO < (TFREZ - 100.0)) cpv%TSNO = cpv%TSNO + TFREZ
            end if
            call check_field_assigned(cpv%WSNO, 'WSNO', 0.0, error_status = z); if (z /= 0) ierr = z
        else
            cpv%ALBS = 0.0
            cpv%RHOS = 0.0
            cpv%TSNO = TFREZ
            cpv%WSNO = 0.0
        end if
        call check_field_assigned(cpv%TBAR, 'TBAR', is_error = .true., error_status = z)
        if (z /= 0) then
            ierr = z
        else if (any(cpv%TBAR < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TBAR' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (cpv%TBAR < (TFREZ - 100.0)) cpv%TBAR = cpv%TBAR + TFREZ
        end if
        call check_field_assigned(cpv%THIC, 'THIC', 0.0, error_status = z)
        if (z == 0) then
            if (any(csfv%SAND <= -3)) then
                call print_remark( &
                    "Setting the volumetric frozen water content 'THIC' of non-soil layers ('SAND' -3 or -4) to zero.")
                where (csfv%SAND <= -3) cpv%THIC = 0.0
            end if
        end if
        call check_field_assigned(cpv%THLQ, 'THLQ', 0.04, error_status = z)
        if (z == 0) then
            if (any(csfv%SAND <= -3)) then
                call print_remark( &
                    "Setting the volumetric liquid water content 'THLQ' of non-soil layers ('SAND' -3 or -4) to zero.")
                where (csfv%SAND <= -3) cpv%THLQ = 0.0
            end if
            if (any(csfv%SAND >= -2 .and. cpv%THLQ < 0.04)) then
                write(line, *) count(csfv%SAND >= -2 .and. cpv%THLQ < 0.04)
                call print_warning( &
                    "The configuration contains " // trim(adjustl(line)) // " tiles where the initial value for the " // &
                    "volumetric liquid water content 'THLQ' is lower than the permissible limit of 0.04.")
            end if
        end if
        call check_field_assigned(cpv%TBAS, 'TBAS', error_status = z)
        if (z /= 0 .and. all(cpv%TBAR /= NO_DATA_REAL)) then
            call print_remark("Assigning the unassigned values of 'TBAS' from the bottom layer value of 'TBAR'.")
            where (cpv%TBAS == NO_DATA_REAL) cpv%TBAS = cpv%TBAR(:, NSL)
        else if (any(cpv%TBAS < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TBAS' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (cpv%TBAS < (TFREZ - 100.0)) cpv%TBAS = cpv%TBAS + TFREZ
        end if
        call check_field_assigned(cpv%TSFS, 'TSFS', error_status = z)
        if (z /= 0 .and. all(cpv%TBAR /= NO_DATA_REAL)) then
            call print_remark("Assigning the unassigned values of 'TSFS' to a background field and from 'TBAR'.")
            do j = 1, 2

                !> Assign snow-covers the freezing point of water.
                where (cpv%TSFS(:, j) == NO_DATA_REAL) cpv%TSFS(:, j) = TFREZ
            end do
            do j = 3, 4

                !> Assign ground covers to the first layer 'TBAR'.
                where (cpv%TSFS(:, j) == NO_DATA_REAL) cpv%TSFS(:, j) = cpv%TBAR(:, 1)
            end do
        else if (any(cpv%TSFS < (TFREZ - 100.0))) then

            !> Implicit conversion.
            call print_remark( &
                "There may be tile values of the variable '" // 'TSFS' // &
                "' in units of 'deg C'. Applying an implicit conversion of these values to units of 'Kelvin'.")
            where (cpv%TSFS < (TFREZ - 100.0)) cpv%TSFS = cpv%TSFS + TFREZ
        end if

        !> Options.
!        IDISP 0: Indistinguishable vegetation displacement height from terrain.
!        IDISP 1: Vegetation displament height is calculated.
!        IZREF 2: Atmospheric bottom is at the local surface roughness length (when coupled with numerical weather outputs).
!        IZREF 1: Atmospheric bottom is not at the local surface roughness length.
!        ISLFD 2: RPN surface layer flux and screen-level routines or if momentum and thermodynamic reference levels are staggered (FLXSURFZ, DIASURFZ).
!        ISLFD 0: CCCma surface layer flux and screen-level routines (DRCOEF).
!        ISLFD 1: CCCma surface layer flux and RPN screen-level routines (DRCOEF, SLDIAG).
!        IPCP 1: Total precipitation rate is partitioned by the freezing point of water.
!        IPCP 2: Total precipitation rate is partitioned by a linear phase change between 0 and 2 degrees Celsius.
!        IPCP 3: Total precipitation rate is partitioned using a polynomial curve between 0 and 6 degrees Celsius.
!        IPCP 4: Rain and snow rate are provided separately.
!        ITC 2: Newton-Raphson iteration for canopy cover (for RPN surface layer flux routines, 'ISLFD 2').
!        ITC 1: Bisection iteration for canopy cover (for CCCma surface layer flux routines, 'ISLFD 0' or 'ISLFD 1').
!        ITCG 2: Newton-Raphson iteration for ground cover under the canopy (for RPN surface layer flux routines, 'ISLFD 2').
!        ITCG 1: Bisection iteration for ground cover under the canopy (for CCCma surface layer flux routines, 'ISLFD 0' or 'ISLFD 1').
!        ITG 2: Newton-Raphson iteration for bare ground cover (for RPN surface layer flux routines, 'ISLFD 2').
!        ITG 1: Bisection iteration for bare ground cover (for CCCma surface layer flux routines, 'ISLFD 0' or 'ISLFD 1').
!        IWF 0: CLASS proper (hard-coded ponding limits, no interflow routines).
!        IWF 1: Sloped overland and interflow routines (WATROF, WATDRN).
!        IWF 2: Prairie pothole statistical parameterization for overland flow (PDMROF) with no interflow routines.
!        IWF 3: Sloped prairie pothole statistical parameterization for overland flow (WATROF, PDMROF; 'LATFLOW') with sloped interflow routines (WATDRN).
!        IPAI 0: No user-provided plant area index per time-step (option not available).
!        IHGT 0: No user-provided vegetation height per time-step (option not available).
!        IALC 0: No user-provided canopy albedo per CLASS time-step (option not available).
!        IALS 0: No user-provided snow albedo per CLASS time-step (option not available).
!        IALG 0: No user-provided ground albedo per CLASS time-step (option not available).
!        ICTEMMOD 0: No coupling with the Canadian Terrestrial Ecosystem Model (CTEM) (option not available).

!        ITC/ITCG/ITG should be 1 if ISLFD is 1
!        ITC/ITCG/ITG should be 2 if ISLFD is 2

!        if no FCLO, 1 when PRE occuring, otherwise a fraction of incoming diffuse radiation between 0.1 and 1.0.
!        if no FSVH/FSIH, 50:50 partitioning between incoming and near-infrared radiation.

        !> Old variables.
        vs%tile%cmas(il1:il2) = cpv%CMAI(il1:il2)
        vs%tile%lqwssno(il1:il2) = cpv%WSNO(il1:il2)
        vs%tile%qacan(il1:il2) = cpv%QAC(il1:il2)
        vs%tile%tcan(il1:il2) = cpv%TCAN(il1:il2)
        vs%tile%tacan(il1:il2) = cpv%TAC(il1:il2)
        vs%tile%tsno(il1:il2) = cpv%TSNO(il1:il2)
        vs%tile%tpnd(il1:il2) = cpv%TPND(il1:il2)
        vs%tile%zpnd(il1:il2) = cpv%ZPND(il1:il2)
        vs%tile%lqwscan(il1:il2) = cpv%RCAN(il1:il2)
        vs%tile%fzwscan(il1:il2) = cpv%SNCAN(il1:il2)
        vs%tile%sno(il1:il2) = cpv%SNO(il1:il2)
        vs%tile%albsno(il1:il2) = cpv%ALBS(il1:il2)
        vs%tile%rhosno(il1:il2) = cpv%RHOS(il1:il2)
        vs%tile%gro(il1:il2) = cpv%GRO(il1:il2)
        vs%tile%tsfs(il1:il2, :) = cpv%TSFS(il1:il2, :)
        vs%tile%tsol(il1:il2, :) = cpv%TBAR(il1:il2, :)
        vs%tile%thlqsol(il1:il2, :) = cpv%THLQ(il1:il2, :)
        vs%tile%thicsol(il1:il2, :) = cpv%THIC(il1:il2, :)
        vs%tile%tbas(il1:il2) = cpv%TBAS(il1:il2)

        !> Derived or other fields.
        cdv%ITCT = 0

        !> FROZENSOILINFILFLAG 1.
        allocate(FRZCGAT(NML))
        allocate(INFILTYPE(NML), SI(NML), TSI(NML), &
                 SNOWMELTD(NML), SNOWMELTD_LAST(NML), SNOWINFIL(NML), &
                 CUMSNOWINFILCS(NML), MELTRUNOFF(NML), CUMSNOWINFILGS(NML))
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
                FRZCGAT(k) = hp%FRZCROW(shd%lc%ILMOS(k), shd%lc%JLMOS(k))
            end do
        end if

        !> *********************************************************************
        !> Call CLASSBG to set more CLASS variables
        !> *********************************************************************
        !> bjd - July 25, 2005: For inputting field measured soil properties.
        call CLASSBG(csfv%THP, csfv%THR, csfv%THM, csfv%BI, csfv%PSIS, csfv%GRKS, &
                     csfv%THRA, csfv%HCPS, csfv%TCS, csfv%THFC, csfv%PSIW, &
                     csfv%DELZW, csfv%ZBTW, csfv%ALGW, csfv%ALGD, &
                     csfv%SAND, csfv%CLAY, csfv%ORGM, shd%lc%sl%DELZ, shd%lc%sl%ZBOT, csfv%SDEP, &
                     csfv%ISND, csfv%IGDR, NML, il1, il2, NSL, ICTEMMOD, &
                     pm%gru%thpor, pm%gru%thlret, pm%gru%thlmin, pm%gru%bi, &
                     pm%gru%psisat, pm%gru%grksat, pm%gru%hcps, pm%gru%tcs, &
                     NA, NTYPE, shd%lc%ILG, shd%lc%ILMOS, shd%lc%JLMOS)
        if (any(csfv%DELZW(:, 2:) < 0.25)) then
            call print_warning( &
                "The configuration contains thin layers at depth which may cause instability in the temperature stepping scheme.")
        end if

        pm%tile%alwet(il1:il2) = csfv%ALGW(il1:il2)
        pm%tile%aldry(il1:il2) = csfv%ALGD(il1:il2)
        pm%tile%thpor(il1:il2, :) = csfv%THP(il1:il2, :)
        pm%tile%thlret(il1:il2, :) = csfv%THR(il1:il2, :)
        pm%tile%thlmin(il1:il2, :) = csfv%THM(il1:il2, :)
        pm%tile%bi(il1:il2, :) = csfv%BI(il1:il2, :)
        pm%tile%psisat(il1:il2, :) = csfv%PSIS(il1:il2, :)
        pm%tile%grksat(il1:il2, :) = csfv%GRKS(il1:il2, :)
        pm%tile%thlrat(il1:il2, :) = csfv%THRA(il1:il2, :)
        pm%tile%hcps(il1:il2, :) = csfv%HCPS(il1:il2, :)
        pm%tile%tcs(il1:il2, :) = csfv%TCS(il1:il2, :)
        pm%tile%thfc(il1:il2, :) = csfv%THFC(il1:il2, :)
        pm%tile%psiwlt(il1:il2, :) = csfv%PSIW(il1:il2, :)
        vs%tile%dzsolhyd(il1:il2, :) = csfv%DELZW(il1:il2, :)
        vs%tile%zsolhyd(il1:il2) = sum(csfv%ZBTW(il1:il2, :), 2)

        !> CLASS output files.
        if (WF_NUM_POINTS > 0) call CLASSOUT_open_files(shd)

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
                catv%RADJ(i) = vs%grid%lat(vs%tile%from_cell(i))*PI/180.0
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
                catv%DLON(i) = vs%grid%lon(vs%tile%from_cell(i))
            end do
        end if

!-        do k = il1, il2
!-            ik = shd%lc%ILMOS(k)
!-            catv%RADJ(k) = shd%ylat(ik)*PI/180.0
!-            catv%DLON(k) = shd%xlng(ik)
!-        end do
        catv%Z0OR = 0.0
        catv%GGEO(:) = GGEOGRD(1)
        catv%ZDM = 10.0
        catv%ZDH = 2.0

        !> Check for errors.
        if (ierr /= 0) then
            call reset_tab()
            call print_error("Variables required by the module are not assigned or were missing during initialization.")
            call program_abort()
        end if

        !> IWF 2 (PDMROF) and IWF 3 (LATFLOW).
        allocate(&
            CMINPDM(NML), CMAXPDM(NML), BPDM(NML), K1PDM(NML), K2PDM(NML), &
            ZPNDPRECS(NML), ZPONDPREC(NML), ZPONDPREG(NML), ZPNDPREGS(NML), &
            UM1CS(NML), UM1C(NML), UM1G(NML), UM1GS(NML), &
            QM1CS(NML), QM1C(NML), QM1G(NML), QM1GS(NML), &
            QM2CS(NML), QM2C(NML), QM2G(NML), QM2GS(NML), &
            UMQ(NML), &
            FSTRCS(NML), FSTRC(NML), FSTRG(NML), FSTRGS(NML))
        ZPNDPRECS = 0.0; ZPONDPREC = 0.0; ZPONDPREG = 0.0; ZPNDPREGS = 0.0
        UM1CS = 0.0; UM1C = 0.0; UM1G = 0.0; UM1GS = 0.0
        QM1CS = 0.0; QM1C = 0.0; QM1G = 0.0; QM1GS = 0.0
        QM2CS = 0.0; QM2C = 0.0; QM2G = 0.0; QM2GS = 0.0
        UMQ = 0.0
        FSTRCS = 0.0; FSTRC = 0.0; FSTRG = 0.0; FSTRGS = 0.0
        if (any(csfv%IWF == 2) .or. any(csfv%IWF == 3)) then
!-            do k = il1, il2
!-                ik = shd%lc%ILMOS(k)
!-                jk = shd%lc%JLMOS(k)
!-                CMINPDM(k) = hp%CMINROW(ik, jk)
!-                CMAXPDM(k) = hp%CMAXROW(ik, jk)
!-                BPDM(k) = hp%BROW(ik, jk)
!-                K1PDM(k) = hp%K1ROW(ik, jk)
!-                K2PDM(k) = hp%K2ROW(ik, jk)
!-            end do
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

    end subroutine

    !> Description:
    !>  Subroutine to create a field list to read the binary sequential file for CLASS states created by 'SAVERESUMEFLAG 3'.
    subroutine read_init_prog_variables_class_row(fls, shd)

!-        use mpi_module
        use model_files_variables
        use sa_mesh_common
!-        use model_dates
        use input_fields
        use mesh_legacy_io

        type(fl_ids) fls
        type(ShedGridParams) shd

        !> For SAVERESUMEFLAG 3
!-        real(kind = 4), dimension(:, :), allocatable :: ALBSROW, CMAIROW, GROROW, QACROW, RCANROW, &
!-            RHOSROW, SCANROW, SNOROW, TACROW, TBASROW, &
!-            TCANROW, TPNDROW, TSNOROW, WSNOROW, ZPNDROW
!-        real(kind = 4), dimension(:, :, :), allocatable :: TBARROW, THICROW, THLQROW, TSFSROW

!-        integer NA, NTYPE, NSL, k, ik, jk, iun, ierr

        !> Local variables.
        integer file_id, error_status
        character(len = SHORT_FIELD_LENGTH) dim_names_cell_gru(2), dim_names_cell_gru_surf(3), dim_names_cell_gru_sol(3)

        !> Return if the process is not active.
        if (.not. RUNCLASS36_flgs%PROCESS_ACTIVE) return

        !> Local indices.
!-        NA = shd%NA
!-        NTYPE = shd%lc%NTYPE
!-        NSL = shd%lc%IGND

        !> Open the resume state file.
!-        iun = fls%fl(mfk%f883)%iun
!-        open(iun, file = trim(adjustl(fls%fl(mfk%f883)%fn)), status = 'old', action = 'read', &
!-             form = 'unformatted', access = 'sequential', iostat = ierr)
!todo: condition for ierr.

        !> Allocate temporary variables.
!-        allocate(ALBSROW(NA, NTYPE), CMAIROW(NA, NTYPE), GROROW(NA, NTYPE), QACROW(NA, NTYPE), RCANROW(NA, NTYPE), &
!-                 RHOSROW(NA, NTYPE), SCANROW(NA, NTYPE), SNOROW(NA, NTYPE), TACROW(NA, NTYPE), TBASROW(NA, NTYPE), &
!-                 TCANROW(NA, NTYPE), TPNDROW(NA, NTYPE), TSNOROW(NA, NTYPE), WSNOROW(NA, NTYPE), ZPNDROW(NA, NTYPE), &
!-                 TBARROW(NA, NTYPE, NSL), THICROW(NA, NTYPE, NSL), THLQROW(NA, NTYPE, NSL), TSFSROW(NA, NTYPE, 4))

        !> Read inital values from the file.
!-        read(iun) ALBSROW
!-        read(iun) CMAIROW
!-        read(iun) GROROW
!-        read(iun) QACROW
!-        read(iun) RCANROW
!-        read(iun) RHOSROW
!-        read(iun) SCANROW
!-        read(iun) SNOROW
!-        read(iun) TACROW
!-        read(iun) TBARROW
!-        read(iun) TBASROW
!-        read(iun) TCANROW
!-        read(iun) THICROW
!-        read(iun) THLQROW
!-        read(iun) TPNDROW
!-        read(iun) TSFSROW
!-        read(iun) TSNOROW
!-        read(iun) WSNOROW
!-        read(iun) ZPNDROW

        !> Close the file to free the unit.
!-        close(iun)

        !> Scatter the temporary variables.
!-        do k = il1, il2

            !> Grab the grid and GRU of the current tile.
!-            ik = shd%lc%ILMOS(k)
!-            jk = shd%lc%JLMOS(k)

            !> Assign values.
!-            vs%tile%albsno(k) = ALBSROW(ik, jk)
!-            vs%tile%cmas(k) = CMAIROW(ik, jk)
!-            vs%tile%gro(k) = GROROW(ik, jk)
!-            vs%tile%qacan(k) = QACROW(ik, jk)
!-            vs%tile%lqwscan(k) = RCANROW(ik, jk)
!-            vs%tile%rhosno(k) = RHOSROW(ik, jk)
!-            vs%tile%fzwscan(k) = SCANROW(ik, jk)
!-            vs%tile%sno(k) = SNOROW(ik, jk)
!-            vs%tile%tacan(k) = TACROW(ik, jk)
!-            vs%tile%tsol(k, :) = TBARROW(ik, jk, :)
!-            vs%tile%tbas(k) = TBASROW(ik, jk)
!-            vs%tile%tcan(k) = TCANROW(ik, jk)
!-            vs%tile%thicsol(k, :) = THICROW(ik, jk, :)
!-            vs%tile%thlqsol(k, :) = THLQROW(ik, jk, :)
!-            vs%tile%tpnd(k) = TPNDROW(ik, jk)
!-            vs%tile%tsfs(k, :) = TSFSROW(ik, jk, :)
!-            vs%tile%tsno(k) = TSNOROW(ik, jk)
!-            vs%tile%lqwssno(k) = WSNOROW(ik, jk)
!-            vs%tile%zpnd(k) = ZPNDROW(ik, jk)

!-        end do

        !> Deallocate temporary variables.
!-        deallocate(ALBSROW, CMAIROW, GROROW, QACROW, RCANROW, &
!-                   RHOSROW, SCANROW, SNOROW, TACROW, TBASROW, &
!-                   TCANROW, TPNDROW, TSNOROW, WSNOROW, ZPNDROW, &
!-                   TBARROW, THICROW, THLQROW, TSFSROW)

        !> Define temporary variables for common dimensions.
        dim_names_cell_gru = (/DIM_NAME_NCELL, DIM_NAME_NGRU/)
        dim_names_cell_gru_surf = (/DIM_NAME_NCELL, DIM_NAME_NGRU, DIM_NAME_NSUBTYPE/)
        dim_names_cell_gru_sol = (/DIM_NAME_NCELL, DIM_NAME_NGRU, DIM_NAME_NSOL/)

        !> Register this as an active input file.
        file_id = check_for_file_in_list(parameter_files, 'resume_runclass36_3seq', fls%fl(mfk%f883)%fn)

        !> Associate and set the properties of the file.
        associate(input_file => parameter_files(file_id))

            !> Set the file type as binary sequential (replacing any existing container).
            if (allocated(input_file%container)) deallocate(input_file%container)
            allocate(input_file%container, source = io_type_seq())

            !> Allocate the field list (replacing any existing list).
            if (allocated(input_file%fields)) deallocate(input_file%fields)
            allocate(input_file%fields(19))

            !> Populate the field list with variable definitions.
            call helper_update_field_in_list( &
                input_file%fields, 1, 'ALBS', 'Snow albedo', '1', dim_names_cell_gru, error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 2, 'CMAI', 'Aggregated mass of vegetation canopy', 'kg m**-2', dim_names_cell_gru, &
                error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 3, 'GRO', 'Vegetation growth index', '1', dim_names_cell_gru, error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 4, 'QAC', 'Specific humidity of air within the vegetation canopy', 'kg kg**-1', &
                dim_names_cell_gru, error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 5, 'RCAN', 'Intercepted liquid water stored in the canopy', 'kg m**-2', dim_names_cell_gru, &
                error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 6, 'RHOS', 'Density of snow', 'kg m**-3', dim_names_cell_gru, error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 7, 'SCAN', 'Intercepted frozen water stored in the canopy', 'kg m**-2', dim_names_cell_gru, &
                error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 8, 'SNO', 'Mass of snow', 'kg m**-2', dim_names_cell_gru, error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 9, 'TAC', 'Temperature of air within the vegetation canopy', 'K', dim_names_cell_gru, &
                error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 10, 'TBAR', 'Temperature of soil layers', 'K', dim_names_cell_gru_sol, &
                error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 11, 'TBAS', 'Temperature of bedrock at depth of impermeable soils', 'K', dim_names_cell_gru, &
                error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 12, 'TCAN', 'Temperature of the vegetation canopy', 'K', dim_names_cell_gru, &
                error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 13, 'THIC', 'Volumetric frozen water content of soil layers', 'm**3 m**-3', &
                dim_names_cell_gru_sol, error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 14, 'THLQ', 'Volumetric liquid water content of soil layers', 'm**3 m**-3', &
                dim_names_cell_gru_sol, error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 15, 'TPND', 'Temperature of ponded water', 'K', dim_names_cell_gru, error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 16, 'TSFS', 'Ground surface temperature over the tile sub-areas', 'K', dim_names_cell_gru_surf, &
                error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 17, 'TSNO', 'Temperature of snow', 'K', dim_names_cell_gru, error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 18, 'WSNO', 'Liquid water stored in the snow', 'kg m**-2', dim_names_cell_gru, &
                error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 19, 'ZPND', 'Depth of ponded water', 'm', dim_names_cell_gru, error_status = error_status)
        end associate

    end subroutine

    subroutine save_init_prog_variables_class_row(fls, shd)

        use mpi_module
        use model_files_variables
        use sa_mesh_common
        use model_dates

        type(fl_ids) fls
        type(ShedGridParams) shd

        !> For SAVERESUMEFLAG 3
        real(kind = 4), dimension(:, :), allocatable :: ALBSROW, CMAIROW, GROROW, QACROW, RCANROW, &
            RHOSROW, SCANROW, SNOROW, TACROW, TBASROW, &
            TCANROW, TPNDROW, TSNOROW, WSNOROW, ZPNDROW
        real(kind = 4), dimension(:, :, :), allocatable :: TBARROW, THICROW, THLQROW, TSFSROW

        integer NA, NTYPE, NSL, k, ik, jk, iun, ierr

        !> Return if not the head node or if the process is not active.
        if (.not. ISHEADNODE .or. .not. RUNCLASS36_flgs%PROCESS_ACTIVE) return

        !> Local indices.
        NA = shd%NA
        NTYPE = shd%lc%NTYPE
        NSL = shd%lc%IGND

        !> Open the resume state file.
        iun = fls%fl(mfk%f883)%iun
        open(iun, file = trim(adjustl(fls%fl(mfk%f883)%fn)), status = 'replace', action = 'write', &
             form = 'unformatted', access = 'sequential', iostat = ierr)
!todo: condition for ierr.

        !> Allocate and initialize temporary variables.
        allocate(ALBSROW(NA, NTYPE), CMAIROW(NA, NTYPE), GROROW(NA, NTYPE), QACROW(NA, NTYPE), RCANROW(NA, NTYPE), &
                 RHOSROW(NA, NTYPE), SCANROW(NA, NTYPE), SNOROW(NA, NTYPE), TACROW(NA, NTYPE), TBASROW(NA, NTYPE), &
                 TCANROW(NA, NTYPE), TPNDROW(NA, NTYPE), TSNOROW(NA, NTYPE), WSNOROW(NA, NTYPE), ZPNDROW(NA, NTYPE), &
                 TBARROW(NA, NTYPE, NSL), THICROW(NA, NTYPE, NSL), THLQROW(NA, NTYPE, NSL), TSFSROW(NA, NTYPE, 4))
        ALBSROW = 0.0; CMAIROW = 0.0; GROROW = 0.0; QACROW = 0.0; RCANROW = 0.0; RHOSROW = 0.0
        SCANROW = 0.0; SNOROW = 0.0; TACROW = 0.0; TBASROW = 0.0; TCANROW = 0.0; TPNDROW = 0.0
        TSNOROW = 0.0; WSNOROW = 0.0; ZPNDROW = 0.0
        TBARROW = 0.0; THICROW = 0.0; THLQROW = 0.0; TSFSROW = 0.0

        !> Gather the temporary variables.
        do k = 1, shd%lc%NML

            !> Grab the grid and GRU of the current tile.
            ik = shd%lc%ILMOS(k)
            jk = shd%lc%JLMOS(k)

            !> Assign values.
            ALBSROW(ik, jk) = real(vs%tile%albsno(k), kind = 4)
            CMAIROW(ik, jk) = real(vs%tile%cmas(k), kind = 4)
            GROROW(ik, jk) = real(vs%tile%gro(k), kind = 4)
            QACROW(ik, jk) = real(vs%tile%qacan(k), kind = 4)
            RCANROW(ik, jk) = real(vs%tile%lqwscan(k), kind = 4)
            RHOSROW(ik, jk) = real(vs%tile%rhosno(k), kind = 4)
            SCANROW(ik, jk) = real(vs%tile%fzwscan(k), kind = 4)
            SNOROW(ik, jk) = real(vs%tile%sno(k), kind = 4)
            TACROW(ik, jk) = real(vs%tile%tacan(k), kind = 4)
            TBARROW(ik, jk, :) = real(vs%tile%tsol(k, :), kind = 4)
            TBASROW(ik, jk) = real(vs%tile%tbas(k), kind = 4)
            TCANROW(ik, jk) = real(vs%tile%tcan(k), kind = 4)
            THICROW(ik, jk, :) = real(vs%tile%thicsol(k, :), kind = 4)
            THLQROW(ik, jk, :) = real(vs%tile%thlqsol(k, :), kind = 4)
            TPNDROW(ik, jk) = real(vs%tile%tpnd(k), kind = 4)
            TSFSROW(ik, jk, :) = real(vs%tile%tsfs(k, :), kind = 4)
            TSNOROW(ik, jk) = real(vs%tile%tsno(k), kind = 4)
            WSNOROW(ik, jk) = real(vs%tile%lqwssno(k), kind = 4)
            ZPNDROW(ik, jk) = real(vs%tile%zpnd(k), kind = 4)

        end do

        !> Read inital values from the file.
        write(iun) ALBSROW
        write(iun) CMAIROW
        write(iun) GROROW
        write(iun) QACROW
        write(iun) RCANROW
        write(iun) RHOSROW
        write(iun) SCANROW
        write(iun) SNOROW
        write(iun) TACROW
        write(iun) TBARROW
        write(iun) TBASROW
        write(iun) TCANROW
        write(iun) THICROW
        write(iun) THLQROW
        write(iun) TPNDROW
        write(iun) TSFSROW
        write(iun) TSNOROW
        write(iun) WSNOROW
        write(iun) ZPNDROW

        !> Close the file to free the unit.
        close(iun)

        !> Deallocate temporary variables.
        deallocate(ALBSROW, CMAIROW, GROROW, QACROW, RCANROW, &
                   RHOSROW, SCANROW, SNOROW, TACROW, TBASROW, &
                   TCANROW, TPNDROW, TSNOROW, WSNOROW, ZPNDROW, &
                   TBARROW, THICROW, THLQROW, TSFSROW)

    end subroutine

    subroutine RUNCLASS36_finalize(fls, shd)

        use mpi_module
        use model_files_variables
        use sa_mesh_common
        use model_dates
        use FLAGS

        type(fl_ids) :: fls
        type(ShedGridParams) :: shd

        !> For GRU-based end of run prognostic variables.
        real, dimension(:, :), allocatable :: tcan, rcan, sncan, gro, zpnd, tpnd, sno, tsno, albs, rhos
        real, dimension(:, :, :), allocatable :: tbar, thlq, thic
        integer, dimension(:), allocatable :: kc
        character cfmt*3, cfmtt*1000
        character(len = DEFAULT_LINE_LENGTH) line_buffer
        integer NA, NTYPE, NML, NSL, m, k, ik, jk, j, i, ignd, iun, ierr

        !> Return if not the head node or if the process is not active.
        if (.not. ISHEADNODE .or. .not. RUNCLASS36_flgs%PROCESS_ACTIVE) return

        !> Local indices.
        NA = shd%NA
        NTYPE = shd%lc%NTYPE
        NML = shd%lc%NML
        NSL = shd%lc%IGND

        !> Write data to the output summary file.
        if (ECHOTXTMODE) then

            !> CLASS states for prognostic variables.
            allocate(tcan(3, NTYPE), rcan(3, NTYPE), sncan(3, NTYPE), gro(3, NTYPE), zpnd(3, NTYPE), tpnd(3, NTYPE), &
                     sno(3, NTYPE), tsno(3, NTYPE), albs(3, NTYPE), rhos(3, NTYPE), &
                     tbar(3, NTYPE, NSL), thlq(3, NTYPE, NSL), thic(3, NTYPE, NSL), kc(NTYPE))
            tcan = 0.0; rcan = 0.0; sncan = 0.0; gro = 0.0; zpnd = 0.0; tpnd = 0.0
            sno = 0.0; tsno = 0.0; albs = 0.0; rhos = 0.0
            tbar = 0.0; thlq = 0.0; thic = 0.0; kc = 0

            !> Loop through the GRUs.
            do m = 1, NTYPE

                !> Cycle if the GRU does not exist.
                if (count(shd%lc%JLMOS(1:NML) == m) == 0) cycle

                !> Canopy.
                tcan(3, m) = maxval(vs%tile%tcan, shd%lc%JLMOS(1:NML) == m)
                if (tcan(3, m) > 0.0) then
                    tcan(1, m) = sum(vs%tile%tcan, shd%lc%JLMOS(1:NML) == m .and. &
                                     vs%tile%tcan /= 0.0)/count(shd%lc%JLMOS(1:NML) == m .and. vs%tile%tcan /= 0.0)
                    tcan(2, m) = minval(vs%tile%tcan, shd%lc%JLMOS(1:NML) == m .and. vs%tile%tcan /= 0.0)
                end if
                where (tcan < 173.16 .or. tcan > 373.16 .or. tcan == 0.0) tcan = 273.16
                rcan(2, m) = minval(vs%tile%lqwscan, shd%lc%JLMOS(1:NML) == m)
                rcan(3, m) = maxval(vs%tile%lqwscan, shd%lc%JLMOS(1:NML) == m)
                rcan(1, m) = sum(vs%tile%lqwscan, shd%lc%JLMOS(1:NML) == m)/count(shd%lc%JLMOS(1:NML) == m)
                sncan(2, m) = minval(vs%tile%fzwscan, shd%lc%JLMOS(1:NML) == m)
                sncan(3, m) = maxval(vs%tile%fzwscan, shd%lc%JLMOS(1:NML) == m)
                sncan(1, m) = sum(vs%tile%fzwscan, shd%lc%JLMOS(1:NML) == m)/count(shd%lc%JLMOS(1:NML) == m)
                gro(2, m) = minval(vs%tile%gro, shd%lc%JLMOS(1:NML) == m)
                gro(3, m) = maxval(vs%tile%gro, shd%lc%JLMOS(1:NML) == m)
                gro(1, m) = sum(vs%tile%gro, shd%lc%JLMOS(1:NML) == m)/count(shd%lc%JLMOS(1:NML) == m)

                !> Ponded water at surface.
                zpnd(2, m) = minval(vs%tile%zpnd, shd%lc%JLMOS(1:NML) == m)
                zpnd(3, m) = maxval(vs%tile%zpnd, shd%lc%JLMOS(1:NML) == m)
                zpnd(1, m) = sum(vs%tile%zpnd, shd%lc%JLMOS(1:NML) == m)/count(shd%lc%JLMOS(1:NML) == m)
                tpnd(3, m) = maxval(vs%tile%tpnd, shd%lc%JLMOS(1:NML) == m)
                if (tpnd(3, m) > 0.0) then
                    tpnd(1, m) = sum(vs%tile%tpnd, shd%lc%JLMOS(1:NML) == m .and. &
                                     vs%tile%tpnd /= 0.0)/count(shd%lc%JLMOS(1:NML) == m .and. vs%tile%tpnd /= 0.0)
                    tpnd(2, m) = minval(vs%tile%tpnd, shd%lc%JLMOS(1:NML) == m .and. vs%tile%tpnd /= 0.0)
                end if
                where (tpnd < 173.16 .or. tpnd > 373.16 .or. tpnd == 0.0) tpnd = 273.16

                !> Snow.
                sno(2, m) = minval(vs%tile%sno, shd%lc%JLMOS(1:NML) == m)
                sno(3, m) = maxval(vs%tile%sno, shd%lc%JLMOS(1:NML) == m)
                sno(1, m) = sum(vs%tile%sno, shd%lc%JLMOS(1:NML) == m)/count(shd%lc%JLMOS(1:NML) == m)
                tsno(3, m) = maxval(vs%tile%tsno, shd%lc%JLMOS(1:NML) == m)
                if (tsno(3, m) > 0.0) then
                    tsno(1, m) = sum(vs%tile%tsno, shd%lc%JLMOS(1:NML) == m .and. &
                                     vs%tile%tsno /= 0.0)/count(shd%lc%JLMOS(1:NML) == m .and. vs%tile%tsno /= 0.0)
                    tsno(2, m) = minval(vs%tile%tsno, shd%lc%JLMOS(1:NML) == m .and. vs%tile%tsno /= 0.0)
                end if
                where (tsno < 173.16 .or. tsno > 373.16 .or. tsno == 0.0) tsno = 273.16
                if (sno(3, m) > 0.0) then
                    albs(1, m) = sum(vs%tile%albsno, shd%lc%JLMOS(1:NML) == m .and. &
                                     vs%tile%sno > 0.0)/count(shd%lc%JLMOS(1:NML) == m .and. vs%tile%sno > 0.0)
                    albs(2, m) = minval(vs%tile%albsno, shd%lc%JLMOS(1:NML) == m .and. vs%tile%sno > 0.0)
                    albs(3, m) = maxval(vs%tile%albsno, shd%lc%JLMOS(1:NML) == m .and. vs%tile%sno > 0.0)
                end if
                rhos(3, m) = maxval(vs%tile%rhosno, shd%lc%JLMOS(1:NML) == m)
                if (rhos(3, m) > 0.0) then
                    rhos(1, m) = sum(vs%tile%rhosno, shd%lc%JLMOS(1:NML) == m .and. &
                                     vs%tile%rhosno /= 0.0)/count(shd%lc%JLMOS(1:NML) == m .and. vs%tile%rhosno /= 0.0)
                    rhos(2, m) = minval(vs%tile%rhosno, shd%lc%JLMOS(1:NML) == m .and. vs%tile%rhosno /= 0.0)
                end if

                !> Soil.
                do j = 1, NSL
                    tbar(1, m, j) = sum(vs%tile%tsol(:, j), shd%lc%JLMOS(1:NML) == m)/count(shd%lc%JLMOS(1:NML) == m)
                    tbar(2, m, j) = minval(vs%tile%tsol(:, j), shd%lc%JLMOS(1:NML) == m)
                    tbar(3, m, j) = maxval(vs%tile%tsol(:, j), shd%lc%JLMOS(1:NML) == m)
                    thlq(1, m, j) = sum(vs%tile%thlqsol(:, j), shd%lc%JLMOS(1:NML) == m)/count(shd%lc%JLMOS(1:NML) == m)
                    thlq(2, m, j) = minval(vs%tile%thlqsol(:, j), shd%lc%JLMOS(1:NML) == m)
                    thlq(3, m, j) = maxval(vs%tile%thlqsol(:, j), shd%lc%JLMOS(1:NML) == m)
                    thic(1, m, j) = sum(vs%tile%thicsol(:, j), shd%lc%JLMOS(1:NML) == m)/count(shd%lc%JLMOS(1:NML) == m)
                    thic(2, m, j) = minval(vs%tile%thicsol(:, j), shd%lc%JLMOS(1:NML) == m)
                    thic(3, m, j) = maxval(vs%tile%thicsol(:, j), shd%lc%JLMOS(1:NML) == m)
                end do
            end do

            !> Write to file.
            if (NRSOILAYEREADFLAG > 3) then
                ignd = min(NRSOILAYEREADFLAG, NSL)
            else if (NRSOILAYEREADFLAG == 1) then
                ignd = NSL
            else
                ignd = 3
            end if
            write(cfmt, '(i3)') ignd
            call print_echo_txt('')
            call print_echo_txt('End of run prognostic states')
            write(line_buffer, '(3x, (a), i4)') 'Number of GRUs: ', NTYPE
            call print_echo_txt(trim(line_buffer))
            do i = 1, 3
                call print_echo_txt('')
                select case (i)
                    case (1); call print_echo_txt('Average values')
                    case (2); call print_echo_txt('Minimum values')
                    case (3); call print_echo_txt('Maximum values')
                end select
                do m = 1, NTYPE
                    write(line_buffer, "(3x, 'GRU ', i3, ':')") m
                    call print_echo_txt(trim(line_buffer))
                    cfmtt = "(" // trim(adjustl(cfmt)) // "(f10.3), 3(f10.3), " // &
                            "2x, '!> TBAR(1:" // trim(adjustl(cfmt)) // ")/TCAN/TSNO/TPND')"
                    write(line_buffer, cfmtt) ((tbar(i, m, j) - 273.16), j = 1, ignd), &
                        (tcan(i, m) - 273.16), (tsno(i, m) - 273.16), (tpnd(i, m) - 273.16)
                    call print_echo_txt(trim(line_buffer))
                    cfmtt = "(" // trim(adjustl(cfmt)) // "(f10.3), " // trim(adjustl(cfmt)) // "(f10.3), f10.3, " // &
                            "2x, '!> THLQ(1:" // trim(adjustl(cfmt)) // ")/THIC(1:" // trim(adjustl(cfmt)) // ")/ZPND')"
                    write(line_buffer, cfmtt) (thlq(i, m, j), j = 1, ignd), (thic(i, m, j), j = 1, ignd), zpnd(i, m)
                    call print_echo_txt(trim(line_buffer))
                    write(line_buffer, "(6(f10.3), 2x, '!> RCAN/SNCAN/SNO/ALBS/RHOS/GRO')") &
                        rcan(i, m), sncan(i, m), sno(i, m), albs(i, m), rhos(i, m), gro(i, m)
                    call print_echo_txt(trim(line_buffer))
                end do
            end do

        end if

    end subroutine

end module
