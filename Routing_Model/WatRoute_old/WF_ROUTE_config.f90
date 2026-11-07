module WF_ROUTE_config

    implicit none

    !> *****************************************************************
    !> Object variables.
    !> *****************************************************************

    !>
    !> Description: WF_RTE_flags contains flags to enable and
    !>              disable the module and control its output.
    !>

    type WF_RTE_flags

        !> Flag used to enable the module.
        logical :: PROCESS_ACTIVE = .true.

        !> Flag to control the reservoir release function used in
        !> WF_ROUTE.f.
        !>  2 = 2-parameter power release function.
        !>  5 = 5-parameter polynomial release function.
!?        integer :: RESVRELSWFB = 2

        !> Channel length 'rl' flag for WF_ROUTE.f.
!todo: Verify units if read from shed file.
        !>  0 = 'rl' is calculated using WF_AL and WF_A1 (default).
        !>  1 = Values are taken from the 'chnllength' attribute from
        !>      the drainage database/r2c shed file.
        integer :: RLFLAG = 0

        !> Bankfull/capacity 'cap' flag for WF_ROUTE.f
!todo: Verify units if read from shed file.
        !>  0 = 'cap' is calculated using WF_DA, WF_A2, WF_A3, and WF_A4
        !>      (default).
        !>  1 = Values are taken from the 'bankfull' attribute from the
        !>      drainage database/r2c shed file.
        integer :: CAPFLAG = 0

        integer :: RTE_TS = 1800

    end type

    type WF_RTE_parameters

        !> Channel roughness coefficients.
        !* r2: River channel roughness coefficient.
        !* r1: Overbank channel roughness coefficient.
        real, dimension(:), allocatable :: r2, r1

        !> Fitting coefficients.
        !* aa1: Channel length coefficient.
        !* aa2: Bankfull area coefficient.
        !* aa3: Bankfull area coefficient.
        !* aa4: Bankfull area coefficient.
        real, dimension(:), allocatable :: aa1, aa2, aa3, aa4

    end type

    type(WF_RTE_parameters), save :: wfp

    !> WF_RTE_flgs: Configuration flags for the module.
    type(WF_RTE_flags), save :: WF_RTE_flgs

    !> *****************************************************************
    !> Local variables.
    !> *****************************************************************

    !>
    !> Description: Variables used by WF_ROUTE. These variables are used
    !> by WF_ROUTE and are only accessible to code that use this module.
    !>

    !> STREAMFLOW VARIABLES
    !* WF_NL: NUMBER OF DATA POINTS
    !* WF_MHRD: NUMBER OF HOURS OF DATA PER MONTH
    !* WF_KT: HOURLY INCREMENT FOR STREAMFLOW INPUT (24 = DAILY)
    !* WF_A1: Channel fitting parameter for average channel length (default: 1.0).
    !* WF_A2: Channel fitting parameter for average bankfull capacity (default: 11.0).
    !* WF_A3: Channel fitting parameter for average bankfull capacity (default: 0.43).
    !* WF_A4: Channel fitting parameter for average bankfull capacity (default: 1.0).
    integer WF_NAA, WF_NL, WF_MHRD, WF_KT

    !> RESERVOIR VARIABLES
    integer, dimension(:), allocatable :: WF_RES
    real, dimension(:), allocatable :: WF_B1, WF_B2, WF_RESSTORE
    real, dimension(:), allocatable :: WF_B3, WF_B4, WF_B5

    !> FOR BASEFLOW INITIALIZATION
    integer JAN

    !* WF_R1: MANNING'S N FOR RIVER CHANNEL
    !* WF_R2: OPTIMIZED RIVER ROUGHNESS FACTOR
    !* WF_QO2: SIMULATED STREAMFLOW VALUE
    real, dimension(:), allocatable :: WF_NHYD, WF_QBASE, WF_QI1, &
        WF_QO1, WF_QR, WF_STORE1

    !> RESERVOIR MEASUREMENTS:
    !* WR_NREL: NUMBER OF DATA POINTS
    !* WF_KTR: HOURLY INCREMENT FOR RESERVOIR INPUR (24 = DAILY)
    integer WF_NREL, WF_KTR, WF_NORESV_CTRL
    integer WF_ROUTETIMESTEP, WF_TIMECOUNT

    contains

    !> *****************************************************************
    !> Subroutines.
    !> *****************************************************************

    !>
    !> Description: Check for the existence of input files, open them,
    !>              print diagnostic information, skip records, and open
    !>              the output files, in preparation for running the
    !>              WF_ROUTE process.
    !>
    subroutine WF_ROUTE_init(fls, shd)

        use model_files_variables
        use sa_mesh_common
        use model_dates
!-        use FLAGS
        use input_fields

        type(fl_ids) :: fls
        type(ShedGridParams), intent(in) :: shd

        !> Temporary variables for frequently accessed terms.
        !* NA: Number of grid cells.
        !* NS: Number of streamflow gauges.
        !* NR: Number of reservoir outlets.
        integer NS, NR, NA

        !> Local variables.
        !* iun: Unit number.
        !* ierr: Error return from external calls.
        integer iun, ierr, i, z

        !> Return if not the head node or if the process is not active.
        if (.not. ISHEADNODE .or. .not. WF_RTE_flgs%PROCESS_ACTIVE) then
            return
        else
            call print_new_section("WF_ROUTE is active.")
            call increase_tab()
        end if

        NA = shd%NA
        WF_NAA = NA - shd%NAA

        allocate(WF_NHYD(NA), WF_QR(NA), &
                 WF_QBASE(NA), WF_QI1(NA), WF_QO1(NA), &
                 WF_STORE1(NA))

        WF_NHYD = 0.0
        WF_QBASE = 0.0
        WF_QI1 = 0.0
        WF_QO1 = 0.0
        WF_QR = 0.0
        WF_STORE1 = 0.0

        !> *************************************************************
        !>  Open and read in values from MESH_input_reservoir.txt file
        !> *************************************************************

        WF_NORESV_CTRL = 0
        NR = fms%rsvr%n
        WF_KTR = fms%rsvr%rlsmeas%dts

        if (NR > 0) then

            !> Allocate and initialize reservoir variables.
            allocate(WF_RES(NR), &
                     WF_B1(NR), WF_B2(NR), WF_B3(NR), WF_B4(NR), WF_B5(NR), &
                     WF_RESSTORE(NR))
            WF_RESSTORE = 0.0
            WF_B1 = fms%rsvr%rls%b1
            WF_B2 = fms%rsvr%rls%b2
            WF_B3 = fms%rsvr%rls%b3
            WF_B4 = fms%rsvr%rls%b4
            WF_B5 = fms%rsvr%rls%b5

            do i = 1, NR
                if (WF_B3(i) > 0.0) then
                    fms%rsvr%rls%cfn(i) = 3
                else if (WF_B1(i) > 0.0) then
                    fms%rsvr%rls%cfn(i) = 2
                else
                    fms%rsvr%rls%cfn(i) = 1
                end if
                if (WF_B1(i) == 0.0) then
                    WF_NORESV_CTRL = WF_NORESV_CTRL + 1
                end if
            end do
        end if

        !> *********************************************************************
        !> Open and read in values from MESH_input_streamflow.txt file
        !> *********************************************************************

        NS = fms%stmg%n

        WF_KT = fms%stmg%qomeas%dts

        WF_ROUTETIMESTEP = 900
        WF_TIMECOUNT = 0

        !* JAN: The first time throught he loop, jan = 1. Jan will equal 2 after that.
        JAN = 1

        !> Allocate and initialize variables to a background 'no data' value.
        ierr = 0
        z = 0
        call allocate_field(wfp%r1, (/shd%NRVR/), error_status = z); if (z /= 0) ierr = z
        call allocate_field(wfp%r2, (/shd%NRVR/), error_status = z); if (z /= 0) ierr = z
        call allocate_field(wfp%aa1, (/shd%NRVR/), error_status = z); if (z /= 0) ierr = z
        call allocate_field(wfp%aa2, (/shd%NRVR/), error_status = z); if (z /= 0) ierr = z
        call allocate_field(wfp%aa3, (/shd%NRVR/), error_status = z); if (z /= 0) ierr = z
        call allocate_field(wfp%aa4, (/shd%NRVR/), error_status = z); if (z /= 0) ierr = z
        JAN = NO_DATA_INT
        vs%grid%qo = NO_DATA_REAL
        vs%grid%stgch = NO_DATA_REAL
        vs%grid%qi = NO_DATA_REAL
        if (ierr /= 0) then
            call print_error("An error occurred allocating input fields.")
            call program_abort()
        end if

        !> Assign input variables from field list.
        do i = 1, size(input_field_list)

            !> Assign recognized fields.
            associate(field => input_field_list(i))
                select case (field%field_name)
                    case ('WF_R2')
                        call get_field_value(field, wfp%r2, DIM_NAME_K)
                    case ('WF_R1')
                        call get_field_value(field, wfp%r1, DIM_NAME_K)
                    case ('AA1')
                        call get_field_value(field, wfp%aa1, DIM_NAME_K)
                    case ('AA2')
                        call get_field_value(field, wfp%aa2, DIM_NAME_K)
                    case ('AA3')
                        call get_field_value(field, wfp%aa3, DIM_NAME_K)
                    case ('AA4')
                        call get_field_value(field, wfp%aa4, DIM_NAME_K)
                    case ('JAN')
                        call get_field_value(field, JAN)
                    case ('QO2', 'QO', 'GRID_QO')
                        call get_field_value(field, vs%grid%qo, DIM_NAME_N)
                    case ('STORE2', 'STGCH', 'GRID_STGCH')
                        call get_field_value(field, vs%grid%stgch, DIM_NAME_N)
                    case ('QI2', 'QI', 'GRID_QI')
                        call get_field_value(field, vs%grid%qi, DIM_NAME_N)
                end select
            end associate
        end do

        !> Assign defaults or derived fields where possible.
        ierr = 0
        z = 0
        call check_field_assigned(wfp%r2, 'WF_R2', is_error = .true., error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(wfp%r1, 'WF_R1', 2.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(wfp%aa1, 'AA1', 1.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(wfp%aa2, 'AA2', 11.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(wfp%aa3, 'AA3', 0.43, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(wfp%aa4, 'AA4', 1.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(JAN, 'JAN', 1, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(vs%grid%qo, 'QO2', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(vs%grid%stgch, 'STORE2', 0.0, error_status = z); if (z /= 0) ierr = z
        call check_field_assigned(vs%grid%qi, 'QI2', 0.0, error_status = z); if (z /= 0) ierr = z

        !> Check for errors.
        if (ierr /= 0) then
            call reset_tab()
            call print_error("Variables required by the module are not assigned or were missing during initialization.")
            call program_abort()
        end if

    end subroutine

    subroutine WF_ROUTE_resume_read(fls, shd)

        use model_files_variables
        use sa_mesh_common
        use input_fields
        use mesh_legacy_io

        type(fl_ids) :: fls
        type(ShedGridParams), intent(in) :: shd

        !> Local variables.
!-        integer(kind = 4) JAN_i4, WF_TIMECOUNT_i4
!-        real(kind = 4), dimension(:), allocatable :: qo_r4, stgch_r4, qi_r4
!-        integer ierr, iun

        !> Local variables.
        integer file_id, error_status
        character(len = SHORT_FIELD_LENGTH) dim_names_cell(1)

        !> Return if not the head node or if the process is not active.
        if (.not. ISHEADNODE .or. .not. WF_RTE_flgs%PROCESS_ACTIVE) return

        !> Open the resume file.
!-        iun = fls%fl(mfk%f883)%iun
!-        open(iun, file = trim(adjustl(fls%fl(mfk%f883)%fn)) // '.wf_route', status = 'old', action = 'read', &
!-             form = 'unformatted', access = 'sequential', iostat = ierr)
!todo: condition for ierr.

        !> Allocate and initialize local variables.
!-        allocate(qo_r4(shd%NA), stgch_r4(shd%NA), qi_r4(shd%NA))
!-        qo_r4 = 0.0
!-        stgch_r4 = 0.0
!-        qi_r4 = 0.0

        !> Read inital values from the file.
!-        read(iun) JAN_i4
!-        read(iun) WF_TIMECOUNT_i4
!-        read(iun)
!-        read(iun)
!-        read(iun)
!-        read(iun)
!-        read(iun)
!-        read(iun)
!-        read(iun) qo_r4
!-        read(iun) stgch_r4
!-        read(iun) qi_r4
!-        read(iun)
!-        read(iun)

        !> Transfer variables.
!-        JAN = int(JAN_i4, kind(JAN))
!-        WF_TIMECOUNT = int(WF_TIMECOUNT_i4, kind(WF_TIMECOUNT))
!-        vs%grid%qo = real(qo_r4, kind(vs%grid%qo))
!-        vs%grid%stgch = real(stgch_r4, kind(vs%grid%stgch))
!-        vs%grid%qi = real(qi_r4, kind(vs%grid%qi))

        !> Close the file to free the unit.
!-        close(iun)

        !> Define temporary variables for common dimensions.
        dim_names_cell = (/DIM_NAME_NCELL/)

        !> Register this as an active input file.
        file_id = check_for_file_in_list(parameter_files, 'resume_wf_route_5seq', trim(fls%fl(mfk%f883)%fn) // '.wf_route')

        !> Associate and set the properties of the file.
        associate(input_file => parameter_files(file_id))

            !> Set the file type as binary sequential (replacing any existing container).
            if (allocated(input_file%container)) deallocate(input_file%container)
            allocate(input_file%container, source = io_type_seq())

            !> Allocate the field list (replacing any existing list).
            if (allocated(input_file%fields)) deallocate(input_file%fields)
            allocate(input_file%fields(11))

            !> Populate the field list with variable definitions.
            call helper_update_field_in_list( &
                input_file%fields, 1, 'JAN', 'WATFLOOD-style initialization marker (JAN)', '1', error_status = error_status)
            allocate(input_file%fields(1)%field, source = model_variable_int(dat = NO_DATA_INT))
            call helper_update_field_in_list( &
                input_file%fields, 9, 'QO2', 'Outflow', 'm**3 s**-1', dim_names_cell, error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 10, 'STORE2', 'Storage', 'm**3', dim_names_cell, error_status = error_status)
            call helper_update_field_in_list( &
                input_file%fields, 11, 'QI2', 'Inflow', 'm**3 s**-1', dim_names_cell, error_status = error_status)
        end associate

    end subroutine

    subroutine WF_ROUTE_resume_read_nots(fls, shd)

        use model_files_variables
        use sa_mesh_common

        type(fl_ids) :: fls
        type(ShedGridParams), intent(in) :: shd

        !> Local variables.
!-        integer(kind = 4) JAN_i4
!-        real(kind = 4), dimension(:), allocatable :: qo_r4, stgch_r4, qi_r4
!-        integer ierr, iun

        !> Return if not the head node or if the process is not active.
        if (.not. ISHEADNODE .or. .not. WF_RTE_flgs%PROCESS_ACTIVE) return

        !> Open the resume file.
!-        iun = fls%fl(mfk%f883)%iun
!-        open(iun, file = trim(adjustl(fls%fl(mfk%f883)%fn)) // '.wf_route', status = 'old', action = 'read', &
!-             form = 'unformatted', access = 'sequential', iostat = ierr)
!todo: condition for ierr.

        !> Allocate and initialize local variables.
!-        allocate(qo_r4(shd%NA), stgch_r4(shd%NA), qi_r4(shd%NA))
!-        qo_r4 = 0.0
!-        stgch_r4 = 0.0
!-        qi_r4 = 0.0

        !> Read inital values from the file.
!-        read(iun) JAN_i4
!-        read(iun)
!-        read(iun)
!-        read(iun)
!-        read(iun)
!-        read(iun)
!-        read(iun)
!-        read(iun)
!-        read(iun) qo_r4
!-        read(iun) stgch_r4
!-        read(iun) qi_r4
!-        read(iun)
!-        read(iun)

        !> Transfer variables.
!-        JAN = int(JAN_i4, kind(JAN))
!-        vs%grid%qo = real(qo_r4, kind(vs%grid%qo))
!-        vs%grid%stgch = real(stgch_r4, kind(vs%grid%stgch))
!-        vs%grid%qi = real(qi_r4, kind(vs%grid%qi))

        !> Close the file to free the unit.
!-        close(iun)

        !> Call base routine.
        call WF_ROUTE_resume_read(fls, shd)

    end subroutine

    subroutine WF_ROUTE_resume_save(fls, shd)

        use mpi_module
        use model_files_variables
        use sa_mesh_common

        type(fl_ids) :: fls
        type(ShedGridParams) :: shd

        !> Local variables.
        integer ierr, iun

        !> Return if not the head node or if the process is not active.
        if (.not. ISHEADNODE .or. .not. WF_RTE_flgs%PROCESS_ACTIVE) return

        !> Return if not the head node.
        if (ipid /= 0) return

        !> Open the resume file.
        iun = fls%fl(mfk%f883)%iun
        open(iun, file = trim(adjustl(fls%fl(mfk%f883)%fn)) // '.wf_route', status = 'replace', action = 'write', &
             form = 'unformatted', access = 'sequential', iostat = ierr)
!todo: condition for ierr.

        !> Write the current state of these variables to the file.
        write(iun) int(JAN, kind = 4)
        write(iun) int(WF_TIMECOUNT, kind = 4)
        write(iun)
        write(iun)
        write(iun)
        write(iun)
        write(iun)
        write(iun)
        write(iun) real(vs%grid%qo, kind = 4)
        write(iun) real(vs%grid%stgch, kind = 4)
        write(iun) real(vs%grid%qi, kind = 4)
        write(iun)
        write(iun)

        !> Close the file to free the unit.
        close(iun)

    end subroutine

    subroutine WF_ROUTE_finalize(fls, shd)

        use mpi_module
        use model_files_variables
        use sa_mesh_common

        type(fl_ids) :: fls
        type(ShedGridParams) :: shd

        !> Return if not the head node or if the process is not active.
        if (.not. ISHEADNODE .or. .not. WF_RTE_flgs%PROCESS_ACTIVE) return

    end subroutine

end module
