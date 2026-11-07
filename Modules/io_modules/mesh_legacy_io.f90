!> Description:
!>  Module containing routines to support legacy I/O formats.
module mesh_legacy_io

    !> 'mesh_io': For I/O field and file types, options and routines.
    use mesh_io

    !> 'txt_io': For text-based I/O utilities.
    use txt_io

    !> 'variable_names': For variable name constants.
    use variable_names

    !> 'strings': For 'uppercase' function.
    use strings, only: uppercase

    implicit none

    interface helper_assign_field
        module procedure helper_assign_field_from_args
        module procedure helper_assign_field_from_line
    end interface

    private helper_create_variable_name
    private helper_print_error_message
    private helper_assign_field
    private helper_assign_field_from_args
    private helper_assign_field_from_line

    contains

    character(len = SHORT_FIELD_LENGTH) function helper_create_variable_name(variable_name, set_id, set_dim_name)

        !> Input/output variables.
        character(len = *), intent(in) :: variable_name
        integer, intent(in), optional :: set_id
        character(len = *), intent(in), optional :: set_dim_name

        !> Local variables.
        character(len = SHORT_FIELD_LENGTH) line

        !> Build variable name.
        helper_create_variable_name = adjustl(variable_name)
        if (present(set_id)) then
            write(line, *) set_id
            if (present(set_dim_name)) line = trim(set_dim_name) // "=" // trim(adjustl(line))
            helper_create_variable_name = trim(helper_create_variable_name) // " " // trim(adjustl(line))
        end if

    end function

    subroutine helper_print_error_message(variable_name, set_id, set_dim_name, line_number, error_status)

        !> Input/output variables.
        character(len = *), intent(in) :: variable_name
        integer, intent(in), optional :: set_id
        character(len = *), intent(in), optional :: set_dim_name
        integer, intent(in), optional :: line_number
        integer, intent(in), optional :: error_status

        !> Local variables.
        character(len = :), allocatable :: message_text
        character(len = SHORT_FIELD_LENGTH) line

        !> Build error message.
        message_text = &
            "An error occurred reading '" // trim(helper_create_variable_name(variable_name, set_id, set_dim_name)) // "'"
        if (present(line_number)) then
            write(line, *) line_number
            message_text = message_text // " at line " // trim(adjustl(line))
        end if
        if (present(error_status)) then
            write(line, *) error_status
            message_text = message_text // " (Code: " // trim(adjustl(line)) // ")"
        end if
        message_text = message_text // "."
        call print_error(message_text)

    end subroutine

    subroutine helper_assign_field_from_args( &
        field_list, field_index, read_status, variable_name, input_field, &
        start_id, stop_id, dim_name, check_numeric, set_id, set_dim_name, line_number, &
        error_status)

        !> Modules.
        use parse_utilities, only: is_numeric

        !> Input/output variables.
        type(io_field) field_list(:)
        integer, intent(in) :: field_index
        integer, intent(in) :: read_status
        character(len = *), intent(in) :: variable_name
        character(len = SHORT_FIELD_LENGTH), intent(in) :: input_field(:)
        integer, intent(in), optional :: start_id
        integer, intent(in), optional :: stop_id
        character(len = *), intent(in), optional :: dim_name
        logical, intent(in), optional :: check_numeric
        integer, intent(in), optional :: set_id
        character(len = *), intent(in), optional :: set_dim_name
        integer, intent(in), optional :: line_number
        integer error_status

        !> Local variables.
        integer id_1, id_2, i
        character(len = SHORT_FIELD_LENGTH) dim_length
        character(len = :), allocatable :: message
        logical size_test, numeric_test

        !> Indices.
        id_1 = 1
        if (present(start_id)) id_1 = start_id
        id_2 = size(input_field)
        if (present(stop_id)) id_2 = stop_id
        size_test = size(input_field) >= id_2

        !> Check numeric.
        if (size_test) then
            numeric_test = .true.
            if (present(check_numeric)) then
                if (check_numeric) then
                    do i = id_1, id_2
                        numeric_test = (numeric_test .and. is_numeric(input_field(i)))
                        if (.not. numeric_test) exit
                    end do
                end if
            end if
        end if

        !> Assign value.
        if (read_status == 0 .and. size_test .and. numeric_test) then
            field_list(field_index)%label = adjustl(helper_create_variable_name(variable_name, set_id, set_dim_name))
            message = "Found the variable '" // trim(field_list(field_index)%label) // "'"
            if (id_1 == id_2) then
                allocate(field_list(field_index)%field, source = model_variable_char(dat = adjustl(input_field(id_1))))
                message = message // " with 1 value"
            else
                allocate(field_list(field_index)%field, source = model_variable_char1d(dat = input_field(id_1:id_2)))
                write(dim_length, *) (id_2 - id_1 + 1)
                message = message // " with " // trim(adjustl(dim_length)) // " values"
            end if
            if (present(dim_name)) then
                allocate(field_list(field_index)%dim_names(1))
                field_list(field_index)%dim_names = dim_name
                message = message // " mapped to the dimension '" // trim(field_list(field_index)%dim_names(1)) // "'"
            end if
            call print_info(message // ".")
        else
            call helper_print_error_message(variable_name, set_id, set_dim_name, line_number, error_status)
            error_status = 1
        end if

    end subroutine

    subroutine helper_assign_field_from_line( &
        field_list, field_index, read_status, variable_name, input_field, set_id, set_dim_name, line_number, error_status)

        !> Input/output variables.
        type(io_field) field_list(:)
        integer, intent(in) :: field_index
        integer, intent(in) :: read_status
        character(len = *), intent(in) :: variable_name
        character(len = *), intent(in) :: input_field
        integer, intent(in), optional :: set_id
        character(len = *), intent(in), optional :: set_dim_name
        integer, intent(in), optional :: line_number
        integer error_status

        !> Assign value.
        if (read_status == 0) then
            field_list(field_index)%label = adjustl(helper_create_variable_name(variable_name, set_id, set_dim_name))
            allocate(field_list(field_index)%field, source = model_variable_char(dat = adjustl(input_field)))
            call print_info("Found the variable '" // trim(field_list(field_index)%label) // "' with 1 value.")
        else
            call helper_print_error_message(variable_name, set_id, set_dim_name, line_number, error_status)
            error_status = 1
        end if

    end subroutine

    subroutine helper_update_field_in_list( &
        field_list, field_index, variable_name, description, units, dim_names, set_id, set_dim_name, error_status)

        !> Input/output variables.
        type(io_field), allocatable :: field_list(:)
        integer, intent(in) :: field_index
        character(len = *), intent(in) :: variable_name
        character(len = *), intent(in), optional :: description
        character(len = *), intent(in), optional :: units
        character(len = *), intent(in), optional :: dim_names(:)
        integer, intent(in), optional :: set_id
        character(len = *), intent(in), optional :: set_dim_name
        integer error_status

        !> Check the size of the field list.
        if (.not. allocated(field_list) .or. size(field_list) < field_index) then

            !> Expand the list.
            call expand_field_list(field_list, increment = 1, error_status = error_status)
        else

            !> Set the return variable.
            error_status = 0
        end if

        !> Update the attributes.
        associate(field => field_list(field_index))

            !> Create the field name.
            field%label = adjustl(helper_create_variable_name(variable_name, set_id, set_dim_name))

            !> Add meta information.
            associate(meta => field%meta)
                meta%label = adjustl(variable_name)
                if (present(description)) meta%description = adjustl(description)
                if (present(units)) meta%units = adjustl(units)
            end associate

            !> Allocate (or reallocate) and assign dimension names.
            if (present(dim_names)) then
                if (allocated(field%dim_names)) deallocate(field%dim_names)
                allocate(field%dim_names(size(dim_names)))
                field%dim_names = dim_names
            end if
        end associate

    end subroutine

    !> Description:
    !>  Subroutine to read the 'MESH_input_soil_levels.txt' to a field list.
    !>
    !> Variables:
    !*  input_file: Input file.
    !*  full_path: Path to replace the default file name (optional).
    !*  global_list: Existing or empty field list of variables read from any files (optional).
    !*  error_status: Return status (0: normal).
    subroutine read_soil_levels_txt_to_field_list(input_file, full_path, global_list, error_status)

        !> Input/output variables.
        type(io_file) input_file
        character(len = *), intent(in), optional :: full_path
        type(io_field), allocatable, optional :: global_list(:)
        integer, intent(out) :: error_status

        !> Local variables.
        integer layer_count, line_count, i, ierr
        character(len = SHORT_FIELD_LENGTH), allocatable :: args(:)

        !> Reset the file structure.
        if (allocated(input_file%container)) deallocate(input_file%container)
        allocate(input_file%container, source = io_type_char_delimited())
        if (allocated(input_file%fields)) deallocate(input_file%fields)

        !> Open the file.
        if (present(full_path)) then
            input_file%full_path = adjustl(full_path)
        else
            input_file%full_path = 'MESH_input_soil_levels.txt'
        end if
        call open_input_file(input_file, error_status = error_status)
        if (error_status /= 0) return

        !> Count the number of levels (lines in the file).
        call count_txt_lines(input_file%iunit, layer_count, error_status)
        if (error_status == 0) then
            if (layer_count == 0) then
                error_status = 1
            else
                call assign_field(layer_count, args, ierr)
            end if
        end if
        if (error_status /= 0) then
            call print_error("The file contains no data.")
        else

            !> Allocate the field list.
            allocate(input_file%fields(1))

            !> Save the number of layers.
            call helper_assign_field(input_file%fields, 1, ierr, DIM_NAME_NSOL, args, 1, error_status = error_status)
        end if

        !> Read the layers from the file.
        if (error_status == 0) then

            !> Rewind file.
            rewind(input_file%iunit)

            !> Initialize line count.
            line_count = 0

            !> Read the file.
            call expand_field_list(input_file%fields, increment = layer_count, error_status = error_status)
            if (error_status == 0) then
                do i = 1, layer_count
                    call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)

                    !> Save the field.
                    call helper_assign_field( &
                        input_file%fields, (i + 1), ierr, VN_DZSOL, args, 1, 1, &
                        DIM_NAME_BASIN, .true., i, DIM_NAME_LAYER, line_count, error_status)
                end do
            end if
        end if

        !> Close unit.
        close(input_file%iunit)

        !> Check for errors.
        if (error_status /= 0) then
            call print_error("Errors occurred reading the file.")
        else

            !> Map fields.
            call create_field_maps_from_list(input_file%fields, ierr)

            !> Update the global list.
            if (present(global_list)) call copy_field_list(global_list, input_file%fields, error_status)
        end if

    end subroutine

    !> Description:
    !>  Subroutine to read the 'MESH_parameters_CLASS.ini' to a field list.
    !>
    !> Variables:
    !*  global_list: Existing field list of variables read from other files (for field dependencies).
    !*  input_file: Input file.
    !*  full_path: Path to replace the default file name (optional).
    !*  error_status: Return status (0: normal).
    subroutine read_parameters_class_ini_to_field_list(global_list, input_file, full_path, error_status)

        !> Input/output variables.
        type(io_field), allocatable :: global_list(:)
        type(io_file) input_file
        character(len = *), intent(in), optional :: full_path
        integer, intent(out) :: error_status

        !> Local variables.
        integer gru_count, read_soil_layer_count, layer_count, line_count, i, ierr
        character(len = SHORT_FIELD_LENGTH), allocatable :: args(:)
        character(len = LONG_FIELD_LENGTH) line

        !> Reset the file structure.
        if (allocated(input_file%container)) deallocate(input_file%container)
        allocate(input_file%container, source = io_type_char_delimited())
        if (allocated(input_file%fields)) deallocate(input_file%fields)

        !> Open the file.
        if (present(full_path)) then
            input_file%full_path = adjustl(full_path)
        else
            input_file%full_path = 'MESH_parameters_CLASS.ini'
        end if
        call open_input_file(input_file, error_status = error_status)
        if (error_status /= 0) return

        !> Look for 'NRSOILAYEREADFLAG' and soil layer count 'NSOL' in the global field list.
        read_soil_layer_count = NO_DATA_INT
        layer_count = NO_DATA_INT
        do i = 1, size(global_list)

            !> Find and assign fields.
            if (global_list(i)%label == 'NRSOILAYEREADFLAG' .and. allocated(global_list(i)%field)) then
                call assign_field(global_list(i)%field, read_soil_layer_count, ierr)
                if (ierr /= 0) read_soil_layer_count = NO_DATA_INT
            else if (global_list(i)%label == DIM_NAME_NSOL .and. allocated(global_list(i)%field)) then
                call assign_field(global_list(i)%field, layer_count, ierr)
                if (ierr /= 0) layer_count = NO_DATA_INT
            else if (read_soil_layer_count /= NO_DATA_INT .and. layer_count /= NO_DATA_INT) then

                !> Exit loop.
                exit
            end if
        end do

        !> Check the state of 'NRSOILAYEREADFLAG'.
        if (read_soil_layer_count /= NO_DATA_INT .and. read_soil_layer_count > 0) then
            if (layer_count == NO_DATA_INT) then

                !> Print a warning that the number of soil layers is not defined.
                call print_warning( &
                    "'NRSOILAYEREADFLAG' is active but the number of soil layers has not been defined. The flag has no effect.")
            else if (layer_count < 3) then

                !> Print a warning if the number of soil layers is too small.
                write(line, *) layer_count
                call print_warning( &
                    "'NRSOILAYEREADFLAG' has been activated with an incompatible number of soil layers (" // &
                    trim(adjustl(line)) // "). The flag has no effect.")
                layer_count = max(layer_count, 3)
            else

                !> The flag specifies to read a user-specified number of layers.
                if (read_soil_layer_count > 3) layer_count = min(read_soil_layer_count, layer_count)

                !> Print a remark about the number of soil layer expected to be read from file.
                write(line, *) layer_count
                call print_remark( &
                    "'NRSOILAYEREADFLAG' has been activated to read " // trim(adjustl(line)) // &
                    " values for all soil-related parameters and initial conditions from the file.")
            end if
        else

            !> Default case is the CLASS legacy of three layers to read from the file.
            layer_count = 3
        end if

        !> Allocate the field list for meta fields and end-of-file dates.
        allocate(input_file%fields(22))

        !> Initialize the line and GRU counts.
        line_count = 0
        gru_count = 0

        !> Read and save meta information.
        if (error_status == 0) then
            call read_txt_line(input_file%iunit, line, preserve_blank = .true., line_count = line_count, ierr = ierr)
            call helper_assign_field( &
                input_file%fields, 1, ierr, 'CLASSOUTTITLE', line, line_number = line_count, error_status = error_status)
            call read_txt_line(input_file%iunit, line, preserve_blank = .true., line_count = line_count, ierr = ierr)
            call helper_assign_field( &
                input_file%fields, 2, ierr, 'CLASSOUTNAME', line, line_number = line_count, error_status = error_status)
            call read_txt_line(input_file%iunit, line, preserve_blank = .true., line_count = line_count, ierr = ierr)
            call helper_assign_field( &
                input_file%fields, 3, ierr, 'CLASSOUTPLACE', line, line_number = line_count, error_status = error_status)
        end if

        !> Read and save configuration information.
        if (error_status == 0) then
            call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
            if (ierr /= 0) then
                error_status = ierr
                call helper_print_error_message(DIM_NAME_DEGLAT, line_number = line_count, error_status = ierr)
            else
                call helper_assign_field( &
                    input_file%fields, 4, ierr, DIM_NAME_DEGLAT, args, 1, 1, line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, 5, ierr, DIM_NAME_DEGLON, args, 2, 2, line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, 6, ierr, VN_ZREFMOMENTUM, args, 3, 3, line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, 7, ierr, VN_ZREFTHERMODYN, args, 4, 4, line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, 8, ierr, 'ZREFBLEND', args, 5, 5, line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, 9, ierr, 'GRIDCOVER', args, 6, 6, line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, 10, ierr, 'ILG', args, 7, 7, line_number = line_count, error_status = error_status)
                if (.not. size(args) > 8) then
                    ierr = 1
                else
                    call assign_field(args(9), gru_count, ierr)
                end if
                if (ierr /= 0) then
                    error_status = 1
                    call helper_print_error_message(DIM_NAME_NGRU, line_number = line_count)
                end if
            end if
        end if
        if (gru_count == 0) then

            !> Print an error if the number of GRUs could not be read from the file.
            call print_error("The number of GRUs could not be read from the file.")
            error_status = 1
        else

            !> Expand field list for the fields of each GRU.
            call expand_field_list(input_file%fields, increment = (39*gru_count), error_status = error_status)
        end if

        !> Read and save the GRU fields.
        do i = 1, gru_count
            if (error_status == 0) then
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message('FCAN', i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 11), ierr, 'FCAN', args, 1, 5, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 12), ierr, 'PAMX', args, 6, 9, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message('LNZ0', i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 13), ierr, 'LNZ0', args, 1, 5, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 14), ierr, 'PAMN', args, 6, 9, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message('ALVC', i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 15), ierr, 'ALVC', args, 1, 5, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 16), ierr, 'CMAS', args, 6, 9, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message('ALIC', i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 17), ierr, 'ALIC', args, 1, 5, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 18), ierr, 'ROOT', args, 6, 9, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message('RSMN', i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 19), ierr, 'RSMN', args, 1, 4, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 20), ierr, 'QA50', args, 5, 8, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message('VPDA', i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 21), ierr, 'VPDA', args, 1, 4, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 22), ierr, 'VPDB', args, 5, 8, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message('PSGA', i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 23), ierr, 'PSGA', args, 1, 4, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 24), ierr, 'PSGB', args, 5, 8, &
                        DIM_NAME_CANOPY, .true., i, DIM_NAME_GRU, line_count, error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message('DRN', i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 25), ierr, 'DRN', args, 1, 1, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 26), ierr, 'SDEP', args, 2, 2, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 27), ierr, 'FARE', args, 3, 3, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 28), ierr, 'DD', args, 4, 4, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    if (error_status == 0) input_file%fields((i - 1)*39 + 28)%meta%units = 'km km**-2'
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message('XSLP', i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 29), ierr, 'XSLP', args, 1, 1, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 30), ierr, 'GRKEF', args, 2, 2, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 31), ierr, 'MANN', args, 3, 3, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 32), ierr, 'KSAT', args, 4, 4, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 33), ierr, 'MOSAICID', args, 5, 5, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message(VN_SOLSAND, i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 34), ierr, VN_SOLSAND, args, 1, layer_count, &
                        DIM_NAME_SOL, .true., i, DIM_NAME_GRU, line_count, error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message(VN_SOLCLAY, i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 35), ierr, VN_SOLCLAY, args, 1, layer_count, &
                        DIM_NAME_SOL, .true., i, DIM_NAME_GRU, line_count, error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message(VN_SOLORGCONTENT, i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 36), ierr, VN_SOLORGCONTENT, args, 1, layer_count, &
                        DIM_NAME_SOL, .true., i, DIM_NAME_GRU, line_count, error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message(VN_TSOL, i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 37), ierr, VN_TSOL, args, 1, layer_count, &
                        DIM_NAME_SOL, .true., i, DIM_NAME_GRU, line_count, error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 38), ierr, VN_TCAN, args, (layer_count + 1), (layer_count + 1), &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 39), ierr, VN_TSNO, args, (layer_count + 2), (layer_count + 2), &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 40), ierr, VN_TPND, args, (layer_count + 3), (layer_count + 3), &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message(VN_THLQSOL, i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 41), ierr, VN_THLQSOL, args, 1, layer_count, &
                        DIM_NAME_SOL, .true., i, DIM_NAME_GRU, line_count, error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 42), ierr, VN_THICSOL, args, (layer_count + 1), layer_count*2, &
                        DIM_NAME_SOL, .true., i, DIM_NAME_GRU, line_count, error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 43), ierr, VN_ZPND, args, (layer_count*2 + 1), (layer_count*2 + 1), &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                end if
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                if (ierr /= 0) then
                    error_status = ierr
                    call helper_print_error_message(VN_LQWSCAN, i, DIM_NAME_GRU, line_count, ierr)
                else
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 44), ierr, VN_LQWSCAN, args, 1, 1, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 45), ierr, VN_FZWSCAN, args, 2, 2, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 46), ierr, VN_SNO, args, 3, 3, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 47), ierr, VN_ALBSNO, args, 4, 4, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 48), ierr, VN_RHOSNO, args, 5, 5, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, ((i - 1)*39 + 49), ierr, 'GRO', args, 6, 6, &
                        check_numeric = .true., set_id = i, set_dim_name = DIM_NAME_GRU, line_number = line_count, &
                        error_status = error_status)
                end if
            end if
        end do

        !> Read and save the dates for CLASS diagnostic output.
        if (error_status == 0) then
            call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
            if (ierr /= 0) then
                error_status = ierr
                call helper_print_error_message('CLASSJOUT1', line_number = line_count, error_status = ierr)
            else
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 11), ierr, 'CLASSJOUT1', args, 1, 1, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 12), ierr, 'CLASSJOUT2', args, 2, 2, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 13), ierr, 'CLASSJAV1', args, 3, 3, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 14), ierr, 'CLASSJAV2', args, 4, 4, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
            end if
            call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
            if (ierr /= 0) then
                error_status = ierr
                call helper_print_error_message('CLASSKOUT1', line_number = line_count, error_status = ierr)
            else
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 15), ierr, 'CLASSKOUT1', args, 1, 1, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 16), ierr, 'CLASSKOUT2', args, 2, 2, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 17), ierr, 'CLASSKAV1', args, 3, 3, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 18), ierr, 'CLASSKAV2', args, 4, 4, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
            end if
        end if

        !> Read and save the forcing data start date.
        if (error_status == 0) then
            call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
            if (ierr /= 0) then
                error_status = ierr
                call helper_print_error_message('METSTARTHOUR', line_number = line_count, error_status = ierr)
            else
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 19), ierr, 'METSTARTHOUR', args, 1, 1, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 20), ierr, 'METSTARTMINUTES', args, 2, 2, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 21), ierr, 'METSTARTJDAY', args, 3, 3, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (gru_count*39 + 22), ierr, 'METSTARTYEAR', args, 4, 4, &
                    check_numeric = .true., line_number = line_count, error_status = error_status)
            end if
        end if

        !> Close unit.
        close(input_file%iunit)

        !> Check for errors.
        if (error_status /= 0) then
            call print_error("Errors occurred reading the file.")
        else

            !> Map fields.
            call create_field_maps_from_list(input_file%fields, ierr)

            !> Update the global list.
            call copy_field_list(global_list, input_file%fields, error_status)
        end if

    end subroutine

    !> Description:
    !>  Subroutine to read the 'MESH_parameters_hydrology.ini' to a field list.
    !>
    !> Variables:
    !*  global_list: Existing field list of variables read from other files (for field dependencies).
    !*  input_file: Input file.
    !*  full_path: Path to replace the default file name (optional).
    !*  error_status: Return status (0: normal).
    subroutine read_parameters_hydrology_ini_to_field_list(global_list, input_file, full_path, error_status)

        !> Input/output variables.
        type(io_field), allocatable :: global_list(:)
        type(io_file) input_file
        character(len = *), intent(in), optional :: full_path
        integer, intent(out) :: error_status

        !> Local variables.
        integer &
            sim_start_year, sim_stop_year, sim_years, frozen_soil_flag, global_iwf_flag, pbsm_flag, baseflow_flag, option_flags, &
            riverclass_variables, gru_indep_variables, gru_dep_variables, line_count, field_count, j, i, ierr
        character(len = SHORT_FIELD_LENGTH), allocatable :: gru_indep_names(:), gru_dep_names(:), args(:)
        character(len = SHORT_FIELD_LENGTH) file_version, option_name
        character(len = LONG_FIELD_LENGTH) option_value, line

        !> Reset the file structure.
        if (allocated(input_file%container)) deallocate(input_file%container)
        allocate(input_file%container, source = io_type_char_delimited())
        if (allocated(input_file%fields)) deallocate(input_file%fields)

        !> Open the file.
        if (present(full_path)) then
            input_file%full_path = adjustl(full_path)
        else
            input_file%full_path = 'MESH_parameters_hydrology.ini'
        end if
        call open_input_file(input_file, error_status = error_status)
        if (error_status /= 0) return

        !> Initialize line and field counts.
        line_count = 0
        field_count = 0

        !> Read first line to get file version.
        call read_txt_line(input_file%iunit, file_version, line_count = line_count, ierr = error_status)
        if (error_status /= 0) then
            call print_error("The file could not be read.")
        else if (index(file_version, ':') > 0) then
            file_version = file_version(1:(index(file_version, ':') - 1))
        else if (index(file_version, ' ') > 0) then
            file_version = file_version(1:(index(file_version, ' ') - 1))
        end if

        !> Look for control flags that may have dependencies in the global field list for the older file format.
        sim_start_year = NO_DATA_INT
        sim_stop_year = NO_DATA_INT
        frozen_soil_flag = NO_DATA_INT
        global_iwf_flag = NO_DATA_INT
        pbsm_flag = NO_DATA_INT
        baseflow_flag = NO_DATA_INT
        if (file_version /= '2.0') then
            do i = 1, size(global_list)

                !> Find and assign fields.
                if (global_list(i)%label == 'SIMSTARTYEAR' .and. allocated(global_list(i)%field)) then
                    call assign_field(global_list(i)%field, sim_start_year, ierr)
                    if (ierr /= 0) sim_start_year = NO_DATA_INT
                else if (global_list(i)%label == 'SIMSTOPYEAR' .and. allocated(global_list(i)%field)) then
                    call assign_field(global_list(i)%field, sim_start_year, ierr)
                    if (ierr /= 0) sim_start_year = NO_DATA_INT
                else if (global_list(i)%label == 'FROZENSOILINFILFLAG' .and. allocated(global_list(i)%field)) then
                    select type (field => global_list(i)%field)
                        type is (model_variable_char1d)
                            do j = 1, size(field%dat)
                                if (field%dat(j) == '1' .or. uppercase(field%dat(j)) == 'ON') frozen_soil_flag = 1
                            end do
                        type is (model_variable_char)
                            if (field%dat == '1' .or. uppercase(field%dat) == 'ON') frozen_soil_flag = 1
                    end select
                else if (global_list(i)%label == 'IWF' .and. allocated(global_list(i)%field)) then
                    select type (field => global_list(i)%field)
                        type is (model_variable_char1d)
                            do j = 1, size(field%dat)
                                select case (uppercase(field%dat(j)))
                                    case ('2', 'PDMROF', '3', 'LATFLOW')
                                        global_iwf_flag = 2
                                end select
                            end do
                        type is (model_variable_char)
                            select case (uppercase(field%dat))
                                case ('2', 'PDMROF', '3', 'LATFLOW')
                                    global_iwf_flag = 2
                            end select
                    end select
                else if (global_list(i)%label == 'PBSMFLAG' .and. allocated(global_list(i)%field)) then
                    select type (field => global_list(i)%field)
                        type is (model_variable_char1d)
                            do j = 1, size(field%dat)
                                if (field%dat(j) == '1' .or. uppercase(field%dat(j)) == 'ON') pbsm_flag = 1
                            end do
                        type is (model_variable_char)
                            if (field%dat == '1' .or. uppercase(field%dat) == 'ON') pbsm_flag = 1
                    end select
                else if (global_list(i)%label == 'BASEFLOWFLAG' .and. allocated(global_list(i)%field)) then
                    select type (field => global_list(i)%field)
                        type is (model_variable_char1d)
                            do j = 1, size(field%dat)
                                select case (uppercase(field%dat(j)))
                                    case ('1', 'LUO_2012')
                                        baseflow_flag = 1
                                    case ('2', 'WF_LZS')
                                        baseflow_flag = 2
                                end select
                            end do
                        type is (model_variable_char)
                            select case (uppercase(field%dat))
                                case ('1', 'LUO_2012')
                                    baseflow_flag = 1
                                case ('2', 'WF_LZS')
                                    baseflow_flag = 2
                            end select
                    end select
                end if
            end do

            !> Check for errors.
            if (frozen_soil_flag /= NO_DATA_INT .and. frozen_soil_flag /= 0) then
                if (sim_start_year == NO_DATA_INT .or. sim_stop_year == NO_DATA_INT) then
                    call print_warning( &
                        "'FROZENSOILINFILFLAG' is active but the simulation start or stop year have not been defined. " // &
                        "An error may occur reading the 'TO_ACC' parameter from this version of the file.")
                else
                    sim_years = sim_stop_year - sim_start_year + 1
                end if
            end if
        end if

        !> Option flags.
        option_flags = 0
        if (error_status == 0) then
            call read_txt_line( &
                input_file%iunit, args, n_skip_rows = 2, preserve_blank = .true., line_count = line_count, ierr = error_status)
            if (error_status == 0) then
                if (.not. size(args) >= 1) then
                    error_status = 1
                else
                    call assign_field(args(1), option_flags, error_status)
                end if
            end if
            if (error_status /= 0) then
                call helper_print_error_message('NOPTIONFLAG', line_number = line_count)
            else

                !> Store values.
                do i = 1, option_flags
                    call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                    if (ierr /= 0) then
                        call helper_print_error_message('OPTIONFLAG', i, line_number = line_count, error_status = ierr)
                    else if (.not. size(args) > 1) then
                        cycle
                    else
                        call expand_field_list(input_file%fields, increment = (size(args) - 1), error_status = error_status)
                        do j = 2, size(args)

                            !> Split the option.
                            if (index(args(j), '=') > 0) then
                                option_name = args(j)(1:(index(args(j), '=') - 1))
                                option_value = args(j)((index(args(j), '=') + 1):)
                            else
                                option_name = args(j)
                                option_value = 'on'
                            end if

                            !> Save the field.
                            call helper_assign_field( &
                                input_file%fields, field_count + (j - 1), ierr, option_name, option_value, &
                                line_number = line_count, error_status = error_status)
                        end do
                        field_count = field_count + (size(args) - 1)
                    end if
                end do
            end if
        end if

        !> River class-based variables.
        if (error_status == 0) then
            call read_txt_line( &
                input_file%iunit, args, n_skip_rows = 2, preserve_blank = .true., line_count = line_count, ierr = error_status)
            select case (file_version)
                case ('2.0')
                    riverclass_variables = 0
                    if (error_status == 0) then
                        if (.not. size(args) >= 1) then
                            error_status = 1
                        else
                            call assign_field(args(1), riverclass_variables, error_status)
                        end if
                    end if
                    if (error_status /= 0) then
                        call helper_print_error_message('NRIVERCLASSPARAMETER', line_number = line_count)
                    else
                        call expand_field_list(input_file%fields, increment = riverclass_variables, error_status = error_status)
                        do i = 1, riverclass_variables
                            call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                            if (ierr /= 0 .or. .not. size(args) > 1) then
                                call helper_print_error_message( &
                                    'RIVERCLASSPARAMETER', i, line_number = line_count, error_status = ierr)
                            else
                                call helper_assign_field( &
                                    input_file%fields, (field_count + i), ierr, args(1), args, 2, &
                                    dim_name = DIM_NAME_RVR, line_number = line_count, error_status = error_status)
                            end if
                        end do
                        field_count = field_count + riverclass_variables
                    end if
                case default
                    riverclass_variables = 1
                    call expand_field_list(input_file%fields, increment = 1, error_status = error_status)
                    call helper_assign_field( &
                        input_file%fields, (field_count + 1), ierr, 'WF_R2', args, &
                        dim_name = DIM_NAME_RVR, line_number = line_count, error_status = error_status)
                    field_count = field_count + 1
            end select
        end if

        !> GRU-independent parameters.
        gru_indep_variables = 0
        if (error_status == 0) then
            call read_txt_line( &
                input_file%iunit, args, n_skip_rows = 2, preserve_blank = .true., line_count = line_count, ierr = error_status)
            if (error_status == 0) then
                if (.not. size(args) >= 1) then
                    error_status = 1
                else
                    call assign_field(args(1), gru_indep_variables, error_status)
                end if
            end if
            if (error_status /= 0) then
                call helper_print_error_message('NGRUINDEPVARS', line_number = line_count)
            else
                call expand_field_list(input_file%fields, increment = gru_indep_variables, error_status = error_status)
                select case (file_version)
                    case ('2.0')
                        do i = 1, gru_indep_variables
                            call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                            if (ierr /= 0 .or. .not. size(args) > 1) then
                                call helper_print_error_message('GRUINDEPVARS', i, line_number = line_count, error_status = ierr)
                            else
                                call helper_assign_field( &
                                    input_file%fields, (i + option_flags + riverclass_variables), ierr, args(1), args, 2, &
                                    dim_name = DIM_NAME_B, line_number = line_count, error_status = error_status)
                            end if
                        end do
                    case default
                        allocate(gru_indep_names(gru_indep_variables))
                        do i = 1, gru_indep_variables
                            write(line, *) i
                            gru_indep_names = 'GRUINDEPVARS_' // trim(adjustl(line))
                        end do
                        j = 0
                        if (frozen_soil_flag == 1) then

                            !> FROZENSOILINFILFLAG == 1 (infiltration into frozen soils).
                            if (gru_indep_variables >= 4 + sim_years) then
                                gru_indep_names(1) = 'SOIL_POR_MAX'
                                gru_indep_names(2) = 'SOIL_DEPTH'
                                gru_indep_names(3) = 'S0'
                                gru_indep_names(4) = 'T_ICE_LENS'
                                do j = 5, (4 + sim_years)
                                    write(line, *) j - 4
                                    gru_indep_names(j) = 'T0_ACC_' // trim(adjustl(line))
                                end do
                            else
                                write(line, *) 4 + sim_years
                                call print_error( &
                                    "The section does not contain the necessary number of lines (" // trim(adjustl(line)) // &
                                    ") to define parameters with 'FROZENSOILINFILFLAG' active.")
                                error_status = 1
                            end if
                        end if
                        if (baseflow_flag == 1 .or. baseflow_flag == 2) then

                            !> BASEFLOWFLAG > 0 (lower zone storage).
                            if (gru_indep_variables >= j + 2) then
                                gru_indep_names(j + 1) = 'WRCHRGINI'
                                gru_indep_names(j + 2) = 'QBINI'
                            else
                                write(line, *) j + 2
                                call print_error( &
                                    "The section does not contain the necessary number of lines (" // trim(adjustl(line)) // &
                                    ") to define parameters with 'BASEFLOWFLAG' active.")
                                error_status = 1
                            end if
                        end if
                        if (error_status == 0) then
                            do i = 1, gru_indep_variables
                                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                                call helper_assign_field( &
                                    input_file%fields, (i + option_flags + riverclass_variables), ierr, gru_indep_names(i), args, &
                                    dim_name = DIM_NAME_B, line_number = line_count, error_status = error_status)
                            end do
                        end if
                end select
            end if
        end if

        !> GRU-dependent variables.
        gru_dep_variables = 0
        if (error_status == 0) then
            j = 2
            if (file_version /= '2.0') j = j + 1
            call read_txt_line( &
                input_file%iunit, args, n_skip_rows = j, preserve_blank = .true., line_count = line_count, ierr = error_status)
            if (error_status == 0) then
                if (.not. size(args) >= 1) then
                    error_status = 1
                else
                    call assign_field(args(1), gru_dep_variables, error_status)
                end if
            end if
            if (error_status /= 0) then
                call helper_print_error_message('NGRUDEPVARS', line_number = line_count)
            else
                call expand_field_list(input_file%fields, increment = gru_dep_variables, error_status = error_status)
                select case (file_version)
                    case ('2.0')
                        do i = 1, gru_dep_variables
                            call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                            if (ierr /= 0 .or. .not. size(args) > 1) then
                                call helper_print_error_message('GRUDEPVARS', i, line_number = line_count, error_status = ierr)
                            else
                                call helper_assign_field( &
                                    input_file%fields, (i + option_flags + riverclass_variables + gru_indep_variables), ierr, &
                                    args(1), args, 2, &
                                    dim_name = DIM_NAME_GRU, line_number = line_count, error_status = error_status)
                            end if
                        end do
                    case default
                        allocate(gru_dep_names(gru_dep_variables))
                        do i = 1, gru_dep_variables
                            write(line, *) i
                            gru_indep_names = 'GRUDEPVARS_' // trim(adjustl(line))
                        end do
                        j = 9
                        if (gru_dep_variables >= 3) then
                            gru_dep_names(1) = 'ZSNL'
                            gru_dep_names(2) = 'ZPLS'
                            gru_dep_names(3) = 'ZPLG'
                        else
                            write(line, *) 3
                            call print_error( &
                                "The section does not contain the necessary number of lines (" // trim(adjustl(line)) // &
                                ") to define 'ZSNL', 'ZPLS' and 'ZPLG').")
                            error_status = 1
                        end if
                        if (frozen_soil_flag == 1) then

                            !> FROZENSOILINFILFLAG == 1 (infiltration into frozen soils).
                            if (gru_dep_variables >= 4) then
                                gru_dep_names(4) = 'FRZC'
                            else
                                write(line, *) 4
                                call print_error( &
                                    "The section does not contain the necessary number of lines (" // trim(adjustl(line)) // &
                                    ") to define parameters with 'FROZENSOILINFILFLAG' active.")
                                error_status = 1
                            end if
                        end if
                        if (global_iwf_flag == 2 .or. global_iwf_flag == 3) then

                            !> IWF == 2 (PDMROF) or IWF == 3 (LATFLOW).
                            if (gru_dep_variables >= 9) then
                                gru_dep_names(5) = 'CMAX'
                                gru_dep_names(6) = 'CMIN'
                                gru_dep_names(7) = 'B'
                                gru_dep_names(8) = 'K1'
                                gru_dep_names(9) = 'K2'
                            else
                                write(line, *) 9
                                call print_error( &
                                    "The section does not contain the necessary number of lines (" // trim(adjustl(line)) // &
                                    ") to define parameters with 'PDMROF' (or 'LATFLOW') active.")
                                error_status = 1
                            end if
                        end if
                        if (pbsm_flag == 1) then

                            !> PBSM (blowing snow).
                            if (gru_dep_variables >= 14) then
                                gru_dep_names(10) = 'FETCH'
                                gru_dep_names(11) = 'HT'
                                gru_dep_names(12) = 'N_S'
                                gru_dep_names(13) = 'A_S'
                                gru_dep_names(14) = 'DISTRIB'
                            else
                                write(line, *) 14
                                call print_error( &
                                    "The section does not contain the necessary number of lines (" // trim(adjustl(line)) // &
                                    ") to define parameters with 'PBSMFLAG' active.")
                                error_status = 1
                            end if
                            j = j + 5
                        end if
                        if (baseflow_flag == 1 .or. baseflow_flag == 2) then
                            if (gru_dep_variables >= j + 2) then
                                if (baseflow_flag == 1) then

                                    !> BASEFLOWFLAG 1 (Luo, 2012).
                                    gru_dep_names(j + 1) = 'DGW'
                                    gru_dep_names(j + 2) = 'AGW'
                                else if (baseflow_flag == 2) then

                                    !> BASEFLOWFLAG == 2 (lower zone storage).
                                    gru_dep_names(j + 1) = 'PWR'
                                    gru_dep_names(j + 2) = 'FLZ'
                                end if
                            else
                                write(line, *) j + 2
                                call print_error( &
                                    "The section does not contain the necessary number of lines (" // trim(adjustl(line)) // &
                                    ") to define parameters with 'BASEFLOWFLAG' active.")
                                error_status = 1
                            end if
                        end if
                        if (error_status == 0) then
                            do i = 1, gru_dep_variables
                                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                                call helper_assign_field( &
                                    input_file%fields, (i + option_flags + riverclass_variables + gru_indep_variables), ierr, &
                                    gru_dep_names(i), args, &
                                    dim_name = DIM_NAME_GRU, line_number = line_count, error_status = error_status)
                            end do
                        end if
                end select
            end if
        end if

        !> Close unit.
        close(input_file%iunit)

        !> Check for errors.
        if (error_status /= 0) then
            call print_error("Errors occurred reading the file.")
        else

            !> Map fields.
            call create_field_maps_from_list(input_file%fields, ierr)

            !> Update the global list.
            call copy_field_list(global_list, input_file%fields, error_status)
        end if

    end subroutine

    !> Description:
    !>  Subroutine to read the 'MESH_input_run_options.ini' to a field list.
    !>
    !> Variables:
    !*  input_file: Input file.
    !*  full_path: Path to replace the default file name (optional).
    !*  global_list: Existing or empty field list of variables read from any files (optional).
    !*  error_status: Return status (0: normal).
    subroutine read_run_options_ini_to_field_list(input_file, full_path, global_list, error_status)

        !> Input/output variables.
        type(io_file) input_file
        character(len = *), intent(in), optional :: full_path
        type(io_field), allocatable, optional :: global_list(:)
        integer, intent(out) :: error_status

        !> Local variables.
        integer control_flags, class_output_points, line_count, field_count, j, i, ierr
        character(len = SHORT_FIELD_LENGTH), allocatable :: args(:)
        character(len = SHORT_FIELD_LENGTH) flag_name
        character(len = LONG_FIELD_LENGTH) line

        !> Reset the file structure.
        if (allocated(input_file%container)) deallocate(input_file%container)
        allocate(input_file%container, source = io_type_char_delimited())
        if (allocated(input_file%fields)) deallocate(input_file%fields)

        !> Open the file.
        if (present(full_path)) then
            input_file%full_path = adjustl(full_path)
        else
            input_file%full_path = 'MESH_input_run_options.ini'
        end if
        call open_input_file(input_file, error_status = error_status)
        if (error_status /= 0) return

        !> Initialize line and field counts.
        line_count = 0
        field_count = 0

        !> Control flags.
        control_flags = 0
        if (error_status == 0) then
            call read_txt_line( &
                input_file%iunit, args, n_skip_rows = 3, preserve_blank = .true., line_count = line_count, ierr = ierr)
            if (ierr == 0) then
                if (.not. size(args) >= 1) then
                    ierr = 1
                else
                    call assign_field(args(1), control_flags, ierr)
                end if
            end if
            if (ierr /= 0) then
                call helper_print_error_message('NCONTROLFLAG', line_number = line_count)
                error_status = ierr
            else

                !> Store values.
                call expand_field_list(input_file%fields, increment = control_flags, error_status = ierr)
                field_count = field_count + control_flags
                do i = 1, control_flags
                    call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                    if (size(args) >= 1) then
                        flag_name = adjustl(args(1))
                    else
                        flag_name = 'CONTROLFLAG'
                    end if
                    call helper_assign_field( &
                        input_file%fields, i, ierr, flag_name, args, 2, line_number = line_count, error_status = error_status)
                end do
            end if
        end if

        !> CLASS output points.
        class_output_points = 0
        j = 0
        if (error_status == 0) then
            call read_txt_line( &
                input_file%iunit, args, n_skip_rows = 2, preserve_blank = .true., line_count = line_count, ierr = ierr)
            if (ierr == 0) then
                if (.not. size(args) >= 1) then
                    ierr = 1
                else
                    call assign_field(args(1), class_output_points, ierr)
                end if
            end if
            if (ierr /= 0) then
                call helper_print_error_message('NCLASSOUTPOINTS', line_number = line_count)
                error_status = ierr
            else if (class_output_points > 0) then
                call expand_field_list(input_file%fields, increment = 3, error_status = ierr)
                field_count = field_count + 3
                call read_txt_line(input_file%iunit, args, n_skip_rows = 1, line_count = line_count, ierr = ierr)
                call helper_assign_field( &
                    input_file%fields, (field_count - 2), ierr, 'CLASSOUTGRIDID', args, line_number = line_count, &
                    error_status = error_status)
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                call helper_assign_field( &
                    input_file%fields, (field_count - 1), ierr, 'CLASSOUTGRUID', args, line_number = line_count, &
                    error_status = error_status)
                call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
                call helper_assign_field( &
                    input_file%fields, field_count, ierr, 'CLASSOUTFOLDER', args, line_number = line_count, &
                    error_status = error_status)
            else
                j = j + 4
            end if
        end if

        !> Generic output folder.
        if (error_status == 0) then
            call read_txt_line( &
                input_file%iunit, args, n_skip_rows = (j + 2), preserve_blank = .true., line_count = line_count, ierr = ierr)
            if (ierr == 0) then
                if (.not. size(args) >= 1) ierr = 1
            end if
            if (ierr /= 0) then
                call helper_print_error_message('OUTPUTDIR', line_number = line_count, error_status = ierr)
                error_status = ierr
            else
                call expand_field_list(input_file%fields, increment = 1, error_status = ierr)
                field_count = field_count + 1
                call helper_assign_field( &
                    input_file%fields, field_count, ierr, 'OUTPUTDIR', args(1), &
                    line_number = line_count, error_status = error_status)
            end if
        end if

        !> Simulation start date.
        if (error_status == 0) then
            call read_txt_line( &
                input_file%iunit, args, n_skip_rows = 2, preserve_blank = .true., line_count = line_count, ierr = ierr)
            if (ierr == 0) then
                if (.not. size(args) >= 4) ierr = 1
            end if
            if (ierr /= 0) then
                call helper_print_error_message('SIMSTARTTIMES', line_number = line_count, error_status = ierr)
                error_status = ierr
            else
                call expand_field_list(input_file%fields, increment = 4, error_status = ierr)
                field_count = field_count + 4
                call helper_assign_field( &
                    input_file%fields, (field_count - 3), ierr, 'SIMSTARTYEAR', args(1), line_number = line_count, &
                    error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (field_count - 2), ierr, 'SIMSTARTJDAY', args(2), line_number = line_count, &
                    error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (field_count - 1), ierr, 'SIMSTARTHOUR', args(3), line_number = line_count, &
                    error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, field_count, ierr, 'SIMSTARTMINUTES', args(4), line_number = line_count, &
                    error_status = error_status)
            end if
        end if

        !> Simulation stop date.
        if (error_status == 0) then
            call read_txt_line(input_file%iunit, args, line_count = line_count, ierr = ierr)
            if (ierr == 0) then
                if (.not. size(args) >= 4) ierr = 1
            end if
            if (ierr /= 0) then
                call helper_print_error_message('SIMSTOPTIMES', line_number = line_count, error_status = ierr)
                error_status = ierr
            else
                call expand_field_list(input_file%fields, increment = 4, error_status = ierr)
                field_count = field_count + 4
                call helper_assign_field( &
                    input_file%fields, (field_count - 3), ierr, 'SIMSTOPYEAR', args(1), line_number = line_count, &
                    error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (field_count - 2), ierr, 'SIMSTOPJDAY', args(2), line_number = line_count, &
                    error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, (field_count - 1), ierr, 'SIMSTOPHOUR', args(3), line_number = line_count, &
                    error_status = error_status)
                call helper_assign_field( &
                    input_file%fields, field_count, ierr, 'SIMSTOPMINUTES', args(4), line_number = line_count, &
                    error_status = error_status)
            end if
        end if

        !> Close unit.
        close(input_file%iunit)

        !> Check for errors.
        if (error_status /= 0) then
            call print_error("Errors occurred reading the file.")
        else

            !> Map fields.
            call create_field_maps_from_list(input_file%fields, ierr)

            !> Update the global list.
            if (present(global_list)) call copy_field_list(global_list, input_file%fields, error_status)
        end if

    end subroutine

!> *********************************************************************
!> Open and read in values from MESH_input_drainage_database.txt file
!>   if new_shd.r2c file was not found
!> *********************************************************************
!        open(UNIT=20, FILE='MESH_input_drainage_database.txt',
!     &    STATUS='OLD', IOSTAT=ierr)
!        if (ierr == 0) then
!          print *, 'Reading Drainage Database from ',
!     &      'MESH_input_drainage_database.txt'
!        else
!          print *, 'MESH_input_drainage_database.txt not found'
!          stop
!        end if
!        read(20, '(i5, 50x, i5)') shd%NA, shd%NAA
!        read(20, '(f10.0, 5x, 2i5)') shd%AL, shd%NRVR, shd%lc%NTYPE
!        GRDN = 0.0
!        GRDE = 0.0

!todo change the code or documentation on these variables.
!ANDY Set ILG from the read-in values
!        ILG = shd%NA*shd%lc%NTYPE

!> Using IOSTAT allows us to try to read input that may or may not exist.
!> If all of the values successfully get read, IOSTAT=VarName will set
!> VarName to 0. If all of the values were not successfully read,
!> VarName would be set to 1 or more. In this case, the VarName that
!> we are using is ierr.
!        read(20, '(12i5, 2f5.0)', IOSTAT=ierr) IYMIN, WF_IYMAX,
!     &    JXMIN, WF_JXMAX, LATDEGMIN, LATMINMIN, LATDEGMAX, LATMINMAX,
!     &    LONDEGMIN, LONMINMIN, LONDEGMAX, LONMINMAX, GRDN, GRDE

!> Condition for Lat/Long by Frank S Sept/1999
!        if (GRDN > 0.0) then
!          IYMIN = LATDEGMIN*60 + LATMINMIN
!          WF_IYMAX = LATDEGMAX*60 + LATMINMAX
!          JXMIN = LONDEGMIN*60 + LONMINMIN
!          WF_JXMAX = LONDEGMAX*60 + LONMINMAX

!        else
!> Define GRDN & GRDE for UTM
!          GRDN = shd%AL/1000.0
!          GRDE = shd%AL/1000.0
!        end if
!        read(20, '(2i5)') YCOUNT, XCOUNT

!        do i = 1, YCOUNT
!          read(20, *)
!        end do

!        do i = 1, shd%NA
!          read(20, '(5x, 2i5, 3f10.5, i7, 5i5, f5.2, 15f5.2)') YYY(i),
!     &      XXX(i), WF_DA(i), WF_BNKFLL(i), WF_CHANNELSLOPE(i),
!     &      WF_ELEV(i), WF_IBN(i), WF_IROUGH(i), WF_ICHNL(i),
!     &      WF_NEXT(i), WF_IREACH(i), FRAC(i),
!     &      (ACLASS(i, j), j = 1, shd%lc%NTYPE)
!> check to make sure land cover areas sum to 100%
!          WF_LAND_COUNT = 1
!          WF_LAND_MAX = 0.0
!          WF_LAND_SUM = 0.0
!          do j = 1, shd%lc%NTYPE
!            WF_LAND_SUM = WF_LAND_SUM + ACLASS(i, j)
!            if (ACLASS(i, j) > WF_LAND_MAX) then
!              WF_LAND_COUNT = j
!              WF_LAND_MAX = ACLASS(i, j)
!            end if
!          end do
!          if (WF_LAND_SUM /= 1.0) THEN
!            ACLASS(i, WF_LAND_COUNT) =
!     &        ACLASS(i, WF_LAND_COUNT) - (WF_LAND_SUM - 1.0)
!          end if
!        end do

!        close(20)

end module
