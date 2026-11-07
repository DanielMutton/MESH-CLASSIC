module input_fields

    !> 'mesh_io': For I/O field and file types, options and routines.
    use mesh_io

    implicit none

    !* input_files: Local list of input files.
    type(io_file), dimension(:), allocatable :: basin_info_files, parameter_files

    !* input_field_list: Local list of fields read from the input files.
    type(io_field), allocatable :: input_field_list(:)

    interface get_field_value
        module procedure get_field_value_real2d
        module procedure get_field_value_int2d
        module procedure get_field_value_real1d
        module procedure get_field_value_int1d
        module procedure get_field_value_char1d
        module procedure get_field_value_scalar
        module procedure get_field_value_from_list_ids_real1d
        module procedure get_field_value_from_list_ids_int1d
        module procedure get_field_value_from_list_ids_char1d
        module procedure get_field_value_from_list_ids_scalar
        module procedure get_field_value_from_list_real1d
        module procedure get_field_value_from_list_int1d
        module procedure get_field_value_from_list_char1d
        module procedure get_field_value_from_list_scalar
    end interface

    interface get_cell_values_from_list
        module procedure get_cell_values_from_list_real1d
        module procedure get_cell_values_from_list_int1d
        module procedure get_cell_values_from_list_char1d
    end interface

    interface get_tile_values_from_list
        module procedure get_tile_values_from_list_real1d
        module procedure get_tile_values_from_list_int1d
        module procedure get_tile_values_from_list_char1d
    end interface

    interface check_field_assigned
        module procedure check_field_assigned_real3d_map_real3d
        module procedure check_field_assigned_real2d_map_real2d
        module procedure check_field_assigned_real2d
        module procedure check_field_assigned_int2d_map_int2d
        module procedure check_field_assigned_int2d
        module procedure check_field_assigned_real1d_map_real1d
        module procedure check_field_assigned_real1d
        module procedure check_field_assigned_int1d_map_int1d
        module procedure check_field_assigned_int1d
        module procedure check_field_assigned_char1d_map_char1d
        module procedure check_field_assigned_char1d
        module procedure check_field_assigned_real
        module procedure check_field_assigned_int
        module procedure check_field_assigned_char
    end interface

    interface check_values
        module procedure check_values_real4d
        module procedure check_values_int4d
        module procedure check_values_real3d
        module procedure check_values_int3d
        module procedure check_values_real2d
        module procedure check_values_int2d
        module procedure check_values_real1d
        module procedure check_values_int1d
        module procedure check_values_char1d
        module procedure check_values_real
        module procedure check_values_int
        module procedure check_values_char
    end interface

    type field_list_search_names
        character(len = SHORT_FIELD_LENGTH), allocatable :: names(:)
    end type

    type field_list_search_result
        character(len = SHORT_FIELD_LENGTH) :: label = ''
        integer, allocatable :: ids(:)
    end type

    contains

    subroutine expand_file_list(file_list, increment, error_status)

        !> Input/output variables.
        type(io_file), allocatable :: file_list(:)
        integer, intent(in), optional :: increment
        integer, intent(out) :: error_status

        !> Local variables.
        type(io_file), allocatable :: temp_list(:)
        integer i

        !> Transfer the increment to a local variable (if present).
        i = 1
        if (present(increment)) i = max(increment, 1)

        !> Status.
        error_status = 0

        !> Transfer fields.
        if (.not. allocated(file_list)) then

            !> No need to transfer if the list is not yet allocated.
            allocate(file_list(i))
        else

            !> Allocate the temporary list.
            allocate(temp_list(size(file_list) + i))

            !> Transfer fields.
            temp_list(1:size(file_list)) = file_list

            !> Deallocate the source list.
            deallocate(file_list)

            !> Transfer the allocation back to the source list.
            call move_alloc(from = temp_list, to = file_list)
        end if

    end subroutine

    function check_for_file_in_list(file_list, label, full_path) result(i)

        !> Input/output variables.
        type(io_file), allocatable :: file_list(:)
        character(len = *), intent(in) :: label
        character(len = *), intent(in), optional :: full_path

        !> Local variables.
        integer j, i, ierr
        logical ltest

        !> Return key.
        i = 0

        !> Check if the file exists in the list by label.
        if (allocated(file_list)) then
            do j = 1, size(file_list)
                ltest = .false.
                if (present(full_path)) then
                    if (file_list(j)%full_path == full_path) ltest = .true.
                end if
                if (file_list(j)%label == label) ltest = .true.
                if (ltest) then

                    !> The file already exists.
                    i = j
                    ltest = .false.
                    exit
                end if
            end do
            if (i == 0) then

                !> File does not exist.
                i = size(file_list) + 1
                ltest = .true.
            end if
        else

            !> The list is not allocated so the file cannot exist.
            i = 1
            ltest = .true.
        end if

        !> Create an entry if the file was added to the list.
        if (ltest) then
            call expand_file_list(file_list, error_status = ierr)
            if (ierr == 0) then
                file_list(i)%label = adjustl(label)
                if (present(full_path)) file_list(i)%full_path = adjustl(full_path)
            else

                !> Unable to update the list.
                i = 0
            end if
        end if

    end function

    subroutine get_dim_map(input_field, desired_map_order, dim_map)

        !> Modules.
        use model_variables, only: vs

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        integer, intent(in) :: desired_map_order(:)
        integer, allocatable, intent(out) :: dim_map(:, :)

        !> Local variables.
        integer i

        !> Check if dimension names are defined.
        if (.not. allocated(input_field%dim_names)) return

        !> Check the desired dimension order.
        associate(mapped_dim_order => input_field%mapping%mapped_dim_order)
            if (size(desired_map_order) == 1) then
                select case (desired_map_order(1))
                    case (MAP_ORDER_N)
                        if (associated(vs%grid)) then

                            !> Allocate the map.
                            allocate(dim_map(size(input_field%dim_names), vs%grid%dim_length))
                            dim_map = 1

                            !> Create map.
                            if (mapped_dim_order(MAP_ORDER_X) > 0 .and. mapped_dim_order(MAP_ORDER_Y) > 0) then
                                dim_map(mapped_dim_order(MAP_ORDER_X), :) = vs%grid%from_grid_xy(1, :)
                                dim_map(mapped_dim_order(MAP_ORDER_Y), :) = vs%grid%from_grid_xy(2, :)
                            else if (mapped_dim_order(MAP_ORDER_M) > 0) then
                                dim_map(mapped_dim_order(MAP_ORDER_M), :) = vs%grid%from_gru
                            else if (mapped_dim_order(MAP_ORDER_K) > 0) then
                                dim_map(mapped_dim_order(MAP_ORDER_K), :) = vs%grid%from_riverclass
                            else if (mapped_dim_order(MAP_ORDER_N) > 0) then
                                do i = 1, vs%grid%dim_length
                                    dim_map(mapped_dim_order(MAP_ORDER_N), i) = i
                                end do
                            else if (mapped_dim_order(MAP_ORDER_B) > 0) then
                                dim_map(mapped_dim_ORDER(MAP_ORDER_B), :) = 1
                            else
                                deallocate(dim_map)
                            end if
                        end if
                    case(MAP_ORDER_G)
                        if (associated(vs%tile)) then

                            !> Allocate the map.
                            allocate(dim_map(size(input_field%dim_names), vs%tile%dim_length))
                            dim_map = 1

                            !> Create map.
                            if (mapped_dim_order(MAP_ORDER_X) > 0 .and. mapped_dim_order(MAP_ORDER_Y) > 0) then
                                dim_map(mapped_dim_order(MAP_ORDER_X), :) = vs%tile%from_grid_xy(1, :)
                                dim_map(mapped_dim_order(MAP_ORDER_Y), :) = vs%tile%from_grid_xy(2, :)
                            else if (mapped_dim_order(MAP_ORDER_M) > 0) then
                                dim_map(mapped_dim_order(MAP_ORDER_M), :) = vs%tile%from_gru
                            else if (mapped_dim_order(MAP_ORDER_N) > 0) then
                                dim_map(mapped_dim_order(MAP_ORDER_N), :) = vs%tile%from_cell
                            else if (mapped_dim_order(MAP_ORDER_B) > 0) then
                                dim_map(mapped_dim_ORDER(MAP_ORDER_B), :) = 1
                            else if (mapped_dim_order(MAP_ORDER_G) > 0) then
                                do i = 1, vs%tile%dim_length
                                    dim_map(mapped_dim_order(MAP_ORDER_G), i) = i
                                end do
                            else
                                deallocate(dim_map)
                            end if
                        end if
                end select
            end if

            !> Check for a level ID.
            if (allocated(dim_map) .and. input_field%level_id > 0 .and. mapped_dim_order(MAP_ORDER_L) > 0) then
                dim_map(mapped_dim_order(MAP_ORDER_L), :) = input_field%level_id
            end if
        end associate

    end subroutine

    subroutine check_values_real4d(input_field, output_values, const_mul, const_add, valid_max, valid_min, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        real, allocatable :: output_values(:, :, :, :)
        real, intent(in), optional :: const_mul
        real, intent(in), optional :: const_add
        real, intent(in), optional :: valid_max
        real, intent(in), optional :: valid_min
        integer, intent(out) :: error_status

        !> Local variables.
        real const_mul_value, const_add_value, valid_max_value, valid_min_value, no_data_value
        integer z
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message

        !> Assign the field.
        call assign_field(input_field, output_values, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            const_mul_value = 1.0
            const_add_value = 0.0
            valid_max_value = NO_DATA_REAL
            valid_min_value = -NO_DATA_REAL
            no_data_value = NO_DATA_REAL
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%const_mul, const_mul_value, z)
                    if (z /= 0) const_mul_value = 1.0
                    call assign_field(field%const_add, const_add_value, z)
                    if (z /= 0) const_add_value = 0.0
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_REAL
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_REAL
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
                class is (model_variable_int_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_REAL
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_REAL
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
            end select

            !> Apply overrides.
            if (present(const_mul)) const_mul_value = const_mul
            if (present(const_add)) const_add_value = const_add
            if (present(valid_max)) valid_max_value = valid_max
            if (present(valid_min)) valid_min_value = valid_min

            !> Apply transformations.
            if (const_mul_value /= 1.0) then
                write(code, *) const_mul_value
                call print_remark( &
                    "Adjusting the values of the variable '" // trim(input_field%label) // &
                    "' by the multiplicative factor " // trim(adjustl(code)) // ".")
                output_values = output_values*const_mul_value
            end if
            if (const_add_value /= 0.0) then
                write(code, *) const_add_value
                call print_remark( &
                    "Adjusting the values of the variable '" // trim(input_field%label) // &
                    "' by the additive factor " // trim(adjustl(code)) // ".")
                output_values = output_values + const_add_value
            end if

            !> Check for 'no data' values.
            if (any(output_values == no_data_value)) then
                write(code, *) count(output_values == no_data_value)
                call print_warning( &
                    "The variable '" // trim(input_field%label) // "' contains (" // &
                    trim(adjustl(code)) // ") 'no data' values.")
            end if

            !> Check range.
            if (any(output_values > valid_max_value)) then
                write(code, *) count(output_values > valid_max_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_max_value
                message = message // ") values greater than the valid maximum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
            if (any(output_values < valid_min_value)) then
                write(code, *) count(output_values < valid_min_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_min_value
                message = message // ") values less than the valid minimum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
        end if

    end subroutine

    subroutine check_values_int4d(input_field, output_values, valid_max, valid_min, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        integer, allocatable :: output_values(:, :, :, :)
        integer, intent(in), optional :: valid_max
        integer, intent(in), optional :: valid_min
        integer, intent(out) :: error_status

        !> Local variables.
        integer valid_max_value, valid_min_value, no_data_value, z
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message

        !> Assign the field.
        call assign_field(input_field, output_values, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            valid_max_value = NO_DATA_INT
            valid_min_value = -NO_DATA_INT
            no_data_value = NO_DATA_INT
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_INT
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_INT
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
                class is (model_variable_int_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_INT
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_INT
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
            end select

            !> Apply overrides.
            if (present(valid_max)) valid_max_value = valid_max
            if (present(valid_min)) valid_min_value = valid_min

            !> Check for 'no data' values.
            if (any(output_values == no_data_value)) then
                write(code, *) count(output_values == no_data_value)
                call print_warning( &
                    "The variable '" // trim(input_field%label) // "' contains (" // &
                    trim(adjustl(code)) // ") 'no data' values.")
            end if

            !> Check range.
            if (any(output_values > valid_max_value)) then
                write(code, *) count(output_values > valid_max_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_max_value
                message = message // ") values greater than the valid maximum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
            if (any(output_values < valid_min_value)) then
                write(code, *) count(output_values < valid_min_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_min_value
                message = message // ") values less than the valid minimum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
        end if

    end subroutine

    subroutine check_values_real3d(input_field, output_values, const_mul, const_add, valid_max, valid_min, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        real, allocatable :: output_values(:, :, :)
        real, intent(in), optional :: const_mul
        real, intent(in), optional :: const_add
        real, intent(in), optional :: valid_max
        real, intent(in), optional :: valid_min
        integer, intent(out) :: error_status

        !> Local variables.
        real const_mul_value, const_add_value, valid_max_value, valid_min_value, no_data_value
        integer z
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message

        !> Assign the field.
        call assign_field(input_field, output_values, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            const_mul_value = 1.0
            const_add_value = 0.0
            valid_max_value = NO_DATA_REAL
            valid_min_value = -NO_DATA_REAL
            no_data_value = NO_DATA_REAL
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%const_mul, const_mul_value, z)
                    if (z /= 0) const_mul_value = 1.0
                    call assign_field(field%const_add, const_add_value, z)
                    if (z /= 0) const_add_value = 0.0
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_REAL
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_REAL
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
                class is (model_variable_int_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_REAL
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_REAL
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
            end select

            !> Apply overrides.
            if (present(const_mul)) const_mul_value = const_mul
            if (present(const_add)) const_add_value = const_add
            if (present(valid_max)) valid_max_value = valid_max
            if (present(valid_min)) valid_min_value = valid_min

            !> Apply transformations.
            if (const_mul_value /= 1.0) then
                write(code, *) const_mul_value
                call print_remark( &
                    "Adjusting the values of the variable '" // trim(input_field%label) // &
                    "' by the multiplicative factor " // trim(adjustl(code)) // ".")
                output_values = output_values*const_mul_value
            end if
            if (const_add_value /= 0.0) then
                write(code, *) const_add_value
                call print_remark( &
                    "Adjusting the values of the variable '" // trim(input_field%label) // &
                    "' by the additive factor " // trim(adjustl(code)) // ".")
                output_values = output_values + const_add_value
            end if

            !> Check for 'no data' values.
            if (any(output_values == no_data_value)) then
                write(code, *) count(output_values == no_data_value)
                call print_warning( &
                    "The variable '" // trim(input_field%label) // "' contains (" // &
                    trim(adjustl(code)) // ") 'no data' values.")
            end if

            !> Check range.
            if (any(output_values > valid_max_value)) then
                write(code, *) count(output_values > valid_max_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_max_value
                message = message // ") values greater than the valid maximum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
            if (any(output_values < valid_min_value)) then
                write(code, *) count(output_values < valid_min_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_min_value
                message = message // ") values less than the valid minimum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
        end if

    end subroutine

    subroutine check_values_int3d(input_field, output_values, valid_max, valid_min, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        integer, allocatable :: output_values(:, :, :)
        integer, intent(in), optional :: valid_max
        integer, intent(in), optional :: valid_min
        integer, intent(out) :: error_status

        !> Local variables.
        integer valid_max_value, valid_min_value, no_data_value, z
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message

        !> Assign the field.
        call assign_field(input_field, output_values, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            valid_max_value = NO_DATA_INT
            valid_min_value = -NO_DATA_INT
            no_data_value = NO_DATA_INT
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_INT
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_INT
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
                class is (model_variable_int_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_INT
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_INT
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
            end select

            !> Apply overrides.
            if (present(valid_max)) valid_max_value = valid_max
            if (present(valid_min)) valid_min_value = valid_min

            !> Check for 'no data' values.
            if (any(output_values == no_data_value)) then
                write(code, *) count(output_values == no_data_value)
                call print_warning( &
                    "The variable '" // trim(input_field%label) // "' contains (" // &
                    trim(adjustl(code)) // ") 'no data' values.")
            end if

            !> Check range.
            if (any(output_values > valid_max_value)) then
                write(code, *) count(output_values > valid_max_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_max_value
                message = message // ") values greater than the valid maximum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
            if (any(output_values < valid_min_value)) then
                write(code, *) count(output_values < valid_min_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_min_value
                message = message // ") values less than the valid minimum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
        end if

    end subroutine

    subroutine check_values_real2d(input_field, output_values, const_mul, const_add, valid_max, valid_min, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        real, allocatable :: output_values(:, :)
        real, intent(in), optional :: const_mul
        real, intent(in), optional :: const_add
        real, intent(in), optional :: valid_max
        real, intent(in), optional :: valid_min
        integer, intent(out) :: error_status

        !> Local variables.
        real const_mul_value, const_add_value, valid_max_value, valid_min_value, no_data_value
        integer z
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message

        !> Assign the field.
        call assign_field(input_field, output_values, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            const_mul_value = 1.0
            const_add_value = 0.0
            valid_max_value = NO_DATA_REAL
            valid_min_value = -NO_DATA_REAL
            no_data_value = NO_DATA_REAL
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%const_mul, const_mul_value, z)
                    if (z /= 0) const_mul_value = 1.0
                    call assign_field(field%const_add, const_add_value, z)
                    if (z /= 0) const_add_value = 0.0
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_REAL
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_REAL
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
                class is (model_variable_int_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_REAL
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_REAL
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
            end select

            !> Apply overrides.
            if (present(const_mul)) const_mul_value = const_mul
            if (present(const_add)) const_add_value = const_add
            if (present(valid_max)) valid_max_value = valid_max
            if (present(valid_min)) valid_min_value = valid_min

            !> Apply transformations.
            if (const_mul_value /= 1.0) then
                write(code, *) const_mul_value
                call print_remark( &
                    "Adjusting the values of the variable '" // trim(input_field%label) // &
                    "' by the multiplicative factor " // trim(adjustl(code)) // ".")
                output_values = output_values*const_mul_value
            end if
            if (const_add_value /= 0.0) then
                write(code, *) const_add_value
                call print_remark( &
                    "Adjusting the values of the variable '" // trim(input_field%label) // &
                    "' by the additive factor " // trim(adjustl(code)) // ".")
                output_values = output_values + const_add_value
            end if

            !> Check for 'no data' values.
            if (any(output_values == no_data_value)) then
                write(code, *) count(output_values == no_data_value)
                call print_warning( &
                    "The variable '" // trim(input_field%label) // "' contains (" // &
                    trim(adjustl(code)) // ") 'no data' values.")
            end if

            !> Check range.
            if (any(output_values > valid_max_value)) then
                write(code, *) count(output_values > valid_max_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_max_value
                message = message // ") values greater than the valid maximum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
            if (any(output_values < valid_min_value)) then
                write(code, *) count(output_values < valid_min_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_min_value
                message = message // ") values less than the valid minimum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
        end if

    end subroutine

    subroutine check_values_int2d(input_field, output_values, valid_max, valid_min, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        integer, allocatable :: output_values(:, :)
        integer, intent(in), optional :: valid_max
        integer, intent(in), optional :: valid_min
        integer, intent(out) :: error_status

        !> Local variables.
        integer valid_max_value, valid_min_value, no_data_value, z
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message

        !> Assign the field.
        call assign_field(input_field, output_values, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            valid_max_value = NO_DATA_INT
            valid_min_value = -NO_DATA_INT
            no_data_value = NO_DATA_INT
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_INT
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_INT
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
                class is (model_variable_int_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_INT
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_INT
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
            end select

            !> Apply overrides.
            if (present(valid_max)) valid_max_value = valid_max
            if (present(valid_min)) valid_min_value = valid_min

            !> Check for 'no data' values.
            if (any(output_values == no_data_value)) then
                write(code, *) count(output_values == no_data_value)
                call print_warning( &
                    "The variable '" // trim(input_field%label) // "' contains (" // &
                    trim(adjustl(code)) // ") 'no data' values.")
            end if

            !> Check range.
            if (any(output_values > valid_max_value)) then
                write(code, *) count(output_values > valid_max_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_max_value
                message = message // ") values greater than the valid maximum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
            if (any(output_values < valid_min_value)) then
                write(code, *) count(output_values < valid_min_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_min_value
                message = message // ") values less than the valid minimum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
        end if

    end subroutine

    subroutine check_values_real1d(input_field, output_values, const_mul, const_add, valid_max, valid_min, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        real, allocatable :: output_values(:)
        real, intent(in), optional :: const_mul
        real, intent(in), optional :: const_add
        real, intent(in), optional :: valid_max
        real, intent(in), optional :: valid_min
        integer, intent(out) :: error_status

        !> Local variables.
        real const_mul_value, const_add_value, valid_max_value, valid_min_value, no_data_value
        integer z
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message

        !> Assign the field.
        call assign_field(input_field, output_values, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            const_mul_value = 1.0
            const_add_value = 0.0
            valid_max_value = NO_DATA_REAL
            valid_min_value = -NO_DATA_REAL
            no_data_value = NO_DATA_REAL
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%const_mul, const_mul_value, z)
                    if (z /= 0) const_mul_value = 1.0
                    call assign_field(field%const_add, const_add_value, z)
                    if (z /= 0) const_add_value = 0.0
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_REAL
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_REAL
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
                class is (model_variable_int_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_REAL
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_REAL
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
            end select

            !> Apply overrides.
            if (present(const_mul)) const_mul_value = const_mul
            if (present(const_add)) const_add_value = const_add
            if (present(valid_max)) valid_max_value = valid_max
            if (present(valid_min)) valid_min_value = valid_min

            !> Apply transformations.
            if (const_mul_value /= 1.0) then
                write(code, *) const_mul_value
                call print_remark( &
                    "Adjusting the values of the variable '" // trim(input_field%label) // &
                    "' by the multiplicative factor " // trim(adjustl(code)) // ".")
                output_values = output_values*const_mul_value
            end if
            if (const_add_value /= 0.0) then
                write(code, *) const_add_value
                call print_remark( &
                    "Adjusting the values of the variable '" // trim(input_field%label) // &
                    "' by the additive factor " // trim(adjustl(code)) // ".")
                output_values = output_values + const_add_value
            end if

            !> Check for 'no data' values.
            if (any(output_values == no_data_value)) then
                write(code, *) count(output_values == no_data_value)
                call print_warning( &
                    "The variable '" // trim(input_field%label) // "' contains (" // &
                    trim(adjustl(code)) // ") 'no data' values.")
            end if

            !> Check range.
            if (any(output_values > valid_max_value)) then
                write(code, *) count(output_values > valid_max_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_max_value
                message = message // ") values greater than the valid maximum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
            if (any(output_values < valid_min_value)) then
                write(code, *) count(output_values < valid_min_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_min_value
                message = message // ") values less than the valid minimum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
        end if

    end subroutine

    subroutine check_values_int1d(input_field, output_values, valid_max, valid_min, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        integer, allocatable :: output_values(:)
        integer, intent(in), optional :: valid_max
        integer, intent(in), optional :: valid_min
        integer, intent(out) :: error_status

        !> Local variables.
        integer valid_max_value, valid_min_value, no_data_value, z
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message

        !> Assign the field.
        call assign_field(input_field, output_values, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            valid_max_value = NO_DATA_INT
            valid_min_value = -NO_DATA_INT
            no_data_value = NO_DATA_INT
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_INT
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_INT
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
                class is (model_variable_int_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_INT
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_INT
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
            end select

            !> Apply overrides.
            if (present(valid_max)) valid_max_value = valid_max
            if (present(valid_min)) valid_min_value = valid_min

            !> Check for 'no data' values.
            if (any(output_values == no_data_value)) then
                write(code, *) count(output_values == no_data_value)
                call print_warning( &
                    "The variable '" // trim(input_field%label) // "' contains (" // &
                    trim(adjustl(code)) // ") 'no data' values.")
            end if

            !> Check range.
            if (any(output_values > valid_max_value)) then
                write(code, *) count(output_values > valid_max_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_max_value
                message = message // ") values greater than the valid maximum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
            if (any(output_values < valid_min_value)) then
                write(code, *) count(output_values < valid_min_value)
                message = "The variable '" // trim(input_field%label) // "' contains (" // trim(adjustl(code))
                write(code, *) valid_min_value
                message = message // ") values less than the valid minimum value of " // trim(adjustl(code)) // "."
                call print_warning(message)
            end if
        end if

    end subroutine

    subroutine check_values_char1d(input_field, output_values, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        character(len = *), allocatable :: output_values(:)
        integer, intent(out) :: error_status

        !> Local variables.
        integer z
        character(len = SHORT_FIELD_LENGTH) no_data_value
        character(len = DEFAULT_FIELD_LENGTH) code

        !> Assign the field.
        call assign_field(input_field, output_values, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            no_data_value = NO_DATA_CHAR
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_CHAR
                class is (model_variable_int_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_CHAR
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_CHAR
            end select

            !> Check for 'no data' values.
            if (any(output_values == no_data_value)) then
                write(code, *) count(output_values == no_data_value)
                call print_warning( &
                    "The variable '" // trim(input_field%label) // "' contains (" // &
                    trim(adjustl(code)) // ") 'no data' values.")
            end if
        end if

    end subroutine

    subroutine check_values_real(input_field, output_value, const_mul, const_add, valid_max, valid_min, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        real output_value
        real, intent(in), optional :: const_mul
        real, intent(in), optional :: const_add
        real, intent(in), optional :: valid_max
        real, intent(in), optional :: valid_min
        integer, intent(out) :: error_status

        !> Local variables.
        real const_mul_value, const_add_value, valid_max_value, valid_min_value, no_data_value
        integer z
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message

        !> Assign the field.
        call assign_field(input_field, output_value, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            const_mul_value = 1.0
            const_add_value = 0.0
            valid_max_value = NO_DATA_REAL
            valid_min_value = -NO_DATA_REAL
            no_data_value = NO_DATA_REAL
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%const_mul, const_mul_value, z)
                    if (z /= 0) const_mul_value = 1.0
                    call assign_field(field%const_add, const_add_value, z)
                    if (z /= 0) const_add_value = 0.0
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_REAL
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_REAL
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
                class is (model_variable_int_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_REAL
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_REAL
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_REAL
            end select

            !> Apply overrides.
            if (present(const_mul)) const_mul_value = const_mul
            if (present(const_add)) const_add_value = const_add
            if (present(valid_max)) valid_max_value = valid_max
            if (present(valid_min)) valid_min_value = valid_min

            !> Apply transformations.
            if (const_mul_value /= 1.0) then
                write(code, *) const_mul_value
                call print_remark( &
                    "Adjusting the values of the variable '" // trim(input_field%label) // &
                    "' by the multiplicative factor " // trim(adjustl(code)) // ".")
                output_value = output_value*const_mul_value
            end if
            if (const_add_value /= 0.0) then
                write(code, *) const_add_value
                call print_remark( &
                    "Adjusting the values of the variable '" // trim(input_field%label) // &
                    "' by the additive factor " // trim(adjustl(code)) // ".")
                output_value = output_value + const_add_value
            end if

            !> Check for 'no data' values.
            if (output_value == no_data_value) then
                call print_warning("The variable '" // trim(input_field%label) // "' value is a 'no data' value.")
            end if

            !> Check range.
            if (output_value > valid_max_value) then
                write(code, *) valid_max_value
                message = &
                    "The variable '" // trim(input_field%label) // "' value is greater than the valid maximum value of " // &
                    trim(adjustl(code)) // "."
                call print_warning(message)
            end if
            if (output_value < valid_min_value) then
                write(code, *) valid_min_value
                message = &
                    "The variable '" // trim(input_field%label) // "' value is less than the valid minimum value of " // &
                    trim(adjustl(code)) // "."
                call print_warning(message)
            end if
        end if

    end subroutine

    subroutine check_values_int(input_field, output_value, valid_max, valid_min, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        integer output_value
        integer, intent(in), optional :: valid_max
        integer, intent(in), optional :: valid_min
        integer, intent(out) :: error_status

        !> Local variables.
        integer valid_max_value, valid_min_value, no_data_value, z
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message

        !> Assign the field.
        call assign_field(input_field, output_value, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            valid_max_value = NO_DATA_INT
            valid_min_value = -NO_DATA_INT
            no_data_value = NO_DATA_INT
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_INT
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_INT
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
                class is (model_variable_int_base)
                    call assign_field(field%valid_max, valid_max_value, z)
                    if (z /= 0) valid_max_value = NO_DATA_INT
                    call assign_field(field%valid_min, valid_min_value, z)
                    if (z /= 0) valid_min_value = -NO_DATA_INT
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_INT
            end select

            !> Apply overrides.
            if (present(valid_max)) valid_max_value = valid_max
            if (present(valid_min)) valid_min_value = valid_min

            !> Check for 'no data' values.
            if (output_value == no_data_value) then
                call print_warning("The variable '" // trim(input_field%label) // "' value is a 'no data' value.")
            end if

            !> Check range.
            if (output_value > valid_max_value) then
                write(code, *) valid_max_value
                message = &
                    "The variable '" // trim(input_field%label) // "' value is greater than the valid maximum value of " // &
                    trim(adjustl(code)) // "."
                call print_warning(message)
            end if
            if (output_value < valid_min_value) then
                write(code, *) valid_min_value
                message = &
                    "The variable '" // trim(input_field%label) // "' value is less than the valid minimum value of " // &
                    trim(adjustl(code)) // "."
                call print_warning(message)
            end if
        end if

    end subroutine

    subroutine check_values_char(input_field, output_value, error_status)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        character(len = *) output_value
        integer, intent(out) :: error_status

        !> Local variables.
        integer z
        character(len = SHORT_FIELD_LENGTH) no_data_value
        character(len = DEFAULT_FIELD_LENGTH) code

        !> Assign the field.
        call assign_field(input_field, output_value, error_status)
        if (error_status /= 0) then
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Check for field attributes.
            no_data_value = NO_DATA_CHAR
            z = 0
            select type (field => input_field%field)
                class is (model_variable_real_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_CHAR
                class is (model_variable_int_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_CHAR
                class is (model_variable_char_base)
                    call assign_field(field%no_data_value, no_data_value, z)
                    if (z /= 0) no_data_value = NO_DATA_CHAR
            end select

            !> Check for 'no data' values.
            if (output_value == no_data_value) then
                call print_warning("The variable '" // trim(input_field%label) // "' value is a 'no data' value.")
            end if
        end if

    end subroutine

    subroutine get_field_value_real2d( &
        input_field, output_values, output_dim_names, level_map, value_mask, &
        const_mul, const_add, valid_max, valid_min)

        !> Modules.
        use model_variables, only: vs
        use parse_utilities, only: to_lowercase

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        real, allocatable :: output_values(:, :)
        character(len = *), intent(in) :: output_dim_names(2)
        character(len = *), intent(in), optional :: level_map(:)
        logical, intent(in), optional :: value_mask(:)
        real, intent(in), optional :: const_mul
        real, intent(in), optional :: const_add
        real, intent(in), optional :: valid_max
        real, intent(in), optional :: valid_min

        !> Local variables.
        real, allocatable :: &
            dat4(:, :, :, :), dat4_mapped(:, :, :, :), dat3(:, :, :), dat3_mapped(:, :, :), dat2(:, :), dat2_mapped(:, :), &
            dat2_temp(:, :), dat1(:)
        real dat
        integer, allocatable :: source_dim1_map(:), source_dim2_map(:), source_dim3_map(:), mapped_dim_order(:)
        integer mappable_dim_count, size1, size2, level_id, j2, j1, j, i_match, i_max, i, z, ierr
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message
        logical, allocatable :: active_value_mask(:)

        !> Check if the dimensions of the field are not defined.
        if (.not. allocated(input_field%dim_names)) then

            !> Try to assign the field.
            call check_values(input_field, output_values, const_mul, const_add, valid_max, valid_min, ierr)
        else if (any(input_field%dim_names == output_dim_names(1)) .and. any(input_field%dim_names == output_dim_names(2))) then

            !> Assign variable to temporary array.
            call check_values(input_field, dat2, const_mul, const_add, valid_max, valid_min, ierr)
            if (ierr == 0) then

                !> Get the order of the dimensions.
                call get_dimension_order(input_field%dim_names, output_dim_names, mapped_dim_order, ierr)
                if (ierr /= 0) then
                    ierr = 3
                else

                    !> Map the field.
                    call map_dimensions(dat2, mapped_dim_order, output_values, ierr)
                end if
            end if
        else if (allocated(input_field%mapping%mapped_dim_order)) then

            !> Check if the field is associated with a specific level ID.
            level_id = 0
            ierr = 0
            if (len_trim(input_field%level) > 0) then
                if (input_field%level_id > 0) then
                    level_id = input_field%level_id
                else if (present(level_map)) then
                    do i = 1, size(level_map)
                        if (to_lowercase(level_map(i)) == to_lowercase(input_field%level)) then
                            level_id = i
                            exit
                        end if
                    end do
                    if (level_id == 0) ierr = 1
                end if
                if (ierr /= 0) then
                    call print_warning( &
                        "Unabled to determine an index for the level '" // trim(input_field%level) // &
                        "' to assign the variable '" // trim(input_field%label) // "'.")
                    return
                end if
            end if

            !> Map the variable.
            associate(mapped_orders => input_field%mapping%mapped_dim_order)

                !> Save the number of mappable dimensions.
                mappable_dim_count = count(mapped_orders > 0)
                if (mapped_orders(MAP_ORDER_T) > 0) mappable_dim_count = mappable_dim_count - 1

                !> Derive dimensions for the output variable.
                size1 = 0
                size2 = 0
                if (any(output_dim_names == DIM_NAME_X) .and. any(output_dim_names == DIM_NAME_Y)) then

                    !> Check if the variable is allocated.
                    if (.not. allocated(output_values)) then
                        if (vs%grid_x > 0 .and. vs%grid_y > 0) then
                            size1 = vs%grid_x
                            size2 = vs%grid_y
                        else if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                            size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_X))
                            size2 = input_field%dim_lengths(mapped_orders(MAP_ORDER_Y))
                        else if (mapped_orders(MAP_ORDER_N) > 0) then
                            size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_N))
                            size2 = 1
                        end if
                    end if
                else

                    !> Derive a size for the level dimension of the variable in case it should be allocated or re-allocated.
                    if (level_id > 0 .and. ( &
                        len_trim(input_field%level_dim_name) == 0 .or. any(DIM_NAMES_OF_L == input_field%level_dim_name))) then
                        size2 = level_id
                    else if (mapped_orders(MAP_ORDER_L) > 0) then
                        size2 = input_field%dim_lengths(mapped_orders(MAP_ORDER_L))
                    else
                        size2 = 1
                    end if

                    !> Derive a size for the elemental dimension.
                    if (any(DIM_NAMES_OF_N == output_dim_names(1))) then
                        if (vs%cell_count > 0) then
                            size1 = vs%cell_count
                        else if (mapped_orders(MAP_ORDER_N) > 0) then
                            size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_N))
                        end if
                    else if (any(DIM_NAMES_OF_G == output_dim_names(1))) then
                        if (vs%landtile_count > 0) then
                            size1 = vs%landtile_count
                        else if (mapped_orders(MAP_ORDER_G) > 0) then
                            size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_G))
                        end if
                    end if
                end if

                !> Check the allocation status of the output variable.
                ierr = 0
                if (size1 > 0 .and. size2 > 0) then
                    if (.not. allocated(output_values)) then

                        !> Allocate the variable.
                        call allocate_field(output_values, (/size1, size2/), ierr)
                    else
                        if (size(output_values, 1) < size1) then

                            !> Expand the first dimension of the variable.
                            call allocate_field(dat2_temp, (/size1, size(output_values, 2)/), ierr)
                            if (ierr == 0) then
                                dat2_temp(1:size(output_values, 1), :) = output_values
                                call move_alloc(dat2_temp, output_values)
                            end if
                        end if
                        if (size(output_values, 2) < size2) then

                            !> Expand the second dimension of the variable.
                            call allocate_field(dat2_temp, (/size(output_values, 1), size2/), ierr)
                            if (ierr == 0) then
                                dat2_temp(:, 1:size(output_values, 2)) = output_values
                                call move_alloc(dat2_temp, output_values)
                            end if
                        end if
                    end if
                end if
                if (ierr /= 0) then
                    call print_warning( &
                        "An allocation error occurred trying to assign the variable '" // trim(input_field%label) // ". " // &
                        "The variable cannot be assigned.")
                else if (.not. allocated(output_values)) then
                    call print_warning( &
                        "Unable to allocate the field to the desired dimensions '" // trim(output_dim_names(1)) // &
                        ", " // trim(output_dim_names(2)) // "'. " // &
                        "The variable '" // trim(input_field%label) // "' cannot be assigned.")
                    ierr = 1
                end if
                if (ierr /= 0) return

                !> Derive maps for the output variable.
                if (any(output_dim_names == DIM_NAME_X) .and. any(output_dim_names == DIM_NAME_Y)) then
                    call allocate_field(source_dim1_map, (/size1/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/size2/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        do i = 1, size1
                            source_dim1_map(i) = i
                        end do
                        do i = 1, size2
                            source_dim2_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        do i = 1, size1
                            source_dim1_map(i) = i
                        end do
                        source_dim2_map = size2
                    else
                        ierr = 3
                    end if
                else if (any(DIM_NAMES_OF_N == output_dim_names(1)) .and. associated(vs%grid)) then
                    call allocate_field(source_dim1_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%grid%from_grid_x
                        source_dim2_map = vs%grid%from_grid_y
                        if (mapped_orders(MAP_ORDER_K) > 0) then
                            call allocate_field(source_dim3_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                            source_dim3_map = vs%tile%from_riverclass
                        else
                            ierr = 3
                        end if
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        do i = 1, vs%cell_count
                            source_dim1_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_K) > 0 .or. any(DIM_NAMES_OF_K == input_field%level_dim_name)) then
                        source_dim1_map = vs%grid%from_riverclass
                    else
                        ierr = 3
                    end if
                else if (any(DIM_NAMES_OF_G == output_dim_names(1)) .and. associated(vs%tile)) then
                    call allocate_field(source_dim1_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%tile%from_grid_x
                        source_dim2_map = vs%tile%from_grid_y
                        if (mapped_orders(MAP_ORDER_M) > 0) then
                            call allocate_field(source_dim3_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                            source_dim3_map = vs%tile%from_gru
                        else
                            ierr = 3
                        end if
                    else if (mapped_orders(MAP_ORDER_G) > 0) then
                        do i = 1, vs%landtile_count
                            source_dim1_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_N) > 0 .and. mapped_orders(MAP_ORDER_M) > 0) then
                        source_dim1_map = vs%tile%from_cell
                        source_dim2_map = vs%tile%from_gru
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        source_dim1_map = vs%tile%from_cell
                    else if (mapped_orders(MAP_ORDER_M) > 0 .or. any(DIM_NAMES_OF_M == input_field%level_dim_name)) then
                        source_dim1_map = vs%tile%from_gru
                    else
                        ierr = 3
                    end if
                else
                    ierr = 3
                end if
                if (ierr /= 0 .and. ierr /= 3) then
                    call print_warning( &
                        "An error occurred allocating the mapping variables for the variable '" // trim(input_field%label) // "'.")
                    return
                end if

                !> Determine indices to assign the output variable.
                i_match = 0
                i_max = size(output_values, 1)
                j1 = 1
                j2 = size(output_values, 2)
                if (level_id > 0) then
                    if (len_trim(input_field%level_dim_name) == 0 .or. any(DIM_NAMES_OF_L == input_field%level_dim_name)) then

                        !> Check the level against the dimension.
                        if (level_id < j1 .or. level_id > j2) then
                            write(code, *) level_id
                            message = &
                                "The index [" // trim(adjustl(code)) // "] to assign for the variable '" // &
                                trim(input_field%label) // "' is outside of the range of active indices ["
                            write(code, *) j1
                            message = message // trim(adjustl(code)) // ":"
                            write(code, *) j2
                            message = &
                                message // trim(adjustl(code)) // "] in the second dimension. The variable cannot be assigned."
                            call print_warning(message)
                            return
                        else
                            j1 = level_id
                            j2 = level_id
                        end if
                    else if (any(source_dim1_map == NO_DATA_INT)) then

                        !> No mapping is active when a non-level mapping index is defined.
                        write(code, *) level_id
                        call print_warning( &
                            "A mapping index [" // trim(adjustl(code)) // "] is defined to assign variable '" // &
                            trim(input_field%label) // "' from the '" // trim(input_field%level_dim_name) // &
                            "' dimension but no mapping could be determined for the output variable to the '" // &
                            trim(output_dim_names(1)) // "' dimension. The variable cannot be assigned.")
                        return
                    else if (.not. any(source_dim1_map == level_id)) then

                        !> The map contains no indices that would map to the mapping index.
                        write(code, *) level_id
                        call print_warning( &
                            "The variable '" // trim(input_field%label) // "' contains a mapping index [" // &
                            trim(adjustl(code)) // "] that does not map to any active element of the output variable in the '" // &
                            trim(output_dim_names(1)) // "' dimension. No values from this variable will be assigned.")
                        return
                    else
                        i_match = level_id
                    end if
                end if

                !> Value mask for distributed fields.
                allocate(active_value_mask(i_max))
                active_value_mask = .true.
                if (present(value_mask)) then
                    if (sum(mapped_orders) == 2 .and. mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        if (size(value_mask) == size(active_value_mask)) active_value_mask = value_mask
                        if (any(.not. active_value_mask)) then
                            write(code, *) count(.not. active_value_mask)
                            call print_info( &
                                "A filter is active to exclude " // trim(adjustl(code)) // " values from being assigned " // &
                                "from the distributed input variable '" // trim(input_field%label) // "'. " // &
                                "The values already assigned from other non-distributed inputs will presist at these locations.")
                        end if
                    end if
                end if

                !> Assign the field.
                if (.not. allocated(input_field%dim_lengths)) then

                    !> Scalar.
                    call check_values(input_field, dat, const_mul, const_add, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        do j = j1, j2
                            if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                where (source_dim1_map == i_match) output_values(:, j) = dat
                            else
                                output_values(:, j) = dat
                            end if
                        end do
                    end if
                else if (mappable_dim_count == 4) then

                    !> Check for compatible dimensions.
                    allocate(mapped_dim_order(4))
                    mapped_dim_order = 0
                    if ( &
                        all(source_dim1_map /= NO_DATA_INT) .and. all(source_dim2_map /= NO_DATA_INT) .and. &
                        all(source_dim3_map /= NO_DATA_INT) .and. &
                        mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0 .and. mapped_orders(MAP_ORDER_L) > 0) &
                        then
                        if (mapped_orders(MAP_ORDER_K) > 0) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_K, DIM_NAME_L/), mapped_dim_order, ierr)
                        else if (mapped_orders(MAP_ORDER_M) > 0) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_M, DIM_NAME_L/), mapped_dim_order, ierr)
                        end if
                    end if
                    if (any(mapped_dim_order == 0)) then

                        !> Can only be mapped from the specified dimension combinations.
                        ierr = 3
                    else

                        !> Assign the values to a temporary variable.
                        call check_values(input_field, dat4, const_mul, const_add, valid_max, valid_min, ierr)
                        if (ierr == 0) then

                            !> Ensure the dimensions of the field are in the expected mapped order.
                            call map_dimensions(dat4, mapped_dim_order, dat4_mapped, ierr)

                            !> Check the dimensions and assign values.
                            if (j1 > size(dat4_mapped, 4)) then

                                !> The values to assign are outside the range of the level dimension.
                                write(code, *) j1
                                message = &
                                    "The first index to assign to the model variable [" // trim(adjustl(code)) // "] " // &
                                    "is higher than the number of values in the level dimension of the variable '" // &
                                    trim(input_field%label) // "' ("
                                write(code, *) size(dat4_mapped, 4)
                                message = message // trim(adjustl(code)) // "). The variable cannot be assigned."
                                call print_warning(message)
                                ierr = 1
                            else if (j2 > size(dat4_mapped, 4)) then

                                !> The source contains lesser values than expected.
                                write(code, *) size(dat4_mapped, 4)
                                message = &
                                    "The variable '" // trim(input_field%label) // "' contains " // trim(adjustl(code)) // &
                                    " values in the level dimension where "
                                write(code, *) j2
                                message = message // trim(adjustl(code)) // " are expected. The value at index ["
                                write(code, *) size(dat4_mapped, 4)
                                message = &
                                    message // trim(adjustl(code)) // &
                                    "] will be used to assign higher indices of the assigned variable."
                                call print_warning(message)

                                !> Assign up to the available values.
                                forall (i = 1:i_max, j = j1:size(dat4_mapped, 4), active_value_mask(i))
                                    output_values(i, j) = &
                                        dat4_mapped(source_dim1_map(i), source_dim2_map(i), source_dim3_map(i), j)
                                end forall

                                !> Assign the last value to higher indices in the dimension.
                                forall (i = 1:i_max, j = (size(dat4_mapped, 4) + 1):j2, active_value_mask(i))
                                    output_values(i, j) = &
                                        dat4_mapped( &
                                            source_dim1_map(i), source_dim2_map(i), source_dim3_map(i), size(dat4_mapped, 4))
                                end forall
                            else
                                if (j2 < size(dat4_mapped, 4)) then

                                    !> The source contains more values than are required.
                                    write(code, *) size(dat4_mapped, 4)
                                    message = &
                                        "The variable '" // trim(input_field%label) // "' contains " // &
                                        trim(adjustl(code)) // " values in the level dimension where "
                                    write(code, *) j2
                                    message = &
                                        message // trim(adjustl(code)) // " are expected. The values of the variable " // &
                                        "at higher indices will not be assigned and will be ignored."
                                    call print_warning(message)
                                end if

                                !> Assign the values.
                                forall (i = 1:i_max, j = j1:j2, active_value_mask(i))
                                    output_values(i, j) = &
                                        dat4_mapped(source_dim1_map(i), source_dim2_map(i), source_dim3_map(i), j)
                                end forall
                            end if
                        end if
                    end if
                else if (mappable_dim_count == 3) then

                    !> Check for compatible dimensions.
                    if (all(source_dim1_map /= NO_DATA_INT) .and. all(source_dim2_map /= NO_DATA_INT)) then
                        if (mapped_orders(MAP_ORDER_K) > 0 .and. all(source_dim3_map /= NO_DATA_INT)) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_K/), mapped_dim_order, ierr)
                        else if (mapped_orders(MAP_ORDER_M) > 0 .and. all(source_dim3_map /= NO_DATA_INT)) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_M/), mapped_dim_order, ierr)
                        else if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_L/), mapped_dim_order, ierr)
                        else if (mapped_orders(MAP_ORDER_N) > 0 .and. mapped_orders(MAP_ORDER_M) > 0) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_N, DIM_NAME_M, DIM_NAME_L/), mapped_dim_order, ierr)
                        else
                            allocate(mapped_dim_order(3))
                            mapped_dim_order = 0
                        end if
                        if (any(mapped_dim_order == 0)) then

                            !> Can only be mapped to 2-D by cartesian grid.
                            ierr = 3
                        else

                            !> Assign the values to a temporary variable.
                            call check_values(input_field, dat3, const_mul, const_add, valid_max, valid_min, ierr)
                            if (ierr == 0) then

                                !> Ensure the dimensions of the field are in the expected mapped order.
                                call map_dimensions(dat3, mapped_dim_order, dat3_mapped, ierr)

                                !> Check the dimensions and assign values.
                                if (j1 > size(dat3_mapped, 3)) then

                                    !> The values to assign are outside the range of the level dimension.
                                    write(code, *) j1
                                    message = &
                                        "The first index to assign to the model variable [" // trim(adjustl(code)) // "] " // &
                                        "is higher than the number of values in the level dimension of the variable '" // &
                                        trim(input_field%label) // "' ("
                                    write(code, *) size(dat3_mapped, 3)
                                    message = message // trim(adjustl(code)) // "). The variable cannot be assigned."
                                    call print_warning(message)
                                    ierr = 1
                                else if (j2 > size(dat3_mapped, 3)) then

                                    !> The source contains lesser values than expected.
                                    write(code, *) size(dat3_mapped, 3)
                                    message = &
                                        "The variable '" // trim(input_field%label) // "' contains " // trim(adjustl(code)) // &
                                        " values in the level dimension where "
                                    write(code, *) j2
                                    message = message // trim(adjustl(code)) // " are expected. The value at index ["
                                    write(code, *) size(dat3_mapped, 3)
                                    message = &
                                        message // trim(adjustl(code)) // &
                                        "] will be used to assign higher indices of the assigned variable."
                                    call print_warning(message)

                                    !> Assign up to the available values.
                                    forall (i = 1:i_max, j = j1:size(dat3_mapped, 3), active_value_mask(i))
                                        output_values(i, j) = dat3_mapped(source_dim1_map(i), source_dim2_map(i), j)
                                    end forall

                                    !> Assign the last value to higher indices in the dimension.
                                    forall (i = 1:i_max, j = (size(dat3_mapped, 3) + 1):j2, active_value_mask(i))
                                        output_values(i, j) = &
                                            dat3_mapped(source_dim1_map(i), source_dim2_map(i), size(dat3_mapped, 3))
                                    end forall
                                else
                                    if (j2 < size(dat3_mapped, 3)) then

                                        !> The source contains more values than are required.
                                        write(code, *) size(dat3_mapped, 3)
                                        message = &
                                            "The variable '" // trim(input_field%label) // "' contains " // &
                                            trim(adjustl(code)) // " values in the level dimension where "
                                        write(code, *) j2
                                        message = &
                                            message // trim(adjustl(code)) // " are expected. The values of the variable " // &
                                            "at higher indices will not be assigned and will be ignored."
                                        call print_warning(message)
                                    end if

                                    !> Assign the values.
                                    forall (i = 1:i_max, j = j1:j2, active_value_mask(i))
                                        output_values(i, j) = dat3_mapped(source_dim1_map(i), source_dim2_map(i), j)
                                    end forall
                                end if
                            end if
                        end if
                    else
                        ierr = 3
                    end if
                else if (mappable_dim_count == 2) then

                    !> Assign the values to a temporary variable.
                    call check_values(input_field, dat2, const_mul, const_add, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT)) then
                            if (all(source_dim2_map /= NO_DATA_INT)) then

                                !> Get the order of the dimensions.
                                if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                                    call get_dimension_order( &
                                        input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y/), mapped_dim_order, ierr)
                                else if (mapped_orders(MAP_ORDER_N) > 0 .and. mapped_orders(MAP_ORDER_M) > 0) then
                                    call get_dimension_order( &
                                        input_field%dim_names, (/DIM_NAME_N, DIM_NAME_M/), mapped_dim_order, ierr)
                                else
                                    allocate(mapped_dim_order(2))
                                    mapped_dim_order = 0
                                end if
                                if (any(mapped_dim_order == 0)) then
                                    ierr = 3
                                else

                                    !> Map the field.
                                    call map_dimensions(dat2, mapped_dim_order, dat2_mapped, ierr)
                                    if (ierr /= 0) then
                                        ierr = 3
                                    else

                                        !> Assign the values.
                                        forall (i = 1:i_max, j = j1:j2, active_value_mask(i))
                                            output_values(i, j) = dat2_mapped(source_dim1_map(i), source_dim2_map(i))
                                        end forall
                                    end if
                                end if
                            else if (mapped_orders(MAP_ORDER_L) > 0) then

                                !> Get the order of the dimensions.
                                z = 0
                                if (mapped_orders(MAP_ORDER_N) > 0) then
                                    call get_dimension_order(input_field%dim_names, (/DIM_NAME_N, DIM_NAME_L/), mapped_dim_order, z)
                                else if (mapped_orders(MAP_ORDER_K) > 0) then
                                    call get_dimension_order(input_field%dim_names, (/DIM_NAME_K, DIM_NAME_L/), mapped_dim_order, z)
                                else if (mapped_orders(MAP_ORDER_G) > 0) then
                                    call get_dimension_order(input_field%dim_names, (/DIM_NAME_G, DIM_NAME_L/), mapped_dim_order, z)
                                else if (mapped_orders(MAP_ORDER_M) > 0) then
                                    call get_dimension_order(input_field%dim_names, (/DIM_NAME_M, DIM_NAME_L/), mapped_dim_order, z)
                                else
                                    ierr = 3
                                end if
                                if (z /= 0) then
                                    ierr = 3
                                else

                                    !> Map the field.
                                    call map_dimensions(dat2, mapped_dim_order, dat2_mapped, ierr)
                                    if (ierr /= 0) then
                                        ierr = 3
                                    else

                                        !> Check the dimensions and assign values.
                                        if (j1 > size(dat2_mapped, 2)) then

                                            !> The values to assign are outside the range of the level dimension.
                                            write(code, *) j1
                                            message = &
                                                "The first index to assign to the model variable [" // trim(adjustl(code)) // &
                                                "] is larger than the size of the level dimension of the variable '" // &
                                                trim(input_field%label) // "' ("
                                            write(code, *) size(dat2_mapped, 2)
                                            message = message // trim(adjustl(code)) // "). The variable cannot be assigned."
                                            call print_warning(message)
                                            ierr = 1
                                        else if (j2 > size(dat2_mapped, 2)) then

                                            !> The source contains lesser values than expected.
                                            write(code, *) size(dat2_mapped, 2)
                                            message = &
                                                "The variable '" // trim(input_field%label) // "' contains " // &
                                                trim(adjustl(code)) // " values in the level dimension where "
                                            write(code, *) j2
                                            message = message // trim(adjustl(code)) // " are expected. The value at index ["
                                            write(code, *) size(dat2_mapped, 2)
                                            message = &
                                                message // trim(adjustl(code)) // &
                                                "] will be used to assign values at higher indices of the assigned variable."
                                            call print_warning(message)

                                            !> Assign up to the available values.
                                            do j = j1, size(dat2_mapped, 2)
                                                if (i_match > 0) then
                                                    forall (i = 1:i_max, source_dim1_map(i) == i_match)
                                                        output_values(i, j) = dat2_mapped(i_match, j)
                                                    end forall
                                                else
                                                    do i = 1, i_max
                                                        output_values(i, j) = dat2_mapped(source_dim1_map(i), j)
                                                    end do
                                                end if
                                            end do

                                            !> Assign the last value to higher indices in the dimension.
                                            do j = (size(dat2_mapped, 2) + 1), j2
                                                if (i_match > 0) then
                                                    forall (i = 1:i_max, source_dim1_map(i) == i_match)
                                                        output_values(i, j) = dat2_mapped(i_match, size(dat2_mapped, 2))
                                                    end forall
                                                else
                                                    do i = 1, i_max
                                                        output_values(i, j) = dat2_mapped(source_dim1_map(i), size(dat2_mapped, 2))
                                                    end do
                                                end if
                                            end do
                                        else
                                            if (j2 < size(dat2_mapped, 2)) then

                                                !> The source contains more values than are required.
                                                write(code, *) size(dat2_mapped, 2)
                                                message = &
                                                    "The variable '" // trim(input_field%label) // "' contains " // &
                                                    trim(adjustl(code)) // " values in the level dimension where "
                                                write(code, *) j2
                                                message = &
                                                    message // trim(adjustl(code)) // " are expected. " // &
                                                    "Values at higher indices in the source variable will be ignored."
                                                call print_warning(message)
                                            end if

                                            !> Assign the values.
                                            do j = j1, j2
                                                if (i_match > 0) then
                                                    forall (i = 1:i_max, source_dim1_map(i) == i_match)
                                                        output_values(i, j) = dat2_mapped(i_match, j)
                                                    end forall
                                                else
                                                    forall (i = 1:i_max) output_values(i, j) = dat2_mapped(source_dim1_map(i), j)
                                                end if
                                            end do
                                        end if
                                    end if
                                end if
                            else
                                ierr = 3
                            end if
                        end if
                    end if
                else if (mappable_dim_count == 1 .or. sum(input_field%dim_lengths) == 1) then

                    !> Assign variable to temporary array.
                    call check_values(input_field, dat1, const_mul, const_add, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT) .and. mapped_orders(MAP_ORDER_L) == 0) then

                            !> Assign mapped values.
                            if (i_match > 0) then
                                forall (i = 1:i_max, j = j1:j2, source_dim1_map(i) == i_match) output_values(i, j) = dat1(i_match)
                            else
                                forall (i = 1:i_max, j = j1:j2) output_values(i, j) = dat1(i)
                            end if
                        else

                            !> Check the dimensions and assign level values.
                            if (j1 > size(dat1)) then

                                !> The values to assign are outside the range of the level dimension.
                                write(code, *) j1
                                message = &
                                    "The first index to assign to the model variable [" // trim(adjustl(code)) // "] " // &
                                    "is higher than the number of values in the level dimension of the variable '" // &
                                    trim(input_field%label) // "' ("
                                write(code, *) size(dat1)
                                message = message // trim(adjustl(code)) // "). The variable cannot be assigned."
                                call print_warning(message)
                                ierr = 1
                            else if (j2 > size(dat1)) then

                                !> The source contains lesser values than expected.
                                write(code, *) size(dat1)
                                message = &
                                    "The variable '" // trim(input_field%label) // "' contains " // trim(adjustl(code)) // &
                                    " values in the level dimension where "
                                write(code, *) j2
                                message = message // trim(adjustl(code)) // " are expected. The value at index ["
                                write(code, *) size(dat1)
                                message = &
                                    message // trim(adjustl(code)) // &
                                    "] will be used to assign higher indices of the assigned variable."
                                call print_warning(message)

                                !> Assign up to the available values.
                                do j = j1, size(dat1)
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i, j) = dat1(j)
                                    else
                                        output_values(:, j) = dat1(j)
                                    end if
                                end do

                                !> Assign the last value to higher indices in the dimension.
                                do j = (size(dat1) + 1), j2
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i, j) = dat1(size(dat1))
                                    else
                                        output_values(:, j) = dat1(size(dat1))
                                    end if
                                end do
                            else
                                if (j2 < size(dat1)) then

                                    !> The source contains more values than are required.
                                    write(code, *) size(dat1)
                                    message = &
                                        "The variable '" // trim(input_field%label) // "' contains " // &
                                        trim(adjustl(code)) // " values in the level dimension where "
                                    write(code, *) j2
                                    message = &
                                        message // trim(adjustl(code)) // " are expected. The values of the variable " // &
                                        "at higher indices will not be assigned and will be ignored."
                                    call print_warning(message)
                                end if

                                !> Assign the values.
                                do j = j1, j2
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i, j) = dat1(j)
                                    else
                                        output_values(:, j) = dat1(j)
                                    end if
                                end do
                            end if
                        end if
                    end if
                else
                    ierr = 3
                end if
            end associate
        end if

        !> Check for errors (case '1' is handled inline and ignored).
        select case (ierr)
            case (3)
                call print_warning("An error occurred mapping the variable '" // trim(input_field%label) // "'.")
            case (0)
                call print_info("Assigned the variable '" // trim(input_field%label) // "'.")
        end select

    end subroutine

    subroutine get_field_value_int2d( &
        input_field, output_values, output_dim_names, level_map, value_mask, &
        valid_max, valid_min)

        !> Modules.
        use model_variables, only: vs
        use parse_utilities, only: to_lowercase

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        integer, allocatable :: output_values(:, :)
        character(len = *), intent(in) :: output_dim_names(2)
        character(len = *), intent(in), optional :: level_map(:)
        logical, intent(in), optional :: value_mask(:)
        integer, intent(in), optional :: valid_max
        integer, intent(in), optional :: valid_min

        !> Local variables.
        integer, allocatable :: &
            dat4(:, :, :, :), dat4_mapped(:, :, :, :), dat3(:, :, :), dat3_mapped(:, :, :), dat2(:, :), dat2_mapped(:, :), &
            dat2_temp(:, :), dat1(:), source_dim1_map(:), source_dim2_map(:), source_dim3_map(:), mapped_dim_order(:)
        integer dat, mappable_dim_count, size1, size2, level_id, j2, j1, j, i_match, i_max, i, z, ierr
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message
        logical, allocatable :: active_value_mask(:)

        !> Check if the dimensions of the field are not defined.
        if (.not. allocated(input_field%dim_names)) then

            !> Try to assign the field.
            call check_values(input_field, output_values, valid_max, valid_min, ierr)
        else if (any(input_field%dim_names == output_dim_names(1)) .and. any(input_field%dim_names == output_dim_names(2))) then

            !> Assign variable to temporary array.
            call check_values(input_field, dat2, valid_max, valid_min, ierr)
            if (ierr == 0) then

                !> Get the order of the dimensions.
                call get_dimension_order(input_field%dim_names, output_dim_names, mapped_dim_order, ierr)
                if (ierr /= 0) then
                    ierr = 3
                else

                    !> Map the field.
                    call map_dimensions(dat2, mapped_dim_order, output_values, ierr)
                end if
            end if
        else if (allocated(input_field%mapping%mapped_dim_order)) then

            !> Check if the field is associated with a specific level ID.
            level_id = 0
            ierr = 0
            if (len_trim(input_field%level) > 0) then
                if (input_field%level_id > 0) then
                    level_id = input_field%level_id
                else if (present(level_map)) then
                    do i = 1, size(level_map)
                        if (to_lowercase(level_map(i)) == to_lowercase(input_field%level)) then
                            level_id = i
                            exit
                        end if
                    end do
                    if (level_id == 0) ierr = 1
                end if
                if (ierr /= 0) then
                    call print_warning( &
                        "Unabled to determine an index for the level '" // trim(input_field%level) // &
                        "' to assign the variable '" // trim(input_field%label) // "'.")
                    return
                end if
            end if

            !> Map the variable.
            associate(mapped_orders => input_field%mapping%mapped_dim_order)

                !> Save the number of mappable dimensions.
                mappable_dim_count = count(mapped_orders > 0)
                if (mapped_orders(MAP_ORDER_T) > 0) mappable_dim_count = mappable_dim_count - 1

                !> Derive dimensions for the output variable.
                size1 = 0
                size2 = 0
                if (any(output_dim_names == DIM_NAME_X) .and. any(output_dim_names == DIM_NAME_Y)) then

                    !> Check if the variable is allocated.
                    if (.not. allocated(output_values)) then
                        if (vs%grid_x > 0 .and. vs%grid_y > 0) then
                            size1 = vs%grid_x
                            size2 = vs%grid_y
                        else if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                            size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_X))
                            size2 = input_field%dim_lengths(mapped_orders(MAP_ORDER_Y))
                        else if (mapped_orders(MAP_ORDER_N) > 0) then
                            size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_N))
                            size2 = 1
                        end if
                    end if
                else

                    !> Derive a size for the level dimension of the variable in case it should be allocated or re-allocated.
                    if (level_id > 0 .and. ( &
                        len_trim(input_field%level_dim_name) == 0 .or. any(DIM_NAMES_OF_L == input_field%level_dim_name))) then
                        size2 = level_id
                    else if (mapped_orders(MAP_ORDER_L) > 0) then
                        size2 = input_field%dim_lengths(mapped_orders(MAP_ORDER_L))
                    else
                        size2 = 1
                    end if

                    !> Derive a size for the elemental dimension.
                    if (any(DIM_NAMES_OF_N == output_dim_names(1))) then
                        if (vs%cell_count > 0) then
                            size1 = vs%cell_count
                        else if (mapped_orders(MAP_ORDER_N) > 0) then
                            size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_N))
                        end if
                    else if (any(DIM_NAMES_OF_G == output_dim_names(1))) then
                        if (vs%landtile_count > 0) then
                            size1 = vs%landtile_count
                        else if (mapped_orders(MAP_ORDER_G) > 0) then
                            size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_G))
                        end if
                    end if
                end if

                !> Check the allocation status of the output variable.
                ierr = 0
                if (size1 > 0 .and. size2 > 0) then
                    if (.not. allocated(output_values)) then

                        !> Allocate the variable.
                        call allocate_field(output_values, (/size1, size2/), ierr)
                    else
                        if (size(output_values, 1) < size1) then

                            !> Expand the first dimension of the variable.
                            call allocate_field(dat2_temp, (/size1, size(output_values, 2)/), ierr)
                            if (ierr == 0) then
                                dat2_temp(1:size(output_values, 1), :) = output_values
                                call move_alloc(dat2_temp, output_values)
                            end if
                        end if
                        if (size(output_values, 2) < size2) then

                            !> Expand the second dimension of the variable.
                            call allocate_field(dat2_temp, (/size(output_values, 1), size2/), ierr)
                            if (ierr == 0) then
                                dat2_temp(:, 1:size(output_values, 2)) = output_values
                                call move_alloc(dat2_temp, output_values)
                            end if
                        end if
                    end if
                end if
                if (ierr /= 0) then
                    call print_warning( &
                        "An allocation error occurred trying to assign the variable '" // trim(input_field%label) // ". " // &
                        "The variable cannot be assigned.")
                else if (.not. allocated(output_values)) then
                    call print_warning( &
                        "Unable to allocate the field to the desired dimensions '" // trim(output_dim_names(1)) // &
                        ", " // trim(output_dim_names(2)) // "'. " // &
                        "The variable '" // trim(input_field%label) // "' cannot be assigned.")
                    ierr = 1
                end if
                if (ierr /= 0) return

                !> Derive maps for the output variable.
                if (any(output_dim_names == DIM_NAME_X) .and. any(output_dim_names == DIM_NAME_Y)) then
                    call allocate_field(source_dim1_map, (/size1/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/size2/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        do i = 1, size1
                            source_dim1_map(i) = i
                        end do
                        do i = 1, size2
                            source_dim2_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        do i = 1, size1
                            source_dim1_map(i) = i
                        end do
                        source_dim2_map = size2
                    else
                        ierr = 3
                    end if
                else if (any(DIM_NAMES_OF_N == output_dim_names(1)) .and. associated(vs%grid)) then
                    call allocate_field(source_dim1_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%grid%from_grid_x
                        source_dim2_map = vs%grid%from_grid_y
                        if (mapped_orders(MAP_ORDER_K) > 0) then
                            call allocate_field(source_dim3_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                            source_dim3_map = vs%tile%from_riverclass
                        else
                            ierr = 3
                        end if
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        do i = 1, vs%cell_count
                            source_dim1_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_K) > 0 .or. any(DIM_NAMES_OF_K == input_field%level_dim_name)) then
                        source_dim1_map = vs%grid%from_riverclass
                    else
                        ierr = 3
                    end if
                else if (any(DIM_NAMES_OF_G == output_dim_names(1)) .and. associated(vs%tile)) then
                    call allocate_field(source_dim1_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%tile%from_grid_x
                        source_dim2_map = vs%tile%from_grid_y
                        if (mapped_orders(MAP_ORDER_M) > 0) then
                            call allocate_field(source_dim3_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                            source_dim3_map = vs%tile%from_gru
                        else
                            ierr = 3
                        end if
                    else if (mapped_orders(MAP_ORDER_G) > 0) then
                        do i = 1, vs%landtile_count
                            source_dim1_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_N) > 0 .and. mapped_orders(MAP_ORDER_M) > 0) then
                        source_dim1_map = vs%tile%from_cell
                        source_dim2_map = vs%tile%from_gru
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        source_dim1_map = vs%tile%from_cell
                    else if (mapped_orders(MAP_ORDER_M) > 0 .or. any(DIM_NAMES_OF_M == input_field%level_dim_name)) then
                        source_dim1_map = vs%tile%from_gru
                    else
                        ierr = 3
                    end if
                else
                    ierr = 3
                end if
                if (ierr /= 0 .and. ierr /= 3) then
                    call print_warning( &
                        "An error occurred allocating the mapping variables for the variable '" // trim(input_field%label) // "'.")
                    return
                end if

                !> Determine indices to assign the output variable.
                i_match = 0
                i_max = size(output_values, 1)
                j1 = 1
                j2 = size(output_values, 2)
                if (level_id > 0) then
                    if (len_trim(input_field%level_dim_name) == 0 .or. any(DIM_NAMES_OF_L == input_field%level_dim_name)) then

                        !> Check the level against the dimension.
                        if (level_id < j1 .or. level_id > j2) then
                            write(code, *) level_id
                            message = &
                                "The index [" // trim(adjustl(code)) // "] to assign for the variable '" // &
                                trim(input_field%label) // "' is outside of the range of active indices ["
                            write(code, *) j1
                            message = message // trim(adjustl(code)) // ":"
                            write(code, *) j2
                            message = &
                                message // trim(adjustl(code)) // "] in the second dimension. The variable cannot be assigned."
                            call print_warning(message)
                            return
                        else
                            j1 = level_id
                            j2 = level_id
                        end if
                    else if (any(source_dim1_map == NO_DATA_INT)) then

                        !> No mapping is active when a non-level mapping index is defined.
                        write(code, *) level_id
                        call print_warning( &
                            "A mapping index [" // trim(adjustl(code)) // "] is defined to assign variable '" // &
                            trim(input_field%label) // "' from the '" // trim(input_field%level_dim_name) // &
                            "' dimension but no mapping could be determined for the output variable to the '" // &
                            trim(output_dim_names(1)) // "' dimension. The variable cannot be assigned.")
                        return
                    else if (.not. any(source_dim1_map == level_id)) then

                        !> The map contains no indices that would map to the mapping index.
                        write(code, *) level_id
                        call print_warning( &
                            "The variable '" // trim(input_field%label) // "' contains a mapping index [" // &
                            trim(adjustl(code)) // "] that does not map to any active element of the output variable in the '" // &
                            trim(output_dim_names(1)) // "' dimension. No values from this variable will be assigned.")
                        return
                    else
                        i_match = level_id
                    end if
                end if

                !> Value mask for distributed fields.
                allocate(active_value_mask(i_max))
                active_value_mask = .true.
                if (present(value_mask)) then
                    if (sum(mapped_orders) == 2 .and. mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        if (size(value_mask) == size(active_value_mask)) active_value_mask = value_mask
                        if (any(.not. active_value_mask)) then
                            write(code, *) count(.not. active_value_mask)
                            call print_info( &
                                "A filter is active to exclude " // trim(adjustl(code)) // " values from being assigned " // &
                                "from the distributed input variable '" // trim(input_field%label) // "'. " // &
                                "The values already assigned from other non-distributed inputs will presist at these locations.")
                        end if
                    end if
                end if

                !> Assign the field.
                if (.not. allocated(input_field%dim_lengths)) then

                    !> Scalar.
                    call check_values(input_field, dat, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        do j = j1, j2
                            if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                where (source_dim1_map == i_match) output_values(:, j) = dat
                            else
                                output_values(:, j) = dat
                            end if
                        end do
                    end if
                else if (mappable_dim_count == 4) then

                    !> Check for compatible dimensions.
                    allocate(mapped_dim_order(4))
                    mapped_dim_order = 0
                    if ( &
                        all(source_dim1_map /= NO_DATA_INT) .and. all(source_dim2_map /= NO_DATA_INT) .and. &
                        all(source_dim3_map /= NO_DATA_INT) .and. &
                        mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0 .and. mapped_orders(MAP_ORDER_L) > 0) &
                        then
                        if (mapped_orders(MAP_ORDER_K) > 0) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_K, DIM_NAME_L/), mapped_dim_order, ierr)
                        else if (mapped_orders(MAP_ORDER_M) > 0) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_M, DIM_NAME_L/), mapped_dim_order, ierr)
                        end if
                    end if
                    if (any(mapped_dim_order == 0)) then

                        !> Can only be mapped from the specified dimension combinations.
                        ierr = 3
                    else

                        !> Assign the values to a temporary variable.
                        call check_values(input_field, dat4, valid_max, valid_min, ierr)
                        if (ierr == 0) then

                            !> Ensure the dimensions of the field are in the expected mapped order.
                            call map_dimensions(dat4, mapped_dim_order, dat4_mapped, ierr)

                            !> Check the dimensions and assign values.
                            if (j1 > size(dat4_mapped, 4)) then

                                !> The values to assign are outside the range of the level dimension.
                                write(code, *) j1
                                message = &
                                    "The first index to assign to the model variable [" // trim(adjustl(code)) // "] " // &
                                    "is higher than the number of values in the level dimension of the variable '" // &
                                    trim(input_field%label) // "' ("
                                write(code, *) size(dat4_mapped, 4)
                                message = message // trim(adjustl(code)) // "). The variable cannot be assigned."
                                call print_warning(message)
                                ierr = 1
                            else if (j2 > size(dat4_mapped, 4)) then

                                !> The source contains lesser values than expected.
                                write(code, *) size(dat4_mapped, 4)
                                message = &
                                    "The variable '" // trim(input_field%label) // "' contains " // trim(adjustl(code)) // &
                                    " values in the level dimension where "
                                write(code, *) j2
                                message = message // trim(adjustl(code)) // " are expected. The value at index ["
                                write(code, *) size(dat4_mapped, 4)
                                message = &
                                    message // trim(adjustl(code)) // &
                                    "] will be used to assign higher indices of the assigned variable."
                                call print_warning(message)

                                !> Assign up to the available values.
                                forall (i = 1:i_max, j = j1:size(dat4_mapped, 4), active_value_mask(i))
                                    output_values(i, j) = &
                                        dat4_mapped(source_dim1_map(i), source_dim2_map(i), source_dim3_map(i), j)
                                end forall

                                !> Assign the last value to higher indices in the dimension.
                                forall (i = 1:i_max, j = (size(dat4_mapped, 4) + 1):j2, active_value_mask(i))
                                    output_values(i, j) = &
                                        dat4_mapped( &
                                            source_dim1_map(i), source_dim2_map(i), source_dim3_map(i), size(dat4_mapped, 4))
                                end forall
                            else
                                if (j2 < size(dat4_mapped, 4)) then

                                    !> The source contains more values than are required.
                                    write(code, *) size(dat4_mapped, 4)
                                    message = &
                                        "The variable '" // trim(input_field%label) // "' contains " // &
                                        trim(adjustl(code)) // " values in the level dimension where "
                                    write(code, *) j2
                                    message = &
                                        message // trim(adjustl(code)) // " are expected. The values of the variable " // &
                                        "at higher indices will not be assigned and will be ignored."
                                    call print_warning(message)
                                end if

                                !> Assign the values.
                                forall (i = 1:i_max, j = j1:j2, active_value_mask(i))
                                    output_values(i, j) = &
                                        dat4_mapped(source_dim1_map(i), source_dim2_map(i), source_dim3_map(i), j)
                                end forall
                            end if
                        end if
                    end if
                else if (mappable_dim_count == 3) then

                    !> Check for compatible dimensions.
                    if (all(source_dim1_map /= NO_DATA_INT) .and. all(source_dim2_map /= NO_DATA_INT)) then
                        if (mapped_orders(MAP_ORDER_K) > 0 .and. all(source_dim3_map /= NO_DATA_INT)) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_K/), mapped_dim_order, ierr)
                        else if (mapped_orders(MAP_ORDER_M) > 0 .and. all(source_dim3_map /= NO_DATA_INT)) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_M/), mapped_dim_order, ierr)
                        else if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_L/), mapped_dim_order, ierr)
                        else if (mapped_orders(MAP_ORDER_N) > 0 .and. mapped_orders(MAP_ORDER_M) > 0) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_N, DIM_NAME_M, DIM_NAME_L/), mapped_dim_order, ierr)
                        else
                            allocate(mapped_dim_order(3))
                            mapped_dim_order = 0
                        end if
                        if (any(mapped_dim_order == 0)) then

                            !> Can only be mapped to 2-D by cartesian grid.
                            ierr = 3
                        else

                            !> Assign the values to a temporary variable.
                            call check_values(input_field, dat3, valid_max, valid_min, ierr)
                            if (ierr == 0) then

                                !> Ensure the dimensions of the field are in the expected mapped order.
                                call map_dimensions(dat3, mapped_dim_order, dat3_mapped, ierr)

                                !> Check the dimensions and assign values.
                                if (j1 > size(dat3_mapped, 3)) then

                                    !> The values to assign are outside the range of the level dimension.
                                    write(code, *) j1
                                    message = &
                                        "The first index to assign to the model variable [" // trim(adjustl(code)) // "] " // &
                                        "is higher than the number of values in the level dimension of the variable '" // &
                                        trim(input_field%label) // "' ("
                                    write(code, *) size(dat3_mapped, 3)
                                    message = message // trim(adjustl(code)) // "). The variable cannot be assigned."
                                    call print_warning(message)
                                    ierr = 1
                                else if (j2 > size(dat3_mapped, 3)) then

                                    !> The source contains lesser values than expected.
                                    write(code, *) size(dat3_mapped, 3)
                                    message = &
                                        "The variable '" // trim(input_field%label) // "' contains " // trim(adjustl(code)) // &
                                        " values in the level dimension where "
                                    write(code, *) j2
                                    message = message // trim(adjustl(code)) // " are expected. The value at index ["
                                    write(code, *) size(dat3_mapped, 3)
                                    message = &
                                        message // trim(adjustl(code)) // &
                                        "] will be used to assign higher indices of the assigned variable."
                                    call print_warning(message)

                                    !> Assign up to the available values.
                                    forall (i = 1:i_max, j = j1:size(dat3_mapped, 3), active_value_mask(i))
                                        output_values(i, j) = dat3_mapped(source_dim1_map(i), source_dim2_map(i), j)
                                    end forall

                                    !> Assign the last value to higher indices in the dimension.
                                    forall (i = 1:i_max, j = (size(dat3_mapped, 3) + 1):j2, active_value_mask(i))
                                        output_values(i, j) = &
                                            dat3_mapped(source_dim1_map(i), source_dim2_map(i), size(dat3_mapped, 3))
                                    end forall
                                else
                                    if (j2 < size(dat3_mapped, 3)) then

                                        !> The source contains more values than are required.
                                        write(code, *) size(dat3_mapped, 3)
                                        message = &
                                            "The variable '" // trim(input_field%label) // "' contains " // &
                                            trim(adjustl(code)) // " values in the level dimension where "
                                        write(code, *) j2
                                        message = &
                                            message // trim(adjustl(code)) // " are expected. The values of the variable " // &
                                            "at higher indices will not be assigned and will be ignored."
                                        call print_warning(message)
                                    end if

                                    !> Assign the values.
                                    forall (i = 1:i_max, j = j1:j2, active_value_mask(i))
                                        output_values(i, j) = dat3_mapped(source_dim1_map(i), source_dim2_map(i), j)
                                    end forall
                                end if
                            end if
                        end if
                    else
                        ierr = 3
                    end if
                else if (mappable_dim_count == 2) then

                    !> Assign the values to a temporary variable.
                    call check_values(input_field, dat2, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT)) then
                            if (all(source_dim2_map /= NO_DATA_INT)) then

                                !> Get the order of the dimensions.
                                if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                                    call get_dimension_order( &
                                        input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y/), mapped_dim_order, ierr)
                                else if (mapped_orders(MAP_ORDER_N) > 0 .and. mapped_orders(MAP_ORDER_M) > 0) then
                                    call get_dimension_order( &
                                        input_field%dim_names, (/DIM_NAME_N, DIM_NAME_M/), mapped_dim_order, ierr)
                                else
                                    allocate(mapped_dim_order(2))
                                    mapped_dim_order = 0
                                end if
                                if (any(mapped_dim_order == 0)) then
                                    ierr = 3
                                else

                                    !> Map the field.
                                    call map_dimensions(dat2, mapped_dim_order, dat2_mapped, ierr)
                                    if (ierr /= 0) then
                                        ierr = 3
                                    else

                                        !> Assign the values.
                                        forall (i = 1:i_max, j = j1:j2, active_value_mask(i))
                                            output_values(i, j) = dat2_mapped(source_dim1_map(i), source_dim2_map(i))
                                        end forall
                                    end if
                                end if
                            else if (mapped_orders(MAP_ORDER_L) > 0) then

                                !> Get the order of the dimensions.
                                z = 0
                                if (mapped_orders(MAP_ORDER_N) > 0) then
                                    call get_dimension_order(input_field%dim_names, (/DIM_NAME_N, DIM_NAME_L/), mapped_dim_order, z)
                                else if (mapped_orders(MAP_ORDER_K) > 0) then
                                    call get_dimension_order(input_field%dim_names, (/DIM_NAME_K, DIM_NAME_L/), mapped_dim_order, z)
                                else if (mapped_orders(MAP_ORDER_G) > 0) then
                                    call get_dimension_order(input_field%dim_names, (/DIM_NAME_G, DIM_NAME_L/), mapped_dim_order, z)
                                else if (mapped_orders(MAP_ORDER_M) > 0) then
                                    call get_dimension_order(input_field%dim_names, (/DIM_NAME_M, DIM_NAME_L/), mapped_dim_order, z)
                                else
                                    ierr = 3
                                end if
                                if (z /= 0) then
                                    ierr = 3
                                else

                                    !> Map the field.
                                    call map_dimensions(dat2, mapped_dim_order, dat2_mapped, ierr)
                                    if (ierr /= 0) then
                                        ierr = 3
                                    else

                                        !> Check the dimensions and assign values.
                                        if (j1 > size(dat2_mapped, 2)) then

                                            !> The values to assign are outside the range of the level dimension.
                                            write(code, *) j1
                                            message = &
                                                "The first index to assign to the model variable [" // trim(adjustl(code)) // &
                                                "] is larger than the size of the level dimension of the variable '" // &
                                                trim(input_field%label) // "' ("
                                            write(code, *) size(dat2_mapped, 2)
                                            message = message // trim(adjustl(code)) // "). The variable cannot be assigned."
                                            call print_warning(message)
                                            ierr = 1
                                        else if (j2 > size(dat2_mapped, 2)) then

                                            !> The source contains lesser values than expected.
                                            write(code, *) size(dat2_mapped, 2)
                                            message = &
                                                "The variable '" // trim(input_field%label) // "' contains " // &
                                                trim(adjustl(code)) // " values in the level dimension where "
                                            write(code, *) j2
                                            message = message // trim(adjustl(code)) // " are expected. The value at index ["
                                            write(code, *) size(dat2_mapped, 2)
                                            message = &
                                                message // trim(adjustl(code)) // &
                                                "] will be used to assign values at higher indices of the assigned variable."
                                            call print_warning(message)

                                            !> Assign up to the available values.
                                            do j = j1, size(dat2_mapped, 2)
                                                if (i_match > 0) then
                                                    forall (i = 1:i_max, source_dim1_map(i) == i_match)
                                                        output_values(i, j) = dat2_mapped(i_match, j)
                                                    end forall
                                                else
                                                    do i = 1, i_max
                                                        output_values(i, j) = dat2_mapped(source_dim1_map(i), j)
                                                    end do
                                                end if
                                            end do

                                            !> Assign the last value to higher indices in the dimension.
                                            do j = (size(dat2_mapped, 2) + 1), j2
                                                if (i_match > 0) then
                                                    forall (i = 1:i_max, source_dim1_map(i) == i_match)
                                                        output_values(i, j) = dat2_mapped(i_match, size(dat2_mapped, 2))
                                                    end forall
                                                else
                                                    do i = 1, i_max
                                                        output_values(i, j) = dat2_mapped(source_dim1_map(i), size(dat2_mapped, 2))
                                                    end do
                                                end if
                                            end do
                                        else
                                            if (j2 < size(dat2_mapped, 2)) then

                                                !> The source contains more values than are required.
                                                write(code, *) size(dat2_mapped, 2)
                                                message = &
                                                    "The variable '" // trim(input_field%label) // "' contains " // &
                                                    trim(adjustl(code)) // " values in the level dimension where "
                                                write(code, *) j2
                                                message = &
                                                    message // trim(adjustl(code)) // " are expected. " // &
                                                    "Values at higher indices in the source variable will be ignored."
                                                call print_warning(message)
                                            end if

                                            !> Assign the values.
                                            do j = j1, j2
                                                if (i_match > 0) then
                                                    forall (i = 1:i_max, source_dim1_map(i) == i_match)
                                                        output_values(i, j) = dat2_mapped(i_match, j)
                                                    end forall
                                                else
                                                    forall (i = 1:i_max) output_values(i, j) = dat2_mapped(source_dim1_map(i), j)
                                                end if
                                            end do
                                        end if
                                    end if
                                end if
                            else
                                ierr = 3
                            end if
                        end if
                    end if
                else if (mappable_dim_count == 1 .or. sum(input_field%dim_lengths) == 1) then

                    !> Assign variable to temporary array.
                    call check_values(input_field, dat1, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT) .and. mapped_orders(MAP_ORDER_L) == 0) then

                            !> Assign mapped values.
                            if (i_match > 0) then
                                forall (i = 1:i_max, j = j1:j2, source_dim1_map(i) == i_match) output_values(i, j) = dat1(i_match)
                            else
                                forall (i = 1:i_max, j = j1:j2) output_values(i, j) = dat1(i)
                            end if
                        else

                            !> Check the dimensions and assign level values.
                            if (j1 > size(dat1)) then

                                !> The values to assign are outside the range of the level dimension.
                                write(code, *) j1
                                message = &
                                    "The first index to assign to the model variable [" // trim(adjustl(code)) // "] " // &
                                    "is higher than the number of values in the level dimension of the variable '" // &
                                    trim(input_field%label) // "' ("
                                write(code, *) size(dat1)
                                message = message // trim(adjustl(code)) // "). The variable cannot be assigned."
                                call print_warning(message)
                                ierr = 1
                            else if (j2 > size(dat1)) then

                                !> The source contains lesser values than expected.
                                write(code, *) size(dat1)
                                message = &
                                    "The variable '" // trim(input_field%label) // "' contains " // trim(adjustl(code)) // &
                                    " values in the level dimension where "
                                write(code, *) j2
                                message = message // trim(adjustl(code)) // " are expected. The value at index ["
                                write(code, *) size(dat1)
                                message = &
                                    message // trim(adjustl(code)) // &
                                    "] will be used to assign higher indices of the assigned variable."
                                call print_warning(message)

                                !> Assign up to the available values.
                                do j = j1, size(dat1)
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i, j) = dat1(j)
                                    else
                                        output_values(:, j) = dat1(j)
                                    end if
                                end do

                                !> Assign the last value to higher indices in the dimension.
                                do j = (size(dat1) + 1), j2
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i, j) = dat1(size(dat1))
                                    else
                                        output_values(:, j) = dat1(size(dat1))
                                    end if
                                end do
                            else
                                if (j2 < size(dat1)) then

                                    !> The source contains more values than are required.
                                    write(code, *) size(dat1)
                                    message = &
                                        "The variable '" // trim(input_field%label) // "' contains " // &
                                        trim(adjustl(code)) // " values in the level dimension where "
                                    write(code, *) j2
                                    message = &
                                        message // trim(adjustl(code)) // " are expected. The values of the variable " // &
                                        "at higher indices will not be assigned and will be ignored."
                                    call print_warning(message)
                                end if

                                !> Assign the values.
                                do j = j1, j2
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i, j) = dat1(j)
                                    else
                                        output_values(:, j) = dat1(j)
                                    end if
                                end do
                            end if
                        end if
                    end if
                else
                    ierr = 3
                end if
            end associate
        end if

        !> Check for errors (case '1' is handled inline and ignored).
        select case (ierr)
            case (3)
                call print_warning("An error occurred mapping the variable '" // trim(input_field%label) // "'.")
            case (0)
                call print_info("Assigned the variable '" // trim(input_field%label) // "'.")
        end select

    end subroutine

    subroutine get_field_value_real1d( &
        input_field, output_values, output_dim_name, level_map, value_mask, &
        const_mul, const_add, valid_max, valid_min)

        !> Modules.
        use model_variables, only: vs
        use parse_utilities, only: to_lowercase

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        real, allocatable :: output_values(:)
        character(len = *), intent(in), optional :: output_dim_name
        character(len = *), intent(in), optional :: level_map(:)
        logical, intent(in), optional :: value_mask(:)
        real, intent(in), optional :: const_mul
        real, intent(in), optional :: const_add
        real, intent(in), optional :: valid_max
        real, intent(in), optional :: valid_min

        !> Local variables.
        real, allocatable :: dat3(:, :, :), dat3_mapped(:, :, :), dat2(:, :), dat2_mapped(:, :), dat1(:), dat1_temp(:)
        real dat
        integer, allocatable :: source_dim1_map(:), source_dim2_map(:), source_dim3_map(:), mapped_dim_order(:)
        integer mappable_dim_count, size1, level_id, j2, j1, j, i_match, i_max, i, z, ierr
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message
        logical, allocatable :: active_value_mask(:)
        logical :: active_mapping

        !> Check if the dimensions of the field are defined for mapping.
        active_mapping = &
            (allocated(input_field%dim_names) .or. len_trim(input_field%level_dim_name) > 0) .and. present(output_dim_name)
!-        if (allocated(input_field%dim_names) .and. present(output_dim_name)) then
!-            active_mapping = .not. any(input_field%dim_names == output_dim_name)
!-        end if
        ierr = 0
        if (.not. active_mapping) then

            !> Try to assign the field.
            call check_values(input_field, output_values, const_mul, const_add, valid_max, valid_min, ierr)
        else if (allocated(input_field%mapping%mapped_dim_order)) then

            !> Check if the field is associated with a specific level ID.
            level_id = 0
            ierr = 0
            if (len_trim(input_field%level) > 0) then
                if (input_field%level_id > 0) then
                    level_id = input_field%level_id
                else if (present(level_map)) then
                    do i = 1, size(level_map)
                        if (to_lowercase(level_map(i)) == to_lowercase(input_field%level)) then
                            level_id = i
                            exit
                        end if
                    end do
                    if (level_id == 0) ierr = 1
                end if
                if (ierr /= 0) then
                    call print_warning( &
                        "Unabled to determine an index for the level '" // trim(input_field%level) // &
                        "' to assign the variable '" // trim(input_field%label) // "'.")
                    return
                end if
            end if

            !> Map the variable.
            associate(mapped_orders => input_field%mapping%mapped_dim_order)

                !> Save the number of mappable dimensions.
                mappable_dim_count = count(mapped_orders > 0)
                if (mapped_orders(MAP_ORDER_T) > 0) mappable_dim_count = mappable_dim_count - 1

                !> Derive dimensions for the output variable.
                size1 = 0
                if (any(DIM_NAMES_OF_N == output_dim_name)) then
                    if (vs%cell_count > 0) then
                        size1 = vs%cell_count
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_N))
                    end if
                else if (any(DIM_NAMES_OF_G == output_dim_name)) then
                    if (vs%landtile_count > 0) then
                        size1 = vs%landtile_count
                    else if (mapped_orders(MAP_ORDER_G) > 0) then
                        size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_G))
                    end if
                end if

                !> Check the allocation status of the output variable.
                ierr = 0
                if (size1 > 0) then
                    if (.not. allocated(output_values)) then

                        !> Allocate the variable.
                        call allocate_field(output_values, (/size1/), ierr)
                    else
                        if (size(output_values) < size1) then

                            !> Expand the first dimension of the variable.
                            call allocate_field(dat1_temp, (/size1/), ierr)
                            if (ierr == 0) then
                                dat1_temp(1:size(output_values)) = output_values
                                call move_alloc(dat1_temp, output_values)
                            end if
                        end if
                    end if
                end if
                if (ierr /= 0) then
                    call print_warning( &
                        "An allocation error occurred trying to assign the variable '" // trim(input_field%label) // ". " // &
                        "The variable cannot be assigned.")
                else if (.not. allocated(output_values)) then
                    call print_warning( &
                        "Unable to allocate the field to the desired dimension '" // trim(output_dim_name) // "'. " // &
                        "The variable '" // trim(input_field%label) // "' cannot be assigned.")
                    ierr = 1
                end if
                if (ierr /= 0) return

                !> Derive maps for the output variable.
                if (any(DIM_NAMES_OF_N == output_dim_name) .and. associated(vs%grid)) then
                    call allocate_field(source_dim1_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%grid%from_grid_x
                        source_dim2_map = vs%grid%from_grid_y
                        if (mapped_orders(MAP_ORDER_K) > 0) then
                            call allocate_field(source_dim3_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                            source_dim3_map = vs%tile%from_riverclass
                        else
                            ierr = 3
                        end if
                    else if (mapped_orders(MAP_ORDER_X) > 0) then
                        source_dim1_map = vs%grid%from_grid_x
                    else if (mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%grid%from_grid_y
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        do i = 1, vs%cell_count
                            source_dim1_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_K) > 0 .or. any(DIM_NAMES_OF_K == input_field%level_dim_name)) then
                        source_dim1_map = vs%grid%from_riverclass
                    else
                        ierr = 3
                    end if
                else if (any(DIM_NAMES_OF_G == output_dim_name) .and. associated(vs%tile)) then
                    call allocate_field(source_dim1_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%tile%from_grid_x
                        source_dim2_map = vs%tile%from_grid_y
                        if (mapped_orders(MAP_ORDER_M) > 0) then
                            call allocate_field(source_dim3_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                            source_dim3_map = vs%tile%from_gru
                        else
                            ierr = 3
                        end if
                    else if (mapped_orders(MAP_ORDER_X) > 0) then
                        source_dim1_map = vs%tile%from_grid_x
                    else if (mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%tile%from_grid_y
                    else if (mapped_orders(MAP_ORDER_G) > 0) then
                        do i = 1, vs%landtile_count
                            source_dim1_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_N) > 0 .and. mapped_orders(MAP_ORDER_M) > 0) then
                        source_dim1_map = vs%tile%from_cell
                        source_dim2_map = vs%tile%from_gru
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        source_dim1_map = vs%tile%from_cell
                    else if (mapped_orders(MAP_ORDER_M) > 0 .or. any(DIM_NAMES_OF_M == input_field%level_dim_name)) then
                        source_dim1_map = vs%tile%from_gru
                    else
                        ierr = 3
                    end if
                else if (allocated(input_field%dim_names) .and. allocated(output_values)) then
                    if (input_field%dim_names(1) == output_dim_name) then
                        call allocate_field(source_dim1_map, (/size(output_values)/), z); if (z /= 0) ierr = z
                        do i = 1, size(output_values)
                            source_dim1_map(i) = i
                        end do
                    end if
                else
                    ierr = 3
                end if
                if (ierr /= 0 .and. ierr /= 3) then
                    call print_warning( &
                        "An error occurred allocating the mapping variables for the variable '" // trim(input_field%label) // "'.")
                    return
                end if

                !> Determine indices to assign the output variable.
                i_match = 0
                i_max = size(output_values)
                j1 = 1
                j2 = 1
                if (level_id > 0) then
                    if (len_trim(input_field%level_dim_name) == 0 .or. any(DIM_NAMES_OF_L == input_field%level_dim_name)) then

                        !> Check the level against the dimension.
                        if (level_id < j1 .or. level_id > j2) then
                            write(code, *) level_id
                            message = &
                                "The index [" // trim(adjustl(code)) // "] to assign for the variable '" // &
                                trim(input_field%label) // "' is outside of the range of active indices ["
                            write(code, *) j1
                            message = message // trim(adjustl(code)) // ":"
                            write(code, *) j2
                            message = &
                                message // trim(adjustl(code)) // "] in the second dimension. The variable cannot be assigned."
                            call print_warning(message)
                            return
                        else
                            j1 = level_id
                            j2 = level_id
                        end if
                    else if (any(source_dim1_map == NO_DATA_INT)) then

                        !> No mapping is active when a non-level mapping index is defined.
                        write(code, *) level_id
                        call print_warning( &
                            "A mapping index [" // trim(adjustl(code)) // "] is defined to assign variable '" // &
                            trim(input_field%label) // "' from the '" // trim(input_field%level_dim_name) // &
                            "' dimension but no mapping could be determined for the output variable to the '" // &
                            trim(output_dim_name) // "' dimension. The variable cannot be assigned.")
                        return
                    else if (.not. any(source_dim1_map == level_id)) then

                        !> The map contains no indices that would map to the mapping index.
                        write(code, *) level_id
                        call print_warning( &
                            "The variable '" // trim(input_field%label) // "' contains a mapping index [" // &
                            trim(adjustl(code)) // "] that does not map to any active element of the output variable in the '" // &
                            trim(output_dim_name) // "' dimension. No values from this variable will be assigned.")
                        return
                    else
                        i_match = level_id
                    end if
                end if

                !> Value mask for distributed fields.
                allocate(active_value_mask(i_max))
                active_value_mask = .true.
                if (present(value_mask)) then
                    if (sum(mapped_orders) == 2 .and. mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        if (size(value_mask) == size(active_value_mask)) active_value_mask = value_mask
                        if (any(.not. active_value_mask)) then
                            write(code, *) count(.not. active_value_mask)
                            call print_info( &
                                "A filter is active to exclude " // trim(adjustl(code)) // " values from being assigned " // &
                                "from the distributed input variable '" // trim(input_field%label) // "'. " // &
                                "The values already assigned from other non-distributed inputs will presist at these locations.")
                        end if
                    end if
                end if

                !> Assign the field.
                if (.not. allocated(input_field%dim_lengths)) then

                    !> Scalar.
                    call check_values(input_field, dat, const_mul, const_add, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                            where (source_dim1_map == i_match) output_values(:) = dat
                        else
                            output_values(:) = dat
                        end if
                    end if
                else if (mappable_dim_count == 3) then

                    !> Check for compatible dimensions.
                    allocate(mapped_dim_order(3))
                    mapped_dim_order = 0
                    if (allocated(source_dim1_map) .and. allocated(source_dim2_map) .and. allocated(source_dim3_map)) then
                        if ( &
                            all(source_dim1_map /= NO_DATA_INT) .and. all(source_dim2_map /= NO_DATA_INT) .and. &
                            all(source_dim3_map /= NO_DATA_INT) .and. &
                            mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) &
                            then
                            if (mapped_orders(MAP_ORDER_K) > 0) then
                                call get_dimension_order( &
                                    input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_K/), mapped_dim_order, ierr)
                            else if (mapped_orders(MAP_ORDER_M) > 0) then
                                call get_dimension_order( &
                                    input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_M/), mapped_dim_order, ierr)
                            end if
                        end if
                    end if
                    if (any(mapped_dim_order == 0)) then

                        !> Can only be mapped from the specified dimension combinations.
                        ierr = 3
                    else

                        !> Assign the values to a temporary variable.
                        call check_values(input_field, dat3, const_mul, const_add, valid_max, valid_min, ierr)
                        if (ierr == 0) then

                            !> Ensure the dimensions of the field are in the expected mapped order.
                            call map_dimensions(dat3, mapped_dim_order, dat3_mapped, ierr)

                            !> Assign the values.
                            forall (i = 1:i_max, active_value_mask(i))
                                output_values(i) = dat3_mapped(source_dim1_map(i), source_dim2_map(i), source_dim3_map(i))
                            end forall
                        end if
                    end if
                else if (mappable_dim_count == 2) then

                    !> Assign the values to a temporary variable.
                    call check_values(input_field, dat2, const_mul, const_add, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT)) then
                            if (all(source_dim2_map /= NO_DATA_INT)) then

                                !> Get the order of the dimensions.
                                if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                                    call get_dimension_order( &
                                        input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y/), mapped_dim_order, ierr)
                                else if (mapped_orders(MAP_ORDER_N) > 0 .and. mapped_orders(MAP_ORDER_M) > 0) then
                                    call get_dimension_order( &
                                        input_field%dim_names, (/DIM_NAME_N, DIM_NAME_M/), mapped_dim_order, ierr)
                                else
                                    allocate(mapped_dim_order(2))
                                    mapped_dim_order = 0
                                end if
                                if (any(mapped_dim_order == 0)) then
                                    ierr = 3
                                else

                                    !> Map the field.
                                    call map_dimensions(dat2, mapped_dim_order, dat2_mapped, ierr)
                                    if (ierr /= 0) then
                                        ierr = 3
                                    else

                                        !> Assign the values.
                                        forall (i = 1:i_max, active_value_mask(i))
                                            output_values(i) = dat2_mapped(source_dim1_map(i), source_dim2_map(i))
                                        end forall
                                    end if
                                end if
                            else
                                ierr = 3
                            end if
                        end if
                    end if
                else if (mappable_dim_count == 1 .or. sum(input_field%dim_lengths) == 1) then

                    !> Assign variable to temporary array.
                    call check_values(input_field, dat1, const_mul, const_add, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT)) then
                            if (i_match > 0) then

                                !> Assign the values.
                                forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(i_match)
                            else
                                do j = 1, size(dat1)
                                    forall (i = 1:i_max, source_dim1_map(i) == j) output_values(i) = dat1(j)
                                end do
                            end if
                        else

                            !> Check the dimensions and assign level values.
                            if (j1 > size(dat1)) then

                                !> The values to assign are outside the range of the level dimension.
                                write(code, *) j1
                                message = &
                                    "The first index to assign to the model variable [" // trim(adjustl(code)) // "] " // &
                                    "is higher than the number of values in the level dimension of the variable '" // &
                                    trim(input_field%label) // "' ("
                                write(code, *) size(dat1)
                                message = message // trim(adjustl(code)) // "). The variable cannot be assigned."
                                call print_warning(message)
                                ierr = 1
                            else if (j2 > size(dat1)) then

                                !> The source contains lesser values than expected.
                                write(code, *) size(dat1)
                                message = &
                                    "The variable '" // trim(input_field%label) // "' contains " // trim(adjustl(code)) // &
                                    " values in the level dimension where "
                                write(code, *) j2
                                message = message // trim(adjustl(code)) // " are expected. The value at index ["
                                write(code, *) size(dat1)
                                message = &
                                    message // trim(adjustl(code)) // &
                                    "] will be used to assign higher indices of the assigned variable."
                                call print_warning(message)

                                !> Assign up to the available values.
                                do j = j1, size(dat1)
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(j)
                                    else
                                        output_values(:) = dat1(j)
                                    end if
                                end do

                                !> Assign the last value to higher indices in the dimension.
                                do j = (size(dat1) + 1), j2
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(size(dat1))
                                    else
                                        output_values(:) = dat1(size(dat1))
                                    end if
                                end do
                            else
                                if (j2 < size(dat1)) then

                                    !> The source contains more values than are required.
                                    write(code, *) size(dat1)
                                    message = &
                                        "The variable '" // trim(input_field%label) // "' contains " // &
                                        trim(adjustl(code)) // " values in the level dimension where "
                                    write(code, *) j2
                                    message = &
                                        message // trim(adjustl(code)) // " are expected. The values of the variable " // &
                                        "at higher indices will not be assigned and will be ignored."
                                    call print_warning(message)
                                end if

                                !> Assign the values.
                                do j = j1, j2
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(j)
                                    else
                                        output_values(:) = dat1(j)
                                    end if
                                end do
                            end if
                        end if
                    end if
                else
                    ierr = 3
                end if
            end associate
        end if

        !> Check for errors (case '1' is handled inline and ignored).
        select case (ierr)
            case (3)
                call print_warning("An error occurred mapping the variable '" // trim(input_field%label) // "'.")
            case (0)
                call print_info("Assigned the variable '" // trim(input_field%label) // "'.")
        end select

    end subroutine

    subroutine get_field_value_int1d( &
        input_field, output_values, output_dim_name, level_map, value_mask, &
        valid_max, valid_min)

        !> Modules.
        use model_variables, only: vs
        use parse_utilities, only: to_lowercase

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        integer, allocatable :: output_values(:)
        character(len = *), intent(in), optional :: output_dim_name
        character(len = *), intent(in), optional :: level_map(:)
        logical, intent(in), optional :: value_mask(:)
        integer, intent(in), optional :: valid_max
        integer, intent(in), optional :: valid_min

        !> Local variables.
        integer, allocatable :: &
            dat3(:, :, :), dat3_mapped(:, :, :), dat2(:, :), dat2_mapped(:, :), dat1(:), dat1_temp(:), &
            source_dim1_map(:), source_dim2_map(:), source_dim3_map(:), mapped_dim_order(:)
        integer dat, mappable_dim_count, size1, level_id, j2, j1, j, i_match, i_max, i, z, ierr
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message
        logical, allocatable :: active_value_mask(:)
        logical :: active_mapping

        !> Check if the dimensions of the field are defined for mapping.
        active_mapping = &
            (allocated(input_field%dim_names) .or. len_trim(input_field%level_dim_name) > 0) .and. present(output_dim_name)
!-        if (allocated(input_field%dim_names) .and. present(output_dim_name)) then
!-            active_mapping = .not. any(input_field%dim_names == output_dim_name)
!-        end if
        ierr = 0
        if (.not. active_mapping) then

            !> Try to assign the field.
            call check_values(input_field, output_values, valid_max, valid_min, ierr)
        else if (allocated(input_field%mapping%mapped_dim_order)) then

            !> Check if the field is associated with a specific level ID.
            level_id = 0
            ierr = 0
            if (len_trim(input_field%level) > 0) then
                if (input_field%level_id > 0) then
                    level_id = input_field%level_id
                else if (present(level_map)) then
                    do i = 1, size(level_map)
                        if (to_lowercase(level_map(i)) == to_lowercase(input_field%level)) then
                            level_id = i
                            exit
                        end if
                    end do
                    if (level_id == 0) ierr = 1
                end if
                if (ierr /= 0) then
                    call print_warning( &
                        "Unabled to determine an index for the level '" // trim(input_field%level) // &
                        "' to assign the variable '" // trim(input_field%label) // "'.")
                    return
                end if
            end if

            !> Map the variable.
            associate(mapped_orders => input_field%mapping%mapped_dim_order)

                !> Save the number of mappable dimensions.
                mappable_dim_count = count(mapped_orders > 0)
                if (mapped_orders(MAP_ORDER_T) > 0) mappable_dim_count = mappable_dim_count - 1

                !> Derive dimensions for the output variable.
                size1 = 0
                if (any(DIM_NAMES_OF_N == output_dim_name)) then
                    if (vs%cell_count > 0) then
                        size1 = vs%cell_count
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_N))
                    end if
                else if (any(DIM_NAMES_OF_G == output_dim_name)) then
                    if (vs%landtile_count > 0) then
                        size1 = vs%landtile_count
                    else if (mapped_orders(MAP_ORDER_G) > 0) then
                        size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_G))
                    end if
                end if

                !> Check the allocation status of the output variable.
                ierr = 0
                if (size1 > 0) then
                    if (.not. allocated(output_values)) then

                        !> Allocate the variable.
                        call allocate_field(output_values, (/size1/), ierr)
                    else
                        if (size(output_values) < size1) then

                            !> Expand the first dimension of the variable.
                            call allocate_field(dat1_temp, (/size1/), ierr)
                            if (ierr == 0) then
                                dat1_temp(1:size(output_values)) = output_values
                                call move_alloc(dat1_temp, output_values)
                            end if
                        end if
                    end if
                end if
                if (ierr /= 0) then
                    call print_warning( &
                        "An allocation error occurred trying to assign the variable '" // trim(input_field%label) // ". " // &
                        "The variable cannot be assigned.")
                else if (.not. allocated(output_values)) then
                    call print_warning( &
                        "Unable to allocate the field to the desired dimension '" // trim(output_dim_name) // "'. " // &
                        "The variable '" // trim(input_field%label) // "' cannot be assigned.")
                    ierr = 1
                end if
                if (ierr /= 0) return

                !> Derive maps for the output variable.
                if (any(DIM_NAMES_OF_N == output_dim_name) .and. associated(vs%grid)) then
                    call allocate_field(source_dim1_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%grid%from_grid_x
                        source_dim2_map = vs%grid%from_grid_y
                        if (mapped_orders(MAP_ORDER_K) > 0) then
                            call allocate_field(source_dim3_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                            source_dim3_map = vs%tile%from_riverclass
                        else
                            ierr = 3
                        end if
                    else if (mapped_orders(MAP_ORDER_X) > 0) then
                        source_dim1_map = vs%grid%from_grid_x
                    else if (mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%grid%from_grid_y
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        do i = 1, vs%cell_count
                            source_dim1_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_K) > 0 .or. any(DIM_NAMES_OF_K == input_field%level_dim_name)) then
                        source_dim1_map = vs%grid%from_riverclass
                    else
                        ierr = 3
                    end if
                else if (any(DIM_NAMES_OF_G == output_dim_name) .and. associated(vs%tile)) then
                    call allocate_field(source_dim1_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%tile%from_grid_x
                        source_dim2_map = vs%tile%from_grid_y
                        if (mapped_orders(MAP_ORDER_M) > 0) then
                            call allocate_field(source_dim3_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                            source_dim3_map = vs%tile%from_gru
                        else
                            ierr = 3
                        end if
                    else if (mapped_orders(MAP_ORDER_X) > 0) then
                        source_dim1_map = vs%tile%from_grid_x
                    else if (mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%tile%from_grid_y
                    else if (mapped_orders(MAP_ORDER_G) > 0) then
                        do i = 1, vs%landtile_count
                            source_dim1_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_N) > 0 .and. mapped_orders(MAP_ORDER_M) > 0) then
                        source_dim1_map = vs%tile%from_cell
                        source_dim2_map = vs%tile%from_gru
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        source_dim1_map = vs%tile%from_cell
                    else if (mapped_orders(MAP_ORDER_M) > 0 .or. any(DIM_NAMES_OF_M == input_field%level_dim_name)) then
                        source_dim1_map = vs%tile%from_gru
                    else
                        ierr = 3
                    end if
                else if (allocated(input_field%dim_names) .and. allocated(output_values)) then
                    if (input_field%dim_names(1) == output_dim_name) then
                        call allocate_field(source_dim1_map, (/size(output_values)/), z); if (z /= 0) ierr = z
                        do i = 1, size(output_values)
                            source_dim1_map(i) = i
                        end do
                    end if
                else
                    ierr = 3
                end if
                if (ierr /= 0 .and. ierr /= 3) then
                    call print_warning( &
                        "An error occurred allocating the mapping variables for the variable '" // trim(input_field%label) // "'.")
                    return
                end if

                !> Determine indices to assign the output variable.
                i_match = 0
                i_max = size(output_values)
                j1 = 1
                j2 = 1
                if (level_id > 0) then
                    if (len_trim(input_field%level_dim_name) == 0 .or. any(DIM_NAMES_OF_L == input_field%level_dim_name)) then

                        !> Check the level against the dimension.
                        if (level_id < j1 .or. level_id > j2) then
                            write(code, *) level_id
                            message = &
                                "The index [" // trim(adjustl(code)) // "] to assign for the variable '" // &
                                trim(input_field%label) // "' is outside of the range of active indices ["
                            write(code, *) j1
                            message = message // trim(adjustl(code)) // ":"
                            write(code, *) j2
                            message = &
                                message // trim(adjustl(code)) // "] in the second dimension. The variable cannot be assigned."
                            call print_warning(message)
                            return
                        else
                            j1 = level_id
                            j2 = level_id
                        end if
                    else if (any(source_dim1_map == NO_DATA_INT)) then

                        !> No mapping is active when a non-level mapping index is defined.
                        write(code, *) level_id
                        call print_warning( &
                            "A mapping index [" // trim(adjustl(code)) // "] is defined to assign variable '" // &
                            trim(input_field%label) // "' from the '" // trim(input_field%level_dim_name) // &
                            "' dimension but no mapping could be determined for the output variable to the '" // &
                            trim(output_dim_name) // "' dimension. The variable cannot be assigned.")
                        return
                    else if (.not. any(source_dim1_map == level_id)) then

                        !> The map contains no indices that would map to the mapping index.
                        write(code, *) level_id
                        call print_warning( &
                            "The variable '" // trim(input_field%label) // "' contains a mapping index [" // &
                            trim(adjustl(code)) // "] that does not map to any active element of the output variable in the '" // &
                            trim(output_dim_name) // "' dimension. No values from this variable will be assigned.")
                        return
                    else
                        i_match = level_id
                    end if
                end if

                !> Value mask for distributed fields.
                allocate(active_value_mask(i_max))
                active_value_mask = .true.
                if (present(value_mask)) then
                    if (sum(mapped_orders) == 2 .and. mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        if (size(value_mask) == size(active_value_mask)) active_value_mask = value_mask
                        if (any(.not. active_value_mask)) then
                            write(code, *) count(.not. active_value_mask)
                            call print_info( &
                                "A filter is active to exclude " // trim(adjustl(code)) // " values from being assigned " // &
                                "from the distributed input variable '" // trim(input_field%label) // "'. " // &
                                "The values already assigned from other non-distributed inputs will presist at these locations.")
                        end if
                    end if
                end if

                !> Assign the field.
                if (.not. allocated(input_field%dim_lengths)) then

                    !> Scalar.
                    call check_values(input_field, dat, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                            where (source_dim1_map == i_match) output_values(:) = dat
                        else
                            output_values(:) = dat
                        end if
                    end if
                else if (mappable_dim_count == 3) then

                    !> Check for compatible dimensions.
                    allocate(mapped_dim_order(3))
                    mapped_dim_order = 0
                    if ( &
                        all(source_dim1_map /= NO_DATA_INT) .and. all(source_dim2_map /= NO_DATA_INT) .and. &
                        all(source_dim3_map /= NO_DATA_INT) .and. &
                        mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) &
                        then
                        if (mapped_orders(MAP_ORDER_K) > 0) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_K/), mapped_dim_order, ierr)
                        else if (mapped_orders(MAP_ORDER_M) > 0) then
                            call get_dimension_order( &
                                input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y, DIM_NAME_M/), mapped_dim_order, ierr)
                        end if
                    end if
                    if (any(mapped_dim_order == 0)) then

                        !> Can only be mapped from the specified dimension combinations.
                        ierr = 3
                    else

                        !> Assign the values to a temporary variable.
                        call check_values(input_field, dat3, valid_max, valid_min, ierr)
                        if (ierr == 0) then

                            !> Ensure the dimensions of the field are in the expected mapped order.
                            call map_dimensions(dat3, mapped_dim_order, dat3_mapped, ierr)

                            !> Assign the values.
                            forall (i = 1:i_max, active_value_mask(i))
                                output_values(i) = dat3_mapped(source_dim1_map(i), source_dim2_map(i), source_dim3_map(i))
                            end forall
                        end if
                    end if
                else if (mappable_dim_count == 2) then

                    !> Assign the values to a temporary variable.
                    call check_values(input_field, dat2, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT)) then
                            if (all(source_dim2_map /= NO_DATA_INT)) then

                                !> Get the order of the dimensions.
                                if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                                    call get_dimension_order( &
                                        input_field%dim_names, (/DIM_NAME_X, DIM_NAME_Y/), mapped_dim_order, ierr)
                                else if (mapped_orders(MAP_ORDER_N) > 0 .and. mapped_orders(MAP_ORDER_M) > 0) then
                                    call get_dimension_order( &
                                        input_field%dim_names, (/DIM_NAME_N, DIM_NAME_M/), mapped_dim_order, ierr)
                                else
                                    allocate(mapped_dim_order(2))
                                    mapped_dim_order = 0
                                end if
                                if (any(mapped_dim_order == 0)) then
                                    ierr = 3
                                else

                                    !> Map the field.
                                    call map_dimensions(dat2, mapped_dim_order, dat2_mapped, ierr)
                                    if (ierr /= 0) then
                                        ierr = 3
                                    else

                                        !> Assign the values.
                                        forall (i = 1:i_max, active_value_mask(i))
                                            output_values(i) = dat2_mapped(source_dim1_map(i), source_dim2_map(i))
                                        end forall
                                    end if
                                end if
                            else
                                ierr = 3
                            end if
                        end if
                    end if
                else if (mappable_dim_count == 1 .or. sum(input_field%dim_lengths) == 1) then

                    !> Assign variable to temporary array.
                    call check_values(input_field, dat1, valid_max, valid_min, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT)) then
                            if (i_match > 0) then

                                !> Assign the values.
                                forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(i_match)
                            else
                                do j = 1, size(dat1)
                                    forall (i = 1:i_max, source_dim1_map(i) == j) output_values(i) = dat1(j)
                                end do
                            end if
                        else

                            !> Check the dimensions and assign level values.
                            if (j1 > size(dat1)) then

                                !> The values to assign are outside the range of the level dimension.
                                write(code, *) j1
                                message = &
                                    "The first index to assign to the model variable [" // trim(adjustl(code)) // "] " // &
                                    "is higher than the number of values in the level dimension of the variable '" // &
                                    trim(input_field%label) // "' ("
                                write(code, *) size(dat1)
                                message = message // trim(adjustl(code)) // "). The variable cannot be assigned."
                                call print_warning(message)
                                ierr = 1
                            else if (j2 > size(dat1)) then

                                !> The source contains lesser values than expected.
                                write(code, *) size(dat1)
                                message = &
                                    "The variable '" // trim(input_field%label) // "' contains " // trim(adjustl(code)) // &
                                    " values in the level dimension where "
                                write(code, *) j2
                                message = message // trim(adjustl(code)) // " are expected. The value at index ["
                                write(code, *) size(dat1)
                                message = &
                                    message // trim(adjustl(code)) // &
                                    "] will be used to assign higher indices of the assigned variable."
                                call print_warning(message)

                                !> Assign up to the available values.
                                do j = j1, size(dat1)
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(j)
                                    else
                                        output_values(:) = dat1(j)
                                    end if
                                end do

                                !> Assign the last value to higher indices in the dimension.
                                do j = (size(dat1) + 1), j2
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(size(dat1))
                                    else
                                        output_values(:) = dat1(size(dat1))
                                    end if
                                end do
                            else
                                if (j2 < size(dat1)) then

                                    !> The source contains more values than are required.
                                    write(code, *) size(dat1)
                                    message = &
                                        "The variable '" // trim(input_field%label) // "' contains " // &
                                        trim(adjustl(code)) // " values in the level dimension where "
                                    write(code, *) j2
                                    message = &
                                        message // trim(adjustl(code)) // " are expected. The values of the variable " // &
                                        "at higher indices will not be assigned and will be ignored."
                                    call print_warning(message)
                                end if

                                !> Assign the values.
                                do j = j1, j2
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(j)
                                    else
                                        output_values(:) = dat1(j)
                                    end if
                                end do
                            end if
                        end if
                    end if
                else
                    ierr = 3
                end if
            end associate
        end if

        !> Check for errors (case '1' is handled inline and ignored).
        select case (ierr)
            case (3)
                call print_warning("An error occurred mapping the variable '" // trim(input_field%label) // "'.")
            case (0)
                call print_info("Assigned the variable '" // trim(input_field%label) // "'.")
        end select

    end subroutine

    subroutine get_field_value_char1d( &
        input_field, output_values, output_dim_name, level_map, value_mask)

        !> Modules.
        use model_variables, only: vs
        use parse_utilities, only: to_lowercase

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        character(len = *), allocatable :: output_values(:)
        character(len = *), intent(in), optional :: output_dim_name
        character(len = *), intent(in), optional :: level_map(:)
        logical, intent(in), optional :: value_mask(:)

        !> Local variables.
        integer, allocatable :: source_dim1_map(:), source_dim2_map(:), mapped_dim_order(:)
        integer mappable_dim_count, size1, level_id, j2, j1, j, i_match, i_max, i, z, ierr
        character(len = len(output_values)), allocatable :: dat1(:), dat1_temp(:)
        character(len = len(output_values)) dat
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: message
        logical, allocatable :: active_value_mask(:)
        logical :: active_mapping

        !> Check if the dimensions of the field are defined for mapping.
        active_mapping = &
            (allocated(input_field%dim_names) .or. len_trim(input_field%level_dim_name) > 0) .and. present(output_dim_name)
!-        if (allocated(input_field%dim_names) .and. present(output_dim_name)) then
!-            active_mapping = .not. any(input_field%dim_names == output_dim_name)
!-        end if
        ierr = 0
        if (.not. active_mapping) then

            !> Try to assign the field.
            call check_values(input_field, output_values, ierr)
        else if (allocated(input_field%mapping%mapped_dim_order)) then

            !> Check if the field is associated with a specific level ID.
            level_id = 0
            ierr = 0
            if (len_trim(input_field%level) > 0) then
                if (input_field%level_id > 0) then
                    level_id = input_field%level_id
                else if (present(level_map)) then
                    do i = 1, size(level_map)
                        if (to_lowercase(level_map(i)) == to_lowercase(input_field%level)) then
                            level_id = i
                            exit
                        end if
                    end do
                    if (level_id == 0) ierr = 1
                end if
                if (ierr /= 0) then
                    call print_warning( &
                        "Unabled to determine an index for the level '" // trim(input_field%level) // &
                        "' to assign the variable '" // trim(input_field%label) // "'.")
                    return
                end if
            end if

            !> Map the variable.
            associate(mapped_orders => input_field%mapping%mapped_dim_order)

                !> Save the number of mappable dimensions.
                mappable_dim_count = count(mapped_orders > 0)
                if (mapped_orders(MAP_ORDER_T) > 0) mappable_dim_count = mappable_dim_count - 1

                !> Derive dimensions for the output variable.
                size1 = 0
                if (any(DIM_NAMES_OF_N == output_dim_name)) then
                    if (vs%cell_count > 0) then
                        size1 = vs%cell_count
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_N))
                    end if
                else if (any(DIM_NAMES_OF_G == output_dim_name)) then
                    if (vs%landtile_count > 0) then
                        size1 = vs%landtile_count
                    else if (mapped_orders(MAP_ORDER_G) > 0) then
                        size1 = input_field%dim_lengths(mapped_orders(MAP_ORDER_G))
                    end if
                end if

                !> Check the allocation status of the output variable.
                ierr = 0
                if (size1 > 0) then
                    if (.not. allocated(output_values)) then

                        !> Allocate the variable.
                        call allocate_field(output_values, (/size1/), ierr)
                    else
                        if (size(output_values) < size1) then

                            !> Expand the first dimension of the variable.
                            call allocate_field(dat1_temp, (/size1/), ierr)
                            if (ierr == 0) then
                                dat1_temp(1:size(output_values)) = output_values
                                call move_alloc(dat1_temp, output_values)
                            end if
                        end if
                    end if
                end if
                if (ierr /= 0) then
                    call print_warning( &
                        "An allocation error occurred trying to assign the variable '" // trim(input_field%label) // ". " // &
                        "The variable cannot be assigned.")
                else if (.not. allocated(output_values)) then
                    call print_warning( &
                        "Unable to allocate the field to the desired dimension '" // trim(output_dim_name) // "'. " // &
                        "The variable '" // trim(input_field%label) // "' cannot be assigned.")
                    ierr = 1
                end if
                if (ierr /= 0) return

                !> Derive maps for the output variable.
                if (any(DIM_NAMES_OF_N == output_dim_name) .and. associated(vs%grid)) then
                    call allocate_field(source_dim1_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/vs%cell_count/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%grid%from_grid_x
                        source_dim2_map = vs%grid%from_grid_y
                    else if (mapped_orders(MAP_ORDER_X) > 0) then
                        source_dim1_map = vs%grid%from_grid_x
                    else if (mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%grid%from_grid_y
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        do i = 1, vs%cell_count
                            source_dim1_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_K) > 0 .or. any(DIM_NAMES_OF_K == input_field%level_dim_name)) then
                        source_dim1_map = vs%grid%from_riverclass
                    else
                        ierr = 3
                    end if
                else if (any(DIM_NAMES_OF_G == output_dim_name) .and. associated(vs%tile)) then
                    call allocate_field(source_dim1_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                    call allocate_field(source_dim2_map, (/vs%landtile_count/), z); if (z /= 0) ierr = z
                    if (mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%tile%from_grid_x
                        source_dim2_map = vs%tile%from_grid_y
                    else if (mapped_orders(MAP_ORDER_X) > 0) then
                        source_dim1_map = vs%tile%from_grid_x
                    else if (mapped_orders(MAP_ORDER_Y) > 0) then
                        source_dim1_map = vs%tile%from_grid_y
                    else if (mapped_orders(MAP_ORDER_G) > 0) then
                        do i = 1, vs%landtile_count
                            source_dim1_map(i) = i
                        end do
                    else if (mapped_orders(MAP_ORDER_N) > 0 .and. mapped_orders(MAP_ORDER_M) > 0) then
                        source_dim1_map = vs%tile%from_cell
                        source_dim2_map = vs%tile%from_gru
                    else if (mapped_orders(MAP_ORDER_N) > 0) then
                        source_dim1_map = vs%tile%from_cell
                    else if (mapped_orders(MAP_ORDER_M) > 0 .or. any(DIM_NAMES_OF_M == input_field%level_dim_name)) then
                        source_dim1_map = vs%tile%from_gru
                    else
                        ierr = 3
                    end if
                else if (allocated(input_field%dim_names) .and. allocated(output_values)) then
                    if (input_field%dim_names(1) == output_dim_name) then
                        call allocate_field(source_dim1_map, (/size(output_values)/), z); if (z /= 0) ierr = z
                        do i = 1, size(output_values)
                            source_dim1_map(i) = i
                        end do
                    end if
                else
                    ierr = 3
                end if
                if (ierr /= 0 .and. ierr /= 3) then
                    call print_warning( &
                        "An error occurred allocating the mapping variables for the variable '" // trim(input_field%label) // "'.")
                    return
                end if

                !> Determine indices to assign the output variable.
                i_match = 0
                i_max = size(output_values)
                j1 = 1
                j2 = 1
                if (level_id > 0) then
                    if (len_trim(input_field%level_dim_name) == 0 .or. any(DIM_NAMES_OF_L == input_field%level_dim_name)) then

                        !> Check the level against the dimension.
                        if (level_id < j1 .or. level_id > j2) then
                            write(code, *) level_id
                            message = &
                                "The index [" // trim(adjustl(code)) // "] to assign for the variable '" // &
                                trim(input_field%label) // "' is outside of the range of active indices ["
                            write(code, *) j1
                            message = message // trim(adjustl(code)) // ":"
                            write(code, *) j2
                            message = &
                                message // trim(adjustl(code)) // "] in the second dimension. The variable cannot be assigned."
                            call print_warning(message)
                            return
                        else
                            j1 = level_id
                            j2 = level_id
                        end if
                    else if (any(source_dim1_map == NO_DATA_INT)) then

                        !> No mapping is active when a non-level mapping index is defined.
                        write(code, *) level_id
                        call print_warning( &
                            "A mapping index [" // trim(adjustl(code)) // "] is defined to assign variable '" // &
                            trim(input_field%label) // "' from the '" // trim(input_field%level_dim_name) // &
                            "' dimension but no mapping could be determined for the output variable to the '" // &
                            trim(output_dim_name) // "' dimension. The variable cannot be assigned.")
                        return
                    else if (.not. any(source_dim1_map == level_id)) then

                        !> The map contains no indices that would map to the mapping index.
                        write(code, *) level_id
                        call print_warning( &
                            "The variable '" // trim(input_field%label) // "' contains a mapping index [" // &
                            trim(adjustl(code)) // "] that does not map to any active element of the output variable in the '" // &
                            trim(output_dim_name) // "' dimension. No values from this variable will be assigned.")
                        return
                    else
                        i_match = level_id
                    end if
                end if

                !> Value mask for distributed fields.
                allocate(active_value_mask(i_max))
                active_value_mask = .true.
                if (present(value_mask)) then
                    if (sum(mapped_orders) == 2 .and. mapped_orders(MAP_ORDER_X) > 0 .and. mapped_orders(MAP_ORDER_Y) > 0) then
                        if (size(value_mask) == size(active_value_mask)) active_value_mask = value_mask
                        if (any(.not. active_value_mask)) then
                            write(code, *) count(.not. active_value_mask)
                            call print_info( &
                                "A filter is active to exclude " // trim(adjustl(code)) // " values from being assigned " // &
                                "from the distributed input variable '" // trim(input_field%label) // "'. " // &
                                "The values already assigned from other non-distributed inputs will presist at these locations.")
                        end if
                    end if
                end if

                !> Assign the field.
                if (.not. allocated(input_field%dim_lengths)) then

                    !> Scalar.
                    call check_values(input_field, dat, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                            where (source_dim1_map == i_match) output_values(:) = dat
                        else
                            output_values(:) = dat
                        end if
                    end if
                else if (mappable_dim_count == 1 .or. sum(input_field%dim_lengths) == 1) then

                    !> Assign variable to temporary array.
                    call check_values(input_field, dat1, ierr)
                    if (ierr == 0) then
                        if (all(source_dim1_map /= NO_DATA_INT)) then
                            if (i_match > 0) then

                                !> Assign the values.
                                forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(i_match)
                            else
                                do j = 1, size(dat1)
                                    forall (i = 1:i_max, source_dim1_map(i) == j) output_values(i) = dat1(j)
                                end do
                            end if
                        else

                            !> Check the dimensions and assign level values.
                            if (j1 > size(dat1)) then

                                !> The values to assign are outside the range of the level dimension.
                                write(code, *) j1
                                message = &
                                    "The first index to assign to the model variable [" // trim(adjustl(code)) // "] " // &
                                    "is higher than the number of values in the level dimension of the variable '" // &
                                    trim(input_field%label) // "' ("
                                write(code, *) size(dat1)
                                message = message // trim(adjustl(code)) // "). The variable cannot be assigned."
                                call print_warning(message)
                                ierr = 1
                            else if (j2 > size(dat1)) then

                                !> The source contains lesser values than expected.
                                write(code, *) size(dat1)
                                message = &
                                    "The variable '" // trim(input_field%label) // "' contains " // trim(adjustl(code)) // &
                                    " values in the level dimension where "
                                write(code, *) j2
                                message = message // trim(adjustl(code)) // " are expected. The value at index ["
                                write(code, *) size(dat1)
                                message = &
                                    message // trim(adjustl(code)) // &
                                    "] will be used to assign higher indices of the assigned variable."
                                call print_warning(message)

                                !> Assign up to the available values.
                                do j = j1, size(dat1)
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(j)
                                    else
                                        output_values(:) = dat1(j)
                                    end if
                                end do

                                !> Assign the last value to higher indices in the dimension.
                                do j = (size(dat1) + 1), j2
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(size(dat1))
                                    else
                                        output_values(:) = dat1(size(dat1))
                                    end if
                                end do
                            else
                                if (j2 < size(dat1)) then

                                    !> The source contains more values than are required.
                                    write(code, *) size(dat1)
                                    message = &
                                        "The variable '" // trim(input_field%label) // "' contains " // &
                                        trim(adjustl(code)) // " values in the level dimension where "
                                    write(code, *) j2
                                    message = &
                                        message // trim(adjustl(code)) // " are expected. The values of the variable " // &
                                        "at higher indices will not be assigned and will be ignored."
                                    call print_warning(message)
                                end if

                                !> Assign the values.
                                do j = j1, j2
                                    if (all(source_dim1_map /= NO_DATA_INT) .and. i_match > 0) then
                                        forall (i = 1:i_max, source_dim1_map(i) == i_match) output_values(i) = dat1(j)
                                    else
                                        output_values(:) = dat1(j)
                                    end if
                                end do
                            end if
                        end if
                    end if
                else
                    ierr = 3
                end if
            end associate
        end if

        !> Check for errors (case '1' is handled inline and ignored).
        select case (ierr)
            case (3)
                call print_warning("An error occurred mapping the variable '" // trim(input_field%label) // "'.")
            case (0)
                call print_info("Assigned the variable '" // trim(input_field%label) // "'.")
        end select

    end subroutine

    subroutine get_field_value_scalar(input_field, output_value)

        !> Input/output variables.
        type(io_field), intent(in) :: input_field
        class(*) output_value

        !> Local variables.
        integer ierr

        !> Assign the field.
        call assign_field(input_field%field, output_value, ierr)
        if (ierr /= 0) then

            !> Print a warning if the function return a non-normal exit.
            call print_warning("An error occurred assigning the variable '" // trim(input_field%label) // "'.")
        else

            !> Print an info message that the variable was assigned.
            call print_info("Assigned the variable '" // trim(input_field%label) // "'.")
        end if

    end subroutine

    subroutine get_field_value_from_list_ids_real1d(field_list, field_ids, output_values)

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        integer, intent(in) :: field_ids(:)
        real, allocatable :: output_values(:)

        !> Local variables.
        integer i

        !> Assign the fields.
        do i = 1, size(field_ids)
            call get_field_value(field_list(field_ids(i)), output_values)
        end do

    end subroutine

    subroutine get_field_value_from_list_ids_int1d(field_list, field_ids, output_values)

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        integer, intent(in) :: field_ids(:)
        integer, allocatable :: output_values(:)

        !> Local variables.
        integer i

        !> Assign the fields.
        do i = 1, size(field_ids)
            call get_field_value(field_list(field_ids(i)), output_values)
        end do

    end subroutine

    subroutine get_field_value_from_list_ids_char1d(field_list, field_ids, output_values)

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        integer, intent(in) :: field_ids(:)
        character(len = *), allocatable :: output_values(:)

        !> Local variables.
        integer i

        !> Assign the fields.
        do i = 1, size(field_ids)
            call get_field_value(field_list(field_ids(i)), output_values)
        end do

    end subroutine

    subroutine get_field_value_from_list_ids_scalar(field_list, field_ids, output_value)

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        integer, intent(in) :: field_ids(:)
        class(*) output_value

        !> Local variables.
        integer i

        !> Assign the fields.
        do i = 1, size(field_ids)
            call get_field_value(field_list(field_ids(i)), output_value)
        end do

    end subroutine

    subroutine get_field_value_from_list_real1d(field_list, field_names, output_values)

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        character(len = *), intent(in) :: field_names(:)
        real, allocatable :: output_values(:)

        !> Local variables.
        integer i

        !> Check for variable in list.
        do i = 1, size(field_list)
            if (any(field_list(i)%field_name == field_names) .and. allocated(field_list(i)%field)) then

                !> Assign the field.
                call get_field_value(field_list(i), output_values)
                exit
            end if
        end do

    end subroutine

    subroutine get_field_value_from_list_int1d(field_list, field_names, output_values)

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        character(len = *), intent(in) :: field_names(:)
        integer, allocatable :: output_values(:)

        !> Local variables.
        integer i

        !> Check for variable in list.
        do i = 1, size(field_list)
            if (any(field_list(i)%field_name == field_names) .and. allocated(field_list(i)%field)) then

                !> Assign the field.
                call get_field_value(field_list(i), output_values)
                exit
            end if
        end do

    end subroutine

    subroutine get_field_value_from_list_char1d(field_list, field_names, output_values)

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        character(len = *), intent(in) :: field_names(:)
        character(len = *), allocatable :: output_values(:)

        !> Local variables.
        integer i

        !> Check for variable in list.
        do i = 1, size(field_list)
            if (any(field_list(i)%field_name == field_names) .and. allocated(field_list(i)%field)) then

                !> Assign the field.
                call get_field_value(field_list(i), output_values)
                exit
            end if
        end do

    end subroutine

    subroutine get_field_value_from_list_scalar(field_list, field_names, output_value)

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        character(len = *), intent(in) :: field_names(:)
        class(*) output_value

        !> Local variables.
        integer i, ierr

        !> Check for variable in list.
        do i = 1, size(field_list)
            if (any(field_list(i)%field_name == field_names) .and. allocated(field_list(i)%field)) then

                !> Assign the field.
                call get_field_value(field_list(i), output_value)
                exit
            end if
        end do

    end subroutine

    function is_value_in_field(field_list, field_names, target_values) result(ltest)

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        character(len = *), intent(in) :: field_names(:)
        class(*), intent(in) :: target_values(:)

        !> Local variables.
        real, allocatable :: dat1_r(:)
        integer, allocatable :: dat1_i(:)
        integer j2, j1, j, i, ierr
        character(len = SHORT_FIELD_LENGTH), allocatable :: dat1_c(:)
        logical ltest

        !> Initialize return variable.
        ltest = .false.

        !> Array bounds.
        j1 = lbound(target_values, 1)
        j2 = ubound(target_values, 1)

        !> Check for variable in list.
        do i = size(field_list), 1, -1
            if (any(field_list(i)%field_name == field_names) .and. allocated(field_list(i)%field)) then

                !> Check for the value in the field.
                select type (target_values)
                    type is (real)
                        call assign_field(field_list(i)%field, dat1_r, ierr)
                        do j = j1, j2
                            if (ierr == 0) ltest = any(dat1_r == target_values(j))
                            if (ltest) exit
                        end do
                    type is (integer)
                        call assign_field(field_list(i)%field, dat1_i, ierr)
                        do j = j1, j2
                            if (ierr == 0) ltest = any(dat1_i == target_values(j))
                            if (ltest) exit
                        end do
                    type is (character(len = *))
                        call assign_field(field_list(i)%field, dat1_c, ierr)
                        do j = j1, j2
                            if (ierr == 0) ltest = any(dat1_c == target_values(j))
                            if (ltest) exit
                        end do
                end select
                exit
            end if
        end do

    end function

    function check_for_field_names_in_list(field_list, field_names) result(field_ids)

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        type(field_list_search_names), intent(in) :: field_names(:)
        type(field_list_search_result), allocatable :: field_ids(:)

        !> Local variables.
        integer j, i

        !> Allocate the return list.
        if (allocated(field_ids)) deallocate(field_ids)
        i = size(field_names)
        allocate(field_ids(i))

        !> Check if the fields exist in the list.
        do i = 1, size(field_list)
            do j = 1, size(field_names)
                if (any(field_names(j)%names == field_list(i)%field_name)) then

                    !> The field matches.
                    if (.not. allocated(field_ids(j)%ids)) then
                        field_ids(j)%ids = [i]
                    else
                        field_ids(j)%ids = [field_ids(j)%ids, i]
                    end if

                    !> Save the field label (only the last if multiple exist).
                    field_ids(j)%label = field_list(i)%label
                end if
            end do
        end do

    end function

    subroutine get_cell_values_from_list_real1d(field_list, field_names, cell_values, mask)

        !> Modules.
        use model_variables, only: vs

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        character(len = *), intent(in) :: field_names(:)
        real, allocatable :: cell_values(:)
        logical, intent(in), optional :: mask(:)

        !> Local variables.
        type(model_variable_real1d) values
        integer i, ierr

        !> Return if grid values are not active.
        if (.not. associated(vs%grid)) return

        !> Check for variable in list.
        do i = 1, size(field_list)
            if (any(field_list(i)%label == field_names) .and. allocated(field_list(i)%field)) then

                !> Assign mapped values to a temporary array.
                if (allocated(field_list(i)%mapping%cell_map)) then
                    call assign_mapped_value(field_list(i), values, field_list(i)%mapping%cell_map, ierr)
                    if (ierr /= 0) then
                        call print_warning("An error occurred assigning the cell-variable for '" // trim(field_names(1)) // "'.")
                    end if
                end if
            end if
        end do

        !> Assign values.
        if (allocated(values%dat)) then

            !> Check if the variable is allocated.
            if (.not. allocated(cell_values)) then
                call allocate_field(cell_values, (/vs%grid%dim_length/), ierr)
                if (ierr /= 0) then
                    call print_warning("An error occurred allocating the cell-variable for '" // trim(field_names(1)) // "'.")
                end if
            end if

            !> Transfer from the field.
            if (present(mask)) then
                where (mask) cell_values = values%dat
            else
                cell_values = values%dat
            end if
        else
            call print_warning("No input field was found to assign '" // trim(field_names(1)) // "'.")
            if (allocated(cell_values)) deallocate(cell_values)
        end if

    end subroutine

    subroutine get_cell_values_from_list_int1d(field_list, field_names, cell_values, mask)

        !> Modules.
        use model_variables, only: vs

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        character(len = *), intent(in) :: field_names(:)
        integer, allocatable :: cell_values(:)
        logical, intent(in), optional :: mask(:)

        !> Local variables.
        type(model_variable_int1d) values
        integer i, ierr

        !> Return if grid values are not active.
        if (.not. associated(vs%grid)) return

        !> Check for variable in list.
        do i = 1, size(field_list)
            if (any(field_list(i)%label == field_names) .and. allocated(field_list(i)%field)) then

                !> Assign mapped values to a temporary array.
                if (allocated(field_list(i)%mapping%cell_map)) then
                    call assign_mapped_value(field_list(i), values, field_list(i)%mapping%cell_map, ierr)
                    if (ierr /= 0) then
                        call print_warning("An error occurred assigning the cell-variable for '" // trim(field_names(1)) // "'.")
                    end if
                end if
            end if
        end do

        !> Assign values.
        if (allocated(values%dat)) then

            !> Check if the variable is allocated.
            if (.not. allocated(cell_values)) then
                call allocate_field(cell_values, (/vs%grid%dim_length/), ierr)
                if (ierr /= 0) then
                    call print_warning("An error occurred allocating the cell-variable for '" // trim(field_names(1)) // "'.")
                end if
            end if

            !> Transfer from the field.
            if (present(mask)) then
                where (mask) cell_values = values%dat
            else
                cell_values = values%dat
            end if
        else
            call print_warning("No input field was found to assign '" // trim(field_names(1)) // "'.")
            if (allocated(cell_values)) deallocate(cell_values)
        end if

    end subroutine

    subroutine get_cell_values_from_list_char1d(field_list, field_names, cell_values, mask)

        !> Modules.
        use model_variables, only: vs

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        character(len = *), intent(in) :: field_names(:)
        character(len = *), allocatable :: cell_values(:)
        logical, intent(in), optional :: mask(:)

        !> Local variables.
        type(model_variable_char1d) values
        integer i, ierr

        !> Return if grid values are not active.
        if (.not. associated(vs%grid)) return

        !> Check for variable in list.
        do i = 1, size(field_list)
            if (any(field_list(i)%label == field_names) .and. allocated(field_list(i)%field)) then

                !> Assign mapped values to a temporary array.
                if (allocated(field_list(i)%mapping%cell_map)) then
                    call assign_mapped_value(field_list(i), values, field_list(i)%mapping%cell_map, ierr)
                    if (ierr /= 0) then
                        call print_warning("An error occurred assigning the cell-variable for '" // trim(field_names(1)) // "'.")
                    end if
                end if
            end if
        end do

        !> Assign values.
        if (allocated(values%dat)) then

            !> Check if the variable is allocated.
            if (.not. allocated(cell_values)) then
                call allocate_field(cell_values, (/vs%grid%dim_length/), ierr)
                if (ierr /= 0) then
                    call print_warning("An error occurred allocating the cell-variable for '" // trim(field_names(1)) // "'.")
                end if
            end if

            !> Transfer from the field.
            if (present(mask)) then
                where (mask) cell_values = values%dat
            else
                cell_values = values%dat
            end if
        else
            call print_warning("No input field was found to assign '" // trim(field_names(1)) // "'.")
            if (allocated(cell_values)) deallocate(cell_values)
        end if

    end subroutine

    subroutine get_tile_values_from_list_real1d(field_list, field_names, tile_values, mask)

        !> Modules.
        use model_variables, only: vs

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        character(len = *), intent(in) :: field_names(:)
        real, allocatable :: tile_values(:)
        logical, intent(in), optional :: mask(:)

        !> Local variables.
        type(model_variable_real1d) values
        integer i, ierr

        !> Return if tile values are not active.
        if (.not. associated(vs%tile)) return

        !> Check for variable in list.
        do i = 1, size(field_list)
            if (any(field_list(i)%label == field_names) .and. allocated(field_list(i)%field)) then

                !> Assign mapped values to a temporary array.
                if (allocated(field_list(i)%mapping%tile_map)) then
                    call assign_mapped_value(field_list(i), values, field_list(i)%mapping%tile_map, ierr)
                    if (ierr /= 0) then
                        call print_warning("An error occurred assigning the tile-variable for '" // trim(field_names(1)) // "'.")
                    end if
                end if
            end if
        end do

        !> Assign values.
        if (allocated(values%dat)) then

            !> Check if the variable is allocated.
            if (.not. allocated(tile_values)) then
                call allocate_field(tile_values, (/vs%tile%dim_length/), ierr)
                if (ierr /= 0) then
                    call print_warning("An error occurred allocating the tile-variable for '" // trim(field_names(1)) // "'.")
                end if
            end if

            !> Transfer from the field.
            if (present(mask)) then
                where (mask) tile_values = values%dat
            else
                tile_values = values%dat
            end if
        else
            call print_warning("No input field was found to assign '" // trim(field_names(1)) // "'.")
            if (allocated(tile_values)) deallocate(tile_values)
        end if

    end subroutine

    subroutine get_tile_values_from_list_int1d(field_list, field_names, tile_values, mask)

        !> Modules.
        use model_variables, only: vs

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        character(len = *), intent(in) :: field_names(:)
        integer, allocatable :: tile_values(:)
        logical, intent(in), optional :: mask(:)

        !> Local variables.
        type(model_variable_int1d) values
        integer i, ierr

        !> Return if tile values are not active.
        if (.not. associated(vs%tile)) return

        !> Check for variable in list.
        do i = 1, size(field_list)
            if (any(field_list(i)%label == field_names) .and. allocated(field_list(i)%field)) then

                !> Assign mapped values to a temporary array.
                if (allocated(field_list(i)%mapping%tile_map)) then
                    call assign_mapped_value(field_list(i), values, field_list(i)%mapping%tile_map, ierr)
                    if (ierr /= 0) then
                        call print_warning("An error occurred assigning the tile-variable for '" // trim(field_names(1)) // "'.")
                    end if
                end if
            end if
        end do

        !> Assign values.
        if (allocated(values%dat)) then

            !> Check if the variable is allocated.
            if (.not. allocated(tile_values)) then
                call allocate_field(tile_values, (/vs%tile%dim_length/), ierr)
                if (ierr /= 0) then
                    call print_warning("An error occurred allocating the tile-variable for '" // trim(field_names(1)) // "'.")
                end if
            end if

            !> Transfer from the field.
            if (present(mask)) then
                where (mask) tile_values = values%dat
            else
                tile_values = values%dat
            end if
        else
            call print_warning("No input field was found to assign '" // trim(field_names(1)) // "'.")
            if (allocated(tile_values)) deallocate(tile_values)
        end if

    end subroutine

    subroutine get_tile_values_from_list_char1d(field_list, field_names, tile_values, mask)

        !> Modules.
        use model_variables, only: vs

        !> Input/output variables.
        type(io_field), intent(in) :: field_list(:)
        character(len = *), intent(in) :: field_names(:)
        character(len = *), allocatable :: tile_values(:)
        logical, intent(in), optional :: mask(:)

        !> Local variables.
        type(model_variable_char1d) values
        integer i, ierr

        !> Return if tile values are not active.
        if (.not. associated(vs%tile)) return

        !> Check for variable in list.
        do i = 1, size(field_list)
            if (any(field_list(i)%label == field_names) .and. allocated(field_list(i)%field)) then

                !> Assign mapped values to a temporary array.
                if (allocated(field_list(i)%mapping%tile_map)) then
                    call assign_mapped_value(field_list(i), values, field_list(i)%mapping%tile_map, ierr)
                    if (ierr /= 0) then
                        call print_warning("An error occurred assigning the tile-variable for '" // trim(field_names(1)) // "'.")
                    end if
                end if
            end if
        end do

        !> Assign values.
        if (allocated(values%dat)) then

            !> Check if the variable is allocated.
            if (.not. allocated(tile_values)) then
                call allocate_field(tile_values, (/vs%tile%dim_length/), ierr)
                if (ierr /= 0) then
                    call print_warning("An error occurred allocating the tile-variable for '" // trim(field_names(1)) // "'.")
                end if
            end if

            !> Transfer from the field.
            if (present(mask)) then
                where (mask) tile_values = values%dat
            else
                tile_values = values%dat
            end if
        else
            call print_warning("No input field was found to assign '" // trim(field_names(1)) // "'.")
            if (allocated(tile_values)) deallocate(tile_values)
        end if

    end subroutine

    subroutine check_field_assigned_real3d_map_real3d(field_values, field_name, default_values, is_error, error_status)

        !> Input/output variables.
        real field_values(:, :, :)
        character(len = *), intent(in) :: field_name
        real, intent(in) :: default_values(:, :, :)
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (any(field_values == NO_DATA_REAL)) then
            if (size(field_values) == size(default_values)) then

                !> Assign default field.
                call print_remark("Assigning the unassigned values of the variable '" // trim(field_name) // "' default values.")
                where (field_values == NO_DATA_REAL) field_values = default_values
            else if (all(field_values == NO_DATA_REAL)) then
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            else
                error_message = "The variable '" // trim(field_name) // "' contains unassigned values."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_real2d_map_real2d(field_values, field_name, default_values, is_error, error_status)

        !> Input/output variables.
        real field_values(:, :)
        character(len = *), intent(in) :: field_name
        real, intent(in) :: default_values(:, :)
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (any(field_values == NO_DATA_REAL)) then
            if (size(field_values) == size(default_values)) then

                !> Assign default field.
                call print_remark("Assigning the unassigned values of the variable '" // trim(field_name) // "' default values.")
                where (field_values == NO_DATA_REAL) field_values = default_values
            else if (all(field_values == NO_DATA_REAL)) then
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            else
                error_message = "The variable '" // trim(field_name) // "' contains unassigned values."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_real2d(field_value, field_name, default_value, is_error, error_status)

        !> Input/output variables.
        real field_value(:, :)
        character(len = *), intent(in) :: field_name
        real, intent(in), optional :: default_value
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        integer i, n
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (any(field_value == NO_DATA_REAL)) then
            if (present(default_value)) then

                !> Assign default field.
                write(code, *) default_value
                call print_remark( &
                    "Assigning the unassigned values of the variable '" // trim(field_name) // "' a default value of '" // &
                    trim(adjustl(code)) // "'.")
                where (field_value == NO_DATA_REAL) field_value = default_value
            else if (all(field_value == NO_DATA_REAL)) then
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            else
                error_message = "The variable '" // trim(field_name) // "' contains unassigned values at index ["
                n = 1
                do i = 1, size(field_value)
                    if (any(field_value(i, :) == NO_DATA_REAL)) then
                        write(code, *) i
                        error_message = error_message // "[" // trim(adjustl(code)) // ", :]"
                        n = n + 1
                        if (n == 3) then
                            exit
                        else if (i < size(field_value)) then
                            error_message = error_message // ", "
                        end if
                    end if
                end do
                error_message = error_message // "]."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_int2d_map_int2d(field_values, field_name, default_values, is_error, error_status)

        !> Input/output variables.
        integer field_values(:, :)
        character(len = *), intent(in) :: field_name
        integer, intent(in) :: default_values(:, :)
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (any(field_values == NO_DATA_INT)) then
            if (size(field_values) == size(default_values)) then

                !> Assign default field.
                call print_remark("Assigning the unassigned values of the variable '" // trim(field_name) // "' default values.")
                where (field_values == NO_DATA_INT) field_values = default_values
            else if (all(field_values == NO_DATA_INT)) then
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            else
                error_message = "The variable '" // trim(field_name) // "' contains unassigned values."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_int2d(field_value, field_name, default_value, is_error, error_status)

        !> Input/output variables.
        integer field_value(:, :)
        character(len = *), intent(in) :: field_name
        integer, intent(in), optional :: default_value
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        integer i, n
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (any(field_value == NO_DATA_INT)) then
            if (present(default_value)) then

                !> Assign default field.
                write(code, *) default_value
                call print_remark( &
                    "Assigning the unassigned values of the variable '" // trim(field_name) // "' a default value of '" // &
                    trim(adjustl(code)) // "'.")
                where (field_value == NO_DATA_INT) field_value = default_value
            else if (all(field_value == NO_DATA_INT)) then
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            else
                error_message = "The variable '" // trim(field_name) // "' contains unassigned values at index ["
                n = 1
                do i = 1, size(field_value)
                    if (any(field_value(i, :) == NO_DATA_INT)) then
                        write(code, *) i
                        error_message = error_message // "[" // trim(adjustl(code)) // ", :]"
                        n = n + 1
                        if (n == 3) then
                            exit
                        else if (i < size(field_value)) then
                            error_message = error_message // ", "
                        end if
                    end if
                end do
                error_message = error_message // "]."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_real1d_map_real1d(field_values, field_name, default_values, is_error, error_status)

        !> Input/output variables.
        real field_values(:)
        character(len = *), intent(in) :: field_name
        real, intent(in) :: default_values(:)
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        integer i, n
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (any(field_values == NO_DATA_REAL)) then
            if (size(field_values) == size(default_values)) then

                !> Assign default field.
                call print_remark("Assigning the unassigned values of the variable '" // trim(field_name) // "' default values.")
                where (field_values == NO_DATA_REAL) field_values = default_values
            else if (all(field_values == NO_DATA_REAL)) then
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            else
                error_message = "The variable '" // trim(field_name) // "' contains unassigned values at index ["
                n = 1
                do i = 1, size(field_values)
                    if (field_values(i) == NO_DATA_REAL) then
                        write(code, *) i
                        error_message = error_message // trim(adjustl(code))
                        n = n + 1
                        if (n == 3) then
                            exit
                        else if (i < size(field_values)) then
                            error_message = error_message // ", "
                        end if
                    end if
                end do
                error_message = error_message // "]."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_real1d(field_value, field_name, default_value, is_error, error_status)

        !> Input/output variables.
        real field_value(:)
        character(len = *), intent(in) :: field_name
        real, intent(in), optional :: default_value
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        integer i, n
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (any(field_value == NO_DATA_REAL)) then
            if (present(default_value)) then

                !> Assign default field.
                write(code, *) default_value
                call print_remark( &
                    "Assigning the unassigned values of the variable '" // trim(field_name) // "' a default value of '" // &
                    trim(adjustl(code)) // "'.")
                where (field_value == NO_DATA_REAL) field_value = default_value
            else if (all(field_value == NO_DATA_REAL)) then
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            else
                error_message = "The variable '" // trim(field_name) // "' contains unassigned values at index ["
                n = 1
                do i = 1, size(field_value)
                    if (field_value(i) == NO_DATA_REAL) then
                        write(code, *) i
                        error_message = error_message // trim(adjustl(code))
                        n = n + 1
                        if (n == 3) then
                            exit
                        else if (i < size(field_value)) then
                            error_message = error_message // ", "
                        end if
                    end if
                end do
                error_message = error_message // "]."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_int1d_map_int1d(field_values, field_name, default_values, is_error, error_status)

        !> Input/output variables.
        integer field_values(:)
        character(len = *), intent(in) :: field_name
        integer, intent(in) :: default_values(:)
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        integer i, n
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (any(field_values == NO_DATA_INT)) then
            if (size(field_values) == size(default_values)) then

                !> Assign default field.
                call print_remark("Assigning the unassigned values of the variable '" // trim(field_name) // "' default values.")
                where (field_values == NO_DATA_INT) field_values = default_values
            else if (all(field_values == NO_DATA_INT)) then
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            else
                error_message = "The variable '" // trim(field_name) // "' contains unassigned values at index ["
                n = 1
                do i = 1, size(field_values)
                    if (field_values(i) == NO_DATA_INT) then
                        write(code, *) i
                        error_message = error_message // trim(adjustl(code))
                        n = n + 1
                        if (n == 3) then
                            exit
                        else if (i < size(field_values)) then
                            error_message = error_message // ", "
                        end if
                    end if
                end do
                error_message = error_message // "]."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_int1d(field_value, field_name, default_value, is_error, error_status)

        !> Input/output variables.
        integer field_value(:)
        character(len = *), intent(in) :: field_name
        integer, intent(in), optional :: default_value
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        integer i, n
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (any(field_value == NO_DATA_INT)) then
            if (present(default_value)) then

                !> Assign default field.
                write(code, *) default_value
                call print_remark( &
                    "Assigning the unassigned values of the variable '" // trim(field_name) // "' a default value of '" // &
                    trim(adjustl(code)) // "'.")
                where (field_value == NO_DATA_INT) field_value = default_value
            else if (all(field_value == NO_DATA_INT)) then
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            else
                error_message = "The variable '" // trim(field_name) // "' contains unassigned values at index ["
                n = 1
                do i = 1, size(field_value)
                    if (field_value(i) == NO_DATA_INT) then
                        write(code, *) i
                        error_message = error_message // trim(adjustl(code))
                        n = n + 1
                        if (n == 3) then
                            exit
                        else if (i < size(field_value)) then
                            error_message = error_message // ", "
                        end if
                    end if
                end do
                error_message = error_message // "]."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_char1d_map_char1d(field_values, field_name, default_values, is_error, error_status)

        !> Input/output variables.
        character(len = *) field_values(:)
        character(len = *), intent(in) :: field_name
        character(len = *), intent(in) :: default_values(:)
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        integer i, n
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (any(field_values == NO_DATA_CHAR)) then
            if (size(field_values) == size(default_values)) then

                !> Assign default field.
                call print_remark("Assigning the unassigned values of the variable '" // trim(field_name) // "' default values.")
                where (field_values == NO_DATA_CHAR) field_values = default_values
            else if (all(field_values == NO_DATA_CHAR)) then
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            else
                error_message = "The variable '" // trim(field_name) // "' contains unassigned values at index ["
                n = 1
                do i = 1, size(field_values)
                    if (field_values(i) == NO_DATA_CHAR) then
                        write(code, *) i
                        error_message = error_message // trim(adjustl(code))
                        n = n + 1
                        if (n == 3) then
                            exit
                        else if (i < size(field_values)) then
                            error_message = error_message // ", "
                        end if
                    end if
                end do
                error_message = error_message // "]."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_char1d(field_value, field_name, default_value, is_error, error_status)

        !> Input/output variables.
        character(len = *) field_value(:)
        character(len = *), intent(in) :: field_name
        character(len = *), intent(in), optional :: default_value
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        integer i, n
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (any(field_value == NO_DATA_CHAR)) then
            if (present(default_value)) then

                !> Assign default field.
                write(code, *) default_value
                call print_remark( &
                    "Assigning the unassigned values of the variable '" // trim(field_name) // "' a default value of '" // &
                    trim(adjustl(code)) // "'.")
                where (field_value == NO_DATA_CHAR) field_value = default_value
            else if (all(field_value == NO_DATA_CHAR)) then
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            else
                error_message = "The variable '" // trim(field_name) // "' contains unassigned values at index ["
                n = 1
                do i = 1, size(field_value)
                    if (field_value(i) == NO_DATA_CHAR) then
                        write(code, *) i
                        error_message = error_message // trim(adjustl(code))
                        n = n + 1
                        if (n == 3) then
                            exit
                        else if (i < size(field_value)) then
                            error_message = error_message // ", "
                        end if
                    end if
                end do
                error_message = error_message // "]."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_real(field_value, field_name, default_value, is_error, error_status)

        !> Input/output variables.
        real field_value
        character(len = *), intent(in) :: field_name
        real, intent(in), optional :: default_value
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (field_value == NO_DATA_REAL) then
            if (present(default_value)) then

                !> Assign default field.
                write(code, *) default_value
                call print_remark( &
                    "Assigning the unassigned value of variable '" // trim(field_name) // "' a default value of '" // &
                    trim(adjustl(code)) // "'.")
                field_value = default_value
            else
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_int(field_value, field_name, default_value, is_error, error_status)

        !> Input/output variables.
        integer field_value
        character(len = *), intent(in) :: field_name
        integer, intent(in), optional :: default_value
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (field_value == NO_DATA_INT) then
            if (present(default_value)) then

                !> Assign default field.
                write(code, *) default_value
                call print_remark( &
                    "Assigning the unassigned value of variable '" // trim(field_name) // "' a default value of '" // &
                    trim(adjustl(code)) // "'.")
                field_value = default_value
            else
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine check_field_assigned_char(field_value, field_name, default_value, is_error, error_status)

        !> Input/output variables.
        character(len = *) field_value
        character(len = *), intent(in) :: field_name
        character(len = *), intent(in), optional :: default_value
        logical, intent(in), optional :: is_error
        integer, intent(out) :: error_status

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code
        character(len = :), allocatable :: error_message
        logical is_warning

        !> Reset the return status.
        error_status = 0

        !> Prepare the message.
        is_warning = .true.
        if (present(is_error)) is_warning = .not. is_error
        error_message = ''

        !> Check if the field is assigned.
        if (field_value == NO_DATA_CHAR) then
            if (present(default_value)) then

                !> Assign default field.
                write(code, *) default_value
                call print_remark( &
                    "Assigning the unassigned value of variable '" // trim(field_name) // "' a default value of '" // &
                    trim(adjustl(code)) // "'.")
                field_value = default_value
            else
                error_message = "The variable '" // trim(field_name) // "' is not assigned."
            end if
        end if

        !> Check for errors.
        if (len_trim(error_message) > 0) then
            if (is_warning) then
                call print_warning(error_message)
            else
                call print_error(error_message)
            end if
            error_status = 1
        end if

    end subroutine

    subroutine parse_basin_info_flag(quiet, error_status)

        !> Modules.
        use parse_utilities, only: parse_line_values
        use strings, only: lowercase

        !> Input/output variables.
        logical, intent(in), optional :: quiet
        integer, intent(out) :: error_status

        !> Local variables.
        real, allocatable :: dat1_r(:)
        integer, allocatable :: dat1_i(:)
        integer j, i, ierr
        character(len = :), allocatable :: flag_name, file_label
        character(len = SHORT_FIELD_LENGTH), allocatable :: dim_names(:), values(:)
        logical v

        !> Check verbosity (override if 'DIAGNOSEMODE' is enabled).
        call reset_tab()
        v = .true.
        if (present(quiet)) v = .not. quiet
        if (DIAGNOSEMODE) v = .true.

        !> Find the flag in the list of input fields.
        flag_name = 'SHDFILEFLAG'
        call get_field_value(input_field_list, (/flag_name/), values)

        !> Check if the flag was found.
        if (.not. allocated(values)) then
            call print_remark("'" // flag_name // "' was not found. Assuming default file formats.")
            allocate(values(1))
            values = 'r2c'
        end if

        !> Identify the file type and assign the file label.
        file_label = ''
        do j = 1, size(values)

            !> Determine the file format by the old-style numeric format or from recognized word options.
            select case (lowercase(values(j)))
                case ('1', 'r2c')
                    file_label = 'r2c'
                case ('2', 'class_ini', 'point_mode')
                    file_label = 'point_mode'
                case ('3', 'netcdf', 'nc', '5', 'nc_subbasin', 'nc_hru')
                    file_label = 'nc'
                case ('4', 'asc')
                    file_label = 'asc'
            end select
        end do

        !> Check for errors.
        if (len_trim(file_label) == 0) then
            call print_error("'" // flag_name // "' specifies an unknown or unsupported file format.")
            error_status = 1
            return
        else

            !> Check if the file already exists or create an entry for the file.
            i = check_for_file_in_list(basin_info_files, 'basin_info_' // file_label, 'MESH_drainage_database.' // file_label)
            if (i == 0) then
                call print_error("An error occurred updating the file list for '" // flag_name // "'.")
                return
            end if

            !> Define the container for new files.
            if (.not. allocated(basin_info_files(i)%container)) then
                select case (file_label)
                    case ('r2c')

                        !> ASCII R2C format.
                        allocate(basin_info_files(i)%container, source = io_type_r2c())
                    case ('point_mode')

                        !> Point mode (location read from parameter files).
!-                        call reset_tab()
!-                        call print_message("READING: (creating basin for point mode)")
                        basin_info_files(i)%full_path = '(creating virtual file for point mode)'
!-                        call increase_tab()
                        allocate(dim_names(1), dat1_r(2), dat1_i(2))
                        dim_names = (/DIM_NAME_CELL/)
                        allocate(basin_info_files(i)%container, source = io_type_virtual())
!-                        allocate(basin_info_files(i)%fields(10))
                        allocate(basin_info_files(i)%fields(8))
                        basin_info_files(i)%fields(1)%label = 'projection'
                        allocate(basin_info_files(i)%fields(1)%field, source = model_variable_char(dat = 'LATLONG'))
                        basin_info_files(i)%fields(2)%label = 'ellipsoid'
                        allocate(basin_info_files(i)%fields(2)%field, source = model_variable_char(dat = 'GRS80'))
                        basin_info_files(i)%fields(3)%label = 'ngru'
                        allocate(basin_info_files(i)%fields(3)%field, source = model_variable_int(dat = 2))
                        dat1_i = (/1, 2/)
                        basin_info_files(i)%fields(4)%label = 'Rank'
                        allocate(basin_info_files(i)%fields(4)%field, source = model_variable_int1d(dat = dat1_i))
                        dat1_i = (/2, 0/)
                        basin_info_files(i)%fields(5)%label = 'Next'
                        allocate(basin_info_files(i)%fields(5)%field, source = model_variable_int1d(dat = dat1_i))
                        dat1_r = (/1.0, 0.0/)
                        basin_info_files(i)%fields(6)%label = 'GridArea'
                        allocate(basin_info_files(i)%fields(6)%field, source = model_variable_real1d(dat = dat1_r))
                        basin_info_files(i)%fields(7)%label = 'gru 1'
                        allocate(basin_info_files(i)%fields(7)%field, source = model_variable_real1d(dat = dat1_r))
                        dat1_r = (/0.0, 0.0/)
                        basin_info_files(i)%fields(8)%label = 'gru 2'
                        allocate(basin_info_files(i)%fields(8)%field, source = model_variable_real1d(dat = dat1_r))
!-                        basin_info_files(i)%fields(9)%label = 'lat'
!-                        allocate(basin_info_files(i)%fields(9)%field, source = model_variable_real1d(dat = dat1_r))
!-                        basin_info_files(i)%fields(10)%label = 'lon'
!-                        allocate(basin_info_files(i)%fields(10)%field, source = model_variable_real1d(dat = dat1_r))
!-                        do j = 4, 10
                        do j = 4, 8
                            allocate(basin_info_files(i)%fields(j)%dim_names(size(dim_names)), source = dim_names)
                        end do
                        deallocate(dim_names)
                    case ('nc')

                        !> NetCDF format.
#ifdef NETCDF
                        allocate(basin_info_files(i)%container, source = io_type_nc())
#else
                        call print_error( &
                            "The format of the input file is specified as NetCDF but the module is not active. " // &
                            "A version of MESH compiled with the NetCDF library must be used to read files in this format.")
                        error_status = 1
                        return
#endif
                    case ('asc')

                        !> Rank-ordered simple text (ASCII) format.
                        allocate(dim_names(1))
                        dim_names = (/DIM_NAME_N/)
                        allocate(basin_info_files(i)%overrides%dim_names(size(dim_names)), source = dim_names)
                        allocate(basin_info_files(i)%container, source = io_type_char_delimited())
                end select
            end if
        end if

        !> Parse the other options to update the file definition.
        do j = 1, size(values)
            if (lowercase(values(j)(1:6)) == 'fpath=') then

                !> Full path including file name and extension.
                basin_info_files(i)%full_path = trim(values(j)(7:))
            else if (lowercase(values(j)) == 'to_map') then

                !> Flag to create the diagnostic WATFLOOD 'map' file (from compatible formats).
                call expand_field_list(input_field_list, error_status = error_status)
                associate (input_field => input_field_list(size(input_field_list)))
                    input_field%label = 'SHDTOMAPFLAG'
                    allocate(input_field%field, source = model_variable_char(dat = 'on'))
                end associate
            end if
        end do

    end subroutine

    subroutine parse_parameter_input_flag(quiet, error_status)

        !> Modules.
        use parse_utilities, only: parse_line_values
        use strings, only: lowercase

        !> Input/output variables.
        logical, intent(in), optional :: quiet
        integer, intent(out) :: error_status

        !> Local variables.
        integer active_files, soil_ini_flag, j, i, ierr
        character(len = :), allocatable :: flag_name
        character(len = SHORT_FIELD_LENGTH), allocatable :: values(:)
        logical no_lss, no_route, v

        !> Check verbosity (override if 'DIAGNOSEMODE' is enabled).
        call reset_tab()
        v = .true.
        if (present(quiet)) v = .not. quiet
        if (DIAGNOSEMODE) v = .true.

        !> Find the flag in the list of input fields and other dependencies.
        flag_name = 'INPUTPARAMSFORMFLAG'
        call get_field_value(input_field_list, (/flag_name/), values)
        no_lss = is_value_in_field(input_field_list, (/'RUNMODE'/), (/'nolss ', 'no_lss'/))
        no_route = is_value_in_field(input_field_list, (/'RUNMODE'/), (/'noroute ', 'no_route'/))
        soil_ini_flag = 0
        call get_field_value(input_field_list, (/'SOILINIFLAG'/), soil_ini_flag)

        !> Check if the flag was found.
        if (.not. allocated(values)) then
            call print_remark("'" // flag_name // "' was not found. Assuming default file formats.")
            allocate(values(1))
            values = 'ini'
        end if

        !> Identify the active file types.
        active_files = 2**0
        do j = 1, size(values)

            !> Determine the file format from recognized word options.
            select case(lowercase(values(j)))
                case ('only')
                    active_files = 0
                case ('none')
                    active_files = 0
                    exit
                case ('r2c')
                    active_files = active_files + 2**1
                case ('csv')
                    active_files = active_files + 2**2
                case ('txt')
                    active_files = active_files + 2**3
                case ('netcdf', 'nc', 'nc_subbasin', 'nc_hru')
                    active_files = active_files + 2**4
            end select
        end do

        !> Check for errors (unless all modules are disabled).
        if (active_files == 0) then
            if (no_lss .and. no_route) then
                call print_remark("No parameter files are active.")
                return
            else
                error_status = 1
                call print_error("'" // flag_name // "' specifies an unknown or unsupported file format.")
                return
            end if
        else

            !> Check if the files already exist or create an entries for the files.
            ierr = 0
            if (btest(active_files, 0)) then

                !> Legacy 'INI' formats.
                if (.not. no_lss) then
                    i = check_for_file_in_list(parameter_files, 'class_ini', 'MESH_parameters_CLASS.ini')
                    if (i == 0) ierr = 1
                end if
                if (.not. no_lss .or. .not. no_route) then
                    i = check_for_file_in_list(parameter_files, 'hydrology_ini', 'MESH_parameters_hydrology.ini')
                    if (i == 0) ierr = 1
                end if
                if (.not. no_lss .and. soil_ini_flag == 5) then
                    i = check_for_file_in_list(parameter_files, 'soil_ini', 'soil.ini')
                    if (i == 0) ierr = 1
                end if
            end if
            if (btest(active_files, 2)) then

                !> CSV format.
                i = check_for_file_in_list(parameter_files, 'input_parameters_csv', 'MESH_parameters.csv')
                if (i == 0) then
                    ierr = 1
                else
                    allocate(parameter_files(i)%container, source = io_type_char_delimited(delimiter = ','))
                end if
            end if
            if (btest(active_files, 3)) then

                !> Simple text format.
                i = check_for_file_in_list(parameter_files, 'input_parameters_txt', 'MESH_parameters.txt')
                if (i == 0) then
                    ierr = 1
                else
                    allocate(parameter_files(i)%container, source = io_type_char_delimited())
                end if
            end if
            if (btest(active_files, 1)) then

                !> ASCII R2C format.
                i = check_for_file_in_list(parameter_files, 'input_parameters_r2c', 'MESH_parameters.r2c')
                if (i == 0) then
                    ierr = 1
                else
                    allocate(parameter_files(i)%container, source = io_type_r2c())
                end if
            end if
            if (btest(active_files, 4)) then

                !> NetCDF format.
#ifdef NETCDF
                i = check_for_file_in_list(parameter_files, 'input_parameters_nc', 'MESH_parameters.nc')
                if (i == 0) then
                    ierr = 1
                else
                    allocate(parameter_files(i)%container, source = io_type_nc())
                end if
#else
                call print_error( &
                    "The format of the input file is specified as NetCDF but the module is not active. " // &
                    "A version of MESH compiled with the NetCDF library must be used to read files in this format.")
                error_status = 1
                return
#endif
            end if

            !> Check for errors.
            if (ierr /= 0) then
                call print_error("An error occurred updating the file list for '" // flag_name // "'.")
                return
            end if
        end if

    end subroutine

    subroutine read_run_options2(full_path, release, error_status)

        !> Modules.
        use mpi_module, only: ipid
        use mesh_legacy_io, only: read_run_options_ini_to_field_list

        !> Input/output variables.
        character(len = *), intent(in) :: full_path
        character(len = *), intent(in) :: release
        integer, intent(out) :: error_status

        !> Local variables.
        type(io_file) run_options
        integer i, ierr
        character(len = SHORT_FIELD_LENGTH), allocatable :: values(:)
        character(len = SHORT_FIELD_LENGTH) line

        !> Create new print section.
        call print_configuring('Simulation configuration...')

        !> Read file.
        call read_run_options_ini_to_field_list(run_options, full_path, input_field_list, error_status)
        if (error_status /= 0) return

        !> Find and set the generic output folder.
        GENERIC_OUTPUT_FOLDER = repeat(' ', LONG_FIELD_LENGTH)
        call get_field_value(input_field_list, (/'OUTPUTDIR'/), GENERIC_OUTPUT_FOLDER)

        !> Find and assign other I/O options.
        call get_field_value(input_field_list, (/PRINTSIMSTATUS_NAME, 'VERBOSEMODE   '/), values)
        if (allocated(values)) then
            call parse_options(PRINTSIMSTATUS_NAME, values)
            deallocate(values)
        end if
        call get_field_value(input_field_list, (/DIAGNOSEMODE_NAME/), values)
        if (allocated(values)) then
            call parse_options(DIAGNOSEMODE_NAME, values)
            deallocate(values)
        end if
        call get_field_value(input_field_list, (/ECHOTXTMODE_NAME // '     ', 'MODELINFOOUTFLAG'/), values)
        if (allocated(values)) then
            call parse_options(ECHOTXTMODE_NAME, values)
            deallocate(values)
        end if

        !> Set the output folder to the current folder if none was defined.
        GENERIC_OUTPUT_FOLDER = trim(GENERIC_OUTPUT_FOLDER)
        if (len(GENERIC_OUTPUT_FOLDER) == 0) GENERIC_OUTPUT_FOLDER = '.'

        !> Check that the output folder exists.
        write(line, *) ipid
        open(100, file = GENERIC_OUTPUT_FOLDER // '/tmp' // trim(adjustl(line)), status = 'unknown', iostat = error_status)
        if (error_status /= 0) then
            ECHOTXTMODE = .false.
            call print_error("The output folder does not exist: " // GENERIC_OUTPUT_FOLDER)
            return
        else

            !> Delete the temporary file.
            close(100, status = 'delete')

            !> Print an output summary.
            call print_output_summary()
        end if

        !> Open the status file.
        call open_echo_txt(GENERIC_OUTPUT_FOLDER // '/MESH_output_echo_print.txt', error_status)
        if (error_status /= 0) return

        !> Write MESH version to file.
        call reset_tab()
        call print_echo_txt(trim(release))
        call print_echo_txt('')

    end subroutine

    subroutine read_basin_info(error_status)

        !> Modules.
        use mesh_legacy_io, only: read_soil_levels_txt_to_field_list
        use variable_names, only: VN_DZSOL

        !> Input/output variables.
        integer, intent(out) :: error_status

        !> Local variables.
        integer j, i, ierr
        logical no_lss

        !> Create new print section.
        call print_configuring('Static characterization and basin information...')

        !> Read file format of basin info file.
        call parse_basin_info_flag(error_status = error_status)
        if (error_status /= 0) return

        !> Check for dependencies for soil levels input file.
        no_lss = is_value_in_field(input_field_list, (/'RUNMODE'/), (/'nolss ', 'no_lss'/))
        if (.not. no_lss) then

            !> Add soil levels input file.
            i = check_for_file_in_list(basin_info_files, 'soil_levels_txt', 'MESH_input_soil_levels.txt')
        else

            !> Create a dummy soil file to prevent segmentation faults.
!-            call expand_field_list(input_field_list, increment = 2, error_status = ierr)
!-            associate (input_field => input_field_list(size(input_field_list) - 1))
!-                input_field%label = DIM_NAME_NSOL
!-                allocate(input_field%field, source = model_variable_int(dat = 1))
!-            end associate
!-            associate (input_field => input_field_list(size(input_field_list)))
!-                input_field%label = VN_DZSOL
!-                allocate(input_field%field, source = model_variable_real1d(dat = (/0.0/)))
!-                allocate(input_field%dim_names(1))
!-                input_field%dim_names = DIM_NAME_SOL
!-            end associate
            i = check_for_file_in_list(basin_info_files, 'soil_levels_virtual', '(creating virtual file for inactive soil profile)')
            allocate(basin_info_files(i)%container, source = io_type_virtual())
            allocate(basin_info_files(i)%fields(2))
            associate (input_field => basin_info_files(i)%fields(1))
                input_field%label = DIM_NAME_NSOL
                allocate(input_field%field, source = model_variable_int(dat = 1))
            end associate
            associate (input_field => basin_info_files(i)%fields(2))
                input_field%label = trim(VN_DZSOL) // ' ' // trim(DIM_NAME_LAYER) // '=1'
                allocate(input_field%field, source = model_variable_real1d(dat = (/0.0/)))
                allocate(input_field%dim_names(1))
                input_field%dim_names = DIM_NAME_BASIN
            end associate
        end if

        !> Read fields from the files.
        do j = 1, size(basin_info_files)

            !> Open and read fields from the file.
            ierr = 0
            select case (basin_info_files(j)%label)

                !> Special cases (i.e., legacy file formats).
                case ('soil_levels_txt')
                    call read_soil_levels_txt_to_field_list( &
                        basin_info_files(j), global_list = input_field_list, error_status = ierr)

                !> Generic case.
                case default

                    !> Open the file.
                    call open_input_file(basin_info_files(j), error_status = ierr)
                    if (ierr == 0) then

                        !> Read fields.
                        call read_file_fields_to_buffer(basin_info_files(j), error_status = ierr)

                        !> Copy fields to the global list.
                        call copy_field_list(input_field_list, basin_info_files(j)%fields, ierr)
                    end if
            end select
            if (ierr /= 0) error_status = ierr
        end do

    end subroutine

    subroutine read_parameter_files(fls, error_status)

        !> Modules.
        use mesh_legacy_io, only: read_parameters_class_ini_to_field_list, read_parameters_hydrology_ini_to_field_list
!temp
use model_files_variabletypes, only: fl_ids
use shd_variables, only: ShedGridParams

        !> Input/output variables.
        integer, intent(out) :: error_status

        !> Local variables.
        integer j, ierr
!temp
type(fl_ids) fls
type(ShedGridParams) shd

        !> Create new print section.
        call print_configuring('Parameter input files...')

        !> Read the parameter input file flag.
        call parse_parameter_input_flag(error_status = error_status)
        if (error_status /= 0) return

        !> Add resume files to the input file list.
        call resumerun_read(fls, shd, error_status)
        if (error_status /= 0) return

        !> Initial values.

        !> Read fields from the files.
        do j = 1, size(parameter_files)

            !> Open and read fields from the file.
            ierr = 0
            select case (parameter_files(j)%label)

                !> Special cases (i.e., legacy file formats).
                case ('class_ini')
                    call read_parameters_class_ini_to_field_list(input_field_list, parameter_files(j), error_status = ierr)
                case ('hydrology_ini')
                    call read_parameters_hydrology_ini_to_field_list( &
                        input_field_list, parameter_files(j), error_status = ierr)
                case ('soil_ini')

                !> Generic case.
                case default

                    !> Open the file.
                    call open_input_file(parameter_files(j), error_status = ierr)
                    if (ierr == 0) then

                        !> Read fields.
                        call read_file_fields_to_buffer(parameter_files(j), error_status = ierr)
                        if (ierr == 0) then

                            !> Copy fields to the global list.
                            call copy_field_list(input_field_list, parameter_files(j)%fields, ierr)
                        end if
                    end if
            end select
            if (ierr /= 0) error_status = ierr
        end do

    end subroutine

end module
