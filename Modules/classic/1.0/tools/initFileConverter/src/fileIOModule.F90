!> The File IO Module contains wrappers for accessing netCDF files.
!! It uses precompile directives to determine if a compilation should be serial, using netCDF or parallel, using MPI/PnetCDF.
module fileIOModule
#if PARALLEL
    use mpi
#endif
    use netcdf
    implicit none
contains

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Creates a new netCDF file and returns the file id.
    integer function ncCreate(fileName, cmode)
        character(*), intent(in)    :: fileName !< File name
        integer, intent(in)         :: cmode    !< Creation mode
#if PARALLEL
        integer                     :: mode, status, info
#endif

#if PARALLEL
        ! we assume MPI_COMM_WORLD and MPI_INFO_NULL are common
        call MPI_INFO_CREATE(info, status)
        call MPI_INFO_SET(info, "IBM_largeblock_io", "true", status)
        mode = 0
        mode = ior(NF90_MPIIO, NF90_CLOBBER)
        mode = ior(mode, NF90_NETCDF4)
        call checkNC(nf90_create_par(trim(fileName), mode, MPI_COMM_WORLD, info, ncCreate), tag = 'ncCreate(' // trim(filename) // ') ')
#else
        call checkNC(nf90_create(trim(fileName), cmode, ncCreate), tag = 'ncCreate(' // trim(filename) // ') ')
#endif
    end function ncCreate

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Opens an existing netCDF file and returns the file id.
    integer function ncOpen(fileName, omode)
        character(*), intent(in)    :: fileName !< File name
        integer, intent(in)         :: omode    !< Open mode
#if PARALLEL
        integer                     :: mode, status, info
#endif

#if PARALLEL
        ! we assume MPI_COMM_WORLD and MPI_INFO_NULL are common
        call MPI_INFO_CREATE(info, status)
        call MPI_INFO_SET(info, "IBM_largeblock_io", "true", status)
        mode = 0
        mode = ior(NF90_MPIIO, omode)
        mode = ior(mode, NF90_NETCDF4)
        call checkNC(nf90_open_par(trim(fileName), mode, MPI_COMM_WORLD, info, ncOpen), tag = 'ncOpen(' // trim(filename) // ') ')
#else
        call checkNC(nf90_open(trim(fileName), omode, ncOpen), tag = 'ncOpen(' // trim(filename) // ') ')
#endif
    end function ncOpen

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns the variable id for a given file id a variable label.
    integer function ncGetVarId(fileId, label)
        integer, intent(in)         :: fileId   !< File id
        character(*), intent(in)    :: label    !< netCDF variable label
        call checkNC(nf90_inq_varid(fileId, label, ncGetVarId), tag = 'ncGetVarId(' // trim(label) // ') ')
    end function ncGetVarId

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns the variable dimensions.
    integer function ncGetVarDimensions(fileId, varId)
        integer, intent(in)         :: fileId   !< File id
        integer, intent(in)         :: varId    !< Variable id
        call checkNC(nf90_inquire_variable(fileId, varId, ndims = ncGetVarDimensions), tag = 'ncGetVarDimensions() ')
    end function ncGetVarDimensions

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns the dimension ID.
    integer function ncGetDimId(fileId, label)
        integer, intent(in)         :: fileId   !< File id
        character(*), intent(in)    :: label    !< Variable label
        call checkNC(nf90_inq_dimid(fileId, label, ncGetDimId), tag = 'ncGetDimId(' // trim(label) // ') ')
    end function ncGetDimId

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns the dimension length for a given dimension
    integer function ncGetDimLen(fileId, label)
        integer, intent(in)             :: fileId   !< File id
        character(*), intent(in)        :: label    !< Dimension label
        integer                         :: dimId    !< Dimension id
        dimId = ncGetDimId(fileId, label)
        call checkNC(nf90_inquire_dimension(fileId, dimId, len=ncGetDimLen), tag = 'ncGetDimLen(' // trim(label) // ') ')
    end function ncGetDimLen

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Defines a new dimension and returns the dimension id
    integer function ncDefDim(fileId, label, length)
        integer, intent(in)                         :: fileId   !< File id
        integer, intent(in)                         :: length   !< Dimension length
        character(*), intent(in)                    :: label    !< Label
#if PARALLEL
        integer                                     :: locaLength
#endif

#if PARALLEL
        if (label == 'time') then
            !locaLength = NFMPI_UNLIMITED !< Special case
            locaLength = 0 !< Special case
        else
            locaLength = length
        end if
#endif
        call checkNC(nf90_def_dim(fileId, label, length, ncDefDim), tag = 'ncDefDim(' // trim(label) // ') ')
    end function ncDefDim

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Defines a new variable and returns the variable id
    integer function ncDefVar(fileId, label, type, dimIds)
        integer, intent(in)                         :: fileId   !< File id
        integer, intent(in)                         :: dimIds(:)!< Dimension ids
        integer, intent(in)                         :: type     !< Variable type
        character(*), intent(in)                    :: label    !< Variable label
        call checkNC(nf90_def_var(fileId, label, type, dimIds, ncDefVar), tag = 'ncDefVar(' // trim(label) // ') ')
    end function ncDefVar

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Ends the definition mode.
    subroutine ncEndDef(fileId)
        integer, intent(in)                 :: fileId   !< File id
        call checkNC(nf90_enddef(fileId))
    end subroutine ncEndDef

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Enters definition mode.
    subroutine ncReDef(fileId)
        integer, intent(in)                 :: fileId   !< File id
        call checkNC(nf90_redef(fileId))
    end subroutine ncReDef

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Closes a given file.
    subroutine ncClose(fileId)
        integer, intent(in)                 :: fileId !< File id
        call checkNC(nf90_close(fileId))
    end subroutine ncClose

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Writes an attribute to a file. Can take in attributes of the char/int/real types. Only one type can be used at a time.
    subroutine ncPutAtt(fileId, varId, label, charValues, intValues, realValues)
        integer, intent(in)     :: fileId   !< File id
        integer, intent(in)     :: varId    !< Variable id

        character(*)            :: label
        character(*), optional  :: charValues
        integer, optional       :: intValues
        real, optional          :: realValues
        integer                 :: counter

        counter = 0
        if (present(charValues)) then
            counter = counter + 1
            call checkNC(nf90_put_att(fileId, varId, label, charValues), tag = 'ncPutAtt(' // trim(label) // ') ')
        else if (present(intValues)) then
            counter = counter + 1
            call checkNC(nf90_put_att(fileId, varId, label, intValues), tag = 'ncPutAtt(' // trim(label) // ') ')
        else if (present(realValues)) then
            counter = counter + 1
            call checkNC(nf90_put_att(fileId, varId, label, realValues), tag = 'ncPutAtt(' // trim(label) // ') ')
        end if
        if (counter /= 1) stop ('In function ncPutAtt, please supply either charValues, intValue or realValues; just one')
    end subroutine ncPutAtt

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns the variable content in the form of a data package
    function ncGetVar(fileId, label, start, count)
        integer, intent(in)                             :: fileId       !< File id
        character(*), intent(in)                        :: label        !< Label
        integer, intent(in)                             :: start(:)     !< Start array
        integer, intent(in)                             :: count(:)     !< Count array
        real, allocatable                               :: ncGetVar(:)  !< Return type
        integer                                         :: varId, ndims
        real, allocatable                               :: temp1D(:), temp2D(:,:), temp3D(:,:,:), temp4D(:,:,:,:), temp5D(:,:,:,:,:)

        varId = ncGetVarId(fileId, label)
        ndims = ncGetVarDimensions(fileId, varId)

        select case(ndims)
        case(1)
            allocate(ncGetVar(count(1)))
            allocate(temp1D(count(1)))
            call checkNC(nf90_get_var(fileId, varId, temp1D, start = start, count = count), tag = 'ncGetVar(' // trim(label) // ') ')
            ncGetVar = temp1D
        case(2)
            allocate(ncGetVar(count(1) * count(2)))
            allocate(temp2D(count(1), count(2)))
            call checkNC(nf90_get_var(fileId, varId, temp2D, start = start, count = count), tag = 'ncGetVar(' // trim(label) // ') ')
            ncGetVar = reshape(temp2D,(/size(temp2D)/))
        case(3)
            allocate(ncGetVar(count(1) * count(2) * count(3)))
            allocate(temp3D(count(1), count(2), count(3)))
            call checkNC(nf90_get_var(fileId, varId, temp3D, start = start, count = count), tag = 'ncGetVar(' // trim(label) // ') ')
            ncGetVar = reshape(temp3D,(/size(temp3D)/))
        case(4)
            allocate(ncGetVar(count(1) * count(2) * count(3) * count(4)))
            allocate(temp4D(count(1), count(2), count(3), count(4)))
            call checkNC(nf90_get_var(fileId, varId, temp4D, start = start, count = count), tag = 'ncGetVar(' // trim(label) // ') ')
            ncGetVar = reshape(temp4D,(/size(temp4D)/))
        case(5)
            allocate(ncGetVar(count(1) * count(2) * count(3) * count(4) * count(5)))
            allocate(temp5D(count(1), count(2), count(3), count(4), count(5)))
            call checkNC(nf90_get_var(fileId, varId, temp5D, start = start, count = count), tag = 'ncGetVar(' // trim(label) // ') ')
            ncGetVar = reshape(temp5D,(/size(temp5D)/))
        case default
            stop ("Only up to 5 dimensions have been implemented!")
        end select

    end function ncGetVar

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Writes the 1D dimension values to a given file
    subroutine ncPutDimValues(fileId, label, realValues, intValues, start, count)
        integer, intent(in)                 :: fileId       !< File id
        character(*), intent(in)            :: label        !< Label
        real, intent(in), optional       :: realValues(:)!< Array of reals
        integer, intent(in), optional    :: intValues(:) !< Array of ints
        integer, intent(in), optional       :: start(1)     !< Start array
        integer, intent(in), optional       :: count(1)     !< Count array
        integer                             :: localStart(1) = [1], localCount(1) = [1], varId, counter

        counter = 0

        if (present(start)) localStart = start
        if (present(count)) localCount = count

        varId = ncGetVarId(fileId, label)

        if (present(realValues)) then
            counter = counter + 1
            call checkNC(nf90_put_var(fileId, varId, realValues, localStart, localCount), tag = 'ncPutDimValues(' // trim(label) // ') ')
        else if (present(intValues)) then
            counter = counter + 1
            call checkNC(nf90_put_var(fileId, varId, intValues, localStart, localCount), tag = 'ncPutDimValues(' // trim(label) // ') ')
        end if

        if (counter /= 1) then
            print*,'In function ncPutDimValues, please supply either intValues or realValues; just one ',trim(label)
            stop
        end if

    end subroutine ncPutDimValues

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Writes a local variable (always 1D input values, either real or int)
    subroutine ncPutVar(fileId, label, realValues, intValues, start, count)
        integer, intent(in)                                     :: fileId       !< File id
        character(*), intent(in)                                :: label        !< Label
        real, intent(in), optional                           :: realValues(:)!< Array of reals
        integer, intent(in), optional                        :: intValues(:) !< Array of ints
        integer, intent(in)                                     :: start(:)     !< Start array
        integer, intent(in)                                     :: count(:)     !< Count array
        integer                                                 :: varId, counter

        counter = 0
        varId = ncGetVarId(fileId, label)

        if (present(realValues)) then
            counter = counter + 1
            call checkNC(nf90_put_var(fileId, varId, realValues, start, count), tag = 'ncPutVar(' // trim(label) // ') ')
        else
            counter = counter + 1
            call checkNC(nf90_put_var(fileId, varId, intValues, start, count), tag = 'ncPutVar(' // trim(label) // ') ')
        end if

        if (counter /= 1) stop ('In function ncPutVar, please supply either intValues or realValues; just one')

    end subroutine ncPutVar

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns the values stored in a dimension (e.g. get Lon, Lat or Time values)
    function ncGetDimValues(fileId, label, start, count)
        integer, intent(in)                         :: fileId   !< File id
        character(*), intent(in)                    :: label    !< Label
        integer, intent(in), optional               :: start(1), count(1)
        real, allocatable                           :: ncGetDimValues(:)
        integer                                     :: localCount(1) = [1], localStart(1) = [1]

        if (present(start)) localStart = start
        if (present(count)) localCount = count
        allocate(ncGetDimValues(count(1)))
        ncGetDimValues = ncGetVar(fileId, label, localStart, localCount)

    end function ncGetDimValues

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns a 1D array from a variable, based on file id, label and coordinates
    function ncGet1DVar(fileId, label, start, count)
        integer, intent(in)                         :: fileId
        character(*), intent(in)                    :: label
        integer, intent(in)                         :: start(:)
        integer, intent(in), optional               :: count(:)
        real, allocatable                           :: ncGet1DVar(:)
        integer, allocatable                        :: localCount(:), localFormat(:)
        integer                                     :: formatSize

        if (present(count)) then
            allocate(localCount(size(count)))
            localCount = count
            formatSize = estimateOnes(count)
            allocate(localFormat(formatSize))
            localFormat = collapseOnes(count)
            if (size(localFormat) /= 1) stop ('Count problem in ncGet1DVar function')
        else
            allocate(localCount(2))
            localCount = [1, 1]
            allocate(localFormat(1))
            localFormat = [1]
        endif

        allocate(ncGet1DVar(localFormat(1)))
        ncGet1DVar = ncGetVar(fileId, label, start, localCount)

    end function ncGet1DVar

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns a 2D array from a variable, based on file id, label and coordinates
    function ncGet2DVar(fileId, label, start, count, format)
        integer, intent(in)                         :: fileId
        character(*), intent(in)                    :: label
        integer, intent(in)                         :: start(:)
        integer, intent(in), optional               :: count(:), format(:)
        real, allocatable                           :: ncGet2DVar(:,:)
        integer, allocatable                        :: localCount(:), localFormat(:)
        integer                                     :: fixedFormat(2), formatSize

        if (present(count)) then
            allocate(localCount(size(count)))
            localCount = count

            if (present(format)) then
                allocate(localFormat(size(format)))
                localFormat = format
            else
                formatSize = estimateOnes(count)
                allocate(localFormat(formatSize))
                localFormat = collapseOnes(count)
                if (size(localFormat) /= 2) stop ('Count and/or format problem found in ncGet2DVar function')
            endif
        else
            allocate(localCount(3))
            localCount = [1, 1, 1]

            if (present(format)) then
                allocate(localFormat(size(format)))
                localFormat = format
            else
                allocate(localFormat(2))
                localFormat = [1, 1]
            endif
        endif

        fixedFormat = localFormat
        allocate(ncGet2DVar(fixedFormat(1), fixedFormat(2)))
        ncGet2DVar = reshape(ncGetVar(fileId, label, start, localCount), fixedFormat)
    end function ncGet2DVar

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns a 3D array from a variable, based on file id, label and coordinates
    function ncGet3DVar(fileId, label, start, count, format)
        integer, intent(in)                         :: fileId
        character(*), intent(in)                    :: label
        integer, intent(in)                         :: start(:)
        integer, intent(in), optional               :: count(:), format(:)
        real, allocatable                           :: ncGet3DVar(:,:,:)
        integer, allocatable                        :: localCount(:), localFormat(:)
        integer                                     :: fixedFormat(3), formatSize

        if (present(count)) then
            allocate(localCount(size(count)))
            localCount = count

            if (present(format)) then
                allocate(localFormat(size(format)))
                localFormat = format
            else
                formatSize = estimateOnes(count)
                allocate(localFormat(formatSize))
                localFormat = collapseOnes(count)
                if (size(localFormat) /= 3) stop ('Count and/or format problem found in ncGet3DVar function')
            endif
        else
            allocate(localCount(4))
            localCount = [1, 1, 1, 1]

            if (present(format)) then
                allocate(localFormat(size(format)))
                localFormat = format
            else
                allocate(localFormat(3))
                localFormat = [1, 1, 1]
            endif
        endif

        fixedFormat = localFormat
        allocate(ncGet3DVar(fixedFormat(1), fixedFormat(2), fixedFormat(3)))
        ncGet3DVar = reshape(ncGetVar(fileId, label, start, localCount), fixedFormat)
    end function ncGet3DVar

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns a 4D array from a variable, based on file id, label and coordinates
    function ncGet4DVar(fileId, label, start, count, format)
        integer, intent(in)                         :: fileId
        character(*), intent(in)                    :: label
        integer, intent(in)                         :: start(:)
        integer, intent(in), optional               :: count(:), format(:)
        real, allocatable                           :: ncGet4DVar(:,:,:,:)
        integer, allocatable                        :: localCount(:), localFormat(:)
        integer                                     :: fixedFormat(4), formatSize

        if (present(count)) then
            allocate(localCount(size(count)))
            localCount = count

            if (present(format)) then
                allocate(localFormat(size(format)))
                localFormat = format
            else
                formatSize = estimateOnes(count)
                allocate(localFormat(formatSize))
                localFormat = collapseOnes(count)
                if (size(localFormat) /= 4) stop ('Count and/or format problem found in ncGet4DVar function')
            endif
        else
            allocate(localCount(5))
            localCount = [1, 1, 1, 1, 1]

            if (present(format)) then
                allocate(localFormat(size(format)))
                localFormat = format
            else
                allocate(localFormat(4))
                localFormat = [1, 1, 1, 1]
            endif
        endif

        fixedFormat = localFormat
        allocate(ncGet4DVar(fixedFormat(1), fixedFormat(2), fixedFormat(3), fixedFormat(4)))
        ncGet4DVar = reshape(ncGetVar(fileId, label, start, localCount), fixedFormat)
    end function ncGet4DVar

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns an estimate count of how many 1s can be collapsed in the array
    integer function estimateOnes(values)
        integer, intent(in)             :: values(:)        !< The input values
        integer                         :: i

        estimateOnes = 0
        do i = 1, size(values)
            if (values(i) /= 1) estimateOnes = estimateOnes + 1
        enddo

        if (estimateOnes == 0) then
            !print*,('The estimateOnes function found no values different than 1 and will default to just one element');
            estimateOnes = 1;
        endif
    end function estimateOnes

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Returns an array that keeps only the values not equal to one.
    function collapseOnes(values)
        integer, intent(in)             :: values(:)        !< The input values
        integer, allocatable            :: collapseOnes(:)  !< The output values
        integer                         :: i, counter, count
        count = estimateOnes(values)
        allocate(collapseOnes(count))

        counter = 0
        do i = 1, size(values)
            if (values(i) /= 1) then
                counter = counter + 1
                collapseOnes(counter) = values(i)
            endif
        enddo
        if (counter == 0) then
            !print*,('The collapseOnes function found no values different than 1 and will default to just one element');
            collapseOnes(1) = 1
        endif
    end function collapseOnes

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    !> Checks for errors in the NetCDF access process
    subroutine checkNC(ncStatus, tag)
        integer, intent(in)         :: ncStatus !< Status variable
        character(*), optional      :: tag  !< Optional tag
        character(100)              :: message
        if (present(tag)) then
            message = tag
        else
            message = 'Unspecified tag'
        endif

        if(ncStatus /= nf90_noerr) then
            print*,'netCDF error with tag ', trim(message), ' : ', trim(nf90_strerror(ncStatus))
            stop
        end if
    end subroutine checkNC

    ! OLD FUNCTIONS TO BE DELETED
    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    ! Write a local 2D variable
#if PARALLEL
    subroutine ncPut2DVar(fileId, label, data, start, count)
        integer, intent(in)                                     :: fileId   !<
        character(*), intent(in)                                :: label
        real, dimension(:,:), intent(in)                        :: data
        integer, dimension(:), intent(in)                       :: start
        integer, dimension(:), optional, intent(in)             :: count
        integer                                                 :: varId
        real, dimension(:,:,:), allocatable                     :: temp3D
        integer, dimension(3)                                   :: localFormat, localCount
        if (present(count)) then
            localCount = count
        else
            localCount = [1, 1, 1]
        endif
        localFormat = count
        varId = ncGetVarId(fileId, label)
        temp3D = reshape(data, localFormat)
        !call checkNC(nf90mpi_put_var_all(fileId, varId, temp3D, int(start,8), int(localCount,8)))
    end subroutine ncPut2DVar
#else
    subroutine ncPut2DVar(fileId, label, data, start, count)
        integer, intent(in)                                     :: fileId
        character(*), intent(in)                                :: label
        real, dimension(:,:), intent(in)                        :: data
        integer, dimension(:), intent(in)                       :: start
        integer, dimension(:), optional, intent(in)             :: count
        integer                                                 :: varId
        real, dimension(:,:,:), allocatable                     :: temp3D
        integer, dimension(3)                                   :: localFormat, localCount
        if (present(count)) then
            localCount = count
        else
            localCount = [1, 1, 1]
        endif
        localFormat = count
        varId = ncGetVarId(fileId, label)
        temp3D = reshape(data, localFormat)
        call checkNC(nf90_put_var(fileId, varId, temp3D, start, localCount))
    end subroutine ncPut2DVar
#endif

    ! Write a local 3D variable (4D in NetCDF file)
#if PARALLEL
    subroutine ncPut3DVar(fileId, label, data, start, count)
        integer, intent(in)                                     :: fileId
        character(*), intent(in)                                :: label
        real, dimension(:,:,:), intent(in)                      :: data
        integer, dimension(:), intent(in)                       :: start, count
        integer                                                 :: varId
        real, dimension(:,:,:,:), allocatable                   :: temp4D
        integer, dimension(4)                                   :: format = 0
        format = count
        varId = ncGetVarId(fileId, label)
        temp4D = reshape(data, format)
        !call checkNC(nf90mpi_put_var_all(fileId, varId, temp4D, int(start,8), int(count,8)))
    end subroutine ncPut3DVar
#else
    subroutine ncPut3DVar(fileId, label, data, start, count)
        integer, intent(in)                                     :: fileId
        character(*), intent(in)                                :: label
        real, dimension(:,:,:), intent(in)                      :: data
        integer, dimension(:), intent(in)                       :: start, count
        integer                                                 :: varId
        real, dimension(:,:,:,:), allocatable                   :: temp4D
        integer, dimension(4)                                   :: format = 0
        format = count
        varId = ncGetVarId(fileId, label)
        temp4D = reshape(data, format)
        call checkNC(nf90_put_var(fileId, varId, temp4D, start, count))
    end subroutine ncPut3DVar
#endif

end module fileIOModule
