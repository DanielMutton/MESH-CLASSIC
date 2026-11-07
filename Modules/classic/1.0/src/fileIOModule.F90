!> \file
!> The File IO Module contains wrappers for accessing netCDF files.
!! It uses precompile directives to determine if a compilation should be serial, using netCDF or parallel, using MPI/PnetCDF.
!!
!> @author E. Wisernig, J. Melton

module fileIOModule
#if PARALLEL
    use mpi
#endif
  use netcdf
  implicit none

  public  :: ncCreate           !> Creates a new netCDF file and returns the file id.
  public  :: ncOpen             !> Opens an existing netCDF file and returns the file id.
  public  :: ncGetVarId         !> Returns the variable id for a given variable label.
  public  :: ncGetVarDimensions !> Returns the variable dimensions.
  public  :: ncGetVarName       !> Returns the variable name.
  public  :: ncGetDimId         !> Returns the dimension ID.
  public  :: ncGetDimLen        !> Returns the dimension length for a given dimension
  public  :: ncDefDim           !> Defines a new dimension and returns the dimension id
  public  :: ncDefVar           !> Defines a new variable and returns the variable id
  public  :: ncEndDef           !> Ends the definition mode.
  public  :: ncReDef            !> Enters definition mode.
  public  :: ncClose            !> Closes a given file.
  public  :: ncPutAtt           !> Writes an attribute to a file. Can take in attributes of the char/int/real :: types. Only one type can be used at a time.
  public  :: ncGetVar           !> Returns the variable content
  public  :: ncGetDimValues     !> Returns the values stored in a dimension (e.g. get Lon, Lat or Time values)
  public  :: ncGet1DVar         !> Returns a 1D array from a variable, based on file id, label and coordinates
  public  :: ncGet2DVar         !> Returns a 2D array from a variable, based on file id, label and coordinates
  public  :: ncGet3DVar         !> Returns a 3D array from a variable, based on file id, label and coordinates
  public  :: ncGet4DVar         !> Writes the 1D dimension values to a given file
  public  :: ncPutDimValues     !> Writes the 1D dimension values to a given file
  public  :: ncPutVar           !> Writes a local variable (1D input values, either real :: or int). Takes in a 1D array (of
  ! either real :: or int elements) and writes it into a netCDF structure, according to start and count
  public  :: ncPut2DVar         !> Writes a local variable of 2D input values. Takes in a 2D array (of either real
  ! or int elements) and writes it into a netCDF structure, according to start and count
  public  :: ncPut3DVar         !> Writes a local variable of 3D input values. Takes in a 3D array (of either real
  ! or int elements) and writes it into a netCDF structure, according to start and count
  public  :: estimateOnes       !> Returns an estimate count of how many 1s can be collapsed in the array
  public  :: collapseOnes       !> Returns an array that keeps only the values not equal to one.
  public  :: checkNC            !> Checks for errors in the NetCDF access process

contains
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileiomodule_ncCreate
  !! @{
  !> Creates a new netCDF file and returns the file id.
  integer function ncCreate(fileName, cmode)
    implicit none
    character( * ), intent(in)    :: fileName !< File name
    integer, intent(in)         :: cmode    !< Creation mode
#if PARALLEL
        integer                     :: mode, status, info
#endif

#if PARALLEL
        ! we assume MPI_COMM_WORLD and MPI_INFO_NULL are common
        call MPI_INFO_CREATE(info, status)
        ! call MPI_INFO_SET(info, "IBM_largeblock_io", "true", status) ! Obsolete and only for specific filesystems. EC, Sep, 2018.
        mode = ior(NF90_MPIIO, NF90_CLOBBER)
        mode = ior(mode, NF90_NETCDF4)
        call checkNC(nf90_create_par(trim(fileName), mode, MPI_COMM_WORLD, info, ncCreate), tag = 'ncCreate(' // trim(filename) // ') ')
#else
        call checkNC(nf90_create(trim(fileName), cmode, ncCreate), tag = 'ncCreate(' // trim(filename) // ') ')
#endif
  end function ncCreate
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileiomodule_ncOpen
  !! @{
  !> Opens an existing netCDF file and returns the file id.
  integer function ncOpen(fileName, omode)
    implicit none
    character( * ), intent(in)    :: fileName !< File name
    integer, intent(in)         :: omode    !< Open mode
#if PARALLEL
        integer                     :: mode, status, info
#endif

#if PARALLEL
        ! we assume MPI_COMM_WORLD and MPI_INFO_NULL are common
        call MPI_INFO_CREATE(info, status)
        ! call MPI_INFO_SET(info, "IBM_largeblock_io", "true", status) ! Obsolete and only for specific filesystems (EC, Sep 2018).
        mode = ior(NF90_MPIIO, omode)
        mode = ior(mode, NF90_NETCDF4)
        call checkNC(nf90_open_par(trim(fileName), mode, MPI_COMM_WORLD, info, ncOpen), tag = 'ncOpen(' // trim(filename) // ') ')
#else
        call checkNC(nf90_open(trim(fileName), omode, ncOpen), tag = 'ncOpen(' // trim(filename) // ') ')
#endif
  end function ncOpen
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileiomodule_ncGetVarId
  !! @{
  !> Returns the variable id for a given variable label.
  integer function ncGetVarId (fileId, label)
    implicit none
    integer, intent(in)         :: fileId   !< File id
    character( * ), intent(in)    :: label    !< netCDF variable label
    call checkNC(nf90_inq_varid(fileId, label, ncGetVarId), tag = 'ncGetVarId(' // trim(label) // ') hint - check if nameInCode is in the xml file, or perhaps you have a duplicate')
  end function ncGetVarId
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileiomodule_ncGetVarDimensions
  !! @{
  !> Returns the variable dimensions.
  integer function ncGetVarDimensions (fileId, varId)
    implicit none
    integer, intent(in)         :: fileId   !< File id
    integer, intent(in)         :: varId    !< Variable id
    call checkNC(nf90_inquire_variable(fileId, varId, ndims = ncGetVarDimensions), tag = 'ncGetVarDimensions() ')
  end function ncGetVarDimensions
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileiomodule_ncGetVarName
  !! @{
  !> Returns the variable name.
  function ncGetVarName (fileId)
    implicit none
    integer, intent(in)         :: fileId   !< File id
    integer                     :: varId    !< Variable id
    character(80)               :: ncGetVarName !< variable name
    logical                     :: notit !< true if not the variable we are seeking

    ! Cycle through the variables in the file skipping over the lon, lat, lev, and time variables. This assumes
    ! that the input has no other variables than the one of interest. lev is in the LUC file so be careful here if the
    ! name is different.
    varId = 1
    notit = .true.
    do while (notit)
      call checkNC(nf90_inquire_variable(fileId, varId, name = ncGetVarName), tag = 'ncGetVarName() ')
      varId = varId + 1
      if (trim(ncGetVarName) == 'lat' .or. trim(ncGetVarName) == 'lon' .or. &
          trim(ncGetVarName) == 'time' .or. trim(ncGetVarName) == 'lev') then
        notit = .true.
      else
        notit = .false.
      end if
    end do

  end function ncGetVarName
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileiomodule_ncGetDimId
  !! @{
  !> Returns the dimension ID.
  integer function ncGetDimId (fileId, label)
    implicit none
    integer, intent(in)         :: fileId   !< File id
    character( * ), intent(in)    :: label    !< Variable label
    call checkNC(nf90_inq_dimid(fileId, label, ncGetDimId), tag = 'ncGetDimId(' // trim(label) // ') ')
  end function ncGetDimId
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileiomodule_ncGetDimLen
  !! @{
  !> Returns the dimension length for a given dimension
  integer function ncGetDimLen (fileId, label)
    implicit none
    integer, intent(in)             :: fileId   !< File id
    character( * ), intent(in)        :: label    !< Dimension label
    integer                         :: dimId    !< Dimension id
    dimId = ncGetDimId(fileId, label)
    call checkNC(nf90_inquire_dimension(fileId, dimId, len = ncGetDimLen), tag = 'ncGetDimLen(' // trim(label) // ') ')
  end function ncGetDimLen
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileiomodule_ncDefDim
  !! @{
  !> Defines a new dimension and returns the dimension id
  integer function ncDefDim (fileId, label, length)
    implicit none
    integer, intent(in)        :: fileId   !< File id
    integer, intent(in)        :: length   !< Dimension length
    character( * ), intent(in) :: label    !< Label
#if PARALLEL
    integer                    :: locaLength
#endif

#if PARALLEL
        if (label == 'time') then
            ! locaLength = NFMPI_UNLIMITED !< Special case
            locaLength = 0 !< Special case
        else
            locaLength = length
        end if
#endif
    call checkNC(nf90_def_dim(fileId, label, length, ncDefDim), tag = 'ncDefDim(' // trim(label) // ') ')
  end function ncDefDim
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileiomodule_ncDefVar
  !! @{
  !> Defines a new variable and returns the variable id
  integer function ncDefVar (fileId, label, type, dimIds)
    implicit none
    integer, intent(in)                         :: fileId   !< File id
    integer, intent(in)                         :: dimIds(:)!< Dimension ids
    integer, intent(in)                         :: type     !< Variable type
    character( * ), intent(in)                    :: label    !< Variable label
    call checkNC(nf90_def_var(fileId, label, type, dimIds, ncDefVar), tag = 'ncDefVar(' // trim(label) // ') ')
  end function ncDefVar
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileiomodule_ncEndDef
  !! @{
  !> Ends the definition mode.
  subroutine ncEndDef (fileId)
    implicit none
    integer, intent(in)                 :: fileId   !< File id
    call checkNC(nf90_enddef(fileId))
  end subroutine ncEndDef
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncReDef
  !! @{
  !> Enters definition mode.
  subroutine ncReDef (fileId)
    implicit none
    integer, intent(in)                 :: fileId   !< File id
    call checkNC(nf90_redef(fileId))
  end subroutine ncReDef
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncClose
  !! @{
  !> Closes a given file.
  subroutine ncClose (fileId)
    implicit none
    integer, intent(in)                 :: fileId !< File id
    call checkNC(nf90_close(fileId))
  end subroutine ncClose
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncPutAtt
  !! @{
  !> Writes an attribute to a file. Can take in attributes of the char/int/real :: types. Only one type can be used at a time.
  subroutine ncPutAtt (fileId, varId, label, charValues, intValues, realValues)
    implicit none
    integer, intent(in)     :: fileId   !< File id
    integer, intent(in)     :: varId    !< Variable id

    character( * )            :: label
    character( * ), optional  :: charValues
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
    if (counter /= 1) stop ('In function ncPutAtt, please supply either charValues, intValue or realValues - just one')
  end subroutine ncPutAtt
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncGetVar
  !! @{
  !> Returns the variable content
  function ncGetVar (fileId, label, start, count)
    implicit none
    integer, intent(in)                             :: fileId       !< File id
    character( * ), intent(in)                        :: label        !< Label
    integer, intent(in)                             :: start(:)     !< Start array
    integer, intent(in)                             :: count(:)     !< Count array
    real, allocatable                               :: ncGetVar(:)  !< Return type
    integer                                         :: varId, ndims
    character(50)                                   :: string
    real, allocatable                               :: temp1D(:), temp2D(:,:), temp3D(:,:,:), temp4D(:,:,:,:), temp5D(:,:,:,:,:)

    varId = ncGetVarId(fileId, label)
    ndims = ncGetVarDimensions(fileId, varId)

    write(string, '(5I10)') start

    ! Currently makes I/O worse in general, but may be useful in the future (EC, Sep 2018).
    ! call checkNC(nf90_var_par_access(fileId, varId, nf90_collective))

    select case (ndims)
    case (1)
      allocate(ncGetVar(count(1)))
      allocate(temp1D(count(1)))
      call checkNC(nf90_get_var(fileId, varId, temp1D, start = start, count = count), tag = 'ncGetVar(' // trim(label) // '; start = ' // trim(string) // ') ')
      ncGetVar = temp1D
    case (2)
      allocate(ncGetVar(count(1) * count(2)))
      allocate(temp2D(count(1), count(2)))
      call checkNC(nf90_get_var(fileId, varId, temp2D, start = start, count = count), tag = 'ncGetVar(' // trim(label) // '; start = ' // trim(string) // ') ')
      ncGetVar = reshape(temp2D, (/size(temp2D)/))
    case (3)
      allocate(ncGetVar(count(1) * count(2) * count(3)))
      allocate(temp3D(count(1), count(2), count(3)))
      call checkNC(nf90_get_var(fileId, varId, temp3D, start = start, count = count), tag = 'ncGetVar(' // trim(label) // '; start = ' // trim(string) // ') ')
      ncGetVar = reshape(temp3D, (/size(temp3D)/))
    case (4)
      allocate(ncGetVar(count(1) * count(2) * count(3) * count(4)))
      allocate(temp4D(count(1), count(2), count(3), count(4)))
      call checkNC(nf90_get_var(fileId, varId, temp4D, start = start, count = count), tag = 'ncGetVar(' // trim(label) // '; start = ' // trim(string) // ') ')
      ncGetVar = reshape(temp4D, (/size(temp4D)/))
    case (5)
      allocate(ncGetVar(count(1) * count(2) * count(3) * count(4) * count(5)))
      allocate(temp5D(count(1), count(2), count(3), count(4), count(5)))
      call checkNC(nf90_get_var(fileId, varId, temp5D, start = start, count = count), tag = 'ncGetVar(' // trim(label) // '; start = ' // trim(string) // ') ')
      ncGetVar = reshape(temp5D, (/size(temp5D)/))
    case default
      stop ("Only up to 5 dimensions have been implemented !")
    end select

  end function ncGetVar
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncGetDimValues
  !! @{
  !> Returns the values stored in a dimension (e.g. get Lon, Lat or Time values)
  function ncGetDimValues (fileId, label, start, count, start2D, count2D)
    implicit none
    integer, intent(in)                         :: fileId   !< File id
    character( * ), intent(in)                    :: label    !< Label
    integer                                     :: varId, ndims
    integer, intent(in), optional               :: start(1), count(1)
    integer, intent(in), optional               :: start2D(2), count2D(2)
    real, allocatable                           :: ncGetDimValues(:)
    integer                                     :: localCount(1) = [1], localStart(1) = [1]
    integer                                     :: localCount2D(2) = [1, 1], localStart2D(2) = [1, 1]

    varId = ncGetVarId(fileId, label)
    ndims = ncGetVarDimensions(fileId, varId)

    select case (ndims)
    case (1)
      if (present(start)) localStart = start
      if (present(count)) localCount = count
      allocate(ncGetDimValues(count(1)))
      ncGetDimValues = ncGetVar(fileId, label, localStart, localCount)
    case (2)
      if (present(start2D)) localStart2D = start2D
      if (present(count2D)) localCount2D = count2D
      allocate(ncGetDimValues(count2D(1) * count2D(2)))
      ncGetDimValues = ncGetVar(fileId, label, localStart2D, localCount2D)
    case default
      stop ("Only up to 2 dimensions have been implemented in ncGetDimValues")
    end select

  end function ncGetDimValues
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncGet1DVar
  !! @{
  !> Returns a 1D array from a variable, based on file id, label and coordinates
  function ncGet1DVar (fileId, label, start, count)
    implicit none
    integer, intent(in)                         :: fileId
    character( * ), intent(in)                    :: label
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
    end if

    allocate(ncGet1DVar(localFormat(1)))
    ncGet1DVar = ncGetVar(fileId, label, start, localCount)

  end function ncGet1DVar
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncGet2DVar
  !! @{
  !> Returns a 2D array from a variable, based on file id, label and coordinates
  function ncGet2DVar (fileId, label, start, count, format)
    implicit none
    integer, intent(in)                         :: fileId
    character( * ), intent(in)                    :: label
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
      end if
    else
      allocate(localCount(3))
      localCount = [1,1,1]

      if (present(format)) then
        allocate(localFormat(size(format)))
        localFormat = format
      else
        allocate(localFormat(2))
        localFormat = [1,1]
      end if
    end if

    fixedFormat = localFormat
    allocate(ncGet2DVar(fixedFormat(1), fixedFormat(2)))
	!print *, "fileId ", fileId
	!print *, "label ", label
	!print *, "start ", start
	!print *, "localCount ", localCount
    ncGet2DVar = reshape(ncGetVar(fileId, label, start, localCount), fixedFormat)
  end function ncGet2DVar
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncGet3DVar
  !! @{
  !> Returns a 3D array from a variable, based on file id, label and coordinates
  function ncGet3DVar (fileId, label, start, count, format)
    implicit none
    integer, intent(in)                         :: fileId
    character( * ), intent(in)                    :: label
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
      end if
    else
      allocate(localCount(4))
      localCount = [1,1,1,1]

      if (present(format)) then
        allocate(localFormat(size(format)))
        localFormat = format
      else
        allocate(localFormat(3))
        localFormat = [1,1,1]
      end if
    end if

    fixedFormat = localFormat
    allocate(ncGet3DVar(fixedFormat(1), fixedFormat(2), fixedFormat(3)))
    ncGet3DVar = reshape(ncGetVar(fileId, label, start, localCount), fixedFormat, order = [1,3,2])

  end function ncGet3DVar
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncGet4DVar
  !! @{
  !> Returns a 4D array from a variable, based on file id, label and coordinates
  function ncGet4DVar (fileId, label, start, count, format)
    implicit none
    integer, intent(in)                         :: fileId
    character( * ), intent(in)                    :: label
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
      end if
    else
      allocate(localCount(5))
      localCount = [1,1,1,1,1]

      if (present(format)) then
        allocate(localFormat(size(format)))
        localFormat = format
      else
        allocate(localFormat(4))
        localFormat = [1,1,1,1]
      end if
    end if

    fixedFormat = localFormat
    allocate(ncGet4DVar(fixedFormat(1), fixedFormat(2), fixedFormat(3), fixedFormat(4)))
    ncGet4DVar = reshape(ncGetVar(fileId, label, start, localCount), fixedFormat)
  end function ncGet4DVar
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncPutDimValues
  !! @{
  !> Writes the 1D dimension values to a given file
  subroutine ncPutDimValues (fileId, label, realValues, intValues, start, count)
    implicit none
    integer, intent(in)                 :: fileId       !< File id
    character( * ), intent(in)            :: label        !< Label
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

    if (counter /= 1) stop ('In function ncPutVar, please supply either intValues or realValues - just one')

  end subroutine ncPutDimValues
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncPutVar
  !! @{
  !> Writes a local variable (1D input values, either real :: or int)
  !! Takes in a 1D array (of either real :: or int elements) and writes it into a netCDF structure, according to start and count
  subroutine ncPutVar (fileId, label, realValues, intValues, start, count)
    implicit none
    integer, intent(in)                                     :: fileId           !< File id
    character( * ), intent(in)                                :: label            !< Label
    real, intent(in), optional                              :: realValues(:)    !< Array of reals
    integer, intent(in), optional                           :: intValues(:)     !< Array of ints
    integer, intent(in)                                     :: start(:)         !< Start array
    integer, intent(in)                                     :: count(:)         !< Count array
    integer                                                 :: varId, counter

    counter = 0
    varId = ncGetVarId(fileId,label)

    if (present(realValues)) then
      counter = counter + 1
      call checkNC(nf90_put_var(fileId, varId, realValues, start, count), tag = 'ncPutVar(' // trim(label) // ') ')
    else if (present(intValues)) then
      counter = counter + 1
      call checkNC(nf90_put_var(fileId, varId, intValues, start, count), tag = 'ncPutVar(' // trim(label) // ') ')
    end if

    if (counter /= 1) stop ('In function ncPutVar, please supply either intValues or realValues - just one')

  end subroutine ncPutVar
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncPut2DVar
  !! @{
  !> Writes a local variable of 2D input values
  !! Takes in a 2D array (of either real :: or int elements) and writes it into a netCDF structure, according to start and count
  subroutine ncPut2DVar (fileId, label, realValues, intValues, start, count)
    implicit none
    integer, intent(in)                                     :: fileId   !< File id
    character( * ), intent(in)                                :: label    !< Label
    real, intent(in), optional                              :: realValues(:,:)    !< Array of reals
    integer, intent(in), optional                           :: intValues(:,:)     !< Array of ints
    integer, dimension(:), intent(in)                       :: start    !< Start array
    integer, dimension(3), intent(in)                       :: count    !< Count array
    integer                                                 :: varId, counter
    varId = ncGetVarId(fileId, label)

    counter = 0
    if (present(realValues)) then
      counter = counter + 1
      call checkNC(nf90_put_var(fileId, varId, reshape(realValues, count), start, count), tag = 'ncPut2DVar(' // trim(label) // ') ')
    else if (present(intValues)) then
      counter = counter + 1
      call checkNC(nf90_put_var(fileId, varId, reshape(intValues, count), start, count), tag = 'ncPut2DVar(' // trim(label) // ') ')
    end if

    if (counter /= 1) stop ('In function ncPut2DVar, please supply either intValues or realValues - just one')
  end subroutine ncPut2DVar
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_ncPut3DVar
  !! @{
  !> Writes a local variable of 3D input values
  !! Takes in a 3D array (of either real :: or int elements) and writes it into a netCDF structure, according to start and count
  subroutine ncPut3DVar (fileId, label, realValues, intValues, start, count)
    implicit none
    integer, intent(in)                                     :: fileId   !< File id
    character( * ), intent(in)                                :: label    !< Label
    real, intent(in), optional                              :: realValues(:,:,:)    !< Array of reals
    integer, intent(in), optional                           :: intValues(:,:,:)     !< Array of ints
    integer, dimension(:), intent(in)                       :: start    !< Start array
    integer, dimension(4), intent(in)                       :: count    !< Count array
    integer                                                 :: varId, counter
    varId = ncGetVarId(fileId, label)

    counter = 0
    if (present(realValues)) then
      counter = counter + 1
      call checkNC(nf90_put_var(fileId, varId, reshape(realValues, count, order = [1,2,4,3]), start, count), tag = 'ncPut3DVar(' // trim(label) // ') ')
    else if (present(intValues)) then
      counter = counter + 1
      call checkNC(nf90_put_var(fileId, varId, reshape(intValues, count, order = [1,2,4,3]), start, count), tag = 'ncPut3DVar(' // trim(label) // ') ')
    end if
    if (counter /= 1) stop ('In function ncPut3DVar, please supply either intValues or realValues - just one')
  end subroutine ncPut3DVar
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_estimateOnes
  !! @{
  !> Returns an estimate count of how many 1s can be collapsed in the array
  integer function estimateOnes (values)
    implicit none
    integer, intent(in)             :: values(:)        !< The input values
    integer                         :: i

    estimateOnes = 0
    do i = 1,size(values)
      if (values(i) /= 1) estimateOnes = estimateOnes + 1
    end do

    if (estimateOnes == 0) then
      ! print*,('The estimateOnes function found no values different than 1 and will default to just one element')
      estimateOnes = 1
    end if
  end function estimateOnes
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_collapseOnes
  !! @{
  !> Returns an array that keeps only the values not equal to one.
  function collapseOnes (values)
    implicit none
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
      end if
    end do
    if (counter == 0) then
      ! print*,('The collapseOnes function found no values different than 1 and will default to just one element')
      collapseOnes(1) = 1
    end if
  end function collapseOnes
  !! @}
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> \ingroup fileIOMOdule_checkNC
  !! @{
  !> Checks for errors in the NetCDF access process
  subroutine checkNC (ncStatus, tag)
    
    use generalUtils, only : abandonCell
    
    implicit none
    integer, intent(in)         :: ncStatus !< Status variable
    character( * ), optional      :: tag  !< Optional tag
    character(150)              :: message

    if (present(tag)) then
      message = tag
    else
      message = 'Unspecified tag'
    end if

    if (ncStatus /= nf90_noerr) then
      print*,'netCDF error with tag ', trim(message), ' : ', trim(nf90_strerror(ncStatus))
      print*,'Stopping the run for this cell'
      call abandonCell
    end if
    
  end subroutine checkNC
  !! @}

  !> \namespace fileiomodule
  !! Contains wrappers for accessing netCDF files thus allowing the same code for MPI and serial
end module fileIOModule
