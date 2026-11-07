!> \file
!> Central module for all netCDF output file operations
!!
!> @author Joe Melton and Ed Wisernig
!ignoreLint(1000)
module outputManager

  use ctemStateVars, only : c_switch

  implicit none

  public :: generateOutputFiles
  private :: generateNetCDFFile
  private :: addVariable
  private :: validGroup
  private :: validTime
  private :: generateFilename
  private :: getDescriptor
  public :: getIdByKey
  public  :: createNetCDF
  private :: determineTime
  private :: identityVector
  public :: writeOutput1D
  public  :: closeNCFiles
  public :: checkForTime

  !> Stores geographic information about the total domain and the region to be simulated
  type simulationDomain
    real, dimension(:), allocatable     :: lonLandCell, latLandCell     !< Long/Lat values of only the land cells in our model domain
    integer, dimension(:), allocatable  :: lonLandIndex, latLandIndex   !< Indexes of only the land cells in our model domain for our resolution
    real, dimension(:), allocatable     :: allLonValues, allLatValues   !< All long/Lat values in our model domain (including ocean/non-land)
    integer, dimension(:), allocatable  :: lonLocalIndex, latLocalIndex  !< The index for only the region that is being simulated
    real, dimension(:), allocatable     :: latUnique, lonUnique          !< The index for only the region that is being simulated with each value only once
    integer                             :: LandCellCount    !> number of land cells that the model will run over
    real, dimension(4) :: domainBounds                      !> Corners of the domain to be simulated (netcdfs)
    integer :: srtx                                         !> starting index for this simulation for longitudes
    integer :: srty                                         !> starting index for this simulation for latitudes
    integer :: cntx                                         !> number of grid cells for this simulation in the longitude direction
    integer :: cnty                                         !> number of grid cells for this simulation in the latitude direction
  end type

  type(simulationDomain) :: myDomain

  !> Used to set up the output netcdf files.
  type outputDescriptor
    character(80)   :: group                = ''
    character(30)   :: shortName            = ''
    character(80)   :: standardName         = ''
    character(400)  :: longName             = ''        !< Long name of the variable
    character(30)   :: timeFreq             = ''        !< Time frequency of variable: half-hourly, daily, monthly, annually
    logical         :: includeBareGround    = .false.   !< If true then expand the PFT list for a final position that is the bare ground.
  end type

  type(outputDescriptor), allocatable     :: outputDescriptors(:)

  !> Contains characteristic information about the output netcdf files and is used in their creation
  type variant
    character(80)   :: nameInCode           = ''
    character(80)   :: timeFrequency        = ''
    character(80)   :: outputForm           = ''
    character(80)   :: shortName            = ''
    character(80)   :: units                = ''        !< Units of the variant
  end type

  type(variant), allocatable              :: variants(:)

  !> Stores the information needed to retrieve the output netcdf files to write to them
  type netcdfVar
    integer         :: ncid
    character(30)   :: key
    character(350)  :: filename
  end type

  integer, parameter  :: maxncVariableNumber = 300        !< Maximum number of netcdf output files to make (can be adjusted)
  type(netcdfVar)     :: netcdfVars(maxncVariableNumber)

  integer :: variableCount = 0, descriptorCount = 0, variantCount = 0 !< Initializes the counter variables

  integer         :: refyr                                     !< Time reference for netcdf output files
  character(30)   :: timestart                                 !< Time reference for netcdf output files
  real   :: fill_value = huge(fill_value) !1.E38                                 !< Default fill value for missing values in netcdf output files
  real, dimension(:), allocatable :: timeVect                  !< Array of the timesteps in days since refyr for this model run and output file

  real :: consecDays                                        !< Consecutive days since refyr in run. Set by setConsecDays

contains

  !---------------------------------------------------------------------------------------
  !> \ingroup outputmanager_generateOutputFiles
  !> @{
  !> Runs through all possible variants and calls for the generation of the required output files.
  subroutine generateOutputFiles

    implicit none

    integer             :: i

    do i = 1, variantCount
      call generateNetCDFFile(variants(i)%nameInCode, variants(i)%timeFrequency, &
                              variants(i)%outputForm, variants(i)%units, variants(i)%shortName)

    end do

    return

  end subroutine generateOutputFiles

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_generateNetCDFFile
  !> @{
  !> Generates a new (netcdf) variable
  subroutine generateNetCDFFile (nameInCode, timeFreq, outputForm, units, descriptorLabel)

    implicit none

    character( * ), intent(in)                :: nameInCode, timeFreq, descriptorLabel, units
    character( * ), intent(inout)             :: outputForm
    type(outputDescriptor)                  :: descriptor
    character(350)                           :: filename = ''
    logical                                 :: isTimeValid, isGroupValid, fileCreatedOk
    integer                                 :: id

    ! Get variable descriptor
    descriptor = getDescriptor(descriptorLabel)

    ! If the requested timeFreq matches the project config, all's good
    isTimeValid = validTime(timeFreq, descriptor)

    ! If the group property of descriptor matches the project config, all's good
    isGroupValid = validGroup(descriptor,outputForm)

    ! If the project config and variable descriptor match the request, process current variable
    if (isTimeValid .and. isGroupValid) then

      ! Generate the filename
      filename = generateFilename(outputForm, descriptor)

      ! Allocate a new variable (ncid, filename, key etc.)
      id = addVariable(nameInCode, filename)

      ! Make the netcdf file for the new variable (mostly definitions)
      call createNetCDF(filename, id, outputForm, descriptor, timeFreq, units, nameInCode)

      ! Now make sure the file was properly created
      fileCreatedOk = checkFileExists(filename)

      if (.not. fileCreatedOk) then
        print * ,'Failed to create',filename
        print * ,'Aborting'
        stop ! can use stop here as not in MPI part of code.
      end if
    end if

  end subroutine generateNetCDFFile

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_checkFileExists
  !> @{
  !> Checks if file exists
  logical function checkFileExists (filename)
    implicit none

    character( * ), intent(in) :: filename

    inquire(file = filename, exist = checkFileExists)

  end function

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_addVariable
  !> @{
  !> Adds the new variable to the list of variables (see the type "netcdfVars")
  integer function addVariable (key, filename)
    use fileIOModule
    implicit none
    character( * ), intent(in)    :: key, filename
    integer                     :: ncid

    ncid = ncCreate(fileName, cmode = NF90_CLOBBER)

    variableCount = variableCount + 1
    netcdfVars(variableCount)%ncid = ncid
    netcdfVars(variableCount)%key = key
    netcdfVars(variableCount)%filename = filename
    addVariable = variableCount

  end function addVariable

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_validGroup
  !> @{
  !> Determines if the current variable matches the project configuration
  logical function validGroup (descriptor,outputForm)

    implicit none

    character( * ), intent(in)              :: outputForm
    type(outputDescriptor), intent(in)    :: descriptor

    ! First weed out the per pft, per tile options if they should not be used.
    if (.not. c_switch%doperpftoutput .and. trim(outputForm) == "pft") then
      validGroup = .false.
      return
    end if

    if (.not. c_switch%dopertileoutput .and. trim(outputForm) == "tile") then
      validGroup = .false.
      return
    end if

    ! Now check over the remaining options
    if (trim(descriptor%group) == "class") then ! CLASS outputs always are valid
      validGroup = .true.
    else if (c_switch%ctem_on .and. trim(descriptor%group) == "ctem") then
      validGroup = .true.
    else if (c_switch%ctem_on) then ! check the CTEM sub-switches
      if (c_switch%dofire .and. trim(descriptor%group) == "fire") then
        validGroup = .true.
      else if (c_switch%lnduseon .and. trim(descriptor%group) == "land") then
        validGroup = .true.
      else if ((c_switch%lnduseon .or. c_switch%dofire) .and. trim(descriptor%group) == "land_fire") then
        validGroup = .true.
      else if (c_switch%PFTCompetition .and. trim(descriptor%group) == "PFTCompetition") then
        validGroup = .true.
      else if (c_switch%useTracer > 0 .and. trim(descriptor%group) == "tracer") then
        validGroup = .true.
      else if (c_switch%doMethane .and. trim(descriptor%group) == "methane") then
        validGroup = .true.
      else
        validGroup = .false.
      end if
    else
      validGroup = .false.
    end if
  end function validGroup

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_validTime
  !> @{
  !> Determines whether the current variable matches the project configuration
  logical function validTime (timeFreq, descriptor)

    implicit none

    type(outputDescriptor), intent(inout) :: descriptor
    character( * ), intent(in)                :: timeFreq
    logical                                 :: valid

    valid = .true.
    if (c_switch%doAnnualOutput .and. timeFreq == 'annually') then
      descriptor%timeFreq = timeFreq
    else if (c_switch%domonthoutput .and. timeFreq == 'monthly') then
      descriptor%timeFreq = timeFreq
    else if (c_switch%dodayoutput .and. timeFreq == 'daily') then
      descriptor%timeFreq = timeFreq
    else if (c_switch%dohhoutput .and. timeFreq == 'halfhourly') then
      descriptor%timeFreq = timeFreq
    else
      valid = .false.
    end if
    validTime = valid
  end function validTime

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_generateFilename
  !> @{
  !> Generates the filename for the current variable
  character(350) function generateFilename (outputForm, descriptor)

    implicit none

    type(outputDescriptor), intent(in)    :: descriptor
    character( * ), intent(in)                :: outputForm
    character(80)                           :: suffix = ''

    select case (trim(outputForm))
    case ("pft")
      suffix = '_perpft'
    case ("tile")
      suffix = '_pertil'
    case default
      suffix = ''
    end select
    generateFilename = trim(c_switch%output_directory) // '/' // &
                       trim(descriptor%shortName) // '_' // &
                       trim(descriptor%timeFreq) // trim(suffix) // '.nc'
  end function generateFilename

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_getDescriptor
  !> @{
  !> Retrieves a variable descriptor based on a given key (e.g. shortName)
  type (outputDescriptor) function getDescriptor (key)

    implicit none

    character(len =* ), intent(in)       :: key
    integer                            :: i

    do i = 1, descriptorCount
      if (outputDescriptors(i)%shortName == key) then
        getDescriptor = outputDescriptors(i)
        return
      end if
    end do
    print * , "Something went awry with the getDescriptor function"
  end function getDescriptor

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_getIdByKey
  !> @{
  !> Finds the id of the variable with the following key
  integer function getIdByKey (key)

    implicit none
    character( * ), intent(in)   :: key
    integer :: i
    do i = 1, variableCount
      if (netcdfVars(i)%key == key) then
        getIdByKey = i
        return
      end if
    end do
    getIdByKey = 0
  end function getIdByKey

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_createNetCDF
  !> @{
  !> Creates the output netcdf files
  subroutine createNetCDF (fileName, id, outputForm, descriptor, timeFreq, units, nameInCode)

    use fileIOModule
    use ctemStateVars,     only : c_switch
    use classicParams,        only : ignd,icc,nmos,iccp1

    implicit none

    character( * ), intent(in)              :: fileName
    character( * ), intent(inout)           :: outputForm
    character( * ), intent(in)              :: timeFreq, units
    type(outputDescriptor), intent(in)    :: descriptor
    integer, intent(in)                   :: id
    character( * ), intent(in)              :: nameInCode

    character(8)  :: today
    character(10) :: now
    character(80) :: format_string
    integer                     :: ncid, varid, suffix,i,timeLength
    integer                     :: DimId,lonDimId,latDimId,tileDimId,pftDimId,layerDimId,timeDimId
    real, dimension(2)          :: xrange, yrange
    integer, dimension(:), allocatable :: intArray

    logical, pointer :: projectedGrid
    integer :: dummyId(0)
    character(:), pointer :: Comment   !< Comment about the run that will be written to the output netcdfs
    logical, pointer :: leap           !< set to true if all/some leap years in the .MET file have data for 366 days
    !< also accounts for leap years in .MET when cycling over meteorology (metLoop > 1)

    ! Point pointers
    projectedGrid => c_switch%projectedGrid
    leap => c_switch%leap
    Comment => c_switch%Comment

    ncid = netcdfVars(variableCount)%ncid

    ! ---------
    call ncPutAtt(ncid, nf90_global, 'title',charvalues = 'CLASSIC output file')
    call date_and_time(today,now)

    call ncPutAtt(ncid,nf90_global,'timestamp',charvalues = today//' '//now(1:4))
    call ncPutAtt(ncid,nf90_global,'Conventions',charvalues = 'COARDS') ! FLAG remove?
    call ncPutAtt(ncid,nf90_global,'node_offset',intvalues = 1)

    !----1 - Longitude and Latitude

    lonDimId = ncDefDim(ncid,'lon',myDomain%cntx)
    latDimId = ncDefDim(ncid,'lat',myDomain%cnty)

    if (.not. projectedGrid) then
      varid = ncDefVar(ncid,'lon',nf90_double,[lonDimId])
    else
      varid = ncDefVar(ncid,'lon',nf90_double,[lonDimId, latDimId])
    end if
    call ncPutAtt(ncid,varid,'standard_name',charvalues = 'Longitude')
    call ncPutAtt(ncid,varid,'long_name',charvalues = 'longitude')
    call ncPutAtt(ncid,varid,'units',charvalues = 'degrees_east')
    ! call ncPutAtt(ncid,varid,'actual_range',xrange) #FLAG need to find the xrange from all_lon.

    if (.not. projectedGrid) then
      varid = ncDefVar(ncid,'lat',nf90_double,[latDimId])
    else
      varid = ncDefVar(ncid,'lat',nf90_double,[londimId, latDimId])
    end if
    call ncPutAtt(ncid,varid,'long_name',charvalues = 'latitude')
    call ncPutAtt(ncid,varid,'standard_name',charvalues = 'Latitude')
    call ncPutAtt(ncid,varid,'units',charvalues = 'degrees_north')
    ! call ncPutAtt(ncid,varid,'actual_range',yrange) #FLAG need to find the xrange from all_lon.

    ! The landCoverFrac, while 'grid' so that it will output in all runs
    ! needs to be switched to 'pft' so that it can be properly given the
    ! dimensions needed.
    if (trim(descriptor%shortName) == 'landCoverFrac') outputForm = 'pft'
    if (trim(descriptor%shortName) == 'landCoverExist') outputForm = 'pft'

    select case (trim(outputForm))

    case ("tile")       ! Per tile outputs

      tileDimId = ncDefDim(ncid,'tile',nmos)
      varid = ncDefVar(ncid,'tile',nf90_short,[tileDimId])
      call ncPutAtt(ncid,varid,'long_name',charvalues = 'tile')
      call ncPutAtt(ncid,varid,'units',charvalues = 'tile number')

    case ("pft")         ! Per PFT outputs

      if (descriptor%includeBareGround) then
        pftDimId = ncDefDim(ncid,'pft',iccp1)
      else
        pftDimId = ncDefDim(ncid,'pft',icc)
      end if

      varid = ncDefVar(ncid,'pft',nf90_short,[pftDimId])
      call ncPutAtt(ncid,varid,'long_name',charvalues = 'Plant Functional Type')
      call ncPutAtt(ncid,varid,'units',charvalues = 'PFT')

    case ("layer")      ! Per layer outputs

      layerDimId = ncDefDim(ncid,'layer',ignd)
      varid = ncDefVar(ncid,'layer',nf90_short,[layerDimId])
      call ncPutAtt(ncid,varid,'long_name',charvalues = 'soil layer')
      call ncPutAtt(ncid,varid,'units',charvalues = 'layer')

    end select

    ! Figure out the total run length, make a time vector and add to file.
    call determineTime(timeFreq)

    timeLength = size(timeVect)

    ! Set up the time dimension
    timeDimId = ncDefDim(ncid,'time',timeLength)

    varid = ncDefVar(ncid,'time',nf90_double,[timeDimId])

    call ncPutAtt(ncid,varid,'long_name',charvalues = 'time')

    ! timestart = "days since "//str(refyr)//"-01-01 00:00"
    format_string = "(A11,I4,A12)"
    write (timestart,format_string) "days since ",refyr,"-01-01 00:00"
    call ncPutAtt(ncid,varid,'units',charvalues = trim(timestart))

    if (leap) then
      call ncPutAtt(ncid,varid,'calendar',charvalues = "standard")
    else
      call ncPutAtt(ncid,varid,'calendar',charvalues = "365_day")
    end if

    call ncEndDef(ncid)

    call ncPutDimValues(ncid, 'time', realValues = timeVect, count = (/timelength/))

    deallocate(timeVect) ! needs to be deallocated so the next file can allocate it.

    ! Fill in the dimension variables and define the model output vars
    if (.not. projectedGrid) then
      call ncPutDimValues(ncid, 'lon', realValues = myDomain%lonUnique, count = (/myDomain%cntx/))
      call ncPutDimValues(ncid, 'lat', realValues = myDomain%latUnique, count = (/myDomain%cnty/))
    else ! projected grid
      ! Since these are flattened arrays we will use the ncPutVar function which will reshape them
      ! call ncPutVar(ncid, 'lon', realValues = myDomain%allLonValues,start = [1, 1], count = [myDomain%cntx,myDomain%cnty])
      ! call ncPutVar(ncid, 'lat', realValues = myDomain%allLatValues,start = [1, 1], count = [myDomain%cntx,myDomain%cnty])
      ! Use arrays that store the lats/lons for the part of the domain being processed.
      call ncPutVar(ncid, 'lon', realValues = myDomain%lonUnique,start = [1,1], count = [myDomain%cntx,myDomain%cnty])
      call ncPutVar(ncid, 'lat', realValues = myDomain%latUnique,start = [1,1], count = [myDomain%cntx,myDomain%cnty])
    end if

    select case (trim(outputForm))
    case ("tile")       ! Per tile outputs
      allocate(intArray(nmos))
      intArray = identityVector(nmos)
      call ncPutDimValues(ncid, 'tile', intValues = intArray, count = (/nmos/))
      call ncReDef(ncid)
      varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_double, [lonDimId,latDimId,tileDimId,timeDimId])

    case ("pft")         ! Per PFT outputs

      if (descriptor%includeBareGround) then
        allocate(intArray(iccp1))
        intArray = identityVector(iccp1)
        call ncPutDimValues(ncid, 'pft', intValues = intArray, count = (/iccp1/))
      else
        allocate(intArray(icc))
        intArray = identityVector(icc)
        call ncPutDimValues(ncid, 'pft', intValues = intArray, count = (/icc/))
      end if


      call ncReDef(ncid)

      if (descriptor%includeBareGround) then
        ! do something for cells that have bare ground
        suffix = suffix + 1
      else
        ! do something for cells that don't have bare ground
        suffix = suffix
      end if
      varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_double, [lonDimId,latDimId,pftDimId,timeDimId])

    case ("layer")      ! Per layer outputs

      allocate(intArray(ignd))
      intArray = identityVector(ignd)
      call ncPutDimValues(ncid, 'layer',intValues = intArray, count = (/ignd/))
      call ncReDef(ncid)
      varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_double, [lonDimId,latDimId,layerDimId,timeDimId])

    case default        ! Grid average outputs

      call ncReDef(ncid)
      varid = ncDefVar(ncid, trim(descriptor%shortName), nf90_double, [lonDimId,latDimId,timeDimId])

    end select

    call ncPutAtt(ncid,varid,'long_name',charvalues = descriptor%longName)
    call ncPutAtt(ncid,varid,'units',charvalues = units)
!    call ncPutAtt(ncid,varid,'_FillValue',realvalues = fill_value)
    call ncPutAtt(ncid,varid,'nameInCode',charvalues = nameInCode)

    if (projectedGrid) then
      ! Identify lat/lon as the coordinate variables (used by cdo) for rotated lat-lon grids.
      call ncPutAtt(ncid,varid,'coordinates',charvalues = 'lat lon')
      ! Add variable with coordinates of rotated pole specified via attributes.
      ! The dimension ID is not optional for the ncDefVar function in fileIOModule.F90, so use a dummy array of size 0.
      varid = ncDefVar(ncid, 'rotated_pole',nf90_byte,dummyId)
      call ncPutAtt(ncid,varid,'grid_mapping_name',charvalues = 'rotated_latitude_longitude')
      call ncPutAtt(ncid,varid,'grid_north_pole_latitude',realvalues = 42.5)
      call ncPutAtt(ncid,varid,'grid_north_pole_longitude',realvalues = 83.)
    end if

    call ncPutAtt(ncid,nf90_global,'Comment',c_switch%Comment)

    call ncEndDef(ncid)

  end subroutine createNetCDF

  !! @}
  !---------------------------------------------------------------------------------------
  !> \ingroup outputmanager_determineTime
  !> @{
  !> Determines the time vector for this run. This implictly
  !! assumes that leap year meteorological forcing is used for runs with metLoop = 1, otherwise
  !! the timing of the leap years will be off in the output files.
  subroutine determineTime (timeFreq)

    use fileIOModule
    use ctemStateVars,     only : c_switch
    use generalUtils,       only : findLeapYears
    use classicParams,        only : monthend,DELT

    implicit none

    character( * ), intent(in)              :: timeFreq
    real, allocatable, dimension(:) :: temptime
    integer :: totsteps, totyrs, i, st, en, j, m, cnt, totyrs2
    integer :: lastDOY, length, styr, endyr,numsteps
    logical :: leapnow

    logical, pointer :: leap           !< set to true if all/some leap years in the .MET file have data for 366 days
    !< also accounts for leap years in .MET when cycling over meteorology (metLoop > 1)
    integer, pointer :: metLoop        !< no. of times the meteorological data is to be looped over. this
    !< option is useful to equilibrate CTEM's C pools
    integer, pointer :: readMetStartYear !< First year of meteorological forcing to read in from the met file
    integer, pointer :: readMetEndYear   !< Last year of meteorological forcing to read in from the met file
    integer, pointer :: jhhstd  !< day of the year to start writing the half-hourly output
    integer, pointer :: jhhendd !< day of the year to stop writing the half-hourly output
    integer, pointer :: jdstd   !< day of the year to start writing the daily output
    integer, pointer :: jdendd  !< day of the year to stop writing the daily output
    integer, pointer :: jhhsty  !< simulation year (iyear) to start writing the half-hourly output
    integer, pointer :: jhhendy !< simulation year (iyear) to stop writing the half-hourly output
    integer, pointer :: jdsty   !< simulation year (iyear) to start writing the daily output
    integer, pointer :: jdendy  !< simulation year (iyear) to stop writing the daily output
    integer, pointer :: jmosty    !< Year to start writing out the monthly output files. If you want to write monthly outputs right

    leap            => c_switch%leap
    metLoop         => c_switch%metLoop
    readMetStartYear=> c_switch%readMetStartYear
    readMetEndYear  => c_switch%readMetEndYear
    jhhstd          => c_switch%jhhstd
    jhhendd         => c_switch%jhhendd
    jdstd           => c_switch%jdstd
    jdendd          => c_switch%jdendd
    jhhsty          => c_switch%jhhsty
    jhhendy         => c_switch%jhhendy
    jdsty           => c_switch%jdsty
    jdendy          => c_switch%jdendy
    jmosty          => c_switch%jmosty

    lastDOY = 365

    !  refyr is set to be the start year of the run.
    refyr = readMetStartYear
    consecDays = 0.

    select case (trim(timeFreq))

    case ("annually")
      totyrs = (readMetEndYear - readMetStartYear + 1) * metLoop
      totsteps = totyrs
      allocate(timeVect(totsteps))
      do i = 1, totsteps
        if (leap) call findLeapYears(readMetStartYear + i - 1,leapnow,lastDOY)
        timeVect(i) = (readMetStartYear + i - 1 - refyr) * lastDOY
      end do

    case ("monthly")
      ! Monthly may start writing later (after jmosty) so make sure to account for that.
      ! Also if leap years are on, it changes the timestamps
      totyrs = (readMetEndYear - readMetStartYear + 1) * metLoop
      styr = readMetStartYear

      if (readMetStartYear < jmosty .and. jmosty <= readMetStartYear + totyrs - 1) then
        totyrs = readMetStartYear + totyrs - jmosty
        styr = jmosty
      else if (jmosty > readMetStartYear + totyrs - 1) then
        print * ,'determineTime says: Warning - jmosty is set to a year beyond the end of the run'
        print * ,'monthly files will still be made but will include all years'
      end if
      totsteps = totyrs * 12
      allocate(timeVect(totsteps))
      do i = 1, totyrs
        ! Find out if this year is a leap year. It adjusts the monthend array.
        if (leap) call findLeapYears(styr + i - 1,leapnow,lastDOY)
        do m = 1, 12
          j = ((i - 1) * 12) + m
          timeVect(j) = (styr + i - 1 - refyr) * lastDOY + monthend(m + 1) - 1
        end do
      end do

    case ("daily")
      ! Daily may start writing later (after jdsty) and end earlier (jdendy) so make sure to account for that.
      ! Also likely doesn't do all days of the year. Lastly if leap years are on, it changes the timestamps

      ! First guess for total years
      totyrs = (readMetEndYear - readMetStartYear + 1) * metLoop

      ! Sanity check on jdsty and jdendy
      if ((readMetStartYear + totyrs - 1) < jdsty .or. readMetStartYear > jdendy) then
        print * ,' ** addTime says: Check your daily output file start and end points, they are outside the range of this run'
        stop
      end if

      ! Take the possible number of years then trim based on jdsty and jdendy
      styr = readMetStartYear

      endyr = readMetStartYear + totyrs - 1

      if ( jdsty > endyr .or. jdendy > endyr  .or. jdsty < readMetStartYear .or. jdendy < readMetStartYear) then
        write( * , * )'It seems you are requesting daily output outside the simulated years of'
        write( * , * )styr,' to ',endyr
        print * ,'Requested was:'
        print * ,jdsty,' to ',jdendy
        stop
      end if

      totyrs2 = jdendy - jdsty + 1

      ! Now determine the total number of timesteps (days) across all years
      totsteps = 0
      allocate(timeVect(0))
      cnt = 0
      ! Create the time vector to write to the file
      !               do i = 1, totyrs2
      do i = jdsty + 1, totyrs2 + jdsty
        if (leap) call findLeapYears(styr + i - 1,leapnow,lastDOY)
        st = max(1, jdstd) - 1  ! minus 1 since we are referencing the 01-01 day of refyr
        en = min(jdendd, lastDOY) - 1 ! minus 1 since we are referencing the 01-01 day of refyr
        totsteps = totsteps + (en - st + 1)
        allocate(temptime(totsteps))
        length = size(timeVect)
        if (i > 1) then
          temptime(1 : length) = timeVect
        end if
        do j = st,en
          cnt = cnt + 1
          temptime(cnt) = (i - 1 - refyr) * lastDOY + j
        end do
        call move_alloc(temptime,timeVect)
      end do

    case ("halfhourly")
      ! Similar to daily in that it may start writing later (after jhhsty) and end earlier (jhhendy) so make sure to account for that.
      ! Also likely doesn't do all days of the year. Lastly if leap years are on, it changes the timestamps

      ! First guess for total years
      totyrs = (readMetEndYear - readMetStartYear + 1) * metLoop

      ! Sanity check on jhhsty and jhhendy
      if ((readMetStartYear + totyrs - 1) < jhhsty .or. readMetStartYear > jhhendy) then
        print * ,' ** addTime says: Check your half - hourly output file start and end points, they are outside the range of this run'
        stop
      end if

      ! Take the possible number of years then trim based on jdsty and jdendy
      styr = readMetStartYear

      endyr = readMetStartYear + totyrs - 1

      if ( jhhsty > endyr .or. jhhendy > endyr  .or. jhhsty < readMetStartYear .or. jhhendy < readMetStartYear) then
        write( * , * )'It seems you are requesting half - hourly output outside the simulated years of'
        write( * , * )styr,' to ',endyr
        print * ,'Requested was:'
        print * ,jhhsty,' to ',jhhendy
        stop
      end if

      totyrs2 = jhhendy - jhhsty + 1

      ! Now determine the total number of timesteps (half-hours) across all years
      totsteps = 0
      allocate(timeVect(0))
      cnt = 0
      numsteps = int(86400./DELT)
      do i = jhhsty + 1, totyrs2 + jhhsty
        if (leap) call findLeapYears(styr + i - 1,leapnow,lastDOY)
        st = max(1, jhhstd) - 1  ! minus 1 since we are referencing the 01-01 day of refyr
        en = min(jhhendd, lastDOY) - 1 ! minus 1 since we are referencing the 01-01 day of refyr
        totsteps = totsteps + (en - st + 1) * numsteps
        allocate(temptime(totsteps))
        length = size(timeVect)
        if (i > 1) then
          temptime(1 : length) = timeVect
        end if
        do j = st,en
          do m = 1,numsteps
            cnt = cnt + 1
            temptime(cnt) = (i - 1 - refyr) * lastDOY + j + (m - 1) / real(numsteps)
          end do
        end do
        call move_alloc(temptime,timeVect)
      end do

    case default
      print * ,'addTime says - Unknown timeFreq: ',timeFreq
      stop
    end select

  end subroutine determineTime

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_writeOutput1D
  !> @{
  !> Write model outputs to already created netcdf files
  subroutine writeOutput1D (lonLocalIndex,latLocalIndex,key,timeStamp,label,data,specStart)

    use fileIOModule



    implicit none

    integer, intent(in) :: lonLocalIndex,latLocalIndex
    character( * ), intent(in) :: key
    integer :: ncid, timeIndex, id, length
    real, dimension(:), intent(in) :: timeStamp
    real, dimension(:), intent(in)  :: data
    character( * ), intent(in)        :: label
    integer,  optional :: specStart
    integer :: start
    integer :: posTimeWanted
    real, dimension(:), allocatable :: localData
    real, dimension(1) :: localStamp
    real, allocatable, dimension(:) :: timeWritten

    allocate(localData(size(data)))
    localStamp = timeStamp
    localData = data

    ! print*,key,timeStamp,label,lonLocalIndex,latLocalIndex
    id = getIdByKey(key)

    if (id == 0) then
      !print * ,'writeOutput1D says: Your requested key does not exist (' // trim(key) // ') in netcdfVars.'
      !print * , 'Possible reasons include '// trim(key) // ' not in xml file so no netcdf created'
      !print * , 'or mismatch between xml group and model switch for this key. Model run will continue'
      !print * , 'without writing this variable.'
      return
    end if

    ncid = netcdfVars(id)%ncid

    length = size(data)
    if (present(specStart)) then
      start = specStart
    else
      start = 1
    end if

    ! Check if the time period has already been added to the file
    timeIndex = ncGetDimLen(ncid, "time")

    allocate(timeWritten(timeIndex))
    timeWritten = ncGetDimValues(ncid, "time", count = (/timeIndex/))
    posTimeWanted = checkForTime(timeIndex,timeWritten,localStamp(1))

    if (posTimeWanted == 0) then
      print * ,'missing timestep in output file ',key,localStamp
      print * ,' ** Did you set leap = true but give CLASSIC inputs that do not have leap years?'
      print * ,'OR did you try and write out a variable that is incorrectly specified in your output xml file?'
      stop
    else
      timeIndex = posTimeWanted
    end if

    if (length > 1) then
      call ncPutVar(ncid, label, localData, start = [lonLocalIndex,latLocalIndex,start,timeIndex], count = [1,1,length,1])
    else
      call ncPutVar(ncid, label, localData, start = [lonLocalIndex,latLocalIndex,timeIndex], count = [1,1,1])
    end if

  end subroutine writeOutput1D

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_checkForTime
  !> @{
  !> Find if a time period is already in the timeIndex of the file
  integer function checkForTime (timeIndex,timeWritten,timeStamp)

    implicit none

    real, intent(in)   :: timeStamp
    integer, intent(in) :: timeIndex
    real, dimension(:), intent(in) :: timeWritten
    integer :: i
    do i = 1, timeIndex
      if (timeWritten(i) == timeStamp) then
        checkForTime = i
        return
      end if
    end do
    checkForTime = 0
  end function checkForTime

  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_closeNCFiles
  !> @{
  !> Close all output netcdfs or just a select file
  subroutine closeNCFiles (incid)

    use fileIOModule

    implicit none

    integer, optional   :: incid
    integer :: i

    if (present(incid)) then
      call ncClose(incid)
    else
      do i = 1, variableCount
        call ncClose(netcdfVars(i)%ncid)
      end do
    end if

  end subroutine closeNCFiles
  !! @}
  !---------------------------------------------------------------------------------------

  !> \ingroup outputmanager_identityVector
  !> @{
  pure function identityVector (n) result(res)
    implicit none
    integer, allocatable :: res(:)
    integer, intent(in) :: n
    integer             :: i
    allocate(res(n))
    forall (i = 1:n)
      res(i) = i
    end forall
  end function identityVector

  !! @}
  !---------------------------------------------------------------------------------------
  !> \namespace outputmanager
  !! Central module for all netCDF output file operations
end module outputManager
