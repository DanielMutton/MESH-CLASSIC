program metASCIILoader
    use fileIOModule
    implicit none
    real, allocatable   :: shortWave(:), longWave(:), precipitation(:), &
        temperature(:), humidity(:), wind(:), pressure(:), time(:)
    integer, allocatable    :: hour(:), minute(:), day(:), year(:)
    character(300)          :: inputFile, charLon, charLat,charContainsLeaps
    logical                 :: containsLeaps
    integer                 :: timesteps, i, currentYear
    integer, dimension(12)  :: daysInMonth = [ 31,28,31,30,31,30,31,31,30,31,30,31 ]
    real, parameter         :: fillValue = 1.E38
    real                    :: lon, lat

    call processArguments
    open(unit = 10, file = inputFile, form = 'formatted', status = 'old', action = 'read')
    call initializeData
    call loadData
    call processTimesteps
    call exportData
    close(unit = 10)
contains
    subroutine processArguments
        if (iargc() .ne. 4) then
            print*,'Usage is: metASCIILoader [input file] [lon] [lat] [containsLeaps]'
            print*,'containsLeaps should be true if the file contains leap years, otherwise false'
            stop
        endif

        call getarg(1, inputFile)
        call getarg(2, charLon)
        call getarg(3, charLat)
        call getarg(4, charContainsLeaps)
        if (trim(charContainsLeaps) == 'true') then
          containsLeaps = .true.
        else
          containsLeaps = .false.
        endif
        lon = charToReal(charLon)
        lat = charToReal(charLat)
    end subroutine processArguments

    subroutine initializeData
        integer                 :: fileSize

        inquire(file = inputFile, size = fileSize)
        timesteps = (fileSize + 1) / 91;
        allocate(hour(timesteps))
        allocate(minute(timesteps))
        allocate(day(timesteps))
        allocate(year(timesteps))
        allocate(shortWave(timesteps))
        allocate(longWave(timesteps))
        allocate(precipitation(timesteps))
        allocate(temperature(timesteps))
        allocate(humidity(timesteps))
        allocate(wind(timesteps))
        allocate(pressure(timesteps))
        allocate(time(timesteps))
    end subroutine initializeData

    subroutine loadData()
        integer                 :: i
        character(*), parameter :: format = '(1X, I2, I3, I5, I6, 2F9.2, E14.4, F9.2, E12.3, F8.2, F12.2, 3F9.2, F9.4)'

        do i = 1, timesteps
            read(unit = 10, fmt = format) hour(i), minute(i), day(i), year(i), &
                shortWave(i), longWave(i), precipitation(i), temperature(i), humidity(i), wind(i), pressure(i)
        enddo
    end subroutine loadData

    logical function isLeapYear(thisYear)
        integer, intent(in) :: thisYear
        if (mod(thisYear,4).ne.0) then
            isLeapYear = .false.
        else if (mod(thisYear,100).ne.0) then
            isLeapYear = .true.
        else if (mod(thisYear,400).ne.0) then
            isLeapYear = .false.
        else
            isLeapYear = .true.
        end if
    end function isLeapYear

    real function buildTimestep(hour, minute, dayIn, year)
        integer, intent(in) :: hour, minute, dayIn, year
        real                :: time, day
        integer             :: i, month

        day = dayIn
        time = (real(hour) * 60 + real(minute)) / 1440  ! 1440 minutes in a day

        month = 0
        do i = 1, 12
            if (day > daysInMonth(i)) then
                day = day - daysInMonth(i)
            else
                month = i
                exit
            endif
        enddo

        buildTimestep = year * 10000 + month * 100 + day + time
    end function buildTimestep

    subroutine processTimesteps
        currentYear = 0
        do i = 1, timesteps
            if (currentYear /= year(i)) then
                currentYear = year(i)
                if (containsLeaps) then
                  if (isLeapYear(year(i))) then
                    daysInMonth(2) = 29
                  else
                    daysInMonth(2) = 28
                  endif
                endif
            endif
            time(i) = buildTimestep(hour(i), minute(i), day(i), year(i))
        enddo
    end subroutine processTimesteps

    real function charToReal(input)
        character(len=*), intent(in)    :: input    !< Char input
        read(input,*) charToReal
    end function charToReal

    subroutine exportData
        call exportVariable('sw', shortWave)
        call exportVariable('lw', longWave)
        call exportVariable('pr', precipitation)
        call exportVariable('ta', temperature)
        call exportVariable('qa', humidity)
        call exportVariable('wi', wind)
        call exportVariable('ap', pressure)
    end subroutine exportData

    subroutine exportVariable(label, variable)
        character(*), intent(in)    :: label
        real, intent(inout)         :: variable(:)
        integer                     :: fileId
        character(200)              :: filename, units, gridType, title, name
        integer                     :: varId, timeDimId, lonDimId, latDimId

        ! Make the filename
        filename = 'metVar_' // label // '.nc'

        ! If the file doesn't already exist, then initialize the file
        if (.not.fileExists(filename)) then
            ! Create the variable file
            fileId = ncCreate(filename, NF90_CLOBBER)

            ! Populate the variable with attribute content
            select case(label)
            case('sw')
                units = 'J/m2'
                gridType = 'gaussian'
                title = 'Incoming short wave radiation'
                name = 'short wave'
            case('ap')
                units = 'Pa'
                gridType = 'gaussian'
                title = 'Atmospheric pressure'
                name = 'athmospheric pressure'
            case('pr')
                units = 'mm/hh'
                gridType = 'gaussian'
                title = 'Total precipitation'
                name = 'precipitation'
            case('qa')
                units = 'g/g'
                gridType = 'gaussian'
                title = 'Air specific humidity'
                name = 'humidity'
            case('ta')
                units = 'K'
                gridType = 'gaussian'
                title = 'Temperature'
                name = 'temperature'
            case('wi')
                units = 'm/s'
                gridType = 'gaussian'
                title = 'Wind'
                name = 'wind'
            case('lw')
                units = ''
                gridType = 'gaussian'
                title = 'Incoming long wave radiation'
                name = 'long wave'
            case default
                stop ('Unrecognized label')
            end select

            ! Define file attributes
            call ncPutAtt(fileId, nf90_global, 'title', charValues = title)
            call ncPutAtt(fileId, nf90_global, 'units', charValues = units)
            call ncPutAtt(fileId, nf90_global, 'grid_type', charValues = gridType)
            call ncPutAtt(fileId, nf90_global, '_FillValue', realValues = fillValue)

            ! Define the longitude dimension
            lonDimId = ncDefDim(fileId, 'lon', 1)
            varid = ncDefVar(fileId, 'lon', nf90_double, [lonDimId])
            call ncPutAtt(fileId, varId, 'standard_name', charValues = 'Longitude')
            call ncPutAtt(fileId, varId, 'units', charValues = 'degrees_east')
            call ncPutAtt(fileId, varId, 'axis', charValues = 'X')

            call ncEndDef(fileId)

            call ncPutDimValues(fileId, 'lon', realValues = [lon], count = [1])

            call ncRedef(fileId)

            ! Define the latitude dimension
            latDimId = ncDefDim(fileId, 'lat', 1)
            varid = ncDefVar(fileId, 'lat', nf90_double, [latDimId])
            call ncPutAtt(fileId, varId, 'standard_name', charValues = 'Latitude')
            call ncPutAtt(fileId, varId, 'units', charValues = 'degrees_north')
            call ncPutAtt(fileId, varId, 'axis', charValues = 'Y')

            call ncEndDef(fileId)

            call ncPutDimValues(fileId, 'lat', [lat], count = [1])

            call ncRedef(fileId)

            ! Define the time dimension
            timeDimId = ncDefDim(fileId, 'time', size(time))
            varid = ncDefVar(fileId, 'time', nf90_double, [timeDimId])
            call ncPutAtt(fileId, varId, 'standard_name', charValues = 'time')
            call ncPutAtt(fileId, varId, 'units', charValues = 'day as YYYYMMDD.FFFF')
            call ncPutAtt(fileId, varId, 'calendar', charValues = 'proleptic_gregorian')
            call ncEndDef(fileId)
            call ncPutDimValues(fileId, 'time', time, count = [size(time)])

            call ncRedef(fileId)

            ! Define the variable
            varid = ncDefVar(fileId, label, nf90_double, [lonDimId, latDimId, timeDimId])
            call ncPutAtt(fileId, varId, '_FillValue', realValues = fillValue)
            call ncEndDef(fileId)

        else
            fileId = ncOpen(filename, nf90_write)
        endif

        ! Put in data
        call ncPutVar(fileId, label, realValues = variable, start = [1, 1, 1], count = [1, 1, size(variable)])

        call ncClose(fileId)
    end subroutine exportVariable

    logical function fileExists(filename)
        character(*), intent(in) :: filename
        inquire(file = filename, exist = fileExists)
    end function fileExists

end program
