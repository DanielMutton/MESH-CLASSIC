!> Converts old format INI and CTM files into netcdf format initialization files.
!! Also takes in namelist format. Useful for tests at the site-level.

program initFileConverter

    use fileIOModule

    implicit none

    real                    :: lon, lat
    integer                 :: NLTEST, NMTEST
    character(300)          :: INIFile, CTMFile,charLon, charLat
    character(3)            :: fileType
    integer                 :: fileId,slopeDimId,iccp2DimId,icctemDimId,icDimId,layerDimId, & !monthsDimId,
                                icp1DimId,lonDimId,latDimId,tileDimId
    character(50)           :: units
    character(350)          :: long_name
    integer, allocatable    :: dimArray(:),start(:),count(:)
    real :: dummy       ! dummyvar

    logical                 :: inclCTEM = .false.   ! True if a .CTM file was provided to the program
    logical                 :: fromnml = .false.    ! True if we read in from a namelist file
    integer                 :: ignd = 3             ! Number of soil layers in INI file.
    integer                 :: ican = 4
    integer                 :: icc = 9
    real, parameter         :: fillValue = -999.
    real                    :: grclarea = 100.     !<area of the grid cell, \f$km^2\f$, kind of meaningless at the point scale but needed for run if CTEM and fire on.

    ! All variable descriptions are listed in exportData
    real                    :: ZRFHROW
    real                    :: ZRFMROW
    real                    :: ZBLDROW
    real                    :: GCROW
    real, allocatable, dimension(:,:):: FCANROT
    real, allocatable, dimension(:)  :: FAREROT
    real, allocatable, dimension(:,:):: RSMNROT
    real, allocatable, dimension(:,:):: QA50ROT
    real, allocatable, dimension(:,:):: VPDAROT
    real, allocatable, dimension(:,:):: VPDBROT
    real, allocatable, dimension(:,:):: PSGAROT
    real, allocatable, dimension(:,:):: PSGBROT
    real, allocatable, dimension(:,:):: ALVCROT
    real, allocatable, dimension(:,:):: ALICROT
    real, allocatable, dimension(:,:):: PAMNROT
    real, allocatable, dimension(:,:):: PAMXROT
    real, allocatable, dimension(:,:):: LNZ0ROT
    real, allocatable, dimension(:,:):: CMASROT
    real, allocatable, dimension(:,:):: ROOTROT
    real, allocatable, dimension(:)  :: DRNROT
    real, allocatable, dimension(:)  :: SDEPROT
    real, allocatable, dimension(:)  :: XSLPROT
    real, allocatable, dimension(:)  :: GRKFROT
    real, allocatable, dimension(:)  :: WFSFROT
    real, allocatable, dimension(:)  :: WFCIROT
    integer, allocatable, dimension(:)  :: MIDROT
    real, allocatable, dimension(:,:):: SANDROT
    real, allocatable, dimension(:,:):: CLAYROT
    real, allocatable, dimension(:,:):: ORGMROT
    real, allocatable, dimension(:,:):: TBARROT
    real, allocatable, dimension(:,:):: THLQROT
    real, allocatable, dimension(:,:):: THICROT
    real, allocatable, dimension(:)   :: DELZ
    real, allocatable, dimension(:,:)   :: ZBOT
    real, allocatable, dimension(:)  :: TCANROT
    real, allocatable, dimension(:)  :: TSNOROT
    real, allocatable, dimension(:)  :: TPNDROT
    real, allocatable, dimension(:)  :: ZPNDROT
    real, allocatable, dimension(:)  :: RCANROT
    real, allocatable, dimension(:)  :: SCANROT
    real, allocatable, dimension(:)  :: SNOROT
    real, allocatable, dimension(:)  :: ALBSROT
    real, allocatable, dimension(:)  :: RHOSROT
    real, allocatable, dimension(:)  :: GROROT
    integer, allocatable, dimension(:)  :: SOCIROT

    !real :: extnprob
    !real :: prbfrhuc
    real, allocatable, dimension(:,:):: ailcminrow
    real, allocatable, dimension(:,:):: ailcmaxrow
    real, allocatable, dimension(:,:):: dvdfcanrow
    real, allocatable, dimension(:,:):: fcancmxrow
    real, allocatable, dimension(:,:):: gleafmasrow
    real, allocatable, dimension(:,:):: bleafmasrow
    real, allocatable, dimension(:,:):: stemmassrow
    real, allocatable, dimension(:,:):: rootmassrow
    real, allocatable, dimension(:,:):: grwtheffrow
    ! the tracer variables:
    ! real, allocatable, dimension(:,:):: tracerGLeafMass
    ! real, allocatable, dimension(:,:):: tracerBLeafMass
    ! real, allocatable, dimension(:,:):: tracerStemMass
    ! real, allocatable, dimension(:,:):: tracerRootMass
    ! real, allocatable, dimension(:,:,:):: tracerLitrMass
    ! real, allocatable, dimension(:,:,:):: tracerSoilCMass
    ! real, allocatable, dimension(:):: tracerMossCMass
    ! real, allocatable, dimension(:):: tracerMossLitrMass
    
    real, allocatable, dimension(:,:):: pstemmassrow
    real, allocatable, dimension(:,:):: pgleafmassrow
    real :: twarmm ! These bioclimatic vars are per grid cell so don't require allocation.
    real :: tcoldm
    real :: gdd5
    real :: aridity
    real :: srplsmon
    real :: defctmon
    real :: anndefct
    real :: annsrpls
    real :: annpcp
    real :: dry_season_length
    real, allocatable, dimension(:,:):: litrmassrow
    real, allocatable, dimension(:,:):: soilcmasrow
    !real, allocatable, dimension(:,:,:):: litrmassrow
    !real, allocatable, dimension(:,:,:):: soilcmasrow
    integer, allocatable, dimension(:,:):: lfstatusrow
    integer, allocatable, dimension(:,:):: pandaysrow
    real, allocatable, dimension(:,:):: slopefrac
    integer, allocatable, dimension(:):: ipeatlandrow   !<Peatland switch: 0 = not a peatland, 1= bog, 2 = fen
    real, allocatable, dimension(:):: Cmossmas          !<C in moss biomass, \f$kg C/m^2\f$
    real, allocatable, dimension(:):: litrmsmoss        !<moss litter mass, \f$kg C/m^2\f$
    real, allocatable, dimension(:):: dmoss             !<depth of living moss (m)
    real, allocatable, dimension(:) :: maxAnnualActLyr  !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).

    !----------

    ! Parse the arguments and determine if the input file is INI format or namelist
    call processArguments

    ! Open the INI file, or the namelist file
    open(unit = 10, file = INIFile, status = 'old', action = 'read')

    ! Read the file headers (INI or nml) so we can allocate arrays
    if (fileType == 'ini') then
        read(10,*) ! Throw out first three lines.
        read(10,*)
        read(10,*)
        READ(10,5020) lat,lon,ZRFMROW,ZRFHROW,ZBLDROW,GCROW,NLTEST,NMTEST
    else ! the file is a namelist
        call readNamelistHeader
    end if

    ! Allocate the arrays
    call setupArrays

    ! Read in the rest of the input from the INI/nml file
    if (fileType == 'ini') then
        call loadINIData
    else
        call readNamelist
        fromnml = .true.
    end if

    ! Read in the CTM file
    if (inclCTEM) then
        call loadCTMData
    end if

    ! Make the file
    call makeNetCDF

    ! Write the data out to netcdf
    call exportData

    ! Close the input and output files
    close(unit = 10)
    call ncClose(fileId)
    if (inclCTEM) close(unit = 20)

    5020  FORMAT(5F10.2,F7.1,3I5)

contains

    ! ------------------------------------------------------------------------------------------------------

    subroutine processArguments

        call getarg(1, INIFile)

        if (iargc() .eq. 2) then
            call getarg(2, CTMFile)
            if (index(CTMFile,'CTM') > 0 .or. index(CTMFile,'ctm') > 0) then
                inclCTEM = .true.
            end if
        elseif (iargc() .gt. 2 .or. iargc() == 0) then
            print*
            print*,'Expecting only one or two arguments: First in the INI or namelist file'
            print*,'the second file is the CTM file and only required if you want to run with'
            print*,'CTEM on. The possible suffixes of the input files are:.INI, .ini, .CTM, .ctm, .nml, .txt'
            print*,'Output netcdf is given the same name as the INI/nml file with the suffix .nc'
        endif

        ! Parse the INIfile to help decide if it is a namelist or INI file
        if (index(INIFile,'INI') > 0 .or. index(INIFile,'ini') > 0) then
            fileType = 'ini'
        elseif (index(INIFile,'txt') > 0 .or. index(INIFile,'nml') > 0) then
            fileType = 'nml'
        else
            print*,'***Unknown file given for first file. Expected suffixes: INI,ini,txt,nml'
            call exit
        end if

    end subroutine processArguments

    ! ------------------------------------------------------------------------------------------------------

    subroutine setupArrays

        allocate(FCANROT(NMTEST,ICAN+1))
        allocate(PAMXROT(NMTEST,ICAN))
        allocate(LNZ0ROT(NMTEST,ICAN+1))
        allocate(PAMNROT(NMTEST,ICAN))
        allocate(ALVCROT(NMTEST,ICAN+1))
        allocate(CMASROT(NMTEST,ICAN))
        allocate(ALICROT(NMTEST,ICAN+1))
        allocate(ROOTROT(NMTEST,ICAN))
        allocate(RSMNROT(NMTEST,ICAN))
        allocate(QA50ROT(NMTEST,ICAN))
        allocate(VPDAROT(NMTEST,ICAN))
        allocate(VPDBROT(NMTEST,ICAN))
        allocate(PSGAROT(NMTEST,ICAN))
        allocate(PSGBROT(NMTEST,ICAN))
        allocate(DRNROT(NMTEST))
        allocate(SDEPROT(NMTEST))
        allocate(FAREROT(NMTEST))
        allocate(XSLPROT(NMTEST))
        allocate(GRKFROT(NMTEST))
        allocate(WFSFROT(NMTEST))
        allocate(WFCIROT(NMTEST))
        allocate(MIDROT(NMTEST))
        allocate(SOCIROT(NMTEST))
        allocate(SANDROT(NMTEST,ignd))
        allocate(CLAYROT(NMTEST,ignd))
        allocate(ORGMROT(NMTEST,ignd))
        allocate(TBARROT(NMTEST,ignd))
        allocate(TCANROT(NMTEST))
        allocate(TSNOROT(NMTEST))
        allocate(TPNDROT(NMTEST))
        allocate(THLQROT(NMTEST,ignd))
        allocate(THICROT(NMTEST,ignd))
        allocate(ZPNDROT(NMTEST))
        allocate(RCANROT(NMTEST))
        allocate(SCANROT(NMTEST))
        allocate(SNOROT(NMTEST))
        allocate(ALBSROT(NMTEST))
        allocate(RHOSROT(NMTEST))
        allocate(GROROT(NMTEST))
        allocate(DELZ(ignd))
        allocate(ZBOT(nmtest,ignd))

        allocate(ailcminrow(NMTEST,icc))
        allocate(ailcmaxrow(NMTEST,icc))
        allocate(dvdfcanrow(NMTEST,icc))
        allocate(gleafmasrow(NMTEST,icc))
        allocate(bleafmasrow(NMTEST,icc))
        allocate(stemmassrow(NMTEST,icc))
        allocate(rootmassrow(NMTEST,icc))
        allocate(grwtheffrow(NMTEST,icc))
        allocate(litrmassrow(NMTEST,icc+2)) ! +2 due to bare ground and land use change products pools
        allocate(soilcmasrow(NMTEST,icc+2)) ! +2 due to bare ground and land use change products pools
        !allocate(litrmassrow(NMTEST,icc+2,ignd)) ! +2 due to bare ground and land use change products pools
        !allocate(soilcmasrow(NMTEST,icc+2,ignd)) ! +2 due to bare ground and land use change products pools
        allocate(lfstatusrow(NMTEST,icc))
        allocate(pandaysrow(NMTEST,icc))
        
        !allocate(tracerGLeafMass(NMTEST,icc))
        !allocate(tracerBLeafMass(NMTEST,icc))
        !allocate(tracerStemMass(NMTEST,icc))
        !allocate(tracerRootMass(NMTEST,icc))
        !allocate(tracerLitrMass(NMTEST,icc+2,ignd))
        !allocate(tracerSoilCMass(NMTEST,icc+2,ignd))
        !allocate(tracerMossCMass(NMTEST))
        !allocate(tracerMossLitrMass(NMTEST))
        
        allocate(slopefrac(NMTEST,8))
        allocate(ipeatlandrow(NMTEST))
        allocate(Cmossmas(NMTEST))
        allocate(litrmsmoss(NMTEST))
        allocate(dmoss(NMTEST))
        allocate(maxAnnualActLyr(NMTEST))

        allocate(fcancmxrow(nmtest,icc))
        !allocate(mlightng(12))

    end subroutine setupArrays

    ! ------------------------------------------------------------------------------------------------------

    subroutine loadINIData()

        integer                 :: m,j

         DO 50 M=1,NMTEST
          READ(10,5040) (FCANROT(M,J),J=1,ICAN+1),(PAMXROT(M,J),J=1,ICAN)
          READ(10,5040) (LNZ0ROT(M,J),J=1,ICAN+1),(PAMNROT(M,J),J=1,ICAN)
          READ(10,5040) (ALVCROT(M,J),J=1,ICAN+1),(CMASROT(M,J),J=1,ICAN)
          READ(10,5040) (ALICROT(M,J),J=1,ICAN+1),(ROOTROT(M,J),J=1,ICAN)
          READ(10,5030) (RSMNROT(M,J),J=1,ICAN),(QA50ROT(M,J),J=1,ICAN)
          READ(10,5030) (VPDAROT(M,J),J=1,ICAN),(VPDBROT(M,J),J=1,ICAN)
          READ(10,5030) (PSGAROT(M,J),J=1,ICAN),(PSGBROT(M,J),J=1,ICAN)
          READ(10,5040) DRNROT(M),SDEPROT(M),FAREROT(M)
          READ(10,5090) XSLPROT(M),GRKFROT(M),WFSFROT(M),WFCIROT(M),MIDROT(M),SOCIROT(M)
          ! The soil colour index is a 'relatively' new input requirement. It is possible some older INI files 
          ! would not have this information in them. If not then throw an error and stop the run.
          if (SOCIROT(M) == 0) THEN
            print*,''
            print*,'*** Error! Your INI/nml file is missing the soil colour index! See CLASSIC Manual for information.***'
            print*,'>>> Enter a soil colour index (integer) or enter -99 to abort the program:'
            read(*,*)SOCIROT(M)
            if (SOCIROT(M) == -99) then
              print*,'Aborting...'
              exit 
            end if
          end if
          READ(10,5080) (SANDROT(M,J),J=1,ignd)
          READ(10,5080) (CLAYROT(M,J),J=1,ignd)
          READ(10,5080) (ORGMROT(M,J),J=1,ignd)
          READ(10,5050) (TBARROT(M,J),J=1,ignd),TCANROT(M),TSNOROT(M),TPNDROT(M)
          READ(10,5060) (THLQROT(M,J),J=1,ignd),(THICROT(M,J),J=1,ignd),ZPNDROT(M)
          READ(10,5070) RCANROT(M),SCANROT(M),SNOROT(M),ALBSROT(M),RHOSROT(M),GROROT(M)
50    CONTINUE
     ! In CLASS 3.6.2, we include this soil info in the INI file.
      DO 25 J=1,IGND
          READ(10,*) DELZ(J),ZBOT(1,J)
 25   CONTINUE

    ! DELZ and ZBOT get put the same for all tiles.
    if (nmtest > 1) then
        do m = 2,nmtest
            zbot(m,:) = zbot(1,:)
        end do
    end if

    ! The peatland stuff should not be in any files so put in dummy vals
    slopefrac(:,:) = 0.
    ipeatlandrow(:) = 0
    Cmossmas(:) = 0.
    litrmsmoss(:) = 0.
    dmoss(:) = 0.

    ! The active layer depth won't be in any init files:
    ! so set to a very deep value.
    maxAnnualActLyr(:) = 999.
    
5040  FORMAT(9F8.3)
5030  FORMAT(4F8.3,8X,4F8.3)
5050  FORMAT(6F10.2)
5060  FORMAT(7F10.3)
5070  FORMAT(2F10.4,F10.2,F10.3,F10.4,F10.3)
5080  FORMAT(3F10.1)
5090  FORMAT(4E8.1,2I8)

    end subroutine loadINIData

    ! ------------------------------------------------------------------------------------------------------

    subroutine loadCTMData()

        integer :: m,j,maxlayer,k

        !>Read from CTEM initialization file (.CTM)
        open(unit = 11, file = CTMFile, form = 'formatted', status = 'old', action = 'read')

        ! dump header
        read(11,*)
        read(11,*)
        read(11,*)

        do m=1,nmtest
            read(11,*) (ailcminrow(m,j),j=1,icc) ! read in but not used.
            read(11,*) (ailcmaxrow(m,j),j=1,icc) ! read in but not used.
            read(11,*) (dvdfcanrow(m,j),j=1,icc)
            read(11,*) (gleafmasrow(m,j),j=1,icc)
            read(11,*) (bleafmasrow(m,j),j=1,icc)
            read(11,*) (stemmassrow(m,j),j=1,icc)
            read(11,*) (rootmassrow(m,j),j=1,icc)
            ! There are no legacy CTM files that use icc+2 (include the luc products pool) so 
            ! assume they only go to icc+1. Put a zero in the icc+2 position later.
            ! Also no CTM files will have per layer values so just put all in layer 1,
            ! read(11,*) (litrmassrow(m,j,1),j=1,icc+1) 
            ! read(11,*) (soilcmasrow(m,j,1),j=1,icc+1)
            ! litrmassrow(m,icc+2,1)=0.
            ! soilcmasrow(m,icc+2,1)=0.
            read(11,*) (litrmassrow(m,j),j=1,icc+1) 
            read(11,*) (soilcmasrow(m,j),j=1,icc+1)

            read(11,*) (lfstatusrow(m,j),j=1,icc)
            read(11,*) (pandaysrow(m,j),j=1,icc)
        
        ! Distribute the litter and soil C across the top layers since we are now per layer
        ! this is just to speed up the time to equilibrium.
        ! maxlayer = min(5,ignd-1)
        ! do k = 1,maxlayer
        !   do j = 1, icc+1
        !     litrmassrow(m,j,k) = litrmassrow(m,j,1) / real(maxlayer)
        !     soilcmasrow(m,j,k) = soilcmasrow(m,j,1) / real(maxlayer)
        !   end do 
        ! end do
        
        end do

        ! None of the tracers will be in CTM files so just set to 0.
        ! tracerGLeafMass = 0.
        ! tracerBLeafMass = 0.
        ! tracerStemMass = 0.
        ! tracerRootMass = 0.
        ! tracerLitrMass = 0.
        ! tracerSoilCMass = 0.
        ! tracerMossCMass = 0.
        ! tracerMossLitrMass = 0. 

            ! Obsolete variables.
            !read(11,*) (mlightng(j),j=1,6)  !mean monthly lightning frequency
            !read(11,*) (mlightng(j),j=7,12) !flashes/km2.year, this is spread over other tiles below
            !read(11,*) extnprob
            !read(11,*) prbfrhuc
            !read(11,*) dummy ! was stdaln but not used so dump.
            
            
            !read in the bioclimatic parameters
            ! If your file does not have these, comment out below.
            read(11,*) twarmm, tcoldm, gdd5, aridity,srplsmon
            read(11,*) defctmon, anndefct, annsrpls, annpcp, dry_season_length
            
            ! The grwtheffrow  won't be in any init files:
            ! so set to a very high value.
            grwtheffrow = 100.

    end subroutine loadCTMData

    ! ------------------------------------------------------------------------------------------------------

    subroutine readNamelistHeader

        namelist /header/ &
        lat,&
        lon,&
        nmtest,&
        ignd,&
        ican,&
        icc

        read(unit=10,nml = header)
        rewind(10)

    end subroutine readNamelistHeader

    ! ------------------------------------------------------------------------------------------------------

    subroutine readNamelist

        namelist /classicvars/ &
            !GCROW,&
            FCANROT,&
            FAREROT,&
            !RSMNROT,&
            !QA50ROT,&
            !VPDAROT,&
            !VPDBROT,&
            !PSGAROT,&
            !PSGBROT,&
            ALVCROT,&
            ALICROT,&
            PAMNROT,&
            PAMXROT,&
            LNZ0ROT,&
            CMASROT,&
            ROOTROT,&
            DRNROT,&
            SDEPROT,&
            !XSLPROT,&
            !GRKFROT,&
            !WFSFROT,&
            !WFCIROT,&
            MIDROT,&
            SANDROT,&
            CLAYROT,&
            ORGMROT,&
            TBARROT,&
            THLQROT,&
            THICROT,&
            DELZ,&
            TCANROT,&
            TSNOROT,&
            TPNDROT,&
            ZPNDROT,&
            RCANROT,&
            SCANROT,&
            SNOROT,&
            ALBSROT,&
            RHOSROT,&
            GROROT,&
            SOCIROT, &
            !ailcminrow,&
            !ailcmaxrow,&
            fcancmxrow,&
            gleafmasrow,&
            bleafmasrow,&
            stemmassrow,&
            rootmassrow,&
            litrmassrow,&
            soilcmasrow,&
            grwtheffrow, &
            lfstatusrow,&
            pandaysrow,&
            !mlightng,&
            !extnprob,&
            !prbfrhuc,&
            slopefrac,&
            ipeatlandrow,&
            Cmossmas,&
            litrmsmoss,&
            dmoss,&
            maxAnnualActLyr,&
            ! tracerGLeafMass, &
            ! tracerBLeafMass, &
            ! tracerStemMass, &
            ! tracerRootMass, &
            ! tracerLitrMass, &
            ! tracerSoilCMass, &
            ! tracerMossCMass, &
            ! tracerMossLitrMass, &
            twarmm, &
            tcoldm, &
            gdd5, &
            aridity, &
            srplsmon, &
            defctmon, &
            anndefct, &
            annsrpls, &
            annpcp, &
            dry_season_length
            
        read(unit=10,nml = classicvars)
        !write(*,nml = classicvars)

    end subroutine readNamelist

    ! ------------------------------------------------------------------------------------------------------

    subroutine makeNetCDF

        character(400) :: filename,title
        character(8)   :: date
        integer :: i,varId
        integer :: tile(nmtest)
        integer :: icp1(ican+1)
        integer :: layer(ignd)
        integer :: ic(ican)
        integer :: icctem(icc)
        integer :: iccp2(icc+2)
        !integer :: months(12)
        integer :: slope(8)
        integer :: indexend

        tile = (/(i, i=1,nmtest, 1)/)
        icp1 = (/(i, i=1,ican+1, 1)/)
        layer = (/(i, i=1,ignd, 1)/)
        ic = (/(i, i=1,ican, 1)/)
        icctem = (/(i, i=1,icc, 1)/)
        iccp2 = (/(i, i=1,icc+2, 1)/)
        !months = (/(i, i=1,12, 1)/)
        slope = (/(i, i=1,8, 1)/)

        ! Filename is going to be the INI file name with .nc as a suffix.
        if (fileType == 'ini') then
            indexend = max(index(INIFile,'INI'),index(INIFile,'ini'))
        else
            indexend = max(index(INIFile,'nml'),index(INIFile,'txt'))
        end if
        filename = INIfile(1:indexend-1)//'nc'
        fileId = ncCreate(filename, NF90_CLOBBER)

        ! Add in the metadata for the file

        title = 'CLASSIC initialization file created from: '//INIfile
        if (inclCTEM) then
            title = title//' and '//CTMFile
        end if

        call ncPutAtt(fileId, nf90_global, 'title', charValues = trim(title))
        call date_and_time(date=date)
        call ncPutAtt(fileId, nf90_global, 'creation_date', charValues = date)

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

        ! Define the tile dimension
        tileDimId = ncDefDim(fileId, 'tile', size(tile))
        varid = ncDefVar(fileId, 'tile', nf90_int, [tileDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'land surface tile')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'tile', intValues=tile, count = [size(tile)])

        call ncRedef(fileId)

        ! Define the icp1 dimension
        icp1DimId = ncDefDim(fileId, 'icp1', size(icp1))
        varid = ncDefVar(fileId, 'icp1', nf90_int, [icp1DimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'physics (CLASS) PFTs plus bareground')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'icp1', intValues=icp1, count = [size(icp1)])

        call ncRedef(fileId)

        ! Define the layer dimension
        layerDimId = ncDefDim(fileId, 'layer', size(layer))
        varid = ncDefVar(fileId, 'layer', nf90_int, [layerDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'ground layers')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'layer', intValues=layer, count = [size(layer)])

        call ncRedef(fileId)

        ! Define the ic dimension
        icDimId = ncDefDim(fileId, 'ic', size(ic))
        varid = ncDefVar(fileId, 'ic', nf90_int, [icDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'physics (CLASS) PFTs')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'ic', intValues=ic, count = [size(ic)])

        call ncRedef(fileId)

        ! Define the icctem dimension
        icctemDimId = ncDefDim(fileId, 'icctem', size(icctem))
        varid = ncDefVar(fileId, 'icctem', nf90_int, [icctemDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'biogeochemical (CTEM) PFTs')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'icctem', intValues=icctem, count = [size(icctem)])

        call ncRedef(fileId)

        ! Define the iccp2 dimension
        iccp2DimId = ncDefDim(fileId, 'iccp2', size(iccp2))
        varid = ncDefVar(fileId, 'iccp2', nf90_int, [iccp2DimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'biogeochemical (CTEM) PFTs plus bareground and land use change product pools')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'iccp2', intValues=iccp2, count = [size(iccp2)])

        call ncRedef(fileId)

        ! Define the months dimension
        ! monthsDimId = ncDefDim(fileId, 'months', size(months))
        ! varid = ncDefVar(fileId, 'months', nf90_int, [monthsDimId])
        ! call ncPutAtt(fileId, varId, 'standard_name', charValues = 'months')
        ! call ncEndDef(fileId)
        ! call ncPutDimValues(fileId, 'months', intValues=months, count = [size(months)])
        ! 
        ! call ncRedef(fileId)

        ! Define the slope dimension
        slopeDimId = ncDefDim(fileId, 'slope', size(slope))
        varid = ncDefVar(fileId, 'slope', nf90_int, [slopeDimId])
        call ncPutAtt(fileId, varId, 'standard_name', charValues = 'wetland slope fractions for 0.025, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3 and 0.35 percent slope threshold')
        call ncEndDef(fileId)
        call ncPutDimValues(fileId, 'slope', intValues=slope, count = [size(slope)])

        ! Add some generic variables that CLASSIC will look for

        ! per grid variables
        allocate(dimArray(2),start(2),count(2))
        dimArray = (/lonDimId,latDimId/)
        start = (/1, 1 /)
        count = (/1, 1 /)        
        call exportVariable('GC',units='-',long_name='GCM surface descriptor - land surfaces (inc. inland water) is -1',intvalues=(/-1/))
        call exportVariable('nmtest',units='-',long_name='Number of tiles in each grid cell',intvalues=(/nmtest/))
        deallocate(dimArray,start,count)

    end subroutine makeNetCDF

    ! ------------------------------------------------------------------------------------------------------

    subroutine exportData

        integer :: m,k

        ! icp1 variables:
        allocate(dimArray(4),start(4),count(4))
        dimArray = (/lonDimId,latDimId,icp1DimId,tileDimId/)
        start = (/1, 1, 1 ,1/)
        count = (/1, 1, ican+1, nmtest/)
        call exportVariable('FCAN',units='-',long_name='Annual maximum fractional coverage of modelled area (read in for CLASS only runs)',values2D=FCANROT)
        call exportVariable('LNZ0',units='-',long_name='Natural logarithm of maximum vegetation roughness length',values2D=LNZ0ROT)
        call exportVariable('ALIC',units='-',long_name='Average near-IR albedo of vegetation category when fully-leafed',values2D=ALICROT)
        call exportVariable('ALVC',units='-',long_name='Average visible albedo of vegetation category when fully-leafed',values2D=ALVCROT)


        ! ic variables:
        dimArray = (/lonDimId,latDimId,icDimId,tileDimId/)
        count = (/1, 1, ican, nmtest/)
        call exportVariable('PAMX',units='m2/m2',long_name='Annual maximum plant area index of vegetation category',values2D=PAMXROT)
        call exportVariable('PAMN',units='m2/m2',long_name='Annual minimum plant area index of vegetation category',values2D=PAMNROT)
        call exportVariable('CMAS',units='$[kg m^{-2} ]$',long_name='Annual maximum canopy mass for vegetation category',values2D=CMASROT)
        call exportVariable('ROOT',units='m',long_name='Annual maximum rooting depth of vegetation category',values2D=ROOTROT)
        ! The following are now read in from the model parameter namelist file. 
        !call exportVariable('RSMN',units='s/m',long_name='Minimum stomatal resistance of vegetation category',values2D=RSMNROT)
        !call exportVariable('QA50',units='W/m2',long_name='Reference value of incoming shortwave radiation (used in stomatal resistance calculation)',values2D=QA50ROT)
        !call exportVariable('VPDA',units='-',long_name='Vapour pressure deficit coefficient (used in stomatal resistance calculation)',values2D=VPDAROT)
        !call exportVariable('VPDB',units='-',long_name='Vapour pressure deficit coefficient (used in stomatal resistance calculation)',values2D=VPDBROT)
        !call exportVariable('PSGA',units='-',long_name='Soil moisture suction coefficient (used in stomatal resistance calculation)',values2D=PSGAROT)
        !call exportVariable('PSGB',units='-',long_name='Soil moisture suction coefficient (used in stomatal resistance calculation)',values2D=PSGBROT)

        ! ignd variables:
        dimArray = (/lonDimId,latDimId,layerDimId,tileDimId/)
        count = (/1, 1, ignd, nmtest/)
        call exportVariable('SAND',units='%',long_name='Percentage sand content',values2D=SANDROT)
        call exportVariable('CLAY',units='%',long_name='Percentage clay content',values2D=CLAYROT)
        call exportVariable('ORGM',units='%',long_name='Percentage organic matter content',values2D=ORGMROT)
        call exportVariable('TBAR',units='C',long_name='Temperature of soil layers',values2D=TBARROT)
        call exportVariable('THIC',units='m3/m3',long_name='Volumetric frozen water content of soil layers',values2D=THICROT)
        call exportVariable('THLQ',units='m3/m3',long_name='Volumetric liquid water content of soil layers',values2D=THLQROT)

        deallocate(dimArray,start,count)

        ! ignd variables:
        allocate(dimArray(1),start(1),count(1))
        dimArray = (/layerDimId/)
        count = (/ignd/)
        start = (/1/)
        call exportVariable('DELZ',units='m',long_name='Ground layer thickness',values=DELZ)
        deallocate(dimArray,start,count)

        
        ! nmtest only variables:
        allocate(dimArray(3),start(3),count(3))
        dimArray = (/lonDimId,latDimId,tileDimId/)
        start = (/1, 1, 1 /)
        count = (/1, 1, nmtest/)
        call exportVariable('DRN',units='-',long_name='Soil drainage index',values=DRNROT)
        call exportVariable('FARE',units='fraction',long_name='Tile fractional area of gridcell',values=FAREROT)
        call exportVariable('MID',units='-',long_name='Mosaic tile type identifier (1 for land surface, 0 for inland lake)',intvalues=MIDROT)
        call exportVariable('SDEP',units='m',long_name='Soil permeable depth',values=SDEPROT)
        !call exportVariable('XSLP',units='-',long_name='Not in Use: parameters lateral movement of soil water',values=XSLPROT)
        !call exportVariable('GRKF',units='-',long_name='Not in Use: parameters lateral movement of soil water',values=GRKFROT)
        !call exportVariable('WFCI',units='-',long_name='Not in Use: parameters lateral movement of soil water',values=WFCIROT)
        !call exportVariable('WFSF',units='-',long_name='Not in Use: parameters lateral movement of soil water',values=WFSFROT)
        call exportVariable('SOCI',units='index',long_name='Soil colour index',intvalues=SOCIROT)
        call exportVariable('TCAN',units='C',long_name='Vegetation canopy temperature',values=TCANROT)
        call exportVariable('ALBS',units='-',long_name='Soil drainage index',values=ALBSROT)
        call exportVariable('GRO',units='-',long_name='Vegetation growth index',values=GROROT)
        call exportVariable('RCAN',units='-',long_name='Intercepted liquid water stored on canopy',values=RCANROT)
        call exportVariable('RHOS',units='kg/m3',long_name='Density of snow',values=RHOSROT)
        call exportVariable('SCAN',units='kg/m2',long_name='Intercepted frozen water stored on canopy',values=SCANROT)
        call exportVariable('SNO',units='kg/m2',long_name='Mass of snow pack',values=SNOROT)
        call exportVariable('TPND',units='C',long_name='Temperature of ponded water',values=TPNDROT)
        call exportVariable('TSNO',units='C',long_name='Snowpack temperature',values=TSNOROT)
        call exportVariable('ZPND',units='m',long_name='Depth of ponded water on surface',values=ZPNDROT)

        deallocate(dimArray,start,count)

        ! per grid variables
        allocate(dimArray(2),start(2),count(2))
        dimArray = (/lonDimId,latDimId/)
        start = (/1, 1 /)
        count = (/1, 1 /)
        ! These are now in the joboptions file used to setup a model run.
        !call exportVariable('ZBLD',units='m',long_name='Atmospheric blending height for surface roughness length averaging',values=(/ZBLDROW/))
        !call exportVariable('ZRFH',units='m',long_name='Reference height associated with forcing air temperature and humidity',values=(/ZRFHROW/))
        !call exportVariable('ZRFM',units='m',long_name='Reference height associated with forcing wind speed',values=(/ZRFMROW/))
        call exportVariable('grclarea',units='km2',long_name='Area of grid cell',values=(/grclarea/))
        call exportVariable('twarmm',units='degree C',long_name='temperature of the warmest month (running averages in an e-folding sense)',values=(/twarmm/))
        call exportVariable('tcoldm',units='degree C',long_name='temperature of the warmest month (running averages in an e-folding sense)',values=(/tcoldm/))
        call exportVariable('gdd5',units='days',long_name='growing degree days above 5 C (running averages in an e-folding sense)',values=(/gdd5/))
        call exportVariable('aridity',units='-',long_name='aridity index, ratio of potential evaporation to precipitation (running averages in an e-folding sense)',values=(/aridity/))
        call exportVariable('srplsmon',units='months',long_name='number of months in a year with surplus water i.e. precipitation more than potential evaporation(running averages in an e-folding sense)',values=(/srplsmon/))
        call exportVariable('defctmon',units='months',long_name='number of months in a year with water deficit i.e. precipitation less than potential evaporation(running averages in an e-folding sense)',values=(/defctmon/))
        call exportVariable('anndefct',units='mm',long_name='annual water deficit(running averages in an e-folding sense)',values=(/anndefct/))
        call exportVariable('annsrpls',units='mm',long_name='annual water surplus (running averages in an e-folding sense)',values=(/annsrpls/))
        call exportVariable('annpcp',units='mm',long_name='annual precipitation (running averages in an e-folding sense)',values=(/annpcp/))
        call exportVariable('dry_season_length',units='months',long_name='annual maximum dry month length (running averages in an e-folding sense)',values=(/dry_season_length/))
        

        deallocate(dimArray,start,count)

        if (inclCTEM .or. fromnml) then

            ! icc variables
            allocate(dimArray(4),start(4),count(4))
            dimArray = (/lonDimId,latDimId,icctemDimId,tileDimId/)
            start = (/1, 1, 1 ,1/)
            count = (/1, 1, icc, nmtest/)
            !call exportVariable('ailcmin',units='m2/m2',long_name='Min. LAI for use with CTEM1 option only. Obsolete',values2D=ailcminrow)
            !call exportVariable('ailcmax',units='m2/m2',long_name='Max. LAI for use with CTEM1 option only. Obsolete',values2D=ailcmaxrow)
            call exportVariable('bleafmas',units='kgC/m2',long_name='Brown leaf mass',values2D=bleafmasrow)
            call exportVariable('gleafmas',units='kgC/m2',long_name='Green leaf mass',values2D=gleafmasrow)
            call exportVariable('stemmass',units='kgC/m2',long_name='Stem mass',values2D=stemmassrow)
            call exportVariable('rootmass',units='kgC/m2',long_name='Root mass',values2D=rootmassrow)
            call exportVariable('grwtheff',units='(kg C/m^2)/(m2/m2)',long_name='Growth efficiency. Change in biomass per year per unit max. LAI,for use in mortality subroutine',values2D=grwtheffrow)
            ! call exportVariable('tracerBLeafMass',units='kgC/m2',long_name='Universal tracer for Brown leaf mass',values2D=tracerBLeafMass)
            ! call exportVariable('tracerGLeafMass',units='kgC/m2',long_name='Universal tracer for Green leaf mass',values2D=tracerGLeafMass)
            ! call exportVariable('tracerStemMass',units='kgC/m2',long_name='Universal tracer for Stem mass',values2D=tracerStemMass)
            ! call exportVariable('tracerRootMass',units='kgC/m2',long_name='Universal tracer for Root mass',values2D=tracerRootMass)
            call exportVariable('lfstatus',units='-',long_name='Leaf status, see Phenology',intvalues2D=lfstatusrow)
            call exportVariable('pandays',units='-',long_name='Days with +ve new photosynthesis, see Phenology',intvalues2D=pandaysrow)
            deallocate(dimArray,start,count)

            ! Peat vars and active layer depths
            allocate(dimArray(3),start(3),count(3))
            dimArray = (/lonDimId,latDimId,tileDimId/)
            start = (/1, 1, 1/)
            count = (/1, 1, nmtest/)
            call exportVariable('ipeatland',units='-',long_name='Peatland flag: 0 = not a peatland, 1= bog, 2 = fen',intvalues=ipeatlandrow)
            call exportVariable('Cmossmas',units='kgC/m2',long_name='C in moss biomass',values=Cmossmas)
            !call exportVariable('tracerMossCMass',units='kgC/m2',long_name='Universal tracer for C in moss biomass',values=tracerMossCMass)
            call exportVariable('litrmsmoss',units='kgC/m2',long_name='Moss litter mass',values=litrmsmoss)
            !call exportVariable('tracerMossLitrMass',units='kgC/m2',long_name='Universal tracer for C in Moss litter mass',values=tracerMossLitrMass)
            call exportVariable('dmoss',units='m',long_name='Depth of living moss',values=dmoss)
            call exportVariable('maxAnnualActLyr',units='m',long_name='Active layer depth maximum over the e-folding period specified by parameter eftime',values=maxAnnualActLyr)
            deallocate(dimArray,start,count)

            allocate(dimArray(4),start(4),count(4))
            dimArray = (/lonDimId,latDimId,slopeDimId,tileDimId/)
            start = (/1, 1, 1 ,1/)
            count = (/1, 1, 8, nmtest/)
            call exportVariable('slopefrac',units='-',long_name='Slope-based fraction for dynamic wetlands',values2D=slopefrac)

            dimArray = (/lonDimId,latDimId,icctemDimId,tileDimId/)
            count = (/1, 1, icc, nmtest/)

            if (fileType == 'ini') then
                if (icc .ne. 9 .and. ican .ne. 4) print*,'Warning - expected ICC =9 and ICAN = 4'
                do m = 1,nmtest
                    fcancmxrow(m,1) = FCANROT(m,1) * dvdfcanrow(m,1)
                    fcancmxrow(m,2) = FCANROT(m,1) * dvdfcanrow(m,2)
                    fcancmxrow(m,3) = FCANROT(m,2) * dvdfcanrow(m,3)
                    fcancmxrow(m,4) = FCANROT(m,2) * dvdfcanrow(m,4)
                    fcancmxrow(m,5) = FCANROT(m,2) * dvdfcanrow(m,5)
                    fcancmxrow(m,6) = FCANROT(m,3) * dvdfcanrow(m,6)
                    fcancmxrow(m,7) = FCANROT(m,3) * dvdfcanrow(m,7)
                    fcancmxrow(m,8) = FCANROT(m,4) * dvdfcanrow(m,8)
                    fcancmxrow(m,9) = FCANROT(m,4) * dvdfcanrow(m,9)
                end do
            !else - namelist reads in the fcancmx directly without this conversion
            end if
            call exportVariable('fcancmx',units='-',long_name='PFT fractional coverage per grid cell',values2D=fcancmxrow)

            ! iccp2 variables
            ! dimArray = (/lonDimId,latDimId,iccp2DimId,layerDimId,tileDimId/)
            ! start = (/1, 1, 1 ,1, 1/)
            ! count = (/1, 1, icc+2, ignd, nmtest/)
            ! call exportVariable('litrmass',units='kgC/m2',long_name='Litter mass per soil layer',values3D=litrmassrow)
            ! call exportVariable('soilcmas',units='kgC/m2',long_name='Soil C mass per soil layer',values3D=soilcmasrow)
            ! call exportVariable('tracerLitrMass',units='kgC/m2',long_name='Universal tracer for Litter mass per soil layer',values3D=tracerLitrMass)
            ! call exportVariable('tracerSoilCMass',units='kgC/m2',long_name='Universal tracer for Soil C mass per soil layer',values3D=tracerSoilCMass)
            dimArray = (/lonDimId,latDimId,iccp2DimId,tileDimId/)
            start = (/1, 1, 1 ,1/)
            count = (/1, 1, icc+2, nmtest/)
            call exportVariable('litrmass',units='kgC/m2',long_name='Litter mass per soil layer',values2D=litrmassrow)
            call exportVariable('soilcmas',units='kgC/m2',long_name='Soil C mass per soil layer',values2D=soilcmasrow)
            
            deallocate(dimArray,start,count)

            end if

    end subroutine exportData

    ! ------------------------------------------------------------------------------------------------------

    subroutine exportVariable(name,units,long_name,values,values2D,values3D,intvalues,intvalues2D)
        character(*), intent(in)    :: name
        real, intent(in),optional   :: values(:),values2D(:,:),values3D(:,:,:)
        integer, intent(in), optional :: intvalues(:),intvalues2D(:,:)
        character(*), intent(in)    :: units
        character(*), intent(in)    :: long_name
        integer                     :: varId, m, k,icnum(1)
        integer, allocatable        :: incstart(:),usecount(:)

        call ncRedef(fileId)

        varid = ncDefVar(fileId, name, nf90_double, dimArray)
        call ncPutAtt(fileId, varid, 'units', charValues = units)
        call ncPutAtt(fileId, varid, '_FillValue', realValues = fillValue)
        call ncPutAtt(fileId, varId, 'long_name', charValues = long_name)
        call ncEndDef(fileId)

        allocate(incstart(size(start)))
        allocate(usecount(size(count)))
        
        ! Put in data
        if (present(values)) then
            call ncPutVar(fileId, name, realValues = values, start= start, count = count)
        else if (present(values2D)) then
          ! Need to change the start and count to respect nmtest here. Count is always by tile so set to 1 whereas
          ! start increments
            usecount = count
            usecount(ubound(usecount))=1
            do m = 1,nmtest
              incstart = start
              incstart(ubound(start)) = m          
                call ncPutVar(fileId, name, realValues = values2D(m,:), start= incstart, count = usecount)
            end do
        else if (present(values3D)) then
          ! Need to change the start and count to respect nmtest here. Count is always by tile so set to 1 whereas
          ! start increments
          !count = (/1, 1, icc+2, ignd, nmtest/)
          ! assume the extra dim is ignd
            usecount = count
            icnum = usecount(ubound(usecount)-1)
            usecount(ubound(usecount))=1
            usecount(ubound(usecount)-1)=1            
            do m = 1,nmtest        
              do k = 1, icnum(1)       
                incstart = start
                incstart(ubound(start)) = m          
                incstart(ubound(start)-1) = k
                call ncPutVar(fileId, name, realValues = values3D(m,:,k), start= incstart, count = usecount)
              end do
            end do
            
        else if (present(intvalues)) then
            call ncPutVar(fileId, name, intValues = intvalues, start= start, count = count)
        else if (present(intvalues2D)) then
            usecount = count
            usecount(ubound(usecount))=1          
            do m = 1,nmtest
                incstart = start
                incstart(ubound(start)) = m
                call ncPutVar(fileId, name, intValues = intvalues2D(m,:), start= incstart, count = usecount)
            end do
        else
            print*,'Problem in exportVariable'
        end if

        deallocate(usecount,incstart)
        
    end subroutine exportVariable

    ! ------------------------------------------------------------------------------------------------------

    logical function fileExists(filename)
        character(*), intent(in) :: filename
        inquire(file = filename, exist = fileExists)
    end function fileExists

    ! ------------------------------------------------------------------------------------------------------

    real function charToReal(input)
        character(len=*), intent(in)    :: input    !< Char input
        read(input,*) charToReal
    end function charToReal
    ! ------------------------------------------------------------------------------------------------------

end program initFileConverter
