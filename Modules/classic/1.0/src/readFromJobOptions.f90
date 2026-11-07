!> \file
!> Parses command line arguments to program and reads in joboptions file.
module readJobOpts

  implicit none

  public :: readFromJobOptions
  public :: parsecoords

contains

  ! ---------------------------------------------------
  !> \ingroup readjobopts_readFromJobOptions
  !! @{
  !> Reads from the joboptions file, assigns the model switches, and determines the geographic domain
  !! of the simulation. All switches are described in the configurationFiles/template_job_options_file.txt file
  !! and also the user manual (housed in type ctem_switches in ctemStateVars.f90.
  !> @author Joe Melton

  subroutine readFromJobOptions

    use outputManager, only : myDomain
    use ctemStateVars, only : c_switch
    use classicParams, only : runParamsFile, PFTCompetitionSwitch, &
                              zbldJobOpt, zrfhJobOpt, zrfmJobOpt

    implicit none

    ! -------------

    logical, pointer :: projectedGrid

    ! ctem model switches

    integer, pointer :: metLoop
    logical, pointer :: ctem_on
    integer, pointer :: readMetStartYear
    integer, pointer :: readMetEndYear
    logical, pointer :: lnduseon
    integer, pointer :: spinfast
    integer, pointer :: useTracer
    character(:), pointer :: tracerCO2file
    logical, pointer :: transientCO2
    character(:), pointer :: CO2File
    integer, pointer :: fixedYearCO2
    logical, pointer :: doMethane
    logical, pointer :: transientCH4
    character(:), pointer :: CH4File
    integer, pointer :: fixedYearCH4
    logical, pointer :: transientPOPD
    character(:), pointer :: POPDFile
    integer, pointer :: fixedYearPOPD
    logical, pointer :: dofire
    integer, pointer :: fixedYearLGHT
    logical, pointer :: transientLGHT
    character(:), pointer :: LGHTFile
    character(:), pointer :: LUCFile
    integer, pointer :: fixedYearLUC
    logical, pointer :: transientOBSWETF
    character(:), pointer :: OBSWETFFile
    integer, pointer  :: fixedYearOBSWETF
    logical, pointer :: PFTCompetition
    logical, pointer :: inibioclim
    logical, pointer :: start_bare
    character(:), pointer :: metFileFss
    character(:), pointer :: metFileFdl
    character(:), pointer :: metFilePre
    character(:), pointer :: metFileTa
    character(:), pointer :: metFileQa
    character(:), pointer :: metFileUv
    character(:), pointer :: metFilePres
    character(:), pointer :: init_file
    character(:), pointer :: rs_file_to_overwrite
    character(:), pointer :: runparams_file
    character(:), pointer :: output_directory
    character(:), pointer :: xmlFile
    logical, pointer :: leap

    ! -------------
    ! class model switches

    integer, pointer :: idisp
    integer, pointer :: izref
    integer, pointer :: islfd
    integer, pointer :: ipcp
    integer, pointer :: iwf                 !modified to match change in iwf in classStateVars - Daniel Mutton
    integer, pointer :: ITC
    integer, pointer :: ITCG
    integer, pointer :: ITG
    integer, pointer :: IPAI
    integer, pointer :: IHGT
    integer, pointer :: IALC
    integer, pointer :: IALS
    integer, pointer :: IALG
    integer, pointer :: isnoalb

    ! -------------
    ! Output switches

    integer, pointer :: jhhstd  !< day of the year to start writing the half-hourly output
    integer, pointer :: jhhendd !< day of the year to stop writing the half-hourly output
    integer, pointer :: jdstd   !< day of the year to start writing the daily output
    integer, pointer :: jdendd  !< day of the year to stop writing the daily output
    integer, pointer :: jhhsty  !< simulation year (iyear) to start writing the half-hourly output
    integer, pointer :: jhhendy !< simulation year (iyear) to stop writing the half-hourly output
    integer, pointer :: jdsty   !< simulation year (iyear) to start writing the daily output
    integer, pointer :: jdendy  !< simulation year (iyear) to stop writing the daily output
    integer, pointer :: jmosty    !< Year to start writing out the monthly output files. If you want to write monthly outputs right
    logical, pointer :: doperpftoutput    !< Switch for making extra output files that are at the per PFT level
    logical, pointer :: dopertileoutput    !< Switch for making extra output files that are at the per tile level
    logical, pointer :: doAnnualOutput    !< Switch for making annual output files
    logical, pointer :: doMonthOutput    !< Switch for making monthly output files
    logical, pointer :: doDayOutput    !< Switch for making daily output files
    logical, pointer :: doHhOutput    !< Switch for making half hourly output files
    logical, pointer :: doChecksums
    character(:), pointer :: Comment   !< Comment about the run that will be written to the output netcdfs

    character(350) :: jobfile
    character(350) :: argbuff
    integer :: argcount, iargc
    real :: ZBLD, ZRFH, ZRFM

    ! Order of the namelist and order in the file don't have to match.

    namelist /joboptions/ &
        projectedGrid, &
        metLoop, &
        readMetStartYear, &
        readMetEndYear, &
        leap, &
        ctem_on, &
        spinfast, &
        useTracer, &
        tracerCO2file, &
        transientCO2, &
        CO2File, &
        fixedYearCO2, &
        doMethane, &
        transientCH4, &
        CH4File, &
        fixedYearCH4, &
        transientPOPD, &
        POPDFile, &
        fixedYearPOPD, &
        lnduseon, &
        LUCFile, &
        fixedYearLUC, &
        PFTCompetition, &
        inibioclim, &
        start_bare, &
        dofire, &
        transientLGHT, &
        fixedYearLGHT, &
        LGHTFile, &
        transientOBSWETF, &
        OBSWETFFile, &
        fixedYearOBSWETF, &
        metFileFss, &
        metFileFdl, &
        metFilePre, &
        metFileTa, &
        metFileQa, &
        metFileUv, &
        metFilePres, &
        init_file, &
        rs_file_to_overwrite, &
        runparams_file, &
        IDISP, &
        IZREF, &
        ZBLD, &
        ZRFH, &
        ZRFM, &
        ISLFD, &
        IPCP, &
        ITC, &
        ITCG, &
        ITG, &
        IWF, &
        IPAI, &
        IHGT, &
        IALC, &
        IALS, &
        IALG, &
        isnoalb, &
        output_directory, &
        xmlFile, &
        doperpftoutput, &
        dopertileoutput, &
        doHhOutput, &
        JHHSTD, &
        JHHENDD, &
        JHHSTY, &
        JHHENDY, &
        doDayOutput, &
        JDSTD, &
        JDENDD, &
        JDSTY, &
        JDENDY, &
        doMonthOutput, &
        JMOSTY, &
        doAnnualOutput, &
        doChecksums, &
        Comment

    ! Point pointers:
    projectedGrid   => c_switch%projectedGrid
    metLoop         => c_switch%metLoop
    readMetStartYear=> c_switch%readMetStartYear
    readMetEndYear  => c_switch%readMetEndYear
    ctem_on         => c_switch%ctem_on
    lnduseon        => c_switch%lnduseon
    LUCFile         => c_switch%LUCFile
    fixedYearLUC    => c_switch%fixedYearLUC
    spinfast        => c_switch%spinfast
    useTracer       => c_switch%useTracer
    tracerCO2file   => c_switch%tracerCO2file
    transientCO2    => c_switch%transientCO2
    CO2File         => c_switch%CO2File
    fixedYearCO2    => c_switch%fixedYearCO2
    doMethane       => c_switch%doMethane
    transientCH4    => c_switch%transientCH4
    CH4File         => c_switch%CH4File
    fixedYearCH4    => c_switch%fixedYearCH4
    transientPOPD   => c_switch%transientPOPD
    POPDFile        => c_switch%POPDFile
    fixedYearPOPD   => c_switch%fixedYearPOPD
    dofire          => c_switch%dofire
    fixedYearLGHT   => c_switch%fixedYearLGHT
    transientLGHT   => c_switch%transientLGHT
    LGHTFile        => c_switch%LGHTFile
    transientOBSWETF=> c_switch%transientOBSWETF
    OBSWETFFile     => c_switch%OBSWETFFile
    fixedYearOBSWETF=> c_switch%fixedYearOBSWETF
    PFTCompetition  => c_switch%PFTCompetition
    inibioclim      => c_switch%inibioclim
    start_bare      => c_switch%start_bare
    rs_file_to_overwrite => c_switch%rs_file_to_overwrite
    output_directory => c_switch%output_directory
    xmlFile         => c_switch%xmlFile
    metFileFss      => c_switch%metFileFss
    metFileFdl      => c_switch%metFileFdl
    metFilePre      => c_switch%metFilePre
    metFileTa       => c_switch%metFileTa
    metFileQa       => c_switch%metFileQa
    metFileUv       => c_switch%metFileUv
    metFilePres     => c_switch%metFilePres
    runparams_file  => c_switch%runparams_file
    init_file       => c_switch%init_file
    IDISP           => c_switch%IDISP
    IZREF           => c_switch%IZREF
    ISLFD           => c_switch%ISLFD
    IPCP            => c_switch%IPCP
    ITC             => c_switch%ITC
    ITCG            => c_switch%ITCG
    ITG             => c_switch%ITG
    IWF             => c_switch%IWF
    IPAI            => c_switch%IPAI
    IHGT            => c_switch%IHGT
    IALC            => c_switch%IALC
    IALS            => c_switch%IALS
    IALG            => c_switch%IALG
    isnoalb         => c_switch%isnoalb
    leap            => c_switch%leap
    jhhstd          => c_switch%jhhstd
    jhhendd         => c_switch%jhhendd
    jdstd           => c_switch%jdstd
    jdendd          => c_switch%jdendd
    jhhsty          => c_switch%jhhsty
    jhhendy         => c_switch%jhhendy
    jdsty           => c_switch%jdsty
    jdendy          => c_switch%jdendy
    jmosty          => c_switch%jmosty
    doAnnualOutput  => c_switch%doAnnualOutput
    doperpftoutput  => c_switch%doperpftoutput
    dopertileoutput => c_switch%dopertileoutput
    doMonthOutput   => c_switch%doMonthOutput
    doDayOutput     => c_switch%doDayOutput
    doHhOutput      => c_switch%doHhOutput
    doChecksums     => c_switch%doChecksums
    Comment         => c_switch%Comment

    !-------------------------
    ! read the joboptions

    argcount = iargc()

    if (argcount /= 2) then
      write( * , * )'Usage is as follows'
      write( * , * )' '
      write( * , * )'bin/CLASSIC joboptions_file longitude/{longitude}/latitude/{latitude}'
      write( * , * )' '
      write( * , * )' - joboptions_file - an example is '
      write( * , * )'  configurationFiles/template_job_options_file.txt.'
      write( * , * )' '
      write( * , * )' - longitude/latitude '
      write( * , * )'  e.g. 105.23/40.91 '
      write( * , * )' '
      write( * , * )' * OR * '
      write( * , * )' if you wish to run a region then you give '
      write( * , * )' the corners of the box you wish to run '
      write( * , * )' '
      write( * , * )' - longitude/longitude/latitude/latitude '
      write( * , * )'  e.g. 90/105/30/45'
      write( * , * )' '
      write( * , * )' ** If you are running a projected grid you must'
      write( * , * )' use the grid cell indices, not coordinates !** '
      write( * , * )' '
      stop
    end if

    !> Argument 1 is the jobfile, which is openned and the namelist is read
    call getarg(1,jobfile)

    open(10,file = jobfile,action = 'read',status = 'old')
	
	print *, "projectedGrid ", projectedGrid
	print *, "metLoop ", metLoop
	print *, "readMetStartYear ", readMetStartYear
	print *, "readMetEndYear ", readMetEndYear
	print *, "leap ", leap
	print *, "ctem_on ", ctem_on
	print *, "spinfast ", spinfast
	print *, "useTracer ", useTracer
	print *, "tracerCO2file ", tracerCO2file
	print *, "transientCO2 ", transientCO2
	print *, "CO2File ", CO2File							! could be the issue
	print *, "fixedYearCO2 ", fixedYearCO2
	print *, "doMethane ", doMethane
	print *, "transientCH4 ", transientCH4
	print *, "CH4File ", CH4File
	print *, "fixedYearCH4 ", fixedYearCH4
	print *, "transientPOPD ", transientPOPD
	print *, "POPDFile ", POPDFile
	print *, "fixedYearPOPD ", fixedYearPOPD
	print *, "lnduseon ", lnduseon
	print *, "LUCFile ", LUCFile
	print *, "fixedYearLUC ", fixedYearLUC
	print *, "PFTCompetition ", PFTCompetition
	print *, "inibioclim ", inibioclim
	print *, "start_bare ", start_bare
	print *, "dofire ", dofire
	print *, "transientLGHT ", transientLGHT
	print *, "fixedYearLGHT ", fixedYearLGHT
	print *, "LGHTFile ", LGHTFile
	print *, "transientOBSWETF ", transientOBSWETF
	print *, "OBSWETFFile ", OBSWETFFile
	print *, "fixedYearOBSWETF ", fixedYearOBSWETF
	print *, "metFileFss ", metFileFss
	print *, "metFileFdl ", metFileFdl
	print *, "metFilePre ", metFilePre
	print *, "metFileTa ", metFileTa
	print *, "metFileQa ", metFileQa
	print *, "metFileUv ", metFileUv
	print *, "metFilePres ", metFilePres
	print *, "init_file ", init_file
	print *, "rs_file_to_overwrite ", rs_file_to_overwrite
	print *, "runparams_file ", runparams_file
	print *, "IDISP ", IDISP
	print *, "IZREF ", IZREF
	print *, "ZBLD ", ZBLD
	print *, "ZRFH ", ZRFH
	print *, "ZRFM ", ZRFM
	print *, "ISLFD ", ISLFD
	print *, "IPCP ", IPCP
	print *, "ITC ", ITC
	print *, "ITCG ", ITCG
	print *, "ITG ", ITG
	print *, "IWF ", IWF
	print *, "IPAI ", IPAI
	print *, "IHGT ", IHGT
	print *, "IALC ", IALC
	print *, "IALS ", IALS
	print *, "IALG ", IALG
	print *, "isnoalb ", isnoalb
	print *, "output_directory ", output_directory
	print *, "xmlFile ", xmlFile
	print *, "doperpftoutput ", doperpftoutput
	print *, "dopertileoutput ", dopertileoutput
	print *, "doHhOutput ", doHhOutput
	print *, "JHHSTD ", JHHSTD
	print *, "JHHENDD ", JHHENDD
	print *, "JHHSTY ", JHHSTY
	print *, "JHHENDY ", JHHENDY
	print *, "doDayOutput ", doDayOutput
	print *, "JDSTD ", JDSTD
	print *, "JDENDD ", JDENDD
	print *, "JDSTY ", JDSTY
	print *, "JDENDY ", JDENDY
	print *, "doMonthOutput ", doMonthOutput
	print *, "JMOSTY ", JMOSTY
	print *, "doAnnualOutput ", doAnnualOutput
	print *, "doChecksums ", doChecksums
	print *, "Comment ", Comment
	
    read(10,nml = joboptions, iostat=spinfast)
	print *, spinfast

    close(10)

    !> Parse the 2nd argument to get the domain that the simulation should be run over
    call getarg(2, argbuff)
    call parsecoords(argbuff, myDomain%domainBounds)

    ! Assign some vars that are passed out
    runParamsFile = runparams_file
    PFTCompetitionSwitch = PFTCompetition
    zbldJobOpt = ZBLD
    zrfhJobOpt = ZRFH
    zrfmJobOpt = ZRFM

  end subroutine readFromJobOptions
  !! @}
  ! ----------------------------------------------------------------------------------

  !> \ingroup readjobopts_parsecoords
  !! @{
  !> Parses a coordinate string
  !> @author Joe Melton

  subroutine parsecoords (coordstring, val)

    implicit none

    character(45), intent(in)  :: coordstring
    real, dimension(4), intent(out) :: val

    character(10), dimension(4) :: cval = '0'

    integer :: i
    integer :: lasti = 1
    integer :: part  = 1

    do i = 1,len_trim(coordstring)
      if (coordstring(i:i) == '/') then
        cval(part) = coordstring(lasti:i - 1)
        lasti = i + 1
        part = part + 1
      end if
    end do

    cval(part) = coordstring(lasti:i - 1)

    read(cval, * )val

    if (part < 4) then
      val(3) = val(2)
      val(4) = val(3)
      val(2) = val(1)
    end if

  end subroutine parsecoords
  !! @}

  !> \namespace readjobopts
  !> Parses command line arguments to program and reads in joboptions file.

end module readJobOpts
