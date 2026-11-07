!> \file
!! Principle driver program to run CLASSIC in stand-alone mode using specified boundary
!! conditions and atmospheric forcing. Depending upon the compiler options chosen, this
!! program will either run in MPI mode for use on parallel computing environments or
!! in serial mode for use in running at single sites.
!!

program CLASSIC

#if PARALLEL
    use mpi
#endif

  use modelStateDrivers,    only : read_modelsetup, rsid
  use xmlManager,           only : loadoutputDescriptor
  use outputManager,        only : generateOutputFiles, closeNCFiles
  use readJobOpts,          only : readFromJobOptions
  use main,                 only : main_driver
  use ctemStateVars,        only : allocCtemVars
  use classStateVars,       only : allocClassVars
  use classicParams,        only : prepareGlobalParams

  implicit none

  double precision  :: time
  integer           :: ierr, rank, size, i, cell, blocks, remainder

  ! MAIN PROGRAM

  ! First initialize the MPI and PnetCDF session. Ignored if run in serial mode.
  call initializeParallelEnvironment

  ! Next load the job options file. This first parses the command line arguments.
  ! Then all model switches are read in from a namelist file. This sets up the
  ! run options and points to input files as needed.
  call readFromJobOptions

  ! Then load the run setup information based on the metadata in the
  ! initialization netcdf file. The bounds given as an argument to
  ! CLASSIC are used to find the start points (srtx and srty)
  ! in the netcdf file, placing the gridcell on the domain of the
  ! input/output netcdfs. In read_modelsetup we use the netcdf to set
  ! the nmos (number of tiles), ignd (number of soil layers), and ilg (number of latitude
  ! points times nmos, which defaults to nmos in offline mode) constants.
  ! It also opens the initial conditions file that is used below in
  ! read_initialstate as well as the restart file that is written to later.
  call read_modelsetup
  
  ! Prepare all of the global parameters in classicParams
  ! which are read from a namelist file.
  call prepareGlobalParams
  
  ! The output files are created based on the model switches in the
  ! joboptions file and the xml file that describes the metadata for
  ! each output file. The loadoutputDescriptor parses the xml file and
  ! creates a data structure to allow us to make all of the netcdf output
  ! files, one per variable per time period (daily, monthly, etc.).
  call loadoutputDescriptor
  
  ! Generate the output files based on options in the joboptions file
  ! and the parameters of the initilization netcdf file.
  call generateOutputFiles
  
  ! Run model over the land grid cells, in parallel or serial
  call processLandCells
  
  ! Close all of the output netcdf files and the restart file
  ! (these were written to so need to ensure buffer is flushed)
  call closeNCFiles
  call closeNCFiles(rsid)
  
#if PARALLEL
    !> Shut down the MPI session
    call MPI_FINALIZE(ierr)
#endif

  ! END MAIN PROGRAM

  !------------------

contains
  !> \ingroup CLASSIC_processLandCells
  !! @{
  !> Runs the model over all of the land cells.
  subroutine processLandCells

    ! PROCESS LAND CELLS
    ! This section runs the model over all of the land cells. There are LandCellCount valid(i.e. land) cells, stored in lonLandCell and latLandCell

    use outputManager,       only : myDomain
    use fileIOModule

    implicit none
        
    ! Since we know the nlat, nmos, ignd, and ilg we can allocate the CLASS and
    ! CTEM variable structures. This has to be done before call to main_driver.
    call allocClassVars()
    call allocCtemVars()

    blocks = myDomain%LandCellCount / size + 1          ! The number of processing blocks
    remainder = mod(myDomain%LandCellCount,size)       ! The number of cells for the last block

    do i = 1,blocks - 1                    ! Go through every block except for the last one
      cell = (i - 1) * size + rank + 1
      print * ,'in process',rank,i,cell,myDomain%lonLandCell(cell),myDomain%latLandCell(cell),myDomain%lonLandIndex(cell),myDomain%latLandIndex(cell), &
                                myDomain%lonLocalIndex(cell),myDomain%latLocalIndex(cell)
      call main_driver(myDomain%lonLandCell(cell), myDomain%latLandCell(cell), myDomain%lonLandIndex(cell), myDomain%latLandIndex(cell), &
                       myDomain%lonLocalIndex(cell), myDomain%latLocalIndex(cell))
    end do

    cell = (blocks - 1) * size + rank + 1   ! In the last block,process only the existing cells

    if (rank < remainder) print * ,'final in process',cell,myDomain%lonLandCell(cell),myDomain%latLandCell(cell), &
        myDomain%lonLandIndex(cell),myDomain%latLandIndex(cell),myDomain%lonLocalIndex(cell),myDomain%latLocalIndex(cell)
    if (rank < remainder) call main_driver(myDomain%lonLandCell(cell),myDomain%latLandCell(cell), &
        myDomain%lonLandIndex(cell),myDomain%latLandIndex(cell),myDomain%lonLocalIndex(cell),myDomain%latLocalIndex(cell))
    
  end subroutine processLandCells
  !! @}
  !----------------------------------------------------------------------------------------

  !> \ingroup CLASSIC_initializeParallelEnvironment
  !! @{
  !> If compiled in parallel mode, initializes MPI variables.
  subroutine initializeParallelEnvironment
    implicit none

    size = 1
    rank = 0
#if PARALLEL
    call MPI_INIT(ierr)
    time = MPI_WTIME()
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
#endif
  end subroutine initializeParallelEnvironment
  !! @}
  ! ---------------------------------------------------------------------------------------

  !> \namespace CLASSIC

  !>
  !! @author J. Melton, E. Wisernig
  !!
  !! The order of events of the main program are summarized here:
  !!
  !! (initializeParallelEnvironment) -
  !! First initialize the MPI and PnetCDF session. Ignored if run in serial mode.
  !!
  !! (readFromJobOptions) -
  !! Next load the job options file. This first parses the command line arguments.
  !! Then all model switches are read in from a namelist file. This sets up the
  !! run options and points to input files as needed.
  !!
  !! (read_modelsetup) -
  !! Then load the run setup information based on the metadata in the
  !! initialization netcdf file. The bounds given as an argument to
  !! CLASSIC are used to find the start points (srtx and srty)
  !! in the netcdf file, placing the gridcell on the domain of the
  !! input/output netcdfs. In read_modelsetup we use the netcdf to set
  !! the nmos (number of tiles), ignd (number of soil layers), and ilg (number of latitude
  !! points times nmos, which defaults to nmos in offline mode) constants.
  !! It also opens the initial conditions file that is used below in
  !! read_initialstate as well as the restart file that is written to later.
  !!
  !! (prepareGlobalParams) -
  !! Prepare all of the global parameters in classicParams
  !! which are read from a namelist file.
  !!
  !! (loadoutputDescriptor) -
  !! The output files are created based on the model switches in the
  !! joboptions file and the xml file that describes the metadata for
  !! each output file. The loadoutputDescriptor parses the xml file and
  !! creates a data structure to allow us to make all of the netcdf output
  !! files, one per variable per time period (daily, monthly, etc.).
  !!
  !! (generateOutputFiles) -
  !! Generate the output files based on options in the joboptions file
  !! and the parameters of the initilization netcdf file.
  !!
  !! (processLandCells) -
  !! Run model over the land grid cells, in parallel or serial
  !!
  !! (closeNCFiles) -
  !! Close all of the output netcdf files and the restart file
  !! (these were written to so need to ensure buffer is flushed)
  !!
  !> \file
end program CLASSIC
