!> \file
!> The XML Manager Module handles the loading of variable descriptors and variable variants from an XML file.
!! Once the variables are generated, corresponding NetCDF output files are created.
!> @author
!> Ed Wisernig

module xmlManager
  use outputManager, only : outputDescriptor, outputDescriptors, descriptorCount, &
                                variant, variants, variantCount, checkFileExists
  use xmlParser, only : xml_process
  implicit none
  public :: loadoutputDescriptor
  public :: startfunc
  public :: datafunc
  public :: endfunc
  private :: charToLogical
  private :: charToInt

  character(len = 80), dimension(2,10)      :: attribs
  character(len = 400), dimension(100)      :: data
  logical                                   :: error, currentDormant = .true.
  character(len = 80)                       :: currentGroup, currentVariableName, variableSetType, variableSetDate, variableSetVersion
  real                                      :: xmlVersion
  character(len=40)                         :: dormant

contains
  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> The loadoutputDescriptor function parses the XML file using the startfunc, datafunc and endfunc for starting tags, data tags and end tags, respectively.
  subroutine loadoutputDescriptor
    use ctemStateVars, only : c_switch
    implicit none
    character(:),pointer :: xmlFile
    logical              :: fileExists
    xmlFile         => c_switch%xmlFile
    ! Check if the xmlFile exists:
    ! Now make sure the file was properly created
    fileExists = checkFileExists(xmlFile)

    if (.not. fileExists) then
      print * ,'Missing xml file: ',xmlFile
      print * ,'Aborting'
      stop ! can use stop here as not in MPI part of code.
    end if

    call xml_process(xmlFile, attribs, data, startfunc, datafunc, endfunc, 0, error)
  end subroutine loadoutputDescriptor

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> The startfunc function parses opening tags or start tags. It gets triggered for each and every opening tag of the XML document.
  subroutine startfunc (tag, attribs, error)
    implicit none
    character(len =* ), intent(in)      :: tag, attribs(:,:)
    character(len = 40)                 :: attribute
    integer                             :: id
    logical, intent(in)                 :: error
    type(outputDescriptor), allocatable :: tempDescriptors(:)
    type(variant), allocatable          :: tempVariants(:)

    ! Examine the tag
    select case (tag)
      ! If the tag is <variableSet>, allocate the variable descriptors and the variants.
    case ( 'variableSet' )
      xmlVersion = 0
      variableSetType = trim(attribs(2,1))
      variableSetVersion = trim(attribs(2,2))
      if (variableSetVersion == '') then
        stop ("The input XML document doesn't feature the required version field of the < variableSet > node")
      else
        xmlVersion = charToReal(variableSetVersion)
        if (xmlVersion < 1.2) stop ('Older XML document found, please upgrade to a more recent version')
      end if

      variableSetDate = trim(attribs(2,3))
      allocate(outputDescriptors(0))
      allocate(variants(0))
      ! If the tag is <group>, remember the current group.
    case ( 'group' )
      currentGroup = trim(attribs(2,1))
      ! If the tag is <variable>, increment the descriptor count and allocate the temporary descriptors.
      ! Then, add a new descriptor to the array, by copying and extending the array, followed by setting some values.
    case ( 'variable' )
      dormant = trim(attribs(2,3))
      if (dormant == 'false') then
        descriptorCount = descriptorCount + 1
        allocate(tempDescriptors(descriptorCount))
        tempDescriptors(1 : descriptorCount - 1) = outputDescriptors(1 : descriptorCount - 1)
        call move_alloc(tempDescriptors,outputDescriptors)
        attribute = attribs(2,2)
        outputDescriptors(descriptorCount)%includeBareGround = charToLogical(attribute)
        outputDescriptors(descriptorCount)%group = currentGroup
        currentDormant = .false.
      else
        currentDormant = .true.
      end if
      ! If the tag is <variant>, then similar to the above section, extend the array and append a new variant.
    case ('variant')
      if (dormant == 'false') then
        variantCount = variantCount + 1
        allocate(tempVariants(variantCount))
        tempVariants(1 : variantCount - 1) = variants(1 : variantCount - 1)
        call move_alloc(tempVariants,variants)
        variants(variantCount)%shortName = currentVariableName
      end if
    end select
  end subroutine

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> The datafunc function parses the content of a tag.
  subroutine datafunc (tag, data, error)
    implicit none
    character(len =* ), intent(in) :: tag, data(:)
    character(len = 400)           :: info
    logical, intent(in)            :: error

    info = trim(data(1))
    ! Examine the tag and store the tag content in the descriptors array.
    if (.not. currentDormant) then
      select case (tag)
      case ( 'shortName' )
        outputDescriptors(descriptorCount)%shortName = info
        currentVariableName = info
      case ( 'longName' )
        outputDescriptors(descriptorCount)%longName = info
      case ( 'standardName' )
        outputDescriptors(descriptorCount)%standardName = info
      case ( 'units' )
        ! outputDescriptors(descriptorCount)%units = info
        variants(variantCount)%units = info
      case ( 'nameInCode' )
        variants(variantCount)%nameInCode = info
      case ( 'timeFrequency' )
        variants(variantCount)%timeFrequency = info
      case ( 'outputForm' )
        variants(variantCount)%outputForm = info
      end select
    end if
  end subroutine

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> The endfunc function gets triggered when the closing XML element is encountered. E.g. if we wanted to deallocate at </variableSet>.
  subroutine endfunc (tag, error)
    implicit none
    character(len =* ), intent(in)   :: tag
    logical, intent(in)              :: error
    ! Place a select case (tag) in here if you need to do anything on a closing tag.
  end subroutine

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> The charToLogical function returns the logical value of a given char input
  logical function charToLogical (input)
    implicit none
    character(len =* ), intent(in)    :: input    !< Char input
    if (input == "true") then
      charToLogical = .true.
    else
      charToLogical = .false.
    end if
  end function charToLogical

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> The charToInt function returns the integer :: value of a given char input
  integer function charToInt (input)
    implicit none
    character(len =* ), intent(in)    :: input    !< Char input
    read(input, * ) charToInt
  end function charToInt

  !-----------------------------------------------------------------------------------------------------------------------------------------------------
  !> The charToInt function returns the integer :: value of a given char input
  real function charToReal (input)
    implicit none
    character(len =* ), intent(in) :: input    !< Char input
    read(input, * ) charToReal
  end function charToReal
  !> \namespace xmlmanager
end module xmlManager
