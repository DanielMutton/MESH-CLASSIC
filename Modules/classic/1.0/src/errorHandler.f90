!> \file
!> Prints the name of the subroutine and an error code when
!! an error condition is encountered.
!! @author J. D. Henderson, E. Chan
!
subroutine errorHandler (NAME, N) ! Formerly XIT

  use generalUtilsGAT, only : abandonCell

  implicit none
  !
  !     * OCT 01/92 - E.CHAN. (CHANGE STOP 1 TO STOP)
  !     * JUN 10/91 - E.CHAN. (TRANSLATE HOLLERITH LITERALS AND
  !     *                      DIMENSION STRINGS)
  !
  !     * OCT 10/78 - J.D.HENDERSON.
  !     * TERMINATES A PROGRAM BY PRINTING THE PROGRAM NAME AND
  !     * A LINE ACROSS THE PAGE FOLLOWED BY A NUMBER N.
  !
  !     * N>=0 IS FOR A NORMAL END. THE LINE IS DASHED.
  !     * NORMAL ENDS TERMINATE WITH   STOP.
  !
  !     * N<0 IS FOR AN ABNORMAL END. THE LINE IS DOTTED.
  !     * IF N IS LESS THAN -100 THE PROGRAM SIMPLY TERMINATES.
  !     * OTHERWISE IF N IS LESS THAN ZERO THE PROGRAM ABORTS.
  !
  character * ( * ), intent(in) :: NAME    !< Name of the subroutine in which the error
  !< was found
  !< N: error code
  character * 8  :: NAME8, DASH, STAR
  integer, intent(in) :: N
  integer             :: i
  !
  !ignoreLint(1)
  DATA DASH /'--------'/, STAR /'********'/
  !---------------------------------------------------------------------
  !
  !>
  !! In CLASS, this subroutine is called when a test of ambient values
  !! of selected variables is performed and an abnormal condition is
  !! encountered. The name of the subroutine in which the condition
  !! arose is passed in, and is printed together with an error code,
  !! flagging the location of the error in the subroutine. A call to
  !! abort is then executed.
  !!
  NAME8 = NAME
  if (N >= 0) write(6,6010) DASH,NAME8,(DASH,I = 1,9),N
  !
  if (N < 0) write(6,6010) STAR,NAME8,(STAR,I = 1,9),N
  !
  if (N >= 0 .or. N < - 100) then
    ! CALL EXIT
    call abandonCell
  else
    ! CALL ABORT
    call abandonCell
  end if
  !
  !---------------------------------------------------------------------
  !ignoreLint(1)
6010 FORMAT('0',A8,'  END  ',A8,9A8,I8)
end subroutine errorHandler
