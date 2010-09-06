MODULE mod_logga
  USE mod_indata
  IMPLICIT  NONE

  INTERFACE LOGGA
    MODULE PROCEDURE LOGGA_I, LOGGA_R, LOGGA_R8
  END INTERFACE

CONTAINS

SUBROUTINE LOGGA_I(lev, event, iarg)
  IMPLICIT NONE
! opens file "statusfile" with unit 44 and logs the event "event"
! if lev > loglevel
  CHARACTER(LEN=*), INTENT(IN) :: event
  INTEGER, INTENT(IN) :: lev    ! 1 minor event, 3 major event
  INTEGER, INTENT(IN) :: iarg
!......................................................Tipi locali
  CHARACTER (LEN=18) :: ora
  CHARACTER (LEN=3) :: mese(12)= (/"GEN","FEB","MAR","APR",       &
       &         "MAG","GIU","LUG","AGO","SET","OTT","NOV","DIC"/)
  CHARACTER (LEN= 10):: chdate, chtime
  INTEGER :: val(8)
  IF ( (lev > loglevel) .AND. (mypenum_mpi == npeio_mpi) ) THEN
    OPEN(44, FILE=statusfile , STATUS="UNKNOWN", POSITION="APPEND")
    CALL DATE_AND_TIME(date= chdate, time= chtime, values= val)
    ora(1:10)= chdate(7:8)//mese(val(2))//chdate(1:4)//" "
    ora(11:18)= chtime(1:2)//":"//chtime(3:4)//":"//chtime(5:6)
    WRITE(44,*) ora//" "// REPEAT("+",lev) //" "//event, iarg
    CLOSE(44)
  END IF
END SUBROUTINE LOGGA_I

SUBROUTINE LOGGA_R(lev, event, rarg)
  IMPLICIT NONE
! opens file "statusfile" with unit 44 and logs the event "event"
! if lev > loglevel
  CHARACTER(LEN=*), INTENT(IN) :: event
  INTEGER, INTENT(IN) :: lev    ! 1 minor event, 3 major event
  REAL, OPTIONAL :: rarg
!......................................................Tipi locali
  CHARACTER (LEN=18) :: ora
  CHARACTER (LEN=3) :: mese(12)= (/"GEN","FEB","MAR","APR",       &
       &         "MAG","GIU","LUG","AGO","SET","OTT","NOV","DIC"/)
  CHARACTER (LEN= 10):: chdate, chtime
  INTEGER :: val(8)
  IF ( (lev > loglevel) .AND. (mypenum_mpi == npeio_mpi) ) THEN
    OPEN(44, FILE=statusfile , STATUS="UNKNOWN", POSITION="APPEND")
    CALL DATE_AND_TIME(date= chdate, time= chtime, values= val)
    ora(1:10)= chdate(7:8)//mese(val(2))//chdate(1:4)//" "
    ora(11:18)= chtime(1:2)//":"//chtime(3:4)//":"//chtime(5:6)
    IF (PRESENT(rarg)) THEN
      WRITE(44,*) ora//" "// REPEAT("+",lev) //" "//event, rarg
    ELSE
      WRITE(44,*) ora//" "// REPEAT("+",lev) //" "//event
    END IF
    CLOSE(44)
  END IF
END SUBROUTINE LOGGA_R

SUBROUTINE LOGGA_R8(lev, event, rarg)
  IMPLICIT NONE
! opens file "statusfile" with unit 44 and logs the event "event"
! if lev > loglevel
  CHARACTER(LEN=*), INTENT(IN) :: event
  INTEGER, INTENT(IN) :: lev    ! 1 minor event, 3 major event
  REAL*8 :: rarg
!......................................................Tipi locali
  CHARACTER (LEN=18) :: ora
  CHARACTER (LEN=3) :: mese(12)= (/"GEN","FEB","MAR","APR",       &
       &         "MAG","GIU","LUG","AGO","SET","OTT","NOV","DIC"/)
  CHARACTER (LEN= 10):: chdate, chtime
  INTEGER :: val(8)
  IF ( (lev > loglevel) .AND. (mypenum_mpi == npeio_mpi) ) THEN
    OPEN(44, FILE=statusfile , STATUS="UNKNOWN", POSITION="APPEND")
    CALL DATE_AND_TIME(date= chdate, time= chtime, values= val)
    ora(1:10)= chdate(7:8)//mese(val(2))//chdate(1:4)//" "
    ora(11:18)= chtime(1:2)//":"//chtime(3:4)//":"//chtime(5:6)
    WRITE(44,*) ora//" "// REPEAT("+",lev) //" "//event, rarg
    CLOSE(44)
  END IF
END SUBROUTINE LOGGA_R8


END MODULE mod_logga

