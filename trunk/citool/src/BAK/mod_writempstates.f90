!$$$$$$$$$$$$$$$$$$$$$$$$$  CItool  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_writempstates
  IMPLICIT  NONE
  SAVE

CONTAINS

SUBROUTINE WRITEMPSTATES_X( nblock, dimblock, ket,  &
     &  namespqn_e, spqn_e, namespqn_h, spqn_h,       &
     &  binfmt_e, binfmt_h, mpenergies, mpstates )
  USE mod_indata
  IMPLICIT NONE
! 
INTEGER, INTENT(IN) :: nblock
INTEGER, INTENT(IN) :: dimblock          ! Hilbert space basis for elecs and holes
INTEGER*8, INTENT(IN) :: ket(dimblock,2)
CHARACTER(LEN=12), INTENT(IN) :: namespqn_e(numspqn_e)
INTEGER, INTENT(IN) :: spqn_e(numspstates_e,numspqn_e)
CHARACTER(LEN=12), INTENT(IN) :: namespqn_h(numspqn_h)
INTEGER, INTENT(IN) :: spqn_h(numspstates_h,numspqn_h)
REAL*8, INTENT(IN) :: mpenergies(nummpenergies)
COMPLEX*16, INTENT(IN) :: mpstates(dimblock,nummpstates)

REAL*8 :: weight
INTEGER*8 :: kete, keth
INTEGER :: nmp, nsd
CHARACTER(LEN=80) :: binfmt_e, binfmt_h


OPEN(21, FILE=TRIM(fileout_mpstates),                    &
     &   ACTION="WRITE", STATUS="UNKNOWN", POSITION="APPEND", FORM="FORMATTED")

WRITE(21,*) "BLOCK", nblock

DO nmp= 1, nummpenergies
  WRITE(21,"(I3,3X)",ADVANCE="NO") nmp
  WRITE(21,*) mpenergies(nmp)

  IF (nmp <= nummpstates) THEN
    DO nsd= 1, dimblock
      weight= ABS(mpstates(nsd,nmp))**2
      IF (weight<1e-6) CYCLE
      kete= ket( nsd, 1 )
      keth= ket( nsd, 2 )

      WRITE(21,"('(',E9.3,',',E9.3,')',3X)",ADVANCE="NO") mpstates(nsd,nmp)
      WRITE(21,"(A,1X)",ADVANCE="NO") "e"
      WRITE(21,binfmt_e) kete
      WRITE(21,"('_______',F6.3,'_%','________ ')",ADVANCE="NO") weight*100
      WRITE(21,"(A,1X)",ADVANCE="NO") "h"
      WRITE(21,binfmt_h) keth

    END DO
  END IF
  WRITE(21,*) 
END DO

CLOSE( 21 )

END SUBROUTINE WRITEMPSTATES_X


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

SUBROUTINE WRITEMPSTATES( nblock, dimblock, ket,  &
     &  namespqn_e, spqn_e, namespqn_h, spqn_h,      &
     &  binfmt_e, binfmt_h, mpenergies, mpstates )
  USE mod_indata
  IMPLICIT NONE
! 
INTEGER, INTENT(IN) :: nblock
INTEGER, INTENT(IN) :: dimblock          ! Hilbert space basis for elecs and holes
INTEGER*8, INTENT(IN) :: ket(dimblock,2)
CHARACTER(LEN=12), INTENT(IN) :: namespqn_e(numspqn_e)
INTEGER, INTENT(IN) :: spqn_e(numspstates_e,numspqn_e)
CHARACTER(LEN=12), INTENT(IN) :: namespqn_h(numspqn_h)
INTEGER, INTENT(IN) :: spqn_h(numspstates_h,numspqn_h)
REAL*8, INTENT(IN) :: mpenergies(nummpenergies)
REAL*8, INTENT(IN) :: mpstates(dimblock,nummpstates)

REAL*8 :: weight
INTEGER*8 :: kete, keth
INTEGER :: nmp, nsd
CHARACTER(LEN=80) :: binfmt_e, binfmt_h


OPEN(21, FILE=TRIM(fileout_mpstates),                    &
     &   ACTION="WRITE", STATUS="UNKNOWN", POSITION="APPEND", FORM="FORMATTED")

WRITE(21,*) "BLOCK", nblock

DO nmp= 1, nummpenergies
  WRITE(21,"(I3,3X)",ADVANCE="NO") nmp
  WRITE(21,*) mpenergies(nmp)

  IF (nmp <= nummpstates) THEN
    DO nsd= 1, dimblock
      weight= ABS(mpstates(nsd,nmp))**2
      IF (weight<1e-6) CYCLE
      kete= ket( nsd, 1 )
      keth= ket( nsd, 2 )

      WRITE(21,"(E9.3,15X)",ADVANCE="NO") mpstates(nsd,nmp)
      WRITE(21,"(A,1X)",ADVANCE="NO") "e"
      WRITE(21,binfmt_e) kete
      WRITE(21,"('_______',F6.3,'_%','________ ')",ADVANCE="NO") weight*100
      WRITE(21,"(A,1X)",ADVANCE="NO") "h"
      WRITE(21,binfmt_h) keth

    END DO
  END IF
  WRITE(21,*) 
END DO

CLOSE( 21 )

END SUBROUTINE WRITEMPSTATES


END MODULE mod_writempstates
