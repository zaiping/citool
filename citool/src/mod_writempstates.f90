!$$$$$$$$$$$$$$$$$$$$$$$$$  CItool  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_writempstates
  IMPLICIT  NONE
  SAVE

CONTAINS

SUBROUTINE WRITEMPSTATES_X( numblock, dimhspacecons, blockstart, ket,           &
     &  nummpenergiescons, mpenergies, blocknummpenergies,                      &
     &  nummpstatescons, mpstates_x, blocknummpstates )
  USE mod_indata
  IMPLICIT NONE
! 
INTEGER, INTENT(IN) :: numblock
INTEGER, INTENT(IN) :: dimhspacecons
INTEGER, INTENT(IN) :: blockstart(dimhspacecons+1)
INTEGER*8, INTENT(IN) :: ket(dimhspacecons,2)
INTEGER, INTENT(IN) :: nummpenergiescons
REAL*8, INTENT(IN) :: mpenergies(nummpenergiescons,numblock)
INTEGER, INTENT(IN) :: blocknummpenergies(dimhspacecons)
INTEGER, INTENT(IN) :: nummpstatescons
COMPLEX*16, INTENT(IN) :: mpstates_x(dimhspacecons,nummpstatescons)
INTEGER, INTENT(IN) :: blocknummpstates(dimhspacecons)

REAL*8 :: weight
INTEGER*8 :: kete, keth
INTEGER :: nb, blockdim, nmp, nsd


IF (fileoutASC_mpstates /= "") THEN
  OPEN(21, FILE=TRIM(fileoutASC_mpstates),                    &
     &   ACTION="WRITE", STATUS="OLD", POSITION="APPEND", FORM="FORMATTED")

  DO nb= 1, numblock 
    blockdim= blockstart(nb+1)-blockstart(nb)

    WRITE(21,"(A6,I4,A9,I8)") "BLOCK ", nb, " ,  dim: ", blockdim

    DO nmp= 1, blocknummpenergies(nb)
      WRITE(21,"(I3,3X)",ADVANCE="NO") nmp
      WRITE(21,*) mpenergies(nmp,nb)
      IF (nmp <= blocknummpstates(nb)) THEN
        DO nsd= blockstart(nb), blockstart(nb+1)-1
          weight= ABS(mpstates_x(nsd,nmp))**2
          IF (weight < cutoff_fileoutASC_mpstates) CYCLE
          kete= ket( nsd, 1 )
          keth= ket( nsd, 2 )
          WRITE(21,"('(',E9.3,',',E9.3,')',3X)",ADVANCE="NO") mpstates_x(nsd,nmp)
          WRITE(21,"(A,1X)",ADVANCE="NO") "e"
          WRITE(21,binfmt_e) kete
          WRITE(21,"('_______',F7.3,'_%','_______ ')",ADVANCE="NO") weight*100
          WRITE(21,"(A,1X)",ADVANCE="NO") "h"
          WRITE(21,binfmt_h) keth
        END DO
      END IF
      WRITE(21,*) 
    END DO

  END DO

  CLOSE(21)

END IF

IF (fileoutBIN_mpstates /= "") THEN
  OPEN(22, FILE=TRIM(fileoutBIN_mpstates),                    &
     &   ACTION="WRITE", STATUS="OLD", POSITION="APPEND", FORM="UNFORMATTED")
  WRITE(22) mpenergies
  WRITE(22) blocknummpenergies
  WRITE(22) mpstates_x
  WRITE(22) blocknummpstates
  CLOSE(22)
END IF

END SUBROUTINE WRITEMPSTATES_X


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

SUBROUTINE WRITEMPSTATES( numblock, dimhspacecons, blockstart, ket,           &
     &  nummpenergiescons, mpenergies, blocknummpenergies,                    &
     &  nummpstatescons, mpstates, blocknummpstates )
  USE mod_indata
  IMPLICIT NONE
! 
INTEGER, INTENT(IN) :: numblock
INTEGER, INTENT(IN) :: dimhspacecons
INTEGER, INTENT(IN) :: blockstart(dimhspacecons+1)
INTEGER*8, INTENT(IN) :: ket(dimhspacecons,2)
INTEGER, INTENT(IN) :: nummpenergiescons
REAL*8, INTENT(IN) :: mpenergies(nummpenergiescons,numblock)
INTEGER, INTENT(IN) :: blocknummpenergies(dimhspacecons)
INTEGER, INTENT(IN) :: nummpstatescons
REAL*8, INTENT(IN) :: mpstates(dimhspacecons,nummpstatescons)
INTEGER, INTENT(IN) :: blocknummpstates(dimhspacecons)

REAL*8 :: weight
INTEGER*8 :: kete, keth
INTEGER :: nb, blockdim, nmp, nsd


IF (fileoutASC_mpstates /= "") THEN
  OPEN(21, FILE=TRIM(fileoutASC_mpstates),                    &
       &   ACTION="WRITE", STATUS="OLD", POSITION="APPEND", FORM="FORMATTED")

  DO nb= 1, numblock 
    blockdim= blockstart(nb+1)-blockstart(nb)

    WRITE(21,"(A6,I4,A9,I8)") "BLOCK ", nb, " ,  dim: ", blockdim

    DO nmp= 1, blocknummpenergies(nb)
      WRITE(21,"(I3,3X)",ADVANCE="NO") nmp
      WRITE(21,*) mpenergies(nmp,nb)
      IF (nmp <= blocknummpstates(nb)) THEN
        DO nsd= blockstart(nb), blockstart(nb+1)-1
          weight= ABS(mpstates(nsd,nmp))**2
          IF (weight < cutoff_fileoutASC_mpstates) CYCLE
          kete= ket( nsd, 1 )
          keth= ket( nsd, 2 )
          WRITE(21,"(E9.3,15X)",ADVANCE="NO") mpstates(nsd,nmp)
          WRITE(21,"(A,1X)",ADVANCE="NO") "e"
          WRITE(21,binfmt_e) kete
          WRITE(21,"('_______',F7.3,'_%','_______ ')",ADVANCE="NO") weight*100
          WRITE(21,"(A,1X)",ADVANCE="NO") "h"
          WRITE(21,binfmt_h) keth
        END DO
      END IF
      WRITE(21,*) 
    END DO

  END DO

  CLOSE(21)

END IF

IF (fileoutBIN_mpstates /= "") THEN
  OPEN(22, FILE=TRIM(fileoutBIN_mpstates),                    &
     &   ACTION="WRITE", STATUS="OLD", POSITION="APPEND", FORM="UNFORMATTED")
  WRITE(22) mpenergies
  WRITE(22) blocknummpenergies
  WRITE(22) mpstates
  WRITE(22) blocknummpstates
  CLOSE(22)
END IF

END SUBROUTINE WRITEMPSTATES


END MODULE mod_writempstates
