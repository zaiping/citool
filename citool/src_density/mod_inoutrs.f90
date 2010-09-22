
MODULE mod_inoutrs
IMPLICIT NONE
SAVE
! for HEXAGONAL domain and grid
! reads the BINARY file for wave functions

! INPUT params and files
CHARACTER(80), PARAMETER :: FILEwavefunction_e= "psi_e.bin"
CHARACTER(80), PARAMETER :: FILEwavefunction_h= ""
! OUTPUT params and files
CHARACTER(80), PARAMETER :: FILEdensTOTe= "densTOTe.hdat"
CHARACTER(80), PARAMETER :: FILEdensUPe= "densUPe.hdat"
CHARACTER(80), PARAMETER :: FILEdensDNe= "densDNe.hdat"
! vars used in this module alone (values inside the psi file)
INTEGER :: numh_inoutrs
REAL*8 :: dh_inoutrs

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE INSPWF (numspwf, numx, psi, filename)
IMPLICIT NONE

INTEGER, INTENT(IN) :: numspwf
INTEGER, INTENT(OUT) :: numx
REAL*8, ALLOCATABLE, INTENT(OUT) :: psi(:,:)
CHARACTER(*), INTENT(IN) :: filename

REAL*8, ALLOCATABLE :: psihex(:,:,:)
INTEGER :: numpsi_read
INTEGER :: np, nq, nx

!.........................................reads single-particle wave functions
OPEN(33, FILE=filename, FORM="UNFORMATTED", ACTION="READ",  &
       &   STATUS="OLD")
READ(33) numh_inoutrs, numpsi_read, dh_inoutrs
IF (2*numpsi_read /= numspwf) STOP "INSPWF: wrong numpsi in psi binary file"
ALLOCATE(psihex(-numh_inoutrs:numh_inoutrs,-numh_inoutrs:numh_inoutrs,numpsi_read))
READ(33) psihex
CLOSE(33)

numx= 3*(numh_inoutrs+1)*numh_inoutrs + 1
ALLOCATE(psi(numx,numspwf))
nx=0
DO nq= -numh_inoutrs, numh_inoutrs
  DO np= -numh_inoutrs, numh_inoutrs
    IF (np <= nq+numh_inoutrs .AND. np >= nq-numh_inoutrs) THEN
      nx= nx + 1
      psi(nx,1:numspwf/2)= psihex(np,nq,1:numpsi_read)
      psi(nx,numspwf/2+1:numspwf)= psihex(np,nq,1:numpsi_read)
    END IF
  END DO
END DO
IF (nx /= numx) STOP "INSPWF: nx /= num nodes"

DEALLOCATE(psihex)

END SUBROUTINE INSPWF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE OUTDENS (numx, dens, filename)
IMPLICIT NONE
INTEGER, INTENT(IN) :: numx
REAL*8, INTENT(IN) :: dens(numx) 
CHARACTER(*), INTENT(IN) :: filename

INTEGER :: np, nq, nx

IF (numx /= 3*(numh_inoutrs+1)*numh_inoutrs + 1) STOP "OUTDENS: numx /= num nodes"

OPEN(15, FILE=filename, FORM="FORMATTED")
WRITE(15,*) numh_inoutrs, dh_inoutrs, numx  ! number of nodes

nx= 0
DO nq= -numh_inoutrs, numh_inoutrs
  DO np= -numh_inoutrs, numh_inoutrs
    IF (np <= nq+numh_inoutrs .AND. np >= nq-numh_inoutrs) THEN
      nx= nx + 1
      WRITE(15,*) np, nq, dens(nx)
    END IF
  END DO
  WRITE(15,*)
END DO
IF (nx /= numx) STOP "OUTDENS: nx /= num nodes"

CLOSE(15)

END SUBROUTINE OUTDENS


END MODULE mod_inoutrs

