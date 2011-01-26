
MODULE mod_inoutrs
!USE mod_indatadensity
IMPLICIT NONE
SAVE
! for HEXAGONAL domain and grid
! reads the BINARY file for wave functions


!INTEGER :: WANTCONSE= 1 , WANTCONSH= 1
!INTEGER :: WANTBLOCK= 0 , WANTRANK= 0

! INPUT params and files
!CHARACTER(80), PARAMETER :: FILEwavefunction_e= "psi_e.bin"
!CHARACTER(80), PARAMETER :: FILEwavefunction_h= ""
! OUTPUT params and files
!CHARACTER(80), PARAMETER :: FILEdensTOTe= "densTOTe.hdat"
!CHARACTER(80), PARAMETER :: FILEdensUPe= "densUPe.hdat"
!CHARACTER(80), PARAMETER :: FILEdensDNe= "densDNe.hdat"


! vars used in this module alone (values inside the psi file)
INTEGER, PRIVATE :: numh_inoutrs
REAL*8, PRIVATE :: dh_inoutrs

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE INSPWF (numspwf, numx, psi, filename)
IMPLICIT NONE

INTEGER, INTENT(OUT) :: numspwf
INTEGER, INTENT(OUT) :: numx
REAL*8, ALLOCATABLE, INTENT(OUT) :: psi(:,:)
CHARACTER(*), INTENT(IN) :: filename

REAL*8, ALLOCATABLE :: psihex(:,:,:)
REAL*8 :: constnorm
INTEGER :: np, nq, nx

!.........................................reads single-particle wave functions
OPEN(33, FILE=filename, FORM="UNFORMATTED", ACTION="READ",  &
       &   STATUS="OLD")
READ(33) numh_inoutrs, numspwf, dh_inoutrs
ALLOCATE(psihex(-numh_inoutrs:numh_inoutrs,-numh_inoutrs:numh_inoutrs,numspwf))
READ(33) psihex
CLOSE(33)

constnorm= SQRT( SQRT(3./4.)*dh_inoutrs*dh_inoutrs )
numx= 3*(numh_inoutrs+1)*numh_inoutrs + 1
ALLOCATE(psi(numx,numspwf))
nx=0
DO nq= -numh_inoutrs, numh_inoutrs
  DO np= -numh_inoutrs, numh_inoutrs
    IF (np <= nq+numh_inoutrs .AND. np >= nq-numh_inoutrs) THEN
      nx= nx + 1
      psi(nx,:)= constnorm * psihex(np,nq,:)
    END IF
  END DO
END DO
IF (nx /= numx) STOP "INSPWF: nx /= num nodes"

DEALLOCATE(psihex)

END SUBROUTINE INSPWF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE OUTDENS (numx, dens, filename, denssum)
IMPLICIT NONE
INTEGER, INTENT(IN) :: numx
REAL*8, INTENT(IN) :: dens(numx) 
CHARACTER(*), INTENT(IN) :: filename
REAL*8, INTENT(OUT) :: denssum

INTEGER :: np, nq, nx

IF (numx /= 3*(numh_inoutrs+1)*numh_inoutrs + 1) STOP "OUTDENS: numx /= num nodes"

OPEN(15, FILE=filename, FORM="FORMATTED")
WRITE(15,*) numh_inoutrs, dh_inoutrs, numx  ! number of nodes

denssum= 0.
nx= 0
DO nq= -numh_inoutrs, numh_inoutrs
  DO np= -numh_inoutrs, numh_inoutrs
    IF (np <= nq+numh_inoutrs .AND. np >= nq-numh_inoutrs) THEN
      nx= nx + 1
      WRITE(15,*) np, nq, dens(nx)
      denssum= denssum + dens(nx)
    END IF
  END DO
  WRITE(15,*)
END DO
IF (nx /= numx) STOP "OUTDENS: nx /= num nodes"

CLOSE(15)

END SUBROUTINE OUTDENS


END MODULE mod_inoutrs

