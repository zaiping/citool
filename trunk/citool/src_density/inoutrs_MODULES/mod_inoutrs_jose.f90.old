
MODULE mod_inoutrs
IMPLICIT NONE
SAVE
! for nanodumbbells

INTEGER :: WANTBLOCK= 0 , WANTRANK= 0

! INPUT params and files
INTEGER, PARAMETER :: NUMRin= 81, NUMZin=961
CHARACTER(80), PARAMETER :: FILEwavefunction_e= "wf_all_mz.dat"
CHARACTER(80), PARAMETER :: FILEwavefunction_h= ""
! OUTPUT params and files
INTEGER, PARAMETER :: NUMRout= 81, NUMZout=961
CHARACTER(80), PARAMETER :: FILEdensTOTe= "densTOTe.dat"
CHARACTER(80), PARAMETER :: FILEdensUPe= "densUPe.dat"
CHARACTER(80), PARAMETER :: FILEdensDNe= "densDNe.dat"

CONTAINS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE INSPWF (numspwf, numx, psi, filename)
IMPLICIT NONE

INTEGER, INTENT(IN) :: numspwf
INTEGER, INTENT(OUT) :: numx
REAL*8, ALLOCATABLE, INTENT(OUT) :: psi(:,:)
CHARACTER(*), INTENT(IN) :: filename

INTEGER :: ni, nj, nk, nnr, nnz
INTEGER :: nn

numx= NUMRin*NUMZin

ALLOCATE(psi(numx,numspwf))

!write(*,*) "numspwf, NUMR, NUMZ", numspwf, NUMRin, NUMZin

OPEN(UNIT=8, FILE=filename)
read (8,*)
read (8,*)
do nk = 1, numspwf/2
  nn= 0
  do ni = 1, NUMRin
    do nj = 1, NUMZin
      nn = nn + 1
      read (8,*) nnr, nnz, psi(nn, nk)
      !write (*,*) psi(n, k)
      if (ni /= nnr) STOP "warning..."
    end do
    ! write(*,*) "i=", i
    ! write(*,*) "psi=", psi(n, k)
  end do
  read(8,*)
!  write(*,*) "k=", nk
end do

psi(:, numspwf/2+1:numspwf)= psi(:, 1:numspwf/2)

CLOSE(8)

END SUBROUTINE INSPWF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE OUTDENS (numx, dens, filename, denssum)
IMPLICIT NONE
INTEGER, INTENT(IN) :: numx
REAL*8, INTENT(IN) :: dens(numx) 
CHARACTER(*), INTENT(IN) :: filename
REAL*8, INTENT(OUT) :: denssum

INTEGER :: ni, nj, nk
INTEGER :: nn

if (numx /= NUMRout*NUMZout) STOP "error here"

OPEN(UNIT=15, FILE=filename)

denssum= 0.
nn= 0
do ni = 1, NUMRout
  do nj = 1, NUMZout
     nn = nn + 1
     write (15,*) ni, nj, dens(nn)
     denssum= denssum + dens(nn)
  end do
  write (15,*)
end do

CLOSE(15)

END SUBROUTINE OUTDENS

END MODULE mod_inoutrs

