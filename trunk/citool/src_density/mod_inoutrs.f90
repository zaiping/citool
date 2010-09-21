
MODULE mod_inoutrs
IMPLICIT NONE
SAVE

! INPUT params and files
INTEGER, PARAMETER :: NUMRin= 81, NUMZin=961
CHARACTER(80), PARAMETER :: FILENAMEwavefunction_e= "wf_all_mz.dat"
CHARACTER(80), PARAMETER :: FILENAMEwavefunction_h= ""
! OUTPUT params and files
INTEGER, PARAMETER :: NUMRout= 81, NUMZout=961
CHARACTER(80), PARAMETER :: FILENAMEdensTOTe= "densTOTe.dat"
CHARACTER(80), PARAMETER :: FILENAMEdensUPe= "densUPe.dat"
CHARACTER(80), PARAMETER :: FILENAMEdensDNe= "densDNe.dat"

CONTAINS

SUBROUTINE INSPWF (numspwf, numx, psi, filename)
IMPLICIT NONE

INTEGER, INTENT(IN) :: numspwf
INTEGER, INTENT(OUT) :: numx
REAL*8, ALLOCATABLE, INTENT(OUT) :: psi(:,:)
CHARACTER(*), INTENT(IN) :: filename

INTEGER :: i, j, k, nnr, nnz
INTEGER :: n

numx= NUMRin*NUMZin

ALLOCATE(psi(numx,numspwf))

write(*,*) "numspwf, NUMR, NUMZ", numspwf, NUMRin, NUMZin

OPEN(UNIT=8, FILE=filename)
read (8,*)
read (8,*)
do k = 1, numspwf/2
  n= 0
  do i = 1, NUMRin
    do j = 1, NUMZin
      n = n + 1
      read (8,*) nnr, nnz, psi(n, k)
      !write (*,*) psi(n, k)
      if (i /= nnr) STOP "warning..."
    end do
    ! write(*,*) "i=", i
    ! write(*,*) "psi=", psi(n, k)
  end do
  read(8,*)
  write(*,*) "k=", k
end do

psi(:, numspwf/2+1:numspwf)= psi(:, 1:numspwf/2)

CLOSE(8)

END SUBROUTINE INSPWF



SUBROUTINE OUTDENS (numx, dens, filename)
IMPLICIT NONE
INTEGER, INTENT(IN) :: numx
REAL*8, INTENT(IN) :: dens(numx) 
CHARACTER(*), INTENT(IN) :: filename

INTEGER :: i, j, k
INTEGER :: n

if (numx /= NUMRout*NUMZout) STOP "error here"

OPEN(UNIT=15, FILE=filename)

n= 0
do i = 1, NUMRout
  do j = 1, NUMZout
     n = n + 1
     write (15,*) i, j, dens(n)
  end do
  write (15,*)
end do

CLOSE(15)

END SUBROUTINE OUTDENS

END MODULE mod_inoutrs

