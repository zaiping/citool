
MODULE mod_inoutrs
IMPLICIT NONE
SAVE

CONTAINS

SUBROUTINE INSPWF (numspwf, numx, psi, filename)
IMPLICIT NONE

INTEGER, INTENT(IN) :: numspwf
INTEGER, INTENT(OUT) :: numx
REAL*8, ALLOCATABLE, INTENT(OUT) :: psi(:,:)
CHARACTER(*), INTENT(IN) :: filename

INTEGER :: i, j, k, nnr, nnz
INTEGER :: n

INTEGER, PARAMETER :: NUMR= 81, NUMZ=961

numx= NUMR*NUMZ

ALLOCATE(psi(numx,numspwf))

write(*,*) numspwf, NUMR, NUMZ

OPEN(UNIT=8, FILE=filename)
read (8,*)
read (8,*)
do k = 1, numspwf/2
  n= 0
  do i = 1, NUMR
    do j = 1, NUMZ
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

write(*,*) NUMR, NUMZ, psi(n, numspwf)

CLOSE(8)

END SUBROUTINE INSPWF



SUBROUTINE OUTDENS (numx, dens, filename)
IMPLICIT NONE
INTEGER, INTENT(IN) :: numx
REAL*8, INTENT(IN) :: dens(numx) 
CHARACTER(*), INTENT(IN) :: filename

INTEGER :: i, j, k
INTEGER :: n

INTEGER, PARAMETER :: NUMR= 81, NUMZ=961

if (numx /= NUMR*NUMZ) STOP "error here"


OPEN(UNIT=15, FILE=filename)

n= 0
do i = 1, NUMR
  do j = 1, NUMZ
     n = n + 1
     write (15,*) i, j, dens(n)
  end do
  write (15,*)
end do

CLOSE(15)

END SUBROUTINE OUTDENS

END MODULE mod_inoutrs

