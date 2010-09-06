MODULE mod_indexx
IMPLICIT  NONE
  SAVE

CONTAINS

SUBROUTINE indexx(nelem, arr, indx)
! Indexes an array arr(1:n), i.e., outputs the array indx(1:n)
! such that arr(indx(j)) is in ascending order for j = 1, 2, . . . , N .
! The input quantities n and arr are not changed.
! from NR cap 8
! code taken from
!   www.phys.uu.nl/DU/num_recipes/fortran.208/f77/recipes/indexx.for
IMPLICIT NONE
INTEGER, INTENT(IN) :: nelem      ! elements in arr and indx
REAL*8, INTENT(IN) :: arr(nelem) 
INTEGER, INTENT(OUT) :: indx(nelem)

INTEGER, PARAMETER :: M= 7 , NSTACK= 50

INTEGER i, indxt, ir, itemp, j, jstack, k, l, istack(NSTACK)
REAL a
DO j= 1, nelem
  indx(j)= j
END DO
jstack= 0
l= 1
ir= nelem

1 IF (ir-l .LT. M) THEN
  DO j= l+1, ir
    indxt= indx(j)
    a= arr(indxt)
    DO i= j-1, l, -1
      IF (arr(indx(i)) .LE. a) GOTO 2
      indx(i+1)= indx(i)
    END DO
    i= l-1
2   indx(i+1)= indxt
  END DO
  IF (jstack .EQ. 0) RETURN
  ir= istack(jstack)
  l= istack(jstack-1)
  jstack= jstack-2
ELSE
  k= (l+ir)/2
  itemp= indx(k)
  indx(k)= indx(l+1)
  indx(l+1)= itemp
  IF (arr(indx(l)) .GT. arr(indx(ir))) THEN
    itemp= indx(l)
    indx(l)= indx(ir)
    indx(ir)= itemp
  ENDIF
  IF (arr(indx(l+1)) .GT. arr(indx(ir))) THEN
    itemp= indx(l+1)
    indx(l+1)= indx(ir)
    indx(ir)= itemp
  ENDIF
  IF (arr(indx(l)) .GT. arr(indx(l+1))) THEN
    itemp= indx(l)
    indx(l)= indx(l+1)
    indx(l+1)= itemp
  ENDIF
  i= l+1
  j= ir
  indxt= indx(l+1)
  a= arr(indxt)
3 CONTINUE
  i= i+1
  IF (arr(indx(i)) .LT. a) GOTO 3
4 CONTINUE
  j=j-1
  IF (arr(indx(j)) .GT. a) GOTO 4
  IF (j .LT. i) GOTO 5
  itemp= indx(i)
  indx(i)= indx(j)
  indx(j)= itemp
  GOTO 3
5 indx(l+1)= indx(j)
  indx(j)= indxt
  jstack= jstack + 2
  IF (jstack .GT. NSTACK) PAUSE 'NSTACK too small in indexx'
  IF (ir-i+1 .GE. j-l) THEN
    istack(jstack)= ir
    istack(jstack-1)= i
    ir= j-1
  ELSE
    istack(jstack)= j-1
    istack(jstack-1)= l
    l= i
  ENDIF
ENDIF
GOTO 1

END SUBROUTINE indexx


END MODULE mod_indexx
