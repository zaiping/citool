MODULE mod_specialf
  IMPLICIT NONE
CONTAINS


!!$!**********************************************************************
!!$FUNCTION BESSELJ(nn,xx)
!!$!  USE IFPORT  ! for intel compiler
!!$  IMPLICIT NONE
!!$!  retuns the value of the Bessel function J using the implicit function
!!$!  dbesjn(n,x), for double precision in/output (PGI compiler)
!!$  INTEGER, INTENT(IN) :: nn
!!$  REAL(8), INTENT(IN) :: xx
!!$  REAL(8) :: BESSELJ
!!$!.................................declares the implicit function
!!$  REAL(8) :: dbesjn
!!$!....................................calls the implicit function
!!$  BESSELJ= dbesjn(nn, xx)
!!$END FUNCTION BESSELJ


!**********************************************************************
FUNCTION LAGUERRE(mm,nn,xx)
  IMPLICIT NONE
!  retuns the value of the Laguerre series L^m_n(x)
  INTEGER, INTENT(IN) :: mm, nn
  REAL(8), INTENT(IN) :: xx
  REAL(8) :: LAGUERRE
!.........................................................Tipi locali
  INTEGER :: ni
  REAL(8) :: lmn, den
  IF (nn < 0) STOP "ERRORE: (nn < 0)  in LAGUERRE"
!....................................computes the series
  lmn= 0.
  DO ni = 0, nn
    den= FACTORIAL(nn-ni) * FACTORIAL(mm+ni) * FACTORIAL(ni)
    lmn= lmn + ( FACTORIAL(nn+mm) *  (-xx)**ni / den )
  END DO
  LAGUERRE= lmn
END FUNCTION LAGUERRE

!**********************************************************************
FUNCTION BINOMIALCO(nn, kk)
  IMPLICIT NONE
!  retuns the binomial coefficient
!  it can handle larger nn and kk than the formula  nn!/(kk!(nn-kk)!)
  INTEGER, INTENT(IN) :: nn, kk
  INTEGER :: BINOMIALCO
!.........................................................Tipi locali
  INTEGER :: ii, k1
  REAL*8 :: prod

  IF (nn < 0 .OR. kk < 0)  STOP "BINOMIALCO: negative arg"
  IF (nn < kk)  STOP "BINOMIALCO: nn < kk"

  k1 = MIN(kk, nn-kk)
  prod= 1.
  DO ii= 1, k1
    prod= prod * REAL(nn-ii+1)/ii
  END DO
  BINOMIALCO= NINT(prod)

END FUNCTION BINOMIALCO

!**********************************************************************
FUNCTION FACTORIAL(nn)
  IMPLICIT NONE
!  retuns the factorial of the input:  nn!
  INTEGER, INTENT(IN) :: nn
  INTEGER :: FACTORIAL
!.........................................................Tipi locali
  INTEGER :: n1, nprod
  IF (nn < 0) STOP "ERRORE: (nn < 0)  in FACTORIAL"
  nprod= 1
  DO n1= 1, nn
    nprod= nprod * n1
  END DO
  FACTORIAL= nprod
END FUNCTION FACTORIAL


SUBROUTINE HERMITE_POLY( n, x, cx )
! from polpak.f90   http://www.csit.fsu.edu/~burkardt/f_src/polpak/polpak.html
!*******************************************************************************
!
!! HERMITE_POLY evaluates the Hermite polynomials at X.
!
!  Differential equation:
!
!    Y'' - 2 X Y' + 2 N Y = 0
!
!  First terms:
!
!      1
!      2 X
!      4 X**2     -  2
!      8 X**3     - 12 X
!     16 X**4     - 48 X**2     + 12
!     32 X**5    - 160 X**3    + 120 X
!     64 X**6    - 480 X**4    + 720 X**2    - 120
!    128 X**7   - 1344 X**5   + 3360 X**3   - 1680 X
!    256 X**8   - 3584 X**6  + 13440 X**4  - 13440 X**2   + 1680
!    512 X**9   - 9216 X**7  + 48384 X**5  - 80640 X**3  + 30240 X
!   1024 X**10 - 23040 X**8 + 161280 X**6 - 403200 X**4 + 302400 X**2 - 30240
!
!  Recursion:
!
!    H(0,X) = 1,
!    H(1,X) = 2*X,
!    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
!
!  Norm:
!
!    Integral ( -Infinity < X < Infinity ) exp ( - X**2 ) * H(N,X)**2 dX
!    = sqrt ( PI ) * 2**N * N!
!
!    H(N,X) = (-1)**N * exp ( X**2 ) * dn/dXn ( exp(-X**2 ) )
!
!  Modified:
!
!    01 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964.
!
!    Larry Andrews,
!    Special Functions of Mathematics for Engineers,
!    Second Edition, 
!    Oxford University Press, 1998.
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials are 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 Hermite
!    polynomials at the point X.
!
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL*8,  INTENT(IN) :: x
  REAL*8,  INTENT(OUT):: cx(0:n)

  INTEGER :: ii

  IF ( n < 0 )  RETURN

  cx(0) = 1.0D+00

  IF ( n == 0 )  RETURN

  cx(1) = 2.0D+00 * x
 
  DO ii = 2, n
    cx(ii) = (2.0D+00 * x * cx(ii-1)) - (2.0D+00 * REAL(ii-1) * cx(ii-2))
  END DO
 
END SUBROUTINE HERMITE_POLY


!**********************************************************************
END MODULE mod_specialf
