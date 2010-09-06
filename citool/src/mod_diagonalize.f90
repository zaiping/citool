!$$$$$$$$$$$$$$$$$$$$$$$$$  CItool  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_diagonalize
  IMPLICIT  NONE
  SAVE

CONTAINS

SUBROUTINE DIAGONALIZE_X( dimham, hami, hamj, ham,          &
     &                    numevals, evals, numevecs, evecs )
  USE mod_indexx
  IMPLICIT NONE
! finds eigenvalues and eigenvectors of the COMPLEX Hamiltonian
! uses ARPACK routines for the diagonalization
! H stored using the Compressed Row Storage format (same as Pardiso)
! see  http://www.cs.utk.edu/~dongarra/etemplates/node371.html
! since H is hermitian, we store only the diag and the upper triangle

INTEGER, INTENT(IN) :: dimham            ! dimension of the H matrix
INTEGER, INTENT(IN) :: hami(:)
INTEGER, INTENT(IN) :: hamj(:)
COMPLEX*16, INTENT(IN) :: ham(:)
INTEGER, INTENT(IN) :: numevals
REAL*8, INTENT(OUT) :: evals(1:)  !(numevals)
INTEGER, INTENT(IN) :: numevecs
COMPLEX*16, INTENT(OUT) :: evecs(1:,1:)  !(dimham,numevecs)

! znaupd args
INTEGER :: ido_ap, n_ap, nev_ap
REAL*8 :: tol_ap
COMPLEX*16, ALLOCATABLE :: resid_ap(:)
INTEGER :: ncv_ap
COMPLEX*16, ALLOCATABLE :: v_ap(:,:)
INTEGER :: iparam_ap(11), ipntr_ap(11)
COMPLEX*16, ALLOCATABLE :: workd_ap(:)
COMPLEX*16, ALLOCATABLE :: workl_ap(:)
INTEGER :: lworkl_ap
REAL*8, ALLOCATABLE :: rwork_ap(:)
INTEGER :: info_ap

! zneupd args
LOGICAL :: rvec_ap
LOGICAL, ALLOCATABLE :: select_ap(:)
COMPLEX*16, ALLOCATABLE :: d_ap(:)
COMPLEX*16 :: sigma_ap
COMPLEX*16, ALLOCATABLE :: workev_ap(:)

! args for computing residuals
INTEGER :: nconv                         ! num of converged eigenvalues
COMPLEX*16, ALLOCATABLE :: ax_ap(:)
REAL*8, ALLOCATABLE :: rd_ap(:,:)
REAL*8, EXTERNAL :: dznrm2, dlapy2

! other vars
REAL*8 :: norm
INTEGER, ALLOCATABLE :: indexd(:)
INTEGER :: nn, kk

IF (numevals>dimham .OR. numevecs>dimham) STOP "DIAGONALIZE: numev>dimham"
IF (numevals<=0 .OR. numevecs<0) STOP "DIAGONALIZE: numev<=0"

IF (dimham==1) THEN
  evals(1)= ABS(ham(1))
  IF (numevecs > 0) THEN
    evecs(1,1)= (1.,0.)
  END IF
  RETURN
END IF

!................................................initializing solver data
ido_ap= 0
n_ap= dimham
nev_ap= numevals  ! number of eigenvalues to be computed
tol_ap= -1.
ncv_ap= MIN(4*(nev_ap+1),n_ap)  ! should be > 2* nev_ap
lworkl_ap= 3*ncv_ap**2 + 5*ncv_ap
info_ap= 0
ipntr_ap(:)= 0
iparam_ap(:)= 0
iparam_ap(1)= 1     ! uses the exact shift strategy
iparam_ap(3)= 600   ! on INPUT: max num Arnoldi update iterations allowed
iparam_ap(4)= 1
iparam_ap(7)= 1     ! mode
rvec_ap= .TRUE.

!................................................allocating solver arrays
ALLOCATE( resid_ap(n_ap) )
ALLOCATE( v_ap(n_ap,ncv_ap) )
ALLOCATE( workd_ap(3*n_ap) )
ALLOCATE( workl_ap(lworkl_ap) )
ALLOCATE( rwork_ap(ncv_ap) )
ALLOCATE( select_ap(ncv_ap) )
ALLOCATE( d_ap(nev_ap+1) )
ALLOCATE( workev_ap(2*ncv_ap) )
ALLOCATE( ax_ap(n_ap) )
ALLOCATE( rd_ap(nev_ap,3) )

ALLOCATE( indexd(numevals) )

!.........................................ARPACK reverse communication loop
DO

  CALL znaupd ( ido_ap, 'I', n_ap, 'SR', nev_ap, tol_ap, resid_ap,    &
       &        ncv_ap, v_ap, n_ap, iparam_ap, ipntr_ap,              &
       &        workd_ap, workl_ap, lworkl_ap, rwork_ap, info_ap )

  IF (ido_ap == 99) EXIT
  IF (ido_ap /= 1 .AND. ido_ap /= -1) STOP "ido_ap /=+-1"

  CALL avmult_x( n_ap, hami, hamj, ham,                        &
       &         workd_ap(ipntr_ap(1)), workd_ap(ipntr_ap(2)) )

  !PRINT*, "ARPACK loop - converged Ritz values: ", iparam_ap(5)

END DO

IF (info_ap < 0 .OR. info_ap > 1) THEN
  PRINT*, "DIAGONALIZE: in ARPACK loop znaupd, info_ap=", info_ap
  PRINT*, n_ap, nev_ap, ncv_ap
  STOP
ELSE IF (info_ap == 1) THEN
  PRINT*, "ARPACK WARNING in DIAGONALIZE: in loop znaupd"
  PRINT*, "info_ap, iparam_ap(5)  =", info_ap, iparam_ap(5)
END IF
!PRINT*, "Ok, number of ARPACK Arnoldi iterations taken: ", iparam_ap(3)

!................................................ARPACK results extraction
CALL zneupd ( rvec_ap, 'A', select_ap, d_ap, v_ap, n_ap, sigma_ap,       &
     &        workev_ap,                                                 &
     &        'I', n_ap, 'SR', nev_ap, tol_ap, resid_ap, ncv_ap, v_ap,   &
     &        n_ap, iparam_ap, ipntr_ap, workd_ap, workl_ap, lworkl_ap,  &
     &        rwork_ap, info_ap )

!PRINT*, "iparam_ap(5), nev_ap, ncv_ap:"
!PRINT*,  iparam_ap(5), nev_ap, ncv_ap
!PRINT*, "ipntr_ap=", ipntr_ap
!PRINT*

IF (info_ap /= 0) THEN
  PRINT*, "DIAGONALIZE: in zneupd, info_ap=", info_ap
  STOP
END IF

!...................................computes and displays the residual norm
nconv= iparam_ap(5)
IF (nconv /= numevals) THEN
  PRINT*, "ARPACK WARNING ! num of computed evals with desired tol:", nconv
  STOP
END IF
DO nn= 1, nconv
  CALL avmult_x( n_ap, hami, hamj, ham,   &
       &         v_ap(:,nn), ax_ap )
  CALL zaxpy(n_ap, -d_ap(nn), v_ap(:,nn), 1, ax_ap, 1)
  rd_ap(nn,1) = REAL(d_ap(nn))
  rd_ap(nn,2) = AIMAG(d_ap(nn))
  rd_ap(nn,3) = dznrm2(n_ap, ax_ap, 1)
  rd_ap(nn,3) = rd_ap(nn,3) / dlapy2(rd_ap(nn,1),rd_ap(nn,2))
END DO
! uncomment to display infos on residual norm
!CALL dmout(6, nconv, 3, rd_ap, nev_ap, -6,                   &
!     &     'Ritz values (Real,Imag) and relative residuals')

!...............................fills evals and evecs in increasing order
CALL indexx(numevals, REAL(d_ap(1:numevals)), indexd(:))

DO nn= 1, numevals
  evals(nn)= REAL( d_ap(indexd(nn)) )
END DO

DO nn= 1, numevecs
  DO kk= 1, dimham
    evecs(kk,nn)= v_ap(kk,indexd(nn))
  END DO
  norm= SUM( ABS(evecs(:,nn))**2 )
  evecs(:,nn)= evecs(:,nn) / norm
END DO

!............................................................finalizations
DEALLOCATE( indexd )

DEALLOCATE( resid_ap )
DEALLOCATE( v_ap )
DEALLOCATE( workd_ap )
DEALLOCATE( workl_ap )
DEALLOCATE( rwork_ap )
DEALLOCATE( select_ap )
DEALLOCATE( d_ap )
DEALLOCATE( workev_ap )
DEALLOCATE( ax_ap)
DEALLOCATE( rd_ap )


CONTAINS

!************************************************************************
  SUBROUTINE avmult_x(dimat, imat, jmat, vmat, vap, wap)
    IMPLICIT NONE
    ! computes  w = M v
    ! for Hermitian matrices M in Compressed Row Storage format
    INTEGER,    INTENT(IN) :: dimat
    INTEGER,    INTENT(IN) :: imat(1:)  !(dimat+1)
    INTEGER,    INTENT(IN) :: jmat(1:)
    COMPLEX*16, INTENT(IN) :: vmat(1:)
    COMPLEX*16, INTENT(IN) :: vap(1:dimat)
    COMPLEX*16, INTENT(OUT) :: wap(1:dimat)
    
    INTEGER :: ni, nj

    wap(1:dimat)= (0.,0.)
    !
    DO ni= 1, dimat
      wap(ni) = wap(ni) + vmat(imat(ni)) * vap(jmat(imat(ni)))
      DO nj = imat(ni)+1, imat(ni+1) - 1
        wap(ni) = wap(ni) + vmat(nj) * vap(jmat(nj))
        wap(jmat(nj)) = wap(jmat(nj)) + CONJG(vmat(nj)) * vap(ni)
      END DO
    END DO

  END SUBROUTINE avmult_x
!************************************************************************

END SUBROUTINE DIAGONALIZE_X


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

SUBROUTINE DIAGONALIZE( dimham, hami, hamj, ham,          &
     &                  numevals, evals, numevecs, evecs )
  USE mod_indexx
  IMPLICIT NONE
! finds eigenvalues and eigenvectors of the REAL Hamiltonian
! uses ARPACK routines for the diagonalization
! H stored using the Compressed Row Storage format (same as Pardiso)
! see  http://www.cs.utk.edu/~dongarra/etemplates/node371.html
! since H is hermitian, we store only the diag and the upper triangle

INTEGER, INTENT(IN) :: dimham            ! dimension of the H matrix
INTEGER, INTENT(IN) :: hami(:)
INTEGER, INTENT(IN) :: hamj(:)
REAL*8,  INTENT(IN) :: ham(:)
INTEGER, INTENT(IN) :: numevals
REAL*8,  INTENT(OUT) :: evals(1:)  !(numevals)
INTEGER, INTENT(IN) :: numevecs
REAL*8,  INTENT(OUT) :: evecs(1:,1:)  !(dimham,numevecs)


! dsaupd args
INTEGER :: ido_ap, n_ap, nev_ap
REAL*8 :: tol_ap
REAL*8, ALLOCATABLE :: resid_ap(:)
INTEGER :: ncv_ap
REAL*8, ALLOCATABLE :: v_ap(:,:)
INTEGER :: iparam_ap(11), ipntr_ap(11)
REAL*8, ALLOCATABLE :: workd_ap(:)
REAL*8, ALLOCATABLE :: workl_ap(:)
INTEGER :: lworkl_ap
!rrr REAL*8, ALLOCATABLE :: rwork_ap(:)
INTEGER :: info_ap

! dseupd args
LOGICAL :: rvec_ap
LOGICAL, ALLOCATABLE :: select_ap(:)
REAL*8, ALLOCATABLE :: d_ap(:)
REAL*8, ALLOCATABLE :: d2_ap(:)   ! rrr
REAL*8 :: sigma_ap
!rrr COMPLEX*16, ALLOCATABLE :: workev_ap(:)

! args for computing residuals
INTEGER :: nconv                         ! num of converged eigenvalues
REAL*8, ALLOCATABLE :: ax_ap(:)
!rrr REAL*8, ALLOCATABLE :: rd_ap(:,:)
!rrr REAL*8, EXTERNAL :: dznrm2, dlapy2
REAL*8, EXTERNAL :: dnrm2  !rrr

! other vars
REAL*8 :: norm
INTEGER, ALLOCATABLE :: indexd(:)
INTEGER :: nn, kk

IF (numevals>dimham .OR. numevecs>dimham) STOP "DIAGONALIZE: numev>dimham"
IF (numevals<=0 .OR. numevecs<0) STOP "DIAGONALIZE: numev<=0"

IF (dimham==1) THEN
  evals(1)= ABS(ham(1))
  IF (numevecs > 0) THEN
    evecs(1,1)= 1.
  END IF
  RETURN
END IF

!................................................initializing solver data
ido_ap= 0
n_ap= dimham
nev_ap= numevals  ! number of eigenvalues to be computed
tol_ap= -1.
ncv_ap= MIN(4*(nev_ap+1),n_ap)  ! should be > 2* nev_ap
!rrr lworkl_ap= 3*ncv_ap**2 + 5*ncv_ap
lworkl_ap= ncv_ap**2 + 8*ncv_ap
info_ap= 0
ipntr_ap(:)= 0
iparam_ap(:)= 0
iparam_ap(1)= 1     ! uses the exact shift strategy
iparam_ap(3)= 600   ! on INPUT: max num Arnoldi update iterations allowed
iparam_ap(4)= 1
iparam_ap(7)= 1     ! mode
rvec_ap= .TRUE.

!................................................allocating solver arrays
ALLOCATE( resid_ap(n_ap) )
ALLOCATE( v_ap(n_ap,ncv_ap) )
ALLOCATE( workd_ap(3*n_ap) )
ALLOCATE( workl_ap(lworkl_ap) )
ALLOCATE( select_ap(ncv_ap) )
ALLOCATE( d_ap(nev_ap) )
ALLOCATE( d2_ap(nev_ap) )
ALLOCATE( ax_ap(n_ap) )

ALLOCATE( indexd(numevals) )

!print*, "ncv, nev, n", ncv_ap, nev_ap, n_ap

!.........................................ARPACK reverse communication loop
DO

  CALL dsaupd ( ido_ap, 'I', n_ap, 'SA', nev_ap, tol_ap, resid_ap,    &
       &        ncv_ap, v_ap, n_ap, iparam_ap, ipntr_ap,              &
       &        workd_ap, workl_ap, lworkl_ap, info_ap )

  IF (ido_ap == 99) EXIT
  IF (ido_ap /= 1 .AND. ido_ap /= -1) STOP "ido_ap /=+-1"

  CALL avmult( n_ap, hami, hamj, ham,                        &
       &       workd_ap(ipntr_ap(1)), workd_ap(ipntr_ap(2)) )

  !PRINT*, "ARPACK loop - converged Ritz values: ", iparam_ap(5)

END DO

IF (info_ap < 0 .OR. info_ap > 1) THEN
  PRINT*, "DIAGONALIZE: in ARPACK loop znaupd, info_ap=", info_ap
  PRINT*, n_ap, nev_ap, ncv_ap
  STOP
ELSE IF (info_ap == 1) THEN
  PRINT*, "ARPACK WARNING in DIAGONALIZE: in loop znaupd"
  PRINT*, "info_ap, iparam_ap(5)  =", info_ap, iparam_ap(5)
END IF
!PRINT*, "Ok, number of ARPACK Arnoldi iterations taken: ", iparam_ap(3)

!................................................ARPACK results extraction
CALL dseupd ( rvec_ap, 'A', select_ap, d_ap, v_ap, n_ap, sigma_ap,       &
     &        'I', n_ap, 'SA', nev_ap, tol_ap, resid_ap, ncv_ap, v_ap,   &
     &        n_ap, iparam_ap, ipntr_ap, workd_ap, workl_ap, lworkl_ap,  &
     &        info_ap )

!PRINT*, "iparam_ap(5), nev_ap, ncv_ap:"
!PRINT*,  iparam_ap(5), nev_ap, ncv_ap
!PRINT*, "ipntr_ap=", ipntr_ap
!PRINT*

IF (info_ap /= 0) THEN
  PRINT*, "DIAGONALIZE: in zneupd, info_ap=", info_ap
  STOP
END IF

!...................................computes and displays the residual norm
nconv= iparam_ap(5)
IF (nconv /= numevals) THEN
  PRINT*, "ARPACK WARNING ! num of computed evals with desired tol:", nconv
  STOP
END IF
DO nn= 1, nconv
  CALL avmult( n_ap, hami, hamj, ham,   &
       &      v_ap(:,nn), ax_ap )
!rrr  CALL zaxpy(n_ap, -d_ap(nn), v_ap(:,nn), 1, ax_ap, 1)
  CALL daxpy(n_ap, -d_ap(nn), v_ap(:,nn), 1, ax_ap, 1)
!rrr  rd_ap(nn,1) = REAL(d_ap(nn))
!rrr  rd_ap(nn,2) = AIMAG(d_ap(nn))
!rrr  rd_ap(nn,3) = dznrm2(n_ap, ax_ap, 1)
!rrr  rd_ap(nn,3) = rd_ap(nn,3) / dlapy2(rd_ap(nn,1),rd_ap(nn,2))
  d2_ap(nn) = dnrm2(n_ap, ax_ap, 1)
  d2_ap(nn) = d2_ap(nn) / ABS(d_ap(nn))
END DO
! uncomment to display infos on residual norm
!rrr CALL dmout(6, nconv, 3, rd_ap, nev_ap, -6,                   &
!     &     'Ritz values (Real,Imag) and relative residuals')
!CALL dmout(6, nconv, 1, d2_ap, nev_ap, -6,                   &
!     &            'Ritz values relative residuals')


!...............................fills evals and evecs in increasing order
CALL indexx(numevals, d_ap(1:numevals), indexd(:))

DO nn= 1, numevals
  evals(nn)= d_ap(indexd(nn))
END DO

DO nn= 1, numevecs
  DO kk= 1, dimham
    evecs(kk,nn)= v_ap(kk,indexd(nn))
  END DO
  norm= SUM( ABS(evecs(:,nn))**2 )
  evecs(:,nn)= evecs(:,nn) / norm
END DO

!............................................................finalizations
DEALLOCATE( indexd )

DEALLOCATE( resid_ap )
DEALLOCATE( v_ap )
DEALLOCATE( workd_ap )
DEALLOCATE( workl_ap )
DEALLOCATE( select_ap )
DEALLOCATE( d_ap )
DEALLOCATE( d2_ap )
DEALLOCATE( ax_ap)

CONTAINS

!!$!************************************************************************
!!$  SUBROUTINE avmult(dimat, maxnonzero, imat, jmat, vmat, vap, wap)
!!$    IMPLICIT NONE
!!$    ! computes  w = M v
!!$    ! for Hermitian matrices M in Compressed Row Storage format
!!$    INTEGER,    INTENT(IN) :: dimat, maxnonzero
!!$    INTEGER,    INTENT(IN) :: imat(dimat+1), jmat(maxnonzero)
!!$    COMPLEX*16, INTENT(IN) :: vmat(maxnonzero)
!!$    COMPLEX*16, INTENT(IN) :: vap(dimat)
!!$    COMPLEX*16, INTENT(OUT) :: wap(dimat)
!!$    
!!$    INTEGER :: ni, nj
!!$
!!$    wap(:)= (0.,0.)
!!$    !
!!$    DO ni= 1, dimat
!!$      wap(ni) = wap(ni) + vmat(imat(ni)) * vap(jmat(imat(ni)))
!!$      DO nj = imat(ni)+1, imat(ni+1) - 1
!!$        wap(ni) = wap(ni) + vmat(nj) * vap(jmat(nj))
!!$        wap(jmat(nj)) = wap(jmat(nj)) + CONJG(vmat(nj)) * vap(ni)
!!$      END DO
!!$    END DO
!!$
!!$  END SUBROUTINE avmult
!!$!************************************************************************

  SUBROUTINE avmult(dimat, imat, jmat, vmat, vap, wap)
    IMPLICIT NONE
    ! computes  w = M v
    ! for Hermitian matrices M in Compressed Row Storage format
    INTEGER,    INTENT(IN) :: dimat
    INTEGER,    INTENT(IN) :: imat(1:) !(dimat+1)
    INTEGER,    INTENT(IN) :: jmat(1:)
    REAL*8, INTENT(IN) :: vmat(1:)
    REAL*8, INTENT(IN) :: vap(1:dimat)
    REAL*8, INTENT(OUT) :: wap(1:dimat)
    
    INTEGER :: ni, nj

    wap(1:dimat)= (0.,0.)
    !
    DO ni= 1, dimat
      wap(ni) = wap(ni) + vmat(imat(ni)) * vap(jmat(imat(ni)))
      DO nj = imat(ni)+1, imat(ni+1) - 1
        wap(ni) = wap(ni) + REAL(vmat(nj),8) * vap(jmat(nj))
        wap(jmat(nj)) = wap(jmat(nj)) + vmat(nj) * vap(ni)
      END DO
    END DO

  END SUBROUTINE avmult
!************************************************************************


END SUBROUTINE DIAGONALIZE


END MODULE mod_diagonalize
