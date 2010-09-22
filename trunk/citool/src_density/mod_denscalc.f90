!$$$$$$$$$$$$$$$$$$$$$$$$$$$  density4CI  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_denscalc
  IMPLICIT NONE
  SAVE

CONTAINS

SUBROUTINE DENSCALC( blockdim, mpstate, numspstates_C, ket_C, numx_C, psi_C, &
     &  mask, numspstates_O, ket_O, dens )

IMPLICIT NONE
!
! "C" stands for e or h i.e. the carrier for which the dens is obtained
! "O" stands for h or e i.e. the other carrier
!
! "S" stands for 0 or 1 i.e. the spin for which the dens is obtained
! "N" stands for 1 or 0 i.e. the other spin
!
INTEGER, INTENT(IN) :: blockdim
REAL*8,  INTENT(IN) :: mpstate(:)   !(blockdim)
INTEGER, INTENT(IN) :: numspstates_C
INTEGER*8, INTENT(IN) :: ket_C(:)   !(blockdim)
INTEGER, INTENT(IN) :: numx_C
REAL*8,  INTENT(IN) :: psi_C(:,:)   !(numx)
LOGICAL, INTENT(IN) :: mask(:)      !(numspstates_C)
INTEGER, INTENT(IN) :: numspstates_O
INTEGER*8, INTENT(IN) :: ket_O(:)   !(blockdim)
REAL*8, INTENT(OUT) :: dens(:)      !(numx)

INTEGER, ALLOCATABLE :: matcc(:,:)
INTEGER*8 :: klef, krig
REAL*8 :: sign, clcl
INTEGER :: nl, nlp, ns, nsp, nb
INTEGER :: nx

ALLOCATE( matcc(numspstates_C, numspstates_C) )

dens= 0d0

DO nl= 1, blockdim
  DO nlp= 1, blockdim
    IF ( ket_O(nl) /= ket_O(nlp) )  CYCLE   
    ! here I suppose that the unused sp states bits are set to zero

    !PRINT*, "nl, nl' =", nl, nlp

    clcl= mpstate(nl)*mpstate(nlp)

    matcc= 0
    DO ns= 1, numspstates_C
      IF ( .NOT. mask(ns) ) CYCLE
      DO nsp= 1, numspstates_C
        IF ( .NOT. mask(nsp) ) CYCLE
        klef= ket_C(nlp)
        krig= ket_C(nl)
        IF ( BTEST(krig,ns-1) .AND. BTEST(klef,nsp-1) ) THEN
          krig= IBCLR(krig,ns-1)  ! destroy ns in right ket
          sign= 1
          DO nb= 1, ns-1
            IF ( BTEST(krig,nb-1) )  sign= -sign
          END DO
          klef= IBCLR(klef,nsp-1)  ! destroy nsp in left ket
          DO nb= 0, nsp-1
            IF ( BTEST(klef,nb-1) )  sign= -sign
          END DO
          IF ( klef == krig ) THEN
            matcc(nsp,ns)= sign
          END IF
        END IF
      END DO
    END DO
    
    DO ns= 1, numspstates_C
      DO nsp= 1, numspstates_C
        IF (matcc(nsp,ns) /= 0) THEN

          DO nx= 1, numx_C
            dens(nx)= dens(nx) + clcl * matcc(nsp,ns) * psi_C(nx,nsp)*psi_C(nx,ns)
          END DO
                       
        END IF
      END DO
    END DO
          
  END DO
END DO

DEALLOCATE( matcc )

END SUBROUTINE DENSCALC


END MODULE mod_denscalc
