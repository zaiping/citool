MODULE mod_composehspace
  IMPLICIT  NONE
  SAVE

CONTAINS

SUBROUTINE composehspace( numspstates_e, numspqn_e, spqn_e, dimhspace_e, ket_e, consline_e,  &
      &                   numspstates_h, numspqn_h, spqn_h, dimhspace_h, ket_h, consline_h,  &
      &                   ket, dimhspacecons )
  IMPLICIT NONE
! composed e and h H spaces in a single H space excluding vectors not meeting
! the constrains

INTEGER, INTENT(IN) :: numspstates_e, numspqn_e
INTEGER, INTENT(IN) :: spqn_e(1:,1:)
INTEGER, INTENT(IN) :: dimhspace_e
INTEGER*8, INTENT(IN) :: ket_e(1:)
INTEGER, INTENT(IN) :: consline_e(1:)
INTEGER, INTENT(IN) :: numspstates_h, numspqn_h
INTEGER, INTENT(IN) :: spqn_h(1:,1:)
INTEGER, INTENT(IN) :: dimhspace_h
INTEGER*8, INTENT(IN) :: ket_h(1:)
INTEGER, INTENT(IN) :: consline_h(1:)
INTEGER*8, INTENT(OUT) :: ket(1:,1:)
INTEGER, INTENT(OUT) :: dimhspacecons

LOGICAL :: incl_e(dimhspace_e), incl_h(dimhspace_h)
INTEGER :: dimhspace, sumqn, nke, nkh
INTEGER :: nn, nqn, ns

dimhspace= dimhspace_e * dimhspace_h
incl_e(:)= .TRUE.
incl_h(:)= .TRUE.

DO nn= 1, dimhspace_e
  DO nqn= 1, numspqn_e
    IF (consline_e(nqn) /= 9999) THEN    
      sumqn= 0
      DO ns= 1, numspstates_e
        IF (BTEST(ket_e(nn),ns-1))  sumqn= sumqn + spqn_e(ns,nqn)
      END DO
      IF (sumqn /= consline_e(nqn))  incl_e(nn)= .FALSE.
    END IF
  END DO
END DO

DO nn= 1, dimhspace_h
  DO nqn= 1, numspqn_h
    IF (consline_h(nqn) /= 9999) THEN    
      sumqn= 0
      DO ns= 1, numspstates_h
        IF (BTEST(ket_h(nn),ns-1))  sumqn= sumqn + spqn_h(ns,nqn)
      END DO
      IF (sumqn /= consline_h(nqn))  incl_h(nn)= .FALSE.
    END IF
  END DO
END DO

ket(:,:)= 0
dimhspacecons= 0
DO nn= 1, dimhspace
  nke= (nn-1)/dimhspace_h + 1
  nkh= MOD(nn-1,dimhspace_h) + 1
  IF (incl_e(nke) .AND. incl_h(nkh)) THEN
    dimhspacecons= dimhspacecons + 1
    ket(dimhspacecons,1)= ket_e(nke)
    ket(dimhspacecons,2)= ket_h(nkh)
  END IF
END DO

END SUBROUTINE composehspace


END MODULE mod_composehspace
