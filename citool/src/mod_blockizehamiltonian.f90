!$$$$$$$$$$$$$$$$$$$$$$$$$  CItool  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_blockizehamiltonian
  IMPLICIT  NONE
  SAVE

CONTAINS

SUBROUTINE BLOCKIZEHAMILTONIAN( dimhspace, ket, cutoff,            &
     &  numci_ee, ci_ee, ciindex_ee, numci_hh, ci_hh, ciindex_hh,  &
     &  numci_eh, ci_eh, ciindex_eh, numblock, blockstart, blocknonzero )
  USE mod_indata
  IMPLICIT NONE
! reorders the basis ket(:,:) in order to have an Hamiltonian
! with explicit blocks
! also computes the nonzero elements

INTEGER, INTENT(IN) :: dimhspace     ! Hilbert space basis for both elecs and holes
INTEGER*8, INTENT(INOUT) :: ket(:,:)  !(dimhspace,2)
REAL*8, INTENT(IN) :: cutoff
INTEGER, INTENT(IN) :: numci_ee                ! ee Coulomb integrals
TYPE( ci_type_real8 ), INTENT(IN) :: ci_ee(:)  !(numci_ee)
INTEGER, INTENT(IN) :: ciindex_ee(:,:,:,:) !(numspstates_e,numspstates_e,numspstates_e,numspstates_e)
INTEGER, INTENT(IN) :: numci_hh                ! eh Coulomb integrals
TYPE( ci_type_real8 ), INTENT(IN) :: ci_hh(:)  !(numci_hh)
INTEGER, INTENT(IN) :: ciindex_hh(:,:,:,:) !(numspstates_h,numspstates_h,numspstates_h,numspstates_h)
INTEGER, INTENT(IN) :: numci_eh                ! hh Coulomb integrals
TYPE( ci_type_real8 ), INTENT(IN) :: ci_eh(:)  !(numci_eh)
INTEGER, INTENT(IN) :: ciindex_eh(:,:,:,:) !(numspstates_h,numspstates_e,numspstates_e,numspstates_h)
INTEGER, INTENT(OUT) :: numblock
INTEGER, INTENT(OUT) :: blockstart(:)  !(dimhspace+1)
INTEGER, INTENT(OUT) :: blocknonzero(:)  !(dimhspace)

INTEGER :: currentblock
INTEGER :: colfixed, numprev0
LOGICAL :: elexists

INTEGER :: row, col
INTEGER*8 :: brae, brah, kete, keth
INTEGER*8 :: braebit, brahbit, ketebit, kethbit, ketbraexor, ketbrahxor
INTEGER :: emoved, hmoved
REAL*8 :: uee, uhh, ueh
INTEGER :: ni, nj, nk, nl   ! counters on single part states


IF (cutoff < 0.) STOP "BLOCKIZEHAMILTONIAN: negative cutoff"
currentblock= 1
blocknonzero(:)= 0
blockstart(1)= 1
colfixed= 1

DO row= 1, dimhspace
  IF (MOD(row,100) == 0)  print*, "row, dimhspace=", row, dimhspace
  numprev0= 0
  ! colfixed= MAX( colfixed, row )
  IF ( colfixed < row ) THEN
    IF (colfixed+1/=row) STOP "BLOCKIZEHAMILTONIAN: colfixed+1 /= row"
    colfixed= colfixed + 1
    currentblock= currentblock + 1
    blockstart(currentblock)= row
  END IF

  ! bra for elects
  !q  brae= ket_e( (row-1)/dimhspace_h + 1 )
  brae= ket( row, 1 )
  ! bra for holes
  !q  brah= ket_h( MOD(row-1,dimhspace_h) + 1 )
  brah= ket( row, 2 )

  ! diagonal element always there
  blocknonzero(currentblock)= blocknonzero(currentblock) + 1

  ! check following nonzero elems in same row in the upper triangle
  DO col= row+1, dimhspace
    elexists= .FALSE.

    ! ket for elects
    ! kete= ket_e( (col-1)/dimhspace_h + 1 )
    kete= ket( col, 1 )
    ! ket for holes
    ! keth= ket_h( MOD(col-1,dimhspace_h) + 1 )
    keth= ket( col, 2 )
    ! PRINT*, row, nbe, nbh
    ! PRINT*, col, nke, nkh
    ! PRINT*
    ketbraexor= IEOR(kete,brae)
    emoved= POPCNT( ketbraexor )
    ketbrahxor= IEOR(keth,brah)
    hmoved= POPCNT( ketbrahxor )

    IF (emoved==0 .AND. hmoved==2) THEN  !   *** 1 hole is moved ***

      kethbit= IAND(ketbrahxor,keth)     ! from ni position
      DO ni= 1, numspstates_h
        IF ( BTEST(kethbit, ni-1) ) EXIT
      END DO
      brahbit= IAND(ketbrahxor,brah)     ! to nj position
      DO nj= 1, numspstates_h
        IF ( BTEST(brahbit, nj-1) ) EXIT
      END DO
      
      uhh= 0.
      DO nk= 1, ni-1
        IF ( BTEST(keth, nk-1) ) THEN
          CALL ADDCI(uhh, ci_hh, ciindex_hh, nk, nj, ni, nk)
          CALL SUBCI(uhh, ci_hh, ciindex_hh, nj, nk, ni, nk)
        END IF
      END DO
      DO nk= ni+1, numspstates_h
        IF ( BTEST(keth, nk-1) ) THEN
          CALL ADDCI(uhh, ci_hh, ciindex_hh, nj, nk, nk, ni)
          CALL SUBCI(uhh, ci_hh, ciindex_hh, nk, nj, nk, ni)
        END IF
      END DO

      ueh= 0.
      DO nk= 1, numspstates_e
        IF ( BTEST(kete, nk-1) ) THEN
          CALL ADDCI(ueh, ci_eh, ciindex_eh, nj, nk, nk, ni)
        END IF
      END DO
      
      CALL  CHECKEMENT(elexists, cutoff, uhh + ueh)

    ELSE IF (emoved==2 .AND. hmoved==0) THEN  ! *** 1 electron is moved ***

      ketebit= IAND(ketbraexor,kete)     ! from ni position
      DO ni= 1, numspstates_e
        IF ( BTEST(ketebit, ni-1) ) EXIT
      END DO
      braebit= IAND(ketbraexor,brae)     ! to nj position
      DO nj= 1, numspstates_e
        IF ( BTEST(braebit, nj-1) ) EXIT
      END DO

      uee= 0.      
      DO nk= 1, ni-1
        IF ( BTEST(kete, nk-1) ) THEN
          CALL ADDCI(uee, ci_ee, ciindex_ee, nk, nj, ni, nk)
          CALL SUBCI(uee, ci_ee, ciindex_ee, nj, nk, ni, nk)
        END IF
      END DO
      DO nk= ni+1, numspstates_e
        IF ( BTEST(kete, nk-1) ) THEN
          CALL ADDCI(uee, ci_ee, ciindex_ee, nj, nk, nk, ni)
          CALL SUBCI(uee, ci_ee, ciindex_ee, nk, nj, nk, ni)
        END IF
      END DO

      ueh= 0.
      DO nk= 1, numspstates_h
        IF ( BTEST(keth, nk-1) ) THEN
          CALL ADDCI(ueh, ci_eh, ciindex_eh, nk, nj, ni, nk)
        END IF
      END DO

      CALL  CHECKEMENT(elexists, cutoff, uee + ueh)

    ELSE IF (emoved==2 .AND. hmoved==2) THEN  !*** 1 hole and 1 elec moved ***
      
      kethbit= IAND(ketbrahxor,keth)     ! hole moved from nl position
      DO nl= 1, numspstates_h
        IF ( BTEST(kethbit, nl-1) ) EXIT
      END DO
      brahbit= IAND(ketbrahxor,brah)     !              to ni position
      DO ni= 1, numspstates_h
        IF ( BTEST(brahbit, ni-1) ) EXIT
      END DO

      ketebit= IAND(ketbraexor,kete)     ! elec moved from nk position
      DO nk= 1, numspstates_e
        IF ( BTEST(ketebit, nk-1) ) EXIT
      END DO
      braebit= IAND(ketbraexor,brae)     ! to nj position
      DO nj= 1, numspstates_e
        IF ( BTEST(braebit, nj-1) ) EXIT
      END DO

      ueh= 0.
      CALL ADDCI(ueh, ci_eh, ciindex_eh, ni, nj, nk, nl)

      CALL  CHECKELEMENT(elexists, cutoff, ueh)

    ELSE IF (emoved==0 .AND. hmoved==4) THEN   ! *** 2 holes are moved ***

      kethbit= IAND(ketbrahxor,keth)     ! 2 holes moved from nl & nk pos
      DO nl= 1, numspstates_h
        IF ( BTEST(kethbit, nl-1) ) EXIT
      END DO
      DO nk= nl+1, numspstates_h
        IF ( BTEST(kethbit, nk-1) ) EXIT
      END DO

      brahbit= IAND(ketbrahxor,brah)     !              to nj & ni pos
      DO nj= 1, numspstates_h
        IF ( BTEST(brahbit, nj-1) ) EXIT
      END DO
      DO ni= nj+1, numspstates_h
        IF ( BTEST(brahbit, ni-1) ) EXIT
      END DO

      uhh= 0.
      CALL ADDCI(uhh, ci_hh, ciindex_hh, ni, nj, nk, nl)
      CALL SUBCI(uhh, ci_hh, ciindex_hh, nj, ni, nk, nl)

      CALL  CHECKELEMENT(elexists, cutoff, uhh)

    ELSE IF (emoved==4 .AND. hmoved==0) THEN   ! *** 2 elecs are moved ***

      ketebit= IAND(ketbraexor,kete)     ! 2 elecs moved from nl & nk pos
      DO nl= 1, numspstates_e
        IF ( BTEST(ketebit, nl-1) ) EXIT
      END DO
      DO nk= nl+1, numspstates_e
        IF ( BTEST(ketebit, nk-1) ) EXIT
      END DO

      braebit= IAND(ketbraexor,brae)     !              to nj & ni pos
      DO nj= 1, numspstates_e
        IF ( BTEST(braebit, nj-1) ) EXIT
      END DO
      DO ni= nj+1, numspstates_e
        IF ( BTEST(braebit, ni-1) ) EXIT
      END DO

      uee= 0.
      CALL CHECKCI(uee, ci_ee, ciindex_ee, ni, nj, nk, nl)
      CALL CHECKCI(uee, ci_ee, ciindex_ee, nj, ni, nk, nl)
      
      CALL  CHECKELEMENT(elexists, cutoff, uee)
      
    END IF

    IF ( elexists ) THEN
      blocknonzero(currentblock)= blocknonzero(currentblock) + 1
    END IF

    IF ( col > colfixed ) THEN

      IF ( elexists ) THEN
        IF ( numprev0 > 0 ) THEN   ! swaps the ket with the first that has 0 coupling
          colfixed= colfixed + 1   ! after clofixed: it sould be in colfixed+1
          kete= ket(col,1)
          keth= ket(col,2)
          ket(col,1)= ket(colfixed,1)
          ket(col,2)= ket(colfixed,2)
          ket(colfixed,1)= kete
          ket(colfixed,2)= keth
          numprev0= numprev0 - 1
        ELSE
          colfixed= colfixed + 1
        END IF
      ELSE
        numprev0= numprev0 + 1
      END IF

    END IF


  END DO ! col

END DO ! row

blockstart(currentblock+1)= dimhspace + 1
numblock= currentblock


!Q hami(dimham+1)= nonzero + 1

!!$PRINT*, "npm00=", npm00 + dimham
!!$PRINT*, "npm01=", npm01, dimham*((numspstates_h-num_h)*num_h)/2
!!$PRINT*, "npm10=", npm10, dimham*((numspstates_e-num_e)*num_e)/2
!!$PRINT*, "npm11=", npm11, dimham*((numspstates_h-num_h)*num_h*(numspstates_e-num_e)*num_e)/2
!!$PRINT*, "npm02=", npm02, dimham*((numspstates_h-num_h)*num_h*(numspstates_h-num_h-1)*(num_h-1)/4)/2
!!$PRINT*, "npm20=", npm20, dimham*((numspstates_e-num_e)*num_e*(numspstates_e-num_e-1)*(num_e-1)/4)/2
!!$PRINT*, "npmx=", npmx
!!$PRINT*, "non zero elems:", npm00+dimham+npm01+npm10+npm11+npm02+npm20

CONTAINS

  SUBROUTINE ADDCI(uu, ci, ciindex, n1, n2, n3, n4)
    IMPLICIT NONE
    REAL*8, INTENT(INOUT) :: uu
    TYPE( ci_type_real8 ), INTENT(IN) :: ci(:)
    INTEGER, INTENT(IN) :: ciindex(:,:,:,:)
    INTEGER, INTENT(IN) :: n1, n2, n3, n4

    IF (ciindex(n1,n2,n3,n4)>0) THEN
      uu= uu + ci(ciindex(n1,n2,n3,n4))%v
    !ELSE IF (ciindex(n1,n2,n3,n4)<0) THEN
    !  uu= uu + ci(-ciindex(n1,n2,n3,n4))%v
    END IF
  END SUBROUTINE ADDCI

  SUBROUTINE SUBCI(uu, ci, ciindex, n1, n2, n3, n4)
    IMPLICIT NONE
    REAL*8, INTENT(INOUT) :: uu
    TYPE( ci_type_real8 ), INTENT(IN) :: ci(:)
    INTEGER, INTENT(IN) :: ciindex(:,:,:,:)
    INTEGER, INTENT(IN) :: n1, n2, n3, n4

    IF (ciindex(n1,n2,n3,n4)>0) THEN
      uu= uu - ci(ciindex(n1,n2,n3,n4))%v
    !ELSE IF (ciindex(n1,n2,n3,n4)<0) THEN
    !  uu= uu - ci(-ciindex(n1,n2,n3,n4))%v
    END IF
  END SUBROUTINE SUBCI
  
  SUBROUTINE CHECKELEMENT(elexists, cutoff, val)
    IMPLICIT NONE
    LOGICAL, INTENT(INOUT) :: elexists
    REAL*8, INTENT(IN) :: cutoff
    REAL*8, INTENT(IN) :: val
    IF ( ABS(val) > cutoff ) THEN
      elexists= .TRUE.
    END IF
  END SUBROUTINE CHECKELEMENT


END SUBROUTINE BLOCKIZEHAMILTONIAN


END MODULE mod_blockizehamiltonian
