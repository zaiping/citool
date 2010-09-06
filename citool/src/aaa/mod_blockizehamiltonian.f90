!$$$$$$$$$$$$$$$$$$$$$$$$$  CItool  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_blockizehamiltonian
  IMPLICIT  NONE
  SAVE

CONTAINS

SUBROUTINE BLOCKIZEHAMILTONIAN( dimhspace, ket,                   &
     &  ciindex_ee, ciindex_hh, ciindex_eh, numblock, blockstart, blocknonzero )
  USE mod_indata
  IMPLICIT NONE
! reorders the basis ket(:,:) in order to have an Hamiltonian
! with explicit blocks
! also computes the nonzero elements

INTEGER, INTENT(IN) :: dimhspace     ! Hilbert space basis for both elecs and holes
INTEGER*8, INTENT(INOUT) :: ket(:,:)  !(dimhspace,2)
INTEGER, INTENT(IN) :: ciindex_ee(:,:,:,:) !(numspstates_e,numspstates_e,numspstates_e,numspstates_e)
INTEGER, INTENT(IN) :: ciindex_hh(:,:,:,:) !(numspstates_h,numspstates_h,numspstates_h,numspstates_h)
INTEGER, INTENT(IN) :: ciindex_eh(:,:,:,:) !(numspstates_h,numspstates_e,numspstates_e,numspstates_h)
INTEGER, INTENT(OUT) :: numblock
INTEGER, INTENT(OUT) :: blockstart(:)  !(dimhspace+1)
INTEGER, INTENT(OUT) :: blocknonzero(:)  !(dimhspace)

INTEGER :: currentblock
INTEGER :: colfixed, numprev0
LOGICAL :: ciexists

INTEGER :: row, col
INTEGER*8 :: brae, brah, kete, keth
INTEGER*8 :: braebit, brahbit, ketebit, kethbit, ketbraexor, ketbrahxor
INTEGER :: emoved, hmoved
REAL*8 :: sig
INTEGER :: ni, nj, nk, nl   ! counters on single part states
REAL*8 :: tinye


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
    ciexists= .FALSE.

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
      
      DO nk= 1, ni-1
        IF ( BTEST(keth, nk-1) ) THEN
          CALL CHECKCI(ciexists, ciindex_hh, nk, nj, ni, nk)
          CALL CHECKCI(ciexists, ciindex_hh, nj, nk, ni, nk)
        END IF
      END DO
      DO nk= ni+1, numspstates_h
        IF ( BTEST(keth, nk-1) ) THEN
          CALL CHECKCI(ciexists, ciindex_hh, nj, nk, nk, ni)
          CALL CHECKCI(ciexists, ciindex_hh, nk, nj, nk, ni)
        END IF
      END DO

      DO nk= 1, numspstates_e
        IF ( BTEST(kete, nk-1) ) THEN
          CALL CHECKCI(ciexists, ciindex_eh, nj, nk, nk, ni)
        END IF
      END DO

    ELSE IF (emoved==2 .AND. hmoved==0) THEN  ! *** 1 electron is moved ***

      ketebit= IAND(ketbraexor,kete)     ! from ni position
      DO ni= 1, numspstates_e
        IF ( BTEST(ketebit, ni-1) ) EXIT
      END DO
      braebit= IAND(ketbraexor,brae)     ! to nj position
      DO nj= 1, numspstates_e
        IF ( BTEST(braebit, nj-1) ) EXIT
      END DO
      
      DO nk= 1, ni-1
        IF ( BTEST(kete, nk-1) ) THEN
          CALL CHECKCI(ciexists, ciindex_ee, nk, nj, ni, nk)
          CALL CHECKCI(ciexists, ciindex_ee, nj, nk, ni, nk)
        END IF
      END DO
      DO nk= ni+1, numspstates_e
        IF ( BTEST(kete, nk-1) ) THEN
          CALL CHECKCI(ciexists, ciindex_ee, nj, nk, nk, ni)
          CALL CHECKCI(ciexists, ciindex_ee, nk, nj, nk, ni)
        END IF
      END DO

      DO nk= 1, numspstates_h
        IF ( BTEST(keth, nk-1) ) THEN
          CALL CHECKCI(ciexists, ciindex_eh, nk, nj, ni, nk)
        END IF
      END DO

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

      CALL CHECKCI(ciexists, ciindex_eh, ni, nj, nk, nl)

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

      CALL CHECKCI(ciexists, ciindex_hh, ni, nj, nk, nl)
      CALL CHECKCI(ciexists, ciindex_hh, nj, ni, nk, nl)

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

      CALL CHECKCI(ciexists, ciindex_ee, ni, nj, nk, nl)
      CALL CHECKCI(ciexists, ciindex_ee, nj, ni, nk, nl)
      
    END IF

    IF ( ciexists ) THEN
      blocknonzero(currentblock)= blocknonzero(currentblock) + 1
    END IF

    IF ( col > colfixed ) THEN

      IF ( ciexists ) THEN
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

  SUBROUTINE CHECKCI(ciexists, ciindex, n1, n2, n3, n4)
    IMPLICIT NONE
    LOGICAL, INTENT(INOUT) :: ciexists
    INTEGER, INTENT(IN) :: ciindex(:,:,:,:)
    INTEGER, INTENT(IN) :: n1, n2, n3, n4

    IF (ciindex(n1,n2,n3,n4) /= 0) THEN
      ciexists= .TRUE.
    END IF

  END SUBROUTINE CHECKCI


END SUBROUTINE BLOCKIZEHAMILTONIAN

END MODULE mod_blockizehamiltonian
