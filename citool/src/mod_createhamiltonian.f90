!$$$$$$$$$$$$$$$$$$$$$$$$$  CItool  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_createhamiltonian
  IMPLICIT  NONE
  SAVE

CONTAINS

SUBROUTINE CREATEHAMILTONIAN_X( dimhspace, ket, cutoff,            &
     &  numci_ee, ci_ee, ciindex_ee, numci_hh, ci_hh, ciindex_hh,  &
     &  numci_eh, ci_eh, ciindex_eh, spenergy_e, spenergy_h,       &
     &  nonzero, hami, hamj, ham )
  USE mod_indata
  IMPLICIT NONE
! creates the COMPLEX Hamiltonian matrix H= Te + Th + Vee + Vhh + Veh
! where the V's are COMPLEX Coulomb interaction terms
! it used the Compressed Row Storage format (same as Pardiso)
! see  http://www.cs.utk.edu/~dongarra/etemplates/node371.html

INTEGER, INTENT(IN) :: dimhspace            ! Hilbert space basis for both elecs and holes
INTEGER*8, INTENT(IN) :: ket(1:,:)  !(dimhspace,2)
REAL*8, INTENT(IN) :: cutoff
INTEGER, INTENT(IN) :: numci_ee                    ! ee Coulomb integrals
TYPE( ci_type_complex16 ), INTENT(IN) :: ci_ee(:)  !(numci_ee)
INTEGER, INTENT(IN) :: ciindex_ee(:,:,:,:) !(numspstates_e,numspstates_e,numspstates_e,numspstates_e)
INTEGER, INTENT(IN) :: numci_hh                    ! eh Coulomb integrals
TYPE( ci_type_complex16 ), INTENT(IN) :: ci_hh(:)  !(numci_hh)
INTEGER, INTENT(IN) :: ciindex_hh(:,:,:,:) !(numspstates_h,numspstates_h,numspstates_h,numspstates_h)
INTEGER, INTENT(IN) :: numci_eh                    ! hh Coulomb integrals
TYPE( ci_type_complex16 ), INTENT(IN) :: ci_eh(:)  !(numci_eh)
INTEGER, INTENT(IN) :: ciindex_eh(:,:,:,:) !(numspstates_h,numspstates_e,numspstates_e,numspstates_h)
REAL*8, INTENT(IN) :: spenergy_e(:)  !(numspstates_e)
REAL*8, INTENT(IN) :: spenergy_h(:)  !(numspstates_h)
INTEGER, INTENT(IN) :: nonzero
INTEGER, INTENT(OUT) :: hami(:)
INTEGER, INTENT(OUT) :: hamj(:)
COMPLEX*16, INTENT(OUT) :: ham(:)

INTEGER :: row, col
INTEGER :: nonzerocount
INTEGER*8 :: brae, brah, kete, keth
INTEGER*8 :: braebit, brahbit, ketebit, kethbit, ketbraexor, ketbrahxor
INTEGER :: emoved, hmoved
COMPLEX*16 :: teth, uee, uhh, ueh
REAL*8 :: sig
INTEGER :: ni, nj, nk, nl   ! counters on single part states

!!$
!!$ ** with blockizehamiltonian this is not necessary: kept for historical reason
!!$  
! the max number of non-zero elements is (BC: binomial coefficients)
!maxnonzero= dimham * ( 1                                                  &
!     & + BC(numspstates_h-num_h,1)*num_h                                  &
!     & + BC(numspstates_e-num_e,1)*num_e                                  &
!     & + BC(numspstates_h-num_h,1)*num_h*BC(numspstates_e-num_e,1)*num_e  &
!     & + BC(numspstates_h-num_h,2)*BC(num_h,2)                            &
!     & + BC(numspstates_e-num_e,2)*BC(num_e,2)
!
! the different terms correspond to the following possibilities:
!-----------------------------------------------------------
!   elecs moved    holes moved    H terms /= 0
!-----------------------------------------------------------
!       0              0          Te + Th + Vee + Vhh + Veh
!       0              1          Vhh + Veh
!       1              0          Vee + Veh
!       1              1          Veh
!       0              2          Vhh
!       2              0          Vee
!
! IF ( MOD(maxnonzero-dimham,2)>0 ) STOP "CREATEHAMILTONIAN: odd nonzeros??!"
! maxnonzero= (maxnonzero-dimham)/2 + dimham
! since H is hermitian, we store only the diag and the upper triangle,
! i.e. maxnonzero only includes the upper triangle elements
!  -- the above calculation is performed in the main parogram --

IF (cutoff < 0d0) STOP "CREATEHAMILTONIAN_X: negative cutoff"
nonzerocount= 0
DO row= 1, dimhspace
  ! before the reordering the total basis was
  ! |col>=|nke,nkh> with col=(nke-1)Nh+nkh
  ! i.e. in the row/col index the hole index runs faster (same for rows)

  ! bra for elects
  brae= ket( row, 1 )
  ! bra for holes
  brah= ket( row, 2 )

  !PPPPPPPPPPPPPPPPPP
  !PRINT*, brae, brah

  ! diagonal element
  teth= 0.
  DO ni= 1, numspstates_e
    IF (BTEST(brae,ni-1))  teth= teth + spenergy_e(ni)
  END DO
  DO nj= 1, numspstates_h
    IF (BTEST(brah,nj-1))  teth= teth + spenergy_h(nj)
  END DO

  uee= 0.
  DO ni= 1, numspstates_e-1
    IF (BTEST(brae,ni-1)) THEN
      DO nj= ni+1, numspstates_e
        IF (BTEST(brae,nj-1)) THEN
          CALL ADDCI_X(uee, ci_ee, ciindex_ee, ni, nj, nj, ni)
          CALL SUBCI_X(uee, ci_ee, ciindex_ee, nj, ni, nj, ni)
        END IF
      END DO
    END IF
  END DO

  uhh= 0.
  DO ni= 1, numspstates_h-1
    IF (BTEST(brah,ni-1)) THEN
      DO nj= ni+1, numspstates_h
        IF (BTEST(brah,nj-1)) THEN
          CALL ADDCI_X(uhh, ci_hh, ciindex_hh, ni, nj, nj, ni)
          CALL SUBCI_X(uhh, ci_hh, ciindex_hh, nj, ni, nj, ni)
        END IF
      END DO
    END IF
  END DO

  ueh= 0.
  DO ni= 1, numspstates_h
    IF (BTEST(brah,ni-1)) THEN
      DO nj= 1, numspstates_e
        IF (BTEST(brae,nj-1)) THEN
          CALL ADDCI_X(ueh, ci_eh, ciindex_eh, ni, nj, nj, ni)
        END IF
      END DO
    END IF
  END DO

  ! first elem of a new row, here row=col
  hami(row)= nonzerocount+1 
  ! diag element always present, also if it is zero
  CALL  ADDELEMENT_X(nonzerocount, row, hamj, ham, -1d0, teth + uee + uhh + ueh)

  !PPPPPPPPPPPPPPPPPPPPP
  !PRINT*, teth, uee, uhh, ueh, "qqq"
  !STOP

  IF ( ABS(AIMAG(ham(nonzerocount)))*1e4 > ABS(REAL(ham(nonzerocount))) ) THEN
    PRINT*, "mod_createhamiltonian: non-real diag element !", row
    STOP
  END IF

  ! calculate following nonzero elems in same row in the upper triangle
  DO col= row+1, dimhspace
    ! ket for elects
    kete= ket( col, 1 )
    ! ket for holes
    keth= ket( col, 2 )
    ! PRINT*, row, nbe, nbh
    ! PRINT*, col, nke, nkh
    ! PRINT*
    ketbraexor= IEOR(kete,brae)
    emoved= POPCNT( ketbraexor )
    ketbrahxor= IEOR(keth,brah)
    hmoved= POPCNT( ketbrahxor )

!!$    IF (row==2 .AND. col==9) THEN
!!$      PRINT*, kete,brae,keth,brah
!!$      PRINT*, ketbraexor, emoved
!!$      PRINT*, ketbrahxor, hmoved
!!$      STOP "sssssssssssssssss"
!!$    END IF

!!$    IF (emoved==0 .AND. hmoved==0) THEN
!!$      npm00= npm00 + 1
!!$    ELSE IF (emoved==0 .AND. hmoved==2) THEN
!!$      npm01= npm01 + 1
!!$    ELSE IF (emoved==2 .AND. hmoved==0) THEN
!!$      npm10= npm10 + 1
!!$    ELSE IF (emoved==2 .AND. hmoved==2) THEN
!!$      npm11= npm11 + 1
!!$    ELSE IF (emoved==0 .AND. hmoved==4) THEN
!!$      npm02= npm02 + 1
!!$    ELSE IF (emoved==4 .AND. hmoved==0) THEN
!!$      npm20= npm20 + 1
!!$    ELSE
!!$      npmx= npmx + 1
!!$    END IF

    IF (emoved==0 .AND. hmoved==2) THEN  !   *** 1 hole is moved ***

      kethbit= IAND(ketbrahxor,keth)     ! from ni position
      DO ni= 1, numspstates_h
        IF ( BTEST(kethbit, ni-1) ) EXIT
      END DO
      brahbit= IAND(ketbrahxor,brah)     ! to nj position
      DO nj= 1, numspstates_h
        IF ( BTEST(brahbit, nj-1) ) EXIT
      END DO
      ! sign according to num of 1's between ni and nj
      sig= (-1)**(POPCNT(IBITS( keth, MIN(ni,nj), ABS(ni-nj)-1 )))

      uhh= 0.
      DO nk= 1, ni-1
        IF ( BTEST(keth, nk-1) ) THEN
          CALL ADDCI_X(uhh, ci_hh, ciindex_hh, nk, nj, ni, nk)
          CALL SUBCI_X(uhh, ci_hh, ciindex_hh, nj, nk, ni, nk)
        END IF
      END DO
      DO nk= ni+1, numspstates_h
        IF ( BTEST(keth, nk-1) ) THEN
          CALL ADDCI_X(uhh, ci_hh, ciindex_hh, nj, nk, nk, ni)
          CALL SUBCI_X(uhh, ci_hh, ciindex_hh, nk, nj, nk, ni)
        END IF
      END DO
      uhh= sig * uhh

      ueh= 0.
      DO nk= 1, numspstates_e
        IF ( BTEST(kete, nk-1) ) THEN
          CALL ADDCI_X(ueh, ci_eh, ciindex_eh, nj, nk, nk, ni)
        END IF
      END DO
      ueh= sig * ueh

      CALL  ADDELEMENT_X(nonzerocount, col, hamj, ham, cutoff, uhh + ueh)

    ELSE IF (emoved==2 .AND. hmoved==0) THEN  ! *** 1 electron is moved ***

      ketebit= IAND(ketbraexor,kete)     ! from ni position
      DO ni= 1, numspstates_e
        IF ( BTEST(ketebit, ni-1) ) EXIT
      END DO
      braebit= IAND(ketbraexor,brae)     ! to nj position
      DO nj= 1, numspstates_e
        IF ( BTEST(braebit, nj-1) ) EXIT
      END DO
      ! sign according to num of 1's between ni and nj
      sig= (-1)**(POPCNT(IBITS( kete, MIN(ni,nj), ABS(ni-nj)-1 )))

      uee= 0.
      DO nk= 1, ni-1
        IF ( BTEST(kete, nk-1) ) THEN
          CALL ADDCI_X(uee, ci_ee, ciindex_ee, nk, nj, ni, nk)
          CALL SUBCI_X(uee, ci_ee, ciindex_ee, nj, nk, ni, nk)
        END IF
      END DO
      DO nk= ni+1, numspstates_e
        IF ( BTEST(kete, nk-1) ) THEN
          CALL ADDCI_X(uee, ci_ee, ciindex_ee, nj, nk, nk, ni)
          CALL SUBCI_X(uee, ci_ee, ciindex_ee, nk, nj, nk, ni)
        END IF
      END DO
      uee= sig * uee

      ueh= 0.
      DO nk= 1, numspstates_h
        IF ( BTEST(keth, nk-1) ) THEN
          CALL ADDCI_X(ueh, ci_eh, ciindex_eh, nk, nj, ni, nk)
        END IF
      END DO
      ueh= sig * ueh

      CALL  ADDELEMENT_X(nonzerocount, col, hamj, ham, cutoff, uee + ueh)

    ELSE IF (emoved==2 .AND. hmoved==2) THEN  !*** 1 hole and 1 elec moved ***
      
      kethbit= IAND(ketbrahxor,keth)     ! hole moved from nl position
      DO nl= 1, numspstates_h
        IF ( BTEST(kethbit, nl-1) ) EXIT
      END DO
      brahbit= IAND(ketbrahxor,brah)     !              to ni position
      DO ni= 1, numspstates_h
        IF ( BTEST(brahbit, ni-1) ) EXIT
      END DO
      ! sign according to num of 1's between nl and ni
      sig= (-1)**(POPCNT(IBITS( keth, MIN(nl,ni), ABS(nl-ni)-1 )))

      ketebit= IAND(ketbraexor,kete)     ! elec moved from nk position
      DO nk= 1, numspstates_e
        IF ( BTEST(ketebit, nk-1) ) EXIT
      END DO
      braebit= IAND(ketbraexor,brae)     ! to nj position
      DO nj= 1, numspstates_e
        IF ( BTEST(braebit, nj-1) ) EXIT
      END DO
      ! sign according to num of 1's between ni and nj
      sig= sig * (-1)**(POPCNT(IBITS( kete, MIN(nk,nj), ABS(nk-nj)-1 )))

      ueh= 0.
      CALL ADDCI_X(ueh, ci_eh, ciindex_eh, ni, nj, nk, nl)
      ueh= sig * ueh

      CALL  ADDELEMENT_X(nonzerocount, col, hamj, ham, cutoff, ueh)

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

      ! sign according to num of 1's between ni & nj + nk & nl 
      ! sign according to num of 1's between ni & nj + nk & nl 
      sig= (-1)**(POPCNT(IBITS( keth, nl, nk-nl-1 )))
      sig= -sig * (-1)**(POPCNT(IBITS( brah, nj, ni-nj-1 )))

      uhh= 0.
      CALL ADDCI_X(uhh, ci_hh, ciindex_hh, ni, nj, nk, nl)
      CALL SUBCI_X(uhh, ci_hh, ciindex_hh, nj, ni, nk, nl)
      uhh= sig * uhh

      CALL  ADDELEMENT_X(nonzerocount, col, hamj, ham, cutoff, uhh)

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

!!$      IF (row==2 .AND. col==9) THEN
!!$        PRINT*, nl, nk
!!$        PRINT*, nj, ni
!!$        STOP "ttttt 1"
!!$      END IF

      ! sign according to num of 1's between ni & nj + nk & nl 
      sig= (-1)**(POPCNT(IBITS( kete, nl, nk-nl-1 )))
      sig= -sig * (-1)**(POPCNT(IBITS( brae, nj, ni-nj-1 )))

      uee= 0.
      CALL ADDCI_X(uee, ci_ee, ciindex_ee, ni, nj, nk, nl)
      CALL SUBCI_X(uee, ci_ee, ciindex_ee, nj, ni, nk, nl)
      uee= sig * uee

      CALL  ADDELEMENT_X(nonzerocount, col, hamj, ham, cutoff, uee)
      
    END IF

  END DO ! col

END DO ! row

hami(dimhspace+1)= nonzerocount + 1

!print*, "aaaaaaaaaaaaaaa", nonzerocount, nonzero
IF ( nonzerocount /= nonzero ) STOP "CREATEHAMILTONIAN_X : nonzerocount /= nonzero"

!!$PRINT*, "npm00=", npm00 + dimham
!!$PRINT*, "npm01=", npm01, dimham*((numspstates_h-num_h)*num_h)/2
!!$PRINT*, "npm10=", npm10, dimham*((numspstates_e-num_e)*num_e)/2
!!$PRINT*, "npm11=", npm11, dimham*((numspstates_h-num_h)*num_h*(numspstates_e-num_e)*num_e)/2
!!$PRINT*, "npm02=", npm02, dimham*((numspstates_h-num_h)*num_h*(numspstates_h-num_h-1)*(num_h-1)/4)/2
!!$PRINT*, "npm20=", npm20, dimham*((numspstates_e-num_e)*num_e*(numspstates_e-num_e-1)*(num_e-1)/4)/2
!!$PRINT*, "npmx=", npmx
!!$PRINT*, "non zero elems:", npm00+dimham+npm01+npm10+npm11+npm02+npm20

CONTAINS

  SUBROUTINE ADDCI_X(uu, ci, ciindex, n1, n2, n3, n4)
    IMPLICIT NONE
    COMPLEX*16, INTENT(INOUT) :: uu
    TYPE( ci_type_complex16 ), INTENT(IN) :: ci(:)
    INTEGER, INTENT(IN) :: ciindex(:,:,:,:)
    INTEGER, INTENT(IN) :: n1, n2, n3, n4

    IF (ciindex(n1,n2,n3,n4)>0) THEN
      uu= uu + ci(ciindex(n1,n2,n3,n4))%v
    ELSE IF (ciindex(n1,n2,n3,n4)<0) THEN
      uu= uu + CONJG(ci(-ciindex(n1,n2,n3,n4))%v)
    END IF
  END SUBROUTINE ADDCI_X

  SUBROUTINE SUBCI_X(uu, ci, ciindex, n1, n2, n3, n4)
    IMPLICIT NONE
    COMPLEX*16, INTENT(INOUT) :: uu
    TYPE( ci_type_complex16 ), INTENT(IN) :: ci(:)
    INTEGER, INTENT(IN) :: ciindex(:,:,:,:)
    INTEGER, INTENT(IN) :: n1, n2, n3, n4

    IF (ciindex(n1,n2,n3,n4)>0) THEN
      uu= uu - ci(ciindex(n1,n2,n3,n4))%v
    ELSE IF (ciindex(n1,n2,n3,n4)<0) THEN
      uu= uu - CONJG(ci(-ciindex(n1,n2,n3,n4))%v)
    END IF
  END SUBROUTINE SUBCI_X

  SUBROUTINE ADDELEMENT_X(nonzerocount, col, hamj, ham, cutoff, val)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: nonzerocount
    INTEGER, INTENT(IN) :: col
    INTEGER, INTENT(INOUT) :: hamj(:)
    COMPLEX*16, INTENT(INOUT) :: ham(:)
    REAL*8, INTENT(IN) :: cutoff
    COMPLEX*16, INTENT(IN) :: val
    IF ( ABS(val) >= cutoff ) THEN
      nonzerocount= nonzerocount + 1
      hamj(nonzerocount)= col
      ham(nonzerocount)= val
    END IF
  END SUBROUTINE ADDELEMENT_X

END SUBROUTINE CREATEHAMILTONIAN_X



!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

SUBROUTINE CREATEHAMILTONIAN( dimhspace, ket, cutoff,              &
     &  numci_ee, ci_ee, ciindex_ee, numci_hh, ci_hh, ciindex_hh,  &
     &  numci_eh, ci_eh, ciindex_eh, spenergy_e, spenergy_h,       &
     &  nonzero, hami, hamj, ham )
  USE mod_indata
  IMPLICIT NONE
! creates the REAL Hamiltonian matrix H= Te + Th + Vee + Vhh + Veh
! where the V's are REAL Coulomb interaction terms
! it used the Compressed Row Storage format (same as Pardiso)
! see  http://www.cs.utk.edu/~dongarra/etemplates/node371.html

INTEGER, INTENT(IN) :: dimhspace            ! Hilbert space basis for both elecs and holes
INTEGER*8, INTENT(IN) :: ket(1:,:)  !(dimhspace,2)
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
REAL*8, INTENT(IN) :: spenergy_e(:)  !(numspstates_e)
REAL*8, INTENT(IN) :: spenergy_h(:)  !(numspstates_h)
INTEGER, INTENT(IN) :: nonzero
INTEGER, INTENT(OUT) :: hami(:)
INTEGER, INTENT(OUT) :: hamj(:)
REAL*8, INTENT(OUT) :: ham(:)

INTEGER :: row, col
INTEGER :: nonzerocount
INTEGER*8 :: brae, brah, kete, keth
INTEGER*8 :: braebit, brahbit, ketebit, kethbit, ketbraexor, ketbrahxor
INTEGER :: emoved, hmoved
REAL*8 :: teth, uee, uhh, ueh
REAL*8 :: sig
INTEGER :: ni, nj, nk, nl   ! counters on single part states


!!$
!!$ ** with blockizehamiltonian this is not necessary: kept for historical reason
!!$  
! the max number of non-zero elements is (BC: binomial coefficients)
!maxnonzero= dimham * ( 1                                                  &
!     & + BC(numspstates_h-num_h,1)*num_h                                  &
!     & + BC(numspstates_e-num_e,1)*num_e                                  &
!     & + BC(numspstates_h-num_h,1)*num_h*BC(numspstates_e-num_e,1)*num_e  &
!     & + BC(numspstates_h-num_h,2)*BC(num_h,2)                            &
!     & + BC(numspstates_e-num_e,2)*BC(num_e,2)
!
! the different terms correspond to the following possibilities:
!-----------------------------------------------------------
!   elecs moved    holes moved    H terms /= 0
!-----------------------------------------------------------
!       0              0          Te + Th + Vee + Vhh + Veh
!       0              1          Vhh + Veh
!       1              0          Vee + Veh
!       1              1          Veh
!       0              2          Vhh
!       2              0          Vee
!
! IF ( MOD(maxnonzero-dimham,2)>0 ) STOP "CREATEHAMILTONIAN: odd nonzeros??!"
! maxnonzero= (maxnonzero-dimham)/2 + dimham
! since H is hermitian, we store only the diag and the upper triangle,
! i.e. maxnonzero only includes the upper triangle elements
!  -- the above calculation is performed in the main parogram --

IF (cutoff < 0d0) STOP "CREATEHAMILTONIAN: negative cutoff"
nonzerocount= 0
DO row= 1, dimhspace
  ! before the reordering the total basis was
  ! the total basis is |col>=|nke,nkh> with col=(nke-1)Nh+nkh
  ! i.e. in the row/col index the hole index runs faster (same for rows)

  ! bra for elects
  brae= ket( row, 1 )
  ! bra for holes
  brah= ket( row, 2 )

  ! diagonal element
  teth= 0.
  DO ni= 1, numspstates_e
    IF (BTEST(brae,ni-1))  teth= teth + spenergy_e(ni)
  END DO
  DO nj= 1, numspstates_h
    IF (BTEST(brah,nj-1))  teth= teth + spenergy_h(nj)
  END DO

  uee= 0.
  DO ni= 1, numspstates_e-1
    IF (BTEST(brae,ni-1)) THEN
      DO nj= ni+1, numspstates_e
        IF (BTEST(brae,nj-1)) THEN
          CALL ADDCI(uee, ci_ee, ciindex_ee, ni, nj, nj, ni)
          CALL SUBCI(uee, ci_ee, ciindex_ee, nj, ni, nj, ni)
        END IF
      END DO
    END IF
  END DO

  uhh= 0.
  DO ni= 1, numspstates_h-1
    IF (BTEST(brah,ni-1)) THEN
      DO nj= ni+1, numspstates_h
        IF (BTEST(brah,nj-1)) THEN
          CALL ADDCI(uhh, ci_hh, ciindex_hh, ni, nj, nj, ni)
          CALL SUBCI(uhh, ci_hh, ciindex_hh, nj, ni, nj, ni)
        END IF
      END DO
    END IF
  END DO

  ueh= 0.
  DO ni= 1, numspstates_h
    IF (BTEST(brah,ni-1)) THEN
      DO nj= 1, numspstates_e
        IF (BTEST(brae,nj-1)) THEN
          CALL ADDCI(ueh, ci_eh, ciindex_eh, ni, nj, nj, ni)
        END IF
      END DO
    END IF
  END DO

  ! first elem of a new row, here row=col
  hami(row)= nonzerocount+1
  ! diag element always present, also if it is zero
  CALL  ADDELEMENT(nonzerocount, row, hamj, ham, -1d0, teth + uee + uhh + ueh)

  ! calculate following nonzero elems in same row in the upper triangle
  DO col= row+1, dimhspace
    ! ket for elects
    kete= ket( col, 1 )
    ! ket for holes
    keth= ket( col, 2 )
    ! PRINT*, row, nbe, nbh
    ! PRINT*, col, nke, nkh
    ! PRINT*
    ketbraexor= IEOR(kete,brae)
    emoved= POPCNT( ketbraexor )
    ketbrahxor= IEOR(keth,brah)
    hmoved= POPCNT( ketbrahxor )

!!$    IF (row==2 .AND. col==9) THEN
!!$      PRINT*, kete,brae,keth,brah
!!$      PRINT*, ketbraexor, emoved
!!$      PRINT*, ketbrahxor, hmoved
!!$      STOP "sssssssssssssssss"
!!$    END IF

!!$    IF (emoved==0 .AND. hmoved==0) THEN
!!$      npm00= npm00 + 1
!!$    ELSE IF (emoved==0 .AND. hmoved==2) THEN
!!$      npm01= npm01 + 1
!!$    ELSE IF (emoved==2 .AND. hmoved==0) THEN
!!$      npm10= npm10 + 1
!!$    ELSE IF (emoved==2 .AND. hmoved==2) THEN
!!$      npm11= npm11 + 1
!!$    ELSE IF (emoved==0 .AND. hmoved==4) THEN
!!$      npm02= npm02 + 1
!!$    ELSE IF (emoved==4 .AND. hmoved==0) THEN
!!$      npm20= npm20 + 1
!!$    ELSE
!!$      npmx= npmx + 1
!!$    END IF

    IF (emoved==0 .AND. hmoved==2) THEN  !   *** 1 hole is moved ***

      kethbit= IAND(ketbrahxor,keth)     ! from ni position
      DO ni= 1, numspstates_h
        IF ( BTEST(kethbit, ni-1) ) EXIT
      END DO
      brahbit= IAND(ketbrahxor,brah)     ! to nj position
      DO nj= 1, numspstates_h
        IF ( BTEST(brahbit, nj-1) ) EXIT
      END DO
      ! sign according to num of 1's between ni and nj
      sig= (-1)**(POPCNT(IBITS( keth, MIN(ni,nj), ABS(ni-nj)-1 )))

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
      uhh= sig * uhh

      ueh= 0.
      DO nk= 1, numspstates_e
        IF ( BTEST(kete, nk-1) ) THEN
          CALL ADDCI(ueh, ci_eh, ciindex_eh, nj, nk, nk, ni)
        END IF
      END DO
      ueh= sig * ueh

      CALL  ADDELEMENT(nonzerocount, col, hamj, ham, cutoff, uhh + ueh)

    ELSE IF (emoved==2 .AND. hmoved==0) THEN  ! *** 1 electron is moved ***

      ketebit= IAND(ketbraexor,kete)     ! from ni position
      DO ni= 1, numspstates_e
        IF ( BTEST(ketebit, ni-1) ) EXIT
      END DO
      braebit= IAND(ketbraexor,brae)     ! to nj position
      DO nj= 1, numspstates_e
        IF ( BTEST(braebit, nj-1) ) EXIT
      END DO
      ! sign according to num of 1's between ni and nj
      sig= (-1)**(POPCNT(IBITS( kete, MIN(ni,nj), ABS(ni-nj)-1 )))

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
      uee= sig * uee

      ueh= 0.
      DO nk= 1, numspstates_h
        IF ( BTEST(keth, nk-1) ) THEN
          CALL ADDCI(ueh, ci_eh, ciindex_eh, nk, nj, ni, nk)
        END IF
      END DO
      ueh= sig * ueh

      CALL  ADDELEMENT(nonzerocount, col, hamj, ham, cutoff, uee + ueh)

    ELSE IF (emoved==2 .AND. hmoved==2) THEN  !*** 1 hole and 1 elec moved ***
      
      kethbit= IAND(ketbrahxor,keth)     ! hole moved from nl position
      DO nl= 1, numspstates_h
        IF ( BTEST(kethbit, nl-1) ) EXIT
      END DO
      brahbit= IAND(ketbrahxor,brah)     !              to ni position
      DO ni= 1, numspstates_h
        IF ( BTEST(brahbit, ni-1) ) EXIT
      END DO
      ! sign according to num of 1's between nl and ni
      sig= (-1)**(POPCNT(IBITS( keth, MIN(nl,ni), ABS(nl-ni)-1 )))

      ketebit= IAND(ketbraexor,kete)     ! elec moved from nk position
      DO nk= 1, numspstates_e
        IF ( BTEST(ketebit, nk-1) ) EXIT
      END DO
      braebit= IAND(ketbraexor,brae)     ! to nj position
      DO nj= 1, numspstates_e
        IF ( BTEST(braebit, nj-1) ) EXIT
      END DO
      ! sign according to num of 1's between ni and nj
      sig= sig * (-1)**(POPCNT(IBITS( kete, MIN(nk,nj), ABS(nk-nj)-1 )))

      ueh= 0.
      CALL ADDCI(ueh, ci_eh, ciindex_eh, ni, nj, nk, nl)
      ueh= sig * ueh

      CALL  ADDELEMENT(nonzerocount, col, hamj, ham, cutoff, ueh)

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

      ! sign according to num of 1's between ni & nj + nk & nl 
      ! sign according to num of 1's between ni & nj + nk & nl 
      sig= (-1)**(POPCNT(IBITS( keth, nl, nk-nl-1 )))
      sig= -sig * (-1)**(POPCNT(IBITS( brah, nj, ni-nj-1 )))

      uhh= 0.
      CALL ADDCI(uhh, ci_hh, ciindex_hh, ni, nj, nk, nl)
      CALL SUBCI(uhh, ci_hh, ciindex_hh, nj, ni, nk, nl)
      uhh= sig * uhh

      CALL  ADDELEMENT(nonzerocount, col, hamj, ham, cutoff, uhh)

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

!!$      IF (row==2 .AND. col==9) THEN
!!$        PRINT*, nl, nk
!!$        PRINT*, nj, ni
!!$        STOP "ttttt 1"
!!$      END IF

      ! sign according to num of 1's between ni & nj + nk & nl 
      sig= (-1)**(POPCNT(IBITS( kete, nl, nk-nl-1 )))
      sig= -sig * (-1)**(POPCNT(IBITS( brae, nj, ni-nj-1 )))

      uee= 0.
      CALL ADDCI(uee, ci_ee, ciindex_ee, ni, nj, nk, nl)
      CALL SUBCI(uee, ci_ee, ciindex_ee, nj, ni, nk, nl)
      uee= sig * uee

      CALL  ADDELEMENT(nonzerocount, col, hamj, ham, cutoff, uee)
      
    END IF

  END DO ! col

END DO ! row

hami(dimhspace+1)= nonzerocount + 1

print*, "aaaaaaaaaaaaaaa", nonzerocount, nonzero
IF ( nonzerocount /= nonzero ) STOP "CREATEHAMILTONIAN : nonzerocount /= nonzero"

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

  SUBROUTINE ADDELEMENT(nonzerocount, col, hamj, ham, cutoff, val)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: nonzerocount
    INTEGER, INTENT(IN) :: col
    INTEGER, INTENT(INOUT) :: hamj(:)
    REAL*8, INTENT(INOUT) :: ham(:)
    REAL*8, INTENT(IN) :: cutoff
    REAL*8, INTENT(IN) :: val
    IF ( ABS(val) >= cutoff ) THEN
      nonzerocount= nonzerocount + 1
      hamj(nonzerocount)= col
      ham(nonzerocount)= val
    END IF
  END SUBROUTINE ADDELEMENT

END SUBROUTINE CREATEHAMILTONIAN



END MODULE mod_createhamiltonian
