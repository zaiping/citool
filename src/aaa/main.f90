!$$$$$$$$$$$$$$$$$$$$$$$$$  CItool  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

PROGRAM MAIN
  USE mod_indata
  USE mod_logga
  USE mod_myaux
  USE mod_specialf
  USE mod_cimakeindex
  USE mod_createhspace
  USE mod_blockizehamiltonian
  USE mod_createhamiltonian
  USE mod_diagonalize
  USE mod_writempstates
  IMPLICIT NONE
!*************************************************************************
!*  myci: FewBodyDiag finds few-body eigenvalues and eigenvectors
!*  by diagonalizing the Hamiltonian on a basis composed by
!*  Slater determinants
!*  Two carrier types: e and h are considered.
!*  An arbitrary number of single-particle quantum numbers can
!*  be used.
!*  The Coulomb integrals and the single-particle energies (and q numbers)
!*  are read from five external files (sps_e, sps_h, U_ee, U_hh, U_eh)
!*  now uses Arpack libaries for the diagonalization
!*
!*  11 sept 2008 - 20 apr 2010 
!*
!*  Intel Fortran compiler                               by Andrea Bertoni
!*************************************************************************

!............................................................Declarations
! single-particle states
CHARACTER(LEN=12), ALLOCATABLE :: namespqn_e(:)
CHARACTER(LEN=12), ALLOCATABLE :: namespqn_h(:)
INTEGER, ALLOCATABLE :: spqn_e(:,:)
INTEGER, ALLOCATABLE :: spqn_h(:,:)
REAL*8, ALLOCATABLE :: spenergy_e(:)
REAL*8, ALLOCATABLE :: spenergy_h(:)
!!$INTEGER, ALLOCATABLE :: sumspqn_e(:)
!!$INTEGER, ALLOCATABLE :: sumspqn_h(:)
!!$INTEGER, ALLOCATABLE :: prevsumspqn_e(:)
!!$INTEGER, ALLOCATABLE :: prevsumspqn_h(:)

! Coulomb integrals
INTEGER :: numci_ee=0, numci_hh=0, numci_eh=0
TYPE( ci_type_complex16 ), ALLOCATABLE :: ci_ee_x(:)
TYPE( ci_type_complex16 ), ALLOCATABLE :: ci_hh_x(:)
TYPE( ci_type_complex16 ), ALLOCATABLE :: ci_eh_x(:)
TYPE( ci_type_real8 ), ALLOCATABLE :: ci_ee(:)
TYPE( ci_type_real8 ), ALLOCATABLE :: ci_hh(:)
TYPE( ci_type_real8 ), ALLOCATABLE :: ci_eh(:)
INTEGER, ALLOCATABLE :: ciindex_ee(:,:,:,:)
INTEGER, ALLOCATABLE :: ciindex_hh(:,:,:,:)
INTEGER, ALLOCATABLE :: ciindex_eh(:,:,:,:)

! multi-electron states
INTEGER :: dimhspace_e, dimhspace_h, dimhspace
INTEGER*8, ALLOCATABLE :: ket_e(:)          ! Slater dets for elecs
INTEGER*8, ALLOCATABLE :: ket_h(:)          ! Slater dets for holes
INTEGER*8, ALLOCATABLE :: ket(:,:)          ! Slater dets for both 

! Hamiltonian in Compressed Row Storage format (same as Pardiso)
! see  http://www.cs.utk.edu/~dongarra/etemplates/node371.html
INTEGER :: numblock
INTEGER :: blockdim, blockfr, blockto
INTEGER :: maxblockdim, maxblocknonzero
INTEGER :: blocknummpenergies, blocknummpstates
INTEGER, ALLOCATABLE :: blockstart(:)
INTEGER, ALLOCATABLE :: blocknonzero(:)
INTEGER, ALLOCATABLE :: hami(:)
INTEGER, ALLOCATABLE :: hamj(:)
COMPLEX*16, ALLOCATABLE :: ham_x(:)
REAL*8, ALLOCATABLE :: ham(:)
COMPLEX*16, ALLOCATABLE :: mpstates_x(:,:)
REAL*8, ALLOCATABLE :: mpstates(:,:)
REAL*8, ALLOCATABLE :: mpenergies(:,:)

! aux vars
INTEGER :: nn, nb, ni, nj


!..........................................init vars definitions and readout
CALL INDATA_GET("citool.nml")
CALL LOGGA(3, " == START ==")

ALLOCATE(namespqn_e(numspqn_e))
ALLOCATE(namespqn_h(numspqn_h))
ALLOCATE(spqn_e(numspstates_e,numspqn_e))
ALLOCATE(spqn_h(numspstates_h,numspqn_h))
ALLOCATE(spenergy_e(numspstates_e))
ALLOCATE(spenergy_h(numspstates_h))
!!$ALLOCATE(sumspqn_e(numspqn_e))
!!$ALLOCATE(sumspqn_h(numspqn_h))
!!$ALLOCATE(prevsumspqn_e(numspqn_e))
!!$ALLOCATE(prevsumspqn_h(numspqn_h))

!...............................................reads single-particle states
CALL LOGGA(2, "number of ELECTRONS: ", num_e)
IF (num_e>0) THEN
  CALL LOGGA(2, "reads single particle elec states")
  CALL INDATA_SPSTATES( "e", numspstates_e, numspqn_e,    &
       &  namespqn_e, spqn_e, spenergy_e )
  !!$PRINT*, numspstates_e, numspqn_e
  !!$PRINT*, namespqn_e, spqn_e
  !!$PRINT*, spenergy_e
  !!$PRINT*
ELSE
  CALL LOGGA(2, "fictious single elec states considered:", numspstates_e)
END IF

CALL LOGGA(2, "number of HOLES: ", num_h)
IF (num_h>0) THEN
  CALL LOGGA(2, "reads single particle hole states")
  CALL INDATA_SPSTATES( "h", numspstates_h, numspqn_h,    &
       &  namespqn_h, spqn_h, spenergy_h )
  !!$PRINT*, numspstates_h, numspqn_h
  !!$PRINT*, namespqn_h, spqn_h
  !!$PRINT*, spenergy_h
  !!$PRINT*
ELSE
  CALL LOGGA(2, "fictious single hole states considered:", numspstates_h)
END IF

!..................reads Coulomb integrals file and builds their indexes
CALL LOGGA(2, "reading and indexing Coulomb integrals...")
! ci and index arrays ALLOCATEd inside the subroutines

IF ( complexrun ) THEN
  IF ( INDEX(fileinformat_coulomb,"real8") /= 0 ) THEN
    CALL LOGGA(2, "WARNING: complex run with real Coulomb integrals")
    CALL LOGGA(2, "         the imag part will be set to zero !")
  END IF
  IF (num_e>1) THEN
    CALL INDATA_COULOMB_X( "ee", numci_ee, ci_ee_x )
    CALL CIMAKEINDEX_X( "ee", numci_ee, ci_ee_x, ciindex_ee )
    CALL LOGGA(2, "  reads complex ee Coulomb integrals")
  ELSE
    numci_ee= 1
    ALLOCATE( ci_ee_x(numci_ee) )
    ALLOCATE( ciindex_ee(numspstates_e,numspstates_e,numspstates_e,numspstates_e) )
    ciindex_ee(:,:,:,:)= 0
  END IF
  IF (num_h>1) THEN
    CALL INDATA_COULOMB_X( "hh", numci_hh, ci_hh_x )
    CALL CIMAKEINDEX_X( "hh", numci_hh, ci_hh_x, ciindex_hh )
    CALL LOGGA(2, "  reads complex hh Coulomb integrals")
  ELSE
    numci_hh= 1
    ALLOCATE( ci_hh_x(numci_hh) )
    ALLOCATE( ciindex_hh(numspstates_h,numspstates_h,numspstates_h,numspstates_h) )
    ciindex_hh(:,:,:,:)= 0
  END IF
  IF (num_e>0 .AND. num_h>0) THEN
    CALL INDATA_COULOMB_X( "eh", numci_eh, ci_eh_x )
    CALL CIMAKEINDEX_X( "eh", numci_eh, ci_eh_x, ciindex_eh )
    CALL LOGGA(2, "  reads complex eh Coulomb integrals")
  ELSE
    numci_eh= 1
    ALLOCATE( ci_eh_x(numci_eh) )
    ALLOCATE( ciindex_eh(numspstates_h,numspstates_e,numspstates_e,numspstates_h) )
    ciindex_eh(:,:,:,:)= 0
  END IF
ELSE
  IF ( INDEX(fileinformat_coulomb,"complex16") /= 0 ) THEN
    CALL LOGGA(2, "WARNING: real run with complex Coulomb integrals")
    CALL LOGGA(2, "         the imag part will be ignored !")
  END IF
  IF (num_e>1) THEN
    CALL INDATA_COULOMB( "ee", numci_ee, ci_ee )
    CALL CIMAKEINDEX( "ee", numci_ee, ci_ee, ciindex_ee )
    CALL LOGGA(2, "  reads real ee Coulomb integrals")
  ELSE
    numci_ee= 1
    ALLOCATE( ci_ee(numci_ee) )
    ALLOCATE( ciindex_ee(numspstates_e,numspstates_e,numspstates_e,numspstates_e) )
    ciindex_ee(:,:,:,:)= 0
  END IF
  IF (num_h>1) THEN
    CALL INDATA_COULOMB( "hh", numci_hh, ci_hh )
    CALL CIMAKEINDEX( "hh", numci_hh, ci_hh, ciindex_hh )
    CALL LOGGA(2, "  reads real hh Coulomb integrals")
  ELSE
    numci_hh= 1
    ALLOCATE( ci_hh(numci_hh) )
    ALLOCATE( ciindex_hh(numspstates_h,numspstates_h,numspstates_h,numspstates_h) )
    ciindex_hh(:,:,:,:)= 0
  END IF
  IF (num_e>0 .AND. num_h>0) THEN
    CALL INDATA_COULOMB( "eh", numci_eh, ci_eh )
    CALL CIMAKEINDEX( "eh", numci_eh, ci_eh, ciindex_eh )
    CALL LOGGA(2, "  reads real eh Coulomb integrals")
  ELSE
    numci_eh= 1
    ALLOCATE( ci_eh(numci_eh) )
    ALLOCATE( ciindex_eh(numspstates_h,numspstates_e,numspstates_e,numspstates_h) )
    ciindex_eh(:,:,:,:)= 0
  END IF
END IF

CALL LOGGA(2, "...done")


!PPPPPPPPPPPPPPPPP
!ci_eh(:)%v= ci_eh(:)%v

!!$PRINT*, "ee"
!!$PRINT*, numci_ee
!!$PRINT*, ci_ee
!!$PRINT*
!!$
!!$PRINT*, "hh"
!!$PRINT*, numci_hh
!!$PRINT*, ci_hh
!!$PRINT*
!!$
!!$PRINT*, "eh"
!!$PRINT*, numci_eh
!!$PRINT*, ci_eh
!!$PRINT*

!.................................................creates the Hilbert spaces
dimhspace_e= BINOMIALCO(numspstates_e,num_e)
CALL LOGGA(2, "elecs Hilbert space dimension:", dimhspace_e)
dimhspace_h= BINOMIALCO(numspstates_h,num_h)
CALL LOGGA(2, "holes Hilbert space dimension:", dimhspace_h)
dimhspace= dimhspace_e * dimhspace_h
CALL LOGGA(2, "total Hilbert space dimension:", dimhspace)

ALLOCATE( ket_e(dimhspace_e) )
ALLOCATE( ket_h(dimhspace_h) )
ALLOCATE( ket(dimhspace,2) )

CALL LOGGA(2, "creating elecs Hilbert space...")
CALL CREATEHSPACE( num_e, numspstates_e, dimhspace_e, ket_e )
CALL LOGGA(2, "...done")

CALL LOGGA(2, "creating holes Hilbert space...")
CALL CREATEHSPACE( num_h, numspstates_h, dimhspace_h, ket_h )
CALL LOGGA(2, "...done")

DO nn= 1, dimhspace
  ket(nn,1)= ket_e( (nn-1)/dimhspace_h + 1 )
  ket(nn,2)= ket_h( MOD(nn-1,dimhspace_h) + 1 )
END DO

!.................................................writes Hilbert spaces
WRITE(*,*) "*** elec Hilbert space ***"
DO nn= 1, dimhspace_e
  IF ( POPCNT(ket_e(nn)) /= num_e ) STOP "MAIN: 1 bits /= num_e"
  WRITE(*,binfmt_e) ket_e(nn)
END DO
WRITE(*,*)

WRITE(*,*) "*** hole Hilbert space ***"
DO nn= 1, dimhspace_h
  IF ( POPCNT(ket_h(nn)) /= num_h ) STOP "MAIN: 1 bits /= num_h"
  WRITE(*,binfmt_h) ket_h(nn)
END DO
WRITE(*,*)

!......................................basis reordering for block Hamiltonian
CALL LOGGA(2, "reordering the basis for block Hamiltonian...")

ALLOCATE( blockstart(dimhspace+1) )
ALLOCATE( blocknonzero(dimhspace) )

CALL BLOCKIZEHAMILTONIAN( dimhspace, ket, ciindex_ee, ciindex_hh, ciindex_eh,  &
     &                    numblock, blockstart, blocknonzero )

!!$numblock= 1
!!$blockstart(1)= 1
!!$blockstart(2)= dimhspace+1
!!$blockstart(3:)= 0
!!$blocknonzero(1)= 475343919
!!$
!!$print*, " ZZZ "
!!$print*, numblock
!!$print*, blockstart(1:7)
!!$print*, " ZZZ "
!!$!stop

!!$WRITE(44,*) numblock
!!$WRITE(44,*) blockstart
!!$WRITE(44,*) blocknonzero
!!$WRITE(44,*) 

CALL LOGGA(2, "...done")
CALL LOGGA(2, "number of blocks found= ", numblock)
CALL LOGGA(2, "total number of non-zero elements= ", SUM(blocknonzero))

!.................................................writes reordered Hilbert spaces
WRITE(*,*) "*** total reordered Hilbert space ***"

DO nb= 1, numblock
  WRITE(*,*) "  * * * * *  BLOCK", nb
  DO nn= blockstart(nb), blockstart(nb+1)-1
!!$    prevsumspqn_e(:)= sumspqn_e(:)
!!$    sumspqn_e(:)= 0
!!$    DO ni= 1, numspstates_e
!!$      IF (BTEST(ket(nn,1),ni-1))  sumspqn_e(:)= sumspqn_e(:) + spqn_e(ni,:)
!!$    END DO
!!$    prevsumspqn_h(:)= sumspqn_h(:)
!!$    sumspqn_h(:)= 0
!!$    DO nj= 1, numspstates_h
!!$      IF (BTEST(ket(nn,2),nj-1))  sumspqn_h(:)= sumspqn_h(:) + spqn_h(nj,:)
!!$    END DO
!!$    IF ( ANY(sumspqn_e(1:2)/=prevsumspqn_e(1:2)) .OR. ANY(sumspqn_h(1:2)/=prevsumspqn_h(1:2)) ) THEN
!!$      WRITE(*,*) sumspqn_e
!!$      WRITE(*,*) sumspqn_h
!!$    END IF
    WRITE(*,binfmt_eh) ket(nn,1), ket(nn,2)

  END DO
END DO
WRITE(*,*)

!......................................finds maximum dim and nonzeros of blocks
maxblockdim= 0
maxblocknonzero= 0
DO nb= 1, numblock
  maxblockdim= MAX( maxblockdim, blockstart(nb+1)-blockstart(nb) )
  maxblocknonzero= MAX( maxblocknonzero, blocknonzero(nb) )
END DO
CALL LOGGA(2, "max dimension of a block= ", maxblockdim)
CALL LOGGA(2, "max number of non-zero elements= ", maxblocknonzero)

!.................................................allocating block Hamiltonian
CALL LOGGA(2, "allocating the block Hamiltonian...")
ALLOCATE( hami(maxblockdim+1) )
ALLOCATE( hamj(maxblocknonzero) )
IF ( complexrun ) THEN
  ALLOCATE( ham_x(maxblocknonzero) )
  CALL LOGGA(2, "complex Hamiltonian storage (Mb)=~", maxblocknonzero*20/1e6)
ELSE
  ALLOCATE( ham(maxblocknonzero) )
  CALL LOGGA(2, "real Hamiltonian storage (Mb)=~", maxblocknonzero*12/1e6)
END IF
CALL LOGGA(2, "...done")

!.................................................allocating solution arrays
CALL LOGGA(2, "allocating eigenvectors/values arrays...")
ALLOCATE( mpenergies(nummpenergies,numblock) )
IF ( complexrun ) THEN
  ALLOCATE( mpstates_x(dimhspace,nummpstates) )
  CALL LOGGA(2, "complex mp states storage (Mb)=~", dimhspace*nummpstates*16/1e6)
ELSE
  ALLOCATE( mpstates(dimhspace,nummpstates) )
  CALL LOGGA(2, "real mp states storage (Mb)=~", dimhspace*nummpstates*8/1e6)
END IF
CALL LOGGA(2, "...done")

! hamj(400000000:)= 5
! ham(400000000:)= 5.2
! mpstates(444:555,:)= 7.3


!!$
!!$ ** with blockizehamiltonian this is not necessary: kept for historical reason
!!$  
!!$dimham= dimhspace_e * dimhspace_h
!!$maxnonzero= dimham * ( 1                                                 &
!!$     & + (numspstates_h-num_h)*num_h                                     &
!!$     & + (numspstates_e-num_e)*num_e                                     &
!!$     & + (numspstates_h-num_h)*num_h*(numspstates_e-num_e)*num_e         &
!!$     & + (numspstates_h-num_h)*num_h*(numspstates_h-num_h-1)*(num_h-1)/4 &
!!$     & + (numspstates_e-num_e)*num_e*(numspstates_e-num_e-1)*(num_e-1)/4 )
!!$!!* this should be equiv to the expression with binomial coefficients
!!$!
!!$!maxnonzero= dimham * ( 1                                                  &
!!$!     & + BC(numspstates_h-num_h,1)*num_h                                  &
!!$!     & + BC(numspstates_e-num_e,1)*num_e                                  &
!!$!     & + BC(numspstates_h-num_h,1)*num_h*BC(numspstates_e-num_e,1)*num_e  &
!!$!     & + BC(numspstates_h-num_h,2)*BC(num_h,2)                            &
!!$!     & + BC(numspstates_e-num_e,2)*BC(num_e,2)
!!$!
!!$! the different terms correspond to the following possibilities:
!!$!-----------------------------------------------------------
!!$!   elecs moved    holes moved    H terms /= 0
!!$!-----------------------------------------------------------
!!$!       0              0          Te + Th + Vee + Vhh + Veh
!!$!       0              1          Vhh + Veh
!!$!       1              0          Vee + Veh
!!$!       1              1          Veh
!!$!       0              2          Vhh
!!$!       2              0          Vee
!!$!
!!$IF ( MOD(maxnonzero-dimham,2)>0 ) STOP "MAIN: odd maxnonzeros??!"
!!$maxnonzero= (maxnonzero-dimham)/2 + dimham
!!$! H stored using the Compressed Row Storage format (same as Pardiso)
!!$! see  http://www.cs.utk.edu/~dongarra/etemplates/node371.html
!!$! since H is hermitian, we store only the diag and the upper triangle,
!!$! i.e. maxnonzero only includes the upper triangle elements
!!$
!!$ALLOCATE( hami(dimham+1) )
!!$ALLOCATE( hamj(maxnonzero) )


!..................................................begin loop over Ham blocks
DO nb= 1, numblock
  CALL LOGGA(2, "  block n. ", nb)

  !............................creates the Hamiltonian matrix for current block
  blockdim= blockstart(nb+1)-blockstart(nb)
  blockfr= blockstart(nb)
  blockto= blockstart(nb+1)-1
  blocknummpenergies= MIN(nummpenergies,blockdim)
  blocknummpstates= MIN(nummpstates,blockdim)

  CALL LOGGA(2, "current block H dimension= ", blockdim)
  CALL LOGGA(2, "current number of non-zero elements= ", blocknonzero(nb))
  CALL LOGGA(2, "creating the block H matrix...")

  IF ( complexrun ) THEN
    CALL CREATEHAMILTONIAN_X( blockdim, ket(blockfr:blockto,1:2),  &
         &  numci_ee, ci_ee_x, ciindex_ee, numci_hh, ci_hh_x, ciindex_hh,      &
         &  numci_eh, ci_eh_x, ciindex_eh, spenergy_e, spenergy_h,             &
         &  blocknonzero(nb), hami, hamj, ham_x )
  ELSE
    CALL CREATEHAMILTONIAN( blockdim, ket(blockfr:blockto,1:2),  &
         &  numci_ee, ci_ee, ciindex_ee, numci_hh, ci_hh, ciindex_hh,        &
         &  numci_eh, ci_eh, ciindex_eh, spenergy_e, spenergy_h,             &
         &  blocknonzero(nb), hami, hamj, ham )
  END IF

  CALL LOGGA(2, "...done")

!!$  OPEN(33, FILE="hamiltonianIJH.bin", POSITION="APPEND", FORM="UNFORMATTED")
!!$  WRITE(33) hami
!!$  WRITE(33) hamj
!!$  WRITE(33) ham
!!$  CLOSE(33)

!!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
!DO nn= 1, dimham !dimham
!  DO mm= hami(nn), hami(nn+1)-1
!    IF ( ABS( ham(mm) ) > 1e-22 ) THEN
!      PRINT*, nn, hamj(mm), ham(mm)
!    END IF
!  END DO
!  PRINT*
!END DO
!stop "ddddddddddddddd"
!!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP

  !..................................finds block Hamiltonian evals and evects
  CALL LOGGA(2, "diagonalizing the block Hamiltonian...")

  IF ( complexrun ) THEN
    CALL DIAGONALIZE_X( blockdim, hami, hamj, ham_x, blocknummpenergies,       &
         &  mpenergies(1:blocknummpenergies,nb), blocknummpstates,             &
         &  mpstates_x(blockfr:blockto,1:blocknummpstates) )
  ELSE
    CALL DIAGONALIZE( blockdim, hami, hamj, ham, blocknummpenergies,         &
         &  mpenergies(1:blocknummpenergies,nb), blocknummpstates,           &
         &  mpstates(blockfr:blockto,1:blocknummpstates) )
  END IF

  CALL LOGGA(2, "...done")

!!$  OPEN(33, FILE="evalsvecs", POSITION="APPEND", FORM="UNFORMATTED")
!!$  WRITE(33) mpenergies(1:blocknummpenergies,nb)
!!$  WRITE(33) mpstates(blockfr:blockto,1:blocknummpstates)
!!$  CLOSE(33)

  !.........................................writes mp energies and states
  CALL LOGGA(2, "writing multiparticle states...")

  print*, "writing mp states of BLOCK", nb
  print*, "blockfr, blockto=", blockfr, blockto
  IF ( complexrun ) THEN
    CALL WRITEMPSTATES_X( nb, blockdim, ket(blockfr:blockto,1:2),         &
         &  namespqn_e, spqn_e, namespqn_h, spqn_h,                     &
         &  blocknummpenergies, mpenergies(1:blocknummpenergies,nb),    &
         &  blocknummpstates, mpstates_x(blockfr:blockto,1:blocknummpstates) )
  ELSE
    CALL WRITEMPSTATES( nb, blockdim, ket(blockfr:blockto,1:2),         &
         &  namespqn_e, spqn_e, namespqn_h, spqn_h,                   &
         &  blocknummpenergies, mpenergies(1:blocknummpenergies,nb),  &
         &  blocknummpstates, mpstates(blockfr:blockto,1:blocknummpstates) )
  END IF
  
  CALL LOGGA(2, "...done")


END DO   ! end loop on the blocks


!.........................................................deallocate arrays
CALL LOGGA(2, "deallocating ci arrays...")
IF ( complexrun ) THEN
  DEALLOCATE(ci_ee_x)
  DEALLOCATE(ciindex_ee)
  DEALLOCATE(ci_hh_x)
  DEALLOCATE(ciindex_hh)
  DEALLOCATE(ci_eh_x)
  DEALLOCATE(ciindex_eh)
ELSE
  DEALLOCATE(ci_ee)
  DEALLOCATE(ciindex_ee)
  DEALLOCATE(ci_hh)
  DEALLOCATE(ciindex_hh)
  DEALLOCATE(ci_eh)
  DEALLOCATE(ciindex_eh)
END IF
CALL LOGGA(2, "...done")

CALL LOGGA(2, "deallocating remaining arrays...")
DEALLOCATE( namespqn_e )
DEALLOCATE( namespqn_h )
DEALLOCATE( spqn_e )
DEALLOCATE( spqn_h )
DEALLOCATE( spenergy_e )
DEALLOCATE( spenergy_h )
DEALLOCATE( ket_e )
DEALLOCATE( ket_h )
DEALLOCATE( ket )
DEALLOCATE( blockstart )
DEALLOCATE( blocknonzero )
DEALLOCATE( hami )
DEALLOCATE( hamj )
DEALLOCATE( mpenergies )
IF ( complexrun ) THEN
  DEALLOCATE( ham_x )
  DEALLOCATE( mpstates_x )
ELSE
  DEALLOCATE( ham )
  DEALLOCATE( mpstates )
END IF
CALL LOGGA(2, "...done")

CALL LOGGA(3, " == END ==")

END PROGRAM MAIN
