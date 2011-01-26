!$$$$$$$$$$$$$$$$$$$$$$$$$  CItool  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

PROGRAM MAIN
  USE mod_indata
  USE mod_logga
  USE mod_myaux
  USE mod_specialf
  USE mod_cimakeindex
  USE mod_createhspace
  USE mod_composehspace
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

! single-particle states (allocated in indata_spstates)
INTEGER :: numspqn_e, numspqn_h
CHARACTER(LEN=12), ALLOCATABLE :: namespqn_e(:)
CHARACTER(LEN=12), ALLOCATABLE :: namespqn_h(:)
CHARACTER(LEN=12), ALLOCATABLE :: typespqn_e(:)
CHARACTER(LEN=12), ALLOCATABLE :: typespqn_h(:)
INTEGER, ALLOCATABLE :: spqn_e(:,:)
INTEGER, ALLOCATABLE :: spqn_h(:,:)
REAL*8, ALLOCATABLE :: spenergy_e(:)
REAL*8, ALLOCATABLE :: spenergy_h(:)
!!$INTEGER, ALLOCATABLE :: sumspqn_e(:)
!!$INTEGER, ALLOCATABLE :: sumspqn_h(:)
!!$INTEGER, ALLOCATABLE :: prevsumspqn_e(:)
!!$INTEGER, ALLOCATABLE :: prevsumspqn_h(:)

! Coulomb integrals (allocated in indata_coulomb)
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

! basis
INTEGER :: dimhspace_e, dimhspace_h, dimhspaceglobal
INTEGER*8, ALLOCATABLE :: ket_e(:)          ! Slater dets for elecs
INTEGER*8, ALLOCATABLE :: ket_h(:)          ! Slater dets for holes
INTEGER*8, ALLOCATABLE :: ket(:,:)          ! Slater dets for both 

! constrains on H (allocated in indata_blockconstrains)
INTEGER :: numhcons_e=0, numhcons_h=0
INTEGER, ALLOCATABLE :: hcons_e(:,:)
INTEGER, ALLOCATABLE :: hcons_h(:,:)
INTEGER :: dimhspacecons
INTEGER :: nummpenergiescons, nummpstatescons

! blocks
INTEGER :: numblock
INTEGER :: blockdim, blockfr, blockto
INTEGER :: maxblockdim, maxblocknonzero
INTEGER, ALLOCATABLE :: blocknummpenergies(:)
INTEGER, ALLOCATABLE :: blocknummpstates(:)
INTEGER, ALLOCATABLE :: blockstart(:)
INTEGER, ALLOCATABLE :: blocknonzero(:)

! Hamiltonian of a block
! stored in Compressed Row Storage format (same as Pardiso)
! see  http://www.cs.utk.edu/~dongarra/etemplates/node371.html
INTEGER, ALLOCATABLE :: hami(:)
INTEGER, ALLOCATABLE :: hamj(:)
COMPLEX*16, ALLOCATABLE :: ham_x(:)
REAL*8, ALLOCATABLE :: ham(:)

! mp states
COMPLEX*16, ALLOCATABLE :: mpstates_x(:,:)
REAL*8, ALLOCATABLE :: mpstates(:,:)
REAL*8, ALLOCATABLE :: mpenergies(:,:)

! aux vars
CHARACTER(80), PARAMETER :: citool_version= "0.9"
REAL*8 :: tinye
! counters
INTEGER :: nn, nb, ncons_e, ncons_h
! debug vars
!INTEGER :: ni, nj
! tmpppppppppppppppppppppppp
!INTEGER :: n1, n2, n3, n4

!..........................................init vars definitions and readout
tinye= TINY(1d1)

CALL INDATA_GET("citool.nml")
CALL LOGGA(2, " ")
CALL LOGGA(3, "  ====  START  CItool v. "//TRIM(citool_version)//"  ====")
IF (TRIM(citoolnml_version) /= TRIM(citool_version)) THEN
  CALL LOGGA(3, "WARNING: input namelist version is "//TRIM(citoolnml_version))
END IF
!!$ALLOCATE(sumspqn_e(numspqn_e))
!!$ALLOCATE(sumspqn_h(numspqn_h))
!!$ALLOCATE(prevsumspqn_e(numspqn_e))
!!$ALLOCATE(prevsumspqn_h(numspqn_h))

!.................................................reads single-particle states
CALL LOGGA(2, "number of ELECTRONS: ", num_e)
IF (num_e>0) THEN
  CALL LOGGA(2, "reads single particle ELEC states", numspstates_e)
  CALL INDATA_SPSTATES( "e", numspstates_e, numspqn_e, namespqn_e, typespqn_e,  &
       &  spqn_e, spenergy_e )
ELSE
  numspstates_e= 1
  numspqn_e= 1
  ALLOCATE(namespqn_e(numspqn_e))
  namespqn_e(1)= "fictiousEqn"
  ALLOCATE(typespqn_e(numspqn_e))
  typespqn_e(1)= "fictiousEty"
  ALLOCATE(spqn_e(numspstates_e,numspqn_e))
  spqn_e(1,1)= 0
  ALLOCATE(spenergy_e(numspstates_e))
  spenergy_e(1)= 0.
  CALL LOGGA(2, "fictious single part ELEC states considered:", numspstates_e)
END IF

CALL LOGGA(2, "number of HOLES: ", num_h)
IF (num_h>0) THEN
  CALL LOGGA(2, "reads single particle HOLE states", numspstates_h)
  CALL INDATA_SPSTATES( "h", numspstates_h, numspqn_h, namespqn_h, typespqn_h,  &
       &  spqn_h, spenergy_h )
ELSE
  numspstates_h= 1
  numspqn_h= 1
  ALLOCATE(namespqn_h(numspqn_h))
  namespqn_h(1)= "fictiousHqn"
  ALLOCATE(typespqn_h(numspqn_h))
  typespqn_h(1)= "fictiousHty"
  ALLOCATE(spqn_h(numspstates_h,numspqn_h))
  spqn_h(1,1)= 0
  ALLOCATE(spenergy_h(numspstates_h))
  spenergy_h(1)= 0.
  CALL LOGGA(2, "fictious single part HOLE states considered:", numspstates_h)
END IF

!..........................................reads possible constrains on H
IF (filein_hconstrains_e /= "" .AND. num_e>0) THEN
  CALL INDATA_HCONSTRAINS("e", numspqn_e, namespqn_e, typespqn_e,   &
       &  numhcons_e, hcons_e)
  CALL LOGGA(2, "num of constraints on H for ELECS", numhcons_e)
ELSE
  numhcons_e= 1
  ALLOCATE(hcons_e(numhcons_e,numspqn_e+2))
  hcons_e(:,1:numspqn_e)= 9999
  hcons_e(:,numspqn_e+1)= nummpenergies
  hcons_e(:,numspqn_e+2)= nummpstates
  CALL LOGGA(2, "no constrain on H for ELECS")
END IF

IF (filein_hconstrains_h /= "" .AND. num_h>0) THEN
  CALL INDATA_HCONSTRAINS("h", numspqn_h, namespqn_h, typespqn_h,   &
       &  numhcons_h, hcons_h)
  CALL LOGGA(2, "num of constraints on H for HOLES:", numhcons_h)
ELSE
  numhcons_h= 1
  ALLOCATE(hcons_h(numhcons_h,numspqn_h+2))
  hcons_h(:,1:numspqn_h)= 9999
  hcons_h(:,numspqn_h+1)= nummpenergies
  hcons_h(:,numspqn_h+2)= nummpstates
  CALL LOGGA(2, "no constrain on H for HOLES")
END IF

!......................reads Coulomb integrals file and builds their indexes
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
    CALL LOGGA(2, "  complex ee Coulomb integrals read", numci_ee)
  ELSE
    numci_ee= 1
    ALLOCATE( ci_ee_x(numci_ee) )
    ALLOCATE( ciindex_ee(numspstates_e,numspstates_e,numspstates_e,numspstates_e) )
    ciindex_ee(:,:,:,:)= 0
  END IF
  IF (num_h>1) THEN
    CALL INDATA_COULOMB_X( "hh", numci_hh, ci_hh_x )
    CALL CIMAKEINDEX_X( "hh", numci_hh, ci_hh_x, ciindex_hh )
    CALL LOGGA(2, "  complex hh Coulomb integrals read", numci_hh)
  ELSE
    numci_hh= 1
    ALLOCATE( ci_hh_x(numci_hh) )
    ALLOCATE( ciindex_hh(numspstates_h,numspstates_h,numspstates_h,numspstates_h) )
    ciindex_hh(:,:,:,:)= 0
  END IF
  IF (num_e>0 .AND. num_h>0) THEN
    CALL INDATA_COULOMB_X( "eh", numci_eh, ci_eh_x )
    CALL CIMAKEINDEX_X( "eh", numci_eh, ci_eh_x, ciindex_eh )
    CALL LOGGA(2, "  complex eh Coulomb integrals read", numci_eh)
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
    CALL LOGGA(2, "  real ee Coulomb integrals read", numci_ee)
  ELSE
    numci_ee= 1
    ALLOCATE( ci_ee(numci_ee) )
    ALLOCATE( ciindex_ee(numspstates_e,numspstates_e,numspstates_e,numspstates_e) )
    ciindex_ee(:,:,:,:)= 0
  END IF
  IF (num_h>1) THEN
    CALL INDATA_COULOMB( "hh", numci_hh, ci_hh )
    CALL CIMAKEINDEX( "hh", numci_hh, ci_hh, ciindex_hh )
    CALL LOGGA(2, "  real hh Coulomb integrals read", numci_hh)
  ELSE
    numci_hh= 1
    ALLOCATE( ci_hh(numci_hh) )
    ALLOCATE( ciindex_hh(numspstates_h,numspstates_h,numspstates_h,numspstates_h) )
    ciindex_hh(:,:,:,:)= 0
  END IF
  IF (num_e>0 .AND. num_h>0) THEN
    CALL INDATA_COULOMB( "eh", numci_eh, ci_eh )
    CALL CIMAKEINDEX( "eh", numci_eh, ci_eh, ciindex_eh )
    CALL LOGGA(2, "  real eh Coulomb integrals read", numci_eh)
  ELSE
    numci_eh= 1
    ALLOCATE( ci_eh(numci_eh) )
    ALLOCATE( ciindex_eh(numspstates_h,numspstates_e,numspstates_e,numspstates_h) )
    ciindex_eh(:,:,:,:)= 0
  END IF
END IF

CALL LOGGA(2, "...done")

!..............................creates separate e and h Hilbert spaces
dimhspace_e= BINOMIALCO(numspstates_e,num_e)
CALL LOGGA(2, "elecs Hilbert space dimension: ", dimhspace_e)
dimhspace_h= BINOMIALCO(numspstates_h,num_h)
CALL LOGGA(2, "holes Hilbert space dimension: ", dimhspace_h)
dimhspaceglobal= dimhspace_e * dimhspace_h
CALL LOGGA(2, "global eh Hilbert space dimension:", dimhspaceglobal)

ALLOCATE( ket_e(dimhspace_e) )
ALLOCATE( ket_h(dimhspace_h) )
!xx ALLOCATE( ket(dimhspace,2) )  ! now in sub composehspace

CALL LOGGA(2, "creating elecs Hilbert space...")
CALL CREATEHSPACE( num_e, numspstates_e, dimhspace_e, ket_e )
CALL LOGGA(2, "...done")

CALL LOGGA(2, "creating holes Hilbert space...")
CALL CREATEHSPACE( num_h, numspstates_h, dimhspace_h, ket_h )
CALL LOGGA(2, "...done")

!DO nn= 1, dimhspace
!  ket(nn,1)= ket_e( (nn-1)/dimhspace_h + 1 )
!  ket(nn,2)= ket_h( MOD(nn-1,dimhspace_h) + 1 )
!END DO

!....................................................prints H space on screen
!WRITE(*,*) "*** elec Hilbert space ***"
!DO nn= 1, dimhspace_e
!  IF ( POPCNT(ket_e(nn)) /= num_e ) STOP "MAIN: 1 bits /= num_e"
!  WRITE(*,binfmt_e) ket_e(nn)
!END DO
!WRITE(*,*)

!WRITE(*,*) "*** hole Hilbert space ***"
!DO nn= 1, dimhspace_h
!  IF ( POPCNT(ket_h(nn)) /= num_h ) STOP "MAIN: 1 bits /= num_h"
!  WRITE(*,binfmt_h) ket_h(nn)
!END DO
!WRITE(*,*)

!...IO..........................................creates file for Hilbert space
IF (fileoutASC_hspace /= "") THEN
  OPEN(21, FILE=TRIM(fileoutASC_hspace), ACTION="WRITE", FORM="FORMATTED")
  WRITE(21,*) "*** reordered Hilbert space(s) for run "//TRIM(ADJUSTL(runname))//" ***"
  WRITE(21,*) "ELEC, HOLE and global H space dim:", dimhspace_e, dimhspace_h, dimhspaceglobal
  WRITE(21,*) "number of constrains for ELECs and HOLEs :", numhcons_e, numhcons_h
  WRITE(21,*)
  CLOSE(21)
END IF
IF (fileoutBIN_hspace /= "") THEN
  OPEN(22, FILE=TRIM(fileoutBIN_hspace), ACTION="WRITE", FORM="UNFORMATTED")
  WRITE(22) citool_version
  WRITE(22) runname
  WRITE(22) dimhspace_e, dimhspace_h
  WRITE(22) numhcons_e, numhcons_h
  CLOSE(22) 
END IF

!...IO.....................................creates file for multi-particle states
IF (fileoutASC_mpstates /= "") THEN
  OPEN(21, FILE=TRIM(fileoutASC_mpstates), ACTION="WRITE", FORM="FORMATTED")
  WRITE(21,*) "CItool version: "//TRIM(citool_version), "  run name: "//TRIM(ADJUSTL(runname)),  &
     &        "  cutoff:", REAL(cutoff_fileoutASC_mpstates)
  WRITE(21,*) "number of constrains for ELECS and HOLES :", numhcons_e, numhcons_h
  WRITE(21,*) 
  CLOSE(21)
END IF
IF (fileoutBIN_mpstates /= "") THEN
  OPEN(22, FILE=TRIM(fileoutBIN_mpstates), ACTION="WRITE", FORM="UNFORMATTED")
  WRITE(22) citool_version
  WRITE(22) runname
  WRITE(22) cutoff_fileoutBIN_mpstates
  WRITE(22) dimhspace_e, dimhspace_h
  WRITE(22) numhcons_e, numhcons_h
  CLOSE(22)
END IF


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC DOUBLE LOOP OVER CONSTRAINS
DO ncons_e= 1, numhcons_e
DO ncons_h= 1, numhcons_h

IF (filein_hconstrains_e /= "" .OR. filein_hconstrains_h /= "") THEN
  CALL LOGGA(2, " ")
  CALL LOGGA(2, "processing H with ELEC constrain", ncons_e)
  CALL LOGGA(2, "              and HOLE constrain", ncons_h)
END IF

nummpenergiescons= MIN(hcons_e(ncons_e,numspqn_e+1),      &
      &  hcons_h(ncons_h,numspqn_h+1), nummpenergies)
nummpstatescons=   MIN(hcons_e(ncons_e,numspqn_e+2),      &
      &  hcons_h(ncons_h,numspqn_h+2), nummpstates)

!...................creates (and allocates) the global eh H space for a given constrain
CALL LOGGA(2, "creating composite eh Hilbert space...")
CALL COMPOSEHSPACE(numspstates_e, numspqn_e, spqn_e, dimhspace_e, ket_e, hcons_e(ncons_e,1:numspqn_e),  &
      &            numspstates_h, numspqn_h, spqn_h, dimhspace_h, ket_h, hcons_h(ncons_h,1:numspqn_h),  &
      &            ket, dimhspacecons)
CALL LOGGA(2, "...done")

IF (filein_hconstrains_e /= "" .OR. filein_hconstrains_h /= "") THEN
  CALL LOGGA(2, "H space dimension for this constrain =", dimhspacecons)
END IF  

ALLOCATE( blockstart(dimhspacecons+1) )
ALLOCATE( blocknonzero(dimhspacecons) )
ALLOCATE( blocknummpenergies(dimhspacecons) )
ALLOCATE( blocknummpstates(dimhspacecons) )

!.......................................... basis reordering for block Hamiltonian
CALL LOGGA(2, "reordering the basis to diagonal block Hamiltonian...")
CALL BLOCKIZEHAMILTONIAN( dimhspacecons, ket, tinye,                   &
     &  numci_ee, ci_ee, ciindex_ee, numci_hh, ci_hh, ciindex_hh,      &
     &  numci_eh, ci_eh, ciindex_eh, numblock, blockstart, blocknonzero )
CALL LOGGA(2, "...done")
CALL LOGGA(2, "number of blocks found= ", numblock)
IF (numblock > 9999)  CALL LOGGA(2, "WARNING ! too many blocks for output format!")
CALL LOGGA(2, "total number of non-zero elements= ", SUM(blocknonzero))

!...IO.................................................writes reordered Hilbert spaces
IF (fileoutASC_hspace /= "") THEN
  OPEN(21, FILE=TRIM(fileoutASC_hspace), ACTION="WRITE", FORM="FORMATTED", POSITION="APPEND")
  IF (filein_hconstrains_e /= "" .OR. filein_hconstrains_h /= "") THEN
    WRITE(21,*) 
    WRITE(21,*) "  * * * * * * ELEC HOLE constrains:", ncons_e, ncons_h
    WRITE(21,*) "              H space dimension for this constrain =", dimhspacecons
  END IF
  DO nb= 1, numblock
    WRITE(21,*) "  * * *  BLOCK", nb
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
      WRITE(21,binfmt_eh) ket(nn,1), ket(nn,2)
    END DO
  END DO
  CLOSE(21)
END IF

IF (fileoutBIN_hspace /= "") THEN
  OPEN(22, FILE=TRIM(fileoutBIN_hspace), ACTION="WRITE", FORM="UNFORMATTED", POSITION="APPEND")
  WRITE(22) ncons_e, ncons_h
  WRITE(22) dimhspacecons
  WRITE(22) numblock
  WRITE(22) blockstart(1:numblock+1)
  WRITE(22) ket
  CLOSE(22)
END IF

!...IO.....................................writes constrain header in multi-particle states file
IF (fileoutASC_mpstates /= "") THEN
  OPEN(21, FILE=TRIM(fileoutASC_mpstates), ACTION="WRITE", FORM="FORMATTED", POSITION="APPEND")
  WRITE(21,*) "  * * * * * ELEC HOLE constrains:", ncons_e, ncons_h
  WRITE(21,*) " H space dimension for this constrain =", dimhspacecons
  WRITE(21,*) " NUMBER OF BLOCKS:", numblock
  WRITE(21,*) " NUMBER OF energies and states:", nummpenergiescons, nummpstatescons
  WRITE(21,*) 
  CLOSE(21)
END IF
IF (fileoutBIN_mpstates /= "") THEN
  OPEN(22, FILE=TRIM(fileoutBIN_mpstates), ACTION="WRITE", FORM="UNFORMATTED", POSITION="APPEND")
  WRITE(22) ncons_e, ncons_h
  WRITE(22) dimhspacecons
  WRITE(22) numblock
  WRITE(22) nummpenergiescons, nummpstatescons
  CLOSE(22)
END IF

!......................................finds maximum dim and nonzeros of blocks
maxblockdim= 0
maxblocknonzero= 0
DO nb= 1, numblock
  maxblockdim= MAX( maxblockdim, blockstart(nb+1)-blockstart(nb) )
  maxblocknonzero= MAX( maxblocknonzero, blocknonzero(nb) )
END DO
CALL LOGGA(2, "max dim of a block in this constrain= ", maxblockdim)
IF (maxblockdim > 99999999)  CALL LOGGA(2, "WARNING ! max dim large for output format!")
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
ALLOCATE( mpenergies(nummpenergiescons,numblock) )
mpenergies= 0.
IF ( complexrun ) THEN
  ALLOCATE( mpstates_x(dimhspacecons,nummpstatescons) )
  mpstates_x= (0.,0.)
  CALL LOGGA(2, "complex mp states storage (Mb)=~", dimhspacecons*nummpstatescons*16/1e6)
ELSE
  ALLOCATE( mpstates(dimhspacecons,nummpstatescons) )
  mpstates= 0.
  CALL LOGGA(2, "real mp states storage (Mb)=~", dimhspacecons*nummpstatescons*8/1e6)
END IF
CALL LOGGA(2, "...done")


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


CALL LOGGA(2, "       - loop over blocks -")
!..................................................begin loop over Ham blocks
DO nb= 1, numblock  !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  CALL LOGGA(2, "  block n. ", nb)

  !............................creates the Hamiltonian matrix for current block
  blockdim= blockstart(nb+1)-blockstart(nb)
  blockfr= blockstart(nb)
  blockto= blockstart(nb+1)-1
  blocknummpenergies(nb)= MIN(nummpenergiescons,MAX(blockdim-2,1))
  blocknummpstates(nb)= MIN(nummpstatescons,MAX(blockdim-2,1))

  CALL LOGGA(2, "current block H dimension= ", blockdim)
  CALL LOGGA(2, "current number of non-zero elements= ", blocknonzero(nb))
  CALL LOGGA(2, "creating the block H matrix...")

  IF ( complexrun ) THEN
    CALL CREATEHAMILTONIAN_X( blockdim, ket(blockfr:blockto,1:2), tinye,       &
         &  numci_ee, ci_ee_x, ciindex_ee, numci_hh, ci_hh_x, ciindex_hh,      &
         &  numci_eh, ci_eh_x, ciindex_eh, spenergy_e, spenergy_h,             &
         &  blocknonzero(nb), hami, hamj, ham_x )
  ELSE
    CALL CREATEHAMILTONIAN( blockdim, ket(blockfr:blockto,1:2), tinye,       &
         &  numci_ee, ci_ee, ciindex_ee, numci_hh, ci_hh, ciindex_hh,        &
         &  numci_eh, ci_eh, ciindex_eh, spenergy_e, spenergy_h,             &
         &  blocknonzero(nb), hami, hamj, ham )
  END IF

  CALL LOGGA(2, "...done")

  !..................................finds block Hamiltonian evals and evects
  CALL LOGGA(2, "diagonalizing the block Hamiltonian...")

  IF ( complexrun ) THEN
    CALL DIAGONALIZE_X( blockdim, hami, hamj, ham_x, blocknummpenergies(nb),       &
         &  mpenergies(1:blocknummpenergies(nb),nb), blocknummpstates(nb),             &
         &  mpstates_x(blockfr:blockto,1:blocknummpstates(nb)) )
  ELSE
    CALL DIAGONALIZE( blockdim, hami, hamj, ham, blocknummpenergies(nb),         &
         &  mpenergies(1:blocknummpenergies(nb),nb), blocknummpstates(nb),           &
         &  mpstates(blockfr:blockto,1:blocknummpstates(nb)) )
  END IF

  CALL LOGGA(2, "...done")

END DO   !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB end loop on the blocks


!...IO.................................writes ASC and BIN mp energies and states
CALL LOGGA(2, "writing multiparticle states...")
!  print*, "writing mp states of BLOCK", nb
!  print*, "blockfr, blockto=", blockfr, blockto
IF ( complexrun ) THEN
  CALL WRITEMPSTATES_X( numblock, dimhspacecons, blockstart, ket,        &
       &  nummpenergiescons, mpenergies, blocknummpenergies,             &
       &  nummpstatescons, mpstates_x, blocknummpstates )
ELSE
  CALL WRITEMPSTATES( numblock, dimhspacecons, blockstart, ket,         &
       &  nummpenergiescons, mpenergies, blocknummpenergies,            &
       &  nummpstatescons, mpstates, blocknummpstates )
END IF
  
CALL LOGGA(2, "...done")

!......................................deallocates arrays of constrains loops
DEALLOCATE( ket )
DEALLOCATE( blockstart )
DEALLOCATE( blocknonzero )
DEALLOCATE( blocknummpenergies )
DEALLOCATE( blocknummpstates )

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


END DO
END DO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  end double loop on constrains

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
DEALLOCATE( typespqn_e )
DEALLOCATE( typespqn_h )
DEALLOCATE( spqn_e )
DEALLOCATE( spqn_h )
DEALLOCATE( spenergy_e )
DEALLOCATE( spenergy_h )
DEALLOCATE( hcons_e )
DEALLOCATE( hcons_h )
DEALLOCATE( ket_e )
DEALLOCATE( ket_h )

CALL LOGGA(2, "...done")

CALL LOGGA(3, " == END ==")

END PROGRAM MAIN
