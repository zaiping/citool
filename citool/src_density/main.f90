!$$$$$$$$$$$$$$$$$$$$$$$$$$  density4CItool  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

PROGRAM MAIN 
  USE mod_indata  !(from CItool program)
  USE mod_logga
  USE mod_myaux
  USE mod_specialf
  USE mod_indexx
  USE mod_specialf
  USE mod_inoutrs
  USE mod_denscalc
  IMPLICIT NONE
!*************************************************************************
!* density4CItool  obtains the real-space charge density of the
!*  neutral or charged (multi-)excitons computed by CItool. 
!* It needs the single-particle wave functions whose particular
!* format (and domain dimensionality) must be put in mod_inspstates
!*  
!*  adapted from dens4CI(19jul2008) that I wrote for donrodrigo
!*  11 jun 2010 - 14 jun 2010     ifort 11 compiler      by Andrea Bertoni
!*************************************************************************

CHARACTER(80), PARAMETER :: citool_version= "0.9"
INTEGER :: dimhspace_e, dimhspace_h

! constrains on H (allocated in indata_blockconstrains)
INTEGER :: numhcons_e, numhcons_h
INTEGER :: numhcons_e_read, numhcons_h_read
INTEGER, ALLOCATABLE :: hcons_e(:,:)
INTEGER, ALLOCATABLE :: hcons_h(:,:)
INTEGER :: dimhspacecons
INTEGER :: dimhspacecons_read
INTEGER :: nummpenergiescons, nummpstatescons
INTEGER :: nummpenergiescons_read, nummpstatescons_read
INTEGER :: numblock
INTEGER :: numblock_read
REAL*8 :: cutoff_fileoutBIN_mpstates_read
INTEGER, ALLOCATABLE :: blockstart(:)
INTEGER, ALLOCATABLE :: blockstart_read(:)
INTEGER, ALLOCATABLE :: blocknummpenergies(:)
INTEGER, ALLOCATABLE :: blocknummpenergies_read(:)
INTEGER, ALLOCATABLE :: blocknummpstates(:)
INTEGER, ALLOCATABLE :: blocknummpstates_read(:)
INTEGER*8, ALLOCATABLE :: ket(:,:)        ! Slater dets for both e and h
INTEGER*8, ALLOCATABLE :: ket_read(:,:)
REAL*8, ALLOCATABLE :: mpenergies(:,:)
REAL*8, ALLOCATABLE :: mpenergies_read(:,:)
REAL*8, ALLOCATABLE :: mpstates(:,:)
REAL*8, ALLOCATABLE :: mpstates_read(:,:)

LOGICAL, ALLOCATABLE :: masksp_e(:)
LOGICAL, ALLOCATABLE :: masksp_h(:)
INTEGER :: numspqn_e, numspqn_h
CHARACTER(LEN=12), ALLOCATABLE :: namespqn_e(:)
CHARACTER(LEN=12), ALLOCATABLE :: namespqn_h(:)
INTEGER, ALLOCATABLE :: spqn_e(:,:)
INTEGER, ALLOCATABLE :: spqn_h(:,:)
REAL*8, ALLOCATABLE :: spenergy_e(:)
REAL*8, ALLOCATABLE :: spenergy_h(:)
INTEGER :: numx_e, numx_h
REAL*8, ALLOCATABLE ::  psi_e(:,:)       ! sp wave functions for ELECs
REAL*8, ALLOCATABLE ::  psi_h(:,:)       ! sp wave functions for HOLEs

!! INTEGER :: blocknummpenergies, nmpene
REAL*8, ALLOCATABLE :: ene_mpene(:)
INTEGER, ALLOCATABLE :: nblock_mpene(:)
INTEGER, ALLOCATABLE :: nrank_mpene(:)
INTEGER, ALLOCATABLE :: indx_mpene(:)

! aux vars
REAL*8, ALLOCATABLE :: dens(:)
INTEGER :: wantblockdim, wantblockfr, wantblockto
CHARACTER(80) :: string80
REAL*8 :: normsum, denssum
LOGICAL :: found= .FALSE.
INTEGER :: ncons_e, ncons_h, ncons_e_read, ncons_h_read
INTEGER :: ns, nx, nb, nmpene, ne

!.........................................init vars definitions and readout
CALL INDATA_GET("citool.nml")
CALL LOGGA(3, "                 ")
CALL LOGGA(3, "  ===  START density4CItool for v. "//TRIM(citool_version)//"  ====")
IF (TRIM(citoolnml_version) /= TRIM(citool_version)) THEN
  CALL LOGGA(3, "WARNING: input namelist version is "//TRIM(citoolnml_version))
END IF

!........reads single-particle states (energies and qn) and wave functions
ALLOCATE(masksp_e(numspstates_e))
ALLOCATE(masksp_h(numspstates_h))

CALL LOGGA(2, "number of ELECTRONS: ", num_e)
IF (num_e>0) THEN
  CALL LOGGA(2, "reads single particle ELEC states", numspstates_e)
  CALL INDATA_SPSTATES( "e", numspstates_e, numspqn_e,    &
       &  namespqn_e, spqn_e, spenergy_e )
  CALL LOGGA(2, "reading single-particle wave functions for ELECs")
  ! psi_e is allocated inside the routine
  CALL INSPWF(numspstates_e, numx_e, psi_e, FILEwavefunction_e)
ELSE
  numspstates_e= 1
  numspqn_e= 1
  ALLOCATE(namespqn_e(numspqn_e))
  namespqn_e(1)= "fictiousEqn"
  ALLOCATE(spqn_e(numspstates_e,numspqn_e))
  spqn_e(1,1)= 0
  ALLOCATE(spenergy_e(numspstates_e))
  spenergy_e(1)= 0.
  CALL LOGGA(2, "fictious single part ELEC states considered:", numspstates_e)
END IF

CALL LOGGA(2, "number of HOLES: ", num_h)
IF (num_h>0) THEN
  CALL LOGGA(2, "reads single particle HOLE states", numspstates_h)
  CALL INDATA_SPSTATES( "h", numspstates_h, numspqn_h,    &
       &  namespqn_h, spqn_h, spenergy_h )
  CALL LOGGA(2, "reading single-particle wave functions for HOLEs")
  ! psi_h is allocated inside the routine
  CALL INSPWF(numspstates_h, numx_h, psi_h, FILEwavefunction_h)
ELSE
  numspstates_h= 1
  numspqn_h= 1
  ALLOCATE(namespqn_h(numspqn_h))
  namespqn_h(1)= "fictiousHqn"
  ALLOCATE(spqn_h(numspstates_h,numspqn_h))
  spqn_h(1,1)= 0
  ALLOCATE(spenergy_h(numspstates_h))
  spenergy_h(1)= 0.
  CALL LOGGA(2, "fictious single part HOLE states considered:", numspstates_h)
END IF

!......................................checks normalization of sp states
IF (num_e>0) THEN
  CALL LOGGA(2, "normalization of ELEC sp states")
  DO ns= 1, numspstates_e
    normsum= 0.
    DO nx= 1, numx_e
      normsum= normsum + ABS( psi_e(nx,ns) )**2
    END DO
    CALL LOGGA(2, " norm of state "//STRING(ns,3), normsum)
  END DO
END IF

IF (num_h>0) THEN
  CALL LOGGA(2, "normalization of HOLE sp states")
  DO ns= 1, numspstates_h
    normsum= 0.
    DO nx= 1, numx_h
      normsum= normsum + ABS( psi_h(nx,ns) )**2
    END DO
    CALL LOGGA(2, " norm of state "//STRING(ns,3), normsum)
  END DO
END IF

!..........................................reads possible constrains on H
IF (filein_hconstrains_e /= "" .AND. num_e>0) THEN
  CALL INDATA_HCONSTRAINS("e", numspqn_e, namespqn_e,   &
       &  numhcons_e, hcons_e)
  CALL LOGGA(2, "num of constraints on H for ELECS", numhcons_e)
ELSE
  numhcons_e= 1
  ALLOCATE(hcons_e(numhcons_e,numspqn_e+2))
  hcons_e(:,1:numspqn_e)= 9999
  hcons_e(:,numspqn_e+1)= nummpenergies
  hcons_e(:,numspqn_e+2)= nummpstates
  CALL LOGGA(2, "no constrain on blocks for ELECS")
END IF

IF (filein_hconstrains_h /= "" .AND. num_h>0) THEN
  CALL INDATA_HCONSTRAINS("h", numspqn_h, namespqn_h,   &
       &  numhcons_h, hcons_h)
  CALL LOGGA(2, "num of constraints on H for HOLES:", numhcons_h)
ELSE
  numhcons_h= 1
  ALLOCATE(hcons_h(numhcons_h,numspqn_h+2))
  hcons_h(:,1:numspqn_h)= 9999
  hcons_h(:,numspqn_h+1)= nummpenergies
  hcons_h(:,numspqn_h+2)= nummpstates
  CALL LOGGA(2, "no constrain on blocks for HOLES")
END IF


!...................................reads the reordered total Hilbert space
CALL LOGGA(3, "reading Hilbert space from "//TRIM(fileoutBIN_hspace))
OPEN(22, FILE=TRIM(fileoutBIN_hspace), ACTION="READ", FORM="UNFORMATTED")
READ(22) string80
IF (TRIM(string80) /= TRIM(citool_version)) STOP "wrong citool_version in bin hspace file"
READ(22) string80
IF (TRIM(string80) /= TRIM(runname)) STOP "wrong runname in bin hspace file"
READ(22) dimhspace_e, dimhspace_h
IF (dimhspace_e /= BINOMIALCO(numspstates_e,num_e)) STOP "main: wrong dimhspace_e"
IF (dimhspace_h /= BINOMIALCO(numspstates_h,num_h)) STOP "main: wrong dimhspace_h"
READ(22) numhcons_e_read, numhcons_h_read
IF (numhcons_e_read /= numhcons_e .OR. numhcons_h_read /= numhcons_h) &
     &   STOP "main: wrong numhcons in fileoutBIN_hspace"

! DOUBLE LOOP OVER CONSTRAINS
found= .FALSE.
DO ncons_e= 1, numhcons_e
  DO ncons_h= 1, numhcons_h
    READ(22) ncons_e_read, ncons_h_read
    IF (ncons_e_read/=ncons_e .OR. ncons_h_read/=ncons_h) STOP "MAIN: ncons_read/=ncons"
    READ(22) dimhspacecons_read
    READ(22) numblock_read
    ALLOCATE( blockstart_read(dimhspacecons_read+1) )
    ALLOCATE( ket_read(dimhspacecons_read,2) )
    READ(22) blockstart_read(1:numblock_read+1)
    READ(22) ket_read
    IF (ncons_e==WANTCONSE .AND. ncons_h==WANTCONSH) THEN
      found= .TRUE.
      dimhspacecons= dimhspacecons_read
      numblock= numblock_read
      ALLOCATE( blockstart(dimhspacecons+1) )
      blockstart= blockstart_read
      ALLOCATE( ket(dimhspacecons,2) )
      ket= ket_read
    END IF
    DEALLOCATE( blockstart_read )
    DEALLOCATE( ket_read )
  END DO
END DO
CLOSE(22)
IF (.NOT. found) STOP "MAIN: WANTCONSE/WANTCONSH not found for H"
CALL LOGGA(2, "Hilbert space dimension:", dimhspacecons)


!..........................................reads the multi-carrier states
CALL LOGGA(2, "reading multi-carrier states from "//TRIM(fileoutBIN_mpstates))
OPEN(22, FILE=TRIM(fileoutBIN_mpstates), ACTION="READ", FORM="UNFORMATTED")
READ(22) string80
IF (TRIM(string80) /= TRIM(citool_version)) STOP "wrong citool_version in bin mpstates file"
READ(22) string80
IF (TRIM(string80) /= TRIM(runname)) STOP "wrong runname in bin mpstates file"
READ(22) cutoff_fileoutBIN_mpstates_read
IF (cutoff_fileoutBIN_mpstates_read/=cutoff_fileoutBIN_mpstates) STOP "main: cutoffBIN"
CALL LOGGA(2, "the cutoff in bin hspace file is", cutoff_fileoutBIN_mpstates_read)
READ(22) dimhspace_e, dimhspace_h
IF (dimhspace_e /= BINOMIALCO(numspstates_e,num_e)) STOP "main: wrong dimhspace_e"
IF (dimhspace_h /= BINOMIALCO(numspstates_h,num_h)) STOP "main: wrong dimhspace_h"
READ(22) numhcons_e_read, numhcons_h_read
IF (numhcons_e_read/=numhcons_e .OR. numhcons_h_read/=numhcons_h)  &
    &  STOP "MAIN: numhcons_read/=numhcons in mp file"

! DOUBLE LOOP OVER CONSTRAINS
found= .FALSE.
DO ncons_e= 1, numhcons_e
  DO ncons_h= 1, numhcons_h
    READ(22) ncons_e_read, ncons_h_read
    IF (ncons_e_read/=ncons_e .OR. ncons_h_read/=ncons_h) &
      &   STOP "MAIN: ncons_read/=ncons in mp file"
    READ(22) dimhspacecons_read
    READ(22) numblock_read
    READ(22) nummpenergiescons_read, nummpstatescons_read
    ALLOCATE( mpenergies_read(nummpenergiescons_read,numblock_read) )
    ALLOCATE( blocknummpenergies_read(dimhspacecons_read) )
    ALLOCATE( mpstates_read(dimhspacecons_read,nummpstatescons_read) )
    ALLOCATE( blocknummpstates_read(dimhspacecons_read) )
    READ(22) mpenergies_read
    READ(22) blocknummpenergies_read
    READ(22) mpstates_read
    READ(22) blocknummpstates_read
    IF (ncons_e==WANTCONSE .AND. ncons_h==WANTCONSH) THEN
      found= .TRUE.
      IF (dimhspacecons/=dimhspacecons_read) STOP "MAIN: dimhspacecons in mp file"
      IF (numblock/=numblock_read) STOP "MAIN: numblock in mp file"
      nummpenergiescons= nummpenergiescons_read
      nummpstatescons= nummpstatescons_read
      ALLOCATE( mpenergies(nummpenergiescons,numblock) )
      mpenergies= mpenergies_read
      ALLOCATE( blocknummpenergies(dimhspacecons) )
      blocknummpenergies= blocknummpenergies_read
      ALLOCATE( mpstates(dimhspacecons,nummpstatescons) )
      mpstates= mpstates_read
      ALLOCATE( blocknummpstates(dimhspacecons) )
      blocknummpstates= blocknummpstates_read
    END IF
    DEALLOCATE( mpenergies_read )
    DEALLOCATE( blocknummpenergies_read )
    DEALLOCATE( mpstates_read )
    DEALLOCATE( blocknummpstates_read )
  END DO
END DO
IF (.NOT. found) STOP "MAIN: WANTCONSE/WANTCONSH not found for mp"


!..........................finds the mp states order with increasing energy
CALL LOGGA(2, "finding the mp states order with increasing energy")
ALLOCATE(ene_mpene(numblock*nummpenergiescons))
ALLOCATE(nblock_mpene(numblock*nummpenergiescons))
ALLOCATE(nrank_mpene(numblock*nummpenergiescons))
ALLOCATE(indx_mpene(numblock*nummpenergiescons))

nmpene= 0
DO nb= 1, numblock
  DO ne= 1, blocknummpenergies(nb)
    nmpene= nmpene+1
    ene_mpene(nmpene)= mpenergies(ne,nb)
    nblock_mpene(nmpene)= nb
    nrank_mpene(nmpene)= ne
  END DO
END DO

CALL indexx(nmpene, ene_mpene, indx_mpene)

DO ne= 1, nmpene
  CALL LOGGA(2, " BLOCK: "//STRING(nblock_mpene(indx_mpene(ne)),4)//     &
     &          "  RANK: "//STRING(nrank_mpene(indx_mpene(ne)),4)//        &
     &          " energy:", REAL(ene_mpene(indx_mpene(ne))))
  !PRINT*, "n , ENERGY, BLOCK, RANK:", ne, REAL(ene_mpene(indx_mpene(ne))),   &
  !     &  nblock_mpene(indx_mpene(ne)), nrank_mpene(indx_mpene(ne))
END DO

!..........................................decides which mp state to consider
IF ( WANTBLOCK == 0 ) THEN
  WANTBLOCK= nblock_mpene(indx_mpene(1))
  WANTRANK= nrank_mpene(indx_mpene(1))
END IF
wantblockdim= blockstart(WANTBLOCK+1)-blockstart(WANTBLOCK)
wantblockfr= blockstart(WANTBLOCK)
wantblockto= blockstart(WANTBLOCK+1)-1


!:::::::::::::::::::::::::::::::::::::::::: ELEC density calculation  :::::::

CALL LOGGA(2, "ELECTRON density of mp state BLOCK", WANTBLOCK)
CALL LOGGA(2, "ELECTRON density of mp state RANK ", WANTRANK)

ALLOCATE(dens(numx_e))

!......................................calculating total ELECTRON density
IF (FILEdensTOTe /= "") THEN
  masksp_e= .TRUE.

  CALL LOGGA(2, "calculating Total ELECTRON density")
  CALL DENSCALC( wantblockdim, mpstates(wantblockfr:wantblockto,WANTRANK),             &
       &  numspstates_e, ket(wantblockfr:wantblockto,1), numx_e, psi_e, masksp_e,  &
       &  numspstates_h, ket(wantblockfr:wantblockto,2), dens )
  CALL OUTDENS( numx_e, dens, FILEdensTOTe, denssum )
  CALL LOGGA(2, "integrated Total density=", denssum)
END IF

!...................................finding SPIN quantum number
IF (FILEdensUPe /= "" .OR. FILEdensDNe /= "") THEN
  ne= 0
  DO nx= 1, numspqn_e
    !print*, namespqn_e(nx)
    IF (INDEX(namespqn_e(nx),"spin") /= 0) THEN
      !print*, "zzzz"
      ne= ne*9999 + nx
    END IF
  END DO
  IF (ne < 1 .OR. ne > 9999) STOP "no or multiple spin sp qn"
  CALL LOGGA(2, " ELEC sp SPIN is QN", ne)
  !print*, namespqn_e(:)
END IF

!...................................calculating SPIN-UP ELECTRON density
IF (FILEdensUPe /= "") THEN
  CALL LOGGA(2, "Computing density of sp ELEC SPIN-UP states:")
  DO ns= 1, numspstates_e
    IF (spqn_e(ns,ne)==1) THEN
      masksp_e(ns)= .TRUE.
      CALL LOGGA(2, "  ", ns)
    ELSE IF (spqn_e(ns,ne)==0) THEN
      masksp_e(ns)= .FALSE.
    ELSE
      STOP "spin QN /= 0 or 1 in sp file"
    END IF
  END DO

  CALL LOGGA(2, "calculating SPIN UP ELECTRON density")
  CALL DENSCALC( wantblockdim, mpstates(wantblockfr:wantblockto,WANTRANK),         &
       &  numspstates_e, ket(wantblockfr:wantblockto,1), numx_e, psi_e, masksp_e,  &
       &  numspstates_h, ket(wantblockfr:wantblockto,2), dens )
  CALL OUTDENS( numx_e, dens, FILEdensUPe, denssum )
  CALL LOGGA(2, "integrated spin Up density=", denssum)
END IF

!...................................calculating SPIN-DN ELECTRON density
IF (FILEdensDNe /= "") THEN
  CALL LOGGA(2, "Computing density of sp ELEC SPIN-DN states:")
  DO ns= 1, numspstates_e
    IF (spqn_e(ns,ne)==0) THEN
      masksp_e(ns)= .TRUE.
      CALL LOGGA(2, "  ", ns)
    ELSE IF (spqn_e(ns,ne)==1) THEN
      masksp_e(ns)= .FALSE.
    ELSE
      STOP "spin QN /= 0 or 1 in sp file"
    END IF
  END DO

  CALL LOGGA(2, "calculating SPIN DN ELECTRON density")
  CALL DENSCALC( wantblockdim, mpstates(wantblockfr:wantblockto,WANTRANK),         &
       &  numspstates_e, ket(wantblockfr:wantblockto,1), numx_e, psi_e, masksp_e,  &
       &  numspstates_h, ket(wantblockfr:wantblockto,2), dens )
  CALL OUTDENS( numx_e, dens, FILEdensDNe, denssum )
  CALL LOGGA(2, "integrated spin Down density=", denssum)
END IF

DEALLOCATE(dens)


!:::::::::::::::::::::::::::::::::::::::::: HOLE density calculation  :::::::

! not implemented



CALL LOGGA(3, "  == END ==")

END PROGRAM MAIN
