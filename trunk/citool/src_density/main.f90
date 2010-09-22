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

INTEGER :: numblock, dimhspace
INTEGER, ALLOCATABLE :: blockstart(:)
INTEGER*8, ALLOCATABLE :: ket(:,:)       ! Slater dets for both e and h
REAL*8, ALLOCATABLE :: mpenergies(:,:)
REAL*8, ALLOCATABLE :: mpstates(:,:)

CHARACTER(LEN=12), ALLOCATABLE :: namespqn_e(:)
CHARACTER(LEN=12), ALLOCATABLE :: namespqn_h(:)
INTEGER, ALLOCATABLE :: spqn_e(:,:)
INTEGER, ALLOCATABLE :: spqn_h(:,:)
REAL*8, ALLOCATABLE :: spenergy_e(:)
REAL*8, ALLOCATABLE :: spenergy_h(:)
LOGICAL, ALLOCATABLE :: masksp_e(:)
LOGICAL, ALLOCATABLE :: masksp_h(:)
INTEGER :: numx_e, numx_h
REAL*8, ALLOCATABLE ::  psi_e(:,:)       ! sp wave functions for ELECs
REAL*8, ALLOCATABLE ::  psi_h(:,:)       ! sp wave functions for HOLEs

INTEGER :: blocknummpenergies, nmpene
REAL*8, ALLOCATABLE :: ene_mpene(:)
INTEGER, ALLOCATABLE :: nblock_mpene(:)
INTEGER, ALLOCATABLE :: nrank_mpene(:)
INTEGER, ALLOCATABLE :: indx_mpene(:)

REAL*8, ALLOCATABLE :: dens(:)
INTEGER :: wantblockdim, wantblockfr, wantblockto
CHARACTER(80) :: string80
REAL*8 :: normsum, denssum
INTEGER :: ns, nx, nb, ne

!.........................................init vars definitions and readout
CALL INDATA_GET("citool.nml")
CALL LOGGA(3, "                 ")
CALL LOGGA(3, "  ===  START density4CItool ===")

!...................................reads the reordered total Hilbert space
CALL LOGGA(3, "reading Hilbert space from "//TRIM(fileoutBIN_hspace))
OPEN(22, FILE=TRIM(fileoutBIN_hspace), ACTION="READ", FORM="UNFORMATTED")

READ(22) string80
IF (TRIM(string80) /= TRIM(runname)) STOP "wrong runname in bin hspace file"

READ(22) numblock
ALLOCATE( blockstart(numblock+1) )
READ(22) blockstart

READ(22) dimhspace
ALLOCATE( ket(dimhspace,2) )
READ(22) ket

CLOSE(22)

IF (dimhspace /= BINOMIALCO(numspstates_e,num_e)*BINOMIALCO(numspstates_h,num_h)) &
     &   STOP "main: wrong dimhspace"
CALL LOGGA(2, "total Hilbert space dimension:", dimhspace)

!..........................................reads the multi-carrier states
CALL LOGGA(3, "reading multi-carrier states from "//TRIM(fileoutBIN_mpstates))
ALLOCATE( mpenergies(nummpenergies,numblock) )
ALLOCATE( mpstates(dimhspace,nummpstates) )

OPEN(22, FILE=TRIM(fileoutBIN_mpstates), ACTION="READ", FORM="UNFORMATTED")

READ(22) string80
IF (TRIM(string80) /= TRIM(runname)) STOP "wrong runname in bin mpstates file"

READ(22) mpenergies
READ(22) mpstates

CLOSE(22)

!........reads single-particle states (energies and qn) and wave functions
ALLOCATE(namespqn_e(numspqn_e))
ALLOCATE(namespqn_h(numspqn_h))
ALLOCATE(spqn_e(numspstates_e,numspqn_e))
ALLOCATE(spqn_h(numspstates_h,numspqn_h))
ALLOCATE(spenergy_e(numspstates_e))
ALLOCATE(spenergy_h(numspstates_h))
ALLOCATE(masksp_e(numspstates_e))
ALLOCATE(masksp_h(numspstates_h))

!print*, "AAAAAAAAAAAAA"

CALL LOGGA(2, "number of ELECTRONS: ", num_e)
IF (num_e>0) THEN
  CALL LOGGA(2, "reads single particle ELEC states", numspstates_e)
  CALL INDATA_SPSTATES( "e", numspstates_e, numspqn_e,    &
       &  namespqn_e, spqn_e, spenergy_e )
  CALL LOGGA(2, "reading single-particle wave functions for ELECs")
  ! psi_e is allocated inside the routine
  CALL INSPWF(numspstates_e, numx_e, psi_e, FILEwavefunction_e)
END IF

CALL LOGGA(2, "number of HOLES: ", num_h)
IF (num_h>0) THEN
  CALL LOGGA(2, "reads single particle HOLE states", numspstates_h)
  CALL INDATA_SPSTATES( "h", numspstates_h, numspqn_h,    &
       &  namespqn_h, spqn_h, spenergy_h )
  CALL LOGGA(2, "reading single-particle wave functions for HOLEs")
  ! psi_h is allocated inside the routine
  CALL INSPWF(numspstates_h, numx_h, psi_h, FILEwavefunction_h)
END IF

!print*, "BBBBBBBBBBBBBBB"

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

!..........................finds the mp states order with increasing energy
CALL LOGGA(2, "finding the mp states order with increasing energy")
ALLOCATE(ene_mpene(numblock*nummpenergies))
ALLOCATE(nblock_mpene(numblock*nummpenergies))
ALLOCATE(nrank_mpene(numblock*nummpenergies))
ALLOCATE(indx_mpene(numblock*nummpenergies))

nmpene= 0
DO nb= 1, numblock
!  blockdim= blockstart(nb+1)-blockstart(nb)
!  blockfr= blockstart(nb)
!  blockto= blockstart(nb+1)-1
  blocknummpenergies= MIN(nummpenergies, blockstart(nb+1)-blockstart(nb))
!  blocknummpstates= MIN(nummpstates,blockdim)

  DO ne= 1, blocknummpenergies
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
!wantblocknummpstates= MIN(nummpstates,wantblockdim)


!:::::::::::::::::::::::::::::::::::::::::: ELEC density calculation  :::::::

CALL LOGGA(2, "ELECTRON density of mp state BLOCK", WANTBLOCK)
CALL LOGGA(2, "ELECTRON density of mp state RANK ", WANTRANK)

ALLOCATE(dens(numx_e))

!......................................calculating total ELECTRON density
IF (FILEdensTOTe /= "") THEN
  masksp_e= .TRUE.

  CALL LOGGA(2, "calculating Total ELECTRON density")
  CALL DENSCALC( wantblockdim, mpstates(wantblockfr:wantblockto,WANTRANK),         &
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
