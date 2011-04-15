!$$$$$$$$$$$$$$$$$$$$$$$$$$  density4CItool  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

PROGRAM MAIN 
  USE mod_indatadensity
  USE mod_denscalc
  USE mod_condenscalc
  USE mod_inoutrs
  USE mod_indata  !(from CItool program)
  USE mod_logga
  USE mod_myaux
  USE mod_specialf
  USE mod_indexx
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

CHARACTER(80), PARAMETER :: density4citool_version= "0.9"

! constrains on H (allocated in indata_blockconstrains)
INTEGER :: dimhspace_e, dimhspace_h
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
REAL*8 :: threshold_mpstate

! SP states
INTEGER :: numspqn_e, numspqn_h
CHARACTER(LEN=12), ALLOCATABLE :: namespqn_e(:)
CHARACTER(LEN=12), ALLOCATABLE :: namespqn_h(:)
CHARACTER(LEN=12), ALLOCATABLE :: typespqn_e(:)
CHARACTER(LEN=12), ALLOCATABLE :: typespqn_h(:)
INTEGER, ALLOCATABLE :: spqn_e(:,:)
INTEGER, ALLOCATABLE :: spqn_h(:,:)
REAL*8, ALLOCATABLE :: spenergy_e(:)
REAL*8, ALLOCATABLE :: spenergy_h(:)

! description of densities to be computed
INTEGER :: numdensdesc_e, numdensdesc_h
INTEGER, ALLOCATABLE :: densdesc_e(:,:)
INTEGER, ALLOCATABLE :: densdesc_h(:,:)
CHARACTER(80), ALLOCATABLE :: densfiles_e(:)
CHARACTER(80), ALLOCATABLE :: densfiles_h(:)
LOGICAL, ALLOCATABLE :: masksp(:)

! description of conditional densities to be computed
INTEGER :: numcondensdesc
INTEGER, ALLOCATABLE :: condensdesc(:,:,:)
CHARACTER(80), ALLOCATABLE :: condensfiles(:)
CHARACTER(80), ALLOCATABLE :: condenspoints(:)
REAL*8, ALLOCATABLE ::  psi_fixed(:)
LOGICAL, ALLOCATABLE :: maskspfixed(:)

! psi wave functions
INTEGER :: numx_e, numx_h
INTEGER :: numspwf_e, numspwf_h
REAL*8, ALLOCATABLE ::  psi_e(:,:)       ! sp wave functions for ELECs
REAL*8, ALLOCATABLE ::  psi_h(:,:)       ! sp wave functions for HOLEs

!! INTEGER :: blocknummpenergies, nmpene
REAL*8, ALLOCATABLE :: ene_mpene(:)
INTEGER, ALLOCATABLE :: nblock_mpene(:)
INTEGER, ALLOCATABLE :: nrank_mpene(:)
INTEGER, ALLOCATABLE :: indx_mpene(:)

! aux vars
REAL*8, ALLOCATABLE :: dens(:)
REAL*8, ALLOCATABLE :: densTOT(:)
INTEGER :: wantblockdim, wantblockfr, wantblockto
CHARACTER(80) :: string80
REAL*8 :: normsum, denssum
REAL*8 :: tinye
LOGICAL :: found= .FALSE.
INTEGER :: ncons_e, ncons_h, ncons_e_read, ncons_h_read
INTEGER :: pospsiqn 
INTEGER :: nwant
INTEGER :: ns, nx, nb, nmpene, ne, ndd


!.........................................init vars definitions and readout
CALL INDATADENSITY_GET("density4CItool.nml")

CALL INDATA_GET("citool.nml")

tinye= TINY(1E1)

CALL LOGGA(3, "                 ")
CALL LOGGA(3, "  ===  START density4CItool v. "//TRIM(density4citool_version)//"  ====")
IF (TRIM(density4citoolnml_version) /= TRIM(density4citool_version)) THEN
  CALL LOGGA(3, "WARNING: input density namelist version is "//TRIM(density4citoolnml_version))
END IF
IF (TRIM(citoolnml_version) /= TRIM(density4citool_version)) THEN
  CALL LOGGA(3, "WARNING: citool input namelist version is "//TRIM(citoolnml_version))
END IF

!.....................................reads single-particle states (energies and qn)
CALL LOGGA(2, "number of ELECTRONS: ", num_e)
IF (num_e>0) THEN
  CALL LOGGA(2, "reads single particle ELEC states", numspstates_e)
  CALL INDATA_SPSTATES( "e", numspstates_e, numspqn_e,    &
       &  namespqn_e, typespqn_e, spqn_e, spenergy_e )
!  CALL LOGGA(2, "reading single-particle wave functions for ELECs")
!  ! psi_e is allocated inside the routine
!  CALL INSPWF(numspstates_e, numx_e, psi_e, filein_wavefunction_e)
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
  CALL INDATA_SPSTATES( "h", numspstates_h, numspqn_h,    &
       &  namespqn_h, typespqn_h, spqn_h, spenergy_h )
!  CALL LOGGA(2, "reading single-particle wave functions for HOLEs")
!  ! psi_h is allocated inside the routine
!  CALL INSPWF(numspstates_h, numx_h, psi_h, filein_wavefunction_h)
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


!..................................reads single-particle wave functions
IF (num_e>0) THEN
  CALL LOGGA(2, "reading single-particle wave functions for ELECs")
  ! psi_e is allocated inside the routine
  CALL INSPWF(numspwf_e, numx_e, psi_e, relthreshold_psi_e, fileinBIN_psi_e)
END IF
!IF (relthreshold_psi_e > tinye) THEN
!  threshold_psi_e= relthreshold_psi_e * SUM(ABS(psi_e)**2)/(numspwf_e*numx_e)
!ELSE
!  threshold_psi_e= tinye
!END IF

IF (num_h>0) THEN
  CALL LOGGA(2, "reading single-particle wave functions for HOLEs")
  ! psi_h is allocated inside the routine
  CALL INSPWF(numspwf_h, numx_h, psi_h, relthreshold_psi_h, fileinBIN_psi_h)
END IF
!IF (relthreshold_psi_h > tinye) THEN
!  threshold_psi_h= relthreshold_psi_h * SUM(ABS(psi_h)**2)/(numspwf_h*numx_h)
!ELSE
!  threshold_psi_h= tinye
!END IF

!!$!......................................checks normalization of sp wave function
!!$IF (num_e>0) THEN
!!$  CALL LOGGA(2, "normalization of ELEC sp states")
!!$  DO ns= 1, numspstates_e
!!$    normsum= 0.
!!$    DO nx= 1, numx_e
!!$      normsum= normsum + ABS( psi_e(nx,ns) )**2
!!$    END DO
!!$    CALL LOGGA(2, " norm of state "//STRING(ns,3), normsum)
!!$  END DO
!!$END IF
!!$
!!$IF (num_h>0) THEN
!!$  CALL LOGGA(2, "normalization of HOLE sp states")
!!$  DO ns= 1, numspstates_h
!!$    normsum= 0.
!!$    DO nx= 1, numx_h
!!$      normsum= normsum + ABS( psi_h(nx,ns) )**2
!!$    END DO
!!$    CALL LOGGA(2, " norm of state "//STRING(ns,3), normsum)
!!$  END DO
!!$END IF

!..........................................reads possible constrains on H
IF (filein_hconstrains_e /= "" .AND. num_e>0) THEN
  CALL INDATA_HCONSTRAINS("e", numspqn_e, namespqn_e, typespqn_e,  &
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
  STOP "main: nonvoid filein_hconstrains_h  not implemented"
  CALL INDATA_HCONSTRAINS("h", numspqn_h, namespqn_h, typespqn_h,  &
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

!................reads the description of the densities to be computed
IF (filein_densdescription_e /= "" .AND. num_e>0) THEN
  CALL INDATADENSITY_DENSDESCRIPTION("e", numspqn_e, namespqn_e,typespqn_e,  &
       &  numdensdesc_e, densdesc_e, densfiles_e )
  CALL LOGGA(2, "num of densities to be computed for ELECS", numdensdesc_e)
ELSE
  numdensdesc_e= 0
  CALL LOGGA(2, "no density to be computed for ELECS")
END IF

IF (filein_densdescription_h /= "" .AND. num_h>0) THEN
  CALL INDATADENSITY_DENSDESCRIPTION("h", numspqn_h, namespqn_h,typespqn_h,  &
       &  numdensdesc_h, densdesc_h, densfiles_h )
  CALL LOGGA(2, "num of densities to be computed for HOLES", numdensdesc_h)
ELSE
  numdensdesc_h= 0
  CALL LOGGA(2, "no density to be computed for HOLES")
END IF

!................reads the description of the conditional densities to be computed
IF (filein_condensdescription /= "") THEN
  CALL INDATADENSITY_CONDENSDESCRIPTION(numspqn_e, namespqn_e, numspqn_h, namespqn_h,  &
       &  numcondensdesc, condensdesc, condensfiles, condenspoints )
  CALL LOGGA(2, "num of conditional densities to be computed", numcondensdesc)
ELSE
  numcondensdesc= 0
  CALL LOGGA(2, "no conditional density to be computed")
END IF



print*, "numcondensdesc, numspqn_e=", numcondensdesc, numspqn_e
do ndd= 1, numcondensdesc
  do nx= -1, numspqn_e
    print*, condensdesc(ndd,nx,1), condensdesc(ndd,nx,2)
  end do
end do
!stop



!...................................reads the reordered total Hilbert space
CALL LOGGA(3, "reading Hilbert space from "//TRIM(fileoutBIN_hspace))
OPEN(22, FILE=TRIM(fileoutBIN_hspace), ACTION="READ", FORM="UNFORMATTED")
READ(22) string80
IF (TRIM(string80) /= TRIM(density4citool_version)) STOP "wrong citool_version in bin hspace file"
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
    IF (ncons_e==want_cons .AND. ncons_h==1) THEN
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
IF (.NOT. found) STOP "MAIN: want_cons not found for H"
CALL LOGGA(2, "Hilbert space dimension:", dimhspacecons)


!..........................................reads the multi-carrier states
CALL LOGGA(2, "reading multi-carrier states from "//TRIM(fileoutBIN_mpstates))
OPEN(22, FILE=TRIM(fileoutBIN_mpstates), ACTION="READ", FORM="UNFORMATTED")
READ(22) string80
IF (TRIM(string80) /= TRIM(density4citool_version)) STOP "wrong citool_version in bin mpstates file"
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
    IF (ncons_e==want_cons .AND. ncons_h==1) THEN
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
IF (.NOT. found) STOP "MAIN: want_cons not found for mp"

print*, "dimhspacecons, numblock = ", dimhspacecons, numblock
print*, "nummpenergiescons, nummpstatescons", nummpenergiescons, nummpstatescons

!..........................finds the mp states order with increasing energy
CALL LOGGA(2, "finding the mp states order with increasing energy")
ALLOCATE(ene_mpene(numblock*nummpenergiescons))
ALLOCATE(nblock_mpene(numblock*nummpenergiescons))
ALLOCATE(nrank_mpene(numblock*nummpenergiescons))
ALLOCATE(indx_mpene(numblock*nummpenergiescons))

nmpene= 0
DO nb= 1, numblock
  print*, "num block : ", nb
  print*, "num block energies : ", blocknummpenergies(nb)
  DO ne= 1, blocknummpenergies(nb)
    nmpene= nmpene+1
    ene_mpene(nmpene)= mpenergies(ne,nb)
    nblock_mpene(nmpene)= nb
    nrank_mpene(nmpene)= ne
  END DO
END DO

print*, "nmpene", nmpene
print*, "ene_mpene", ene_mpene

CALL indexx(nmpene, ene_mpene, indx_mpene)

print*, "after call indexx"
print*, "indx_mpene", indx_mpene

DO ne= 1, nmpene
  CALL LOGGA(2, " BLOCK: "//STRING(nblock_mpene(indx_mpene(ne)),4)//     &
     &          "  RANK: "//STRING(nrank_mpene(indx_mpene(ne)),4)//        &
     &          " energy:", REAL(ene_mpene(indx_mpene(ne))))
  !PRINT*, "n , ENERGY, BLOCK, RANK:", ne, REAL(ene_mpene(indx_mpene(ne))),   &
  !     &  nblock_mpene(indx_mpene(ne)), nrank_mpene(indx_mpene(ne))
END DO

!..........................................decides which mp state(s) to consider
IF ( want_energylevel > nmpene ) STOP "MAIN: want_energylevel > nmpene"
IF ( want_energylevel /= 0 ) THEN
  want_block(1)= nblock_mpene(indx_mpene(want_energylevel))
  want_rank(1)= nrank_mpene(indx_mpene(want_energylevel))
  numwantmpstates= 1
END IF

CALL LOGGA(2, "number of MP states included=", numwantmpstates)

DO nwant= 1, numwantmpstates
  print*, "zz", nwant, want_block(nwant), want_rank(nwant), blocknummpenergies(want_block(nwant))
  IF (want_block(nwant)>numblock .OR.  want_rank(nwant)>blocknummpenergies(want_block(nwant)))  &
       &  STOP "MAIN: want_block/rank too large"
END DO

!:::::::::::::::::::::::::::::::::::::::::::::: ELEC density calculation  :::::::
ALLOCATE(dens(numx_e))
ALLOCATE(densTOT(numx_e))
ALLOCATE(masksp(numspstates_e))

!.........................................finding psi-type SP quantum number
pospsiqn= 0
DO nx= 1, numspqn_e
  IF (INDEX(typespqn_e(nx),"psi") /= 0 .OR. INDEX(typespqn_e(nx),"PSI") /= 0) THEN
    pospsiqn= pospsiqn * 9999 + nx
  END IF
END DO
IF (pospsiqn < 1 .OR. pospsiqn > 9999) STOP "no or multiple psi qn type"
CALL LOGGA(2, " ELEC SP psi-type qn is n.", pospsiqn)

!............................................. loop on density descriptions
DO ndd= 1, numdensdesc_e
  densTOT= 0.

  CALL LOGGA(2, "computing ELEC density n. ", ndd)

  !..................................selecting elec SP states to be included
  masksp= .TRUE.
  DO ns= 1, numspstates_e
    DO nx= 1, numspqn_e
      IF (densdesc_e(ndd,nx) /= 9999 .AND. densdesc_e(ndd,nx) /= spqn_e(ns,nx)) THEN
        masksp(ns)= .FALSE.
      END IF
    END DO
  END DO

  DO nwant= 1, numwantmpstates
    dens= 0.

    CALL LOGGA(2, " ELEC density of mp state ", nwant)
    CALL LOGGA(2, " BLOCK", want_block(nwant))
    CALL LOGGA(2, " RANK ", want_rank(nwant))

    wantblockdim= blockstart(want_block(nwant)+1)-blockstart(want_block(nwant))
    wantblockfr= blockstart(want_block(nwant))
    wantblockto= blockstart(want_block(nwant)+1)-1

    IF (relthreshold_mpstate > 2*tinye) THEN
      threshold_mpstate= relthreshold_mpstate *    &
           &  SUM(ABS( mpstates(wantblockfr:wantblockto,want_rank(nwant)) )) / wantblockdim
    ELSE
      threshold_mpstate= 0.
    END IF

    CALL DENSCALC( wantblockdim, mpstates(wantblockfr:wantblockto,want_rank(nwant)),     &
         &  threshold_mpstate,                                                           &
         &  numspstates_e, ket(wantblockfr:wantblockto,1), masksp,                       &
         &  numspwf_e, numx_e, psi_e, spqn_e(1:numspstates_e,pospsiqn),                  &
         &  numspstates_h, ket(wantblockfr:wantblockto,2), dens )
    densTOT= densTOT + dens

  END DO

  !.....................................................writing density to files
  densTOT= densTOT/numwantmpstates
  CALL OUTDENS( numx_e, densTOT, densfiles_e(ndd), denssum )
  CALL LOGGA(2, "integrated density=", denssum)

END DO

!......................................................finalizations
DEALLOCATE(dens)
DEALLOCATE(densTOT)
DEALLOCATE(masksp)


!:::::::::::::::::::::::::::::::::::::::::::::: HOLE density calculation  :::::::

! NOT IMPLEMENTED JET


!::::::::::::::::::::::::::::::::::::::::: conditional density calculation  :::::::
! WARNING: only for ELEC for the moment
DO ndd= 1, numcondensdesc
  IF ( condensdesc(ndd,-1,1) /= 0 ) STOP "MAIN: only ELEC conditional dens for now"
  IF ( condensdesc(ndd,-1,2) /= 0 ) STOP "MAIN: only fixed ELEC for now"
END DO

ALLOCATE(dens(numx_e))
ALLOCATE(densTOT(numx_e))
ALLOCATE(masksp(numspstates_e))
ALLOCATE(psi_fixed(numspwf_e))
ALLOCATE(maskspfixed(numspstates_e))

!.........................................finding psi-type SP quantum number
pospsiqn= 0
DO nx= 1, numspqn_e
  IF (INDEX(typespqn_e(nx),"psi") /= 0 .OR. INDEX(typespqn_e(nx),"PSI") /= 0) THEN
    pospsiqn= pospsiqn * 9999 + nx
  END IF
END DO
IF (pospsiqn < 1 .OR. pospsiqn > 9999) STOP "no or multiple psi qn type"
CALL LOGGA(2, " ELEC SP psi-type qn is n.", pospsiqn)

!........................................ loop on conditional density descriptions
DO ndd= 1, numcondensdesc
  densTOT= 0.

  CALL LOGGA(2, "computing conditional ELEC density n. ", ndd)

  !..................................selecting elec SP states to be included
  masksp= .TRUE.
  DO ns= 1, numspstates_e
    DO nx= 1, numspqn_e  ! condensdesc(ndd,0,1)
      IF (condensdesc(ndd,nx,1) /= 9999 .AND. condensdesc(ndd,nx,1) /= spqn_e(ns,nx)) THEN
        masksp(ns)= .FALSE.
      END IF
    END DO
  END DO
  print*, "masksp", masksp

  !..................................selecting elec SP states for fixed electron
  maskspfixed= .TRUE.
  DO ns= 1, numspstates_e
    DO nx= 1, numspqn_e  ! condensdesc(ndd,0,2)
      IF (condensdesc(ndd,nx,2) /= 9999 .AND. condensdesc(ndd,nx,2) /= spqn_e(ns,nx)) THEN
        maskspfixed(ns)= .FALSE.
      END IF
    END DO
  END DO
  print*, "maskspfixed", maskspfixed

  CALL FINDPSIFIXED( condenspoints(ndd), numx_e, numspwf_e, psi_e, psi_fixed )
  !print*, "psi_fixed", psi_fixed

  DO nwant= 1, numwantmpstates
    dens= 0.

    CALL LOGGA(2, " conditional ELEC density of mp state ", nwant)
    CALL LOGGA(2, " BLOCK", want_block(nwant))
    CALL LOGGA(2, " RANK ", want_rank(nwant))

    wantblockdim= blockstart(want_block(nwant)+1)-blockstart(want_block(nwant))
    wantblockfr= blockstart(want_block(nwant))
    wantblockto= blockstart(want_block(nwant)+1)-1

    IF (relthreshold_mpstate > 2*tinye) THEN
      threshold_mpstate= relthreshold_mpstate *    &
           &  SUM(ABS( mpstates(wantblockfr:wantblockto,want_rank(nwant)) )) / wantblockdim
    ELSE
      threshold_mpstate= 0.
    END IF

    print*, "wantblockdim, wantblockfr, wantblockto", wantblockdim, wantblockfr, wantblockto

    CALL CONDENSCALC( wantblockdim, mpstates(wantblockfr:wantblockto,want_rank(nwant)),  &
         &  threshold_mpstate,                                                           &
         &  numspstates_e, ket(wantblockfr:wantblockto,1),                               &
         &  masksp, numspwf_e, numx_e, psi_e, spqn_e(1:numspstates_e,pospsiqn),          &
         &  numspstates_h, ket(wantblockfr:wantblockto,2),                               &
         &  numspstates_e, ket(wantblockfr:wantblockto,1),                               &
         &  maskspfixed, numspwf_e, psi_fixed, spqn_e(1:numspstates_e,pospsiqn),         &
         &  numspstates_h, ket(wantblockfr:wantblockto,2),                        dens )

    densTOT= densTOT + dens

  END DO

  !.....................................................writing density to files
  densTOT= densTOT/numwantmpstates
  CALL OUTDENS( numx_e, densTOT, condensfiles(ndd), denssum )
  CALL LOGGA(2, "integrated conditional density=", denssum)

  print*

END DO

!......................................................finalizations
DEALLOCATE(dens)
DEALLOCATE(densTOT)
DEALLOCATE(masksp)
DEALLOCATE(psi_fixed)
DEALLOCATE(maskspfixed)

CALL LOGGA(3, "  == END density4CItool ==")

END PROGRAM MAIN
