!$$$$$$$$$$$$$$$$$$$$$$$$$$$$  CItool  $$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_indata
  USE mod_staticdata
  IMPLICIT  NONE
  SAVE

! indata_singleparticle
  INTEGER :: numspstates_e
  INTEGER :: numspqn_e
  INTEGER :: numspstates_h
  INTEGER :: numspqn_h
! derived
  CHARACTER(120) :: binfmt_e, binfmt_h, binfmt_eh

! indata_multiparticle
  INTEGER :: num_e
  INTEGER :: num_h
  INTEGER :: nummpenergies
  INTEGER :: nummpstates
  LOGICAL :: complexrun

! indata_inoutput
  CHARACTER(80) :: filein_spstates_e 
  CHARACTER(80) :: filein_spstates_h
  CHARACTER(30) :: fileinformat_coulomb  ! Coulomb file (dat|cit real|complex)
  CHARACTER(80) :: filein_coulomb_ee 
  CHARACTER(80) :: filein_coulomb_hh
  CHARACTER(80) :: filein_coulomb_eh
  CHARACTER(80) :: fileout_mpstates
  INTEGER :: loglevel             ! with 0 everything is logged
  CHARACTER(80) :: statusfile     ! name of status file

! indata_parallel  (reserved for future parallelization)
  INTEGER :: npeio_mpi= 0   ! proces for i/o
  INTEGER :: mypenum_mpi= 0

  NAMELIST /indata_singleparticle/             &
       &  numspstates_e,                       &
       &  numspqn_e,                           &
       &  numspstates_h,                       &
       &  numspqn_h

  NAMELIST /indata_multiparticle/              &
       &  num_e,                               &
       &  num_h,                               &
       &  nummpenergies,                       &
       &  nummpstates,                         &
       &  complexrun

  NAMELIST /indata_inoutput/                   &
       &  filein_spstates_e,                   &
       &  filein_spstates_h,                   &
       &  fileinformat_coulomb,                &
       &  filein_coulomb_ee,                   &
       &  filein_coulomb_hh,                   &
       &  filein_coulomb_eh,                   &
       &  fileout_mpstates,                    &
       &  loglevel,                            &
       &  statusfile

CONTAINS


!===================================================================
SUBROUTINE INDATA_GET(nmlfile)
  USE mod_myaux
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: nmlfile
  INTEGER :: nx

  OPEN(33, FILE=TRIM(nmlfile), FORM="FORMATTED", ACTION="READ",  &
       &   STATUS="OLD")

  READ(33,NML=indata_singleparticle)
  IF ( numspstates_e < 0 .OR. numspstates_e > 64 )  &
       &  STOP "INDATA_GET: numspstates_e out of bounds"
  IF ( numspstates_h < 0 .OR. numspstates_h > 64 )  &
       &  STOP "INDATA_GET: numspstates_h out of bounds"
  binfmt_e= "(B"//STRING(numspstates_e,2)//"."//STRING(numspstates_e,2)//")"
  binfmt_h= "(B"//STRING(numspstates_h,2)//"."//STRING(numspstates_h,2)//")"
  binfmt_eh= "(B"//STRING(numspstates_e,2)//"."//STRING(numspstates_e,2)//"' ' "//  &
       &  "B"//STRING(numspstates_h,2)//"."//STRING(numspstates_h,2)//")"

  READ(33,NML=indata_multiparticle)
  IF ( num_e < 0 .OR. num_e > numspstates_e )  &
       &  STOP "INDATA_GET: num_e < 0  or  > numspstates_e"
  IF ( num_h < 0 .OR. num_h > numspstates_h )  &
       &  STOP "INDATA_GET: num_h < 0  or  > numspstates_h"
  IF ( num_e == 0 .AND. num_h == 0  ) &
       &  STOP "INDATA_GET: no particles at all ?? too easy !"
  IF ( nummpenergies < 1 )  &
       &  STOP "INDATA_GET: nummpenergies < 1"
  IF ( nummpstates < 0 .OR. nummpstates > nummpenergies )  &
       &  STOP "INDATA_GET: nummpstates < 0 or > nummpenergies"
  IF ( num_e == 0 ) THEN
    numspstates_e= 1
    numspqn_e= 1
  END IF
  IF ( num_h == 0 ) THEN
    numspstates_h= 1
    numspqn_h= 1
  END IF

  fileinformat_coulomb= REPEAT(" ",30) 
  READ(33,NML=indata_inoutput)
  IF (   INDEX(fileinformat_coulomb,"dat")==0 .EQV.     &
       & INDEX(fileinformat_coulomb,"cit")==0        )  &
       & STOP "INDATA_GET:  Unknown or redundant fileinformat_coulomb 1 !"
  IF (   INDEX(fileinformat_coulomb,"real8")==0 .EQV.        &
       & INDEX(fileinformat_coulomb,"complex16")==0        )  &
       & STOP "INDATA_GET:  Unknown or redundant fileinformat_coulomb 2 !"

  CLOSE(33)

END SUBROUTINE INDATA_GET

!===================================================================
SUBROUTINE INDATA_SPSTATES( partype, numspstates, numspqn,    &
     &  namespqn, spqn, spenergy )
  IMPLICIT NONE
! reads the single-particle states and energies
  CHARACTER(*), INTENT(IN) :: partype
  INTEGER, INTENT(IN) :: numspstates
  INTEGER, INTENT(IN) :: numspqn
  CHARACTER(LEN=12), INTENT(OUT) :: namespqn(numspqn)
  INTEGER, INTENT(OUT) :: spqn(numspstates,numspqn)
  REAL*8, INTENT(OUT) :: spenergy(numspstates)

  CHARACTER(80) :: filein_spstates
  CHARACTER(1) :: partype_arg, partype_read
  INTEGER :: numspstates_read
  INTEGER :: numspqn_read
  INTEGER :: ns_read
  INTEGER :: ns, nqn

  IF ( partype=="e" .OR. partype=="E" ) THEN
    partype_arg= "e"
    filein_spstates= filein_spstates_e
  ELSE IF ( partype=="h" .OR. partype=="H" ) THEN
    partype_arg= "h"
    filein_spstates= filein_spstates_h
  ELSE
    STOP "INDATA_SPSTATES: unknown partype"
  END IF

  OPEN(31, FILE=TRIM(filein_spstates),                    &
       &   ACTION="READ", STATUS="OLD", FORM="FORMATTED")

  READ(31,*) partype_read

  IF ( IACHAR(partype_read) < 97 ) THEN
    partype_read= ACHAR(IACHAR(partype_read)+32)  ! to lower case
  END IF

  IF ( partype_read /= partype_arg ) THEN
    STOP "INDATA_SPSTATES: partype_read does not match"
  END IF

  READ(31,*) numspstates_read, numspqn_read

  IF ( numspstates_read < numspstates ) THEN
    STOP "INDATA_SPSTATES: not enough sp states in the file"
  ELSE IF ( numspstates_read /= numspstates ) THEN
    PRINT*, "WARNING: numspstates < sp states available in the file"
  END IF

  IF ( numspqn_read /= numspqn ) THEN
    STOP "INDATA_SPSTATES: numspqn_read does not match"
  END IF

  READ(31,"(XXXX)",ADVANCE="NO")
  DO nqn= 1, numspqn
    READ(31,"(A12)",ADVANCE="NO") namespqn(nqn)
  END DO
  READ(31,*)

  spenergy(:)= 0d0
  DO ns= 1, numspstates
    READ(31,"(I3,X)",ADVANCE="NO") ns_read
    IF ( spenergy(ns_read) /= 0d0 ) THEN
      STOP "INDATA_SPSTATES: double-defined sp state"
    END IF
    DO nqn= 1, numspqn
      READ(31,"(I12)",ADVANCE="NO") spqn(ns_read,nqn)
    END DO
    READ(31,*) spenergy(ns_read)
  END DO

END SUBROUTINE INDATA_SPSTATES

!===================================================================
SUBROUTINE INDATA_COULOMB_X( citype, numci, ci_x )
  IMPLICIT NONE
! reads COMPLEX Coulomb integrals
!
!  the convention for Coulomb Integrals indexes is:
!
!     ci_kl(n1,n2,n3,n4) =  <n1,n2| U |n3,n4>  =  /int dr /int ds 
!     /phi*^l_n1(r) /phi*^k_n2(s) U(|r-s|) /phi^k_n3(s) /phi^l_n4(r)
!
!  where kl = ee ; hh ; eh
!  n1 and n4 refer to l, and n2 and n3 to k
!
!...................................................................Types
!  TYPE ci_type_complex
!    INTEGER*1, :: n1, n2, n3, n4
!    COMPLEX*16 :: v
!  END TYPE ci_type_complex

  CHARACTER(*), INTENT(IN) :: citype
  INTEGER, INTENT(OUT) :: numci
  TYPE( ci_type_complex16 ), ALLOCATABLE, INTENT(OUT) :: ci_x(:)

  TYPE( ci_type_real8 ), ALLOCATABLE :: ci(:)
  CHARACTER(80) :: filein_coulomb
  REAL*8 :: auxr8
  INTEGER :: numspstates14, numspstates23
  INTEGER :: nc

  IF ( citype=="ee" .OR. citype=="EE" ) THEN
    filein_coulomb= filein_coulomb_ee
    numspstates14= numspstates_e
    numspstates23= numspstates_e
  ELSE IF ( citype=="hh" .OR. citype=="HH" ) THEN
    filein_coulomb= filein_coulomb_hh
    numspstates14= numspstates_h
    numspstates23= numspstates_h
  ELSE IF ( citype=="eh" .OR. citype=="EH" ) THEN
    filein_coulomb= filein_coulomb_eh
    numspstates14= numspstates_h
    numspstates23= numspstates_e
  ELSE
     STOP "INDATA_COULOMB_X: unknown citype"
  END IF

  IF ( INDEX(fileinformat_coulomb,"dat") /= 0 ) THEN

    OPEN(31, FILE=TRIM(filein_coulomb),                     &
         &   ACTION="READ", STATUS="OLD", FORM="FORMATTED")
    READ(31,*) numci
    ALLOCATE( ci_x(numci) )
    IF ( INDEX(fileinformat_coulomb,"complex16") /= 0 ) THEN
      DO nc= 1, numci
        READ(31,"(1X,2(I3,1X),2X,2(I3,1X),2X)",ADVANCE="NO")   &
             &  ci_x(nc)%n1, ci_x(nc)%n2, ci_x(nc)%n3, ci_x(nc)%n4
        READ(31,*)  ci_x(nc)%v
      END DO
    ELSE IF ( INDEX(fileinformat_coulomb,"real8") /= 0 ) THEN
      DO nc= 1, numci
        READ(31,"(1X,2(I3,1X),2X,2(I3,1X),2X)",ADVANCE="NO")   &
             &  ci_x(nc)%n1, ci_x(nc)%n2, ci_x(nc)%n3, ci_x(nc)%n4
        READ(31,*)  auxr8
        ci_x(nc)%v= auxr8
      END DO
    ELSE
      STOP "INDATA_COULOMB: fileinformat_coulomb not complex16 | real8"
    END IF
    CLOSE(31)

  ELSE IF ( fileinformat_coulomb=="cit" ) THEN
    
    OPEN(31, FILE=TRIM(filein_coulomb),                     &
         &   ACTION="READ", STATUS="OLD", FORM="UNFORMATTED")
    READ(31) numci
    ALLOCATE( ci_x(numci) )
    IF ( INDEX(fileinformat_coulomb,"complex16") /= 0 ) THEN
      READ(31) ci_x
    ELSE IF ( INDEX(fileinformat_coulomb,"real8") /= 0 ) THEN
      ALLOCATE( ci(numci) )
      READ(31) ci
      ci_x(:)%n1= ci(:)%n1
      ci_x(:)%n2= ci(:)%n2
      ci_x(:)%n3= ci(:)%n3
      ci_x(:)%n4= ci(:)%n4
      ci_x(:)%v = ci(:)%v
      DEALLOCATE( ci )
    ELSE
      STOP "INDATA_COULOMB: fileinformat_coulomb not complex16 | real8"
    END IF
    CLOSE(31)

  ELSE
    STOP "INDATA_COULOMB: unknown fileinformat_coulomb"
  END IF

  !PRINT*, citype, MAXVAL(ci(:)%n2), numspstates23

  IF (MINVAL(ci_x(:)%n1) < 1 .OR. MAXVAL(ci_x(:)%n1) > numspstates14) &
       & STOP "INDATA_COULOMB: ci n1 index out of bounds"
  IF (MINVAL(ci_x(:)%n2) < 1 .OR. MAXVAL(ci_x(:)%n2) > numspstates23) &
       & STOP "INDATA_COULOMB: ci n2 index out of bounds"
  IF (MINVAL(ci_x(:)%n3) < 1 .OR. MAXVAL(ci_x(:)%n3) > numspstates23) &
       & STOP "INDATA_COULOMB: ci n3 index out of bounds"
  IF (MINVAL(ci_x(:)%n4) < 1 .OR. MAXVAL(ci_x(:)%n4) > numspstates14) &
       & STOP "INDATA_COULOMB: ci n4 index out of bounds"

END SUBROUTINE INDATA_COULOMB_X

!===================================================================
SUBROUTINE INDATA_COULOMB( citype, numci, ci )
  IMPLICIT NONE
! reads REAL Coulomb integrals
!
!  the convention for Coulomb Integrals indexes is:
!
!     ci_kl(n1,n2,n3,n4) =  <n1,n2| U |n3,n4>  =  /int dr /int ds 
!     /phi*^l_n1(r) /phi*^k_n2(s) U(|r-s|) /phi^k_n3(s) /phi^l_n4(r)
!
!  where kl = ee ; hh ; eh
!  n1 and n4 refer to l, and n2 and n3 to k
!
!...................................................................Types
!  TYPE ci_type_real
!    INTEGER*1, :: n1, n2, n3, n4
!    REAL*8 :: v
!  END TYPE ci_type_real

  CHARACTER(*), INTENT(IN) :: citype
  INTEGER, INTENT(OUT) :: numci
  TYPE( ci_type_real8 ), ALLOCATABLE, INTENT(OUT) :: ci(:)

  TYPE( ci_type_complex16 ), ALLOCATABLE :: ci_x(:)
  CHARACTER(80) :: filein_coulomb
  COMPLEX*16 :: auxc16
  INTEGER :: numspstates14, numspstates23
  INTEGER :: nc


  IF ( citype=="ee" .OR. citype=="EE" ) THEN
    filein_coulomb= filein_coulomb_ee
    numspstates14= numspstates_e
    numspstates23= numspstates_e
  ELSE IF ( citype=="hh" .OR. citype=="HH" ) THEN
    filein_coulomb= filein_coulomb_hh
    numspstates14= numspstates_h
    numspstates23= numspstates_h
  ELSE IF ( citype=="eh" .OR. citype=="EH" ) THEN
    filein_coulomb= filein_coulomb_eh
    numspstates14= numspstates_h
    numspstates23= numspstates_e
  ELSE
     STOP "INDATA_COULOMB: unknown citype"
  END IF

  IF ( INDEX(fileinformat_coulomb,"dat") /= 0 ) THEN

    OPEN(31, FILE=TRIM(filein_coulomb),                     &
         &   ACTION="READ", STATUS="OLD", FORM="FORMATTED")
    READ(31,*) numci
    ALLOCATE( ci(numci) )
    IF ( INDEX(fileinformat_coulomb,"real8") /= 0 ) THEN
      DO nc= 1, numci
        READ(31,"(1X,2(I3,1X),2X,2(I3,1X),2X)",ADVANCE="NO")   &
             &  ci(nc)%n1, ci(nc)%n2, ci(nc)%n3, ci(nc)%n4
        READ(31,*)  ci(nc)%v
      END DO
      CLOSE(31)
    ELSE IF ( INDEX(fileinformat_coulomb,"complex16") /= 0 ) THEN
      DO nc= 1, numci
        READ(31,"(1X,2(I3,1X),2X,2(I3,1X),2X)",ADVANCE="NO")   &
             &  ci(nc)%n1, ci(nc)%n2, ci(nc)%n3, ci(nc)%n4
        READ(31,*)  auxc16
        ci(nc)%v= REAL(auxc16)
      END DO
    ELSE
      STOP "INDATA_COULOMB: fileinformat_coulomb not complex16 | real8"
    END IF

  ELSE IF ( fileinformat_coulomb=="cit" ) THEN
    
    OPEN(31, FILE=TRIM(filein_coulomb),                     &
         &   ACTION="READ", STATUS="OLD", FORM="UNFORMATTED")
    READ(31) numci
    ALLOCATE( ci(numci) )
    IF ( INDEX(fileinformat_coulomb,"real8") /= 0 ) THEN
      READ(31) ci
    ELSE IF ( INDEX(fileinformat_coulomb,"complex16") /= 0 ) THEN
      ALLOCATE( ci_x(numci) )
      READ(31) ci_x
      ci(:)%n1= ci_x(:)%n1
      ci(:)%n2= ci_x(:)%n2
      ci(:)%n3= ci_x(:)%n3
      ci(:)%n4= ci_x(:)%n4
      ci(:)%v = ci_x(:)%v
      DEALLOCATE( ci_x )
    ELSE
      STOP "INDATA_COULOMB: fileinformat_coulomb not complex16 | real8"
    END IF
    CLOSE(31)

  ELSE
    STOP "INDATA_COULOMB: unknown fileinformat_coulomb"
  END IF

  !PRINT*, citype, MAXVAL(ci(:)%n2), numspstates23

  IF (MINVAL(ci(:)%n1) < 1 .OR. MAXVAL(ci(:)%n1) > numspstates14) &
       & STOP "INDATA_COULOMB: ci n1 index out of bounds"
  IF (MINVAL(ci(:)%n2) < 1 .OR. MAXVAL(ci(:)%n2) > numspstates23) &
       & STOP "INDATA_COULOMB: ci n2 index out of bounds"
  IF (MINVAL(ci(:)%n3) < 1 .OR. MAXVAL(ci(:)%n3) > numspstates23) &
       & STOP "INDATA_COULOMB: ci n3 index out of bounds"
  IF (MINVAL(ci(:)%n4) < 1 .OR. MAXVAL(ci(:)%n4) > numspstates14) &
       & STOP "INDATA_COULOMB: ci n4 index out of bounds"

END SUBROUTINE INDATA_COULOMB


END MODULE mod_indata
