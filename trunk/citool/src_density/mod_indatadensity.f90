MODULE mod_indatadensity
  IMPLICIT  NONE
  SAVE

! indatadensity_wantmpstate
INTEGER :: want_cons
INTEGER :: want_energylevel
INTEGER :: want_block(21)= 0
INTEGER :: want_rank(21)= 0
! aux vars
INTEGER :: numwantmpstates


! indatadensity_inoutput
CHARACTER(80) :: density4citoolnml_version = ""
CHARACTER(80) :: fileinBIN_psi_e= "psi_e.bin"
CHARACTER(80) :: fileinBIN_psi_h= ""
CHARACTER(80) :: filein_densdescription_e= "densdescription_e.dat"
CHARACTER(80) :: filein_densdescription_h= ""

NAMELIST /indatadensity_wantmpstate/         &
     &  want_cons,                           &
     &  want_energylevel,                    &
     &  want_block,                          &
     &  want_rank

NAMELIST /indatadensity_inoutput/            &
     &  density4citoolnml_version,           &
     &  fileinBIN_psi_e,                     &
     &  fileinBIN_psi_h,                     &
     &  filein_densdescription_e,            &
     &  filein_densdescription_h

CONTAINS


!===================================================================
SUBROUTINE INDATADENSITY_GET(nmlfile)
!  USE mod_myaux
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: nmlfile
  INTEGER :: nn

  OPEN(33, FILE=TRIM(nmlfile), FORM="FORMATTED", ACTION="READ",  &
       &   STATUS="OLD")

  READ(33,NML=indatadensity_wantmpstate)
  IF ( want_energylevel < 0 ) STOP "INDATADENSITY_GET:  want_energylevel < 0"
  IF ( want_energylevel == 0 ) THEN
    IF ( want_block(21) /= 0 .OR. want_rank(21) /= 0 )  &
         &  STOP "INDATADENSITY_GET: only 20 want_block/rank values please"
    !IF ( want_block(1) == 0 .NEQV. want_rank(1) == 0 )  &
    !     &  STOP "INDATADENSITY_GET: want_block AND want_rank = 0 for ground state"
    DO nn= 1, 20
      IF ( want_block(nn) /= 0 .NEQV. want_rank(nn) /= 0 )  &
           &  STOP "INDATADENSITY_GET: want_block/rank different n of values "
      IF ( want_block(nn) == 0 .AND. want_rank(nn) == 0 ) EXIT
    END DO
    numwantmpstates= nn-1
    IF ( numwantmpstates < 1) &
       &  STOP "INDATADENSITY_GET: numwantmpstates < 1"
  END IF

  READ(33,NML=indatadensity_inoutput)

  CLOSE(33)

END SUBROUTINE INDATADENSITY_GET

!===================================================================
SUBROUTINE INDATADENSITY_DENSDESCRIPTION( partype, numspqn, namespqn, typespqn,  &
     &  numdensdesc, densdesc, densfiles )
  IMPLICIT NONE
! reads the file with the description of the densities to be computed
  CHARACTER(*), INTENT(IN) :: partype
  INTEGER, INTENT(IN) :: numspqn
  CHARACTER(LEN=12), INTENT(IN) :: namespqn(numspqn)
  CHARACTER(LEN=12), INTENT(IN) :: typespqn(numspqn)
  INTEGER, INTENT(OUT) :: numdensdesc
  INTEGER, ALLOCATABLE, INTENT(OUT) :: densdesc(:,:)
  CHARACTER(80), ALLOCATABLE, INTENT(OUT) :: densfiles(:)

  CHARACTER(80) :: filein_densdescription
  CHARACTER(1) :: partype_arg, partype_read
  INTEGER :: numspqn_read
  CHARACTER(6) :: string6
  CHARACTER(LEN=12) :: namespqn_read, typespqn_read, dede_read

  INTEGER :: nd_read, nd, nqn

  IF ( partype=="e" .OR. partype=="E" ) THEN
    partype_arg= "e"
    filein_densdescription= filein_densdescription_e
  ELSE IF ( partype=="h" .OR. partype=="H" ) THEN
    partype_arg= "h"
    filein_densdescription= filein_densdescription_h
  ELSE
    STOP "INDATADENSITY_DENSDESCRIPTION: unknown partype"
  END IF

  OPEN(31, FILE=TRIM(filein_densdescription),             &
       &   ACTION="READ", STATUS="OLD", FORM="FORMATTED")

  READ(31,*) partype_read

  IF ( IACHAR(partype_read) < 97 ) THEN
    partype_read= ACHAR(IACHAR(partype_read)+32)  ! to lower case
  END IF

  IF ( partype_read /= partype_arg ) THEN
    STOP "INDATADENSITY_DENSDESCRIPTION: partype_read does not match"
  END IF

  READ(31,*) numdensdesc, numspqn_read

  IF ( numdensdesc < 1 ) THEN
    STOP "INDATADENSITY_DENSDESCRIPTION: numdensdesc < 1"
  END IF
  IF ( numspqn_read /= numspqn ) THEN
    STOP "INDATADENSITY_DENSDESCRIPTION: numspqn_read does not match"
  END IF

  ALLOCATE(densdesc(numdensdesc,numspqn))
  ALLOCATE(densfiles(numdensdesc))
  densdesc(:,:)= 9999
  densfiles(:)= ""

  READ(31,"(A6)",ADVANCE="NO") string6
  IF (string6/="NAME: " .AND. string6/="name: ") STOP "DENSDESCRIPTION: 3rd line error"
  DO nqn= 1, numspqn
    READ(31,"(A12)",ADVANCE="NO") namespqn_read
    IF ( namespqn_read /= namespqn(nqn) ) THEN
      STOP "DENSDESCRIPTION: namespqn mismatch"
    END IF
  END DO
  READ(31,*)

  READ(31,"(A6)",ADVANCE="NO") string6
  IF (string6/="TYPE: " .AND. string6/="type: ") STOP "DENSDESCRIPTION: 4th line error"
  DO nqn= 1, numspqn
    READ(31,"(A12)",ADVANCE="NO") typespqn_read
    IF ( typespqn_read /= typespqn(nqn) ) THEN
      STOP "DENSDESCRIPTION: typespqn mismatch"
    END IF
  END DO
  READ(31,*)

  DO nd= 1, numdensdesc
    READ(31,"(I4,XX)",ADVANCE="NO") nd_read
    IF ( nd_read /= nd ) THEN
      STOP "DENSDESCRIPTION: densdesc wrong #"
    END IF
    DO nqn= 1, numspqn
      READ(31,"(A12)",ADVANCE="NO") dede_read
      IF (TRIM(ADJUSTL(dede_read)) /= "*") THEN
        READ(dede_read,*) densdesc(nd,nqn)
      END IF
    END DO
    READ(31,*) densfiles(nd)
    IF (TRIM(ADJUSTL(densfiles(nd))) == "") STOP "DENSDESCRIPTION: void densfiles(nd)"
  END DO

END SUBROUTINE INDATADENSITY_DENSDESCRIPTION


END MODULE mod_indatadensity
