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

! indatadensity_compute
REAL*8 :: relthreshold_mpstate
REAL*8 :: relthreshold_psi_e
REAL*8 :: relthreshold_psi_h

! indatadensity_inoutput
CHARACTER(80) :: density4citoolnml_version = ""
CHARACTER(80) :: fileinBIN_psi_e= "psi_e.bin"
CHARACTER(80) :: fileinBIN_psi_h= ""
CHARACTER(80) :: filein_densdescription_e= "densdescription_e.dat"
CHARACTER(80) :: filein_densdescription_h= ""
CHARACTER(80) :: filein_condensdescription= ""

NAMELIST /indatadensity_wantmpstate/         &
     &  want_cons,                           &
     &  want_energylevel,                    &
     &  want_block,                          &
     &  want_rank

NAMELIST /indatadensity_compute/             &
     &  relthreshold_mpstate,                &
     &  relthreshold_psi_e,                  &
     &  relthreshold_psi_h

NAMELIST /indatadensity_inoutput/            &
     &  density4citoolnml_version,           &
     &  fileinBIN_psi_e,                     &
     &  fileinBIN_psi_h,                     &
     &  filein_densdescription_e,            &
     &  filein_densdescription_h,            &
     &  filein_condensdescription

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

  READ(33,NML=indatadensity_compute)

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

!===================================================================
SUBROUTINE INDATADENSITY_CONDENSDESCRIPTION( numspqn_e, namespqn_e, numspqn_h, namespqn_h, &
     &  numcondensdesc, condensdesc, condensfiles, condenspoints )
  IMPLICIT NONE
! reads the file with the description of the conditional densities to be computed
! The 3D array condensdesc contains the descriptions, in particular
!   the first index labels the # of conditional density
!   the second index indicates:
!     (#,-1,1) = 0 an E density is requested   r=E 
!     (#,-1,1) = 1 an H density is requested   r=H
!     (#,-1,2) = 0 an E is fixed               f=E
!     (#,-1,2) = 1 an H is fixed               f=H
!     (#,0,1) = numspqn_r   numspqn of requested-dens particles
!     (#,0,2) = numspqn_f   numspqn of fixed particle
!     (#,1:numspqn_r,1) description of the density computed, in terms of spqn
!     (#,1:numspqn_f,2) description of the fixed particle, in terms of spqn
! the element condensfiles(#) contains the filename
! the element condenspoints(#) contains the string with the point of fixed particle

  INTEGER, INTENT(IN) :: numspqn_e, numspqn_h
  CHARACTER(LEN=12), INTENT(IN) :: namespqn_e(numspqn_e), namespqn_h(numspqn_h)
  INTEGER, INTENT(OUT) :: numcondensdesc
  INTEGER, ALLOCATABLE, INTENT(OUT) :: condensdesc(:,:,:)
  CHARACTER(80), ALLOCATABLE, INTENT(OUT) :: condensfiles(:)
  CHARACTER(80), ALLOCATABLE, INTENT(OUT) :: condenspoints(:)

  CHARACTER(LEN=1) :: partype_read
  INTEGER :: numspqn_read
  CHARACTER(LEN=6) :: string6
  CHARACTER(LEN=12) :: string12

  INTEGER :: nd_read, nd, nqn


  OPEN(31, FILE=TRIM(filein_condensdescription),             &
       &   ACTION="READ", STATUS="OLD", FORM="FORMATTED")

  READ(31,*) partype_read, numspqn_read

  IF ( partype_read /= "e" .AND. partype_read /= "E" ) THEN
    STOP "INDATADENSITY_CONDENSDESCRIPTION: 1st partype_read is not an E"
  END IF
  IF ( numspqn_e/=1 .OR. namespqn_e(1)/= "fictiousEqn" ) THEN  ! if num_e/=0
    IF ( numspqn_read /= numspqn_e ) THEN
      STOP "INDATADENSITY_CONDENSDESCRIPTION: numspqn_e does not match"
    END IF
  END IF

!print*, partype_read, numspqn_read

  READ(31,"(A5,X)",ADVANCE="NO") string6

  IF (string6/="NAME:" .AND. string6/="name:") STOP "CONDENSDESCRIPTION: 2nd line error"
  IF ( numspqn_e/=1 .OR. namespqn_e(1)/= "fictiousEqn" ) THEN  ! if num_e/=0
    DO nqn= 1, numspqn_e
      READ(31,"(A12)",ADVANCE="NO") string12
      IF ( string12 /= namespqn_e(nqn) ) THEN
        STOP "CONDENSDESCRIPTION: namespqn_e mismatch"
      END IF
    END DO
  END IF
  READ(31,*)

!print*, string12

  READ(31,*) partype_read, numspqn_read

  IF ( partype_read /= "h" .AND. partype_read /= "H" ) THEN
    STOP "INDATADENSITY_CONDENSDESCRIPTION: 2nd partype_read is not an H"
  END IF
  IF ( numspqn_h/=1 .OR. namespqn_h(1)/= "fictiousHqn" ) THEN  ! if num_h/=0
    IF ( numspqn_read /= numspqn_h ) THEN
      STOP "INDATADENSITY_CONDENSDESCRIPTION: numspqn_h does not match"
    END IF
  END IF

!print*, partype_read, numspqn_read

  READ(31,"(A5,X)",ADVANCE="NO") string6

  IF (string6/="NAME:" .AND. string6/="name:") STOP "CONDENSDESCRIPTION: 4th line error"
  IF ( numspqn_h/=1 .OR. namespqn_h(1)/= "fictiousHqn" ) THEN  ! if num_h/=0
    DO nqn= 1, numspqn_h
      READ(31,"(A12)",ADVANCE="NO") string12
      IF ( string12 /= namespqn_h(nqn) ) THEN
        STOP "CONDENSDESCRIPTION: namespqn_h mismatch"
      END IF
    END DO
  END IF
  READ(31,*)

!print*, string6

  READ(31,*) numcondensdesc

  IF ( numcondensdesc < 1 ) THEN
    STOP "INDATADENSITY_CONDENSDESCRIPTION: numcondensdesc < 1"
  END IF

  ALLOCATE(condensdesc( numcondensdesc,-1:MAX(numspqn_e,numspqn_h),2 ))
  ALLOCATE(condensfiles(numcondensdesc))
  ALLOCATE(condenspoints(numcondensdesc))
  condensdesc(:,:,:)= 9999
  condensfiles(:)= ""
  condenspoints(:)= ""

  DO nd= 1, numcondensdesc
    READ(31,"(I4,X)",ADVANCE="NO") nd_read
    IF ( nd_read /= nd )  STOP "CONDENSDESCRIPTION: condensdesc wrong #"

    ! reads the description of the density to be computed
    READ(31,"(A1)",ADVANCE="NO") partype_read
    IF ( partype_read == "E" .OR. partype_read /= "e" ) THEN
      condensdesc(nd,-1,1)= 0
      condensdesc(nd,0,1)= numspqn_e
    ELSE IF ( partype_read == "H" .OR. partype_read /= "h" ) THEN
      condensdesc(nd,-1,1)= 1
      condensdesc(nd,0,1)= numspqn_h
    ELSE
      STOP "CONDENSDESCRIPTION: wrong partype in a description"
    END IF
    DO nqn= 1, condensdesc(nd,0,1)
      READ(31,"(A12)",ADVANCE="NO") string12
      IF (TRIM(ADJUSTL(string12)) /= "*") THEN
        READ(string12,*) condensdesc(nd,nqn,1)
      END IF
    END DO
    READ(31,*) condensfiles(nd)
    IF (TRIM(ADJUSTL(condensfiles(nd))) == "") STOP "CONDENSDESCRIPTION: void densfiles(nd)"

    ! reads the description of the fixed particle
    READ(31,"(XXXXX,A1)",ADVANCE="NO") partype_read
    IF ( partype_read == "E" .OR. partype_read /= "e" ) THEN
      condensdesc(nd,-1,2)= 0
      condensdesc(nd,0,2)= numspqn_e
    ELSE IF ( partype_read == "H" .OR. partype_read /= "h" ) THEN
      condensdesc(nd,-1,2)= 1
      condensdesc(nd,0,2)= numspqn_h
    ELSE
      STOP "CONDENSDESCRIPTION: wrong partype in a description"
    END IF
    DO nqn= 1, condensdesc(nd,0,2)
      READ(31,"(A12)",ADVANCE="NO") string12
      IF (TRIM(ADJUSTL(string12)) /= "*") THEN
        READ(string12,*) condensdesc(nd,nqn,2)
      END IF
    END DO
    READ(31,*) condenspoints(nd)
    IF (TRIM(ADJUSTL(condenspoints(nd))) == "") STOP "CONDENSDESCRIPTION: void densfiles(nd)"

  END DO

END SUBROUTINE INDATADENSITY_CONDENSDESCRIPTION


END MODULE mod_indatadensity
