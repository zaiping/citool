!$$$$$$$$$$$$$$$$$$$$$$$$$  CItool  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_cimakeindex
  IMPLICIT  NONE
  SAVE

CONTAINS

SUBROUTINE CIMAKEINDEX_X( citype, numci, ci_x, ciindex )
  USE mod_indata
  IMPLICIT NONE
! creates the 4D index array for Coulomb integrals
  CHARACTER(*), INTENT(IN) :: citype
  INTEGER, INTENT(IN) :: numci
  TYPE( ci_type_complex16 ), INTENT(IN) :: ci_x(numci)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: ciindex(:,:,:,:)

  INTEGER :: numspstates14, numspstates23
  INTEGER :: nc
  REAL*8 :: tinye
  
  IF ( citype=="ee" .OR. citype=="EE" ) THEN
    numspstates14= numspstates_e
    numspstates23= numspstates_e
  ELSE IF ( citype=="hh" .OR. citype=="HH" ) THEN
    numspstates14= numspstates_h
    numspstates23= numspstates_h
  ELSE IF ( citype=="eh" .OR. citype=="EH" ) THEN
    numspstates14= numspstates_h
    numspstates23= numspstates_e
  ELSE
     STOP "CIMAKEINDEX: unknown citype"
  END IF

  IF ( numspstates14 > 256 .OR. numspstates23 > 256 )  &
       &  STOP "CIMAKEINDEX: numspstates too large"

  tinye= TINY(1E1)

  ALLOCATE( ciindex(numspstates14,numspstates23,numspstates23,numspstates14) )
  
  ciindex(:,:,:,:)= 0

  DO nc= 1, numci
    IF ( ABS(ci_x(nc)%v) > tinye ) THEN
      IF (   ciindex(ci_x(nc)%n1,ci_x(nc)%n2,ci_x(nc)%n3,ci_x(nc)%n4)==0 .AND.  &
           & ciindex(ci_x(nc)%n4,ci_x(nc)%n3,ci_x(nc)%n2,ci_x(nc)%n1)==0 ) THEN
        ciindex(ci_x(nc)%n4,ci_x(nc)%n3,ci_x(nc)%n2,ci_x(nc)%n1) = -nc  !2be conjg
        ciindex(ci_x(nc)%n1,ci_x(nc)%n2,ci_x(nc)%n3,ci_x(nc)%n4) = nc 
      ELSE
        PRINT*, ci_x(nc)%n1,ci_x(nc)%n2,ci_x(nc)%n3,ci_x(nc)%n4
        STOP "CIMAKEINDEX: double-defined Coulomb integral (1234=4321*)"
      END IF
    END IF
  END DO

  ! for ee and hh there is an additional symmetry
  IF ( citype/="eh" .AND. citype/="EH" ) THEN
    DO nc= 1, numci
      IF ( ABS(ci_x(nc)%v) > tinye ) THEN
        IF ( ci_x(nc)%n1/=ci_x(nc)%n3 .OR. ci_x(nc)%n2/=ci_x(nc)%n4 ) THEN

        IF (   ciindex(ci_x(nc)%n2,ci_x(nc)%n1,ci_x(nc)%n4,ci_x(nc)%n3)==0 .AND.  &
             & ciindex(ci_x(nc)%n3,ci_x(nc)%n4,ci_x(nc)%n1,ci_x(nc)%n2)==0 ) THEN
          ciindex(ci_x(nc)%n3,ci_x(nc)%n4,ci_x(nc)%n1,ci_x(nc)%n2) = -nc  !2be conjg
          ciindex(ci_x(nc)%n2,ci_x(nc)%n1,ci_x(nc)%n4,ci_x(nc)%n3) = nc
        ELSE
          PRINT*, ci_x(nc)%n1,ci_x(nc)%n2,ci_x(nc)%n3,ci_x(nc)%n4
          STOP "CIMAKEINDEX: double-defined Coulomb integral (2143=3412*)"
        END IF

        END IF
      END IF
    END DO
  END IF

END SUBROUTINE CIMAKEINDEX_X


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE CIMAKEINDEX( citype, numci, ci, ciindex )
  USE mod_indata
  IMPLICIT NONE
! creates the 4D index array for Coulomb integrals
  CHARACTER(*), INTENT(IN) :: citype
  INTEGER, INTENT(IN) :: numci
  TYPE( ci_type_real8 ), INTENT(IN) :: ci(numci)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: ciindex(:,:,:,:)

  INTEGER :: numspstates14, numspstates23
  INTEGER :: nc
  REAL*8 :: tinye
  
  IF ( citype=="ee" .OR. citype=="EE" ) THEN
    numspstates14= numspstates_e
    numspstates23= numspstates_e
  ELSE IF ( citype=="hh" .OR. citype=="HH" ) THEN
    numspstates14= numspstates_h
    numspstates23= numspstates_h
  ELSE IF ( citype=="eh" .OR. citype=="EH" ) THEN
    numspstates14= numspstates_h
    numspstates23= numspstates_e
  ELSE
     STOP "CIMAKEINDEX: unknown citype"
  END IF

  IF ( numspstates14 > 256 .OR. numspstates23 > 256 )  &
       &  STOP "CIMAKEINDEX: numspstates too large"

  tinye= TINY(1E1)

  ALLOCATE( ciindex(numspstates14,numspstates23,numspstates23,numspstates14) )
  
  ciindex(:,:,:,:)= 0

  DO nc= 1, numci
    IF ( ABS(ci(nc)%v) > tinye ) THEN
      IF (   ciindex(ci(nc)%n1,ci(nc)%n2,ci(nc)%n3,ci(nc)%n4)==0 .AND.  &
           & ciindex(ci(nc)%n4,ci(nc)%n3,ci(nc)%n2,ci(nc)%n1)==0 ) THEN
        ciindex(ci(nc)%n4,ci(nc)%n3,ci(nc)%n2,ci(nc)%n1) = nc
        ciindex(ci(nc)%n1,ci(nc)%n2,ci(nc)%n3,ci(nc)%n4) = nc 
      ELSE
        PRINT*, ci(nc)%n1,ci(nc)%n2,ci(nc)%n3,ci(nc)%n4
        STOP "CIMAKEINDEX: double-defined Coulomb integral (1234=4321*)"
      END IF
    END IF
  END DO

  ! for ee and hh there is an additional symmetry
  IF ( citype/="eh" .AND. citype/="EH" ) THEN
    DO nc= 1, numci
      IF ( ABS(ci(nc)%v) > tinye ) THEN
        IF ( ci(nc)%n1/=ci(nc)%n3 .OR. ci(nc)%n2/=ci(nc)%n4 ) THEN

        IF (   ciindex(ci(nc)%n2,ci(nc)%n1,ci(nc)%n4,ci(nc)%n3)==0 .AND.  &
             & ciindex(ci(nc)%n3,ci(nc)%n4,ci(nc)%n1,ci(nc)%n2)==0 ) THEN
          ciindex(ci(nc)%n3,ci(nc)%n4,ci(nc)%n1,ci(nc)%n2) = nc
          ciindex(ci(nc)%n2,ci(nc)%n1,ci(nc)%n4,ci(nc)%n3) = nc
        ELSE
          PRINT*, ci(nc)%n1,ci(nc)%n2,ci(nc)%n3,ci(nc)%n4
          STOP "CIMAKEINDEX: double-defined Coulomb integral (2143=3412*)"
        END IF

        END IF
      END IF
    END DO
  END IF



!!$  IF (citype=="eh" .OR. citype=="EH") THEN
!!$
!!$    DO nc= 1, numci
!!$      IF ( ABS(ci(nc)%v) > tinye ) THEN
!!$        IF (   ciindex(ci(nc)%n1,ci(nc)%n2,ci(nc)%n3,ci(nc)%n4)==0 .AND.  &
!!$             & ciindex(ci(nc)%n4,ci(nc)%n3,ci(nc)%n2,ci(nc)%n1)==0      ) THEN
!!$          ciindex(ci(nc)%n1,ci(nc)%n2,ci(nc)%n3,ci(nc)%n4) = nc 
!!$          ciindex(ci(nc)%n4,ci(nc)%n3,ci(nc)%n2,ci(nc)%n1) = nc
!!$        ELSE
!!$          PRINT*, ci(nc)%n1,ci(nc)%n2,ci(nc)%n3,ci(nc)%n4
!!$          STOP "CIMAKEINDEX: double-defined EH Coulomb integral (1234=4231=1324=4321)"
!!$        END IF
!!$      END IF
!!$    END DO
!!$
!!$  ! for real Ci, in ee and hh cases there are additional symmetries
!!$  ELSE IF (citype=="ee" .OR. citype=="EE" .OR. citype=="hh" .OR. citype=="HH") THEN
!!$
!!$    DO nc= 1, numci
!!$      IF ( ABS(ci(nc)%v) > tinye ) THEN
!!$        IF (   ciindex(ci(nc)%n1,ci(nc)%n2,ci(nc)%n3,ci(nc)%n4)==0 .AND.  &
!!$             & ciindex(ci(nc)%n4,ci(nc)%n3,ci(nc)%n2,ci(nc)%n1)==0 .AND.  &
!!$             & ciindex(ci(nc)%n2,ci(nc)%n1,ci(nc)%n4,ci(nc)%n3)==0 .AND.  &
!!$             & ciindex(ci(nc)%n3,ci(nc)%n4,ci(nc)%n1,ci(nc)%n2)==0     ) THEN
!!$          IF ( ci(nc)%n1/=ci(nc)%n2 .AND. ci(nc)%n3/=ci(nc)%n4 ) THEN
!!$            ciindex(ci(nc)%n1,ci(nc)%n2,ci(nc)%n3,ci(nc)%n4) = nc
!!$            ciindex(ci(nc)%n4,ci(nc)%n3,ci(nc)%n2,ci(nc)%n1) = nc
!!$            ciindex(ci(nc)%n2,ci(nc)%n1,ci(nc)%n4,ci(nc)%n3) = nc
!!$            ciindex(ci(nc)%n3,ci(nc)%n4,ci(nc)%n1,ci(nc)%n2) = nc
!!$          ELSE
!!$            PRINT*, ci(nc)%n1,ci(nc)%n2,ci(nc)%n3,ci(nc)%n4
!!$            STOP "CIMAKEINDEX: EE or HH Coulomb integral has index 1=2 or 3=4"
!!$          END IF
!!$        ELSE
!!$          PRINT*, ci(nc)%n1,ci(nc)%n2,ci(nc)%n3,ci(nc)%n4
!!$          STOP "CIMAKEINDEX: double-defined EE/HH Ci (1234=4231=1324=4321=2143=3142=2413=3412)"
!!$        END IF
!!$      END IF
!!$    END DO
!!$
!!$  END IF

END SUBROUTINE CIMAKEINDEX


END MODULE mod_cimakeindex
