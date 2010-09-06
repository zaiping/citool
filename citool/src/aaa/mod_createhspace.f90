!$$$$$$$$$$$$$$$$$$$$$$$$$  myci  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_createhspace
  IMPLICIT  NONE
  SAVE

CONTAINS

SUBROUTINE createhspace( numpart, numspstates, dimhspace, ket )
  IMPLICIT NONE
! creates the array of kets representing the Slater determinants
! in binary notation

INTEGER, INTENT(IN) :: numpart            ! n of particles (elecs or holes)
INTEGER, INTENT(IN) :: numspstates        ! n single-particle states
INTEGER, INTENT(IN) :: dimhspace          ! dim of Hilbert space
INTEGER*8, INTENT(OUT) :: ket(dimhspace)

INTEGER :: index

IF ( numpart < 0 ) THEN
  STOP "ERROR mod_hilbertspace: numpart < 0"
ELSE IF ( numpart > numspstates ) THEN
  STOP "ERROR mod_hilbertspace: numpart > numspstates"
ELSE IF ( dimhspace < 1 ) THEN
  STOP "ERROR mod_hilbertspace: dimhspace < 1"
END IF

IF ( dimhspace == 1 ) THEN
  ket(1)= 2**numpart - 1
  RETURN
END IF

index= 1
ket(1)= 2**numpart - 1

CALL movebits(1, numspstates-numpart, ket, index)

IF (dimhspace /= index) THEN
  PRINT*, "dimhspace, index : ", dimhspace, index
  STOP  "ERROR mod_hilbertspace: wrong dimhspace"
END IF

DO index= 1, dimhspace
  IF ( POPCNT(ket(index)) /= numpart ) THEN
    PRINT*, "ERROR mod_hilbertspace: number of 1 bits /= numpart"
    PRINT*, "                        for index", index
    STOP
  END IF
END DO

RETURN


CONTAINS

  RECURSIVE SUBROUTINE movebits(fr, to, ketx, mm)
    INTEGER, INTENT(IN) :: fr, to
    INTEGER*8, INTENT(INOUT) :: ketx(:)
    INTEGER, INTENT(INOUT) :: mm
    !
    INTEGER*8 :: kettmp1, kettmp2
    INTEGER :: n1
    !
    DO n1= fr, to
      kettmp1= ketx(mm)
      IF ( BTEST(ketx(mm),n1) ) THEN
        CALL movebits(n1+1, to+1, ketx, mm)
      END IF
      mm= mm + 1
      ! PRINT*, BIT_SIZE(kettmp)
      ! STOP
      kettmp2= ISHFT(kettmp1, 1)
      ketx(mm)= ISHFTC(kettmp2, -1, n1)

      ! IF ( POPCNT(ketx(mm)) /= 3 ) PRINT*, "dddddddddddddd"
      ! WRITE(*,*) n1, fr, to
      ! WRITE(*,"(B64.64)") ketx(mm)
      
    END DO
  END SUBROUTINE movebits


END SUBROUTINE createhspace


END MODULE mod_createhspace
