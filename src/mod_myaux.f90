MODULE mod_myaux
IMPLICIT  NONE
  SAVE

CONTAINS

FUNCTION STRING(inn,npos)
IMPLICIT NONE
!  converts INTEGER "inn" in a text STRING with npos characters
!   by Andrea Bertoni
INTEGER, INTENT(IN) :: inn, npos
CHARACTER(LEN=npos) :: STRING
!.........................................................Tipi locali
INTEGER :: cifra, np, mm, num
IF (npos <= 0)  STOP "ERROR: npos <= 0  in STRING"
IF (inn > (10**npos)-1) STOP "ERROR: (inn > (10**nposin)-1)  in STRING"
num= inn
DO np= 1, npos
  mm= npos-np
  cifra= num/(10**mm)            ! divisione fra interi
  STRING(np:np)= ACHAR(48+cifra)
  num= num - cifra*(10**mm)
END DO
END FUNCTION STRING

END MODULE mod_myaux

