MODULE mod_staticdata

  IMPLICIT  NONE
  SAVE

! physics constants   (in MKS)
  REAL*8, PARAMETER    :: ELMASS0= 9.1095e-31,                     &
       &                ELCH= 1.602176e-19,                        &
       &                HBAR= 1.05459e-34,                         &
       &                BOLTZ= 1.3807e-23,                         &
       &                PIG=3.14159265,                            &
       &                EPSILON0= 8.854e-12

! def of type for Coulomb Integrals

  TYPE ci_type_complex16
    INTEGER*1  ::  n1, n2, n3, n4
    COMPLEX*16 ::  v
  END TYPE ci_type_complex16

  TYPE ci_type_real8
    INTEGER*1  ::  n1, n2, n3, n4
    REAL*8 ::  v
  END TYPE ci_type_real8

END MODULE mod_staticdata
