PROGRAM GENERATECITEST
  USE mod_indata
  IMPLICIT NONE
!*************************************************************************
!*  created dummy single-perticle and Coulomb integral files to test myci
!*************************************************************************

!............................................................Declarations
! single-particle states
CHARACTER(LEN=12), ALLOCATABLE :: namespqn_e(:)
CHARACTER(LEN=12), ALLOCATABLE :: namespqn_h(:)
INTEGER, ALLOCATABLE :: spqn_e(:,:)
INTEGER, ALLOCATABLE :: spqn_h(:,:)
REAL*8, ALLOCATABLE :: spenergy_e(:)
REAL*8, ALLOCATABLE :: spenergy_h(:)

! Coulomb integrals
INTEGER :: numci_ee, numci_hh, numci_eh
TYPE( ci_type ), ALLOCATABLE :: ci_ee(:)
TYPE( ci_type ), ALLOCATABLE :: ci_hh(:)
TYPE( ci_type ), ALLOCATABLE :: ci_eh(:)
INTEGER, ALLOCATABLE :: ciindex_ee(:,:,:,:)
INTEGER, ALLOCATABLE :: ciindex_hh(:,:,:,:)
INTEGER, ALLOCATABLE :: ciindex_eh(:,:,:,:)

! aux vars
INTEGER :: nqn, ns, nc
INTEGER :: n1, n2, n3, n4
CHARACTER(LEN=80) :: binfmt_e, binfmt_h

!..........................................init vars definitions and readout
CALL INDATA_GET("myci.nml")
IF (numspqn_e/=3) STOP "only works with numspqn_e=3"
IF (numspqn_h/=3) STOP "only works with numspqn_h=3"

ALLOCATE(namespqn_e(numspqn_e))
ALLOCATE(namespqn_h(numspqn_h))
ALLOCATE(spqn_e(numspstates_e,numspqn_e))
ALLOCATE(spqn_h(numspstates_h,numspqn_h))
ALLOCATE(spenergy_e(numspstates_e))
ALLOCATE(spenergy_h(numspstates_h))

!...............fills the array with dummy data
DO n1= 1, numspstates_e
  spqn_e(n1,1)= 1
  spqn_e(n1,2)= n1/10
  spqn_e(n1,3)= n1
  spenergy_e(n1)= 0.001*n1
END DO

DO n1= 1, numspstates_h
  spqn_h(n1,1)= 1
  spqn_h(n1,2)= n1/10
  spqn_h(n1,3)= n1
  spenergy_h(n1)= 0.0001*n1
END DO

numci_ee= 81
ALLOCATE( ci_ee(numci_ee) )
nc= 0
DO n1= 1, 3
  DO n2= 1, 3
    DO n3= 1, 3
      DO n4= 1, 3
        nc= nc + 1
        ci_ee(nc)%n1= n1
        ci_ee(nc)%n2= n2
        ci_ee(nc)%n3= n3
        ci_ee(nc)%n4= n4
        ci_ee(nc)%v= n1*n2*(1e-3,0.0)
      END DO
    END DO
  END DO
END DO

numci_hh= 256
ALLOCATE( ci_hh(numci_hh) )
nc= 0
DO n1= 1, 4
  DO n2= 1, 4
    DO n3= 1, 4
      DO n4= 1, 4
        nc= nc + 1
        ci_hh(nc)%n1= n1
        ci_hh(nc)%n2= n2
        ci_hh(nc)%n3= n3
        ci_hh(nc)%n4= n4
        ci_hh(nc)%v= n1*n2*(1e-3,0.0)
      END DO
    END DO
  END DO
END DO

numci_eh= 144
ALLOCATE( ci_eh(numci_eh) )
nc= 0
DO n1= 1, 4
  DO n2= 1, 3
    DO n3= 1, 3
      DO n4= 1, 4
        nc= nc + 1
        ci_eh(nc)%n1= n1
        ci_eh(nc)%n2= n2
        ci_eh(nc)%n3= n3
        ci_eh(nc)%n4= n4
        ci_eh(nc)%v= n1*n2*(1e-3,0.0)
      END DO
    END DO
  END DO
END DO

!........................................writes ELECS single-particle states
OPEN(21, FILE=TRIM(filein_spstates_e),                    &
       &   ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")

! I'm describing electron s-p states
WRITE(21,*) "E"
! the number of states and the number of quantum numbers
WRITE(21,*) numspstates_e, numspqn_e
! writes the name of quantum numbers
IF (numspqn_e/=3) STOP "only works with numspqn=3"
WRITE(21,"(XXXX,3A12)") "qne1", "qne2", "qne3"
! writes single-particle quantum numbers and energy (in eV !) for each state
DO ns= 1, numspstates_e
  WRITE(21,"(I3,X)",ADVANCE="NO") ns
  DO nqn= 1, numspqn_e
    WRITE(21,"(I12)",ADVANCE="NO") spqn_e(ns,nqn)
  END DO
  WRITE(21,*) spenergy_e(ns)   ! in eV
END DO

CLOSE(21)

!........................................writes HOLES single-particle states
OPEN(21, FILE=TRIM(filein_spstates_h),                    &
       &   ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")

! I'm describing electron s-p states
WRITE(21,*) "H"
! the number of states and the number of quantum numbers
WRITE(21,*) numspstates_h, numspqn_h
! writes the name of quantum numbers
IF (numspqn_h/=3) STOP "only works with numspqn=3"
WRITE(21,"(3A12)") "qnh1", "qnh2", "qnh3"
! writes single-particle quantum numbers and energy (in eV !) for each state
DO ns= 1, numspstates_h
  DO nqn= 1, numspqn_h
    WRITE(21,"(I12)",ADVANCE="NO") spqn_h(ns,nqn)
  END DO
  WRITE(21,*) spenergy_h(ns)   ! in eV
END DO

CLOSE(21)

!..........................................writes Vee Coulomb integral
IF ( fileinformat_coulomb=="dat" ) THEN

  OPEN(21, FILE=TRIM(filein_coulomb_ee),                     &
       &   ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  WRITE(21,*) numci_ee
  DO nc= 1, numci_ee
    WRITE(21,"(1X,2(I3,1X),2X,2(I3,1X),2X)",ADVANCE="NO")   &
         &  ci_ee(nc)%n1, ci_ee(nc)%n2, ci_ee(nc)%n3, ci_ee(nc)%n4
    WRITE(21,*)  ci_ee(nc)%v
  END DO
  CLOSE(21)

ELSE IF ( fileinformat_coulomb=="mci" ) THEN
    
  OPEN(21, FILE=TRIM(filein_coulomb_ee),                     &
       &   ACTION="WRITE", STATUS="UNKNOWN", FORM="UNFORMATTED")
  WRITE(21) numci_ee
  WRITE(21) ci_ee
  CLOSE(21)

ELSE
  STOP "unknown fileinformat_coulomb"
END IF

!..........................................writes Vhh Coulomb integral
IF ( fileinformat_coulomb=="dat" ) THEN

  OPEN(21, FILE=TRIM(filein_coulomb_hh),                     &
       &   ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  WRITE(21,*) numci_hh
  DO nc= 1, numci_hh
    WRITE(21,"(1X,2(I3,1X),2X,2(I3,1X),2X)",ADVANCE="NO")   &
         &  ci_hh(nc)%n1, ci_hh(nc)%n2, ci_hh(nc)%n3, ci_hh(nc)%n4
    WRITE(21,*)  ci_hh(nc)%v
  END DO
  CLOSE(21)

ELSE IF ( fileinformat_coulomb=="mci" ) THEN
    
  OPEN(21, FILE=TRIM(filein_coulomb_hh),                     &
       &   ACTION="WRITE", STATUS="UNKNOWN", FORM="UNFORMATTED")
  WRITE(21) numci_hh
  WRITE(21) ci_hh
  CLOSE(21)

ELSE
  STOP "unknown fileinformat_coulomb"
END IF

!..........................................writes Veh Coulomb integral
IF ( fileinformat_coulomb=="dat" ) THEN

  OPEN(21, FILE=TRIM(filein_coulomb_eh),                     &
       &   ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  WRITE(21,*) numci_eh
  DO nc= 1, numci_eh
    WRITE(21,"(1X,2(I3,1X),2X,2(I3,1X),2X)",ADVANCE="NO")   &
         &  ci_eh(nc)%n1, ci_eh(nc)%n2, ci_eh(nc)%n3, ci_eh(nc)%n4
    WRITE(21,*)  ci_eh(nc)%v
  END DO
  CLOSE(21)

ELSE IF ( fileinformat_coulomb=="mci" ) THEN
    
  OPEN(21, FILE=TRIM(filein_coulomb_eh),                     &
       &   ACTION="WRITE", STATUS="UNKNOWN", FORM="UNFORMATTED")
  WRITE(21) numci_eh
  WRITE(21) ci_eh
  CLOSE(21)

ELSE
  STOP "unknown fileinformat_coulomb"
END IF


END PROGRAM GENERATECITEST

!  the convention for Coulomb Integrals indexes is:
!
!     ci_kl(n1,n2,n3,n4) =  <n1,n2| U |n3,n4>  =  /int dr /int ds 
!     /phi*^l_n1(r) /phi*^k_n2(s) U(|r-s|) /phi^k_n3(s) /phi^l_n4(r)
!
!  where kl = ee ; hh ; eh
!  n1 and n4 refer to l, and n2 and n3 to k
!
! due to U simmetry:  <n1,n2| U |n3,n4>  =  <n2,n1| U |n4,n3>
!...................................................................Types
!  TYPE ci_type
!    INTEGER,          DIMENSION( 4 ) ::  n
!    COMPLEX                          ::  v
!  END TYPE ci_type
!
