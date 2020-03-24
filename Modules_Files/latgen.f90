!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
SUBROUTINE env_latgen(ibrav,celldm,a1,a2,a3,omega)
  !-----------------------------------------------------------------------
  !     sets up the crystallographic vectors a1, a2, and a3.
  !
  !     ibrav is the structure index:
  !       1  cubic P (sc)                8  orthorhombic P
  !       2  cubic F (fcc)               9  1-face (C) centered orthorhombic
  !       3  cubic I (bcc)              10  all face centered orthorhombic
  !       4  hexagonal and trigonal P   11  body centered orthorhombic
  !       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
  !       6  tetragonal P (st)          13  one face (base) centered monoclinic
  !       7  tetragonal I (bct)         14  triclinic P
  !     Also accepted:
  !       0  "free" structure          -12  monoclinic P (unique axis: b)
  !      -3  cubic bcc with a more symmetric choice of axis
  !      -5  trigonal R, threefold axis along (111)
  !      -9  alternate description for base centered orthorhombic
  !     -13  one face (base) centered monoclinic (unique axis: b)
  !      91  1-face (A) centered orthorombic
  !
  !     celldm are parameters which fix the shape of the unit cell
  !     omega is the unit-cell volume
  !
  !     NOTA BENE: all axis sets are right-handed
  !     Boxes for US PPs do not work properly with left-handed axis
  !
  USE env_kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ibrav
  real(DP), INTENT(inout) :: celldm(6)
  real(DP), INTENT(inout) :: a1(3), a2(3), a3(3)
  real(DP), INTENT(out) :: omega
  !
  real(DP), PARAMETER:: sr2 = 1.414213562373d0, &
                        sr3 = 1.732050807569d0
  INTEGER :: i,j,k,l,iperm,ir
  real(DP) :: term, cbya, s, term1, term2, singam, sen
  !
  !  user-supplied lattice vectors
  !
  IF (ibrav == 0) THEN
     IF (sqrt( a1(1)**2 + a1(2)**2 + a1(3)**2 ) == 0 )  &
         CALL env_errore ('latgen', 'wrong at for ibrav=0', 1)
     IF (sqrt( a2(1)**2 + a2(2)**2 + a2(3)**2 ) == 0 )  &
         CALL env_errore ('latgen', 'wrong at for ibrav=0', 2)
     IF (sqrt( a3(1)**2 + a3(2)**2 + a3(3)**2 ) == 0 )  &
         CALL env_errore ('latgen', 'wrong at for ibrav=0', 3)

     IF ( celldm(1) /= 0.D0 ) THEN
     !
     ! ... input at are in units of alat => convert them to a.u.
     !
         a1(:) = a1(:) * celldm(1)
         a2(:) = a2(:) * celldm(1)
         a3(:) = a3(:) * celldm(1)
     ELSE
     !
     ! ... input at are in atomic units: define celldm(1) from a1
     !
         celldm(1) = sqrt( a1(1)**2 + a1(2)**2 + a1(3)**2 )
     ENDIF
     !
  ELSE
     a1(:) = 0.d0
     a2(:) = 0.d0
     a3(:) = 0.d0
  ENDIF
  !
  IF (celldm (1) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(1)', abs(ibrav) )
  !
  !  index of bravais lattice supplied
  !
  IF (ibrav == 1) THEN
     !
     !     simple cubic lattice
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)
     !
  ELSEIF (ibrav == 2) THEN
     !
     !     fcc lattice
     !
     term=celldm(1)/2.d0
     a1(1)=-term
     a1(3)=term
     a2(2)=term
     a2(3)=term
     a3(1)=-term
     a3(2)=term
     !
  ELSEIF (abs(ibrav) == 3) THEN
     !
     !     bcc lattice
     !
     term=celldm(1)/2.d0
     DO ir=1,3
        a1(ir)=term
        a2(ir)=term
        a3(ir)=term
     ENDDO
     IF ( ibrav < 0 ) THEN
        a1(1)=-a1(1)
        a2(2)=-a2(2)
        a3(3)=-a3(3)
     ELSE
        a2(1)=-a2(1)
        a3(1)=-a3(1)
        a3(2)=-a3(2)
     ENDIF
     !
  ELSEIF (ibrav == 4) THEN
     !
     !     hexagonal lattice
     !
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(1)=-celldm(1)/2.d0
     a2(2)=celldm(1)*sr3/2.d0
     a3(3)=celldm(1)*cbya
     !
  ELSEIF (abs(ibrav) == 5) THEN
     !
     !     trigonal lattice
     !
     IF (celldm (4) <= -0.5_dp .or. celldm (4) >= 1.0_dp) &
          CALL env_errore ('latgen', 'wrong celldm(4)', abs(ibrav))
     !
     term1=sqrt(1.0_dp + 2.0_dp*celldm(4))
     term2=sqrt(1.0_dp - celldm(4))
     !
     IF ( ibrav == 5) THEN
        !     threefold axis along c (001)
        a2(2)=sr2*celldm(1)*term2/sr3
        a2(3)=celldm(1)*term1/sr3
        a1(1)=celldm(1)*term2/sr2
        a1(2)=-a1(1)/sr3
        a1(3)= a2(3)
        a3(1)=-a1(1)
        a3(2)= a1(2)
        a3(3)= a2(3)
     ELSEIF ( ibrav == -5) THEN
        !     threefold axis along (111)
        ! Notice that in the cubic limit (alpha=90, celldm(4)=0, term1=term2=1)
        ! does not yield the x,y,z axis, but an equivalent rotated triplet:
        !    a/3 (-1,2,2), a/3 (2,-1,2), a/3 (2,2,-1)
        ! If you prefer the x,y,z axis as cubic limit, you should modify the
        ! definitions of a1(1) and a1(2) as follows:'
        !    a1(1) = celldm(1)*(term1+2.0_dp*term2)/3.0_dp
        !    a1(2) = celldm(1)*(term1-term2)/3.0_dp
        ! (info by G. Pizzi and A. Cepellotti)
        !
        a1(1) = celldm(1)*(term1-2.0_dp*term2)/3.0_dp
        a1(2) = celldm(1)*(term1+term2)/3.0_dp
        a1(3) = a1(2)
        a2(1) = a1(3)
        a2(2) = a1(1)
        a2(3) = a1(2)
        a3(1) = a1(2)
        a3(2) = a1(3)
        a3(3) = a1(1)
     ENDIF
  ELSEIF (ibrav == 6) THEN
     !
     !     tetragonal lattice
     !
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)*cbya
     !
  ELSEIF (ibrav == 7) THEN
     !
     !     body centered tetragonal lattice
     !
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     cbya=celldm(3)
     a2(1)=celldm(1)/2.d0
     a2(2)=a2(1)
     a2(3)=cbya*celldm(1)/2.d0
     a1(1)= a2(1)
     a1(2)=-a2(1)
     a1(3)= a2(3)
     a3(1)=-a2(1)
     a3(2)=-a2(1)
     a3(3)= a2(3)
     !
  ELSEIF (ibrav == 8) THEN
     !
     !     Simple orthorhombic lattice
     !
     IF (celldm (2) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(2)', ibrav)
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(3)=celldm(1)*celldm(3)
     !
  ELSEIF ( abs(ibrav) == 9) THEN
     !
     !     One face (base) centered orthorhombic lattice  (C type)
     !
     IF (celldm (2) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(2)', &
                                                                 abs(ibrav))
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)', &
                                                                 abs(ibrav))
     !
     IF ( ibrav == 9 ) THEN
        !   old PWscf description
        a1(1) = 0.5d0 * celldm(1)
        a1(2) = a1(1) * celldm(2)
        a2(1) = - a1(1)
        a2(2) = a1(2)
     ELSE
        !   alternate description
        a1(1) = 0.5d0 * celldm(1)
        a1(2) =-a1(1) * celldm(2)
        a2(1) = a1(1)
        a2(2) =-a1(2)
     ENDIF
     a3(3) = celldm(1) * celldm(3)
     !
  ELSEIF ( ibrav == 91 ) THEN
     !
     !     One face (base) centered orthorhombic lattice  (A type)
     !
     IF (celldm (2) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(2)', ibrav)
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     a1(1) = celldm(1)
     a2(2) = celldm(1) * celldm(2) * 0.5_DP
     a2(3) = - celldm(1) * celldm(3) * 0.5_DP
     a3(2) = a2(2)
     a3(3) = - a2(3)
     !
  ELSEIF (ibrav == 10) THEN
     !
     !     All face centered orthorhombic lattice
     !
     IF (celldm (2) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(2)', ibrav)
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     a2(1) = 0.5d0 * celldm(1)
     a2(2) = a2(1) * celldm(2)
     a1(1) = a2(1)
     a1(3) = a2(1) * celldm(3)
     a3(2) = a2(1) * celldm(2)
     a3(3) = a1(3)
     !
  ELSEIF (ibrav == 11) THEN
     !
     !     Body centered orthorhombic lattice
     !
     IF (celldm (2) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(2)', ibrav)
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)', ibrav)
     !
     a1(1) = 0.5d0 * celldm(1)
     a1(2) = a1(1) * celldm(2)
     a1(3) = a1(1) * celldm(3)
     a2(1) = - a1(1)
     a2(2) = a1(2)
     a2(3) = a1(3)
     a3(1) = - a1(1)
     a3(2) = - a1(2)
     a3(3) = a1(3)
     !
  ELSEIF (ibrav == 12) THEN
     !
     !     Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
     !
     IF (celldm (2) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(2)', ibrav)
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)', ibrav)
     IF (abs(celldm(4))>=1.d0) CALL env_errore ('latgen', 'wrong celldm(4)', ibrav)
     !
     sen=sqrt(1.d0-celldm(4)**2)
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(4)
     a2(2)=celldm(1)*celldm(2)*sen
     a3(3)=celldm(1)*celldm(3)
     !
  ELSEIF (ibrav ==-12) THEN
     !
     !     Simple monoclinic lattice, unique axis: b (more common)
     !
     IF (celldm (2) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(2)',-ibrav)
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)',-ibrav)
     IF (abs(celldm(5))>=1.d0) CALL env_errore ('latgen', 'wrong celldm(5)',-ibrav)
     !
     sen=sqrt(1.d0-celldm(5)**2)
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(3)=celldm(1)*celldm(3)*sen
     !
  ELSEIF (ibrav == 13) THEN
     !
     !     One face centered monoclinic lattice unique axis c
     !
     IF (celldm (2) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(2)', ibrav)
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)', ibrav)
     IF (abs(celldm(4))>=1.d0) CALL env_errore ('latgen', 'wrong celldm(4)', ibrav)
     !
     sen = sqrt( 1.d0 - celldm(4) ** 2 )
     a1(1) = 0.5d0 * celldm(1)
     a1(3) =-a1(1) * celldm(3)
     a2(1) = celldm(1) * celldm(2) * celldm(4)
     a2(2) = celldm(1) * celldm(2) * sen
     a3(1) = a1(1)
     a3(3) =-a1(3)
  ELSEIF (ibrav == -13) THEN
     !
     !     One face centered monoclinic lattice unique axis b
     !
     IF (celldm (2) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(2)',-ibrav)
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)',-ibrav)
     IF (abs(celldm(5))>=1.d0) CALL env_errore ('latgen', 'wrong celldm(5)',-ibrav)
     !
     sen = sqrt( 1.d0 - celldm(5) ** 2 )
     a1(1) = 0.5d0 * celldm(1)
     a1(2) =-a1(1) * celldm(2)
     a2(1) = a1(1)
     a2(2) =-a1(2)
     a3(1) = celldm(1) * celldm(3) * celldm(5)
     a3(3) = celldm(1) * celldm(3) * sen
     !
  ELSEIF (ibrav == 14) THEN
     !
     !     Triclinic lattice
     !
     IF (celldm (2) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(2)', ibrav)
     IF (celldm (3) <= 0.d0) CALL env_errore ('latgen', 'wrong celldm(3)', ibrav)
     IF (abs(celldm(4))>=1.d0) CALL env_errore ('latgen', 'wrong celldm(4)', ibrav)
     IF (abs(celldm(5))>=1.d0) CALL env_errore ('latgen', 'wrong celldm(5)', ibrav)
     IF (abs(celldm(6))>=1.d0) CALL env_errore ('latgen', 'wrong celldm(6)', ibrav)
     !
     singam=sqrt(1.d0-celldm(6)**2)
     term= (1.d0+2.d0*celldm(4)*celldm(5)*celldm(6)             &
          -celldm(4)**2-celldm(5)**2-celldm(6)**2)
     IF (term < 0.d0) CALL env_errore &
        ('latgen', 'celldm do not make sense, check your data', ibrav)
     term= sqrt(term/(1.d0-celldm(6)**2))
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(6)
     a2(2)=celldm(1)*celldm(2)*singam
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
     a3(3)=celldm(1)*celldm(3)*term
     !
  ELSE
     !
     CALL env_errore('latgen',' nonexistent bravais lattice',ibrav)
     !
  ENDIF
  !
  !  calculate unit-cell volume omega
  !
  CALL env_volume (1.0_dp, a1, a2, a3, omega)
  !
  RETURN
  !
END SUBROUTINE env_latgen
!
SUBROUTINE env_abc2celldm ( ibrav, a,b,c,cosab,cosac,cosbc, celldm )
  !
  !  returns internal parameters celldm from crystallographics ones
  !
  USE env_kinds,     ONLY: dp
  USE env_constants, ONLY: bohr_radius_angs
  IMPLICIT NONE
  !
  INTEGER,  INTENT (in) :: ibrav
  REAL(DP), INTENT (in) :: a,b,c, cosab, cosac, cosbc
  REAL(DP), INTENT (out) :: celldm(6)
  !
  IF (a <= 0.0_dp) CALL env_errore('abc2celldm','incorrect lattice parameter (a)',1)
  IF (b <  0.0_dp) CALL env_errore('abc2celldm','incorrect lattice parameter (b)',1)
  IF (c <  0.0_dp) CALL env_errore('abc2celldm','incorrect lattice parameter (c)',1)
  IF ( abs (cosab) > 1.0_dp) CALL env_errore('abc2celldm', &
                   'incorrect lattice parameter (cosab)',1)
  IF ( abs (cosac) > 1.0_dp) CALL env_errore('abc2celldm', &
                   'incorrect lattice parameter (cosac)',1)
  IF ( abs (cosbc) > 1.0_dp) CALL env_errore('abc2celldm', &
       'incorrect lattice parameter (cosbc)',1)
  !
  celldm(1) = a / bohr_radius_angs
  celldm(2) = b / a
  celldm(3) = c / a
  !
  IF ( ibrav == 14 .or. ibrav == 0 ) THEN
     !
     ! ... triclinic lattice
     !
     celldm(4) = cosbc
     celldm(5) = cosac
     celldm(6) = cosab
     !
  ELSEIF ( ibrav ==-12 .or. ibrav ==-13 ) THEN
     !
     ! ... monoclinic P or base centered lattice, unique axis b
     !
     celldm(4) = 0.0_dp
     celldm(5) = cosac
     celldm(6) = 0.0_dp
     !
  ELSEIF ( ibrav ==-5 .or. ibrav ==5 .or. ibrav ==12 .or. ibrav ==13 ) THEN
     !
     ! ... trigonal and monoclinic lattices, unique axis c
     !
     celldm(4) = cosab
     celldm(5) = 0.0_dp
     celldm(6) = 0.0_dp
     !
  ELSE
     !
     celldm(4) = 0.0_dp
     celldm(5) = 0.0_dp
     celldm(6) = 0.0_dp
     !
  ENDIF
  !
END SUBROUTINE env_abc2celldm
!
SUBROUTINE env_celldm2abc ( ibrav, celldm, a,b,c,cosab,cosac,cosbc )
  !
  !  returns crystallographic parameters a,b,c from celldm
  !
  USE env_kinds,     ONLY: dp
  USE env_constants, ONLY: bohr_radius_angs
  IMPLICIT NONE
  !
  INTEGER,  INTENT (in) :: ibrav
  REAL(DP), INTENT (in) :: celldm(6)
  REAL(DP), INTENT (out) :: a,b,c, cosab, cosac, cosbc
  !
  !
  a = celldm(1) * bohr_radius_angs
  b = celldm(1)*celldm(2) * bohr_radius_angs
  c = celldm(1)*celldm(3) * bohr_radius_angs
  !
  IF ( ibrav == 14 .or. ibrav == 0 ) THEN
     !
     ! ... triclinic lattice
     !
     cosbc = celldm(4)
     cosac = celldm(5)
     cosab = celldm(6)
     !
  ELSEIF ( ibrav ==-12 .or. ibrav ==-13 ) THEN
     !
     ! ... monoclinic P or base centered lattice, unique axis b
     !
     cosab = 0.0_dp
     cosac = celldm(5)
     cosbc = 0.0_dp
     !
  ELSEIF ( ibrav ==-5 .or. ibrav ==5 .or. ibrav ==12 .or. ibrav ==13 ) THEN
     !
     ! ... trigonal and monoclinic lattices, unique axis c
     !
     cosab = celldm(4)
     cosac = 0.0_dp
     cosbc = 0.0_dp
     !
  ELSE
     cosab = 0.0_dp
     cosac = 0.0_dp
     cosbc = 0.0_dp
  ENDIF
  !
END SUBROUTINE env_celldm2abc