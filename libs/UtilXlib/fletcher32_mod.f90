
! Copyright (C) 2002-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module contains interfaces to call the c-routine clib/fletcher32.c
! implementing the Fletcher-32 checksum algorithm as reported in
! https://en.wikipedia.org/wiki/Fletcher%27s_checksum#Optimizations
!
! SdG September 3rd 2017
!
!------------------------------------------------------------------------------!
    MODULE env_fletcher32_mod
!------------------------------------------------------------------------------!
    USE env_util_param,     ONLY : DP
    !
    IMPLICIT NONE
    PRIVATE
    integer(2) :: dat(1)

    PUBLIC :: env_fletcher32_cksum, env_fletcher32
!
    INTERFACE env_fletcher32_cksum
       MODULE PROCEDURE env_fletcher32_i1, env_fletcher32_r1, env_fletcher32_c1, env_fletcher32_z,  env_fletcher32_l,  &
                        env_fletcher32_iv, env_fletcher32_rv, env_fletcher32_cv, env_fletcher32_zv, env_fletcher32_lv, &
                        env_fletcher32_im, env_fletcher32_rm, env_fletcher32_cm,                env_fletcher32_lm, &
                        env_fletcher32_it, env_fletcher32_rt, env_fletcher32_ct, &
                        env_fletcher32_i4, env_fletcher32_r4, env_fletcher32_c4, &
                        env_fletcher32_r5, env_fletcher32_c5
    END INTERFACE

    INTERFACE
       FUNCTION env_fletcher32( dat, dat_size ) BIND(C,name="fletcher32") RESULT(t)
          USE ISO_C_BINDING
          integer(kind=c_int16_t) :: dat(*)
          integer(kind=c_int32_t) :: dat_size
          integer(kind=c_int32_t) :: t
       END FUNCTION env_fletcher32
    END INTERFACE
!
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!

!..fletcher32_cksum
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_i1(msg, cksum)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: msg
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_i1
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_iv(msg, cksum)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: msg(:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_iv
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_im( msg, cksum )
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: msg(:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_im
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_it( msg, cksum )
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: msg(:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_it
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_i4(msg, cksum )
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: msg(:,:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_i4
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_r1( msg, cksum  )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_r1
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_rv(msg, cksum )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg(:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_rv
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_rm(msg, cksum )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg(:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_rm
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_rt(msg, cksum )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg(:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_rt
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_r4(msg, cksum )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg(:,:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_r4
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_r5(msg, cksum )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg(:,:,:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_r5
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_c1(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_c1
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_cv(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg(:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_cv
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_cm(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg(:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_cm
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_ct(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg(:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_ct
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_c4(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg(:,:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_c4
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_c5(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg(:,:,:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_c5
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_l(msg, cksum )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: msg
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_l
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_lv(msg, cksum )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: msg(:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_lv
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_lm(msg, cksum )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: msg(:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_lm
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_z(msg, cksum )
         IMPLICIT NONE
         CHARACTER(len=*), INTENT(IN) :: msg
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_z
!
!------------------------------------------------------------------------------!
      SUBROUTINE env_fletcher32_zv(msg, cksum )
         IMPLICIT NONE
         CHARACTER(len=*), INTENT(IN) :: msg(:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = env_fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE env_fletcher32_zv
!
!------------------------------------------------------------------------------!
    END MODULE env_fletcher32_mod
!------------------------------------------------------------------------------!
