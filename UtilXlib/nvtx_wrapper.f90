!MIT License

!Copyright (c) 2019 maxcuda

!This module has been downloaded and adapted from
!   https://github.com/maxcuda/NVTX_example     
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.


! ----
! env_nvtx
! ----

module env_nvtx
  use iso_c_binding
#ifdef __CUDA
  use cudafor
#endif
  implicit none
#ifdef __PROFILE_NVTX
  integer,private :: col(7) = [ Z'0000ff00', Z'000000ff', Z'00ffff00',Z'00ff00ff',Z'0000ffff', &
                                Z'00ff0000', Z'00ffffff']
  character(len=256),private :: tempName
!  logical, save :: __PROFILE_NVTX=.false.
  type, bind(C):: env_nvtxEventAttributes
     integer(C_INT16_T):: version=1
     integer(C_INT16_T):: size=48 !
     integer(C_INT):: category=0
     integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
     integer(C_INT):: color
     integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
     integer(C_INT):: reserved0
     integer(C_INT64_T):: payload   ! union uint,int,double
     integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1 
     type(C_PTR):: message  ! ascii char
  end type env_nvtxEventAttributes

  interface env_nvtxRangePush
     ! push range with custom label and standard color
     subroutine env_nvtxRangePushA(name) bind(C, name='env_nvtxRangePushA')
       use iso_c_binding
       character(kind=C_CHAR,len=*) :: name
     end subroutine env_nvtxRangePushA

     ! push range with custom label and custom color
     subroutine env_nvtxRangePushEx(event) bind(C, name='env_nvtxRangePushEx')
       use iso_c_binding
       import:: env_nvtxEventAttributes
       type(env_nvtxEventAttributes):: event
     end subroutine env_nvtxRangePushEx
  end interface env_nvtxRangePush

  interface env_nvtxRangePop
     subroutine env_nvtxRangePop() bind(C, name='env_nvtxRangePop')
     end subroutine env_nvtxRangePop
  end interface env_nvtxRangePop
#endif

contains

  subroutine env_nvtxStartRange(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
#ifdef __PROFILE_NVTX
    type(env_nvtxEventAttributes):: event
#if defined(__CUDA) && defined(__SYNC_NVPROF)
    integer :: istat
    istat = cudaDeviceSynchronize()
#endif

    tempName=trim(name)//c_null_char

    if ( .not. present(id)) then
       call env_nvtxRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call env_nvtxRangePushEx(event)
    end if
#endif
  end subroutine env_nvtxStartRange

  subroutine env_nvtxStartRangeAsync(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
#ifdef __PROFILE_NVTX
    type(env_nvtxEventAttributes):: event

    tempName=trim(name)//c_null_char

    if ( .not. present(id)) then
       call env_nvtxRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call env_nvtxRangePushEx(event)
    end if
#endif
  end subroutine env_nvtxStartRangeAsync


  subroutine env_nvtxEndRange
#ifdef __PROFILE_NVTX
#if defined(__CUDA) && defined(__SYNC_NVPROF)
    integer :: istat
    istat = cudaDeviceSynchronize()
#endif
    call env_nvtxRangePop
#endif
  end subroutine env_nvtxEndRange

  subroutine env_nvtxEndRangeAsync
#ifdef __PROFILE_NVTX
    call env_nvtxRangePop
#endif
  end subroutine env_nvtxEndRangeAsync

end module env_nvtx
