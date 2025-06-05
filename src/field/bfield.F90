#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine bfield(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   B F I E L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   subroutine for update of b-field for half time step    .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xu,yu,zu
  integer(kind=4) :: ps


!-------------------- test particle simulation
      if(juncan.ge.1000) return


!-------------------- 
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- loop for update of the b-field
      do k=-1,zu
      do j=-1,yu
      do i=-1,xu
        eb(BX,i,j,k,ps) = eb(BX,i  ,j  ,k  ,ps) &
       &                + eb(EY,i  ,j  ,k+1,ps) - eb(EY,i,j,k,ps) &
       &                - eb(EZ,i  ,j+1,k  ,ps) + eb(EZ,i,j,k,ps)
        eb(BY,i,j,k,ps) = eb(BY,i  ,j  ,k  ,ps) &
       &                + eb(EZ,i+1,j  ,k  ,ps) - eb(EZ,i,j,k,ps) &
       &                - eb(EX,i  ,j  ,k+1,ps) + eb(EX,i,j,k,ps)
        eb(BZ,i,j,k,ps) = eb(BZ,i  ,j  ,k  ,ps) &
       &                + eb(EX,i  ,j+1,k  ,ps) - eb(EX,i,j,k,ps) &
       &                - eb(EY,i+1,j  ,k  ,ps) + eb(EY,i,j,k,ps)
      end do
      end do
      end do


!-------------------- boundary treatment
!      call fsmask(1)
!      call fbound(1)


  return
  end subroutine
