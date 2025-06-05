#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine rhoupd(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   R H O U P D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a charge distribution to      .
!   .      grid points from particle locations by the          .
!   .      "area-sharing" or "area-weighting" scheme.          .
!   ............................................................
!
!     *"SUBROUTINE CHARGP" is optimized for scalar parallel
!       processors (thread-parallel capabilities). 
!
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xu,yu,zu
  integer(kind=4) :: ps


!-------------------- 
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- continuity equation for charge
      do k=0,zu
      do j=0,yu
      do i=0,xu
        rho(1,i,j,k,ps) = rho(1,i,j,k,ps) &
       &                - ( aj(JX,i,j,k,ps) - aj(JX,i-1,j,k,ps) &
       &                  + aj(JY,i,j,k,ps) - aj(JY,i,j-1,k,ps) &
       &                  + aj(JZ,i,j,k,ps) - aj(JZ,i,j,k-1,ps) )*2.0d0
      end do
      end do
      end do


  return
  end subroutine
