!
#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine efield(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   E F I E L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   subroutine to update e-field by current contribution   .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i,j,k, ig
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ps
  real(kind=8) :: dlt
  real(kind=8) :: capc1


!-------------------- test particle simulation
      if(nflag_testp.ne.1) then
        dlt = 2.0d0
      else
        dlt = 0.0d0
      end if


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- 
      do ig=1,ngap
        if(i_gap(ig).ge.xl.and.i_gap(ig).le.xu.and. &
       &   j_gap(ig).ge.yl.and.j_gap(ig).le.yu.and. &
       &   k_gap(ig).ge.zl.and.k_gap(ig).le.zu) then
          ezgap(ig) = eb(EZ,i_gap(ig)-xl,j_gap(ig)-yl,k_gap(ig)-zl,ps)
        end if
      end do


!-------------------- loop for update of the e-field
      do k=0,zu-zl
      do j=0,yu-yl
      do i=0,xu-xl
        eb(EX,i,j,k,ps) = eb(EX,i,j,k,ps) &
       &                + tcs*( eb(BZ,i,j,k,ps) - eb(BZ,i  ,j-1,k  ,ps) &
       &                      - eb(BY,i,j,k,ps) + eb(BY,i  ,j  ,k-1,ps) ) &
       &                - dlt*aj(JX,i,j,k,ps)
        eb(EY,i,j,k,ps) = eb(EY,i,j,k,ps) &
       &                + tcs*( eb(BX,i,j,k,ps) - eb(BX,i  ,j  ,k-1,ps) &
       &                      - eb(BZ,i,j,k,ps) + eb(BZ,i-1,j  ,k  ,ps) ) &
       &                - dlt*aj(JY,i,j,k,ps)
        eb(EZ,i,j,k,ps) = eb(EZ,i,j,k,ps) &
       &                + tcs*( eb(BY,i,j,k,ps) - eb(BY,i-1,j  ,k  ,ps) &
       &                      - eb(BX,i,j,k,ps) + eb(BX,i  ,j-1,k  ,ps) ) &
       &                - dlt*aj(JZ,i,j,k,ps)
      end do
      end do
      end do


!-------------------- modeling of lumped elements
      if(mode_dipole.eq.3) then
        do ig=1,ngap
          i = i_gap(ig) - xl
          j = j_gap(ig) - yl
          k = k_gap(ig) - zl
          if(i.ge.0.and.i.le.xu-xl.and. &
         &   j.ge.0.and.j.le.yu-yl.and. &
         &   k.ge.0.and.k.le.zu-zl) then
            capc1 = capc(ig) + 1.0d0
            if(resc(ig).lt.0) then
              eb(EZ,i,j,k,ps) = ezgap(ig) &
             &                + tcs*( eb(BY,i,j,k,ps) - eb(BY,i-1,j,k,ps) &
             &                - eb(BX,i,j,k,ps) + eb(BX,i,j-1,k,ps) )/capc1 &
             &                - dlt*aj(JZ,i,j,k,ps)/capc1
            else
              eb(EZ,i,j,k,ps) = ezgap(ig)*(resc(ig)*capc1 - 1.0d0) &
             &                 /(resc(ig)*capc1 + 1.0d0) &
             &                + tcs*resc(ig) &
             &                 *( eb(BY,i,j,k,ps) - eb(BY,i-1,j,k,ps) &
             &                  - eb(BX,i,j,k,ps) + eb(BX,i,j-1,k,ps) ) &
             &                 /(resc(ig)*capc1 + 1.0d0) &
             &                -dlt*resc(ig)*aj(JZ,i,j,k,ps) &
             &                 /(resc(ig)*capc1 + 1.0d0)
            end if
          end if
        end do
      end if


  return
  end subroutine

