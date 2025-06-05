#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine fsmask(ps,func)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   F S M A S K
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .          if func = 0 then b-field                        .
!   .          if func = 1 then e-field (transverse)           .
!   .          if func = 2 then e-field (longitudinal)         .
!   .          if func = 3 then current                        .
!   .          if func = 4 then charge                         .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xl,yl,zl, xu,yu,zu
  integer(kind=4) :: xlmask,ylmask,zlmask, xumask,yumask,zumask
  integer(kind=4) :: func
  integer(kind=4) :: ps
  real(kind=8),external :: mskfnx,mskfny,mskfnz


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)
      xumask = medges(1,1,sdid(ps)+1); xlmask = medges(2,1,sdid(ps)+1)
      yumask = medges(1,2,sdid(ps)+1); ylmask = medges(2,2,sdid(ps)+1)
      zumask = medges(1,3,sdid(ps)+1); zlmask = medges(2,3,sdid(ps)+1)


!-------------------- magnetic field
      if(func.eq.0) then
        if(nfbnd(1).eq.2) then
          do k=-1,zu-zl
          do j=-1,yu-yl
          do i=-1,xumask-xl
!            eb(BX,i,j,k,ps) = eb(BX,i,j,k,ps)*mskfnx(i+xl+0.0d0)
            eb(BY,i,j,k,ps) = eb(BY,i,j,k,ps)*mskfnx(i+xl+0.5d0)
            eb(BZ,i,j,k,ps) = eb(BZ,i,j,k,ps)*mskfnx(i+xl+0.5d0)
          end do
          end do
          end do
          do k=-1,zu-zl
          do j=-1,yu-yl
          do i=xlmask-1-xl,xu-xl
!            eb(BX,i,j,k,ps) = eb(BX,i,j,k,ps)*mskfnx(i+xl+0.0d0)
            eb(BY,i,j,k,ps) = eb(BY,i,j,k,ps)*mskfnx(i+xl+0.5d0)
            eb(BZ,i,j,k,ps) = eb(BZ,i,j,k,ps)*mskfnx(i+xl+0.5d0)
          end do
          end do
          end do
        end if
!
        if(nfbnd(2).eq.2) then
          do k=-1,zu-zl
          do j=-1,yumask-yl
          do i=-1,xu-xl
            eb(BX,i,j,k,ps) = eb(BX,i,j,k,ps)*mskfny(j+yl+0.5d0)
!            eb(BY,i,j,k,ps) = eb(BY,i,j,k,ps)*mskfny(j+yl+0.0d0)
            eb(BZ,i,j,k,ps) = eb(BZ,i,j,k,ps)*mskfny(j+yl+0.5d0)
          end do
          end do
          end do
          do k=-1,zu-zl
          do j=ylmask-1-yl,yu-yl
          do i=-1,xu-xl
            eb(BX,i,j,k,ps) = eb(BX,i,j,k,ps)*mskfny(j+yl+0.5d0)
!            eb(BY,i,j,k,ps) = eb(BY,i,j,k,ps)*mskfny(j+yl+0.0d0)
            eb(BZ,i,j,k,ps) = eb(BZ,i,j,k,ps)*mskfny(j+yl+0.5d0)
          end do
          end do
          end do
        end if
!
        if(nfbnd(3).eq.2) then
          do k=-1,zumask-zl
          do j=-1,yu-yl
          do i=-1,xu-xl
            eb(BX,i,j,k,ps) = eb(BX,i,j,k,ps)*mskfnz(k+zl+0.5d0)
            eb(BY,i,j,k,ps) = eb(BY,i,j,k,ps)*mskfnz(k+zl+0.5d0)
!            eb(BZ,i,j,k,ps) = eb(BZ,i,j,k,ps)*mskfnz(k+zl+0.0d0)
          end do
          end do
          end do
          do k=zlmask-1-zl,zu-zl
          do j=-1,yu-yl
          do i=-1,xu-xl
            eb(BX,i,j,k,ps) = eb(BX,i,j,k,ps)*mskfnz(k+zl+0.5d0)
            eb(BY,i,j,k,ps) = eb(BY,i,j,k,ps)*mskfnz(k+zl+0.5d0)
!            eb(BZ,i,j,k,ps) = eb(BZ,i,j,k,ps)*mskfnz(k+zl+0.0d0)
          end do
          end do
          end do
        end if


!-------------------- electric field
      else if(func.eq.1) then
        if(nfbnd(1).eq.2) then
          do k=0,zu-zl
          do j=0,yu-yl
          do i=0,xumask-xl
            eb(EY,i,j,k,ps) = eb(EY,i,j,k,ps)*mskfnx(i+xl+0.0d0)
            eb(EZ,i,j,k,ps) = eb(EZ,i,j,k,ps)*mskfnx(i+xl+0.0d0)
          end do
          end do
          end do
          do k=0,zu-zl
          do j=0,yu-yl
          do i=xlmask-xl,xu-xl
            eb(EY,i,j,k,ps) = eb(EY,i,j,k,ps)*mskfnx(i+xl+0.0d0)
            eb(EZ,i,j,k,ps) = eb(EZ,i,j,k,ps)*mskfnx(i+xl+0.0d0)
          end do
          end do
          end do
        end if
!
        if(nfbnd(2).eq.2) then
          do k=0,zu-zl
          do j=0,yumask-yl
          do i=0,xu-xl
            eb(EX,i,j,k,ps) = eb(EX,i,j,k,ps)*mskfny(j+yl+0.0d0)
            eb(EZ,i,j,k,ps) = eb(EZ,i,j,k,ps)*mskfny(j+yl+0.0d0)
          end do
          end do
          end do
          do k=0,zu-zl
          do j=ylmask-yl,yu-yl
          do i=0,xu-xl
            eb(EX,i,j,k,ps) = eb(EX,i,j,k,ps)*mskfny(j+yl+0.0d0)
            eb(EZ,i,j,k,ps) = eb(EZ,i,j,k,ps)*mskfny(j+yl+0.0d0)
          end do
          end do
          end do
        end if
!
        if(nfbnd(3).eq.2) then
          do k=0,zumask-zl
          do j=0,yu-yl
          do i=0,xu-xl
            eb(EX,i,j,k,ps) = eb(EX,i,j,k,ps)*mskfnz(k+zl+0.0d0)
            eb(EY,i,j,k,ps) = eb(EY,i,j,k,ps)*mskfnz(k+zl+0.0d0)
          end do
          end do
          end do
          do k=zlmask-zl,zu-zl
          do j=0,yu-yl
          do i=0,xu-xl
            eb(EX,i,j,k,ps) = eb(EX,i,j,k,ps)*mskfnz(k+zl+0.0d0)
            eb(EY,i,j,k,ps) = eb(EY,i,j,k,ps)*mskfnz(k+zl+0.0d0)
          end do
          end do
          end do
        end if


!-------------------- electric field (longitudinal)
      else if(func.eq.2) then
        if(nfbnd(1).eq.2) then
          do k=0,zu-zl
          do j=0,yu-yl
          do i=0,xumask-xl
            eb(EX,i,j,k,ps) = eb(EX,i,j,k,ps)*mskfnx(i+xl+0.5d0)
          end do
          end do
          end do
          do k=0,zu-zl
          do j=0,yu-yl
          do i=xlmask-xl,xu-xl
            eb(EX,i,j,k,ps) = eb(EX,i,j,k,ps)*mskfnx(i+xl+0.5d0)
          end do
          end do
          end do
        end if
!
        if(nfbnd(2).eq.2) then
          do k=0,zu-zl
          do j=0,yumask-yl
          do i=0,xu-xl
            eb(EY,i,j,k,ps) = eb(EY,i,j,k,ps)*mskfny(j+yl+0.5d0)
          end do
          end do
          end do
          do k=0,zu-zl
          do j=ylmask-yl,yu-yl
          do i=0,xu-xl
            eb(EY,i,j,k,ps) = eb(EY,i,j,k,ps)*mskfny(j+yl+0.5d0)
          end do
          end do
          end do
        end if
!
        if(nfbnd(3).eq.2) then
          do k=0,zumask-zl
          do j=0,yu-yl
          do i=0,xu-xl
            eb(EZ,i,j,k,ps) = eb(EZ,i,j,k,ps)*mskfnz(k+zl+0.5d0)
          end do
          end do
          end do
          do k=zlmask-zl,zu-zl
          do j=0,yu-yl
          do i=0,xu-xl
            eb(EZ,i,j,k,ps) = eb(EZ,i,j,k,ps)*mskfnz(k+zl+0.5d0)
          end do
          end do
          end do
        end if


!-------------------- current density
      else if(func.eq.3) then
        if(nfbnd(1).eq.2) then
          do k=0,zu-zl
          do j=0,yu-yl
          do i=0,xumask-xl
            aj(JX,i,j,k,ps) = aj(JX,i,j,k,ps)*mskfnx(i+xl+0.5d0)
          end do
          end do
          end do
          do k=0,zu-zl
          do j=0,yu-yl
          do i=xlmask-xl,xu-xl
            aj(JX,i,j,k,ps) = aj(JX,i,j,k,ps)*mskfnx(i+xl+0.5d0)
          end do
          end do
          end do
        end if
!
        if(nfbnd(2).eq.2) then
          do k=0,zu-zl
          do j=0,yumask-yl
          do i=0,xu-xl
            aj(JY,i,j,k,ps) = aj(JY,i,j,k,ps)*mskfny(j+yl+0.5d0)
          end do
          end do
          end do
          do k=0,zu-zl
          do j=ylmask-yl,yu-yl
          do i=0,xu-xl
            aj(JY,i,j,k,ps) = aj(JY,i,j,k,ps)*mskfny(j+yl+0.5d0)
          end do
          end do
          end do
        end if
!
        if(nfbnd(3).eq.2) then
          do k=0,zumask-zl
          do j=0,yu-yl
          do i=0,xu-xl
            aj(JZ,i,j,k,ps) = aj(JZ,i,j,k,ps)*mskfnz(k+zl+0.5d0)
          end do
          end do
          end do
          do k=zlmask-zl,zu-zl
          do j=0,yu-yl
          do i=0,xu-xl
            aj(JZ,i,j,k,ps) = aj(JZ,i,j,k,ps)*mskfnz(k+zl+0.5d0)
          end do
          end do
          end do
        end if


!-------------------- charge density
      else if(func.eq.4) then
        if(nfbnd(1).eq.2) then
          do k=0,zu-zl
          do j=0,yu-yl
          do i=0,xumask-xl
            rho(1,i,j,k,ps) = rho(1,i,j,k,ps)*mskfnx(i+xl+0.0d0)
            rhodg(:,i,j,k,ps) = rhodg(:,i,j,k,ps)*mskfnx(i+xl+0.0d0)
          end do
          end do
          end do
          do k=0,zu-zl
          do j=0,yu-yl
          do i=xlmask-xl,xu-xl
            rho(1,i,j,k,ps) = rho(1,i,j,k,ps)*mskfnx(i+xl+0.0d0)
            rhodg(:,i,j,k,ps) = rhodg(:,i,j,k,ps)*mskfnx(i+xl+0.0d0)
          end do
          end do
          end do
        end if
!
        if(nfbnd(2).eq.2) then
          do k=0,zu-zl
          do j=0,yumask-yl
          do i=0,xu-xl
            rho(1,i,j,k,ps) = rho(1,i,j,k,ps)*mskfny(j+yl+0.0d0)
            rhodg(:,i,j,k,ps) = rhodg(:,i,j,k,ps)*mskfny(j+yl+0.0d0)
          end do
          end do
          end do
          do k=0,zu-zl
          do j=ylmask-yl,yu-yl
          do i=0,xu-xl
            rho(1,i,j,k,ps) = rho(1,i,j,k,ps)*mskfny(j+yl+0.0d0)
            rhodg(:,i,j,k,ps) = rhodg(:,i,j,k,ps)*mskfny(j+yl+0.0d0)
          end do
          end do
          end do
        end if
!
        if(nfbnd(3).eq.2) then
          do k=0,zumask -zl
          do j=0,yu-yl
          do i=0,xu-xl
            rho(1,i,j,k,ps) = rho(1,i,j,k,ps)*mskfnz(k+zl+0.0d0)
            rhodg(:,i,j,k,ps) = rhodg(:,i,j,k,ps)*mskfnz(k+zl+0.0d0)
          end do
          end do
          end do
          do k=zlmask-zl,zu-zl
          do j=0,yu-yl
          do i=0,xu-xl
            rho(1,i,j,k,ps) = rho(1,i,j,k,ps)*mskfnz(k+zl+0.0d0)
            rhodg(:,i,j,k,ps) = rhodg(:,i,j,k,ps)*mskfnz(k+zl+0.0d0)
          end do
          end do
          end do
        end if


!-------------------- end of masking
      end if


  return
  end subroutine



!_._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._
!
! ============================================
  double precision function mskfnx(pos)
! ============================================
!
!-------------------- 
  use oh_type
  use paramt
  use allcom
  implicit none
!
  real(kind=8) :: pos
!  real(kind=8) :: mskfnx


!-------------------- 
      if(pos.le.0.0d0) then
        mskfnx = 0.0d0
      else if(pos.ge.slx) then
        mskfnx = 0.0d0
      else if(pos.lt.nxl0) then
        mskfnx = 1.0d0-((nxl0-pos)/nxl)**2
      else if(pos.gt.nxr0) then
        mskfnx = 1.0d0-((pos-nxr0)/nxr)**2
      else
        mskfnx = 1.0d0
      end if


  end function


!
! ============================================
  double precision function mskfny(pos)
! ============================================
!
!-------------------- 
  use oh_type
  use paramt
  use allcom
  implicit none
!
  real(kind=8) :: pos
!  real(kind=8) :: mskfny


!-------------------- 
      if(pos.le.0.0d0) then
        mskfny = 0.0d0
      else if(pos.ge.sly) then
        mskfny = 0.0d0
      else if(pos.lt.nyl0) then
        mskfny = 1.0d0-((nyl0-pos)/nyl)**2
      else if(pos.gt.nyr0) then
        mskfny = 1.0d0-((pos-nyr0)/nyr)**2
      else
        mskfny = 1.0d0
      end if


  end function


!
! ============================================
  double precision function mskfnz(pos)
! ============================================
!
!-------------------- 
  use oh_type
  use paramt
  use allcom
  implicit none
!
  real(kind=8) :: pos
!  real(kind=8) :: mskfnz


!-------------------- 
      if(pos.le.0.0d0) then
        mskfnz = 0.0d0
      else if(pos.ge.slz) then
        mskfnz = 0.0d0
      else if(pos.lt.nzl0) then
        mskfnz = 1.0d0-((nzl0-pos)/nzl)**2
      else if(pos.gt.nzr0) then
        mskfnz = 1.0d0-((pos-nzr0)/nzr)**2
      else
        mskfnz = 1.0d0
      end if


  end function
