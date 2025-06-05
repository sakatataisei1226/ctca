#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine medium
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   M E D I U M
!   ____________________________________________________________
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  interface rfold
    integer(kind=4) function rfold_i(r,nr,dir)
      integer(kind=4),intent(in) :: r, nr
      integer(kind=4),intent(in) :: dir
    end function rfold_i
!
    real(kind=4) function rfold_r(r,nr,dir)
      real(kind=4),intent(in) :: r, nr
      integer(kind=4),intent(in) :: dir
    end function rfold_r
!
    real(kind=8) function rfold_d(r,nr,dir)
      real(kind=8),intent(in) :: r, nr
      integer(kind=4),intent(in) :: dir
    end function rfold_d
  end interface rfold
!
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ipc
  real(kind=8) :: x,y,z
  real(kind=8) :: xlbd,xubd,ylbd,yubd,zlbd,zubd
  real(kind=8) :: cylaln, cylrad, cylaxs1,cylaxs2, cyledg1,cyledg2
  real(kind=8) :: sphrad, sphcnt1,sphcnt2,sphcnt3
  real(kind=8) :: radsq,radsqin


!--------------------
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)


!-------------------- 
      do ipc=npc,1,-1
        xlbd = xlpc(ipc); xubd = xupc(ipc)
        ylbd = ylpc(ipc); yubd = yupc(ipc)
        zlbd = zlpc(ipc); zubd = zupc(ipc)
!
        cylaln = cylinder(ipc)%align
        cylrad = cylinder(ipc)%radius
        cylaxs1 = cylinder(ipc)%axis(1)
        cylaxs2 = cylinder(ipc)%axis(2)
        cyledg1 = cylinder(ipc)%edge(1)
        cyledg2 = cylinder(ipc)%edge(2)
!
        sphrad = sphere(ipc)%radius
        sphcnt1 = sphere(ipc)%center(1)
        sphcnt2 = sphere(ipc)%center(2)
        sphcnt3 = sphere(ipc)%center(3)
!
!GEOTYPE0 or 1
        if(geotype(ipc).eq.0.or.geotype(ipc).eq.1) then
          do k=-1,zu-zl+1
          do j=-1,yu-zl+1
          do i=-1,xu-xl+1
            x = (i + xl + 0.5d0)*dr
            y = (j + yl)*dr
            z = (k + zl)*dr
            if(x.gt.xlbd.and.x.lt.xubd.and. &
           &   y.ge.ylbd.and.y.le.yubd.and. &
           &   z.ge.zlbd.and.z.le.zubd) then
              mp(EX,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl + 0.5d0)*dr
            y = (rfold(j,ny,2) + yl)*dr
            z = (rfold(k,nz,3) + zl)*dr
            if(x.gt.xlbd.and.x.lt.xubd.and. &
           &   y.ge.ylbd.and.y.le.yubd.and. &
           &   z.ge.zlbd.and.z.le.zubd) then
              mp(EX,i,j,k,1) = 0.0d0
            end if
!
            x = (i + xl)*dr
            y = (j + yl + 0.5d0)*dr
            z = (k + zl)*dr
            if(x.ge.xlbd.and.x.le.xubd.and. &
           &   y.gt.ylbd.and.y.lt.yubd.and. &
           &   z.ge.zlbd.and.z.le.zubd) then
              mp(EY,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl)*dr
            y = (rfold(j,ny,2) + yl + 0.5d0)*dr
            z = (rfold(k,nz,3) + zl)*dr
            if(x.ge.xlbd.and.x.le.xubd.and. &
           &   y.gt.ylbd.and.y.lt.yubd.and. &
           &   z.ge.zlbd.and.z.le.zubd) then
              mp(EY,i,j,k,1) = 0.0d0
            end if
!
            x = (i + xl)*dr
            y = (j + yl)*dr
            z = (k + zl + 0.5d0)*dr
            if(x.ge.xlbd.and.x.le.xubd.and. &
           &   y.ge.ylbd.and.y.le.yubd.and. &
           &   z.gt.zlbd.and.z.lt.zubd) then
              mp(EZ,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl)*dr
            y = (rfold(j,ny,2) + yl)*dr
            z = (rfold(k,nz,3) + zl + 0.5d0)*dr
            if(x.ge.xlbd.and.x.le.xubd.and. &
           &   y.ge.ylbd.and.y.le.yubd.and. &
           &   z.gt.zlbd.and.z.lt.zubd) then
              mp(EZ,i,j,k,1) = 0.0d0
            end if
!
            x = (i + xl)*dr
            y = (j + yl + 0.5d0)*dr
            z = (k + zl + 0.5d0)*dr
            if(x.ge.xlbd.and.x.le.xubd.and. &
           &   y.gt.ylbd.and.y.lt.yubd.and. &
           &   z.gt.zlbd.and.z.lt.zubd) then
              mp(BX,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl)*dr
            y = (rfold(j,ny,2) + yl + 0.5d0)*dr
            z = (rfold(k,nz,3) + zl + 0.5d0)*dr
            if(x.ge.xlbd.and.x.le.xubd.and. &
           &   y.gt.ylbd.and.y.lt.yubd.and. &
           &   z.gt.zlbd.and.z.lt.zubd) then
              mp(BX,i,j,k,1) = 0.0d0
            end if
!
            x = (i + xl + 0.5d0)*dr
            y = (j + yl)*dr
            z = (k + zl + 0.5d0)*dr
            if(x.gt.xlbd.and.x.lt.xubd.and. &
           &   y.ge.ylbd.and.y.le.yubd.and. &
           &   z.gt.zlbd.and.z.lt.zubd) then
              mp(BY,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl + 0.5d0)*dr
            y = (rfold(j,ny,2) + yl)*dr
            z = (rfold(k,nz,3) + zl + 0.5d0)*dr
            if(x.gt.xlbd.and.x.lt.xubd.and. &
           &   y.ge.ylbd.and.y.le.yubd.and. &
           &   z.gt.zlbd.and.z.lt.zubd) then
              mp(BY,i,j,k,1) = 0.0d0
            end if
!
            x = (i + xl + 0.5d0)*dr
            y = (j + yl + 0.5d0)*dr
            z = (k + zl)*dr
            if(x.gt.xlbd.and.x.lt.xubd.and. &
           &   y.gt.ylbd.and.y.lt.yubd.and. &
           &   z.ge.zlbd.and.z.le.zubd) then
              mp(BZ,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl + 0.5d0)*dr
            y = (rfold(j,ny,2) + yl + 0.5d0)*dr
            z = (rfold(k,nz,3) + zl)*dr
            if(x.gt.xlbd.and.x.lt.xubd.and. &
           &   y.gt.ylbd.and.y.lt.yubd.and. &
           &   z.ge.zlbd.and.z.le.zubd) then
              mp(BZ,i,j,k,1) = 0.0d0
            end if
          end do
          end do
          end do
!
!GEOTYPE2
        else if(geotype(ipc).eq.2) then
          radsq = cylrad**2
          radsqin = (cylrad - 1.0d0)**2
!GEOTYPE2-1
          if(cylaln.eq.1) then
            do k=-1,zu-zl+1
            do j=-1,yu-zl+1
            do i=-1,xu-xl+1
              x = (i + xl + 0.5d0)*dr
              y = (j + yl)*dr
              z = (k + zl)*dr
              if(x.gt.cyledg1.and.x.lt.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl + 0.5d0)*dr
              y = (rfold(j,ny,2) + yl)*dr
              z = (rfold(k,nz,3) + zl)*dr
              if(x.gt.cyledg1.and.x.lt.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl)*dr
              y = (j + yl + 0.5d0)*dr
              z = (k + zl)*dr
              if(x.ge.cyledg1.and.x.le.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl)*dr
              y = (rfold(j,ny,2) + yl + 0.5d0)*dr
              z = (rfold(k,nz,3) + zl)*dr
              if(x.ge.cyledg1.and.x.le.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl)*dr
              y = (j + yl)*dr
              z = (k + zl + 0.5d0)*dr
              if(x.ge.cyledg1.and.x.le.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl)*dr
              y = (rfold(j,ny,2) + yl)*dr
              z = (rfold(k,nz,3) + zl + 0.5d0)*dr
              if(x.ge.cyledg1.and.x.le.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl)*dr
              y = (j + yl + 0.5d0)*dr
              z = (k + zl + 0.5d0)*dr
              if(x.ge.cyledg1.and.x.le.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.lt.radsq) then
                mp(BX,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl)*dr
              y = (rfold(j,ny,2) + yl + 0.5d0)*dr
              z = (rfold(k,nz,3) + zl + 0.5d0)*dr
              if(x.ge.cyledg1.and.x.le.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.lt.radsq) then
                mp(BX,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl + 0.5d0)*dr
              y = (j + yl)*dr
              z = (k + zl + 0.5d0)*dr
              if(x.gt.cyledg1.and.x.lt.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.lt.radsq) then
                mp(BY,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl + 0.5d0)*dr
              y = (rfold(j,ny,2) + yl)*dr
              z = (rfold(k,nz,3) + zl + 0.5d0)*dr
              if(x.gt.cyledg1.and.x.lt.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.lt.radsq) then
                mp(BY,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl + 0.5d0)*dr
              y = (j + yl + 0.5d0)*dr
              z = (k + zl)*dr
              if(x.gt.cyledg1.and.x.lt.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.lt.radsq) then
                mp(BZ,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl + 0.5d0)*dr
              y = (rfold(j,ny,2) + yl + 0.5d0)*dr
              z = (rfold(k,nz,3) + zl)*dr
              if(x.gt.cyledg1.and.x.lt.cyledg2.and. &
             &   (y-cylaxs1)**2+(z-cylaxs2)**2.lt.radsq) then          
               mp(BZ,i,j,k,1) = 0.0d0
              end if
            end do
            end do
            end do
!GEOTYPE2-2
          else if(cylaln.eq.2) then
            do k=-1,zu-zl+1
            do j=-1,yu-zl+1
            do i=-1,xu-xl+1
              x = (i + xl + 0.5d0)*dr
              y = (j + yl)*dr
              z = (k + zl)*dr
              if(y.ge.cyledg1.and.y.le.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl + 0.5d0)*dr
              y = (rfold(j,ny,2) + yl)*dr
              z = (rfold(k,nz,3) + zl)*dr
              if(y.ge.cyledg1.and.y.le.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl)*dr
              y = (j + yl + 0.5d0)*dr
              z = (k + zl)*dr
              if(y.gt.cyledg1.and.y.lt.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl)*dr
              y = (rfold(j,ny,2) + yl + 0.5d0)*dr
              z = (rfold(k,nz,3) + zl)*dr
              if(y.gt.cyledg1.and.y.lt.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl)*dr
              y = (j + yl)*dr
              z = (k + zl + 0.5d0)*dr
              if(y.ge.cyledg1.and.y.le.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl)*dr
              y = (rfold(j,ny,2) + yl)*dr
              z = (rfold(k,nz,3) + zl + 0.5d0)*dr
              if(y.ge.cyledg1.and.y.le.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl)*dr
              y = (j + yl + 0.5d0)*dr
              z = (k + zl + 0.5d0)*dr
              if(y.gt.cyledg1.and.y.lt.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.lt.radsq) then
                mp(BX,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl)*dr
              y = (rfold(j,ny,2) + yl + 0.5d0)*dr
              z = (rfold(k,nz,3) + zl + 0.5d0)*dr
              if(y.gt.cyledg1.and.y.lt.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.lt.radsq) then
                mp(BX,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl + 0.5d0)*dr
              y = (j + yl)*dr
              z = (k + zl + 0.5d0)*dr
              if(y.ge.cyledg1.and.y.le.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.lt.radsq) then
                mp(BY,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl + 0.5d0)*dr
              y = (rfold(j,ny,2) + yl)*dr
              z = (rfold(k,nz,3) + zl + 0.5d0)*dr
              if(y.ge.cyledg1.and.y.le.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.lt.radsq) then
                mp(BY,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl + 0.5d0)*dr
              y = (j + yl + 0.5d0)*dr
              z = (k + zl)*dr
              if(y.gt.cyledg1.and.y.lt.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.lt.radsq) then
                mp(BZ,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl + 0.5d0)*dr
              y = (rfold(j,ny,2) + yl + 0.5d0)*dr
              z = (rfold(k,nz,3) + zl)*dr
              if(y.gt.cyledg1.and.y.lt.cyledg2.and. &
             &   (z-cylaxs1)**2+(x-cylaxs2)**2.lt.radsq) then
               mp(BZ,i,j,k,1) = 0.0d0
              end if
            end do
            end do
            end do
!GEOTYPE2-3
          else if(cylaln.eq.3) then
            do k=-1,zu-zl+1
            do j=-1,yu-zl+1
            do i=-1,xu-xl+1
              x = (i + xl + 0.5d0)*dr
              y = (j + yl)*dr
              z = (k + zl)*dr
              if(z.ge.cyledg1.and.z.le.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl + 0.5d0)*dr
              y = (rfold(j,ny,2) + yl)*dr
              z = (rfold(k,nz,3) + zl)*dr
              if(z.ge.cyledg1.and.z.le.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl)*dr
              y = (j + yl + 0.5d0)*dr
              z = (k + zl)*dr
              if(z.ge.cyledg1.and.z.le.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl)*dr
              y = (rfold(j,ny,2) + yl + 0.5d0)*dr
              z = (rfold(k,nz,3) + zl)*dr
              if(z.ge.cyledg1.and.z.le.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl)*dr
              y = (j + yl)*dr
              z = (k + zl + 0.5d0)*dr
              if(z.gt.cyledg1.and.z.lt.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl)*dr
              y = (rfold(j,ny,2) + yl)*dr
              z = (rfold(k,nz,3) + zl + 0.5d0)*dr
              if(z.gt.cyledg1.and.z.lt.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl)*dr
              y = (j + yl + 0.5d0)*dr
              z = (k + zl + 0.5d0)*dr
              if(z.gt.cyledg1.and.z.lt.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.lt.radsq) then
                mp(BX,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl)*dr
              y = (rfold(j,ny,2) + yl + 0.5d0)*dr
              z = (rfold(k,nz,3) + zl + 0.5d0)*dr
              if(z.gt.cyledg1.and.z.lt.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.lt.radsq) then
                mp(BX,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl + 0.5d0)*dr
              y = (j + yl)*dr
              z = (k + zl + 0.5d0)*dr
              if(z.gt.cyledg1.and.z.lt.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.lt.radsq) then
                mp(BY,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl + 0.5d0)*dr
              y = (rfold(j,ny,2) + yl)*dr
              z = (rfold(k,nz,3) + zl + 0.5d0)*dr
              if(z.gt.cyledg1.and.z.lt.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.lt.radsq) then
                mp(BY,i,j,k,1) = 0.0d0
              end if
!
              x = (i + xl + 0.5d0)*dr
              y = (j + yl + 0.5d0)*dr
              z = (k + zl)*dr
              if(z.ge.cyledg1.and.z.le.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.lt.radsq) then
                mp(BZ,i,j,k,1) = 0.0d0
              end if
              x = (rfold(i,nx,1) + xl + 0.5d0)*dr
              y = (rfold(j,ny,2) + yl + 0.5d0)*dr
              z = (rfold(k,nz,3) + zl)*dr
              if(z.ge.cyledg1.and.z.le.cyledg2.and. &
             &   (x-cylaxs1)**2+(y-cylaxs2)**2.lt.radsq) then
               mp(BZ,i,j,k,1) = 0.0d0
              end if
            end do
            end do
            end do
          end if
!
!GEOTYPE3
        else if(geotype(ipc).eq.3) then
          radsq = sphrad**2
          radsqin = (sphrad - 1.0d0)**2
!
          do k=-1,zu-zl+1
          do j=-1,yu-zl+1
          do i=-1,xu-xl+1
            x = (i + xl + 0.5d0)*dr
            y = (j + yl)*dr
            z = (k + zl)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.le.radsq) then
              mp(EX,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl + 0.5d0)*dr
            y = (rfold(j,ny,2) + yl)*dr
            z = (rfold(k,nz,3) + zl)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.le.radsq) then
              mp(EX,i,j,k,1) = 0.0d0
            end if
!
            x = (i + xl)*dr
            y = (j + yl + 0.5d0)*dr
            z = (k + zl)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.le.radsq) then
              mp(EY,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl)*dr
            y = (rfold(j,ny,2) + yl + 0.5d0)*dr
            z = (rfold(k,nz,3) + zl)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.le.radsq) then
              mp(EY,i,j,k,1) = 0.0d0
            end if
!
            x = (i + xl)*dr
            y = (j + yl)*dr
            z = (k + zl + 0.5d0)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.le.radsq) then
              mp(EZ,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl)*dr
            y = (rfold(j,ny,2) + yl)*dr
            z = (rfold(k,nz,3) + zl + 0.5d0)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.le.radsq) then
              mp(EZ,i,j,k,1) = 0.0d0
            end if
!
            x = (i + xl)*dr
            y = (j + yl + 0.5d0)*dr
            z = (k + zl + 0.5d0)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.lt.radsq) then
              mp(BX,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl)*dr
            y = (rfold(j,ny,2) + yl + 0.5d0)*dr
            z = (rfold(k,nz,3) + zl + 0.5d0)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.lt.radsq) then
              mp(BX,i,j,k,1) = 0.0d0
            end if
!
            x = (i + xl + 0.5d0)*dr
            y = (j + yl)*dr
            z = (k + zl + 0.5d0)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.lt.radsq) then
              mp(BY,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl + 0.5d0)*dr
            y = (rfold(j,ny,2) + yl)*dr
            z = (rfold(k,nz,3) + zl + 0.5d0)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.lt.radsq) then
              mp(BY,i,j,k,1) = 0.0d0
            end if
!
            x = (i + xl + 0.5d0)*dr
            y = (j + yl + 0.5d0)*dr
            z = (k + zl)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.lt.radsq) then
              mp(BZ,i,j,k,1) = 0.0d0
            end if
            x = (rfold(i,nx,1) + xl + 0.5d0)*dr
            y = (rfold(j,ny,2) + yl + 0.5d0)*dr
            z = (rfold(k,nz,3) + zl)*dr
            if((x-sphcnt1)**2+(y-sphcnt2)**2+(z-sphcnt3)**2.lt.radsq) then
              mp(BZ,i,j,k,1) = 0.0d0
            end if
          end do
          end do
          end do
!
        end if
      end do


    return
  end subroutine medium

