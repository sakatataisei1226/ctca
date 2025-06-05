#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine reflect_field(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!                    R E F L E C T _ F I E L D
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
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ps


!-------------------- test particle simulation
      if(juncan.ge.1000) return


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- 
      if(bounds(1,3,sdid(ps)+1).ne.1) then
        do j=-1,yu-yl+1
        do i=-1,xu-xl+1
          if(-1.ge.zl-1.and.-1.le.zu+1) then
            eb(EX,i,j,-1  -zl,ps) = -eb(EX,i,j,+1  -zl,ps)
            eb(EY,i,j,-1  -zl,ps) = -eb(EY,i,j,+1  -zl,ps)
            eb(BZ,i,j,-1  -zl,ps) = -eb(BZ,i,j,+1  -zl,ps)
!
            eb(EZ,i,j,-1  -zl,ps) = +eb(EZ,i,j, 0  -zl,ps)
            eb(BX,i,j,-1  -zl,ps) = +eb(BX,i,j, 0  -zl,ps)
            eb(BY,i,j,-1  -zl,ps) = +eb(BY,i,j, 0  -zl,ps)
          end if
!
          if( 0.ge.zl-1.and. 0.le.zu+1) then
            eb(EX,i,j, 0  -zl,ps) = 0.0d0
            eb(EY,i,j, 0  -zl,ps) = 0.0d0
            eb(BZ,i,j, 0  -zl,ps) = 0.0d0
          end if
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,3,sdid(ps)+1).ne.1) then
        do j=-1,yu-yl+1
        do i=-1,xu-xl+1
          if(nz.ge.zl-1.and.nz.le.zu+1) then
            eb(EX,i,j,nz  -zl,ps) = 0.0d0
            eb(EY,i,j,nz  -zl,ps) = 0.0d0
            eb(BZ,i,j,nz  -zl,ps) = 0.0d0
!
            eb(EZ,i,j,nz  -zl,ps) = +eb(EZ,i,j,nz-1-zl,ps)
            eb(BX,i,j,nz  -zl,ps) = +eb(BX,i,j,nz-1-zl,ps)
            eb(BY,i,j,nz  -zl,ps) = +eb(BY,i,j,nz-1-zl,ps)
          end if
!
          if(nz+1.ge.zl-1.and.nz+1.le.zu+1) then
            eb(EX,i,j,nz+1-zl,ps) = -eb(EX,i,j,nz-1-zl,ps)
            eb(EY,i,j,nz+1-zl,ps) = -eb(EY,i,j,nz-1-zl,ps)
            eb(BZ,i,j,nz+1-zl,ps) = -eb(BZ,i,j,nz-1-zl,ps)
!
            eb(EZ,i,j,nz+1-zl,ps) = +eb(EZ,i,j,nz-2-zl,ps)
            eb(BX,i,j,nz+1-zl,ps) = +eb(BX,i,j,nz-2-zl,ps)
            eb(BY,i,j,nz+1-zl,ps) = +eb(BY,i,j,nz-2-zl,ps)
          end if
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,2,sdid(ps)+1).ne.1) then
        do k=-1,zu-zl+1
        do i=-1,xu-xl+1
          if(-1.ge.yl-1.and.-1.le.yu+1) then
            eb(EX,i,-1  -yl,k,ps) = -eb(EX,i,+1  -yl,k,ps)
            eb(EZ,i,-1  -yl,k,ps) = -eb(EZ,i,+1  -yl,k,ps)
            eb(BY,i,-1  -yl,k,ps) = -eb(BY,i,+1  -yl,k,ps)
!
            eb(EY,i,-1  -yl,k,ps) = +eb(EY,i, 0  -yl,k,ps)
            eb(BX,i,-1  -yl,k,ps) = +eb(BX,i, 0  -yl,k,ps)
            eb(BZ,i,-1  -yl,k,ps) = +eb(BZ,i, 0  -yl,k,ps)
          end if
!
          if( 0.ge.yl-1.and. 0.le.yu+1) then
            eb(EX,i, 0  -yl,k,ps) = 0.0d0
            eb(EZ,i, 0  -yl,k,ps) = 0.0d0
            eb(BY,i, 0  -yl,k,ps) = 0.0d0
          end if
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,2,sdid(ps)+1).ne.1) then
        do k=-1,zu-zl+1
        do i=-1,xu-xl+1
          if(ny.ge.yl-1.and.ny.le.yu+1) then
            eb(EX,i,ny  -yl,k,ps) = 0.0d0
            eb(EZ,i,ny  -yl,k,ps) = 0.0d0
            eb(BY,i,ny  -yl,k,ps) = 0.0d0
!
            eb(EY,i,ny  -yl,k,ps) = +eb(EY,i,ny-1-yl,k,ps)
            eb(BX,i,ny  -yl,k,ps) = +eb(BX,i,ny-1-yl,k,ps)
            eb(BZ,i,ny  -yl,k,ps) = +eb(BZ,i,ny-1-yl,k,ps)
          end if
!
          if(ny+1.ge.yl-1.and.ny+1.le.yu+1) then
            eb(EX,i,ny+1-yl,k,ps) = -eb(EX,i,ny-1-yl,k,ps)
            eb(EZ,i,ny+1-yl,k,ps) = -eb(EZ,i,ny-1-yl,k,ps)
            eb(BY,i,ny+1-yl,k,ps) = -eb(BY,i,ny-1-yl,k,ps)
!
            eb(EY,i,ny+1-yl,k,ps) = +eb(EY,i,ny-2-yl,k,ps)
            eb(BX,i,ny+1-yl,k,ps) = +eb(BX,i,ny-2-yl,k,ps)
            eb(BZ,i,ny+1-yl,k,ps) = +eb(BZ,i,ny-2-yl,k,ps)
          end if
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,1,sdid(ps)+1).ne.1) then
        do k=-1,zu-zl+1
        do j=-1,yu-yl+1
          if(-1.ge.xl-1.and.-1.le.xu+1) then
            eb(EY,-1  -xl,j,k,ps) = -eb(EY,+1  -xl,j,k,ps)
            eb(EZ,-1  -xl,j,k,ps) = -eb(EZ,+1  -xl,j,k,ps)
            eb(BX,-1  -xl,j,k,ps) = -eb(BX,+1  -xl,j,k,ps)
!
            eb(EX,-1  -xl,j,k,ps) = +eb(EX, 0  -xl,j,k,ps)
            eb(BY,-1  -xl,j,k,ps) = +eb(BY, 0  -xl,j,k,ps)
            eb(BZ,-1  -xl,j,k,ps) = +eb(BZ, 0  -xl,j,k,ps)
          end if
!
          if( 0.ge.xl-1.and. 0.le.xu+1) then
            eb(EY, 0  -xl,j,k,ps) = 0.0d0
            eb(EZ, 0  -xl,j,k,ps) = 0.0d0
            eb(BX, 0  -xl,j,k,ps) = 0.0d0
          end if
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,1,sdid(ps)+1).ne.1) then
        do k=-1,zu-zl+1
        do j=-1,yu-yl+1
          if(nx.ge.xl-1.and.nx.le.xu+1) then
            eb(EY,nx  -xl,j,k,ps) = 0.0d0
            eb(EZ,nx  -xl,j,k,ps) = 0.0d0
            eb(BX,nx  -xl,j,k,ps) = 0.0d0
!
            eb(EX,nx  -xl,j,k,ps) = +eb(EX,nx-1-xl,j,k,ps)
            eb(BY,nx  -xl,j,k,ps) = +eb(BY,nx-1-xl,j,k,ps)
            eb(BZ,nx  -xl,j,k,ps) = +eb(BZ,nx-1-xl,j,k,ps)
          end if
!
          if(nx+1.ge.xl-1.and.nx+1.le.xu+1) then
            eb(EY,nx+1-xl,j,k,ps) = -eb(EY,nx-1-xl,j,k,ps)
            eb(EZ,nx+1-xl,j,k,ps) = -eb(EZ,nx-1-xl,j,k,ps)
            eb(BX,nx+1-xl,j,k,ps) = -eb(BX,nx-1-xl,j,k,ps)
!
            eb(EX,nx+1-xl,j,k,ps) = +eb(EX,nx-2-xl,j,k,ps)
            eb(BY,nx+1-xl,j,k,ps) = +eb(BY,nx-2-xl,j,k,ps)
            eb(BZ,nx+1-xl,j,k,ps) = +eb(BZ,nx-2-xl,j,k,ps)
          end if
        end do
        end do
      end if


  return
  end subroutine



!
  subroutine replace_boundary_field(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!              A D D _ B O U N D A R Y _ C U R R E N T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a charge distribution to      .
!   .      grid points from particle locations by the          .
!   .      "area-sharing" or "area-weighting" scheme.          .
!   ............................................................
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
  integer(kind=4) :: sl,su, dl,du, nl,nu
  integer(kind=4) :: ps


!-------------------- 
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- 
      sl = ctypes(2,2,1,CAJ)	!=-2
      su = ctypes(2,1,1,CAJ)	!=+2
!
      nl = ctypes(3,2,1,CAJ)	!= 1
      nu = ctypes(3,1,1,CAJ)	!= 2
!
      dl = sl + nl		!=-1
      du = su - nu		!= 0

!-------------------- 
      if(bounds(1,3,sdid(ps)+1).eq.1) then
        do k=0,nl-1		!(i.e., do k=0,0)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+5)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+5)
          eb(EX,sl+i,sl+j,dl+k,ps) = eb(EX,sl+i,sl+j,sl+k,ps)
          eb(EY,sl+i,sl+j,dl+k,ps) = eb(EY,sl+i,sl+j,sl+k,ps)
          eb(EZ,sl+i,sl+j,dl+k,ps) = eb(EZ,sl+i,sl+j,sl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,3,sdid(ps)+1).eq.1) then
        do k=0,nu-1		!(i.e., do k=0,1)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+5)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+5)
          eb(EX,sl+i,sl+j,zu+du+k,ps) = eb(EX,sl+i,sl+j,zu+su+k,ps)
          eb(EY,sl+i,sl+j,zu+du+k,ps) = eb(EY,sl+i,sl+j,zu+su+k,ps)
          eb(EZ,sl+i,sl+j,zu+du+k,ps) = eb(EZ,sl+i,sl+j,zu+su+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,2,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nl-1		!(i.e., do j=0,0)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+5)
          eb(EX,sl+i,dl+j,dl+k,ps) = eb(EX,sl+i,sl+j,dl+k,ps)
          eb(EY,sl+i,dl+j,dl+k,ps) = eb(EY,sl+i,sl+j,dl+k,ps)
          eb(EZ,sl+i,dl+j,dl+k,ps) = eb(EZ,sl+i,sl+j,dl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,2,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nu-1		!(i.e., do j=0,1)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+5)
          eb(EX,sl+i,yu+du+j,dl+k,ps) = eb(EX,sl+i,yu+su+j,dl+k,ps)
          eb(EY,sl+i,yu+du+j,dl+k,ps) = eb(EY,sl+i,yu+su+j,dl+k,ps)
          eb(EZ,sl+i,yu+du+j,dl+k,ps) = eb(EZ,sl+i,yu+su+j,dl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,1,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nl-1		!(i.e., do i=0,0)
          eb(EX,dl+i,dl+j,dl+k,ps) = eb(EX,sl+i,dl+j,dl+k,ps)
          eb(EY,dl+i,dl+j,dl+k,ps) = eb(EY,sl+i,dl+j,dl+k,ps)
          eb(EZ,dl+i,dl+j,dl+k,ps) = eb(EZ,sl+i,dl+j,dl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,1,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nu-1		!(i.e., do i=0,1)
          eb(EX,xu+du+i,dl+j,dl+k,ps) = eb(EX,xu+su+i,dl+j,dl+k,ps)
          eb(EY,xu+du+i,dl+j,dl+k,ps) = eb(EY,xu+su+i,dl+j,dl+k,ps)
          eb(EZ,xu+du+i,dl+j,dl+k,ps) = eb(EZ,xu+su+i,dl+j,dl+k,ps)
        end do
        end do
        end do
      end if


  return
  end subroutine
