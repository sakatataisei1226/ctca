#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine energy
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   E N E R G Y
!   ____________________________________________________________
!
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
#define MCW local_comm
  implicit none
!
  integer(kind=8) :: m, ns,ne, nee
  integer(kind=4) :: i,j,k
  integer(kind=4) :: is
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ps
  integer(kind=4) :: ierr
  real(kind=8) :: tkx,tky,tkz, te,tb


!-------------------- 
      engk(1:3,1:nspec) = 0.0d0
      do ps=1,2
        if(sdid(ps).lt.0) cycle
        nee = pbase(ps)
        do is=1,nspec
          ns = nee + 1
          ne = nee + totalp(is,ps)
          nee = ne
          tkx = 0.0d0; tky = 0.0d0; tkz = 0.0d0
          do m=ns,ne
            if(pbuf(m)%nid.eq.-1) cycle
            tkx = tkx + pbuf(m)%vx*pbuf(m)%vx
            tky = tky + pbuf(m)%vy*pbuf(m)%vy
            tkz = tkz + pbuf(m)%vz*pbuf(m)%vz
          end do
          engk(1,is) = engk(1,is) + 0.5d0*rm(is)*tkx/(slx*sly*slz)
          engk(2,is) = engk(2,is) + 0.5d0*rm(is)*tky/(slx*sly*slz)
          engk(3,is) = engk(3,is) + 0.5d0*rm(is)*tkz/(slx*sly*slz)
        end do
      end do


!-------------------- 
      te = 0.0d0; tb = 0.0d0
      do ps=1,1
        xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
        yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
        zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)
        do k=0,zu-zl-1
        do j=0,yu-yl-1
        do i=0,xu-xl-1
          te = te + eb(EX,i,j,k,ps)*eb(EX,i,j,k,ps) &
         &        + eb(EY,i,j,k,ps)*eb(EY,i,j,k,ps) &
         &        + eb(EZ,i,j,k,ps)*eb(EZ,i,j,k,ps)
          tb = tb + eb(BX,i,j,k,ps)*eb(BX,i,j,k,ps) &
         &        + eb(BY,i,j,k,ps)*eb(BY,i,j,k,ps) &
         &        + eb(BZ,i,j,k,ps)*eb(BZ,i,j,k,ps)
        end do
        end do
        end do
      end do
      engeb(1) = 0.5d0*te/(slx*sly*slz)
      engeb(2) = 0.5d0*cs*tb/(slx*sly*slz)


!-------------------- 
      call MPI_Reduce(engk,engkg,3*nspec,MPI_REAL8,MPI_SUM,0,MCW,ierr)
      call MPI_Reduce(engeb,engebg,2,MPI_REAL8,MPI_SUM,0,MCW,ierr)


  return
  end subroutine
