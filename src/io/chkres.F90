#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine chkres(ps,level,fcomp,message)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   P O S I T N
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .     this subroutine updates particle position.           .
!   .     the leap-frog scheme is used.                        .
!   ............................................................
!

!-------------------- parameters and variables
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=8) :: m, ns,ne
  integer(kind=4) :: is
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ps, level, fcomp, rid
  character(len=10) :: message


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!============================== species loop
      ne = pbase(ps)
 ISL: do is=1,nspec
!-------------------- inner loop
        ns = ne + 1
        ne = ne + totalp(is,ps)
        do m=ns,ne
          if(pbuf(m)%nid.eq.-1) cycle
!
          if(level.eq.0) then
            if(pbuf(m)%x.lt.xl.or.pbuf(m)%x.ge.xu.or. &
           &   pbuf(m)%y.lt.yl.or.pbuf(m)%y.ge.yu.or. &
           &   pbuf(m)%z.lt.zl.or.pbuf(m)%z.ge.zu) then
              print*, message, "A:myid,is,m,x,y,z,nid =", myid,is,m, &
             &  pbuf(m)%x,pbuf(m)%y,pbuf(m)%z,pbuf(m)%nid
            end if
            rid = oh3_map_particle_to_subdomain &
           &        (pbuf(m)%x,pbuf(m)%y,pbuf(m)%z)
            if(rid.ne.pbuf(m)%nid) then
              print*, message, "B:myid,is,m,x,y,z,nid,rid =", myid,is,m, &
             &  pbuf(m)%x,pbuf(m)%y,pbuf(m)%z,pbuf(m)%nid,rid
            end if
          else if(level.eq.1) then
            if(floor(pbuf(m)%x-xl).lt.fsizes(1,1,fcomp).or. &
           &   floor(pbuf(m)%x-xl).gt.fsizes(2,1,fcomp)-1.or. &
           &   floor(pbuf(m)%y-yl).lt.fsizes(1,2,fcomp).or. &
           &   floor(pbuf(m)%y-yl).gt.fsizes(2,2,fcomp)-1.or. &
           &   floor(pbuf(m)%z-zl).lt.fsizes(1,3,fcomp).or. &
           &   floor(pbuf(m)%z-zl).gt.fsizes(2,3,fcomp)-1) then
              print*, message, "myid,is,m,x,y,z,nid =", myid,is,m, &
             &  pbuf(m)%x,pbuf(m)%y,pbuf(m)%z,pbuf(m)%nid
            end if
          else if(level.eq.2) then
            if(pbuf(m)%x.eq.128.0d0.or. &
           &   pbuf(m)%y.eq.128.0d0.or. &
           &   pbuf(m)%z.eq.128.0d0) then
              print*, message, "myid,is,m,x,y,z,nid =", myid,is,m, &
             &  pbuf(m)%x,pbuf(m)%y,pbuf(m)%z,pbuf(m)%nid
            end if
          end if
        end do
!
      end do ISL
!============================== species loop end


  return
  end subroutine
