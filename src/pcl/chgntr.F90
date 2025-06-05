#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine chgntr
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   C H G N T R
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine gives charge density so as to keep      .
!   .  charge neutrality                                       .
!   .                                                          .
!   .  option for background fixed ions                        .
!   .   ionchg = 0 : uniform neutralization in the system      .
!   .          = 1 : local newtralization at each grid point   .
!   ............................................................
!
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
#define MCW local_comm
#define MSS MPI_STATUS_SIZE
  implicit none
!
  integer(kind=4) :: i,j,k, ii,jj,kk, i1,j1,k1
  integer(kind=4) :: is
  integer(kind=4) :: ipc
  integer(kind=4) :: itch
  integer(kind=4) :: ibdy,vol
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  real(kind=8) :: x,y,z, xxc,yyc,zzc
  real(kind=8) :: id0s, d1s
  real(kind=8) :: radsq(inpc)
  real(kind=8) :: rion, qele
  real(kind=8) :: xchg,ychg,zchg
  real(kind=8) :: disp1,disp2,disp3, displ

  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8

!-------------------- 
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)


!--------------------
      do ipc=1,npc
        if(geotype(ipc).eq.2) then
          radsq(ipc) = cylinder(ipc)%radius*cylinder(ipc)%radius
        else if(geotype(ipc).eq.3) then
          radsq(ipc) = sphere(ipc)%radius*sphere(ipc)%radius
        end if
      end do


    if(jobnum(1).eq.0.and.nflag_testp.ne.1) then
!-------------------- uniform neutralization in the system
      if(ionchg.eq.0) then
        if(myid.eq.0) then
          print*,'kempo: chgntr: uniform neutralization in system.'
        end if
!
        rion = 0.0d0
        do is=1,nspec
          if(npin(is).ge.0) then
            rion = rion - q(is)*npin(is)
          end if
        end do
!
        rion = rion/(nx*ny*nz)*1.0d-3
!
        if(abs(rion).lt.1.0d-10) rion = 0.0d0
!
        if(myid.eq.0) print*,'rion=',rion/renrho
!
!        do k=zl,zu
!        do j=yl,yu
!        do i=xl,xu
!          vol = 0
!          do ipc=1,npc
!            if(k.gt.nzpc1(ipc).and.k.lt.nzpc2(ipc).and. &
!           &   j.gt.nypc1(ipc).and.j.lt.nypc2(ipc).and. &
!           &   i.gt.nxpc1(ipc).and.i.lt.nxpc2(ipc)) then
!              vol = 8
!            end if
!          end do
!          if(vol.eq.0) then
!            do ibdy=1,nbdsf
!              if(i.ge.bdyvoxel(1,ibdy).and.i.le.bdyvoxel(1,ibdy)+1.and. &
!             &   j.ge.bdyvoxel(2,ibdy).and.j.le.bdyvoxel(2,ibdy)+1.and. &
!             &   k.ge.bdyvoxel(3,ibdy).and.k.le.bdyvoxel(3,ibdy)+1) then
!                vol = vol + 1
!              end if
!            end do
!          end if
!          rhobk(1,i-xl,j-yl,k-zl,3) = rion*(8-vol)*0.125d0
!        end do
!        end do
!        end do
        do k=-1,(zu-zl)+1
        do j=-1,(yu-yl)+1
        do i=-1,(xu-xl)+1
          if(i.ne.(xu-xl)+1.and.j.ne.(yu-yl)+1.and.k.ne.(zu-zl)+1) then
          do kk=0,9
          do jj=0,9
    iilp: do ii=0,9
            qele = rion
!
            xchg = i + ii*0.1d0 + xl
            ychg = j + jj*0.1d0 + yl
            zchg = k + kk*0.1d0 + zl
!
            do ipc=1,npc
              if((geotype(ipc).eq.0.or.geotype(ipc).eq.1).and. &
             &   (xchg.ge.xlpc(ipc).and.xchg.le.xupc(ipc).and. &
             &    ychg.ge.ylpc(ipc).and.ychg.le.yupc(ipc).and. &
             &    zchg.ge.zlpc(ipc).and.zchg.le.zupc(ipc))) then
                cycle iilp
              else if(geotype(ipc).eq.2) then
                if(cylinder(ipc)%align.eq.1) then
                  disp1 = ychg - cylinder(ipc)%axis(1)
                  disp2 = zchg - cylinder(ipc)%axis(2)
                  if(xchg.ge.cylinder(ipc)%edge(1).and.xchg.le.cylinder(ipc)%edge(2).and. &
                 &   disp1*disp1+disp2*disp2.le.radsq(ipc)) then
                    cycle iilp
                  end if
                else if(cylinder(ipc)%align.eq.2) then
                  disp1 = zchg - cylinder(ipc)%axis(1)
                  disp2 = xchg - cylinder(ipc)%axis(2)
                  if(ychg.ge.cylinder(ipc)%edge(1).and.ychg.le.cylinder(ipc)%edge(2).and. &
                 &   disp1*disp1+disp2*disp2.le.radsq(ipc)) then
                    cycle iilp
                  end if
                else if(cylinder(ipc)%align.eq.3) then
                  disp1 = xchg - cylinder(ipc)%axis(1)
                  disp2 = ychg - cylinder(ipc)%axis(2)
                  if(zchg.ge.cylinder(ipc)%edge(1).and.zchg.le.cylinder(ipc)%edge(2).and. &
                 &   disp1*disp1+disp2*disp2.le.radsq(ipc)) then
                    cycle iilp
                  end if
                end if
              else if(geotype(ipc).eq.3) then
                disp1 = xchg - sphere(ipc)%center(1)
                disp2 = ychg - sphere(ipc)%center(2)
                disp3 = zchg - sphere(ipc)%center(3)
                if(disp1*disp1+disp2*disp2+disp3*disp3.le.radsq(ipc)) then
                  cycle iilp
                end if
              end if
            end do
!
            i1 = i + 1
            j1 = j + 1
            k1 = k + 1
!
            x1 = ii*0.1d0
            y1 = jj*0.1d0
            z1 = kk*0.1d0
!
            xy1 = x1 * y1
            xz1 = x1 * z1
            yz1 = y1 * z1
            z2 = 1.0d0 - z1
            xz2 = x1 * z2
            yz2 = y1 * z2
!
            v3 = xy1 * z1
            v2 = xz1 - v3
            v4 = yz1 - v3
            v1 = z1 - xz1 - v4
            v7 = xy1 * z2
            v6 = xz2 - v7
            v8 = yz2 - v7
            v5 = z2 - xz2 - v8
!
            rhobk(1,i ,j ,k1,3) = rhobk(1,i ,j ,k1,3) + v1*qele
            rhobk(1,i1,j ,k1,3) = rhobk(1,i1,j ,k1,3) + v2*qele
            rhobk(1,i1,j1,k1,3) = rhobk(1,i1,j1,k1,3) + v3*qele
            rhobk(1,i ,j1,k1,3) = rhobk(1,i ,j1,k1,3) + v4*qele
            rhobk(1,i ,j ,k ,3) = rhobk(1,i ,j ,k ,3) + v5*qele
            rhobk(1,i1,j ,k ,3) = rhobk(1,i1,j ,k ,3) + v6*qele
            rhobk(1,i1,j1,k ,3) = rhobk(1,i1,j1,k ,3) + v7*qele
            rhobk(1,i ,j1,k ,3) = rhobk(1,i ,j1,k ,3) + v8*qele
          end do iilp
          end do
          end do
          end if
          if(i.lt.0.or.i.gt.(xu-xl).or. &
             j.lt.0.or.j.gt.(yu-yl).or. &
             k.lt.0.or.k.gt.(zu-zl)) then
            rhobk(1,i,j,k,3) = 0.0d0
          end if
        end do
        end do
        end do



!-------------------- local newtralization at each grid point
      else
        if(myid.eq.0) then
          print*,'kempo: chgntr: local newtralization at each grid point.'
        end if
!       ---------------
        call charge(1,0)
        if(sdid(2).ge.0) call charge(2,0)
        if(currmode.ne.0) &
       &  call oh3_reduce_field(rho(1,0,0,0,1),rho(1,0,0,0,2),FRH)
        call oh3_exchange_borders(rho(1,0,0,0,1),rho(1,0,0,0,2),CRH,0)
        call add_boundary_charge(rho(:,:,:,:,1), &
       &                         1,1,sdoms(:,:,sdid(1)+1), &
       &                         bounds(:,:,sdid(1)+1),ctypes(:,:,1,CRH), &
       &                         fsizes(:,:,FRH),myid)
!       ---------------
        do k=0,zu-zl
        do j=0,yu-yl
        do i=0,xu-xl
          rhobk(1,i,j,k,3) = -rho(1,i,j,k,1)
        end do
        end do
        end do
!
      end if
!
!
    end if


      do itch=1,nqclst
        xxc = rclst(1,itch); yyc = rclst(2,itch); zzc = rclst(3,itch)
        id0s = 0.5d0/(rclst(4,itch)*rclst(4,itch))
        d1s = rclst(5,itch)*rclst(5,itch)
        if(myid.eq.0) print*, "chgntr: ", rclst(:,itch),rhopeak(itch),xxc,yyc,zzc,id0s,d1s
        do k=0,zu-zl
        do j=0,yu-yl
        do i=0,xu-xl
          x = i + xl; y = j + yl; z = k + zl
          displ = dimclst(1,itch)*(x - xxc)*(x - xxc) &
         &      + dimclst(2,itch)*(y - yyc)*(y - yyc) &
         &      + dimclst(3,itch)*(z - zzc)*(z - zzc)
          if(displ.lt.d1s) then
            rhobk(1,i,j,k,3) = rhobk(1,i,j,k,3) &
           &                 + rhopeak(itch)*exp(-displ*id0s)
          end if
        end do
        end do
        end do
      end do


  return
  end subroutine
