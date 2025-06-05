#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine surchg
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   S U R C H G
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .    calculate surface charge on internal boundaries by    .
!   .    3-dimemsional poisson equation solver included        .
!   .    capcity matrix method.                                .
!   ............................................................
!
!
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
#define MCW local_comm
#define MSS MPI_STATUS_SIZE
#define MIP MPI_IN_PLACE
  implicit none
!
  integer(kind=4) :: i,j,k, i1,j1,k1
  integer(kind=4) :: is
  integer(kind=4) :: ipc,jpc, ipcg,jpcg, icap,jcap, mcap
  integer(kind=4) :: ipp(ixw)
  integer(kind=4) :: icon
  integer(kind=4) :: src1, src2
  integer(kind=4) :: icomm, ierr
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: siu,sju,sku
  integer(kind=4) :: swstep
  real(kind=8) :: ipswper
  real(kind=8) :: amx(inpc,inpc)
  real(kind=8) :: wok(ixw)
  real(kind=8) :: wrelaxinst
  real(kind=8) :: xlocal,ylocal,zlocal
  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8
  real(kind=8) :: sendw(2,ixw),recvw(2,ixw)
  real(kind=8) :: vfactor, disp, r1, r2
  real(kind=8) :: tew


!-------------------- 
      if(abs(modeww).lt.2) return


!-------------------- 
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)
      if(pftmode.eq.3.and.myid.ne.snode-1) then
        siu = sxu - sxl - 1
      else
        siu = sxu - sxl
      end if
      if(pftmode.eq.2.and.myid.ne.snode-1) then
        sju = syu - syl - 1
      else
        sju = syu - syl
      end if
      if((pftmode.eq.1.or.pftmode.eq.4).and.myid.ne.snode-1) then
        sku = szu - szl - 1
      else
        sku = szu - szl
      end if


!-------------------- zero clear
      phi(:,:,:,:,:) = 0.0d0
      poi(:,:,:,:,:) = 0.0d0


!-------------------- 
      if(pswper.ge.4) then
        swstep = istep - pswstr
        ipswper = 1.0d0/pswper
        do ipcg=1,npcg
          do ipc=pcgs(ipcg)+1,pcgs(ipcg+1)
            if(ipc.eq.pcgs(ipcg)+1.and.abs(mtd_vchg(ipcg)).eq.1.and. &
           &   pswspn(ipc).lt.0.0d0) then
              if(swstep.ge.0) then
                pfixed(ipcg) = sign(pswspn(ipc)*(1.0d0-abs(cos(pi2*(swstep*ipswper+pswini)))), &
               &                    -sin(pi2*(swstep*ipswper+pswini)))
!                pfixed(ipcg) = pswspn(ipc)*sin(pi2*(swstep*ipswper+pswini))
              else
                pfixed(ipcg) = 0.0d0
              end if
            else if(ipc.ne.pcgs(ipcg)+1.and.pswspn(ipc).lt.0.0d0) then
              if(swstep.ge.0) then
                biasp(ipc) = sign(pswspn(ipc)*(1.0d0-abs(cos(pi2*(swstep*ipswper+pswini)))), &
               &                    -sin(pi2*(swstep*ipswper+pswini)))
!                biasp(ipc) = pswspn(ipc)*sin(pi2*(swstep*ipswper+pswini))
              else
                biasp(ipc) = 0.0d0
              end if
            end if
          end do
        end do
      end if


!-------------------- Poisson Solver
      call poisson(1)


      if(abs(modeww).ge.2) then
!---- difference of solution and reference 
!-------------------- get potential on bodies (not modified yet)
!        sfrho(:) = 0.0d0
!        mcap = 0
        diff(:) = 0.0d0
        do icap=1,ncpmx
          xlocal = bdygrid(1,icap) - sxl
          ylocal = bdygrid(2,icap) - syl
          zlocal = bdygrid(3,icap) - szl
!
          i = floor(xlocal)
          j = floor(ylocal)
          k = floor(zlocal)
!
          i1 = i + 1
          j1 = j + 1
          k1 = k + 1
!
          x1 = xlocal - i
          y1 = ylocal - j
          z1 = zlocal - k
!
!          i = floor(bdygrid(1,icap))
!          j = floor(bdygrid(2,icap)) - syl
!          k = floor(bdygrid(3,icap)) - szl
!
!          i1 = i + 1
!          j1 = j + 1
!          k1 = k + 1
!
!          x1 = bdygrid(1,icap) - floor(bdygrid(1,icap))
!          y1 = bdygrid(2,icap) - floor(bdygrid(2,icap))
!          z1 = bdygrid(3,icap) - floor(bdygrid(3,icap))
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
!
          v7 = xy1 * z2
          v6 = xz2 - v7
          v8 = yz2 - v7
          v5 = z2 - xz2 - v8
!
          if(i.ge.0.and.i.le.siu.and. &
         &   j.ge.0.and.j.le.sju.and. &
         &   k.ge.0.and.k.le.sku) then
!            mcap = mcap + 1
!            sfrho(mcap) = poi(1,i ,j ,k1,1)*v1 &
!           &            + poi(1,i1,j ,k1,1)*v2 &
!           &            + poi(1,i1,j1,k1,1)*v3 &
!           &            + poi(1,i ,j1,k1,1)*v4 &
!           &            + poi(1,i ,j ,k ,1)*v5 &
!           &            + poi(1,i1,j ,k ,1)*v6 &
!           &            + poi(1,i1,j1,k ,1)*v7 &
!           &            + poi(1,i ,j1,k ,1)*v8
            diff(icap) = poi(1,i ,j ,k1,1)*v1 &
           &           + poi(1,i1,j ,k1,1)*v2 &
           &           + poi(1,i1,j1,k1,1)*v3 &
           &           + poi(1,i ,j1,k1,1)*v4 &
           &           + poi(1,i ,j ,k ,1)*v5 &
           &           + poi(1,i1,j ,k ,1)*v6 &
           &           + poi(1,i1,j1,k ,1)*v7 &
           &           + poi(1,i ,j1,k ,1)*v8
          end if
        end do
!        call MPI_Allgatherv(sfrho,mcap,MPI_REAL8, &
!       &                    diff,ncpmxs,scpmxs,MPI_REAL8,MCW,ierr)
        call MPI_Allreduce(MIP,diff,ncpmx,MPI_REAL8,MPI_SUM,MCW,ierr)
        if(ewmodel.eq.1.or.ewmodel.eq.3) then
          tew = t - dt*nretard
          diff(1:ncpmx) = diff(1:ncpmx) &
         &              - (bdygrid(1,1:ncpmx) - prefcrd(1))* &
         &                (e0x + Ew(1,0,1)*cos(omegaw(1)*tew)) &
         &              - (bdygrid(2,1:ncpmx) - prefcrd(2))* &
         &                (e0y + Ew(2,0,1)*sin(omegaw(1)*tew)) &
         &              - (bdygrid(3,1:ncpmx) - prefcrd(3))*e0z
        else
          diff(1:ncpmx) = diff(1:ncpmx) &
         &              - (bdygrid(1,1:ncpmx) - prefcrd(1))*e0x &
         &              - (bdygrid(2,1:ncpmx) - prefcrd(2))*e0y &
         &              - (bdygrid(3,1:ncpmx) - prefcrd(3))*e0z
        end if


!-------------------- sequential computation from here
!      if(myid.eq.cproot) then
!-------------------- get charge on bodies
!  (sfrho) = [cpmx](diff)
!
        sfrho(:) = 0.0d0
        if(abs(modeww).eq.2) then
          call dmav(cpmx,size(cpmx,1),ncpmx,ncpmx,diff,sfrho,icon)
        else if(abs(modeww).eq.3.and.myid.lt.snode) then
          call dmav(cpmx(nsdcpmx(myid+1):nedcpmx(myid+1),1:ncpmx), &
         &          size(cpmx,1),ndcpmx(myid+1),ncpmx, &
         &          diff,sfrho(nsdcpmx(myid+1):nedcpmx(myid+1)),icon)
          call MPI_Allreduce(MIP,sfrho,ncpmx,MPI_REAL8,MPI_SUM,subcomm,ierr)
        end if
        if(icon.ne.0) print*,'[surchg2].WARN01: icon@mav=',icon


!-------------------- sum (cpmx*phi)
!        if(t.lt.dt*512) then
!          wrelaxinst = 0.0d0
!        else if(t.le.dt*2048) then
!        if(t.le.dt*2048) then
!          wrelaxinst = (tanh(-4+(t-dt*0)/(dt*2048)*8)+1.0d0)*0.5d0*wrelax
!        else
          wrelaxinst = wrelax
!        end if

!       ------------- accumulated charge
        if(modeww.lt.0.and.(emflag.eq.0.or.istep.eq.0)) then
          do ipcg=1,npcg
            groupp(ipcg) = sum(gcount(2)%chgacm(2,ccgs(ipcg)+1:ccgs(ipcg+1)))*wrelaxinst
            groupp(ipcg) = groupp(ipcg) &
           &  + sum(sfrho(nscpmx(ccgs(ipcg)+1):necpmx(ccgs(ipcg+1))))
            do jpc=ccgs(ipcg)+1,ccgs(ipcg+1)
              do ipc=1,npc
                groupp(ipcg) = groupp(ipcg) - biasp(ipc)*bpmx(jpc,ipc)*dscaled(ipc)
              end do
            end do
            do jpcg=1,npcg
              if(mtd_vchg(jpcg).eq.1) then
                groupp(ipcg) = groupp(ipcg) - pfixed(jpcg)*apmx(ipcg,jpcg)
              end if
            end do
          end do
        else
          do ipcg=1,npcg
            groupp(ipcg) = &
           &  + sum(sfrho(nscpmx(ccgs(ipcg)+1):necpmx(ccgs(ipcg+1))))
            do jpc=ccgs(ipcg)+1,ccgs(ipcg+1)
              do ipc=1,npc
                groupp(ipcg) = groupp(ipcg) - biasp(ipc)*bpmx(jpc,ipc)*dscaled(ipc)
              end do
            end do
            do jpcg=1,npcg
              if(mtd_vchg(jpcg).eq.1) then
                groupp(ipcg) = groupp(ipcg) - pfixed(jpcg)*apmx(ipcg,jpcg)
              end if
            end do
          end do
        end if


!-------------------- corrected potential of each body
        do jpcg=1,npcg
          if(mtd_vchg(jpcg).ne.1) then
            do ipcg=1,npcg
              amx(ipcg,jpcg) = apmx(ipcg,jpcg)
            end do
          else
            do ipcg=1,npcg
              amx(ipcg,jpcg) = bcmx(ipcg,jpcg)
            end do
          end if
        end do
        if(npcg.eq.1) then
          groupp(1) = groupp(1)/amx(1,1)
        else
          call dlax(amx,inpc,npcg,groupp,0.0d0,1,is,wok,ipp,icon)
          if(icon.ne.0) print*,'[surchg2].WARN02: icon@lax=',icon
        end if


!-------------------- potential value expansion
        do ipcg=1,npcg
          if(abs(mtd_vchg(ipcg)).ne.1) then
            do ipc=pcgs(ipcg)+1,pcgs(ipcg+1)
              selfp(ipc) = groupp(ipcg) + biasp(ipc)
            end do
          else
            do ipc=pcgs(ipcg)+1,pcgs(ipcg+1)
              selfp(ipc) = pfixed(ipcg) + biasp(ipc)
            end do
            i = 1
            do while(biasc(i)%to.ne.0)
              if(abs(biasc(i)%to).ge.ccgs(ipcg)+1.and. &
           &     abs(biasc(i)%to).le.ccgs(ipcg+1)) then
                biasc(i)%val = groupp(ipcg) - groupph(ipcg)
                groupph(ipcg) = groupp(ipcg)
              end if
              i = i + 1
            end do
          end if
        end do


!-------------------- get potential modification
        do ipc=1,npc
          if(boom(ipc)%align.eq.1) then
            vfactor = 0.5d0/log(2.0d0*boom(ipc)%hlength/boom(ipc)%rradius)
            do icap=nscpmx(ipc),necpmx(ipc)
              disp = bdygrid(1,icap) - boom(ipc)%origin(1)
              r1 = sqrt((disp - boom(ipc)%hlength)*(disp - boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              r2 = sqrt((disp + boom(ipc)%hlength)*(disp + boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              diff(icap) = &
             &  selfp(ipc)*vfactor &
             &            *log((+boom(ipc)%hlength - disp + r1) &
             &                /(-boom(ipc)%hlength - disp + r2)) &
             &  - diff(icap)
            end do
          else if(boom(ipc)%align.eq.2) then
            vfactor = 0.5d0/log(2.0d0*boom(ipc)%hlength/boom(ipc)%rradius)
            do icap=nscpmx(ipc),necpmx(ipc)
              disp = bdygrid(2,icap) - boom(ipc)%origin(2)
              r1 = sqrt((disp - boom(ipc)%hlength)*(disp - boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              r2 = sqrt((disp + boom(ipc)%hlength)*(disp + boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              diff(icap) = &
             &  selfp(ipc)*vfactor &
             &            *log((+boom(ipc)%hlength - disp + r1) &
             &                /(-boom(ipc)%hlength - disp + r2)) &
             &  - diff(icap)
            end do
          else if(boom(ipc)%align.eq.3) then
            vfactor = 0.5d0/log(2.0d0*boom(ipc)%hlength/boom(ipc)%rradius)
            do icap=nscpmx(ipc),necpmx(ipc)
              disp = bdygrid(3,icap) - boom(ipc)%origin(3)
              r1 = sqrt((disp - boom(ipc)%hlength)*(disp - boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              r2 = sqrt((disp + boom(ipc)%hlength)*(disp + boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              diff(icap) = &
             &  selfp(ipc)*vfactor &
             &            *log((+boom(ipc)%hlength - disp + r1) &
             &                /(-boom(ipc)%hlength - disp + r2)) &
             &  - diff(icap)
            end do
          else
            do icap=nscpmx(ipc),necpmx(ipc)
              diff(icap) = selfp(ipc) - diff(icap)
            end do
          end if
        end do


!-------------------- surface charge modification
        sfrho(1:ncpmx) = 0.0d0
        if(abs(modeww).eq.2) then
          call dmav(cpmx,size(cpmx,1),ncpmx,ncpmx,diff,sfrho,icon)
        else if(abs(modeww).eq.3) then
          if(myid.lt.snode) then
            call dmav(cpmx(nsdcpmx(myid+1):nedcpmx(myid+1),1:ncpmx), &
           &          size(cpmx,1),ndcpmx(myid+1),ncpmx, &
           &          diff,sfrho(nsdcpmx(myid+1):nedcpmx(myid+1)),icon)
          end if
          call MPI_Allreduce(MIP,sfrho,ncpmx,MPI_REAL8,MPI_SUM,MCW,ierr)
        end if
        do ipc=1,npc
          sfrho(ncpmx+ipc) = selfp(ipc)
        end do


!-------------------- body surface charge
        do ipc=1,npc
          rhoindh(ipc) = rhoind(ipc)
          rhoind(ipc) = 0.0d0
!
          do icap=nscpmx(ipc),necpmx(ipc)
            xlocal = bdygrid(1,icap) - xl
            ylocal = bdygrid(2,icap) - yl
            zlocal = bdygrid(3,icap) - zl
!
            i = floor(xlocal)
            j = floor(ylocal)
            k = floor(zlocal)
!
            i1 = i + 1
            j1 = j + 1
            k1 = k + 1
!
            if(i1.ge.0.and.i.le.xu-xl.and. &
           &   j1.ge.0.and.j.le.yu-yl.and. &
           &   k1.ge.0.and.k.le.zu-zl) then
              x1 = xlocal - i
              y1 = ylocal - j
              z1 = zlocal - k
!             ------- temporary value
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
!
              v7 = xy1 * z2
              v6 = xz2 - v7
              v8 = yz2 - v7
              v5 = z2 - xz2 - v8
!             ------- distribute charge with weight
              rho(1,i ,j ,k1,1) = rho(1,i ,j ,k1,1) + v1*sfrho(icap)
              rho(1,i1,j ,k1,1) = rho(1,i1,j ,k1,1) + v2*sfrho(icap)
              rho(1,i1,j1,k1,1) = rho(1,i1,j1,k1,1) + v3*sfrho(icap)
              rho(1,i ,j1,k1,1) = rho(1,i ,j1,k1,1) + v4*sfrho(icap)
              rho(1,i ,j ,k ,1) = rho(1,i ,j ,k ,1) + v5*sfrho(icap)
              rho(1,i1,j ,k ,1) = rho(1,i1,j ,k ,1) + v6*sfrho(icap)
              rho(1,i1,j1,k ,1) = rho(1,i1,j1,k ,1) + v7*sfrho(icap)
              rho(1,i ,j1,k ,1) = rho(1,i ,j1,k ,1) + v8*sfrho(icap)
              rhobk(1,i ,j ,k1,3) = rhobk(1,i ,j ,k1,3) + v1*sfrho(icap)
              rhobk(1,i1,j ,k1,3) = rhobk(1,i1,j ,k1,3) + v2*sfrho(icap)
              rhobk(1,i1,j1,k1,3) = rhobk(1,i1,j1,k1,3) + v3*sfrho(icap)
              rhobk(1,i ,j1,k1,3) = rhobk(1,i ,j1,k1,3) + v4*sfrho(icap)
              rhobk(1,i ,j ,k ,3) = rhobk(1,i ,j ,k ,3) + v5*sfrho(icap)
              rhobk(1,i1,j ,k ,3) = rhobk(1,i1,j ,k ,3) + v6*sfrho(icap)
              rhobk(1,i1,j1,k ,3) = rhobk(1,i1,j1,k ,3) + v7*sfrho(icap)
              rhobk(1,i ,j1,k ,3) = rhobk(1,i ,j1,k ,3) + v8*sfrho(icap)
            end if
            sfrhoh(icap) = sfrhoh(icap) + sfrho(icap)
            rhoind(ipc) = rhoind(ipc) + sfrho(icap)
          end do
!
          selfp(ipc) = sfrho(ncpmx+ipc)
        end do
!
      end if


  return
  end subroutine surchg
