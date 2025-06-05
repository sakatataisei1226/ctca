#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine capmtx
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   C A P M T X
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   calculate 'capacity matrix'                            .
!   .   this matrix(cpmx) is to evaluate surface charge on     .
!   .   internal boundary                                      .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  use hdf
#define MCW local_comm
#define MSS MPI_STATUS_SIZE
#define MIP MPI_IN_PLACE
  implicit none
!
  integer(kind=4) :: m, ns,ne
  integer(kind=4) :: i,j,k, i1,j1,k1
  integer(kind=4) :: ipc,jpc, ipcg,jpcg, icap,jcap, mcap
  integer(kind=4) :: icon
  integer(kind=4) :: ipp(ixw), is
  integer(kind=4) :: src1, src2
  integer(kind=4) :: icomm, ierr
  integer(kind=4) :: xl,xu, yl,yu, zl,zu, siu,sju,sku
  integer(kind=4) :: ps
  integer(kind=4) :: nodeid
  real(kind=8) :: aamax, cpmsum
  real(kind=8) :: wok(ixw)
  real(kind=8) :: xlocal,ylocal,zlocal
  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8
  real(kind=8) :: phipeak(8,ixw)
!
  integer(kind=HID_T) :: fileid
  integer(kind=4) :: stats0,stats1
  integer(kind=8) :: dim(1)
  character(len=20) :: filename,dsname


!-------------------- 
      if(npc.ge.1) then
        ncpmx = necpmx(npc)
        nbdsf = nbdsf1(npc)
      else
        ncpmx = 0
        nbdsf = 0
      end if


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


!-------------------- 
      if(abs(modeww).lt.2.or.npc.eq.0) then
        return
!-------------------- 
      else if(abs(modeww).eq.2) then
        if(myid.eq.0) print*,"--bdygrid check--"
        call MPI_Barrier(MCW,ierr)
!
        allocate(cpmx(ncpmx,ncpmx))
        cpmx(:,:) = 0.0d0
!
JCAPL:  do jcap=1,ncpmx
          if(myid.eq.0) print '(A,i6,3f8.3)','jcap, i, j, k:', &
         &              jcap, bdygrid(1,jcap), bdygrid(2,jcap), bdygrid(3,jcap)
          poi(:,:,:,:,:) = 0.0d0
!
          xlocal = bdygrid(1,jcap) - sxl
          ylocal = bdygrid(2,jcap) - syl
          zlocal = bdygrid(3,jcap) - szl
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
!         ----------- distribute charge with weight
          if(i .ge.0.and.i .le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            poi(1,i ,j ,k1,1) = v1
          end if
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            poi(1,i1,j ,k1,1) = v2
          end if
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            poi(1,i1,j1,k1,1) = v3
          end if
          if(i .ge.0.and.i .le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            poi(1,i ,j1,k1,1) = v4
          end if
          if(i .ge.0.and.i .le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            poi(1,i ,j ,k ,1) = v5
          end if
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            poi(1,i1,j ,k ,1) = v6
          end if
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            poi(1,i1,j1,k ,1) = v7
          end if
          if(i .ge.0.and.i .le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            poi(1,i ,j1,k ,1) = v8
          end if
          call poisson(0)


!-------------------- store solution to array
!          mcap = 0
          diff(:) = 0.0d0
 ICAPL:   do icap=1,ncpmx
            xlocal = bdygrid(1,icap) - sxl
            ylocal = bdygrid(2,icap) - syl
            zlocal = bdygrid(3,icap) - szl
  
            i = floor(xlocal)
            j = floor(ylocal)
            k = floor(zlocal)
  
            i1 = i + 1
            j1 = j + 1
            k1 = k + 1
  
            x1 = xlocal - i
            y1 = ylocal - j
            z1 = zlocal - k
!
!            i = floor(bdygrid(1,icap))
!            j = floor(bdygrid(2,icap)) - syl
!            k = floor(bdygrid(3,icap)) - szl
!
!            i1 = i + 1
!            j1 = j + 1
!            k1 = k + 1
!
!            x1 = bdygrid(1,icap) - floor(bdygrid(1,icap))
!            y1 = bdygrid(2,icap) - floor(bdygrid(2,icap))
!            z1 = bdygrid(3,icap) - floor(bdygrid(3,icap))
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
!              if(jcap.eq.1) print*, "icap,myid,bgz",icap,myid,bdygrid(3,icap)
!              mcap = mcap + 1
!              sfrho(mcap) = poi(1,i ,j ,k1,1)*v1 &
!             &            + poi(1,i1,j ,k1,1)*v2 &
!             &            + poi(1,i1,j1,k1,1)*v3 &
!             &            + poi(1,i ,j1,k1,1)*v4 &
!             &            + poi(1,i ,j ,k ,1)*v5 &
!             &            + poi(1,i1,j ,k ,1)*v6 &
!             &            + poi(1,i1,j1,k ,1)*v7 &
!             &            + poi(1,i ,j1,k ,1)*v8
              diff(icap) = poi(1,i ,j ,k1,1)*v1 &
             &           + poi(1,i1,j ,k1,1)*v2 &
             &           + poi(1,i1,j1,k1,1)*v3 &
             &           + poi(1,i ,j1,k1,1)*v4 &
             &           + poi(1,i ,j ,k ,1)*v5 &
             &           + poi(1,i1,j ,k ,1)*v6 &
             &           + poi(1,i1,j1,k ,1)*v7 &
             &           + poi(1,i ,j1,k ,1)*v8
            end if
          end do ICAPL
!          call MPI_Allgatherv(sfrho,mcap,MPI_REAL8, &
!         &                    cpmx(:,jcap),ncpmxs,scpmxs,MPI_REAL8,MCW,ierr)
          call MPI_Allreduce(diff,cpmx(:,jcap),ncpmx,MPI_REAL8,MPI_SUM,MCW,ierr)
        end do JCAPL
        if(myid.eq.0) then
          open(5245,file='cpmxdiag1')
          do icap=1,ncpmx
            write(5245,*) cpmx(icap,icap)
          end do
          close(5245)
        end if


!-------------------- 
        phipeak(:,:) = 0.0d0
        do ipc=1,npc
        do icap=nmxcpmx(ipc),necpmx(ipc)
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

!-------------------- 
          poi(:,:,:,:,:) = 0.0d0
          if(i .ge.0.and.i .le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            poi(1,i,j,k1,1) = 1.0d0
          end if
          call poisson(0)
          if(i .ge.0.and.i .le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            phipeak(1,icap) = poi(1,i,j,k1,1)
          end if
!
          poi(:,:,:,:,:) = 0.0d0
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            poi(1,i1,j,k1,1) = 1.0d0
          end if
          call poisson(0)
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            phipeak(2,icap) = poi(1,i1,j,k1,1)
          end if
!
          poi(:,:,:,:,:) = 0.0d0
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            poi(1,i1,j1,k1,1) = 1.0d0
          end if
          call poisson(0)
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            phipeak(3,icap) = poi(1,i1,j1,k1,1)
          end if
!
          poi(:,:,:,:,:) = 0.0d0
          if(i .ge.0.and.i .le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            poi(1,i,j1,k1,1) = 1.0d0
          end if
          call poisson(0)
          if(i .ge.0.and.i .le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k1.ge.0.and.k1.le.sku) then
            phipeak(4,icap) = poi(1,i,j1,k1,1)
          end if
!
          poi(:,:,:,:,:) = 0.0d0
          if(i .ge.0.and.i .le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            poi(1,i,j,k,1) = 1.0d0
          end if
          call poisson(0)
          if(i .ge.0.and.i .le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            phipeak(5,icap) = poi(1,i,j,k,1)
          end if
!
          poi(:,:,:,:,:) = 0.0d0
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            poi(1,i1,j,k,1) = 1.0d0
          end if
          call poisson(0)
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j .ge.0.and.j .le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            phipeak(6,icap) = poi(1,i1,j,k,1)
          end if
!
          poi(:,:,:,:,:) = 0.0d0
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            poi(1,i1,j1,k,1) = 1.0d0
          end if
          call poisson(0)
          if(i1.ge.0.and.i1.le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            phipeak(7,icap) = poi(1,i1,j1,k,1)
          end if
!
          poi(:,:,:,:,:) = 0.0d0
          if(i .ge.0.and.i .le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            poi(1,i,j1,k,1) = 1.0d0
          end if
          call poisson(0)
          if(i .ge.0.and.i .le.siu.and. &
         &   j1.ge.0.and.j1.le.sju.and. &
         &   k .ge.0.and.k .le.sku) then
            phipeak(8,icap) = poi(1,i,j1,k,1)
          end if
        end do
        end do
!
        call MPI_Allreduce(MIP,phipeak,8*ncpmx,MPI_REAL8,MPI_SUM,MCW,ierr)
!
        do ipc=1,npc
        do icap=nmxcpmx(ipc),necpmx(ipc)
          x1 = bdygrid(1,icap) - floor(bdygrid(1,icap))
          y1 = bdygrid(2,icap) - floor(bdygrid(2,icap))
          z1 = bdygrid(3,icap) - floor(bdygrid(3,icap))
!         ----------- temporary value
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
          cpmx(icap,icap) = phipeak(1,icap)*v1 &
         &                + phipeak(2,icap)*v2 &
         &                + phipeak(3,icap)*v3 &
         &                + phipeak(4,icap)*v4 &
         &                + phipeak(5,icap)*v5 &
         &                + phipeak(6,icap)*v6 &
         &                + phipeak(7,icap)*v7 &
         &                + phipeak(8,icap)*v8
        end do
        end do
        if(myid.eq.0) then
          open(5245,file='cpmxdiag2')
          do icap=1,ncpmx
            write(5245,*) cpmx(icap,icap)
          end do
          close(5245)
        end if


!-------------------- check capmtx
        if(myid.eq.0) print*,'capacity matrix size =',ncpmx
        do jcap=1,ncpmx
          aamax = 0.0d0
          do icap=1,ncpmx
            if(abs(cpmx(icap,jcap)).gt.aamax) aamax = abs(cpmx(icap,jcap))
          end do
          if(myid.eq.0.and.aamax.eq.0.0d0) print *, 'Singular matrix. ', icap,jcap
        end do


!-------------------- calculate inverse 'cpmx'
        call dalu(cpmx,size(cpmx,1),ncpmx,0.0d0,ipp,is,wok,icon)
        if(myid.eq.0) print*,'icon(alu)  =',icon 
        call dluiv(cpmx,size(cpmx,1),ncpmx,ipp,icon)
        if(myid.eq.0) print*,'icon(luiv) =',icon
!       ------------- cpmx=D(i,j)
        cpmsum = 0.0d0
        do jcap=1,ncpmx
        do icap=1,ncpmx
          cpmsum = cpmsum + cpmx(icap,jcap)
        end do
        end do
        if(myid.eq.0) print*,'  cpmsum = ',cpmsum


!-------------------- 
        do jpc=1,npc
        do ipc=1,npc
          bpmx(ipc,jpc) = 0.0d0
          do jcap=nscpmx(jpc),necpmx(jpc)
          do icap=nscpmx(ipc),necpmx(ipc)
            bpmx(ipc,jpc) = bpmx(ipc,jpc) + cpmx(icap,jcap)
          end do
          end do
        end do
        end do


!-------------------- 
        apmx(:,:) = 0.0d0
        do jpcg=1,npcg
        do jpc=pcgs(jpcg)+1,pcgs(jpcg+1)
          do ipcg=1,npcg
          do ipc=ccgs(ipcg)+1,ccgs(ipcg+1)
            apmx(ipcg,jpcg) = apmx(ipcg,jpcg) + bpmx(ipc,jpc)
          end do
          end do
        end do
        end do


!-------------------- 
        do jpcg=1,npcg
        do ipcg=1,npcg
          if(mtd_vchg(ipcg).eq.1) then
            if(ipcg.eq.jpcg) then
              bcmx(ipcg,jpcg) = -1.0d0
            else if(mtd_vchg(jpcg).ne.1) then
              bcmx(ipcg,jpcg) =  apmx(ipcg,jpcg)
            else
              bcmx(ipcg,jpcg) =  0.0d0
            end if
          else
            if(mtd_vchg(jpcg).ne.1) then
              bcmx(ipcg,jpcg) =  apmx(ipcg,jpcg)              
            else
              bcmx(ipcg,jpcg) =  0.0d0
              i = 1
              do while(biasc(i)%to.ne.0)
                if(abs(biasc(i)%from).ge.ccgs(ipcg)+1.and. &
             &     abs(biasc(i)%from).le.ccgs(ipcg+1).and. &
                   abs(biasc(i)%to).ge.ccgs(jpcg)+1.and. &
             &     abs(biasc(i)%to).le.ccgs(jpcg+1)) then
                  bcmx(ipcg,jpcg) = 1.0d0
                end if
                i = i + 1
              end do
            end if
          end if
        end do
        end do
!
!
      else if(abs(modeww).eq.3) then
        filename = "./capmtx.h5"
        dim(1) = ncpmx
        call hdfopen(filename,fileid,DFACC_READ)
!
        if(myid.eq.0) then
          cpmsum = 0.0d0
          do jcap=1,ncpmx
            diff(:) = 0.0d0
            write(dsname,'(i8.8)') jcap
            if(ncpmx.gt.0) then
              call read1d(fileid,dsname,dim,diff(1:ncpmx),stats0,stats1)
            else
              diff(1:ncpmx) = 0.0d0
            end if
            cpmsum = cpmsum + sum(diff(1:ncpmx))
          end do
          print*,'  cpmsum = ',cpmsum
        end if
!
        do jpc=1,npc
        do ipc=1,npc
          bpmx(ipc,jpc) = 0.0d0
          do jcap=nscpmx(jpc),necpmx(jpc)
            diff(:) = 0.0d0
            write(dsname,'(i8.8)') jcap
            if(ncpmx.gt.0) then
              call read1d(fileid,dsname,dim,diff(1:ncpmx),stats0,stats1)
            else
              diff(1:ncpmx) = 0.0d0
            end if
            do icap=nscpmx(ipc),necpmx(ipc)
              bpmx(ipc,jpc) = bpmx(ipc,jpc) + diff(icap)
            end do
          end do
        end do
        end do
        diff(:) = 0.0d0
!
        apmx(:,:) = 0.0d0
        do jpcg=1,npcg
        do jpc=pcgs(jpcg)+1,pcgs(jpcg+1)
          do ipcg=1,npcg
          do ipc=ccgs(ipcg)+1,ccgs(ipcg+1)
            apmx(ipcg,jpcg) = apmx(ipcg,jpcg) + bpmx(ipc,jpc)
          end do
          end do
        end do
        end do
!
        do jpcg=1,npcg
        do ipcg=1,npcg
          if(mtd_vchg(ipcg).eq.1) then
            if(ipcg.eq.jpcg) then
              bcmx(ipcg,jpcg) = -1.0d0
            else if(mtd_vchg(jpcg).ne.1) then
              bcmx(ipcg,jpcg) =  apmx(ipcg,jpcg)
            else
              bcmx(ipcg,jpcg) =  0.0d0
            end if
          else
            if(mtd_vchg(jpcg).ne.1) then
              bcmx(ipcg,jpcg) =  apmx(ipcg,jpcg)              
            else
              bcmx(ipcg,jpcg) =  0.0d0
              i = 1
              do while(biasc(i)%to.ne.0)
                if(abs(biasc(i)%from).ge.ccgs(ipcg)+1.and. &
             &     abs(biasc(i)%from).le.ccgs(ipcg+1).and. &
                   abs(biasc(i)%to).ge.ccgs(jpcg)+1.and. &
             &     abs(biasc(i)%to).le.ccgs(jpcg+1)) then
                  bcmx(ipcg,jpcg) = 1.0d0
                end if
                i = i + 1
              end do
            end if
          end if
        end do
        end do
!
        nsdcpmx(1) = 1
        if(0.lt.mod(ncpmx,snode)) then
          nedcpmx(1) = ncpmx/snode + 1
        else
          nedcpmx(1) = ncpmx/snode
        end if
        ndcpmx(1) = nedcpmx(1) - nsdcpmx(1) + 1
        do nodeid=1,nnode-1
          nsdcpmx(nodeid+1) = nedcpmx(nodeid) + 1
          if(nodeid.lt.mod(ncpmx,snode)) then
            nedcpmx(nodeid+1) = nedcpmx(nodeid) + ncpmx/snode + 1
          else if(nodeid.lt.snode) then
            nedcpmx(nodeid+1) = nedcpmx(nodeid) + ncpmx/snode
          else
            nedcpmx(nodeid+1) = nedcpmx(nodeid)
          end if
          ndcpmx(nodeid+1) = nedcpmx(nodeid+1) - nsdcpmx(nodeid+1) + 1
        end do
        ndcpmx2(:) = ndcpmx(:)
        ndcpmx2(snode) = ndcpmx2(snode) + npc
        nsdcpmx2(:) = nsdcpmx(:)
        nsdcpmx2(snode+1:nnode) = nsdcpmx2(snode+1:nnode) + npc
        if(myid.eq.0.and.nedcpmx(snode).ne.ncpmx) then
          print*, "Warning!! nedcpmx(snode) & ncpmx inconsistent!!"
        end if
!
        allocate(cpmx(nsdcpmx(myid+1):nedcpmx(myid+1),ncpmx))
        cpmx(:,:) = 0.0d0
        do jcap=1,ncpmx
          write(dsname,'(i8.8)') jcap
          if(ncpmx.gt.0) then
            call read1d(fileid,dsname,dim,diff(1:ncpmx),stats0,stats1)
          else
            diff(1:ncpmx) = 0.0d0
          end if
          cpmx(nsdcpmx(myid+1):nedcpmx(myid+1),jcap) = &
         &  diff(nsdcpmx(myid+1):nedcpmx(myid+1))
        end do
        diff(:) = 0.0d0
!
        call hdfclose(fileid,stats0)
!
!      if(myid.eq.0) then
!        do jcap=1,5
!        do icap=1,5
!          print*, "check: myid,cpmx", myid,icap,jcap, cpmx(icap,jcap)
!        end do
!        end do
!      else if(myid.eq.1) then
!        do jcap=1,5
!        do icap=31,35
!          print*, "check: myid,cpmx", myid,icap,jcap, cpmx(icap,jcap)
!        end do
!        end do
!      end if
      end if


  return
  end subroutine capmtx
