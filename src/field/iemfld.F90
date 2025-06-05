#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine ibfield1
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   I B F I E L D 1
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   subroutine for solving linear equations for delta B    .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
#define MCW local_comm
#define MIP MPI_IN_PLACE
  implicit none
!
  integer(kind=4) :: i,j,k, ii,jj,kk, iii, ifloor, islice
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: iu,ju,ku
  integer(kind=4) :: siu,sju,sku, stiu
  integer(kind=4) :: ps
  integer(kind=4) :: sd, ifm, from, to
  real(kind=8) :: thet
  real(kind=8),parameter :: dlt=2.0d0, dts=4.0d0, ieps0=1.0d0
  real(kind=8) :: capc1
  real(kind=8) :: divisor

  integer(kind=4) :: lfloor,ufloor
  integer(kind=4) :: nxsbar,nysbar,nzsbar
  integer(kind=4) :: count, icomreq
  integer(kind=4) :: src,dst
  integer(kind=4) :: stag=0,rtag=0
  integer(kind=4) :: mpierr

  double precision,parameter :: freq=0.001d0
  double precision,parameter :: j0=100.0d0


!-------------------- 
      thet = gfactor


!-------------------- 
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)
      if(idxsd.ne.nodes(1)-1) then
        iu = xu - xl - 1
      else
        iu = xu - xl
      end if
      if(idysd.ne.nodes(2)-1) then
        ju = yu - yl - 1
      else
        ju = yu - yl
      end if
      if(idzsd.ne.nodes(3)-1) then
        ku = zu - zl - 1
      else
        ku = zu - zl
      end if
      if(pftmode.eq.3.and.myid.ne.snode-1) then
        siu = sxu - sxl - 1
      else
        siu = sxu - sxl
      end if
      if((pftmode.eq.1.or.pftmode.eq.2).and.myid.ne.snode-1) then
        stiu = stxu - stxl - 1
      else
        stiu = stxu - stxl
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
      rcnts(1:nnode-1) = slcs(1:nnode-1,3)
      rcnts(nnode) = slcs(nnode,3) + 1


!-------------------- 
      db(:,:,:,:,:) = 0.0d0
      dbs(:,:,:,:,:) = 0.0d0


!-------------------- 
      do k=0,ku
      do j=0,ju
      do i=0,iu
        db(CX,i,j,k,1) = thet*dts*(cs*(eb(BX,i+1,j,k,1) + eb(BX,i-1,j,k,1) &
       &                             + eb(BX,i,j+1,k,1) + eb(BX,i,j-1,k,1) &
       &                             + eb(BX,i,j,k+1,1) + eb(BX,i,j,k-1,1) &
       &                             - 6.0d0*eb(BX,i,j,k,1)) &
       &                         + ieps0*(-aj(JY,i,j,k+1,1) + aj(JY,i,j,k,1) &
       &                                 + aj(JZ,i,j+1,k,1) - aj(JZ,i,j,k,1))) &
       &               - dlt*(-eb(EY,i,j,k+1,1) + eb(EY,i,j,k,1) &
       &                     + eb(EZ,i,j+1,k,1) - eb(EZ,i,j,k,1))
        db(CY,i,j,k,1) = thet*dts*(cs*(eb(BY,i+1,j,k,1) + eb(BY,i-1,j,k,1) &
       &                             + eb(BY,i,j+1,k,1) + eb(BY,i,j-1,k,1) &
       &                             + eb(BY,i,j,k+1,1) + eb(BY,i,j,k-1,1) &
       &                             - 6.0d0*eb(BY,i,j,k,1)) &
       &                         + ieps0*(-aj(JZ,i+1,j,k,1) + aj(JZ,i,j,k,1) &
       &                                 + aj(JX,i,j,k+1,1) - aj(JX,i,j,k,1))) &
       &               - dlt*(-eb(EZ,i+1,j,k,1) + eb(EZ,i,j,k,1) &
       &                     + eb(EX,i,j,k+1,1) - eb(EX,i,j,k,1))
        db(CZ,i,j,k,1) = thet*dts*(cs*(eb(BZ,i+1,j,k,1) + eb(BZ,i-1,j,k,1) &
       &                             + eb(BZ,i,j+1,k,1) + eb(BZ,i,j-1,k,1) &
       &                             + eb(BZ,i,j,k+1,1) + eb(BZ,i,j,k-1,1) &
       &                             - 6.0d0*eb(BZ,i,j,k,1)) &
       &                         + ieps0*(-aj(JX,i,j+1,k,1) + aj(JX,i,j,k,1) &
       &                                 + aj(JY,i+1,j,k,1) - aj(JY,i,j,k,1))) &
       &               - dlt*(-eb(EX,i,j+1,k,1) + eb(EX,i,j,k,1) &
       &                     + eb(EY,i+1,j,k,1) - eb(EY,i,j,k,1))
      end do
      end do
      end do


!-------------------- 
      if(pftmode.eq.1) then
        icomreq = 0
        if(myid.lt.snode) then
!        lfloor = floor(dble(szl/nzsd)); ufloor = floor(dble((szu-1)/nzsd))
          lfloor = int(szl/nzsd); ufloor = int((szu-1)/nzsd)
          do ifloor=lfloor,ufloor
            count = min(szu,(ifloor+1)*nzsd) - max(szl,ifloor*nzsd)
            if(myid.eq.snode-1.and.ifloor.eq.ufloor) &
           &   count = count + 1
            do j=0,nodes(2)-1
            do i=0,nodes(1)-1
              src = i + j*nodes(1) + ifloor*nodes(1)*nodes(2)
              icomreq = icomreq + 1
              call MPI_Irecv(dbs(1,i*nxsd,j*nysd,max(0,ifloor*nzsd-szl),1), &
             &               count, mptype_rs(xtype(i),ytype(j),1,1), &
             &               src, 0, MCW, ireqs(icomreq), mpierr)
            end do
            end do
          end do
        end if
!
!      do dst=floor(dble(zl/nzslc)),floor(dble((zu-1)/nzslc))
        do dst=int(zl/nzslc),int((zu-1)/nzslc)
          count = 0
          do islice=dst*nzslc,(dst+1)*nzslc-1
            if(zl.le.islice.and.islice.lt.zu) &
           &  count = count + 1
          end do
          if(dst.eq.snode-1) count = count + 1
          icomreq = icomreq + 1
          call MPI_Isend(db(1,0,0,max(dst*nzslc-zl,0),1), &
         &               count, mptype_sc(xtype(idxsd),ytype(idysd),1,1), &
         &               dst, 0, MCW, ireqs(icomreq), mpierr)
        end do
!
        call MPI_Waitall(icomreq,ireqs,istatus,mpierr)
!
!     --------------- 
      else if(pftmode.eq.2) then
        icomreq = 0
        if(myid.lt.snode) then
!        lfloor = floor(dble(syl/nysd)); ufloor = floor(dble((syu-1)/nysd))
          lfloor = int(syl/nysd); ufloor = int((syu-1)/nysd)
          do ifloor=lfloor,ufloor
            count = min(syu,(ifloor+1)*nysd) - max(syl,ifloor*nysd)
            if(myid.eq.snode-1.and.ifloor.eq.ufloor) &
           &   count = count + 1
            do k=0,nodes(3)-1
            do i=0,nodes(1)-1
              src = i + ifloor*nodes(1) + k*nodes(1)*nodes(2)
              icomreq = icomreq + 1
              call MPI_Irecv(dbs(1,i*nxsd,max(0,ifloor*nysd-syl),k*nzsd,1), &
             &               count, mptype_rs(xtype(i),ztype(k),1,2), &
             &               src, 0, MCW, ireqs(icomreq), mpierr)
            end do
            end do
          end do
        end if
!
!      do dst=floor(dble(yl/nyslc)),floor(dble((yu-1)/nyslc))
        do dst=int(yl/nyslc),int((yu-1)/nyslc)
          count = 0
          do islice=dst*nyslc,(dst+1)*nyslc-1
            if(yl.le.islice.and.islice.lt.yu) &
           &  count = count + 1
          end do
          if(dst.eq.snode-1) count = count + 1
          icomreq = icomreq + 1
          call MPI_Isend(db(1,0,max(dst*nyslc-yl,0),0,1), &
         &               count, mptype_sc(xtype(idxsd),ztype(idzsd),1,2), &
         &               dst, 0, MCW, ireqs(icomreq), mpierr)
        end do
!
        call MPI_Waitall(icomreq,ireqs,istatus,mpierr)
!
!     --------------- 
      else if(pftmode.eq.3) then
        dbs(CX:CZ,0:iu,0:ny,0:nz,1) = db(CX:CZ,0:iu,0:ny,0:nz,1)
!
!     --------------- 
      else if(pftmode.eq.4) then
        dbs(CX:CZ,0:nx,0:ny,0:ku,1) = db(CX:CZ,0:nx,0:ny,0:ku,1)
      end if


!-------------------- 
      if(myid.lt.snode) then
      if(pftmode.eq.1) then
!       ------------- FFT-xy
        do k=0-szl,nz-szl
          if(k.ge.0.and.k.le.sku) then
            call dfftw_execute_r2r(fftplan(CX,1,1), &
           &  dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1), &
           &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1))
            dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1) = &
           &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1)
            call dfftw_execute_r2r(fftplan(CY,1,1), &
           &  dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1), &
           &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1))
            dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1) = &
           &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1)
            call dfftw_execute_r2r(fftplan(CZ,1,1), &
           &  dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1), &
           &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1))
            dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1) = &
           &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1)
        end if
        end do
!       ------------- transpose x-z
        call MPI_Alltoall(dbs(1,0,0,0,1),1,mptype_aa(0,0,1,1), &
       &                  dbs(1,0,0,0,2),1,mptype_aa(1,0,1,1), &
       &                  subcomm,mpierr)
        iii = 0
        do ii=0,(nxslc+1)*(snode-1),nxslc+1
          nxsbar = stiu
          if(ii.eq.(nxslc+1)*(snode-1)) then
            nzsbar = nzslc
          else
            nzsbar = nzslc - 1
          end if
          do k=0,nxsbar
            do i=0,nzsbar
              dbs(CX:CZ,iii+i,0:ny,k,1) = dbs(CX:CZ,ii+k,0:ny,i,2)
            end do
          end do
          iii = iii + nzslc
        end do
!       ------------- FFT-z
        do k=0-stxl,nx-stxl
        do j=0,ny
          if(k.ge.0.and.k.le.stiu) then
            call dfftw_execute_r2r(fftplan(CX,2,1), &
           &  dbs(CX,lzfft(CX):uzfft(CX),j,k,1), &
           &  warray1d(lzfft(CX):uzfft(CX),1))
            dbs(CX,lzfft(CX):uzfft(CX),j,k,1) = &
           &  warray1d(lzfft(CX):uzfft(CX),1)
            call dfftw_execute_r2r(fftplan(CY,2,1), &
           &  dbs(CY,lzfft(CY):uzfft(CY),j,k,1), &
           &  warray1d(lzfft(CY):uzfft(CY),1))
            dbs(CY,lzfft(CY):uzfft(CY),j,k,1) = &
           &  warray1d(lzfft(CY):uzfft(CY),1)
            call dfftw_execute_r2r(fftplan(CZ,2,1), &
           &  dbs(CZ,lzfft(CZ):uzfft(CZ),j,k,1), &
           &  warray1d(lzfft(CZ):uzfft(CZ),1))
            dbs(CZ,lzfft(CZ):uzfft(CZ),j,k,1) = &
           &  warray1d(lzfft(CZ):uzfft(CZ),1)
          end if
        end do
        end do
!
!     --------------- 
      else if(pftmode.eq.2) then
!       ------------- FFT-xz
        do j=0-syl,ny-syl
          if(j.ge.0.and.j.le.sju) then
            call dfftw_execute_r2r(fftplan(CX,1,1), &
           &  dbs(CX,lxfft(CX):uxfft(CX),j,lzfft(CX):uzfft(CX),1), &
           &  warray2d(lxfft(CX):uxfft(CX),lzfft(CX):uzfft(CX),1))
            dbs(CX,lxfft(CX):uxfft(CX),j,lzfft(CX):uzfft(CX),1) = &
           &  warray2d(lxfft(CX):uxfft(CX),lzfft(CX):uzfft(CX),1)
            call dfftw_execute_r2r(fftplan(CY,1,1), &
           &  dbs(CY,lxfft(CY):uxfft(CY),j,lzfft(CY):uzfft(CY),1), &
           &  warray2d(lxfft(CY):uxfft(CY),lzfft(CY):uzfft(CY),1))
            dbs(CY,lxfft(CY):uxfft(CY),j,lzfft(CY):uzfft(CY),1) = &
           &  warray2d(lxfft(CY):uxfft(CY),lzfft(CY):uzfft(CY),1)
            call dfftw_execute_r2r(fftplan(CZ,1,1), &
           &  dbs(CZ,lxfft(CZ):uxfft(CZ),j,lzfft(CZ):uzfft(CZ),1), &
           &  warray2d(lxfft(CZ):uxfft(CZ),lzfft(CZ):uzfft(CZ),1))
            dbs(CZ,lxfft(CZ):uxfft(CZ),j,lzfft(CZ):uzfft(CZ),1) = &
           &  warray2d(lxfft(CZ):uxfft(CZ),lzfft(CZ):uzfft(CZ),1)
        end if
        end do
!       ------------- transpose x-y
        call MPI_Alltoall(dbs(1,0,0,0,1),1,mptype_aa(0,0,1,2), &
       &                  dbs(1,0,0,0,2),1,mptype_aa(1,0,1,2), &
       &                  subcomm,mpierr)
        iii = 0
        do ii=0,(nxslc+1)*(snode-1),nxslc+1
          nxsbar = stiu
          if(ii.eq.(nxslc+1)*(snode-1)) then
            nysbar = nyslc
          else
            nysbar = nyslc - 1
          end if
          do j=0,nxsbar
            do i=0,nysbar
              dbs(CX:CZ,iii+i,j,0:nz,1) = dbs(CX:CZ,ii+j,i,0:nz,2)
            end do
          end do
          iii = iii + nyslc
        end do
!       ------------- FFT-y
        do k=0,nz
        do j=0-stxl,nx-stxl
          if(j.ge.0.and.j.le.stiu) then
            call dfftw_execute_r2r(fftplan(CX,2,1), &
           &  dbs(CX,lyfft(CX):uyfft(CX),j,k,1), &
           &  warray1d(lyfft(CX):uyfft(CX),1))
            dbs(CX,lyfft(CX):uyfft(CX),j,k,1) = &
           &  warray1d(lyfft(CX):uyfft(CX),1)
            call dfftw_execute_r2r(fftplan(CY,2,1), &
           &  dbs(CY,lyfft(CY):uyfft(CY),j,k,1), &
           &  warray1d(lyfft(CY):uyfft(CY),1))
            dbs(CY,lyfft(CY):uyfft(CY),j,k,1) = &
           &  warray1d(lyfft(CY):uyfft(CY),1)
            call dfftw_execute_r2r(fftplan(CZ,2,1), &
           &  dbs(CZ,lyfft(CZ):uyfft(CZ),j,k,1), &
           &  warray1d(lyfft(CZ):uyfft(CZ),1))
            dbs(CZ,lyfft(CZ):uyfft(CZ),j,k,1) = &
           &  warray1d(lyfft(CZ):uyfft(CZ),1)
          end if
        end do
        end do
!
!     --------------- 
      else if(pftmode.eq.3) then
!       ------------- FFT-yz
        do i=0-sxl,nx-sxl
        if(i.ge.0.and.i.le.siu) then
          call dfftw_execute_r2r(fftplan(CX,1,1), &
         &  dbs(CX,i,lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1), &
         &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1))
          dbs(CX,i,lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1) = &
         &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,1,1), &
         &  dbs(CY,i,lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1), &
         &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1))
          dbs(CY,i,lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1) = &
         &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,1,1), &
         &  dbs(CZ,i,lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1), &
         &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1))
          dbs(CZ,i,lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1) = &
         &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1)
        end if
        end do
!       ------------- gather in x-dir
        call MPI_Allgatherv(dbs(1,0,-1,-1,1),siu+1,mptype_yz(1), &
       &                    dbs(1,0,-1,-1,2),rcnts,slcs(:,1),mptype_yz(1), &
       &                    subcomm,mpierr)
        dbs(:,:,:,:,1) = 0.0d0
!       ------------- FFT-x
        do k=0,nz
        do j=0,ny
          call dfftw_execute_r2r(fftplan(CX,2,1), &
         &  dbs(CX,lxfft(CX):uxfft(CX),j,k,2), &
         &  warray1d(lxfft(CX):uxfft(CX),1))
          dbs(CX,lxfft(CX):uxfft(CX),j,k,1) = &
         &  warray1d(lxfft(CX):uxfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,2,1), &
         &  dbs(CY,lxfft(CY):uxfft(CY),j,k,2), &
         &  warray1d(lxfft(CY):uxfft(CY),1))
          dbs(CY,lxfft(CY):uxfft(CY),j,k,1) = &
         &  warray1d(lxfft(CY):uxfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,2,1), &
         &  dbs(CZ,lxfft(CZ):uxfft(CZ),j,k,2), &
         &  warray1d(lxfft(CZ):uxfft(CZ),1))
          dbs(CZ,lxfft(CZ):uxfft(CZ),j,k,1) = &
         &  warray1d(lxfft(CZ):uxfft(CZ),1)
        end do
        end do
!
!     --------------- 
      else if(pftmode.eq.4) then
!       ------------- FFT-xy
        do k=0-szl,nz-szl
        if(k.ge.0.and.k.le.sku) then
          call dfftw_execute_r2r(fftplan(CX,1,1), &
         &  dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1), &
         &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1))
          dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1) = &
         &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,1,1), &
         &  dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1), &
         &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1))
          dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1) = &
         &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,1,1), &
         &  dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1), &
         &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1))
          dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1) = &
         &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1)
        end if
        end do
!       ------------- gather in z-dir
        call MPI_Allgatherv(dbs(1,-1,-1,0,1),sku+1,mptype_xy(1), &
       &                    dbs(1,-1,-1,0,2),rcnts,slcs(:,1),mptype_xy(1), &
       &                    subcomm,mpierr)
!       ------------- FFT-z
        do j=0,ny
        do i=0,nx
          call dfftw_execute_r2r(fftplan(CX,2,1), &
         &  dbs(CX,i,j,lzfft(CX):uzfft(CX),2), &
         &  warray1d(lzfft(CX):uzfft(CX),1))
          dbs(CX,i,j,lzfft(CX):uzfft(CX),1) = &
         &  warray1d(lzfft(CX):uzfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,2,1), &
         &  dbs(CY,i,j,lzfft(CY):uzfft(CY),2), &
         &  warray1d(lzfft(CY):uzfft(CY),1))
          dbs(CY,i,j,lzfft(CY):uzfft(CY),1) = &
         &  warray1d(lzfft(CY):uzfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,2,1), &
         &  dbs(CZ,i,j,lzfft(CZ):uzfft(CZ),2), &
         &  warray1d(lzfft(CZ):uzfft(CZ),1))
          dbs(CZ,i,j,lzfft(CZ):uzfft(CZ),1) = &
         &  warray1d(lzfft(CZ):uzfft(CZ),1)
        end do
        end do
      end if
      end if


!-------------------- 
      if(myid.lt.snode) then
      if(pftmode.eq.1) then
        do k=0-stxl,nx-stxl
        do j=0,ny
        do i=0,nz
          if(k.ge.0.and.k.le.stiu) then
            kk = k + stxl
            divisor = (1.0d0 + kmod(1,kk,1) + kmod(2,j,1) + kmod(3,i,1))*mfactor(1)
            if(divisor.ne.0.0d0) then
              dbs(CX,i,j,k,1) = dbs(CX,i,j,k,1)/divisor
              dbs(CY,i,j,k,1) = dbs(CY,i,j,k,1)/divisor
              dbs(CZ,i,j,k,1) = dbs(CZ,i,j,k,1)/divisor
            else
              dbs(CX,i,j,k,1) = 0.0d0
              dbs(CY,i,j,k,1) = 0.0d0
              dbs(CZ,i,j,k,1) = 0.0d0
            end if
          end if
        end do
        end do
        end do
      else if(pftmode.eq.2) then
        do k=0,nz
        do j=0-stxl,nx-stxl
        do i=0,ny
          if(j.ge.0.and.j.le.stiu) then
            jj = j + stxl
            divisor = (1.0d0 + kmod(1,jj,1) + kmod(2,i,1) + kmod(3,k,1))*mfactor(1)
            if(divisor.ne.0.0d0) then
              dbs(CX,i,j,k,1) = dbs(CX,i,j,k,1)/divisor
              dbs(CY,i,j,k,1) = dbs(CY,i,j,k,1)/divisor
              dbs(CZ,i,j,k,1) = dbs(CZ,i,j,k,1)/divisor
            else
              dbs(CX,i,j,k,1) = 0.0d0
              dbs(CY,i,j,k,1) = 0.0d0
              dbs(CZ,i,j,k,1) = 0.0d0
            end if
          end if
        end do
        end do
        end do
      else if(pftmode.eq.3.or.pftmode.eq.4) then
        do k=0,nz
        do j=0,ny
        do i=0,nx
          divisor = (1.0d0 + kmod(1,i,1) + kmod(2,j,1) + kmod(3,k,1))*mfactor(1)
          if(divisor.ne.0.0d0) then
            dbs(CX,i,j,k,1) = dbs(CX,i,j,k,1)/divisor
            dbs(CY,i,j,k,1) = dbs(CY,i,j,k,1)/divisor
            dbs(CZ,i,j,k,1) = dbs(CZ,i,j,k,1)/divisor
          else
            dbs(CX,i,j,k,1) = 0.0d0
            dbs(CY,i,j,k,1) = 0.0d0
            dbs(CZ,i,j,k,1) = 0.0d0
          end if
        end do
        end do
        end do
      end if
      end if


!-------------------- 
      if(myid.lt.snode) then
      if(pftmode.eq.1) then
!       ------------- IFFT-z
        do k=0-stxl,nx-stxl
        do j=0,ny
          if(k.ge.0.and.k.le.stiu) then
            call dfftw_execute_r2r(fftplan(CX,2,2), &
           &  dbs(CX,lzfft(CX):uzfft(CX),j,k,1), &
           &  warray1d(lzfft(CX):uzfft(CX),1))
            dbs(CX,lzfft(CX):uzfft(CX),j,k,1) = &
           &  warray1d(lzfft(CX):uzfft(CX),1)
            call dfftw_execute_r2r(fftplan(CY,2,2), &
           &  dbs(CY,lzfft(CY):uzfft(CY),j,k,1), &
           &  warray1d(lzfft(CY):uzfft(CY),1))
            dbs(CY,lzfft(CY):uzfft(CY),j,k,1) = &
           &  warray1d(lzfft(CY):uzfft(CY),1)
            call dfftw_execute_r2r(fftplan(CZ,2,2), &
           &  dbs(CZ,lzfft(CZ):uzfft(CZ),j,k,1), &
           &  warray1d(lzfft(CZ):uzfft(CZ),1))
            dbs(CZ,lzfft(CZ):uzfft(CZ),j,k,1) = &
           &  warray1d(lzfft(CZ):uzfft(CZ),1)
          end if
        end do
        end do
!       ------------- transpose z-x
        call MPI_Alltoall(dbs(1,0,0,0,1),1,mptype_aa(0,1,1,1), &
       &                  dbs(1,0,0,0,2),1,mptype_aa(1,1,1,1), &
       &                  subcomm,mpierr)
        iii = 0
        do ii=0,(nzslc+1)*(snode-1),nzslc+1
          if(ii.eq.(nzslc+1)*(snode-1)) then
            nxsbar = nxslc
          else
            nxsbar = nxslc - 1
          end if
          nzsbar = sku
          do k=0,nzsbar
            do i=0,nxsbar
              dbs(CX:CZ,iii+i,0:ny,k,1) = dbs(CX:CZ,ii+k,0:ny,i,2)
            end do
          end do
          iii = iii + nxslc
        end do
!       ------------- IFFT-xy
        do k=0-szl,nz-szl
          if(k.ge.0.and.k.le.sku) then
            call dfftw_execute_r2r(fftplan(CX,1,2), &
           &  dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1), &
           &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1))
            dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1) = &
           &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1)
            call dfftw_execute_r2r(fftplan(CY,1,2), &
           &  dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1), &
           &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1))
            dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1) = &
           &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1)
            call dfftw_execute_r2r(fftplan(CZ,1,2), &
           &  dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1), &
           &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1))
            dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1) = &
           &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1)
          end if
        end do
!
!     --------------- 
      else if(pftmode.eq.2) then
!       ------------- IFFT-y
        do k=0,nz
        do j=0-stxl,nx-stxl
          if(j.ge.0.and.j.le.stiu) then
            call dfftw_execute_r2r(fftplan(CX,2,2), &
           &  dbs(CX,lyfft(CX):uyfft(CX),j,k,1), &
           &  warray1d(lyfft(CX):uyfft(CX),1))
            dbs(CX,lyfft(CX):uyfft(CX),j,k,1) = &
           &  warray1d(lyfft(CX):uyfft(CX),1)
            call dfftw_execute_r2r(fftplan(CY,2,2), &
           &  dbs(CY,lyfft(CY):uyfft(CY),j,k,1), &
           &  warray1d(lyfft(CY):uyfft(CY),1))
            dbs(CY,lyfft(CY):uyfft(CY),j,k,1) = &
           &  warray1d(lyfft(CY):uyfft(CY),1)
            call dfftw_execute_r2r(fftplan(CZ,2,2), &
           &  dbs(CZ,lyfft(CZ):uyfft(CZ),j,k,1), &
           &  warray1d(lyfft(CZ):uyfft(CZ),1))
            dbs(CZ,lyfft(CZ):uyfft(CZ),j,k,1) = &
           &  warray1d(lyfft(CZ):uyfft(CZ),1)
          end if
        end do
        end do
!       ------------- transpose y-x
        call MPI_Alltoall(dbs(1,0,0,0,1),1,mptype_aa(0,1,1,2), &
       &                  dbs(1,0,0,0,2),1,mptype_aa(1,1,1,2), &
       &                  subcomm,mpierr)
        iii = 0
        do ii=0,(nyslc+1)*(snode-1),nyslc+1
          if(ii.eq.(nyslc+1)*(snode-1)) then
            nxsbar = nxslc
          else
            nxsbar = nxslc - 1
          end if
          nysbar = sju
          do j=0,nysbar
            do i=0,nxsbar
              dbs(CX:CZ,iii+i,j,0:nz,1) = dbs(CX:CZ,ii+j,i,0:nz,2)
            end do
          end do
          iii = iii + nxslc
        end do
!       ------------- IFFT-xz
        do j=0-syl,ny-syl
          if(j.ge.0.and.j.le.sju) then
            call dfftw_execute_r2r(fftplan(CX,1,2), &
           &  dbs(CX,lxfft(CX):uxfft(CX),j,lzfft(CX):uzfft(CX),1), &
           &  warray2d(lxfft(CX):uxfft(CX),lzfft(CX):uzfft(CX),1))
            dbs(CX,lxfft(CX):uxfft(CX),j,lzfft(CX):uzfft(CX),1) = &
           &  warray2d(lxfft(CX):uxfft(CX),lzfft(CX):uzfft(CX),1)
            call dfftw_execute_r2r(fftplan(CY,1,2), &
           &  dbs(CY,lxfft(CY):uxfft(CY),j,lzfft(CY):uzfft(CY),1), &
           &  warray2d(lxfft(CY):uxfft(CY),lzfft(CY):uzfft(CY),1))
            dbs(CY,lxfft(CY):uxfft(CY),j,lzfft(CY):uzfft(CY),1) = &
           &  warray2d(lxfft(CY):uxfft(CY),lzfft(CY):uzfft(CY),1)
            call dfftw_execute_r2r(fftplan(CZ,1,2), &
           &  dbs(CZ,lxfft(CZ):uxfft(CZ),j,lzfft(CZ):uzfft(CZ),1), &
           &  warray2d(lxfft(CZ):uxfft(CZ),lzfft(CZ):uzfft(CZ),1))
            dbs(CZ,lxfft(CZ):uxfft(CZ),j,lzfft(CZ):uzfft(CZ),1) = &
           &  warray2d(lxfft(CZ):uxfft(CZ),lzfft(CZ):uzfft(CZ),1)
          end if
        end do
!
!     --------------- 
      else if(pftmode.eq.3) then
!       ------------- IFFT-x
        do k=0,nz
        do j=0,ny
          call dfftw_execute_r2r(fftplan(CX,2,2), &
         &  dbs(CX,lxfft(CX):uxfft(CX),j,k,1), &
         &  warray1d(lxfft(CX):uxfft(CX),1))
          dbs(CX,lxfft(CX):uxfft(CX),j,k,1) = &
         &  warray1d(lxfft(CX):uxfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,2,2), &
         &  dbs(CY,lxfft(CY):uxfft(CY),j,k,1), &
         &  warray1d(lxfft(CY):uxfft(CY),1))
          dbs(CY,lxfft(CY):uxfft(CY),j,k,1) = &
         &  warray1d(lxfft(CY):uxfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,2,2), &
         &  dbs(CZ,lxfft(CZ):uxfft(CZ),j,k,1), &
         &  warray1d(lxfft(CZ):uxfft(CZ),1))
          dbs(CZ,lxfft(CZ):uxfft(CZ),j,k,1) = &
         &  warray1d(lxfft(CZ):uxfft(CZ),1)
        end do
        end do
!       ------------- IFFT-yz
        if(nfbnd(1).eq.0.and.myid.eq.snode-1) then
        do i=0,1
          call dfftw_execute_r2r(fftplan(CX,1,2), &
         &  dbs(CX,i,lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1), &
         &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1))
          dbs(CX,i,lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1) = &
         &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,1,2), &
         &  dbs(CY,i,lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1), &
         &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1))
          dbs(CY,i,lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1) = &
         &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,1,2), &
         &  dbs(CZ,i,lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1), &
         &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1))
          dbs(CZ,i,lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1) = &
         &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1)
        end do
        end if
        do i=max(sxl-1,0),min(sxu+1,nx)
          call dfftw_execute_r2r(fftplan(CX,1,2), &
         &  dbs(CX,i,lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1), &
         &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1))
          dbs(CX,i,lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1) = &
         &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,1,2), &
         &  dbs(CY,i,lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1), &
         &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1))
          dbs(CY,i,lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1) = &
         &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,1,2), &
         &  dbs(CZ,i,lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1), &
         &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1))
          dbs(CZ,i,lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1) = &
         &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1)
        end do
        if(nfbnd(1).eq.0.and.myid.eq.0) then
          i = nx - 1
          call dfftw_execute_r2r(fftplan(CX,1,2), &
         &  dbs(CX,i,lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1), &
         &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1))
          dbs(CX,i,lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1) = &
         &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,1,2), &
         &  dbs(CY,i,lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1), &
         &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1))
          dbs(CY,i,lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1) = &
         &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,1,2), &
         &  dbs(CZ,i,lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1), &
         &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1))
          dbs(CZ,i,lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1) = &
         &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1)
        end if
!
!     --------------- 
      else if(pftmode.eq.4) then
!       ------------- IFFT-z
        do j=0,ny
        do i=0,nx
          call dfftw_execute_r2r(fftplan(CX,2,2), &
         &  dbs(CX,i,j,lzfft(CX):uzfft(CX),1), &
         &  warray1d(lzfft(CX):uzfft(CX),1))
          dbs(CX,i,j,lzfft(CX):uzfft(CX),1) = &
         &  warray1d(lzfft(CX):uzfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,2,2), &
         &  dbs(CY,i,j,lzfft(CY):uzfft(CY),1), &
         &  warray1d(lzfft(CY):uzfft(CY),1))
          dbs(CY,i,j,lzfft(CY):uzfft(CY),1) = &
         &  warray1d(lzfft(CY):uzfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,2,2), &
         &  dbs(CZ,i,j,lzfft(CZ):uzfft(CZ),1), &
         &  warray1d(lzfft(CZ):uzfft(CZ),1))
          dbs(CZ,i,j,lzfft(CZ):uzfft(CZ),1) = &
         &  warray1d(lzfft(CZ):uzfft(CZ),1)
        end do
        end do
!       ------------- IFFT-xy
        if(nfbnd(3).eq.0.and.myid.eq.snode-1) then
        do k=0,1
          call dfftw_execute_r2r(fftplan(CX,1,2), &
         &  dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1), &
         &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1))
          dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1) = &
         &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,1,2), &
         &  dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1), &
         &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1))
          dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1) = &
         &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,1,2), &
         &  dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1), &
         &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1))
          dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1) = &
         &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1)
        end do
        end if
        do k=max(szl-1,0),min(szu+1,nz)
          call dfftw_execute_r2r(fftplan(CX,1,2), &
         &  dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1), &
         &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1))
          dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1) = &
         &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,1,2), &
         &  dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1), &
         &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1))
          dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1) = &
         &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,1,2), &
         &  dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1), &
         &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1))
          dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1) = &
         &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1)
        end do
        if(nfbnd(3).eq.0.and.myid.eq.0) then
          k = nz - 1
          call dfftw_execute_r2r(fftplan(CX,1,2), &
         &  dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1), &
         &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1))
          dbs(CX,lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),k,1) = &
         &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1)
          call dfftw_execute_r2r(fftplan(CY,1,2), &
         &  dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1), &
         &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1))
          dbs(CY,lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),k,1) = &
         &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1)
          call dfftw_execute_r2r(fftplan(CZ,1,2), &
         &  dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1), &
         &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1))
          dbs(CZ,lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),k,1) = &
         &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1)
        end if
      end if
      end if


!-------------------- 
      if(myid.lt.snode) then
      if(pftmode.eq.1) then
        dbs(CX,         -1:lxfft(CX)-1,:,:,1) = 0.0d0
        dbs(CY,         -1:lxfft(CY)-1,:,:,1) = 0.0d0
        dbs(CZ,         -1:lxfft(CZ)-1,:,:,1) = 0.0d0
!
        dbs(CX,uxfft(CX)+1:nx       +1,:,:,1) = 0.0d0
        dbs(CY,uxfft(CY)+1:nx       +1,:,:,1) = 0.0d0
        dbs(CZ,uxfft(CZ)+1:nx       +1,:,:,1) = 0.0d0
!
        dbs(CX,:,         -1:lyfft(CX)-1,:,1) = 0.0d0
        dbs(CY,:,         -1:lyfft(CY)-1,:,1) = 0.0d0
        dbs(CZ,:,         -1:lyfft(CZ)-1,:,1) = 0.0d0
!
        dbs(CX,:,uyfft(CX)+1:ny       +1,:,1) = 0.0d0
        dbs(CY,:,uyfft(CY)+1:ny       +1,:,1) = 0.0d0
        dbs(CZ,:,uyfft(CZ)+1:ny       +1,:,1) = 0.0d0
!
        if(myid.eq.0) then
          dbs(CX,:,:,         -1    :lzfft(CX)-1,1) = 0.0d0
          dbs(CY,:,:,         -1    :lzfft(CY)-1,1) = 0.0d0
          dbs(CZ,:,:,         -1    :lzfft(CZ)-1,1) = 0.0d0
        end if
        if(myid.eq.snode-1) then
          dbs(CX,:,:,uzfft(CX)+1-szl:nhslc    +1,1) = 0.0d0
          dbs(CY,:,:,uzfft(CY)+1-szl:nhslc    +1,1) = 0.0d0
          dbs(CZ,:,:,uzfft(CZ)+1-szl:nhslc    +1,1) = 0.0d0
        end if
!
        if(nfbnd(1).eq.0) then
          dbs(CX:CZ,  -1:  -1,:,:,1) = dbs(CX:CZ,nx-1:nx-1,:,:,1)
          dbs(CX:CZ,nx  :nx  ,:,:,1) = dbs(CX:CZ,   0:   0,:,:,1)
          dbs(CX:CZ,nx+1:nx+1,:,:,1) = dbs(CX:CZ,  +1:  +1,:,:,1)
        end if
!
        if(nfbnd(2).eq.0) then
          dbs(CX:CZ,:,  -1:  -1,:,1) = dbs(CX:CZ,:,ny-1:ny-1,:,1)
          dbs(CX:CZ,:,ny  :ny  ,:,1) = dbs(CX:CZ,:,   0:   0,:,1)
          dbs(CX:CZ,:,ny+1:ny+1,:,1) = dbs(CX:CZ,:,  +1:  +1,:,1)
        end if
!
        call MPI_sendrecv(dbs(1,-1,-1,szu-szl-1,1),1, &
       &                  mptype_xy(1),uslice(1),stag, &
       &                  dbs(1,-1,-1,-1,1),1, &
       &                  mptype_xy(1),lslice(1),rtag, &
       &                  subcomm, MPI_STATUS_IGNORE, mpierr)
        call MPI_sendrecv(dbs(1,-1,-1,0,1),1, &
       &                  mptype_xy(1),lslice(1),stag, &
       &                  dbs(1,-1,-1,szu-szl,1),1, &
       &                  mptype_xy(1),uslice(1),rtag, &
       &                  subcomm, MPI_STATUS_IGNORE, mpierr)
!
!     --------------- 
      else if(pftmode.eq.2) then
        dbs(CX,         -1:lxfft(CX)-1,:,:,1) = 0.0d0
        dbs(CY,         -1:lxfft(CY)-1,:,:,1) = 0.0d0
        dbs(CZ,         -1:lxfft(CZ)-1,:,:,1) = 0.0d0
!
        dbs(CX,uxfft(CX)+1:nx       +1,:,:,1) = 0.0d0
        dbs(CY,uxfft(CY)+1:nx       +1,:,:,1) = 0.0d0
        dbs(CZ,uxfft(CZ)+1:nx       +1,:,:,1) = 0.0d0
!
        dbs(CX,:,:,         -1:lzfft(CX)-1,1) = 0.0d0
        dbs(CY,:,:,         -1:lzfft(CY)-1,1) = 0.0d0
        dbs(CZ,:,:,         -1:lzfft(CZ)-1,1) = 0.0d0
!
        dbs(CX,:,:,uzfft(CX)+1:nz       +1,1) = 0.0d0
        dbs(CY,:,:,uzfft(CY)+1:nz       +1,1) = 0.0d0
        dbs(CZ,:,:,uzfft(CZ)+1:nz       +1,1) = 0.0d0
!
        if(myid.eq.0) then
          dbs(CX,:,         -1    :lyfft(CX)-1,:,1) = 0.0d0
          dbs(CY,:,         -1    :lyfft(CY)-1,:,1) = 0.0d0
          dbs(CZ,:,         -1    :lyfft(CZ)-1,:,1) = 0.0d0
        end if
        if(myid.eq.snode-1) then
          dbs(CX,:,uyfft(CX)+1-syl:ndslc    +1,:,1) = 0.0d0
          dbs(CY,:,uyfft(CY)+1-syl:ndslc    +1,:,1) = 0.0d0
          dbs(CZ,:,uyfft(CZ)+1-syl:ndslc    +1,:,1) = 0.0d0
        end if
!
        if(nfbnd(1).eq.0) then
          dbs(CX:CZ,  -1:  -1,:,:,1) = dbs(CX:CZ,nx-1:nx-1,:,:,1)
          dbs(CX:CZ,nx  :nx  ,:,:,1) = dbs(CX:CZ,   0:   0,:,:,1)
          dbs(CX:CZ,nx+1:nx+1,:,:,1) = dbs(CX:CZ,  +1:  +1,:,:,1)
        end if
!
        if(nfbnd(3).eq.0) then
          dbs(CX:CZ,:,:,  -1:  -1,1) = dbs(CX:CZ,:,:,nz-1:nz-1,1)
          dbs(CX:CZ,:,:,nz  :nz  ,1) = dbs(CX:CZ,:,:,   0:   0,1)
          dbs(CX:CZ,:,:,nz+1:nz+1,1) = dbs(CX:CZ,:,:,  +1:  +1,1)
        end if
!
        call MPI_sendrecv(dbs(1,-1,syu-syl-1,-1,1),1, &
       &                  mptype_xz(1),uslice(1),stag, &
       &                  dbs(1,-1,-1,-1,1),1, &
       &                  mptype_xz(1),lslice(1),rtag, &
       &                  subcomm, MPI_STATUS_IGNORE, mpierr)
        call MPI_sendrecv(dbs(1,-1,0,-1,1),1, &
       &                  mptype_xz(1),lslice(1),stag, &
       &                  dbs(1,-1,syu-syl,-1,1),1, &
       &                  mptype_xz(1),uslice(1),rtag, &
       &                  subcomm, MPI_STATUS_IGNORE, mpierr)
!
!     --------------- 
      else if(pftmode.eq.3) then
        dbs(CX,:,         -1:lyfft(CX)-1,:,1) = 0.0d0
        dbs(CY,:,         -1:lyfft(CY)-1,:,1) = 0.0d0
        dbs(CZ,:,         -1:lyfft(CZ)-1,:,1) = 0.0d0
!
        dbs(CX,:,uyfft(CX)+1:ny       +1,:,1) = 0.0d0
        dbs(CY,:,uyfft(CY)+1:ny       +1,:,1) = 0.0d0
        dbs(CZ,:,uyfft(CZ)+1:ny       +1,:,1) = 0.0d0
!
        dbs(CX,:,:,         -1:lzfft(CX)-1,1) = 0.0d0
        dbs(CY,:,:,         -1:lzfft(CY)-1,1) = 0.0d0
        dbs(CZ,:,:,         -1:lzfft(CZ)-1,1) = 0.0d0
!
        dbs(CX,:,:,uzfft(CX)+1:nz       +1,1) = 0.0d0
        dbs(CY,:,:,uzfft(CY)+1:nz       +1,1) = 0.0d0
        dbs(CZ,:,:,uzfft(CZ)+1:nz       +1,1) = 0.0d0
!
        if(myid.eq.0) then
          dbs(CX,         -1:lxfft(CX)-1,:,:,1) = 0.0d0
          dbs(CY,         -1:lxfft(CY)-1,:,:,1) = 0.0d0
          dbs(CZ,         -1:lxfft(CZ)-1,:,:,1) = 0.0d0
        end if
        if(myid.eq.snode-1) then
          dbs(CX,uxfft(CX)+1:nwslc    +1,:,:,1) = 0.0d0
          dbs(CY,uxfft(CY)+1:nwslc    +1,:,:,1) = 0.0d0
          dbs(CZ,uxfft(CZ)+1:nwslc    +1,:,:,1) = 0.0d0
        end if
!
        if(nfbnd(1).eq.0) then
          if(myid.eq.0) &
         &  dbs(CX:CZ,  -1:  -1,:,:,1) = dbs(CX:CZ,nx-1:nx-1,:,:,1)
          if(myid.eq.snode-1) then
            dbs(CX:CZ,nx  :nx  ,:,:,1) = dbs(CX:CZ,   0:   0,:,:,1)
            dbs(CX:CZ,nx+1:nx+1,:,:,1) = dbs(CX:CZ,  +1:  +1,:,:,1)
          end if
        end if
!
        if(nfbnd(2).eq.0) then
          dbs(CX:CZ,:,  -1:  -1,:,1) = dbs(CX:CZ,:,ny-1:ny-1,:,1)
          dbs(CX:CZ,:,ny  :ny  ,:,1) = dbs(CX:CZ,:,   0:   0,:,1)
          dbs(CX:CZ,:,ny+1:ny+1,:,1) = dbs(CX:CZ,:,  +1:  +1,:,1)
        end if
!
        if(nfbnd(3).eq.0) then
          dbs(CX:CZ,:,:,  -1:  -1,1) = dbs(CX:CZ,:,:,nz-1:nz-1,1)
          dbs(CX:CZ,:,:,nz  :nz  ,1) = dbs(CX:CZ,:,:,   0:   0,1)
          dbs(CX:CZ,:,:,nz+1:nz+1,1) = dbs(CX:CZ,:,:,  +1:  +1,1)
        end if
!
!     --------------- 
      else if(pftmode.eq.4) then
        dbs(CX,         -1:lxfft(CX)-1,:,:,1) = 0.0d0
        dbs(CY,         -1:lxfft(CY)-1,:,:,1) = 0.0d0
        dbs(CZ,         -1:lxfft(CZ)-1,:,:,1) = 0.0d0
!
        dbs(CX,uxfft(CX)+1:nx       +1,:,:,1) = 0.0d0
        dbs(CY,uxfft(CY)+1:nx       +1,:,:,1) = 0.0d0
        dbs(CZ,uxfft(CZ)+1:nx       +1,:,:,1) = 0.0d0
!
        dbs(CX,:,         -1:lyfft(CX)-1,:,1) = 0.0d0
        dbs(CY,:,         -1:lyfft(CY)-1,:,1) = 0.0d0
        dbs(CZ,:,         -1:lyfft(CZ)-1,:,1) = 0.0d0
!
        dbs(CX,:,uyfft(CX)+1:ny       +1,:,1) = 0.0d0
        dbs(CY,:,uyfft(CY)+1:ny       +1,:,1) = 0.0d0
        dbs(CZ,:,uyfft(CZ)+1:ny       +1,:,1) = 0.0d0
!
        if(myid.eq.0) then
          dbs(CX,:,:,         -1:lzfft(CX)-1,1) = 0.0d0
          dbs(CY,:,:,         -1:lzfft(CY)-1,1) = 0.0d0
          dbs(CZ,:,:,         -1:lzfft(CZ)-1,1) = 0.0d0
        end if
        if(myid.eq.snode-1) then
          dbs(CX,:,:,uzfft(CX)+1:nhslc    +1,1) = 0.0d0
          dbs(CY,:,:,uzfft(CY)+1:nhslc    +1,1) = 0.0d0
          dbs(CZ,:,:,uzfft(CZ)+1:nhslc    +1,1) = 0.0d0
        end if
!
        if(nfbnd(1).eq.0) then
          dbs(CX:CZ,  -1:  -1,:,:,1) = dbs(CX:CZ,nx-1:nx-1,:,:,1)
          dbs(CX:CZ,nx  :nx  ,:,:,1) = dbs(CX:CZ,   0:   0,:,:,1)
          dbs(CX:CZ,nx+1:nx+1,:,:,1) = dbs(CX:CZ,  +1:  +1,:,:,1)
        end if
!
        if(nfbnd(2).eq.0) then
          dbs(CX:CZ,:,  -1:  -1,:,1) = dbs(CX:CZ,:,ny-1:ny-1,:,1)
          dbs(CX:CZ,:,ny  :ny  ,:,1) = dbs(CX:CZ,:,   0:   0,:,1)
          dbs(CX:CZ,:,ny+1:ny+1,:,1) = dbs(CX:CZ,:,  +1:  +1,:,1)
        end if
!
        if(nfbnd(3).eq.0) then
          if(myid.eq.0) &
         &  dbs(CX:CZ,:,:,  -1:  -1,1) = dbs(CX:CZ,:,:,nz-1:nz-1,1)
          if(myid.eq.snode-1) then
            dbs(CX:CZ,:,:,nz  :nz  ,1) = dbs(CX:CZ,:,:,   0:   0,1)
            dbs(CX:CZ,:,:,nz+1:nz+1,1) = dbs(CX:CZ,:,:,  +1:  +1,1)
          end if
        end if
      end if
      end if


!-------------------- 
      if(pftmode.eq.1) then
        db(:,:,:,:,1) = 0.0d0
        icomreq = 0
        do ps=1,2
          if(sdid(ps).eq.-1) cycle
          xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
          yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
          zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)
          do src=int(zl/nzslc),int((zu-1)/nzslc)
            count = 0 
            do islice=src*nzslc,(src+1)*nzslc-1
              if(zl.le.islice.and.islice.lt.zu) then
                count = count + 1
              end if
            end do
!            if(count.ne.nzslc) print*, "nzslc1", count, nzslc
            ztype(1) = 0
            if(src.eq.int(zl/nzslc)) then
              count = count + 1
              ztype(1) = ztype(1) - 1
            end if
            if(src.eq.int((zu-1)/nzslc)) then
              count = count + 2
            end if
            icomreq = icomreq + 1
            call MPI_Irecv(db(1,0,0,max(src*nzslc-zl,0),1), &
           &               count, mptype_rc(ztype(1),1,1), &
           &               src, 0, MCW, ireqs(icomreq), mpierr)
          end do
        end do
!
        if(myid.lt.snode) then
          lfloor = int(szl/nzsd); ufloor = int((szu-1)/nzsd)
          do ifloor=lfloor,ufloor
            count = min(szu,(ifloor+1)*nzsd) - max(szl,ifloor*nzsd)
            if(count.ne.nzslc) print*, "nzslc2", count, nzslc
            ztype(1) = 0
            if(mod(max(szl,ifloor*nzsd),nzsd).eq.0) then
              count = count + 1
              ztype(1) = ztype(1) - 1
            end if
            if(mod(min(szu,(ifloor+1)*nzsd),nzsd).eq.0) then
              count = count + 2
            end if
            do j=0,nodes(2)-1
            do i=0,nodes(1)-1
              sd = i + j*nodes(1) + ifloor*nodes(1)*nodes(2)
              from = famind(sd+1) + 1
              to = famind(sd+2)
              do ifm=from,to
                dst = fammbr(ifm)
                icomreq = icomreq + 1
                call MPI_Isend(dbs(1,i*nxsd,j*nysd,max(ifloor*nzsd-szl,0),1), &
               &               count, mptype_ss(ztype(1),1,1), &
               &               dst, 0, MCW, ireqs(icomreq), mpierr)
              end do
            end do
            end do
          end do
        end if
!
        call MPI_Waitall(icomreq,ireqs,istatus,mpierr)
!
!     --------------- 
      else if(pftmode.eq.2) then
        db(:,:,:,:,1) = 0.0d0
        icomreq = 0
        do ps=1,2
          if(sdid(ps).eq.-1) cycle
          xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
          yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
          zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)
          do src=int(yl/nyslc),int((yu-1)/nyslc)
            count = 0 
            do islice=src*nyslc,(src+1)*nyslc-1
              if(yl.le.islice.and.islice.lt.yu) then
                count = count + 1
              end if
            end do
!            if(count.ne.nyslc) print*, "nyslc1", count, nyslc
            ytype(1) = 0
            if(src.eq.int(yl/nyslc)) then
              count = count + 1
              ytype(1) = ytype(1) - 1
            end if
            if(src.eq.int((yu-1)/nyslc)) then
              count = count + 2
            end if
            icomreq = icomreq + 1
            call MPI_Irecv(db(1,0,max(src*nyslc-yl,0),0,1), &
           &               count, mptype_rc(ytype(1),1,2), &
           &               src, 0, MCW, ireqs(icomreq), mpierr)
          end do
        end do
!
        if(myid.lt.snode) then
          lfloor = int(syl/nysd); ufloor = int((syu-1)/nysd)
          do ifloor=lfloor,ufloor
            count = min(syu,(ifloor+1)*nysd) - max(syl,ifloor*nysd)
            if(count.ne.nyslc) print*, "nyslc2", count, nyslc
            ytype(1) = 0
            if(mod(max(syl,ifloor*nysd),nysd).eq.0) then
              count = count + 1
              ytype(1) = ytype(1) - 1
            end if
            if(mod(min(syu,(ifloor+1)*nysd),nysd).eq.0) then
              count = count + 2
            end if
            do k=0,nodes(3)-1
            do i=0,nodes(1)-1
              sd = i + ifloor*nodes(1) + k*nodes(1)*nodes(2)
              from = famind(sd+1) + 1
              to = famind(sd+2)
              do ifm=from,to
                dst = fammbr(ifm)
                icomreq = icomreq + 1
                call MPI_Isend(dbs(1,i*nxsd,max(ifloor*nysd-syl,0),k*nzsd,1), &
               &               count, mptype_ss(ytype(1),1,2), &
               &               dst, 0, MCW, ireqs(icomreq), mpierr)
              end do
            end do
            end do
          end do
        end if
!
        call MPI_Waitall(icomreq,ireqs,istatus,mpierr)
!
!     --------------- 
      else if(pftmode.eq.3) then
        db(CX:CZ,-1:xu-xl+1,-1:ny+1,-1:nz+1,1) = &
       &  dbs(CX:CZ,xl-1:xu+1,-1:ny+1,-1:nz+1,1)
!
!     --------------- 
      else if(pftmode.eq.4) then
        db(CX:CZ,-1:nx+1,-1:ny+1,-1:zu-zl+1,1) = &
       &  dbs(CX:CZ,-1:nx+1,-1:ny+1,zl-1:zu+1,1)
      end if


  return
  end subroutine ibfield1


!
  subroutine ibfield2(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   I B F I E L D 2
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   subroutine for update of b-field for full time step    .
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


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- loop for update of the b-field
      do k=0,zu-zl
      do j=0,yu-yl
      do i=0,xu-xl
        eb(BX,i,j,k,ps) = eb(BX,i,j,k,ps) + db(CX,i,j,k,ps)
        eb(BY,i,j,k,ps) = eb(BY,i,j,k,ps) + db(CY,i,j,k,ps)
        eb(BZ,i,j,k,ps) = eb(BZ,i,j,k,ps) + db(CZ,i,j,k,ps)
      end do
      end do
      end do


!-------------------- boundary treatment
!      call fsmask(1)
!      call fbound(1)


  return
  end subroutine ibfield2


!
  subroutine iefield(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   I E F I E L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   subroutine for update of e-field for full time step    .
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
  real(kind=8) :: thet
  real(kind=8),parameter :: dlt=2.0d0, ieps0=1.0d0
  real(kind=8) :: capc1


!-------------------- 
      thet = gfactor


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- loop for update of the e-field
      do k=0,zu-zl
      do j=0,yu-yl
      do i=0,xu-xl
        eb(EX,i,j,k,ps) = eb(EX,i,j,k,ps) &
       &                + thet*cs*dlt*( db(CZ,i,j,k,ps) - db(CZ,i,j-1,k,ps) &
       &                              - db(CY,i,j,k,ps) + db(CY,i,j,k-1,ps) ) &
       &                + cs*dlt*( eb(BZ,i,j,k,ps) - eb(BZ,i,j-1,k,ps) &
       &                         - eb(BY,i,j,k,ps) + eb(BY,i,j,k-1,ps) ) &
       &                - dlt*ieps0*aj(JX,i,j,k,ps)
        eb(EY,i,j,k,ps) = eb(EY,i,j,k,ps) &
       &                + thet*cs*dlt*( db(CX,i,j,k,ps) - db(CX,i,j,k-1,ps) &
       &                              - db(CZ,i,j,k,ps) + db(CZ,i-1,j,k,ps) ) &
       &                + cs*dlt*( eb(BX,i,j,k,ps) - eb(BX,i,j,k-1,ps) &
       &                         - eb(BZ,i,j,k,ps) + eb(BZ,i-1,j,k,ps) ) &
       &                - dlt*ieps0*aj(JY,i,j,k,ps)
        eb(EZ,i,j,k,ps) = eb(EZ,i,j,k,ps) &
       &                + thet*cs*dlt*( db(CY,i,j,k,ps) - db(CY,i-1,j,k,ps) &
       &                              - db(CX,i,j,k,ps) + db(CX,i,j-1,k,ps) ) &
       &                + cs*dlt*( eb(BY,i,j,k,ps) - eb(BY,i-1,j,k,ps) &
       &                         - eb(BX,i,j,k,ps) + eb(BX,i,j-1,k,ps) ) &
       &                - dlt*ieps0*aj(JZ,i,j,k,ps)
      end do
      end do
      end do


  return
  end subroutine iefield
