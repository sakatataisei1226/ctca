#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine inifft
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   I N I F F T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  use hdf
#define MCW local_comm
#define MSS MPI_STATUS_SIZE
  implicit none
  integer(kind=4) :: i,j,k, m


!-------------------- kmod(:,:,1): -(theta*cv*dt)^2*laplacian
      mfactor(1) = 1.0d0
!
      if(nfbnd(1).eq.0) then
        m = 0
        do i=0,int(nx/2)
          kmod(1,m,1) = 2.0d0*cfactor(1)*sin(PI*i*rnx)
          kmod(1,m,1) = kmod(1,m,1)*kmod(1,m,1)
          m = m + 1
        end do
        do i=int((nx+1)/2)-1,1,-1
          kmod(1,m,1) = 2.0d0*cfactor(1)*sin(PI*i*rnx)
          kmod(1,m,1) = kmod(1,m,1)*kmod(1,m,1)
          m = m + 1
        end do
        if(m.ne.nx.and.myid.eq.0) print*, "Warning@kmod11: nx,m", nx, m
        mfactor(1) = mfactor(1)*nx
        lxfft(CX:CZ) = 0; uxfft(CX:CZ) = nx - 1
        fftw_type(CX:CZ,1,1) = FFTW_R2HC
        fftw_type(CX:CZ,1,2) = FFTW_HC2R
      else if(nfbnd(1).ge.1) then
        m = 0
        do i=0,nx
          kmod(1,m,1) = 2.0d0*(1.0d0 - cos(i*PI/nx))*cfactor(1)*cfactor(1)
          m = m + 1
        end do
        if(m.ne.nx+1.and.myid.eq.0) print*, "Warning@kmod11: nx,m", nx, m
        mfactor(1) = mfactor(1)*nx*2.0d0
        lxfft(CX:CX) = 1; uxfft(CX:CX) = nx - 1
        lxfft(CY:CZ) = 0; uxfft(CY:CZ) = nx - 1
        fftw_type(CX:CX,1,1) = FFTW_RODFT00
        fftw_type(CX:CX,1,2) = FFTW_RODFT00
        fftw_type(CY:CZ,1,1) = FFTW_RODFT10 ! FFTW_REDFT10
        fftw_type(CY:CZ,1,2) = FFTW_RODFT01 ! FFTW_REDFT01
      end if
!
      if(nfbnd(2).eq.0) then
        m = 0
        do j=0,int(ny/2)
          kmod(2,m,1) = 2.0d0*cfactor(2)*sin(PI*j*rny)
          kmod(2,m,1) = kmod(2,m,1)*kmod(2,m,1)
          m = m + 1
        end do
        do j=int((ny+1)/2)-1,1,-1
          kmod(2,m,1) = 2.0d0*cfactor(2)*sin(PI*j*rny)
          kmod(2,m,1) = kmod(2,m,1)*kmod(2,m,1)
          m = m + 1
        end do
        if(m.ne.ny.and.myid.eq.0) print*, "Warning@kmod21: ny,m", ny, m
        mfactor(1) = mfactor(1)*ny
        lyfft(CX:CZ) = 0; uyfft(CX:CZ) = ny - 1
        fftw_type(CX:CZ,2,1) = FFTW_R2HC
        fftw_type(CX:CZ,2,2) = FFTW_HC2R
      else if(nfbnd(2).ge.1) then
        m = 0
        do j=0,ny
          kmod(2,m,1) = 2.0d0*(1.0d0 - cos(j*PI/ny))*cfactor(2)*cfactor(2)
          m = m + 1
        end do
        if(m.ne.ny+1.and.myid.eq.0) print*, "Warning@kmod21: ny,m", ny, m
        mfactor(1) = mfactor(1)*ny*2.0d0
        lyfft(CX) = 0; uyfft(CX) = ny - 1
        lyfft(CY) = 1; uyfft(CY) = ny - 1
        lyfft(CZ) = 0; uyfft(CZ) = ny - 1
        fftw_type(CX,2,1) = FFTW_RODFT10 ! FFTW_REDFT10
        fftw_type(CX,2,2) = FFTW_RODFT01 ! FFTW_REDFT01
        fftw_type(CY,2,1) = FFTW_RODFT00
        fftw_type(CY,2,2) = FFTW_RODFT00
        fftw_type(CZ,2,1) = FFTW_RODFT10 ! FFTW_REDFT10
        fftw_type(CZ,2,2) = FFTW_RODFT01 ! FFTW_REDFT01
      end if
!
      if(nfbnd(3).eq.0) then
        m = 0
        do k=0,int(nz/2)
          kmod(3,m,1) = 2.0d0*cfactor(3)*sin(PI*k*rnz)
          kmod(3,m,1) = kmod(3,m,1)*kmod(3,m,1)
          m = m + 1
        end do
        do k=int((nz+1)/2)-1,1,-1
          kmod(3,m,1) = 2.0d0*cfactor(3)*sin(PI*k*rnz)
          kmod(3,m,1) = kmod(3,m,1)*kmod(3,m,1)
          m = m + 1
        end do
        if(m.ne.nz.and.myid.eq.0) print*, "Warning@kmod31: nz,m", nz, m
        mfactor(1) = mfactor(1)*nz
        lzfft(CX:CZ) = 0; uzfft(CX:CZ) = nz - 1
        fftw_type(CX:CZ,3,1) = FFTW_R2HC
        fftw_type(CX:CZ,3,2) = FFTW_HC2R
      else if(nfbnd(3).eq.1) then
        m = 0
        do k=0,nz
          kmod(3,m,1) = 2.0d0*(1.0d0 - cos(k*PI/nz))*cfactor(3)*cfactor(3)
          m = m + 1
        end do
        if(m.ne.nz+1.and.myid.eq.0) print*, "Warning@kmod31: nz,m", nz, m
        mfactor(1) = mfactor(1)*nz*2.0d0
        lzfft(CX:CY) = 0; uzfft(CX:CY) = nz - 1
        lzfft(CZ:CZ) = 1; uzfft(CZ:CZ) = nz - 1
        fftw_type(CX:CY,3,1) = FFTW_RODFT10 ! FFTW_REDFT10
        fftw_type(CX:CY,3,2) = FFTW_RODFT01 ! FFTW_REDFT01
        fftw_type(CZ:CZ,3,1) = FFTW_RODFT00
        fftw_type(CZ:CZ,3,2) = FFTW_RODFT00
      end if
!
!     --------------- pftmode==1
      if(pftmode.eq.1) then
        call dfftw_plan_r2r_2d(fftplan(CX,1,1), &
       &  uxfft(CX)-lxfft(CX)+1,uyfft(CX)-lyfft(CX)+1, &
       &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1), &
       &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),2), &
       &  fftw_type(CX,1,1),fftw_type(CX,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CY,1,1), &
       &  uxfft(CY)-lxfft(CY)+1,uyfft(CY)-lyfft(CY)+1, &
       &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1), &
       &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),2), &
       &  fftw_type(CY,1,1),fftw_type(CY,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CZ,1,1), &
       &  uxfft(CZ)-lxfft(CZ)+1,uyfft(CZ)-lyfft(CZ)+1, &
       &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1), &
       &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),2), &
       &  fftw_type(CZ,1,1),fftw_type(CZ,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CX,1,2), &
       &  uxfft(CX)-lxfft(CX)+1,uyfft(CX)-lyfft(CX)+1, &
       &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1), &
       &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),2), &
       &  fftw_type(CX,1,2),fftw_type(CX,2,2),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CY,1,2), &
       &  uxfft(CY)-lxfft(CY)+1,uyfft(CY)-lyfft(CY)+1, &
       &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1), &
       &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),2), &
       &  fftw_type(CY,1,2),fftw_type(CY,2,2),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CZ,1,2), &
       &  uxfft(CZ)-lxfft(CZ)+1,uyfft(CZ)-lyfft(CZ)+1, &
       &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1), &
       &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),2), &
       &  fftw_type(CZ,1,2),fftw_type(CZ,2,2),FFTW_MEASURE)
!
        call dfftw_plan_r2r_1d(fftplan(CX,2,1), &
       &  uzfft(CX)-lzfft(CX)+1, &
       &  warray1d(lzfft(CX):uzfft(CX),1), &
       &  warray1d(lzfft(CX):uzfft(CX),2), &
       &  fftw_type(CX,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CY,2,1), &
       &  uzfft(CY)-lzfft(CY)+1, &
       &  warray1d(lzfft(CY):uzfft(CY),1), &
       &  warray1d(lzfft(CY):uzfft(CY),2), &
       &  fftw_type(CY,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CZ,2,1), &
       &  uzfft(CZ)-lzfft(CZ)+1, &
       &  warray1d(lzfft(CZ):uzfft(CZ),1), &
       &  warray1d(lzfft(CZ):uzfft(CZ),2), &
       &  fftw_type(CZ,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CX,2,2), &
       &  uzfft(CX)-lzfft(CX)+1, &
       &  warray1d(lzfft(CX):uzfft(CX),1), &
       &  warray1d(lzfft(CX):uzfft(CX),2), &
       &  fftw_type(CX,3,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CY,2,2), &
       &  uzfft(CY)-lzfft(CY)+1, &
       &  warray1d(lzfft(CY):uzfft(CY),1), &
       &  warray1d(lzfft(CY):uzfft(CY),2), &
       &  fftw_type(CY,3,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CZ,2,2), &
       &  uzfft(CZ)-lzfft(CZ)+1, &
       &  warray1d(lzfft(CZ):uzfft(CZ),1), &
       &  warray1d(lzfft(CZ):uzfft(CZ),2), &
       &  fftw_type(CZ,3,2),FFTW_MEASURE)
!
!     --------------- 
      else if(pftmode.eq.2) then
        call dfftw_plan_r2r_2d(fftplan(CX,1,1), &
       &  uxfft(CX)-lxfft(CX)+1,uzfft(CX)-lzfft(CX)+1, &
       &  warray2d(lxfft(CX):uxfft(CX),lzfft(CX):uzfft(CX),1), &
       &  warray2d(lxfft(CX):uxfft(CX),lzfft(CX):uzfft(CX),2), &
       &  fftw_type(CX,1,1),fftw_type(CX,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CY,1,1), &
       &  uxfft(CY)-lxfft(CY)+1,uzfft(CY)-lzfft(CY)+1, &
       &  warray2d(lxfft(CY):uxfft(CY),lzfft(CY):uzfft(CY),1), &
       &  warray2d(lxfft(CY):uxfft(CY),lzfft(CY):uzfft(CY),2), &
       &  fftw_type(CY,1,1),fftw_type(CY,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CZ,1,1), &
       &  uxfft(CZ)-lxfft(CZ)+1,uzfft(CZ)-lzfft(CZ)+1, &
       &  warray2d(lxfft(CZ):uxfft(CZ),lzfft(CZ):uzfft(CZ),1), &
       &  warray2d(lxfft(CZ):uxfft(CZ),lzfft(CZ):uzfft(CZ),2), &
       &  fftw_type(CZ,1,1),fftw_type(CZ,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CX,1,2), &
       &  uxfft(CX)-lxfft(CX)+1,uzfft(CX)-lzfft(CX)+1, &
       &  warray2d(lxfft(CX):uxfft(CX),lzfft(CX):uzfft(CX),1), &
       &  warray2d(lxfft(CX):uxfft(CX),lzfft(CX):uzfft(CX),2), &
       &  fftw_type(CX,1,2),fftw_type(CX,3,2),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CY,1,2), &
       &  uxfft(CY)-lxfft(CY)+1,uzfft(CY)-lzfft(CY)+1, &
       &  warray2d(lxfft(CY):uxfft(CY),lzfft(CY):uzfft(CY),1), &
       &  warray2d(lxfft(CY):uxfft(CY),lzfft(CY):uzfft(CY),2), &
       &  fftw_type(CY,1,2),fftw_type(CY,3,2),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CZ,1,2), &
       &  uxfft(CZ)-lxfft(CZ)+1,uzfft(CZ)-lzfft(CZ)+1, &
       &  warray2d(lxfft(CZ):uxfft(CZ),lzfft(CZ):uzfft(CZ),1), &
       &  warray2d(lxfft(CZ):uxfft(CZ),lzfft(CZ):uzfft(CZ),2), &
       &  fftw_type(CZ,1,2),fftw_type(CZ,3,2),FFTW_MEASURE)
!
        call dfftw_plan_r2r_1d(fftplan(CX,2,1), &
       &  uyfft(CX)-lyfft(CX)+1, &
       &  warray1d(lyfft(CX):uyfft(CX),1), &
       &  warray1d(lyfft(CX):uyfft(CX),2), &
       &  fftw_type(CX,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CY,2,1), &
       &  uyfft(CY)-lyfft(CY)+1, &
       &  warray1d(lyfft(CY):uyfft(CY),1), &
       &  warray1d(lyfft(CY):uyfft(CY),2), &
       &  fftw_type(CY,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CZ,2,1), &
       &  uyfft(CZ)-lyfft(CZ)+1, &
       &  warray1d(lyfft(CZ):uyfft(CZ),1), &
       &  warray1d(lyfft(CZ):uyfft(CZ),2), &
       &  fftw_type(CZ,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CX,2,2), &
       &  uyfft(CX)-lyfft(CX)+1, &
       &  warray1d(lyfft(CX):uyfft(CX),1), &
       &  warray1d(lyfft(CX):uyfft(CX),2), &
       &  fftw_type(CX,2,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CY,2,2), &
       &  uyfft(CY)-lyfft(CY)+1, &
       &  warray1d(lyfft(CY):uyfft(CY),1), &
       &  warray1d(lyfft(CY):uyfft(CY),2), &
       &  fftw_type(CY,2,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CZ,2,2), &
       &  uyfft(CZ)-lyfft(CZ)+1, &
       &  warray1d(lyfft(CZ):uyfft(CZ),1), &
       &  warray1d(lyfft(CZ):uyfft(CZ),2), &
       &  fftw_type(CZ,2,2),FFTW_MEASURE)
!
!     --------------- pftmode==3
      else if(pftmode.eq.3) then
        call dfftw_plan_r2r_2d(fftplan(CX,1,1), &
       &  uyfft(CX)-lyfft(CX)+1,uzfft(CX)-lzfft(CX)+1, &
       &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1), &
       &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),2), &
       &  fftw_type(CX,2,1),fftw_type(CX,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CY,1,1), &
       &  uyfft(CY)-lyfft(CY)+1,uzfft(CY)-lzfft(CY)+1, &
       &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1), &
       &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),2), &
       &  fftw_type(CY,2,1),fftw_type(CY,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CZ,1,1), &
       &  uyfft(CZ)-lyfft(CZ)+1,uzfft(CZ)-lzfft(CZ)+1, &
       &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1), &
       &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),2), &
       &  fftw_type(CZ,2,1),fftw_type(CZ,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CX,1,2), &
       &  uyfft(CX)-lyfft(CX)+1,uzfft(CX)-lzfft(CX)+1, &
       &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),1), &
       &  warray2d(lyfft(CX):uyfft(CX),lzfft(CX):uzfft(CX),2), &
       &  fftw_type(CX,2,2),fftw_type(CX,3,2),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CY,1,2), &
       &  uyfft(CY)-lyfft(CY)+1,uzfft(CY)-lzfft(CY)+1, &
       &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),1), &
       &  warray2d(lyfft(CY):uyfft(CY),lzfft(CY):uzfft(CY),2), &
       &  fftw_type(CY,2,2),fftw_type(CY,3,2),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CZ,1,2), &
       &  uyfft(CZ)-lyfft(CZ)+1,uzfft(CZ)-lzfft(CZ)+1, &
       &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),1), &
       &  warray2d(lyfft(CZ):uyfft(CZ),lzfft(CZ):uzfft(CZ),2), &
       &  fftw_type(CZ,2,2),fftw_type(CZ,3,2),FFTW_MEASURE)
!
        call dfftw_plan_r2r_1d(fftplan(CX,2,1), &
       &  uxfft(CX)-lxfft(CX)+1, &
       &  warray1d(lxfft(CX):uxfft(CX),1), &
       &  warray1d(lxfft(CX):uxfft(CX),2), &
       &  fftw_type(CX,1,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CY,2,1), &
       &  uxfft(CY)-lxfft(CY)+1, &
       &  warray1d(lxfft(CY):uxfft(CY),1), &
       &  warray1d(lxfft(CY):uxfft(CY),2), &
       &  fftw_type(CY,1,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CZ,2,1), &
       &  uxfft(CZ)-lxfft(CZ)+1, &
       &  warray1d(lxfft(CZ):uxfft(CZ),1), &
       &  warray1d(lxfft(CZ):uxfft(CZ),2), &
       &  fftw_type(CZ,1,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CX,2,2), &
       &  uxfft(CX)-lxfft(CX)+1, &
       &  warray1d(lxfft(CX):uxfft(CX),1), &
       &  warray1d(lxfft(CX):uxfft(CX),2), &
       &  fftw_type(CX,1,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CY,2,2), &
       &  uxfft(CY)-lxfft(CY)+1, &
       &  warray1d(lxfft(CY):uxfft(CY),1), &
       &  warray1d(lxfft(CY):uxfft(CY),2), &
       &  fftw_type(CY,1,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CZ,2,2), &
       &  uxfft(CZ)-lxfft(CZ)+1, &
       &  warray1d(lxfft(CZ):uxfft(CZ),1), &
       &  warray1d(lxfft(CZ):uxfft(CZ),2), &
       &  fftw_type(CZ,1,2),FFTW_MEASURE)
!
!     --------------- pftmode==4
      else if(pftmode.eq.4) then
        call dfftw_plan_r2r_2d(fftplan(CX,1,1), &
       &  uxfft(CX)-lxfft(CX)+1,uyfft(CX)-lyfft(CX)+1, &
       &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1), &
       &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),2), &
       &  fftw_type(CX,1,1),fftw_type(CX,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CY,1,1), &
       &  uxfft(CY)-lxfft(CY)+1,uyfft(CY)-lyfft(CY)+1, &
       &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1), &
       &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),2), &
       &  fftw_type(CY,1,1),fftw_type(CY,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CZ,1,1), &
       &  uxfft(CZ)-lxfft(CZ)+1,uyfft(CZ)-lyfft(CZ)+1, &
       &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1), &
       &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),2), &
       &  fftw_type(CZ,1,1),fftw_type(CZ,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CX,1,2), &
       &  uxfft(CX)-lxfft(CX)+1,uyfft(CX)-lyfft(CX)+1, &
       &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),1), &
       &  warray2d(lxfft(CX):uxfft(CX),lyfft(CX):uyfft(CX),2), &
       &  fftw_type(CX,1,2),fftw_type(CX,2,2),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CY,1,2), &
       &  uxfft(CY)-lxfft(CY)+1,uyfft(CY)-lyfft(CY)+1, &
       &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),1), &
       &  warray2d(lxfft(CY):uxfft(CY),lyfft(CY):uyfft(CY),2), &
       &  fftw_type(CY,1,2),fftw_type(CY,2,2),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CZ,1,2), &
       &  uxfft(CZ)-lxfft(CZ)+1,uyfft(CZ)-lyfft(CZ)+1, &
       &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),1), &
       &  warray2d(lxfft(CZ):uxfft(CZ),lyfft(CZ):uyfft(CZ),2), &
       &  fftw_type(CZ,1,2),fftw_type(CZ,2,2),FFTW_MEASURE)
!
        call dfftw_plan_r2r_1d(fftplan(CX,2,1), &
       &  uzfft(CX)-lzfft(CX)+1, &
       &  warray1d(lzfft(CX):uzfft(CX),1), &
       &  warray1d(lzfft(CX):uzfft(CX),2), &
       &  fftw_type(CX,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CY,2,1), &
       &  uzfft(CY)-lzfft(CY)+1, &
       &  warray1d(lzfft(CY):uzfft(CY),1), &
       &  warray1d(lzfft(CY):uzfft(CY),2), &
       &  fftw_type(CY,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CZ,2,1), &
       &  uzfft(CZ)-lzfft(CZ)+1, &
       &  warray1d(lzfft(CZ):uzfft(CZ),1), &
       &  warray1d(lzfft(CZ):uzfft(CZ),2), &
       &  fftw_type(CZ,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CX,2,2), &
       &  uzfft(CX)-lzfft(CX)+1, &
       &  warray1d(lzfft(CX):uzfft(CX),1), &
       &  warray1d(lzfft(CX):uzfft(CX),2), &
       &  fftw_type(CX,3,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CY,2,2), &
       &  uzfft(CY)-lzfft(CY)+1, &
       &  warray1d(lzfft(CY):uzfft(CY),1), &
       &  warray1d(lzfft(CY):uzfft(CY),2), &
       &  fftw_type(CY,3,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CZ,2,2), &
       &  uzfft(CZ)-lzfft(CZ)+1, &
       &  warray1d(lzfft(CZ):uzfft(CZ),1), &
       &  warray1d(lzfft(CZ):uzfft(CZ),2), &
       &  fftw_type(CZ,3,2),FFTW_MEASURE)
      end if


!-------------------- kmod(:,:,2): -laplacian
      mfactor(2) = 1.0d0
!
      if(mtd_vbnd(1).eq.0) then
        m = 0
        do i=0,int(nx/2)
          kmod(1,m,2) = 2.0d0*sin(PI*i*rnx)*dri
          kmod(1,m,2) = kmod(1,m,2)*kmod(1,m,2)
          m = m + 1
        end do
        do i=int((nx+1)/2)-1,1,-1
          kmod(1,m,2) = 2.0d0*sin(PI*i*rnx)*dri
          kmod(1,m,2) = kmod(1,m,2)*kmod(1,m,2)
          m = m + 1
        end do
        if(m.ne.nx.and.myid.eq.0) print*, "Warning@kmod12: nx,m", nx, m
        mfactor(2) = mfactor(2)*nx
        lxfft(CP) = 0; uxfft(CP) = nx - 1
        fftw_type(CP,1,1) = FFTW_R2HC
        fftw_type(CP,1,2) = FFTW_HC2R
      else if(mtd_vbnd(1).ge.1) then
        m = 0
        do i=0,nx
          kmod(1,m,2) = 2.0d0*(1.0d0 - cos(i*PI/nx))*dri*dri
          m = m + 1
        end do
        if(m.ne.nx+1.and.myid.eq.0) print*, "Warning@kmod12: nx,m", nx, m
        mfactor(2) = mfactor(2)*nx*2.0d0
        if(mtd_vbnd(1).eq.1) then
          lxfft(CP) = 1; uxfft(CP) = nx - 1
          fftw_type(CP,1,1) = FFTW_RODFT00
          fftw_type(CP,1,2) = FFTW_RODFT00
        else if(mtd_vbnd(1).eq.2) then
          lxfft(CP) = 0; uxfft(CP) = nx
          fftw_type(CP,1,1) = FFTW_REDFT00
          fftw_type(CP,1,2) = FFTW_REDFT00
        end if
      end if
!
      if(mtd_vbnd(2).eq.0) then
        m = 0
        do j=0,int(ny/2)
          kmod(2,m,2) = 2.0d0*sin(PI*j*rny)*dri
          kmod(2,m,2) = kmod(2,m,2)*kmod(2,m,2)
          m = m + 1
        end do
        do j=int((ny+1)/2)-1,1,-1
          kmod(2,m,2) = 2.0d0*sin(PI*j*rny)*dri
          kmod(2,m,2) = kmod(2,m,2)*kmod(2,m,2)
          m = m + 1
        end do
        if(m.ne.ny.and.myid.eq.0) print*, "Warning@kmod22: ny,m", ny, m
        mfactor(2) = mfactor(2)*ny
        lyfft(CP) = 0; uyfft(CP) = ny - 1
        fftw_type(CP,2,1) = FFTW_R2HC
        fftw_type(CP,2,2) = FFTW_HC2R
      else if(mtd_vbnd(2).ge.1) then
        m = 0
        do j=0,ny
          kmod(2,m,2) = 2.0d0*(1.0d0 - cos(j*PI/ny))*dri*dri
          m = m + 1
        end do
        if(m.ne.ny+1.and.myid.eq.0) print*, "Warning@kmod22: ny,m", ny, m
        mfactor(2) = mfactor(2)*ny*2.0d0
        if(mtd_vbnd(2).eq.1) then
          lyfft(CP) = 1; uyfft(CP) = ny - 1
          fftw_type(CP,2,1) = FFTW_RODFT00
          fftw_type(CP,2,2) = FFTW_RODFT00
        else if(mtd_vbnd(2).eq.2) then
          lyfft(CP) = 0; uyfft(CP) = ny
          fftw_type(CP,2,1) = FFTW_REDFT00
          fftw_type(CP,2,2) = FFTW_REDFT00
        end if
      end if
!
      if(mtd_vbnd(3).eq.0) then
        m = 0
        do k=0,int(nz/2)
          kmod(3,m,2) = 2.0d0*sin(PI*k*rnz)*dri
          kmod(3,m,2) = kmod(3,m,2)*kmod(3,m,2)
          m = m + 1
        end do
        do k=int((nz+1)/2)-1,1,-1
          kmod(3,m,2) = 2.0d0*sin(PI*k*rnz)*dri
          kmod(3,m,2) = kmod(3,m,2)*kmod(3,m,2)
          m = m + 1
        end do
        if(m.ne.nz.and.myid.eq.0) print*, "Warning@kmod32: nz,m", nz, m
        mfactor(2) = mfactor(2)*nz
        lzfft(CP) = 0; uzfft(CP) = nz - 1
        fftw_type(CP,3,1) = FFTW_R2HC
        fftw_type(CP,3,2) = FFTW_HC2R
      else if(mtd_vbnd(3).ge.1) then
        m = 0
        do k=0,nz
          kmod(3,m,2) = 2.0d0*(1.0d0 - cos(k*PI/nz))*dri*dri
          m = m + 1
        end do
        if(m.ne.nz+1.and.myid.eq.0) print*, "Warning@kmod32: nz,m", nz, m
        mfactor(2) = mfactor(2)*nz*2.0d0
        if(mtd_vbnd(3).eq.1) then
          lzfft(CP) = 1; uzfft(CP) = nz - 1
          fftw_type(CP,3,1) = FFTW_RODFT00
          fftw_type(CP,3,2) = FFTW_RODFT00
        else if(mtd_vbnd(3).eq.2) then
          lzfft(CP) = 0; uzfft(CP) = nz
          fftw_type(CP,3,1) = FFTW_REDFT00
          fftw_type(CP,3,2) = FFTW_REDFT00
        end if
      end if
!
!     --------------- 
      if(pftmode.eq.1) then
        call dfftw_plan_r2r_2d(fftplan(CP,1,1), &
       &  uxfft(CP)-lxfft(CP)+1,uyfft(CP)-lyfft(CP)+1, &
       &  warray2d(lxfft(CP):uxfft(CP),lyfft(CP):uyfft(CP),1), &
       &  warray2d(lxfft(CP):uxfft(CP),lyfft(CP):uyfft(CP),2), &
       &  fftw_type(CP,1,1),fftw_type(CP,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CP,1,2), &
       &  uxfft(CP)-lxfft(CP)+1,uyfft(CP)-lyfft(CP)+1, &
       &  warray2d(lxfft(CP):uxfft(CP),lyfft(CP):uyfft(CP),1), &
       &  warray2d(lxfft(CP):uxfft(CP),lyfft(CP):uyfft(CP),2), &
       &  fftw_type(CP,1,2),fftw_type(CP,2,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CP,2,1), &
       &  uzfft(CP)-lzfft(CP)+1, &
       &  warray1d(lzfft(CP):uzfft(CP),1), &
       &  warray1d(lzfft(CP):uzfft(CP),2), &
       &  fftw_type(CP,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CP,2,2), &
       &  uzfft(CP)-lzfft(CP)+1, &
       &  warray1d(lzfft(CP):uzfft(CP),1), &
       &  warray1d(lzfft(CP):uzfft(CP),2), &
       &  fftw_type(CP,3,2),FFTW_MEASURE)
!
!     --------------- 
      else if(pftmode.eq.2) then
        call dfftw_plan_r2r_2d(fftplan(CP,1,1), &
       &  uxfft(CP)-lxfft(CP)+1,uzfft(CP)-lzfft(CP)+1, &
       &  warray2d(lxfft(CP):uxfft(CP),lzfft(CP):uzfft(CP),1), &
       &  warray2d(lxfft(CP):uxfft(CP),lzfft(CP):uzfft(CP),2), &
       &  fftw_type(CP,1,1),fftw_type(CP,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CP,1,2), &
       &  uxfft(CP)-lxfft(CP)+1,uzfft(CP)-lzfft(CP)+1, &
       &  warray2d(lxfft(CP):uxfft(CP),lzfft(CP):uzfft(CP),1), &
       &  warray2d(lxfft(CP):uxfft(CP),lzfft(CP):uzfft(CP),2), &
       &  fftw_type(CP,1,2),fftw_type(CP,3,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CP,2,1), &
       &  uyfft(CP)-lyfft(CP)+1, &
       &  warray1d(lyfft(CP):uyfft(CP),1), &
       &  warray1d(lyfft(CP):uyfft(CP),2), &
       &  fftw_type(CP,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CP,2,2), &
       &  uyfft(CP)-lyfft(CP)+1, &
       &  warray1d(lyfft(CP):uyfft(CP),1), &
       &  warray1d(lyfft(CP):uyfft(CP),2), &
       &  fftw_type(CP,2,2),FFTW_MEASURE)
!
!     --------------- 
      else if(pftmode.eq.3) then
        call dfftw_plan_r2r_2d(fftplan(CP,1,1), &
       &  uyfft(CP)-lyfft(CP)+1,uzfft(CP)-lzfft(CP)+1, &
       &  warray2d(lyfft(CP):uyfft(CP),lzfft(CP):uzfft(CP),1), &
       &  warray2d(lyfft(CP):uyfft(CP),lzfft(CP):uzfft(CP),2), &
       &  fftw_type(CP,2,1),fftw_type(CP,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CP,1,2), &
       &  uyfft(CP)-lyfft(CP)+1,uzfft(CP)-lzfft(CP)+1, &
       &  warray2d(lyfft(CP):uyfft(CP),lzfft(CP):uzfft(CP),1), &
       &  warray2d(lyfft(CP):uyfft(CP),lzfft(CP):uzfft(CP),2), &
       &  fftw_type(CP,2,2),fftw_type(CP,3,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CP,2,1), &
       &  uxfft(CP)-lxfft(CP)+1, &
       &  warray1d(lxfft(CP):uxfft(CP),1), &
       &  warray1d(lxfft(CP):uxfft(CP),2), &
       &  fftw_type(CP,1,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CP,2,2), &
       &  uxfft(CP)-lxfft(CP)+1, &
       &  warray1d(lxfft(CP):uxfft(CP),1), &
       &  warray1d(lxfft(CP):uxfft(CP),2), &
       &  fftw_type(CP,1,2),FFTW_MEASURE)
!
!     --------------- 
      else if(pftmode.eq.4) then
        call dfftw_plan_r2r_2d(fftplan(CP,1,1), &
       &  uxfft(CP)-lxfft(CP)+1,uyfft(CP)-lyfft(CP)+1, &
       &  warray2d(lxfft(CP):uxfft(CP),lyfft(CP):uyfft(CP),1), &
       &  warray2d(lxfft(CP):uxfft(CP),lyfft(CP):uyfft(CP),2), &
       &  fftw_type(CP,1,1),fftw_type(CP,2,1),FFTW_MEASURE)
        call dfftw_plan_r2r_2d(fftplan(CP,1,2), &
       &  uxfft(CP)-lxfft(CP)+1,uyfft(CP)-lyfft(CP)+1, &
       &  warray2d(lxfft(CP):uxfft(CP),lyfft(CP):uyfft(CP),1), &
       &  warray2d(lxfft(CP):uxfft(CP),lyfft(CP):uyfft(CP),2), &
       &  fftw_type(CP,1,2),fftw_type(CP,2,2),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CP,2,1), &
       &  uzfft(CP)-lzfft(CP)+1, &
       &  warray1d(lzfft(CP):uzfft(CP),1), &
       &  warray1d(lzfft(CP):uzfft(CP),2), &
       &  fftw_type(CP,3,1),FFTW_MEASURE)
        call dfftw_plan_r2r_1d(fftplan(CP,2,2), &
       &  uzfft(CP)-lzfft(CP)+1, &
       &  warray1d(lzfft(CP):uzfft(CP),1), &
       &  warray1d(lzfft(CP):uzfft(CP),2), &
       &  fftw_type(CP,3,2),FFTW_MEASURE)
      end if


!-------------------- 
      if(myid.eq.0) print*, "mfactor(1:2) =", mfactor(1:2)


  return
  end subroutine inifft
