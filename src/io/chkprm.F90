#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine chkprm
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   C H K P R M
!   ____________________________________________________________
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i, j, is
  real(kind=8) :: sss, fff, omega, dts, theta, damp, attn
  real(kind=8) :: qqq, qmin, qmax, totalc


!-------------------- dimension size check.
!      if(max(jx,jy,jz).gt.jxy) then
!        print*,'kempo: chkprm: ========== ERROR !!! =========='
!        print*,'Parameter IXY must be >= max(ix,iy,jz)'
!        print*,'Parameter ix,iy,iz,ixy',jx,jy,jz,jxy
!        call kestop("CHKPRM:1  ",0) 
!      end if
!      if(nx+3.gt.jx) then
!        print*,'kempo: chkprm: ========== ERROR !!! =========='
!        print*,'Parameter IX must be >= (nx+3)'
!        print*,'Parameter ix,nx',ix,nx
!        call kestop("CHKPRM:2  ",0)
!      end if
!      if(ny+3.gt.jy) then
!        print*,'kempo: chkprm: ========== ERROR !!! =========='
!        print*,'Parameter IY must be >= (ny+3)'
!        print*,'Parameter iy,ny',iy,ny
!        call kestop("CHKPRM:3  ",0)
!      end if
!      if(nz+3.gt.jz) then
!        print*,'kempo: chkprm: ========== ERROR !!! =========='
!        print*,'Parameter IZ must be >= (nz+3)'
!        print*,'Parameter iz,nz',iz,nz
!        call kestop("CHKPRM:4  ",0)
!      end if
!      if(npsum.gt.in) then
!        print*,'kempo: chkprm: ========== ERROR !!! =========='
!        print*,'Parameter IN must be >= np(1)+np(2)+np(3)...'
!        call kestop("CHKPRM:5  ",0)
!      end if
!      if(nprsum.gt.inr) then
!        print*,'kempo: chkprm: ========== ERROR !!! =========='
!        print*, 'Parameter INR must be >= npr(1)+npr(2)+npr(3)+...''
!        call kestop("CHKPRM:1  ",0)
!      end if
!-------------------- diagnotics frequency check
      do is=nspec+1,ispec
        imdig(is) = 0
        ipadig(is) = 0
        ildig(is) = 0
      end do
!      if(ifdiag.gt.0.and.ifdiag.lt.mltstpf) ifdiag = mltstpf
      if(iediag.gt.0.and.iediag.lt.mltstpf) iediag = mltstpf
      if(isdiag.gt.0.and.isdiag.lt.mltstpf) isdiag = mltstpf
      if(iddiag.gt.0.and.iddiag.lt.mltstpf) iddiag = mltstpf
      do is=1,nspec
        if(ipadig(is).gt.0.and.ipadig(is).lt.mltstpf) ipadig(is) = mltstpf
!        if(imdig(is).gt.0.and.imdig(is).lt.mltstpf) imdig(is) = mltstpf
!        if(ildig(is).gt.0.and.ildig(is).lt.mltstpf) ildig(is) = mltstpf
      end do
!      if(ijdiag.gt.0.and.ijdiag.lt.mltstpf) ijdiag = mltstpf
!      if(iadiag.gt.0.and.iadiag.lt.mltstpf) iadiag = mltstpf
      if(itchck.gt.0.and.itchck.lt.mltstpf) itchck = mltstpf
      if(ivdiag.gt.0.and.ivdiag.lt.mltstpf) ivdiag = mltstpf
!-------------------- multiple step count check
      if(mod(nstep,mltstpf).ne.0) then
        print*,'kempo: chkprm:1 ----- WARNING! -----'
        print*,'NSTEP is not consistent with MLTSTPf'
        nstep = int(nstep/mltstpf)*mltstpf
        print*,'NSTEP has been changed.  NSTEP = ',nstep
      end if
!-------------------- courant condition check
      sss = dsqrt(3.0d0)
      if(nx.eq.1.or.ny.eq.1) sss = 1.0d0
      fff = sss*cv/renv*dt/dr
      if(juncan.ge.1000) then
        print*,'****** Test particle simulation ******'
        if(mltstp.ne.1) then
          mltstp = 1
          print*,'=== MLTSTP must be 1 for test particle simulaton'
          print*,'=== MLTSTP is set equal to 1 at kempo:chkprm'
        end if
        if(mltstpf.ne.1) then
          mltstpf = 1
          print*,'=== MLTSTPf must be 1 for test particle simulaton'
          print*,'=== MLTSTPf is set equal to 1 at kempo:chkprm'
        end if
      else if(fff.gt.1.) then
        if(emflag.eq.1) then
          if(myid.eq.0) then
            print*,'kempo: chkprm:2---------- WARNING!!!! ----------'
            print*,'--- Courant condition is not satisfied.'
            print*,'--- Check the parameters  DT DR CV.'
            print*,'DT=',dt,'DR=',dr,'CV=',cv,'RENV=',renv
            print*,'--- execution terminated at chkprm '
          end if
          stop
!          call kestop("CHKPRM:6  ",0)
        else
          if(myid.eq.0) then
            print*,'kempo: chkprm:2---------- information ----------'
            print*,'--- Courant condition is not satisfied,'
            print*,'--- ... but execution continued because this is ES/IM sim.'
          end if
        end if
      end if
!-------------------- multiple time step mode check
      if(mltstp.gt.1) then
        print*,'********** Multiple time step mode **********'
        print*,'       MLTSTP = ',mltstp
        print*,'Frequency: attenuation rate: attenuation factor'
        print*,'          at t=',dt*nstep
        do i=1,10
          omega = wp(1)*0.2d0*i
          dts = dt
!          call mtsatn(mltstp,dt,nstep,omega,attn,damp)
          if(mod(mltstp,2).eq.0) then
            attn = 0.0d0
            do j=1,mltstp/2
              theta = omega*dt*(dble(j) - 0.5d0)
              attn = attn + dcos(theta)
            end do
            attn = attn*2.0d0/dble(mltstp)
          else
            attn = 1.0d0
            do j=1,(mltstp-1)/2
              theta = omega*dt*dble(j)
              attn = attn + 2.0d0*dcos(theta)
            end do
            attn = attn/dble(mltstp)
          end if
          damp = attn**(nstep/mltstp)
          attn = dlog(attn)/(dt*mltstp)
          write(6,333) omega, attn, damp
  333     format(2x,1p,e13.3e3,8x,e13.3e3,8x,e13.3e3)
        end do
      end if
      if(mltstpf.gt.1) then
        print*,'********** Multiple time step mode **********'
        print*,'       MLTSTPf = ',mltstpf
        print*,'Frequency: attenuation rate: attenuation factor'
        print*,'          at t=',dt*nstep
        do i=1,10
          omega = wp(1)*0.2*i
          dts = dt
!          call mtsatn(mltstpf,dt,nstep,omega,attn,damp)
          if(mod(mltstpf,2).eq.0) then
            attn = 0.0d0
            do j=1,mltstpf/2
              theta = omega*dt*(dble(j) - 0.5d0)
              attn = attn + dcos(theta)
            end do
            attn = attn*2.0d0/dble(mltstpf)
          else
            attn = 1.0d0
            do j=1,(mltstpf-1)/2
              theta = omega*dt*dble(j)
              attn = attn + 2.0d0*dcos(theta)
            end do
            attn = attn/dble(mltstp)
          end if
          damp = attn**(nstep/mltstp)
          attn = dlog(attn)/(dt*mltstp)
          write(6,444) omega, attn, damp
  444     format(2x,1p,e13.3e3,8x,e13.3e3,8x,e13.3e3)
        end do
      end if
!-------------------- check of charge neutrality
!      if(nspec.eq.1) qqq = 1.0d0
!      if(nspec.ge.2) then
!        qmin = +1.0d10
!        qmax = -1.0d10
!        do is=1,nspec
!          if(qm(is).gt.qmax) qmax = qm(is)
!          if(qm(is).lt.qmin) qmin = qm(is)
!        end do
!        qqq = qmax*qmin
!      end if
!      if(qqq.lt.0.0d0) then
!        totalc = 0.0d0
!        do is=1,nspec
!          totalc = totalc + wp(is)**2/qm(is)
!        end do
!        if(totalc.ge.1.0d-3) then
!          print *,'kempo: chkprm:3---------- WARNING!!! ----------'
!          print *,'--- Charge neutrality is not satisfied.'
!          print *,'--- Total Charge:',totalc
!          print *,'--- Check the parameters,wp(nspec),qm(nspec) etc.'
!          call kestop("CHKPRM:7  ",0)
!        end if
!      end if
!-------------------- check of particle initialization
!      do is=1,nspec
!        IF(MOD(np(is),nphi(is)*ndst(is)).NE.0) THEN
!          print*,'np(is),nphi(is)*ndst(is)',np(is),nphi(is)*ndst(is)
!          print*,'MOD(np(is),nphi(is)*ndst(is))', &
!     &            MOD(np(is),nphi(is)*ndst(is))
!          print*,'kempo: chkprm:4---------- WARNING!!! ----------'
!          print*,'Species',is,': NPHI*NDST must be a factor of np'
!        END IF
!      end do


  return
  end subroutine
