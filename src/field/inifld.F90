#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine inifld
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   I N I F L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   this subroutine gives an initial setting of fields     .
!   ............................................................

!-------------------- parameter common blocks
  use oh_type
  use paramt
  use allcom
  use hdf
  implicit none
!
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: iw
!
  integer(kind=HID_T) :: fileid
  integer(kind=4) :: stats0,stats1
  integer(kind=8) :: dims(3)
  real(kind=8) :: x,y,z
  real(kind=8) :: orge(3), orgb(3)
  real(kind=8) :: xlmcol,xumcol,ylmcol,yumcol,zlmcol,zumcol
  real(kind=8) :: xlsgain,xusgain,ylsgain,yusgain,zlsgain,zusgain
  real(kind=8) :: sigmxl,sigmxu,sigmyl,sigmyu,sigmzl,sigmzu
  character(len=30) :: filename,dsname


!-------------------- 
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)


!-------------------- zero clear
      eb(:,:,:,:,:) = 0.0d0
      ebav(:,:,:,:) = 0.0d0
      mp(:,:,:,:,:) = 0.0d0
      aj(:,:,:,:,:) = 0.0d0
      rho(:,:,:,:,:) = 0.0d0
      phi(:,:,:,:,:) = 0.0d0
      phiav(:,:,:,:) = 0.0d0
      wrk(:,:,:,:) = 0.0d0
      rhobk(:,:,:,:,:) = 0.0d0
      rhodg(:,:,:,:,:) = 0.0d0
      rhobk(:,:,:,:,:) = 0.0d0
      rhoav(:,:,:,:) = 0.0d0
      ajdg(:,:,:,:,:) = 0.0d0
      ajav(:,:,:,:) = 0.0d0
      colf(:,:,:,:,:) = 0.0d0
      if(jobnum(1).gt.0) then
        if(jobnum(1).eq.1) then
          write(filename,'(a,i4.4,a)') './SNAPSHOT0/esdat', myid, '.h5'
        else
          write(filename,'(a,i4.4,a)') './SNAPSHOT0/emdat', myid, '.h5'
        end if
        call hdfopen(filename,fileid,DFACC_READ)
!
        dsname = 'rhobk'
        dims(1) = xu - xl + 1; dims(2) = yu - yl + 1; dims(3) = zu - zl + 1
        call read3d(fileid,dsname,dims(1:3), &
       &  rhobk(1,0:xu-xl,0:yu-yl,0:zu-zl,3),stats0,stats1)
!
        call hdfclose(fileid,stats0)
      end if
      rhodg(:,:,:,:,:) = 0.0d0
      colf(:,:,:,:,:) = 1.0d0
      if(xlcol(2).lt.xlcol(1).and.xucol(2).gt.xucol(1).and. &
     &   ylcol(2).lt.ylcol(1).and.yucol(2).gt.yucol(1).and. &
     &   zlcol(2).lt.zlcol(1).and.zucol(2).gt.zucol(1)) then
        xlmcol = 0.5d0*(xlcol(1) + xlcol(2))
        xumcol = 0.5d0*(xucol(1) + xucol(2))
        ylmcol = 0.5d0*(ylcol(1) + ylcol(2))
        yumcol = 0.5d0*(yucol(1) + yucol(2))
        zlmcol = 0.5d0*(zlcol(1) + zlcol(2))
        zumcol = 0.5d0*(zucol(1) + zucol(2))
        xlsgain = 15.0d0/abs(xlcol(2) - xlcol(1))
        xusgain = 15.0d0/abs(xucol(2) - xucol(1))
        ylsgain = 15.0d0/abs(ylcol(2) - ylcol(1))
        yusgain = 15.0d0/abs(yucol(2) - yucol(1))
        zlsgain = 15.0d0/abs(zlcol(2) - zlcol(1))
        zusgain = 15.0d0/abs(zucol(2) - zucol(1))
        do k=-1,zu-zl
        do j=-1,yu-yl
        do i=-1,xu-xl
          sigmxl = 0.5d0 + 0.5d0*tanh(-0.5*xlsgain*(i+xl-xlmcol))
          sigmxu = 0.5d0 + 0.5d0*tanh(+0.5*xusgain*(i+xl-xumcol))
          sigmyl = 0.5d0 + 0.5d0*tanh(-0.5*ylsgain*(j+yl-ylmcol))
          sigmyu = 0.5d0 + 0.5d0*tanh(+0.5*yusgain*(j+yl-yumcol))
          sigmzl = 0.5d0 + 0.5d0*tanh(-0.5*zlsgain*(k+zl-zlmcol))
          sigmzu = 0.5d0 + 0.5d0*tanh(+0.5*zusgain*(k+zl-zumcol))
          colf(1,i,j,k,1) = 1.0d0 - sigmxl*sigmxu*sigmyl*sigmyu*sigmzl*sigmzu
        end do
        end do
        end do
      end if


!-------------------- background field
!      do k=-1,zu-zl
!      do j=-1,yu-yl
!      do i=-1,xu-xl
!        eb(EX,i,j,k,3) = 
!        eb(EY,i,j,k,3) = 
!        eb(EZ,i,j,k,3) = 
!        eb(BX,i,j,k,3) = 
!        eb(BY,i,j,k,3) = 
!        eb(BZ,i,j,k,3) = 
!      end do
!      end do
!      end do


!-------------------- medium property (standard)
      do k=-1,zu-zl+1
      do j=-1,yu-yl+1
      do i=-1,xu-xl+1
        mp(EX,i,j,k,1) = 1.0d0
        mp(EY,i,j,k,1) = 1.0d0
        mp(EZ,i,j,k,1) = 1.0d0
        mp(BX,i,j,k,1) = 1.0d0
        mp(BY,i,j,k,1) = 1.0d0
        mp(BZ,i,j,k,1) = 1.0d0
      end do
      end do
      end do


!-------------------- initial wave field
       do iw=1,nwave
         if(twave(iw).eq.1) then
           angwv(:,iw) = angwv(:,iw)*pi/180.0d0
           orge(1:3) = orgwv(1:3,iw)
           orgb(1) = orgwv(1,iw) - cv*0.5d0*dt*sin(angwv(1,iw))*cos(angwv(2,iw))
           orgb(2) = orgwv(2,iw) - cv*0.5d0*dt*sin(angwv(1,iw))*sin(angwv(2,iw))
           orgb(3) = orgwv(3,iw) - cv*0.5d0*dt*cos(angwv(1,iw))
           do k=-1,zu-zl
           do j=-1,yu-yl
           do i=-1,xu-xl
!            -------- E-components
             x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(EX,i,j,k,1) = &
            &  ebamp(EX,iw)*sinfld(x,y,z,ldwv(iw),orge,angwv(:,iw))
             x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(EY,i,j,k,1) = &
            &  ebamp(EY,iw)*sinfld(x,y,z,ldwv(iw),orge,angwv(:,iw))
             x = (i + xl)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(EZ,i,j,k,1) = &
            &  ebamp(EZ,iw)*sinfld(x,y,z,ldwv(iw),orge,angwv(:,iw))
!            -------- B-components
             x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl + 0.5d0)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(BX,i,j,k,1) = &
            &  ebamp(BX,iw)*sinfld(x,y,z,ldwv(iw),orgb,angwv(:,iw))
             x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(BY,i,j,k,1) = &
            &  ebamp(BY,iw)*sinfld(x,y,z,ldwv(iw),orgb,angwv(:,iw))
             x = (i + xl + 0.5d0)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(BZ,i,j,k,1) = &
            &  ebamp(BZ,iw)*sinfld(x,y,z,ldwv(iw),orgb,angwv(:,iw))
           end do
           end do
           end do
         else if(twave(iw).eq.2) then
           angwv(:,iw) = angwv(:,iw)*pi/180.0d0
           orge(1:3) = orgwv(1:3,iw)
           orgb(1) = orgwv(1,iw) - cv*0.5d0*dt*sin(angwv(1,iw))*cos(angwv(2,iw))
           orgb(2) = orgwv(2,iw) - cv*0.5d0*dt*sin(angwv(1,iw))*sin(angwv(2,iw))
           orgb(3) = orgwv(3,iw) - cv*0.5d0*dt*cos(angwv(1,iw))
           do k=-1,zu-zl
           do j=-1,yu-yl
           do i=-1,xu-xl
!            -------- E-components
             x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(EX,i,j,k,1) = &
            &  ebamp(EX,iw)*gaus1dA(x,y,z,ldwv(iw),orge,angwv(:,iw))
             x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(EY,i,j,k,1) = &
            &  ebamp(EY,iw)*gaus1dA(x,y,z,ldwv(iw),orge,angwv(:,iw))
             x = (i + xl)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(EZ,i,j,k,1) = &
            &  ebamp(EZ,iw)*gaus1dA(x,y,z,ldwv(iw),orge,angwv(:,iw))
!            -------- B-components
             x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl + 0.5d0)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(BX,i,j,k,1) = &
            &  ebamp(BX,iw)*gaus1dA(x,y,z,ldwv(iw),orgb,angwv(:,iw))
             x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(BY,i,j,k,1) = &
            &  ebamp(BY,iw)*gaus1dA(x,y,z,ldwv(iw),orgb,angwv(:,iw))
             x = (i + xl + 0.5d0)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
             x = mod(x+slx,slx); y = mod(y+sly,sly); z = mod(z+slz,slz)
             eb(BZ,i,j,k,1) = &
            &  ebamp(BZ,iw)*gaus1dA(x,y,z,ldwv(iw),orgb,angwv(:,iw))
           end do
           end do
           end do

!         else if

         end if
       end do

       ! Bフィールドの直接代入
       ! if (bfield_path /= "") then
       !    eb(BX, i, j, k) = ***
       !    eb(BY, i, j, k) = ***
       !    eb(BZ, i, j, k) = ***
       ! end if


!-------------------- initial values of particle pushing fields
!!!      do k=1,nzm
!!!      do j=1,nym
!!!      do i=1,nxm
!!!        pex(i,j,k)=eb(EX,i,j,k)*mltstp
!!!        pey(i,j,k)=eb(EY,i,j,k)*mltstp
!!!        pez(i,j,k)=eb(EX,i,j,k)*mltstp
!!!        pbx(i,j,k)=eb(BX,i,j,k)*mltstp
!!!        pby(i,j,k)=eb(BY,i,j,k)*mltstp
!!!        pbz(i,j,k)=eb(BZ,i,j,k)*mltstp
!!!      end do
!!!      end do
!!!      end do


!-------------------- masking of transverse field
!!!      call fsmask(1)
!!!      call fbound(1)
!!!      call fsmask(2)
!!!      call fbound(2)
!!!      call fsmask(3)
!!!      call fbound(3)
!!!      call fsmask(4)
!!!      call fbound(4)


  return


  contains


  function sinfld(x,y,z,ld,org,ang)
    implicit none
    real(kind=8) :: sinfld
    real(kind=8),intent(in) :: x,y,z,ld,org(:),ang(:)
    real(kind=8) :: displ

    displ = (x - org(1))*sin(ang(1))*cos(ang(2)) &
   &      + (y - org(2))*sin(ang(1))*sin(ang(2)) &
   &      + (z - org(3))*cos(ang(1))

    sinfld = sin(2*pi*displ/ld)

    return
  end function sinfld


  function gaus1dA(x,y,z,pw,org,ang)
    implicit none
    real(kind=8) :: gaus1dA
    real(kind=8),intent(in) :: x,y,z,pw,org(:),ang(:)
    real(kind=8) :: displ

    displ = (x - org(1))*sin(ang(1))*cos(ang(2)) &
   &      + (y - org(2))*sin(ang(1))*sin(ang(2)) &
   &      + (z - org(3))*cos(ang(1))

    gaus1dA = exp(-displ*displ/pw/pw)

    return
  end function gaus1dA


  end subroutine inifld
