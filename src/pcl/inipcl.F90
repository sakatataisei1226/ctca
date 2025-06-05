#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine inipcl
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   I N I P C L
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine gives an initial setting of particles.  .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  use hdf
  use m_ohhinfo
#define MCW local_comm
  implicit none
!
  integer(kind=4) :: type_rdist = 1

  integer(kind=8) :: m,mm, ns,ne
  integer(kind=4) :: is
  integer(kind=4) :: ipc
  integer(kind=4) :: icon
  integer(kind=4) :: rid
  integer(kind=4) :: ierr
  real(kind=8) :: vpesq
  real(kind=8) :: vxtemp, vytemp, vztemp
  real(kind=8) :: disp1,disp2,disp3
  real(kind=8) :: radsq(inpc), rbwlsq, rdomsq
  real(kind=8) :: xsepa,ysepa,zsepa
  real(kind=8) :: vxyz, vxyztmp(ispec)
  real(kind=8) :: tfrac
  type(particle3) :: ptmp
!
  integer(kind=4) :: inttmp(ispec)
  integer(kind=HID_T) :: fileid
  integer(kind=4) :: stats0,stats1
  integer(kind=8) :: dims(1)
  real(kind=8) :: rens(4)
  character(len=30) :: filename,dsname


      call ohhinfo_update(sdoms(:, :, sdid(1) + 1))

      do ipc=1,npc
        if(geotype(ipc).eq.2) then
          radsq(ipc) = cylinder(ipc)%radius*cylinder(ipc)%radius
        else if(geotype(ipc).eq.3) then
          radsq(ipc) = sphere(ipc)%radius*sphere(ipc)%radius
        end if
      end do
      if(rbowl.gt.0.0d0) then
        rbwlsq = rbowl*rbowl
      else
        rbwlsq = 0.0d0
      end if
      if(rdome.gt.0.0d0) then
        rdomsq = rdome*rdome
      else
        rdomsq = 0.0d0
      end if


    if(jobnum(1).eq.0) then
      ! top of species loop 1
      nphgram(:,:,:) = 0
      ne = 0
ISL1: do is=1,nspec
        totalp(is,1) = npin(is)/nnode
        if(myid.lt.mod(npin(is),nnode)) then
          totalp(is,1) = totalp(is,1) + 1
        end if
        ns = ne + 1
        ne = ne + totalp(is,1)
!
        call RANU0(dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
        do m=ns,ne
          pbuf(m)%vx = dranu(m-ns+1)
        end do
        call RANU0(dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
        do m=ns,ne
          pbuf(m)%vy = dranu(m-ns+1)
        end do
        call RANU0(dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
        do m=ns,ne
          pbuf(m)%vz = dranu(m-ns+1)
        end do
!
        m = ns
        if(type_rdist.eq.1) then
          do mm=0,totalp(is,1)-1
            pbuf(m)%x = ngx*pbuf(ns+mm)%vx + xl
            pbuf(m)%y = ngy*pbuf(ns+mm)%vy + yl
            pbuf(m)%z = ngz*pbuf(ns+mm)%vz + zl
            ptmp%x = pbuf(m)%x
            ptmp%y = pbuf(m)%y
            ptmp%z = pbuf(m)%z
            rid = oh3_map_particle_to_subdomain &
           &        (ptmp%x,ptmp%y,ptmp%z)
            pbuf(m)%nid = rid
            do ipc=1,npc
              if((geotype(ipc).eq.0.or.geotype(ipc).eq.1).and. &
             &   (ptmp%x.ge.xlpc(ipc).and.ptmp%x.le.xupc(ipc).and. &
             &    ptmp%y.ge.ylpc(ipc).and.ptmp%y.le.yupc(ipc).and. &
             &    ptmp%z.ge.zlpc(ipc).and.ptmp%z.le.zupc(ipc))) then
                pbuf(m)%x = 0.0d0
                pbuf(m)%y = 0.0d0
                pbuf(m)%z = 0.0d0
                pbuf(m)%nid = -1
              else if(geotype(ipc).eq.2) then
                if(cylinder(ipc)%align.eq.1) then
                  disp1 = ptmp%y - cylinder(ipc)%axis(1)
                  disp2 = ptmp%z - cylinder(ipc)%axis(2)
                  if(ptmp%x.ge.xlpc(ipc).and.ptmp%x.le.xupc(ipc).and. &
                 &   disp1*disp1+disp2*disp2.lt.radsq(ipc)) then
                    pbuf(m)%x = 0.0d0
                    pbuf(m)%y = 0.0d0
                    pbuf(m)%z = 0.0d0
                    pbuf(m)%nid = -1
                  end if
                else if(cylinder(ipc)%align.eq.2) then
                  disp1 = ptmp%z - cylinder(ipc)%axis(1)
                  disp2 = ptmp%x - cylinder(ipc)%axis(2)
                  if(ptmp%y.ge.ylpc(ipc).and.ptmp%y.le.yupc(ipc).and. &
                 &   disp1*disp1+disp2*disp2.lt.radsq(ipc)) then
                    pbuf(m)%x = 0.0d0
                    pbuf(m)%y = 0.0d0
                    pbuf(m)%z = 0.0d0
                    pbuf(m)%nid = -1
                  end if
                else if(cylinder(ipc)%align.eq.3) then
                  disp1 = ptmp%x - cylinder(ipc)%axis(1)
                  disp2 = ptmp%y - cylinder(ipc)%axis(2)
                  if(ptmp%z.ge.zlpc(ipc).and.ptmp%z.le.zupc(ipc).and. &
                 &   disp1*disp1+disp2*disp2.lt.radsq(ipc)) then
                    pbuf(m)%x = 0.0d0
                    pbuf(m)%y = 0.0d0
                    pbuf(m)%z = 0.0d0
                    pbuf(m)%nid = -1
                  end if
                end if
              else if(geotype(ipc).eq.3) then
                disp1 = ptmp%x - sphere(ipc)%center(1)
                disp2 = ptmp%y - sphere(ipc)%center(2)
                disp3 = ptmp%z - sphere(ipc)%center(3)
                if(disp1*disp1+disp2*disp2+disp3*disp3.lt.radsq(ipc)) then
                  pbuf(m)%x = 0.0d0
                  pbuf(m)%y = 0.0d0
                  pbuf(m)%z = 0.0d0
                  pbuf(m)%nid = -1
                end if
              end if
            end do
!
#include "defsurf_init.fnc"
!
            if(pbuf(m)%preside == OH_PCL_TO_BE_ACCUMULATED) then
              pbuf(m)%x = 0.0d0
              pbuf(m)%y = 0.0d0
              pbuf(m)%z = 0.0d0
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            end if
!
            if(pbuf(m)%nid.ge.0) then
              nphgram(pbuf(m)%nid+1,is,1) = nphgram(pbuf(m)%nid+1,is,1) + 1
              m = m + 1
            end if
          end do
        else
          do mm=0,totalp(is,1)-1
            pbuf(m)%x = slx*pbuf(ns+mm)%vx
            pbuf(m)%y = sly*pbuf(ns+mm)%vy
            pbuf(m)%z = slz*pbuf(ns+mm)%vz
            ptmp%x = pbuf(m)%x/dr
            ptmp%y = pbuf(m)%y/dr
            ptmp%z = pbuf(m)%z/dr
            rid = oh3_map_particle_to_subdomain &
           &        (ptmp%x,ptmp%y,ptmp%z)
            pbuf(m)%nid = rid
            do ipc=1,npc
              if((geotype(ipc).eq.0.or.geotype(ipc).eq.1).and. &
             &   (ptmp%x.ge.xlpc(ipc).and.ptmp%x.le.xupc(ipc).and. &
             &    ptmp%y.ge.ylpc(ipc).and.ptmp%y.le.yupc(ipc).and. &
             &    ptmp%z.ge.zlpc(ipc).and.ptmp%z.le.zupc(ipc))) then
                pbuf(m)%x = 0.0d0
                pbuf(m)%y = 0.0d0
                pbuf(m)%z = 0.0d0
                pbuf(m)%nid = -1
              else if(geotype(ipc).eq.2) then
                if(cylinder(ipc)%align.eq.1) then
                  disp1 = ptmp%y - cylinder(ipc)%axis(1)
                  disp2 = ptmp%z - cylinder(ipc)%axis(2)
                  if(ptmp%x.ge.xlpc(ipc).and.ptmp%x.le.xupc(ipc).and. &
                 &   disp1*disp1+disp2*disp2.lt.radsq(ipc)) then
                    pbuf(m)%x = 0.0d0
                    pbuf(m)%y = 0.0d0
                    pbuf(m)%z = 0.0d0
                    pbuf(m)%nid = -1
                  end if
                else if(cylinder(ipc)%align.eq.2) then
                  disp1 = ptmp%z - cylinder(ipc)%axis(1)
                  disp2 = ptmp%x - cylinder(ipc)%axis(2)
                  if(ptmp%y.ge.ylpc(ipc).and.ptmp%y.le.yupc(ipc).and. &
                 &   disp1*disp1+disp2*disp2.lt.radsq(ipc)) then
                    pbuf(m)%x = 0.0d0
                    pbuf(m)%y = 0.0d0
                    pbuf(m)%z = 0.0d0
                    pbuf(m)%nid = -1
                  end if
                else if(cylinder(ipc)%align.eq.3) then
                  disp1 = ptmp%x - cylinder(ipc)%axis(1)
                  disp2 = ptmp%y - cylinder(ipc)%axis(2)
                  if(ptmp%z.ge.zlpc(ipc).and.ptmp%z.le.zupc(ipc).and. &
                 &   disp1*disp1+disp2*disp2.lt.radsq(ipc)) then
                    pbuf(m)%x = 0.0d0
                    pbuf(m)%y = 0.0d0
                    pbuf(m)%z = 0.0d0
                    pbuf(m)%nid = -1
                  end if
                end if
              else if(geotype(ipc).eq.3) then
                disp1 = ptmp%x - sphere(ipc)%center(1)
                disp2 = ptmp%y - sphere(ipc)%center(2)
                disp3 = ptmp%z - sphere(ipc)%center(3)
                if(disp1*disp1+disp2*disp2+disp3*disp3.lt.radsq(ipc)) then
                  pbuf(m)%x = 0.0d0
                  pbuf(m)%y = 0.0d0
                  pbuf(m)%z = 0.0d0
                  pbuf(m)%nid = -1
                end if
              end if
            end do
!
#include "defsurf_init.fnc"
!
            if(pbuf(m)%preside == OH_PCL_TO_BE_ACCUMULATED) then
              pbuf(m)%x = 0.0d0
              pbuf(m)%y = 0.0d0
              pbuf(m)%z = 0.0d0
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            end if
!
            if(pbuf(m)%nid.ge.0) then
              nphgram(pbuf(m)%nid+1,is,1) = nphgram(pbuf(m)%nid+1,is,1) + 1
              m = m + 1
            end if
          end do
        end if
        totalp(is,1) = m - ns
        ne = m - 1
!
        if(lcgamma(is).ge.1.0d0.or.lcgamma(is).lt.0.0d0.or. &
       &   lcbeta(is).le.0.0d0) then
          call RANN0(spe(is)*dcos(speth(is)/180.0d0*pi),peth(is), &
         &            dranu,totalp(is,1),icon)
          if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
          do m=ns,ne
            pbuf(m)%vx = dranu(m-ns+1)
          end do
          call RANN0(spe(is)*dsin(speth(is)/180.0d0*pi),peth(is), &
         &            dranu,totalp(is,1),icon)
          if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
          do m=ns,ne
            pbuf(m)%vy = dranu(m-ns+1)
          end do
        else
          print*, "WARNING!!"
          call RANN0(spe(is)*dcos(speth(is)/180.0d0*pi),peth(is), &
         &            dranu,int(totalp(is,1)*lcgamma(is)),icon)
          if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon
          do m=ns,ns-1+int(totalp(is,1)*lcgamma(is))
            pbuf(m)%vx = dranu(m-ns+1)
          end do
          call RANN0(spe(is)*dsin(speth(is)/180.0d0*pi),peth(is), &
         &            dranu,int(totalp(is,1)*lcgamma(is)),icon)
          if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon
          do m=ns,ns-1+int(totalp(is,1)*lcgamma(is))
            pbuf(m)%vy = dranu(m-ns+1)
          end do
!
          m = ns + int(totalp(is,1)*lcgamma(is))
          do while(m.le.ne)
            call RANN0(spe(is),peth(is),dranu(1:1),1,icon)
            call RANN0(spe(is),peth(is),dranu(2:2),1,icon)
            vpesq = dranu(1)*dranu(1) + dranu(2)*dranu(2)
            call RANU0(dranu(3:3),1,icon)
            if(exp(-vpesq/lcbeta(is)/peth(is)/peth(is)).le.dranu(3)) then
              pbuf(m)%vx = dranu(1) + spe(is)*dcos(speth(is)/180.0d0*pi)
              pbuf(m)%vy = dranu(2) + spe(is)*dsin(speth(is)/180.0d0*pi)
              m = m + 1
            end if
          end do
        end if

        call RANN0(spa(is),path(is),dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
        do m=ns,ne
          pbuf(m)%vz = dranu(m-ns+1)
        end do
!
        do m=ns,ne
          vxtemp = pbuf(m)%vx*t11 + pbuf(m)%vy*t12 + pbuf(m)%vz*t13
          vytemp = pbuf(m)%vx*t21 + pbuf(m)%vy*t22 + pbuf(m)%vz*t23
          vztemp = pbuf(m)%vx*t31 + pbuf(m)%vy*t32 + pbuf(m)%vz*t33
          vxtemp = vxtemp &
         &       + vdri(is)*sin(vdthz(is)/180.0d0*pi)*cos(vdthxy(is)/180.0d0*pi)
          vytemp = vytemp &
         &       + vdri(is)*sin(vdthz(is)/180.0d0*pi)*sin(vdthxy(is)/180.0d0*pi)
          vztemp = vztemp &
         &       + vdri(is)*cos(vdthz(is)/180.0d0*pi)
          pbuf(m)%vx = vxtemp + vdx(is)
          pbuf(m)%vy = vytemp + vdy(is)
          pbuf(m)%vz = vztemp + vdz(is)
        end do
!
!        if((ewmodel.eq.1.or.ewmodel.eq.2).and.is.eq.1.and.phiz.eq.0.0d0) then
!          do m=ns,ne
!            pbuf(m)%vx = pbuf(m)%vx + Ew(1,is,1)*qm(is)*omegaw(1)*id2omega*sin(0.0d0) &
!           &                        + Ew(2,is,1)*qm(is)*wc*id2omega*sin(0.0d0)
!            pbuf(m)%vy = pbuf(m)%vy - Ew(1,is,1)*qm(is)*wc*id2omega*cos(0.0d0) &
!           &                        - Ew(2,is,1)*qm(is)*omegaw(1)*id2omega*cos(0.0d0)
!            pbuf(m)%vz = pbuf(m)%vz
!          end do
!        end if
!
        do m=ns,ne
          pbuf(m)%spec = is
          pbuf(m)%preside = OH_PCL_ALIVE
          pbuf(m)%pid = 0
        end do
!
        do m=ns,ne
          vxyz = sqrt((pbuf(m)%vx-vdtx(is))*(pbuf(m)%vx-vdtx(is)) &
         &          + (pbuf(m)%vy-vdty(is))*(pbuf(m)%vy-vdty(is)) &
         &          + (pbuf(m)%vz-vdtz(is))*(pbuf(m)%vz-vdtz(is)))
          if(vxyzmax(is).lt.vxyz) vxyztmp(is) = vxyz
        end do
      end do ISL1
!
      call MPI_Allreduce(vxyztmp,vxyzmax,minspec,MPI_REAL8,MPI_MAX,MCW,ierr)
      if(myid.eq.0) print*, "vxyzmax =",vxyzmax(1:minspec)
!
      do is=1,nspec
        vxyzmax(is) = vxyzmax(is)*1.25d0
      end do
!
    else
!****************************** top of species loop 2
      nphgram(:,:,:) = 0
      if(jobnum(1).eq.1) then
        write(filename,'(a,i4.4,a)') './SNAPSHOT0/esdat', myid, '.h5'
      else
        write(filename,'(a,i4.4,a)') './SNAPSHOT0/emdat', myid, '.h5'
      end if
      call hdfopen(filename,fileid,DFACC_READ)
!
      dsname = 'np'
      dims(1) = nspec
      call read1i(fileid,dsname,dims(1:1),inttmp(1:nspec),stats0,stats1)
!
      dsname = 'rens'
      dims(1) = 4
      call read1d(fileid,dsname,dims(1:1),rens(1:4),stats0,stats1)
!
      dsname = 'px'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%x,stats0,stats1)
!
      dsname = 'py'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%y,stats0,stats1)
!
      dsname = 'pz'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%z,stats0,stats1)
!
      dsname = 'pvx'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%vx,stats0,stats1)
!
      dsname = 'pvy'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%vy,stats0,stats1)
!
      dsname = 'pvz'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%vz,stats0,stats1)
!
      dsname = 'pres'
      dims(1) = sum(inttmp(1:nspec))
      call read1i(fileid,dsname,dims(1:1),pbuf(1:dims(1))%preside,stats0,stats1)
!
      dsname = 'pis'
      dims(1) = sum(inttmp(1:nspec))
      call read1i(fileid,dsname,dims(1:1),pbuf(1:dims(1))%spec,stats0,stats1)
!
      totalp(1:nspec,1) = inttmp(1:nspec)
!
      call specsort
!
      do m=1,sum(inttmp(1:nspec))
        ptmp%x = pbuf(m)%x/rens(1)/dr
        ptmp%y = pbuf(m)%y/rens(1)/dr
        ptmp%z = pbuf(m)%z/rens(1)/dr
        rid = oh3_map_particle_to_subdomain &
       &        (ptmp%x,ptmp%y,ptmp%z)
        pbuf(m)%nid = rid
        if(rid.ge.0) &
       &  nphgram(rid+1,pbuf(m)%spec,1) = nphgram(rid+1,pbuf(m)%spec,1) + 1
        pbuf(m)%pid = 0
      end do
!
      call hdfclose(fileid,stats0)
    end if


  return
  end subroutine
