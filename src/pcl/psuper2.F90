#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine psuper2(ustep)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   P S U P E R
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine supervises particle injection from      .
!   .  the outer/inner boundaries of the simulation box.       .
!   ............................................................
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
#define MCW local_comm
#define MSS MPI_STATUS_SIZE
  implicit none
!
  integer(kind=4) :: m, l
  integer(kind=8) :: nprs,npre
  integer(kind=8) :: nofluxt(inpc)
  integer(kind=8) :: nofluxsum
  integer(kind=4) :: i,j,k, i1,j1,k1
  integer(kind=4) :: ib,jb,kb
  integer(kind=4) :: ia,ja,ka, i1a,j1a,k1a
  integer(kind=4) :: nns,nne, injctisw(ispec)
  integer(kind=4) :: ie, ipcplane, nemitw
  integer(kind=4) :: iran, iepl, ipc, is, iss, ustep
  integer(kind=4) :: addr0,addr1,addr2
  integer(kind=4) :: icon
  integer(kind=4) :: nemit(inepl), omni
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: idim, rid
  integer(kind=4) :: ierr
!  integer(kind=4) :: oh3_map_region_to_node
  real(kind=8) :: fnemit
  real(kind=8) :: align, radius, axis(2)
  real(kind=8) :: dxl,dxu, dyl,dyu, dzl,dzu
  real(kind=8) :: xlocalb,ylocalb,zlocalb
  real(kind=8) :: xb,yb,zb
  real(kind=8) :: xlocala,ylocala,zlocala
  real(kind=8) :: xa,ya,za
  real(kind=8) :: xr,yr,zr
  real(kind=8) :: xd2,yd2,zd2
  real(kind=8) :: tslx,tsly, rustep, qs
  real(kind=8) :: vvx2,vvy2,vvz2
  real(kind=8) :: vx2w1,vx2w2,vx2w3,vx2w4
  real(kind=8) :: vy2w1,vy2w2,vy2w3,vy2w4
  real(kind=8) :: vz2w1,vz2w2,vz2w3,vz2w4
  real(kind=8) :: csz,snz, csxy,snxy, xew,yew,zew
  real(kind=8) :: te11,te12,te13,te21,te22,te23,te31,te32,te33
  real(kind=8) :: arearatio(12)
  real(kind=8) :: tanthetaz,costhetaxy,sinthetaxy
  real(kind=8) :: tantzcostxy,tantzsintxy
  real(kind=8) :: betav, zdepth, rtmp
  real(kind=8) :: xsepa1,ysepa1, xsepa2,ysepa2
  real(kind=8) :: xlocal,ylocal,zlocal
  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8
!  type(oh_particle) :: pinj


!-------------------- 
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)
      dxl = xl; dxu = xu; dyl = yl; dyu = yu; dzl = zl; dzu = zu


!-------------------- 
!     call MPI_AllReduce(sey,seygl,15*ispec,MPI_DOUBLE_PRECISION,MPI_SUM,MCW,ierr)
!     seygl(15,:) = seygl(12,:)*0.25d0
!     seygl(14,:) = seygl(12,:)*0.25d0
!     seygl(13,:) = seygl(12,:)*0.25d0
!     seygl(12,:) = seygl(12,:)*0.25d0


!-------------------- 
      tanthetaz = tan(thetaz)
      costhetaxy = cos(thetaxy)
      sinthetaxy = sin(thetaxy)
      tantzcostxy = tanthetaz*costhetaxy
      tantzsintxy = tanthetaz*sinthetaxy


!-------------------- supervision 1
!
!     injection of particles from inner/outer boundaries
!
      do is=1,nspec
        injctisw(is) = injct(is)*ustep
      end do
!
!****************************** top of species loop 1
      npre = 0
      nne = 0
      rustep = 1.0d0/ustep
      tslx = 2.0d0*slx
      tsly = 2.0d0*sly
ISL1: do is=1,nspec
        nprs = npre + 1
        npre = npre + npr(is)
        nns = nne + 1
        nne = nne + nepl(is)
!
!============================== from outer boundary
        if((nflag_emit(is).eq.0).and. &
       &   (npbnd(1,is).eq.2.or.npbnd(2,is).eq.2.or.npbnd(3,is).eq.2)) then
!-------------------- x-bottom
          if(npbnd(1,is).eq.2.and.bared(1,1,sdid(1)+1).eq.1) then
            narrxf(is) = (int((istep + lstep)*arrxf(is),8) &
           &            - int((istep + lstep - injct(is))*arrxf(is),8))*ustep
            if(narrxf(is).gt.0) then
              call RANU0(dranu,narrxf(is)*5,icon)
              iran = 0
              do m=1,narrxf(is)
                iran = iran + 5
!               ----- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ----- assign position and velocity
                pinj(1)%x = vxf(addr1)*injctisw(is)*dranu(iran-2)
                pinj(1)%y = yl + (yu - yl)*dranu(iran-1)
                pinj(1)%z = zl + (zu - zl)*dranu(iran  )
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxf(addr1)
                pinj(1)%vy = vyr(addr2)
                pinj(1)%vz = vzr(addr2)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%preside = 0
                pinj(1)%pid = 0
!
                call oh2_inject_particle(pinj(1))
              end do
            end if
          end if

!-------------------- x-top
          if(npbnd(1,is).eq.2.and.bared(2,1,sdid(1)+1).eq.1) then
            narrxb(is) = (int((istep + lstep)*arrxb(is),8) &
           &            - int((istep + lstep - injct(is))*arrxb(is),8))*ustep
            if(narrxb(is).gt.0) then
              call RANU0(dranu,narrxb(is)*5,icon)
              iran = 0
              do m=1,narrxb(is)
                iran = iran + 5
!               ----- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ----- assign position and velocity
                pinj(1)%x = xu + vxb(addr1)*injctisw(is)*dranu(iran-2)
                pinj(1)%y = yl + (yu - yl)*dranu(iran-1)
                pinj(1)%z = zl + (zu - zl)*dranu(iran  )
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxb(addr1)
                pinj(1)%vy = vyr(addr2)
                pinj(1)%vz = vzr(addr2)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%preside = 0
                pinj(1)%pid = 0
!
                call oh2_inject_particle(pinj(1))
              end do
            end if
          end if

!-------------------- y-bottom
          if(npbnd(2,is).eq.2.and.bared(1,2,sdid(1)+1).eq.1) then
            narryf(is) = (int((istep + lstep)*arryf(is),8) &
           &            - int((istep + lstep - injct(is))*arryf(is),8))*ustep
            if(narryf(is).gt.0) then
              call RANU0(dranu,narryf(is)*5,icon)
              iran = 0
              do m=1,narryf(is)
                iran = iran + 5
!               ----- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ----- assign position and velocity
                pinj(1)%x = xl + (xu - xl)*dranu(iran-2)
                pinj(1)%y = vyf(addr1)*injctisw(is)*dranu(iran-1)
                pinj(1)%z = zl + (zu - zl)*dranu(iran  )
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxr(addr2)
                pinj(1)%vy = vyf(addr1)
                pinj(1)%vz = vzr(addr2)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%preside = 0
                pinj(1)%pid = 0
!
                call oh2_inject_particle(pinj(1))
              end do
            end if
          end if

!-------------------- y-top
          if(npbnd(2,is).eq.2.and.bared(2,2,sdid(1)+1).eq.1) then
            narryb(is) = (int((istep + lstep)*arryb(is),8) &
           &            - int((istep + lstep - injct(is))*arryb(is),8))*ustep
            if(narryb(is).gt.0) then
              call RANU0(dranu,narryb(is)*5,icon)
              iran = 0
              do m=1,narryb(is)
                iran = iran + 5
!               ----- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ----- assign position and velocity
                pinj(1)%x = xl + (xu - xl)*dranu(iran-2)
                pinj(1)%y = yu + vyb(addr1)*injctisw(is)*dranu(iran-1)
                pinj(1)%z = zl + (zu - zl)*dranu(iran  )
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxr(addr2)
                pinj(1)%vy = vyb(addr1)
                pinj(1)%vz = vzr(addr2)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%preside = 0
                pinj(1)%pid = 0
!
                call oh2_inject_particle(pinj(1))
              end do
            end if
          end if


!-------------------- z-bottom
          if(npbnd(3,is).eq.2.and.bared(1,3,sdid(1)+1).eq.1) then
            narrzf(is) = (int((istep + lstep)*arrzf(is),8) &
           &            - int((istep + lstep - injct(is))*arrzf(is),8))*ustep
            if(narrzf(is).gt.0) then
              call RANU0(dranu,narrzf(is)*5,icon)
              iran = 0
              do m=1,narrzf(is)
                iran = iran + 5
!               ----- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ----- assign position and velocity
                pinj(1)%x = xl + (xu - xl)*dranu(iran-2)
                pinj(1)%y = yl + (yu - yl)*dranu(iran-1)
                pinj(1)%z = vzf(addr1)*injctisw(is)*dranu(iran  )
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxr(addr2)
                pinj(1)%vy = vyr(addr2)
                pinj(1)%vz = vzf(addr1)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%preside = 0
                pinj(1)%pid = 0
!
                call oh2_inject_particle(pinj(1))
              end do
            end if
          end if

!-------------------- z-top
          if(npbnd(3,is).eq.2.and.bared(2,3,sdid(1)+1).eq.1) then
            narrzb(is) = (int((istep + lstep)*arrzb(is),8) &
           &            - int((istep + lstep - injct(is))*arrzb(is),8))*ustep
            if(narrzb(is).gt.0) then
              call RANU0(dranu,narrzb(is)*5,icon)
              iran = 0
              do m=1,narrzb(is)
                iran = iran + 5
!               ----- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ----- assign position and velocity
                pinj(1)%x = xl + (xu - xl)*dranu(iran-2)
                pinj(1)%y = yl + (yu - yl)*dranu(iran-1)
                pinj(1)%z = zu + vzb(addr1)*injctisw(is)*dranu(iran  )
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxr(addr2)
                pinj(1)%vy = vyr(addr2)
                pinj(1)%vz = vzb(addr1)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%preside = 0
                pinj(1)%pid = 0
!
                call oh2_inject_particle(pinj(1))
              end do
            end if
          end if

!============================== from inner boundary
        else if(nflag_emit(is).ge.1.and.nflag_emit(is).le.2) then
          if(injct(is).eq.0) cycle
          ie = mod(istep,injct(is))
          if(ie.ne.0.or.inpf(is).eq.0) cycle
!
!------------------------------ plane loop
          do iepl=nns,0
!            if(nflag_emit(is).eq.1) then
              ! fluxf = (curfs(is)/abs(q(is)))*xl*yl*dt*0.5d0
              ! emit = fluxf * ustep
              nemit(iepl) = (int((istep + lstep)*fluxf(iepl),8) &
             &            -  int((istep + lstep - injct(is))*fluxf(iepl),8))*ustep
!            else if(nflag_emit(is).eq.2) then
!              do iss=1,nspec
!                fnemit = fnemit + seygl(iepl-nns+1,iss)*q(iss)/q(is)
!              end do
!              fnemit = fnemit*peject(iepl)%area/peject(iepl)%tarea
!              call RANU0(dranu,1,icon)
!              if(fnemit-int(fnemit).le.dranu(1)) then
!                nemit(iepl) = int(fnemit)
!              else
!                nemit(iepl) = int(fnemit) + 1
!              end if
!            end if
            ipcplane = ipcpl(iepl)
            nemitw = nemit(iepl)
            omni = omniemit(iepl)
            if(nemitw.gt.0) then
              call RANU0(dranu,nemitw*6,icon)
              if(icon.ne.0) print*,"[psuper] icon=",icon
            end if
            iran = 0
!
!           ------------------- pos. & vel. assignment
            if(geotype(ipcplane).eq.0.or.geotype(ipcplane).eq.1) then
              do m=1,nemit(iepl)
                if(nemd(iepl).eq. 1) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz,pinj(1)%x,pinj(1)%y,pinj(1)%z,+1, &
                 &   peject(iepl)%grd,peject(iepl)%yl,peject(iepl)%yu, &
                 &   peject(iepl)%zl,peject(iepl)%zu,m*6-5)
                else if(nemd(iepl).eq.-1) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz,pinj(1)%x,pinj(1)%y,pinj(1)%z,-1, &
                 &   peject(iepl)%grd,peject(iepl)%yl,peject(iepl)%yu, &
                 &   peject(iepl)%zl,peject(iepl)%zu,m*6-5)
                else if(nemd(iepl).eq. 2) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vy,pinj(1)%vz,pinj(1)%vx,pinj(1)%y,pinj(1)%z,pinj(1)%x,+1, &
                 &   peject(iepl)%grd,peject(iepl)%zl,peject(iepl)%zu, &
                 &   peject(iepl)%xl,peject(iepl)%xu,m*6-5)
                else if(nemd(iepl).eq.-2) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vy,pinj(1)%vz,pinj(1)%vx,pinj(1)%y,pinj(1)%z,pinj(1)%x,-1, &
                 &   peject(iepl)%grd,peject(iepl)%zl,peject(iepl)%zu, &
                 &   peject(iepl)%xl,peject(iepl)%xu,m*6-5)
                else if(nemd(iepl).eq. 3) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy,pinj(1)%z,pinj(1)%x,pinj(1)%y,+1, &
                 &   peject(iepl)%grd,peject(iepl)%xl,peject(iepl)%xu, &
                 &   peject(iepl)%yl,peject(iepl)%yu,m*6-5)
                else if(nemd(iepl).eq.-3) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy,pinj(1)%z,pinj(1)%x,pinj(1)%y,-1, &
                 &   peject(iepl)%grd,peject(iepl)%xl,peject(iepl)%xu, &
                 &   peject(iepl)%yl,peject(iepl)%yu,m*6-5)
                end if
!
                if(pinj(1)%x.lt.xl.or.pinj(1)%x.ge.xu.or. &
               &   pinj(1)%y.lt.yl.or.pinj(1)%y.ge.yu.or. &
               &   pinj(1)%z.lt.zl.or.pinj(1)%z.ge.zu) then
                  rid = oh3_map_particle_to_subdomain(pinj(1)%x,pinj(1)%y,pinj(1)%z)
                  pinj(1)%nid = rid
                else
                  pinj(1)%nid = sdid(1)
                end if
                pinj(1)%spec = is
                pinj(1)%pid = ipcplane
!
                call oh2_inject_particle(pinj(1))
              end do
            else if(geotype(ipcplane).eq.2) then
              align = cylinder(ipcplane)%align
              radius = cylinder(ipcplane)%radius
              if(abs(nemd(iepl)).eq.align) then
                do m=1,nemit(iepl)
                  if     (nemd(iepl).eq. 1) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz,pinj(1)%x,pinj(1)%y,pinj(1)%z,+1, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%grd,m*6-5)
                  else if(nemd(iepl).eq.-1) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz,pinj(1)%x,pinj(1)%y,pinj(1)%z,-1, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%grd,m*6-5)
                  else if(nemd(iepl).eq. 2) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vy,pinj(1)%vz,pinj(1)%vx,pinj(1)%y,pinj(1)%z,pinj(1)%x,+1, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%grd,m*6-5)
                  else if(nemd(iepl).eq.-2) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vy,pinj(1)%vz,pinj(1)%vx,pinj(1)%y,pinj(1)%z,pinj(1)%x,-1, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%grd,m*6-5)
                  else if(nemd(iepl).eq. 3) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy,pinj(1)%z,pinj(1)%x,pinj(1)%y,+1, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%grd,m*6-5)
                  else if(nemd(iepl).eq.-3) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy,pinj(1)%z,pinj(1)%x,pinj(1)%y,-1, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%grd,m*6-5)
                  end if
!
                  if(pinj(1)%x.lt.xl.or.pinj(1)%x.ge.xu.or. &
                 &   pinj(1)%y.lt.yl.or.pinj(1)%y.ge.yu.or. &
                 &   pinj(1)%z.lt.zl.or.pinj(1)%z.ge.zu) then
                    rid = oh3_map_particle_to_subdomain(pinj(1)%x,pinj(1)%y,pinj(1)%z)
                    pinj(1)%nid = rid
                  else
                    pinj(1)%nid = sdid(1)
                  end if
                  pinj(1)%spec = is
                  pinj(1)%pid = ipcplane
!
                  call oh2_inject_particle(pinj(1))
                end do
              else
                do m=1,nemit(iepl)
                  if     (align.eq.1.and.omni.eq.0.and.radius.ge.1.0d0) then
                    call emit_from_cylindrical_surf1 &
                   &  (pinj(1)%vz,pinj(1)%vy,pinj(1)%vx,pinj(1)%z,pinj(1)%y,pinj(1)%x,psizy, &
                   &   cylinder(ipcplane)%axis(2),cylinder(ipcplane)%axis(1), &
                   &   peject(iepl)%xl,peject(iepl)%xu,m*6-5)
                  else if(align.eq.1.and.omni.eq.0.and.radius.lt.1.0d0) then
                    call emit_from_cylindrical_surf2 &
                   &  (pinj(1)%vz,pinj(1)%vy,pinj(1)%vx,pinj(1)%z,pinj(1)%y,pinj(1)%x,psizy, &
                   &   cylinder(ipcplane)%axis(2),cylinder(ipcplane)%axis(1), &
                   &   peject(iepl)%xl,peject(iepl)%xu,m*5-4)
                  else if(align.eq.1.and.omni.eq.1) then
                    call emit_from_cylindrical_surf3 &
                   &  (pinj(1)%vz,pinj(1)%vy,pinj(1)%vx,pinj(1)%z,pinj(1)%y,pinj(1)%x,psizy, &
                   &   cylinder(ipcplane)%axis(2),cylinder(ipcplane)%axis(1), &
                   &   peject(iepl)%xl,peject(iepl)%xu,m*6-5)
                  else if(align.eq.2.and.omni.eq.0.and.radius.ge.1.0d0) then
                    call emit_from_cylindrical_surf1 &
                   &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy,pinj(1)%z,pinj(1)%x,pinj(1)%y,psizx, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%yl,peject(iepl)%yu,m*6-5)
                  else if(align.eq.2.and.omni.eq.0.and.radius.lt.1.0d0) then
                    call emit_from_cylindrical_surf2 &
                   &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy,pinj(1)%z,pinj(1)%x,pinj(1)%y,psizx, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%yl,peject(iepl)%yu,m*5-4)
                  else if(align.eq.2.and.omni.eq.1) then
                    call emit_from_cylindrical_surf3 &
                   &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy,pinj(1)%z,pinj(1)%x,pinj(1)%y,psizx, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%yl,peject(iepl)%yu,m*6-5)
                  else if(align.eq.3.and.omni.eq.0.and.radius.ge.1.0d0) then
                    call emit_from_cylindrical_surf1 &
                   &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz,pinj(1)%x,pinj(1)%y,pinj(1)%z,psixy, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%zl,peject(iepl)%zu,m*6-5)
                  else if(align.eq.3.and.omni.eq.0.and.radius.lt.1.0d0) then
                    call emit_from_cylindrical_surf2 &
                   &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz,pinj(1)%x,pinj(1)%y,pinj(1)%z,psixy, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%zl,peject(iepl)%zu,m*5-4)
                  else if(align.eq.3.and.omni.eq.1) then
                    call emit_from_cylindrical_surf3 &
                   &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz,pinj(1)%x,pinj(1)%y,pinj(1)%z,psixy, &
                   &   cylinder(ipcplane)%axis(1),cylinder(ipcplane)%axis(2), &
                   &   peject(iepl)%zl,peject(iepl)%zu,m*6-5)
                  end if
!
                  if(pinj(1)%x.lt.xl.or.pinj(1)%x.ge.xu.or. &
                 &   pinj(1)%y.lt.yl.or.pinj(1)%y.ge.yu.or. &
                 &   pinj(1)%z.lt.zl.or.pinj(1)%z.ge.zu) then
                    rid = oh3_map_particle_to_subdomain(pinj(1)%x,pinj(1)%y,pinj(1)%z)
                    pinj(1)%nid = rid
                  else
                    pinj(1)%nid = sdid(1)
                  end if
                  pinj(1)%spec = is
                  pinj(1)%pid = ipcplane
!
                  call oh2_inject_particle(pinj(1))
                end do
              end if
!            else if(geotype(ipcplane).eq.3) then
!              !under construction!
            end if
!
!           --------- calculate accumulated charge
            gcount(1)%chgacm(:,ipcplane) = gcount(1)%chgacm(:,ipcplane) &
           &                             - q(is)*nemitw*sqdscaled(ipcplane)
            gcount(1)%outflux(ipcplane,is) = gcount(1)%outflux(ipcplane,is) &
           &                               + nemitw
          end do

!-------------------- from land surface


!-------------------- from hole basement


!-------------------- from hole flank


        end if
      end do ISL1
!****************************** end of species loop 1


!-------------------- other treatments
!     --------------- charge transfer simulating bias current
      if(myid.eq.0) then
        i = 1
        do while(biasc(i)%to.ne.0)
          if(biasc(i)%to.gt.0) then
            if(biasc(i)%from.gt.0) &
           &  gcount(1)%chgacm(:,biasc(i)%from) = &
           &  gcount(1)%chgacm(:,biasc(i)%from) - biasc(i)%val*ustep
            gcount(1)%chgacm(:,biasc(i)%to) = &
           &  gcount(1)%chgacm(:,biasc(i)%to) + biasc(i)%val*ustep
          end if
          i = i + 1
        end do
      end if


!-------------------- diagnostics
!    if(myid.eq.0) then
!      if(intfoc.ne.0) then
!        nofluxsum = 0
!        do ipc=1,npc
!          nofluxsum = nofluxsum + nofluxt(ipc)
!        end do
!        if(istep.eq.0.or.mod(istep-1,intfoc).eq.0) &
!       &  open(89,file='noflux',position='append')
!        write(89,*) t, (nofluxt(ipc),ipc=1,npc), nofluxsum
!        if(mod(istep,intfoc).eq.0.or.istep.eq.nstep) close(89)
!      end if
!    end if


    return


    contains


    subroutine emit_from_rectangular_surf &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,ud,surfp,hl1,hu1,hl2,hu2,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ud,ir
      real(kind=8),intent(in)    :: surfp,hl1,hu1,hl2,hu2
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)

!-------------------- assign position and velocity
      vnorm  = ud*vnormf(addr0)
      vtang1 = vtangf(addr1)
      vtang2 = vtangf(addr2)
      rnorm  = surfp + vnorm*injctisw(is)*dranu(ir+3)
      rtang1 = hl1 + (hu1 - hl1)*dranu(ir+4) + vtang1*injctisw(is)*dranu(ir+3)
      rtang2 = hl2 + (hu2 - hl2)*dranu(ir+5) + vtang2*injctisw(is)*dranu(ir+3)

    return
    end subroutine emit_from_rectangular_surf


    subroutine emit_from_rectangular_surf_with_hole &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,ud,surfp,hl1,hu1,hl2,hu2, &
   &   holel1,holeu1,holel2,holeu2,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ud,ir
      real(kind=8),intent(in)    :: surfp,hl1,hu1,hl2,hu2
      real(kind=8),intent(in)    :: holel1,holeu1,holel2,holeu2
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2, icon
      real(kind=8) :: dran(2)

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)

!-------------------- assign position and velocity
      vnorm  = ud*vnormf(addr0)
      vtang1 = vtangf(addr1)
      vtang2 = vtangf(addr2)
      rnorm  = surfp + vnorm*injctisw(is)*dranu(ir+3)
      do!while(.true.)
        call RANU0(dran,2,icon)
        rtang1 = hl1 + (hu1 - hl1)*dran(1)
        rtang2 = hl2 + (hu2 - hl2)*dran(2)
        if(rtang1.lt.holel1.or.rtang1.gt.holeu1.or. &
       &   rtang2.lt.holel2.or.rtang2.gt.holeu2) then
          rtang1 = rtang1 + vtang1*injctisw(is)*dranu(ir+3)
          rtang2 = rtang2 + vtang2*injctisw(is)*dranu(ir+3)
          exit
        end if
      end do

    return
    end subroutine emit_from_rectangular_surf_with_hole


    subroutine emit_from_circular_surf &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,ud,axis1,axis2,surfp,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ud,ir
      real(kind=8),intent(in)    :: surfp,axis1,axis2
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2,j3
      real(kind=8) :: radial, angular

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)
      j3 = nprs + npr(is)*dranu(ir+3)

!-------------------- assign position and velocity
      vnorm  = ud*vnormf(addr0)
      vtang1 = vtangf(addr1)
      vtang2 = vtangf(addr2)
      radial = radius*runiform(j3)
      angular = pi2*dranu(ir+4)
      rnorm  = surfp + vnorm*injctisw(is)*dranu(ir+5)
      rtang1 = axis1 + radial*dcos(angular) + vtang1*injctisw(is)*dranu(ir+5)
      rtang2 = axis2 + radial*dsin(angular) + vtang2*injctisw(is)*dranu(ir+5)

    return
    end subroutine emit_from_circular_surf


    subroutine emit_from_cylindrical_surf1 &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,psi,axis1,axis2,hl,hu,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ir
      real(kind=8),intent(in)    :: psi,axis1,axis2,hl,hu
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2,j3

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)
      j3 = nprs + npr(is)*dranu(ir+3)

!-------------------- assign position and velocity
      vnorm  = vnormc(addr0)*dcos(psicos(addr2)+psi) &
     &       - vtangc(addr1)*dsin(psicos(addr2)+psi)
      vtang1 = vnormc(addr0)*dsin(psicos(addr2)+psi) &
     &       + vtangc(addr1)*dcos(psicos(addr2)+psi)
      vtang2 = vtangf(j3)
      rnorm  = axis1 &
     &       + radius*dcos(psicos(addr2)+psi) &
     &       + vnorm*injctisw(is)*dranu(ir+4)
      rtang1 = axis2 &
     &       + radius*dsin(psicos(addr2)+psi) &
     &       + vtang1*injctisw(is)*dranu(ir+4)
      rtang2 = hl + (hu - hl)*dranu(ir+5) + vtang2*injctisw(is)*dranu(ir+4)
     
    return
    end subroutine emit_from_cylindrical_surf1


    subroutine emit_from_cylindrical_surf2 &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,psi,axis1,axis2,hl,hu,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ir
      real(kind=8),intent(in)    :: psi,axis1,axis2,hl,hu
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2
      real(kind=8) :: v2normi

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)

!-------------------- assign position and velocity
      vnorm  = vnormc(addr0)*dcos(psi) &
     &       - vtangc(addr1)*dsin(psi)
      vtang1 = vnormc(addr0)*dsin(psi) &
     &       + vtangc(addr1)*dcos(psi)
      vtang2 = vtangf(addr2)
      v2normi = radius/dsqrt(vnorm*vnorm + vtang1*vtang1)
      rnorm  = axis1 &
     &       + vnorm*(v2normi + injctisw(is)*dranu(ir+3))
      rtang1 = axis2 &
     &       + vtang1*(v2normi + injctisw(is)*dranu(ir+3))
      rtang2 = hl + (hu - hl)*dranu(ir+4) + vtang2*injctisw(is)*dranu(ir+3)
     
    return
    end subroutine emit_from_cylindrical_surf2


    subroutine emit_from_cylindrical_surf3 &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,psi,axis1,axis2,hl,hu,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ir
      real(kind=8),intent(in)    :: psi,axis1,axis2,hl,hu
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2,j3
      real(kind=8) :: eazimuth

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      eazimuth = pi2*dranu(ir+2)
      j3 = nprs + npr(is)*dranu(ir+3)

!-------------------- assign position and velocity
      vnorm  = vnormc(addr0)*dcos(eazimuth) &
     &       - vtangc(addr1)*dsin(eazimuth)
      vtang1 = vnormc(addr0)*dsin(eazimuth) &
     &       + vtangc(addr1)*dcos(eazimuth)
      vtang2 = vtangf(j3)
      rnorm  = axis1 &
     &       + radius*dcos(eazimuth) &
     &       + vnorm*injctisw(is)*dranu(ir+4)
      rtang1 = axis2 &
     &       + radius*dsin(eazimuth) &
     &       + vtang1*injctisw(is)*dranu(ir+4)
      rtang2 = hl + (hu - hl)*dranu(ir+5) + vtang2*injctisw(is)*dranu(ir+4)

    return
    end subroutine emit_from_cylindrical_surf3


    subroutine emit_from_spherical_surf1 &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,psi,axis1,axis2,hl,hu,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ir
      real(kind=8),intent(in)    :: psi,axis1,axis2,hl,hu
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2
      real(kind=8) :: v2normi

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)

!-------------------- assign position and velocity
      vnorm  = vnormc(addr0)*dcos(psi) &
     &       - vtangc(addr1)*dsin(psi)
      vtang1 = vnormc(addr0)*dsin(psi) &
     &       + vtangc(addr1)*dcos(psi)
      vtang2 = vtangf(addr2)
      v2normi = radius/dsqrt(vnorm*vnorm + vtang1*vtang1)
      rnorm  = axis1 &
     &       + vnorm*(v2normi + injctisw(is)*dranu(ir+3))
      rtang1 = axis2 &
     &       + vtang1*(v2normi + injctisw(is)*dranu(ir+3))
      rtang2 = hl + (hu - hl)*dranu(ir+4)
     
    return
    end subroutine emit_from_spherical_surf1


  end subroutine
