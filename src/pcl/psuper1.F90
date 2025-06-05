#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"

subroutine psuper1(ustep, func)
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

    integer(kind=4) :: m, l
    integer(kind=8) :: nprs, npre
    integer(kind=8) :: nofluxt(inpc)
    integer(kind=8) :: nofluxsum
    integer(kind=4) :: i, j, k, i1, j1, k1
    integer(kind=4) :: ib, jb, kb, i1b, j1b, k1b
    integer(kind=4) :: ia, ja, ka, i1a, j1a, k1a
    integer(kind=4) :: nns, nne, isw, injctisw(ispec)
    integer(kind=4) :: ie, ipcplane, nemitw
    integer(kind=4) :: iran, iepl, ipc, is, iss, ustep
    integer(kind=4) :: addr0, addr1, addr2
    integer(kind=4) :: icon
    integer(kind=4) :: nemit(inepl), omni
    integer(kind=4) :: xl, xu, yl, yu, zl, zu
    integer(kind=4) :: func
    integer(kind=4) :: idim, rid
    integer(kind=4) :: ierr
!  integer(kind=4) :: oh3_map_region_to_node
    real(kind=8) :: fnemit
    real(kind=8) :: align, radius, axis(2)
    real(kind=8) :: dxl, dxu, dyl, dyu, dzl, dzu
    real(kind=8) :: xlocalb, ylocalb, zlocalb
    real(kind=8) :: xb, yb, zb
    real(kind=8) :: xlocala, ylocala, zlocala
    real(kind=8) :: xa, ya, za
    real(kind=8) :: xr, yr, zr
    real(kind=8) :: xd1, yd1, zd1
    real(kind=8) :: xd2, yd2, zd2
    real(kind=8) :: tslx, tsly, rustep, qs
    real(kind=8) :: vvx1, vvy1, vvz1
    real(kind=8) :: vvx2, vvy2, vvz2
    real(kind=8) :: vx1w1, vx1w2, vx1w3, vx1w4
    real(kind=8) :: vy1w1, vy1w2, vy1w3, vy1w4
    real(kind=8) :: vz1w1, vz1w2, vz1w3, vz1w4
    real(kind=8) :: vx2w1, vx2w2, vx2w3, vx2w4
    real(kind=8) :: vy2w1, vy2w2, vy2w3, vy2w4
    real(kind=8) :: vz2w1, vz2w2, vz2w3, vz2w4
    real(kind=8) :: csz, snz, csxy, snxy, xew, yew, zew
    real(kind=8) :: te11, te12, te13, te21, te22, te23, te31, te32, te33
    real(kind=8) :: arearatio(12)
    real(kind=8) :: tanthetaz, costhetaxy, sinthetaxy
    real(kind=8) :: tantzcostxy, tantzsintxy
    real(kind=8) :: betav, zdepth, rtmp
    real(kind=8) :: xsepa1, ysepa1, xsepa2, ysepa2
    real(kind=8) :: xlocal, ylocal, zlocal
    real(kind=8) :: x1, y1, z1, z2, xy1, xz1, yz1, xz2, yz2
    real(kind=8) :: v1, v2, v3, v4, v5, v6, v7, v8
    ! type(oh_particle) :: pinj
    ! integer(kind=4) :: ninjct

    xl = sdoms(1, 1, sdid(1) + 1); xu = sdoms(2, 1, sdid(1) + 1)
    yl = sdoms(1, 2, sdid(1) + 1); yu = sdoms(2, 2, sdid(1) + 1)
    zl = sdoms(1, 3, sdid(1) + 1); zu = sdoms(2, 3, sdid(1) + 1)
    dxl = xl; dxu = xu; dyl = yl; dyu = yu; dzl = zl; dzu = zu

    tanthetaz = tan(thetaz)
    costhetaxy = cos(thetaxy)
    sinthetaxy = sin(thetaxy)
    tantzcostxy = tanthetaz*costhetaxy
    tantzsintxy = tanthetaz*sinthetaxy

    ! for three dimensional system
    ! zero clear
    aj(:, :, :, :, 1:2) = 0.0d0
    if (func .eq. 1) then
        ajdg(:, :, :, :, 1:2) = 0.0d0
    end if

    ! supervision 1

    ! injection of particles from inner/outer boundaries
    do is = 1, nspec
        injctisw(is) = injct(is)*ustep
    end do

    ! top of species loop 1
    npre = 0
    nne = 0
    rustep = 1.0d0/ustep
    tslx = 2.0d0*slx
    tsly = 2.0d0*sly
    do is = 1, nspec
        nprs = npre + 1
        npre = npre + npr(is)
        nns = nne + 1
        nne = nne + nepl(is)
        ! zero clear work array
        wrk(TJX:TJZ, :, :, :) = 0.0d0

        ! from outer boundary
        if ((nflag_emit(is) .eq. 0) .and. &
            (npbnd(1, is) .eq. 2 .or. &
             npbnd(2, is) .eq. 2 .or. &
             npbnd(3, is) .eq. 2)) then

            ! x-bottom
            if (npbnd(1, is) .eq. 2 .and. bared(1, 1, sdid(1) + 1) .eq. 1) then
                ! number of inject particles
                narrxf(is) = (int((istep + lstep)*arrxf(is), 8) - &
                              int((istep + lstep - injct(is))*arrxf(is), 8))*ustep
                if (narrxf(is) .gt. 0) then
                    call RANU0(dranu, narrxf(is)*5, icon)
                    iran = 0
                    do m = 1, narrxf(is)
                        iran = iran + 5
                        ! pick up from velocity distribution
                        addr1 = nprs + npr(is)*dranu(iran - 4)
                        addr2 = nprs + npr(is)*dranu(iran - 3)
                        ! assign position and velocity
                        pinj(1)%x = vxf(addr1)*injctisw(is)*dranu(iran - 2)
                        pinj(1)%y = yl + (yu - yl)*dranu(iran - 1)
                        pinj(1)%z = zl + (zu - zl)*dranu(iran)
                        if (pinj(1)%z .le. zssurf .or. pinj(1)%z .le. ubhole) cycle
                        pinj(1)%vx = vxf(addr1)
                        pinj(1)%vy = vyr(addr2)
                        pinj(1)%vz = vzr(addr2)
                        xlocalb = -1.0d0
                        ylocalb = pinj(1)%y - pinj(1)%vy*ustep
                        zlocalb = pinj(1)%z - pinj(1)%vz*ustep
                        pinj(1)%nid = sdid(1)
                        pinj(1)%spec = is
                        pinj(1)%preside = 0
                        pinj(1)%pid = 0

                        call oh2_inject_particle(pinj(1))
                        call current_deposition
                    end do
                end if
            end if

            ! x-top
            if (npbnd(1, is) .eq. 2 .and. bared(2, 1, sdid(1) + 1) .eq. 1) then
                narrxb(is) = (int((istep + lstep)*arrxb(is), 8) &
               &            - int((istep + lstep - injct(is))*arrxb(is), 8))*ustep
                if (narrxb(is) .gt. 0) then
                    call RANU0(dranu, narrxb(is)*5, icon)
                    iran = 0
                    do m = 1, narrxb(is)
                        iran = iran + 5
                        ! pick up from velocity distribution
                        addr1 = nprs + npr(is)*dranu(iran - 4)
                        addr2 = nprs + npr(is)*dranu(iran - 3)
                        ! assign position and velocity
                        pinj(1)%x = xu + vxb(addr1)*injctisw(is)*dranu(iran - 2)
                        pinj(1)%y = yl + (yu - yl)*dranu(iran - 1)
                        pinj(1)%z = zl + (zu - zl)*dranu(iran)
                        if (pinj(1)%z .le. zssurf .or. pinj(1)%z .le. ubhole) cycle
                        pinj(1)%vx = vxb(addr1)
                        pinj(1)%vy = vyr(addr2)
                        pinj(1)%vz = vzr(addr2)
                        xlocalb = xu + 0.99999999d0
                        ylocalb = pinj(1)%y - pinj(1)%vy*ustep
                        zlocalb = pinj(1)%z - pinj(1)%vz*ustep
                        pinj(1)%nid = sdid(1)
                        pinj(1)%spec = is
                        pinj(1)%preside = 0
                        pinj(1)%pid = 0

                        call oh2_inject_particle(pinj(1))
                        call current_deposition
                    end do
                end if
            end if

!-------------------- y-bottom
            if (npbnd(2, is) .eq. 2 .and. bared(1, 2, sdid(1) + 1) .eq. 1) then
                narryf(is) = (int((istep + lstep)*arryf(is), 8) &
               &            - int((istep + lstep - injct(is))*arryf(is), 8))*ustep
                if (narryf(is) .gt. 0) then
                    call RANU0(dranu, narryf(is)*5, icon)
                    iran = 0
                    do m = 1, narryf(is)
                        iran = iran + 5
!               ----- pick up from velocity distribution
                        addr1 = nprs + npr(is)*dranu(iran - 4)
                        addr2 = nprs + npr(is)*dranu(iran - 3)
!               ----- assign position and velocity
                        pinj(1)%x = xl + (xu - xl)*dranu(iran - 2)
                        pinj(1)%y = vyf(addr1)*injctisw(is)*dranu(iran - 1)
                        pinj(1)%z = zl + (zu - zl)*dranu(iran)
                        if (pinj(1)%z .le. zssurf .or. pinj(1)%z .le. ubhole) cycle
                        pinj(1)%vx = vxr(addr2)
                        pinj(1)%vy = vyf(addr1)
                        pinj(1)%vz = vzr(addr2)
                        xlocalb = pinj(1)%x - pinj(1)%vx*ustep
                        ylocalb = -1.0d0
                        zlocalb = pinj(1)%z - pinj(1)%vz*ustep
                        pinj(1)%nid = sdid(1)
                        pinj(1)%spec = is
                        pinj(1)%preside = 0
                        pinj(1)%pid = 0
!
                        call oh2_inject_particle(pinj(1))
                        call current_deposition
                    end do
                end if
            end if

!-------------------- y-top
            if (npbnd(2, is) .eq. 2 .and. bared(2, 2, sdid(1) + 1) .eq. 1) then
                narryb(is) = (int((istep + lstep)*arryb(is), 8) &
               &            - int((istep + lstep - injct(is))*arryb(is), 8))*ustep
                if (narryb(is) .gt. 0) then
                    call RANU0(dranu, narryb(is)*5, icon)
                    iran = 0
                    do m = 1, narryb(is)
                        iran = iran + 5
!               ----- pick up from velocity distribution
                        addr1 = nprs + npr(is)*dranu(iran - 4)
                        addr2 = nprs + npr(is)*dranu(iran - 3)
!               ----- assign position and velocity
                        pinj(1)%x = xl + (xu - xl)*dranu(iran - 2)
                        pinj(1)%y = yu + vyb(addr1)*injctisw(is)*dranu(iran - 1)
                        pinj(1)%z = zl + (zu - zl)*dranu(iran)
                        if (pinj(1)%z .le. zssurf .or. pinj(1)%z .le. ubhole) cycle
                        pinj(1)%vx = vxr(addr2)
                        pinj(1)%vy = vyb(addr1)
                        pinj(1)%vz = vzr(addr2)
                        xlocalb = pinj(1)%x - pinj(1)%vx*ustep
                        ylocalb = yu + 0.99999999d0
                        zlocalb = pinj(1)%z - pinj(1)%vz*ustep
                        pinj(1)%nid = sdid(1)
                        pinj(1)%spec = is
                        pinj(1)%preside = 0
                        pinj(1)%pid = 0
!
                        call oh2_inject_particle(pinj(1))
                        call current_deposition
                    end do
                end if
            end if

!-------------------- z-bottom
            if (npbnd(3, is) .eq. 2 .and. bared(1, 3, sdid(1) + 1) .eq. 1) then
                narrzf(is) = (int((istep + lstep)*arrzf(is), 8) &
               &            - int((istep + lstep - injct(is))*arrzf(is), 8))*ustep
                if (narrzf(is) .gt. 0) then
                    call RANU0(dranu, narrzf(is)*5, icon)
                    iran = 0
                    do m = 1, narrzf(is)
                        iran = iran + 5
!               ----- pick up from velocity distribution
                        addr1 = nprs + npr(is)*dranu(iran - 4)
                        addr2 = nprs + npr(is)*dranu(iran - 3)
!               ----- assign position and velocity
                        pinj(1)%x = xl + (xu - xl)*dranu(iran - 2)
                        pinj(1)%y = yl + (yu - yl)*dranu(iran - 1)
                        pinj(1)%z = vzf(addr1)*injctisw(is)*dranu(iran)
                        if (pinj(1)%z .le. zssurf .or. pinj(1)%z .le. ubhole) cycle
                        pinj(1)%vx = vxr(addr2)
                        pinj(1)%vy = vyr(addr2)
                        pinj(1)%vz = vzf(addr1)
                        xlocalb = pinj(1)%x - pinj(1)%vx*ustep
                        ylocalb = pinj(1)%y - pinj(1)%vy*ustep
                        zlocalb = -1.0d0
                        pinj(1)%nid = sdid(1)
                        pinj(1)%spec = is
                        pinj(1)%preside = 0
                        pinj(1)%pid = 0
!
                        call oh2_inject_particle(pinj(1))
                        call current_deposition
                    end do
                end if
            end if

!-------------------- z-top
            if (npbnd(3, is) .eq. 2 .and. bared(2, 3, sdid(1) + 1) .eq. 1) then
                narrzb(is) = (int((istep + lstep)*arrzb(is), 8) &
               &            - int((istep + lstep - injct(is))*arrzb(is), 8))*ustep
                if (narrzb(is) .gt. 0) then
                    call RANU0(dranu, narrzb(is)*5, icon)
                    iran = 0
                    do m = 1, narrzb(is)
                        iran = iran + 5
!               ------- pick up from velocity distribution
                        addr1 = nprs + npr(is)*dranu(iran - 4)
                        addr2 = nprs + npr(is)*dranu(iran - 3)
!               ------- assign position and velocity
                        pinj(1)%x = xl + (xu - xl)*dranu(iran - 2)
                        pinj(1)%y = yl + (yu - yl)*dranu(iran - 1)
                        pinj(1)%z = zu + vzb(addr1)*injctisw(is)*dranu(iran)
                        if (pinj(1)%z .le. zssurf .or. pinj(1)%z .le. ubhole) cycle
                        pinj(1)%vx = vxr(addr2)
                        pinj(1)%vy = vyr(addr2)
                        pinj(1)%vz = vzb(addr1)
                        xlocalb = pinj(1)%x - pinj(1)%vx*ustep
                        ylocalb = pinj(1)%y - pinj(1)%vy*ustep
                        zlocalb = zu + 0.99999999d0
                        pinj(1)%nid = sdid(1)
                        pinj(1)%spec = is
                        pinj(1)%preside = 0
                        pinj(1)%pid = 0
!
                        call oh2_inject_particle(pinj(1))
                        call current_deposition
                    end do
                end if
            end if

!============================== from inner boundary
        else if (nflag_emit(is) .ge. 1 .and. nflag_emit(is) .le. 2) then
            if (injct(is) .eq. 0) cycle
            ie = mod(istep, injct(is))
            if (ie .ne. 0 .or. inpf(is) .eq. 0) cycle
!
!------------------------------ top plane loop
            do iepl = nns, 0
!            if(nflag_emit(is).eq.1) then
                nemit(iepl) = (int((istep + lstep)*fluxf(iepl), 8) &
               &            - int((istep + lstep - injct(is))*fluxf(iepl), 8))*ustep
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
                if (nemitw .gt. 0) then
                    call RANU0(dranu, nemitw*6, icon)
                    if (icon .ne. 0) print *, "[psuper] icon=", icon
                end if
                iran = 0
!
!           ------------------- pos. & vel. assignment
                if (geotype(ipcplane) .eq. 0 .or. geotype(ipcplane) .eq. 1) then
                    do m = 1, nemit(iepl)
                        if (nemd(iepl) .eq. 1) then
                            call emit_from_rectangular_surf &
                           &  (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, &
                           &   xlocalb, ylocalb, zlocalb, +1, &
                           &   peject(iepl)%grd, peject(iepl)%yl, peject(iepl)%yu, &
                           &   peject(iepl)%zl, peject(iepl)%zu, m*6 - 5)
                        else if (nemd(iepl) .eq. -1) then
                            call emit_from_rectangular_surf &
                           &  (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, &
                           &   xlocalb, ylocalb, zlocalb, -1, &
                           &   peject(iepl)%grd, peject(iepl)%yl, peject(iepl)%yu, &
                           &   peject(iepl)%zl, peject(iepl)%zu, m*6 - 5)
                        else if (nemd(iepl) .eq. 2) then
                            call emit_from_rectangular_surf &
                           &  (pinj(1)%vy, pinj(1)%vz, pinj(1)%vx, pinj(1)%y, pinj(1)%z, pinj(1)%x, &
                           &   ylocalb, zlocalb, xlocalb, +1, &
                           &   peject(iepl)%grd, peject(iepl)%zl, peject(iepl)%zu, &
                           &   peject(iepl)%xl, peject(iepl)%xu, m*6 - 5)
                        else if (nemd(iepl) .eq. -2) then
                            call emit_from_rectangular_surf &
                           &  (pinj(1)%vy, pinj(1)%vz, pinj(1)%vx, pinj(1)%y, pinj(1)%z, pinj(1)%x, &
                           &   ylocalb, zlocalb, xlocalb, -1, &
                           &   peject(iepl)%grd, peject(iepl)%zl, peject(iepl)%zu, &
                           &   peject(iepl)%xl, peject(iepl)%xu, m*6 - 5)
                        else if (nemd(iepl) .eq. 3) then
                            call emit_from_rectangular_surf &
                           &  (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, &
                           &   zlocalb, xlocalb, ylocalb, +1, &
                           &   peject(iepl)%grd, peject(iepl)%xl, peject(iepl)%xu, &
                           &   peject(iepl)%yl, peject(iepl)%yu, m*6 - 5)
                        else if (nemd(iepl) .eq. -3) then
                            call emit_from_rectangular_surf &
                           &  (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, &
                           &   zlocalb, xlocalb, ylocalb, -1, &
                           &   peject(iepl)%grd, peject(iepl)%xl, peject(iepl)%xu, &
                           &   peject(iepl)%yl, peject(iepl)%yu, m*6 - 5)
                        end if
!
                        if (pinj(1)%x .lt. xl .or. pinj(1)%x .ge. xu .or. &
                       &   pinj(1)%y .lt. yl .or. pinj(1)%y .ge. yu .or. &
                       &   pinj(1)%z .lt. zl .or. pinj(1)%z .ge. zu) then
                            rid = oh3_map_particle_to_subdomain(pinj(1)%x, pinj(1)%y, pinj(1)%z)
                            pinj(1)%nid = rid
                        else
                            pinj(1)%nid = sdid(1)
                        end if
                        pinj(1)%spec = is
                        pinj(1)%pid = ipcplane
!
                        call oh2_inject_particle(pinj(1))
                        call current_deposition
                    end do
                else if (geotype(ipcplane) .eq. 2) then
                    align = cylinder(ipcplane)%align
                    radius = cylinder(ipcplane)%radius
                    if (abs(nemd(iepl)) .eq. align) then
                        do m = 1, nemit(iepl)
                            if (nemd(iepl) .eq. 1) then
                                call emit_from_circular_surf &
                               &  (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, &
                               &   xlocalb, ylocalb, zlocalb, +1, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%grd, m*6 - 5)
                            else if (nemd(iepl) .eq. -1) then
                                call emit_from_circular_surf &
                               &  (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, &
                               &   xlocalb, ylocalb, zlocalb, -1, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%grd, m*6 - 5)
                            else if (nemd(iepl) .eq. 2) then
                                call emit_from_circular_surf &
                               &  (pinj(1)%vy, pinj(1)%vz, pinj(1)%vx, pinj(1)%y, pinj(1)%z, pinj(1)%x, &
                               &   ylocalb, zlocalb, xlocalb, +1, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%grd, m*6 - 5)
                            else if (nemd(iepl) .eq. -2) then
                                call emit_from_circular_surf &
                               &  (pinj(1)%vy, pinj(1)%vz, pinj(1)%vx, pinj(1)%y, pinj(1)%z, pinj(1)%x, &
                               &   ylocalb, zlocalb, xlocalb, -1, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%grd, m*6 - 5)
                            else if (nemd(iepl) .eq. 3) then
                                call emit_from_circular_surf &
                               &  (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, &
                               &   zlocalb, xlocalb, ylocalb, +1, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%grd, m*6 - 5)
                            else if (nemd(iepl) .eq. -3) then
                                call emit_from_circular_surf &
                               &  (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, &
                               &   zlocalb, xlocalb, ylocalb, -1, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%grd, m*6 - 5)
                            end if
!
                            if (pinj(1)%x .lt. xl .or. pinj(1)%x .ge. xu .or. &
                           &   pinj(1)%y .lt. yl .or. pinj(1)%y .ge. yu .or. &
                           &   pinj(1)%z .lt. zl .or. pinj(1)%z .ge. zu) then
                                rid = oh3_map_particle_to_subdomain(pinj(1)%x, pinj(1)%y, pinj(1)%z)
                                pinj(1)%nid = rid
                            else
                                pinj(1)%nid = sdid(1)
                            end if
                            pinj(1)%spec = is
                            pinj(1)%pid = ipcplane
!
                            call oh2_inject_particle(pinj(1))
                            call current_deposition
                        end do
                    else
                        do m = 1, nemit(iepl)
                            if (align .eq. 1 .and. omni .eq. 0 .and. radius .ge. 1.0d0) then
                                call emit_from_cylindrical_surf1 &
                               &  (pinj(1)%vz, pinj(1)%vy, pinj(1)%vx, pinj(1)%z, pinj(1)%y, pinj(1)%x, &
                               &   zlocalb, ylocalb, xlocalb, psizy, &
                               &   cylinder(ipcplane)%axis(2), cylinder(ipcplane)%axis(1), &
                               &   peject(iepl)%xl, peject(iepl)%xu, m*6 - 5)
                            else if (align .eq. 1 .and. omni .eq. 0 .and. radius .lt. 1.0d0) then
                                call emit_from_cylindrical_surf2 &
                               &  (pinj(1)%vz, pinj(1)%vy, pinj(1)%vx, pinj(1)%z, pinj(1)%y, pinj(1)%x, &
                               &   zlocalb, ylocalb, xlocalb, psizy, &
                               &   cylinder(ipcplane)%axis(2), cylinder(ipcplane)%axis(1), &
                               &   peject(iepl)%xl, peject(iepl)%xu, m*5 - 4)
                            else if (align .eq. 1 .and. omni .eq. 1) then
                                call emit_from_cylindrical_surf3 &
                               &  (pinj(1)%vz, pinj(1)%vy, pinj(1)%vx, pinj(1)%z, pinj(1)%y, pinj(1)%x, &
                               &   zlocalb, ylocalb, xlocalb, psizy, &
                               &   cylinder(ipcplane)%axis(2), cylinder(ipcplane)%axis(1), &
                               &   peject(iepl)%xl, peject(iepl)%xu, m*6 - 5)
                            else if (align .eq. 2 .and. omni .eq. 0 .and. radius .ge. 1.0d0) then
                                call emit_from_cylindrical_surf1 &
                               &  (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, &
                               &   zlocalb, xlocalb, ylocalb, psizx, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%yl, peject(iepl)%yu, m*6 - 5)
                            else if (align .eq. 2 .and. omni .eq. 0 .and. radius .lt. 1.0d0) then
                                call emit_from_cylindrical_surf2 &
                               &  (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, &
                               &   zlocalb, xlocalb, ylocalb, psizx, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%yl, peject(iepl)%yu, m*5 - 4)
                            else if (align .eq. 2 .and. omni .eq. 1) then
                                call emit_from_cylindrical_surf3 &
                               &  (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, &
                               &   zlocalb, xlocalb, ylocalb, psizx, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%yl, peject(iepl)%yu, m*6 - 5)
                            else if (align .eq. 3 .and. omni .eq. 0 .and. radius .ge. 1.0d0) then
                                call emit_from_cylindrical_surf1 &
                               &  (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, &
                               &   xlocalb, ylocalb, zlocalb, psixy, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%zl, peject(iepl)%zu, m*6 - 5)
                            else if (align .eq. 3 .and. omni .eq. 0 .and. radius .lt. 1.0d0) then
                                call emit_from_cylindrical_surf2 &
                               &  (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, &
                               &   xlocalb, ylocalb, zlocalb, psixy, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%zl, peject(iepl)%zu, m*5 - 4)
                            else if (align .eq. 3 .and. omni .eq. 1) then
                                call emit_from_cylindrical_surf3 &
                               &  (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, &
                               &   xlocalb, ylocalb, zlocalb, psixy, &
                               &   cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                               &   peject(iepl)%zl, peject(iepl)%zu, m*6 - 5)
                            end if
!
                            if (pinj(1)%x .lt. xl .or. pinj(1)%x .ge. xu .or. &
                           &   pinj(1)%y .lt. yl .or. pinj(1)%y .ge. yu .or. &
                           &   pinj(1)%z .lt. zl .or. pinj(1)%z .ge. zu) then
                                rid = oh3_map_particle_to_subdomain(pinj(1)%x, pinj(1)%y, pinj(1)%z)
                                pinj(1)%nid = rid
                            else
                                pinj(1)%nid = sdid(1)
                            end if
                            pinj(1)%spec = is
                            pinj(1)%pid = ipcplane
!
                            call oh2_inject_particle(pinj(1))
                            call current_deposition
                        end do
                    end if
!            else if(geotype(ipcplane).eq.3) then
!              !under construction!
                end if
!
!           --------- calculate accumulated charge
                gcount(1)%chgacm(:, ipcplane) = gcount(1)%chgacm(:, ipcplane) &
               &                             - q(is)*nemitw*sqdscaled(ipcplane)
                gcount(1)%outflux(ipcplane, is) = gcount(1)%outflux(ipcplane, is) &
               &                               + nemitw

            end do

!-------------------- from land surface

!-------------------- from hole basement

!-------------------- from hole flank

        end if
!============================== end of plane loop

!---------------- store current value from work to ajx
        qs = q(is)
        do k = -1, zu - zl + 1
            do j = -1, yu - yl + 1
                do i = -1, xu - xl + 1
                    aj(JX, i, j, k, 1) = aj(JX, i, j, k, 1) + wrk(TJX, i, j, k)*qs
                    aj(JY, i, j, k, 1) = aj(JY, i, j, k, 1) + wrk(TJY, i, j, k)*qs
                    aj(JZ, i, j, k, 1) = aj(JZ, i, j, k, 1) + wrk(TJZ, i, j, k)*qs
                end do
            end do
        end do
        if (func .eq. 1) then
            iss = (is - 1)*3
            do k = -1, zu - zl + 1
                do j = -1, yu - yl + 1
                    do i = -1, xu - xl + 1
                        ajdg(iss + JX, i, j, k, 1) = ajdg(iss + JX, i, j, k, 1) &
                      &                     + wrk(TJX, i, j, k)*qs
                        ajdg(iss + JY, i, j, k, 1) = ajdg(iss + JY, i, j, k, 1) &
                      &                     + wrk(TJY, i, j, k)*qs
                        ajdg(iss + JZ, i, j, k, 1) = ajdg(iss + JZ, i, j, k, 1) &
                      &                     + wrk(TJZ, i, j, k)*qs
                    end do
                end do
            end do
        end if

    end do

    ! end of species loop 1

    ! other treatments
    ! charge transfer simulating bias current
    if (myid .eq. 0) then
        i = 1
        do while (biasc(i)%to .ne. 0)
            if (biasc(i)%to .gt. 0) then
                if (biasc(i)%from .gt. 0) &
               &  gcount(1)%chgacm(:, biasc(i)%from) = &
               &  gcount(1)%chgacm(:, biasc(i)%from) - biasc(i)%val*ustep
                gcount(1)%chgacm(:, biasc(i)%to) = &
               &  gcount(1)%chgacm(:, biasc(i)%to) + biasc(i)%val*ustep
            end if
            i = i + 1
        end do
    end if

    return

contains

    subroutine emit_from_rectangular_surf &
   &  (vnorm, vtang1, vtang2, rnorm, rtang1, rtang2, rnorm_o, rtang1_o, rtang2_o, &
   &   ud, surfp, hl1, hu1, hl2, hu2, ir)
!-------------------- args & vars
        integer(kind=4), intent(in) :: ud, ir
        real(kind=8), intent(in)    :: surfp, hl1, hu1, hl2, hu2
        real(kind=8), intent(out)   :: vnorm, vtang1, vtang2
        real(kind=8), intent(out)   :: rnorm, rtang1, rtang2
        real(kind=8), intent(out)   :: rnorm_o, rtang1_o, rtang2_o
        integer(kind=4) :: addr0, addr1, addr2

!-------------------- pick up from velocity distribution
        addr0 = nprs + npr(is)*dranu(ir)
        addr1 = nprs + npr(is)*dranu(ir + 1)
        addr2 = nprs + npr(is)*dranu(ir + 2)

!-------------------- assign position and velocity
        vnorm = ud*vnormf(addr0)
        vtang1 = vtangf(addr1)
        vtang2 = vtangf(addr2)
        rnorm_o = surfp + ud*1e-3
        rtang1_o = hl1 + (hu1 - hl1)*dranu(ir + 4)
        rtang2_o = hl2 + (hu2 - hl2)*dranu(ir + 5)
        rnorm = rnorm_o + vnorm*injctisw(is)*dranu(ir + 3)
        rtang1 = rtang1_o + vtang1*injctisw(is)*dranu(ir + 3)
        rtang2 = rtang2_o + vtang2*injctisw(is)*dranu(ir + 3)

        return
    end subroutine emit_from_rectangular_surf

    subroutine emit_from_rectangular_surf_with_hole &
   &  (vnorm, vtang1, vtang2, rnorm, rtang1, rtang2, rnorm_o, rtang1_o, rtang2_o, &
   &   ud, surfp, hl1, hu1, hl2, hu2, &
   &   holel1, holeu1, holel2, holeu2, ir)
!-------------------- args & vars
        integer(kind=4), intent(in) :: ud, ir
        real(kind=8), intent(in)    :: surfp, hl1, hu1, hl2, hu2
        real(kind=8), intent(in)    :: holel1, holeu1, holel2, holeu2
        real(kind=8), intent(out)   :: vnorm, vtang1, vtang2
        real(kind=8), intent(out)   :: rnorm, rtang1, rtang2
        real(kind=8), intent(out)   :: rnorm_o, rtang1_o, rtang2_o
        integer(kind=4) :: addr0, addr1, addr2, icon
        real(kind=8) :: dran(2)

!-------------------- pick up from velocity distribution
        addr0 = nprs + npr(is)*dranu(ir)
        addr1 = nprs + npr(is)*dranu(ir + 1)
        addr2 = nprs + npr(is)*dranu(ir + 2)

!-------------------- assign position and velocity
        vnorm = ud*vnormf(addr0)
        vtang1 = vtangf(addr1)
        vtang2 = vtangf(addr2)
        rnorm_o = surfp
        rnorm = rnorm_o + vnorm*injctisw(is)*dranu(ir + 3)
        do!while(.true.)
            call RANU0(dran, 2, icon)
            rtang1_o = hl1 + (hu1 - hl1)*dran(1)
            rtang2_o = hl2 + (hu2 - hl2)*dran(2)
            if (rtang1_o .lt. holel1 .or. rtang1_o .gt. holeu1 .or. &
           &   rtang2_o .lt. holel2 .or. rtang2_o .gt. holeu2) then
                rtang1 = rtang1_o + vtang1*injctisw(is)*dranu(ir + 3)
                rtang2 = rtang2_o + vtang2*injctisw(is)*dranu(ir + 3)
                exit
            end if
        end do

        return
    end subroutine emit_from_rectangular_surf_with_hole

    subroutine emit_from_circular_surf &
   &  (vnorm, vtang1, vtang2, rnorm, rtang1, rtang2, rnorm_o, rtang1_o, rtang2_o, &
   &   ud, axis1, axis2, surfp, ir)
!-------------------- args & vars
        integer(kind=4), intent(in) :: ud, ir
        real(kind=8), intent(in)    :: surfp, axis1, axis2
        real(kind=8), intent(out)   :: vnorm, vtang1, vtang2
        real(kind=8), intent(out)   :: rnorm, rtang1, rtang2
        real(kind=8), intent(out)   :: rnorm_o, rtang1_o, rtang2_o
        integer(kind=4) :: addr0, addr1, addr2, j3
        real(kind=8) :: radial, angular

!-------------------- pick up from velocity distribution
        addr0 = nprs + npr(is)*dranu(ir)
        addr1 = nprs + npr(is)*dranu(ir + 1)
        addr2 = nprs + npr(is)*dranu(ir + 2)
        j3 = nprs + npr(is)*dranu(ir + 3)

!-------------------- assign position and velocity
        vnorm = ud*vnormf(addr0)
        vtang1 = vtangf(addr1)
        vtang2 = vtangf(addr2)
        radial = radius*runiform(j3)
        angular = pi2*dranu(ir + 4)
        rnorm_o = surfp
        rtang1_o = axis1 + radial*cos(angular)
        rtang2_o = axis2 + radial*sin(angular)
        rnorm = rnorm_o + vnorm*injctisw(is)*dranu(ir + 5)
        rtang1 = rtang1_o + vtang1*injctisw(is)*dranu(ir + 5)
        rtang2 = rtang2_o + vtang2*injctisw(is)*dranu(ir + 5)

        return
    end subroutine emit_from_circular_surf

    subroutine emit_from_cylindrical_surf1 &
   &  (vnorm, vtang1, vtang2, rnorm, rtang1, rtang2, rnorm_o, rtang1_o, rtang2_o, &
   &   psi, axis1, axis2, hl, hu, ir)
!-------------------- args & vars
        integer(kind=4), intent(in) :: ir
        real(kind=8), intent(in)    :: psi, axis1, axis2, hl, hu
        real(kind=8), intent(out)   :: vnorm, vtang1, vtang2
        real(kind=8), intent(out)   :: rnorm, rtang1, rtang2
        real(kind=8), intent(out)   :: rnorm_o, rtang1_o, rtang2_o
        integer(kind=4) :: addr0, addr1, addr2, j3

!-------------------- pick up from velocity distribution
        addr0 = nprs + npr(is)*dranu(ir)
        addr1 = nprs + npr(is)*dranu(ir + 1)
        addr2 = nprs + npr(is)*dranu(ir + 2)
        j3 = nprs + npr(is)*dranu(ir + 3)

!-------------------- assign position and velocity
        vnorm = vnormc(addr0)*cos(psicos(addr2) + psi) &
       &         - vtangc(addr1)*sin(psicos(addr2) + psi)
        vtang1 = vnormc(addr0)*sin(psicos(addr2) + psi) &
       &         + vtangc(addr1)*cos(psicos(addr2) + psi)
        vtang2 = vtangf(j3)
        rnorm_o = axis1 &
       &         + radius*cos(psicos(addr2) + psi)
        rtang1_o = axis2 &
       &         + radius*sin(psicos(addr2) + psi)
        rtang2_o = hl + (hu - hl)*dranu(ir + 5)
        rnorm = rnorm_o + vnorm*injctisw(is)*dranu(ir + 4)
        rtang1 = rtang1_o + vtang1*injctisw(is)*dranu(ir + 4)
        rtang2 = rtang2_o + vtang2*injctisw(is)*dranu(ir + 4)

        return
    end subroutine emit_from_cylindrical_surf1

    subroutine emit_from_cylindrical_surf2 &
   &  (vnorm, vtang1, vtang2, rnorm, rtang1, rtang2, rnorm_o, rtang1_o, rtang2_o, &
   &   psi, axis1, axis2, hl, hu, ir)
!-------------------- args & vars
        integer(kind=4), intent(in) :: ir
        real(kind=8), intent(in)    :: psi, axis1, axis2, hl, hu
        real(kind=8), intent(out)   :: vnorm, vtang1, vtang2
        real(kind=8), intent(out)   :: rnorm, rtang1, rtang2
        real(kind=8), intent(out)   :: rnorm_o, rtang1_o, rtang2_o
        integer(kind=4) :: addr0, addr1, addr2
        real(kind=8) :: v2normi

!-------------------- pick up from velocity distribution
        addr0 = nprs + npr(is)*dranu(ir)
        addr1 = nprs + npr(is)*dranu(ir + 1)
        addr2 = nprs + npr(is)*dranu(ir + 2)

!-------------------- assign position and velocity
        vnorm = vnormc(addr0)*cos(psi) &
       &         - vtangc(addr1)*sin(psi)
        vtang1 = vnormc(addr0)*sin(psi) &
       &         + vtangc(addr1)*cos(psi)
        vtang2 = vtangf(addr2)
        v2normi = radius/dsqrt(vnorm*vnorm + vtang1*vtang1)
        rnorm_o = axis1 + vnorm*v2normi
        rtang1_o = axis2 + vtang1*v2normi
        rtang2_o = hl + (hu - hl)*dranu(ir + 4)
        rnorm = rnorm_o + vnorm*injctisw(is)*dranu(ir + 3)
        rtang1 = rtang1_o + vtang1*injctisw(is)*dranu(ir + 3)
        rtang2 = rtang2_o + vtang2*injctisw(is)*dranu(ir + 3)

        return
    end subroutine emit_from_cylindrical_surf2

    subroutine emit_from_cylindrical_surf3 &
   &  (vnorm, vtang1, vtang2, rnorm, rtang1, rtang2, rnorm_o, rtang1_o, rtang2_o, &
   &   psi, axis1, axis2, hl, hu, ir)
!-------------------- args & vars
        integer(kind=4), intent(in) :: ir
        real(kind=8), intent(in)    :: psi, axis1, axis2, hl, hu
        real(kind=8), intent(out)   :: vnorm, vtang1, vtang2
        real(kind=8), intent(out)   :: rnorm, rtang1, rtang2
        real(kind=8), intent(out)   :: rnorm_o, rtang1_o, rtang2_o
        integer(kind=4) :: addr0, addr1, addr2, j3
        real(kind=8) :: eazimuth

!-------------------- pick up from velocity distribution
        addr0 = nprs + npr(is)*dranu(ir)
        addr1 = nprs + npr(is)*dranu(ir + 1)
        eazimuth = pi2*dranu(ir + 2)
        j3 = nprs + npr(is)*dranu(ir + 3)

!-------------------- assign position and velocity
        vnorm = vnormc(addr0)*dcos(eazimuth) &
       &         - vtangc(addr1)*dsin(eazimuth)
        vtang1 = vnormc(addr0)*dsin(eazimuth) &
       &         + vtangc(addr1)*dcos(eazimuth)
        vtang2 = vtangf(j3)
        rnorm_o = axis1 &
       &         + radius*dcos(eazimuth)
        rtang1_o = axis2 &
       &         + radius*dsin(eazimuth)
        rtang2_o = hl + (hu - hl)*dranu(ir + 5)
        rnorm = rnorm_o &
       &       + vnorm*injctisw(is)*dranu(ir + 4)
        rtang1 = rtang1_o &
       &       + vtang1*injctisw(is)*dranu(ir + 4)
        rtang2 = rtang2_o + vtang2*injctisw(is)*dranu(ir + 4)

        return
    end subroutine emit_from_cylindrical_surf3

    subroutine emit_from_spherical_surf1 &
   &  (vnorm, vtang1, vtang2, rnorm, rtang1, rtang2, psi, axis1, axis2, hl, hu, ir)
!-------------------- args & vars
        integer(kind=4), intent(in) :: ir
        real(kind=8), intent(in)    :: psi, axis1, axis2, hl, hu
        real(kind=8), intent(out)   :: vnorm, vtang1, vtang2, rnorm, rtang1, rtang2
        integer(kind=4) :: addr0, addr1, addr2
        real(kind=8) :: v2normi

        ! pick up from velocity distribution
        addr0 = nprs + npr(is)*dranu(ir)
        addr1 = nprs + npr(is)*dranu(ir + 1)
        addr2 = nprs + npr(is)*dranu(ir + 2)

!-------------------- assign position and velocity
        vnorm = vnormc(addr0)*dcos(psi) &
       &       - vtangc(addr1)*dsin(psi)
        vtang1 = vnormc(addr0)*dsin(psi) &
       &       + vtangc(addr1)*dcos(psi)
        vtang2 = vtangf(addr2)
        v2normi = radius/dsqrt(vnorm*vnorm + vtang1*vtang1)
        rnorm = axis1 &
       &       + vnorm*(v2normi + injctisw(is)*dranu(ir + 3))
        rtang1 = axis2 &
       &       + vtang1*(v2normi + injctisw(is)*dranu(ir + 3))
        rtang2 = hl + (hu - hl)*dranu(ir + 4)
        rnorm = axis1 &
       &       + vnorm*(v2normi + injctisw(is)*dranu(ir + 3))
        rtang1 = axis2 &
       &       + vtang1*(v2normi + injctisw(is)*dranu(ir + 3))
        rtang2 = hl + (hu - hl)*dranu(ir + 4)

        return
    end subroutine emit_from_spherical_surf1

    subroutine current_deposition
        implicit none

        xlocalb = xlocalb - xl
        ylocalb = ylocalb - yl
        zlocalb = zlocalb - zl

        ib = floor(xlocalb)
        jb = floor(ylocalb)
        kb = floor(zlocalb)

        xb = ib
        yb = jb
        zb = kb

        i1b = ib + 1
        j1b = jb + 1
        k1b = kb + 1

        xlocala = pinj(1)%x - xl
        ylocala = pinj(1)%y - yl
        zlocala = pinj(1)%z - zl

        ia = floor(xlocala)
        ja = floor(ylocala)
        ka = floor(zlocala)

        xa = ia
        ya = ja
        za = ka

        i1a = ia + 1
        j1a = ja + 1
        k1a = ka + 1

        if (ib .eq. ia) then
            xr = (xlocalb + xlocala)*0.5d0
        else
            xr = max(xb, xa)
        end if

        if (jb .eq. ja) then
            yr = (ylocalb + ylocala)*0.5d0
        else
            yr = max(yb, ya)
        end if

        if (kb .eq. ka) then
            zr = (zlocalb + zlocala)*0.5d0
        else
            zr = max(zb, za)
        end if

        xd1 = (xr + xlocalb)*0.5d0 - xb
        yd1 = (yr + ylocalb)*0.5d0 - yb
        zd1 = (zr + zlocalb)*0.5d0 - zb

        vvx1 = (xr - xlocalb)*rustep
        vvy1 = (yr - ylocalb)*rustep
        vvz1 = (zr - zlocalb)*rustep

        vx1w4 = yd1*zd1
        vx1w2 = yd1 - vx1w4
        vx1w3 = zd1 - vx1w4
        vx1w1 = 1.0d0 - yd1 - vx1w3
        vy1w4 = xd1*zd1
        vy1w2 = xd1 - vy1w4
        vy1w3 = zd1 - vy1w4
        vy1w1 = 1.0d0 - xd1 - vy1w3
        vz1w4 = xd1*yd1
        vz1w2 = xd1 - vz1w4
        vz1w3 = yd1 - vz1w4
        vz1w1 = 1.0d0 - xd1 - vz1w3

        xd2 = (xlocala + xr)*0.5d0 - xa
        yd2 = (ylocala + yr)*0.5d0 - ya
        zd2 = (zlocala + zr)*0.5d0 - za

        vvx2 = (xlocala - xr)*rustep
        vvy2 = (ylocala - yr)*rustep
        vvz2 = (zlocala - zr)*rustep

        vx2w4 = yd2*zd2
        vx2w2 = yd2 - vx2w4
        vx2w3 = zd2 - vx2w4
        vx2w1 = 1.0d0 - yd2 - vx2w3
        vy2w4 = xd2*zd2
        vy2w2 = xd2 - vy2w4
        vy2w3 = zd2 - vy2w4
        vy2w1 = 1.0d0 - xd2 - vy2w3
        vz2w4 = xd2*yd2
        vz2w2 = xd2 - vz2w4
        vz2w3 = yd2 - vz2w4
        vz2w1 = 1.0d0 - xd2 - vz2w3

        wrk(TJX, ib, jb, kb) = wrk(TJX, ib, jb, kb) + vvx1*vx1w1
        wrk(TJY, ib, jb, kb) = wrk(TJY, ib, jb, kb) + vvy1*vy1w1
        wrk(TJZ, ib, jb, kb) = wrk(TJZ, ib, jb, kb) + vvz1*vz1w1

        wrk(TJX, ib, j1b, kb) = wrk(TJX, ib, j1b, kb) + vvx1*vx1w2
        wrk(TJY, i1b, jb, kb) = wrk(TJY, i1b, jb, kb) + vvy1*vy1w2
        wrk(TJZ, i1b, jb, kb) = wrk(TJZ, i1b, jb, kb) + vvz1*vz1w2

        wrk(TJX, ib, jb, k1b) = wrk(TJX, ib, jb, k1b) + vvx1*vx1w3
        wrk(TJY, ib, jb, k1b) = wrk(TJY, ib, jb, k1b) + vvy1*vy1w3
        wrk(TJZ, ib, j1b, kb) = wrk(TJZ, ib, j1b, kb) + vvz1*vz1w3

        wrk(TJX, ib, j1b, k1b) = wrk(TJX, ib, j1b, k1b) + vvx1*vx1w4
        wrk(TJY, i1b, jb, k1b) = wrk(TJY, i1b, jb, k1b) + vvy1*vy1w4
        wrk(TJZ, i1b, j1b, kb) = wrk(TJZ, i1b, j1b, kb) + vvz1*vz1w4

        wrk(TJX, ia, ja, ka) = wrk(TJX, ia, ja, ka) + vvx2*vx2w1
        wrk(TJY, ia, ja, ka) = wrk(TJY, ia, ja, ka) + vvy2*vy2w1
        wrk(TJZ, ia, ja, ka) = wrk(TJZ, ia, ja, ka) + vvz2*vz2w1

        wrk(TJX, ia, j1a, ka) = wrk(TJX, ia, j1a, ka) + vvx2*vx2w2
        wrk(TJY, i1a, ja, ka) = wrk(TJY, i1a, ja, ka) + vvy2*vy2w2
        wrk(TJZ, i1a, ja, ka) = wrk(TJZ, i1a, ja, ka) + vvz2*vz2w2

        wrk(TJX, ia, ja, k1a) = wrk(TJX, ia, ja, k1a) + vvx2*vx2w3
        wrk(TJY, ia, ja, k1a) = wrk(TJY, ia, ja, k1a) + vvy2*vy2w3
        wrk(TJZ, ia, j1a, ka) = wrk(TJZ, ia, j1a, ka) + vvz2*vz2w3

        wrk(TJX, ia, j1a, k1a) = wrk(TJX, ia, j1a, k1a) + vvx2*vx2w4
        wrk(TJY, i1a, ja, k1a) = wrk(TJY, i1a, ja, k1a) + vvy2*vy2w4
        wrk(TJZ, i1a, j1a, ka) = wrk(TJZ, i1a, j1a, ka) + vvz2*vz2w4

        return
    end subroutine current_deposition

end subroutine psuper1
