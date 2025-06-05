#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"

subroutine psolve1(ps, ustep, func)
!   ____________________________________________________________
!
!               S U B R O U T I N E   P S O L V E
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a current distribution to     .
!   .      grid points from particle locations by the          .
!   .      charge conservation method based on the "zigzag"    .
!   .      scheme. [Umeda, T., 2003]                           .
!   ............................................................
!
!
!-------------------- parameter and common block
    use oh_type
    use paramt
    use allcom
    ! use detect_collision
    use m_boundary_base
    use particle_collision
    use m_ohhinfo
    use m_str
    use m_logging
    implicit none

    !> 1: Primary subdomain
    !> 2: Secondary subdomain
    integer, intent(in) :: ps
    !> Simulation time width
    integer, intent(in) :: ustep
    !> True if you want to create data for output
    integer, intent(in) :: func

    integer(kind=8) :: m, ns, ne, nee
    integer(kind=8) :: isfluxt(inpc, ispec)
    integer(kind=8) :: influxsum
    integer(kind=4) :: i, j, k
    integer(kind=4) :: ib, jb, kb, i1b, j1b, k1b
    integer(kind=4) :: ia, ja, ka, i1a, j1a, k1a
    integer(kind=4) :: is, iss, ibdy, ipc, itch
    integer(kind=4) :: nb(3)
    integer(kind=4) :: rid, nud
    integer(kind=4) :: icon
    real(kind=8) :: xlc(inpc), ylc(inpc), zlc(inpc)
    real(kind=8) :: xuc(inpc), yuc(inpc), zuc(inpc)
    real(kind=8) :: xlocalb, ylocalb, zlocalb
    real(kind=8) :: xb, yb, zb
    real(kind=8) :: xlocala, ylocala, zlocala
    real(kind=8) :: xa, ya, za
    real(kind=8) :: xr, yr, zr
    real(kind=8) :: xd1, yd1, zd1, xd2, yd2, zd2
    real(kind=8) :: x1, y1, z1, z2, xy1, xz1, yz1, xz2, yz2
    real(kind=8) :: tslx, tsly, tslz, gustep, rustep, qs, qmp
    real(kind=8) :: v1, v2, v3, v4, v5, v6, v7, v8
    real(kind=8) :: eex, eey, eez, bbx, bby, bbz, boris, vxt, vyt, vzt, vxyz
    real(kind=8) :: vvx1, vvy1, vvz1, vvx2, vvy2, vvz2
    real(kind=8) :: vx1w1, vx1w2, vx1w3, vx1w4, vx2w1, vx2w2, vx2w3, vx2w4
    real(kind=8) :: vy1w1, vy1w2, vy1w3, vy1w4, vy2w1, vy2w2, vy2w3, vy2w4
    real(kind=8) :: vz1w1, vz1w2, vz1w3, vz1w4, vz2w1, vz2w2, vz2w3, vz2w4

    real(kind=8) :: displx, disply, displz, displ3i

    real(kind=8) :: disp1, disp2, disp3
    real(kind=8) :: radsq(inpc), rbwlsq, rdomsq
    real(kind=8) :: xpast0, ypast0, zpast0, xpast1, ypast1, zpast1
    real(kind=8) :: xmove, ymove, zmove, xsepa, ysepa, zsepa, termA, termB
    real(kind=8) :: tfrac, tfrac1, tfrac2
    real(kind=8) :: xpast1a, ypast1a, xpast1b, ypast1b
    real(kind=8) :: xpast2a, ypast2a, xpast2b, ypast2b
    real(kind=8) :: nsigmadt, colfp, uppb
    real(kind=8) :: ewx, ewy, tew

    real(kind=8) :: wrk_e(EX:EZ, 8)

    type(t_CollisionRecord) :: record
    type(t_BoundaryList) :: boundaries

    logical :: is_boundary_condition_applied

    call ohhinfo_update(sdoms(:, :, sdid(ps) + 1))
    do ipc = 1, npc
        xlc(ipc) = xlpc(ipc) - dxl; xuc(ipc) = xupc(ipc) - dxl
        ylc(ipc) = ylpc(ipc) - dyl; yuc(ipc) = yupc(ipc) - dyl
        zlc(ipc) = zlpc(ipc) - dzl; zuc(ipc) = zupc(ipc) - dzl
    end do

    boundaries = create_simple_collision_boundaries(sdoms(1:2, 1:3, sdid(ps) + 1))

    if (ps .eq. 1) sey(:, :) = 0.0d0

    do ipc = 1, npc
        if (geotype(ipc) .eq. 2) then
            radsq(ipc) = cylinder(ipc)%radius*cylinder(ipc)%radius
        else if (geotype(ipc) .eq. 3) then
            radsq(ipc) = sphere(ipc)%radius*sphere(ipc)%radius
        end if
    end do
    if (rbowl .gt. 0.0d0) then
        rbwlsq = rbowl*rbowl
    else
        rbwlsq = 0.0d0
    end if
    if (rdome .gt. 0.0d0) then
        rdomsq = rdome*rdome
    else
        rdomsq = 0.0d0
    end if

    ! species loop
    nee = pbase(ps)
    rustep = 1.0d0/ustep
    tslx = 2.0d0*slx
    tsly = 2.0d0*sly
    tslz = 2.0d0*slz
    do is = 1, nspec
        call logging_debuglog('Psolve species loop: '//str(is))
        ns = nee + 1
        ne = nee + totalp(is, ps)
        nee = ne
        if (ewmodel .eq. 1 .or. ewmodel .eq. 2) then
            tew = t - dt*nretard
            ewx = Ew(1, is, 1)*cos(omegaw(1)*tew)
            ewy = Ew(2, is, 1)*sin(omegaw(1)*tew)
        else
            ewx = 0.0d0
            ewy = 0.0d0
        end if
!       -------------
        if (mfpath(is) > 1.0d0) then
            nsigmadt = ustep/mfpath(is)
        else if (mfpath(is) < -1.0d0) then
            nsigmadt = ustep/mfpath(is)
        else
            nsigmadt = 0.0d0
        end if
        uppb = 1.0d0 - exp(-vxyzmax(is)*abs(nsigmadt))
        ! zero clear work array
        wrk(TJX:TJZ, :, :, :) = 0.0d0

        ! inner loop
        do m = ns, ne
            if (pbuf(m)%preside < 0) then
                nphgram(pbuf(m)%nid + 1, is, ps) = &
                    nphgram(pbuf(m)%nid + 1, is, ps) - 1
                pbuf(m)%nid = -1
                pbuf(m)%preside = 0
            else if (pbuf(m)%preside == OH_PCL_INJECTED) then
                call RANU0(dranu, 1, icon)
                gustep = ustep*dranu(1)
                qmp = qm(is)*gustep*0.5d0
                pbuf(m)%preside = 0
            else ! if pbuf(m)%preside == OH_PCL_ALIVE
                gustep = ustep
                qmp = qm(is)*gustep*0.5d0
            end if

            if (pbuf(m)%nid == -1) cycle

            is_boundary_condition_applied = .false.

            xlocalb = pbuf(m)%x - dxl
            ylocalb = pbuf(m)%y - dyl
            zlocalb = pbuf(m)%z - dzl

            ib = floor(xlocalb)
            jb = floor(ylocalb)
            kb = floor(zlocalb)

            xb = ib
            yb = jb
            zb = kb

            i1b = ib + 1
            j1b = jb + 1
            k1b = kb + 1

            x1 = xlocalb - xb
            y1 = ylocalb - yb
            z1 = zlocalb - zb

            xy1 = x1*y1
            xz1 = x1*z1
            yz1 = y1*z1
            z2 = 1.0d0 - z1
            xz2 = x1*z2
            yz2 = y1*z2
            v3 = xy1*z1
            v2 = xz1 - v3
            v4 = yz1 - v3
            v1 = z1 - xz1 - v4
            v7 = xy1*z2
            v6 = xz2 - v7
            v8 = yz2 - v7
            v5 = z2 - xz2 - v8

            wrk_e(EX:EZ, 1) = wrk(EX:EZ, ib, jb, k1b)
            wrk_e(EX:EZ, 2) = wrk(EX:EZ, i1b, jb, k1b)
            wrk_e(EX:EZ, 3) = wrk(EX:EZ, i1b, j1b, k1b)
            wrk_e(EX:EZ, 4) = wrk(EX:EZ, ib, j1b, k1b)
            wrk_e(EX:EZ, 5) = wrk(EX:EZ, ib, jb, kb)
            wrk_e(EX:EZ, 6) = wrk(EX:EZ, i1b, jb, kb)
            wrk_e(EX:EZ, 7) = wrk(EX:EZ, i1b, j1b, kb)
            wrk_e(EX:EZ, 8) = wrk(EX:EZ, ib, j1b, kb)

            ! eb-fields interpolation
            eex = wrk_e(EX, 1)*v1 + wrk_e(EX, 2)*v2 &
                  + wrk_e(EX, 3)*v3 + wrk_e(EX, 4)*v4 &
                  + wrk_e(EX, 5)*v5 + wrk_e(EX, 6)*v6 &
                  + wrk_e(EX, 7)*v7 + wrk_e(EX, 8)*v8
            eey = wrk_e(EY, 1)*v1 + wrk_e(EY, 2)*v2 &
                  + wrk_e(EY, 3)*v3 + wrk_e(EY, 4)*v4 &
                  + wrk_e(EY, 5)*v5 + wrk_e(EY, 6)*v6 &
                  + wrk_e(EY, 7)*v7 + wrk_e(EY, 8)*v8
            eez = wrk_e(EZ, 1)*v1 + wrk_e(EZ, 2)*v2 &
                  + wrk_e(EZ, 3)*v3 + wrk_e(EZ, 4)*v4 &
                  + wrk_e(EZ, 5)*v5 + wrk_e(EZ, 6)*v6 &
                  + wrk_e(EZ, 7)*v7 + wrk_e(EZ, 8)*v8
            bbx = wrk(BX, ib, jb, k1b)*v1 + wrk(BX, i1b, jb, k1b)*v2 &
                  + wrk(BX, i1b, j1b, k1b)*v3 + wrk(BX, ib, j1b, k1b)*v4 &
                  + wrk(BX, ib, jb, kb)*v5 + wrk(BX, i1b, jb, kb)*v6 &
                  + wrk(BX, i1b, j1b, kb)*v7 + wrk(BX, ib, j1b, kb)*v8
            bby = wrk(BY, ib, jb, k1b)*v1 + wrk(BY, i1b, jb, k1b)*v2 &
                  + wrk(BY, i1b, j1b, k1b)*v3 + wrk(BY, ib, j1b, k1b)*v4 &
                  + wrk(BY, ib, jb, kb)*v5 + wrk(BY, i1b, jb, kb)*v6 &
                  + wrk(BY, i1b, j1b, kb)*v7 + wrk(BY, ib, j1b, kb)*v8
            bbz = wrk(BZ, ib, jb, k1b)*v1 + wrk(BZ, i1b, jb, k1b)*v2 &
                  + wrk(BZ, i1b, j1b, k1b)*v3 + wrk(BZ, ib, j1b, k1b)*v4 &
                  + wrk(BZ, ib, jb, kb)*v5 + wrk(BZ, i1b, jb, kb)*v6 &
                  + wrk(BZ, i1b, j1b, kb)*v7 + wrk(BZ, ib, j1b, kb)*v8

            do itch = 1, ntch
                displx = pbuf(m)%x - rtch(1, itch)
                disply = pbuf(m)%y - rtch(2, itch)
                displz = pbuf(m)%z - rtch(3, itch)
                displ3i = displx*displx + disply*disply + displz*displz
                displ3i = 1.0d0/((displ3i + r2cutoff(itch))*sqrt(displ3i))
                eex = eex + e1tch(itch)*displx*displ3i
                eey = eey + e1tch(itch)*disply*displ3i
                eez = eez + e1tch(itch)*displz*displ3i
            end do

            ! charge-to-mass ratio is taken into consideration
            eex = (eex + ewx)*qmp
            eey = (eey + ewy)*qmp
            eez = eez*qmp
            bbx = bbx*qmp
            bby = bby*qmp
            bbz = bbz*qmp

            ! update particle velocities (Buneman-Boris method)
            boris = 2.0d0/(1.0d0 + bbx*bbx + bby*bby + bbz*bbz)

            pbuf(m)%vx = pbuf(m)%vx + eex
            pbuf(m)%vy = pbuf(m)%vy + eey
            pbuf(m)%vz = pbuf(m)%vz + eez

            vxt = pbuf(m)%vx + pbuf(m)%vy*bbz - pbuf(m)%vz*bby
            vyt = pbuf(m)%vy + pbuf(m)%vz*bbx - pbuf(m)%vx*bbz
            vzt = pbuf(m)%vz + pbuf(m)%vx*bby - pbuf(m)%vy*bbx

            pbuf(m)%vx = pbuf(m)%vx + boris*(vyt*bbz - vzt*bby)
            pbuf(m)%vy = pbuf(m)%vy + boris*(vzt*bbx - vxt*bbz)
            pbuf(m)%vz = pbuf(m)%vz + boris*(vxt*bby - vyt*bbx)

            pbuf(m)%vx = pbuf(m)%vx + eex
            pbuf(m)%vy = pbuf(m)%vy + eey
            pbuf(m)%vz = pbuf(m)%vz + eez

            ! velocity modification associated with collisions
            colfp = colf(1, ib, jb, k1b, ps)*v1 + colf(1, i1b, jb, k1b, ps)*v2 &
                    + colf(1, i1b, j1b, k1b, ps)*v3 + colf(1, ib, j1b, k1b, ps)*v4 &
                    + colf(1, ib, jb, kb, ps)*v5 + colf(1, i1b, jb, kb, ps)*v6 &
                    + colf(1, i1b, j1b, kb, ps)*v7 + colf(1, ib, j1b, kb, ps)*v8
            if ((nsigmadt .gt. 0.0d0) .and. &
                (pbuf(m)%x .le. xlcol(1) .or. pbuf(m)%x .ge. xucol(1) .or. &
                 pbuf(m)%y .le. ylcol(1) .or. pbuf(m)%y .ge. yucol(1) .or. &
                 pbuf(m)%z .le. zlcol(1) .or. pbuf(m)%z .ge. zucol(1))) then
                vxt = pbuf(m)%vx - vdtx(is)
                vyt = pbuf(m)%vy - vdty(is)
                vzt = pbuf(m)%vz - vdtz(is)
                vxyz = sqrt(vxt*vxt + vyt*vyt + vzt*vzt)
                call RANU0(dranu, 1, icon)
                if (exp(-vxyz*nsigmadt) .lt. dranu(1)*colfp) then
                    vxyz = 0.0d0
                    call RANU0(dranu, 1, icon)
                    do while (exp(-vxyz*nsigmadt) .ge. 1.0d0 - dranu(1)*uppb)
                        call RANN0(0.0d0, peth(is), dranu, 1, icon)
                        vxt = dranu(1)
                        call RANN0(0.0d0, peth(is), dranu, 1, icon)
                        vyt = dranu(1)
                        call RANN0(0.0d0, path(is), dranu, 1, icon)
                        vzt = dranu(1)
                        vxyz = sqrt(vxt*vxt + vyt*vyt + vzt*vzt)
                        call RANU0(dranu, 1, icon)
                    end do
                    pbuf(m)%vx = vxt + vdtx(is)
                    pbuf(m)%vy = vyt + vdty(is)
                    pbuf(m)%vz = vzt + vdtz(is)
                end if
            else if ((nsigmadt .lt. 0.0d0) .and. &
                     (pbuf(m)%x .le. xlcol(1) .or. pbuf(m)%x .ge. xucol(1) .or. &
                      pbuf(m)%y .le. ylcol(1) .or. pbuf(m)%y .ge. yucol(1) .or. &
                      pbuf(m)%z .le. zlcol(1) .or. pbuf(m)%z .ge. zucol(1))) then
                vxt = peth(is)
                vyt = peth(is)
                vzt = path(is)
                vxyz = sqrt(vxt*vxt + vyt*vyt + vzt*vzt)
                call RANU0(dranu, 1, icon)
                if (exp(vxyz*nsigmadt) .lt. dranu(1)*colfp) then
                    call RANN0(0.0d0, peth(is), dranu, 1, icon)
                    vxt = dranu(1)
                    call RANN0(0.0d0, peth(is), dranu, 1, icon)
                    vyt = dranu(1)
                    call RANN0(0.0d0, path(is), dranu, 1, icon)
                    vzt = dranu(1)
                    pbuf(m)%vx = vxt + vdtx(is)
                    pbuf(m)%vy = vyt + vdty(is)
                    pbuf(m)%vz = vzt + vdtz(is)
                end if
            end if

            ! update particle positions
            xmove = pbuf(m)%vx*gustep
            ymove = pbuf(m)%vy*gustep
            zmove = pbuf(m)%vz*gustep
            pbuf(m)%x = pbuf(m)%x + xmove
            pbuf(m)%y = pbuf(m)%y + ymove
            pbuf(m)%z = pbuf(m)%z + zmove

            ! internal boundary treatment 0
            block
                double precision :: p1(3), p2(3), v(3)
                integer :: iboundary
                type(t_CollisionRecord) :: tmp_record

                v(:) = [pbuf(m)%vx, pbuf(m)%vy, pbuf(m)%vz]
                p2(:) = [pbuf(m)%x, pbuf(m)%y, pbuf(m)%z]
                p1(:) = p2(:) - v(:)*gustep
                
                record%is_collided = .false.
                record%t = 100.0d0  ! "0 < t < 1" should be established when collided.
                do iboundary = 1, boundaries%nboundaries
                    tmp_record = boundaries%boundaries(iboundary)%check_collision(p1, p2)
                    if (tmp_record%is_collided .and. tmp_record%t < record%t) then
                        record = tmp_record
                    end if
                end do

                ! record = boundaries%check_collision(p1(:), p2(:))

                if (record%is_collided) then
                    pbuf(m)%x = record%position(1)
                    pbuf(m)%y = record%position(2)
                    pbuf(m)%z = record%position(3)
                    pbuf(m)%preside = OH_PCL_TO_BE_ACCUMULATED
                end if
            end block
! #include "defsurf.fnc"

            ! external boundary treatment
            nud = 14
            if (pbuf(m)%x .lt. dxl) then
                nud = nud - 1
            else if (pbuf(m)%x .ge. dxu) then
                nud = nud + 1
            end if
            if (pbuf(m)%y .lt. dyl) then
                nud = nud - 3
            else if (pbuf(m)%y .ge. dyu) then
                nud = nud + 3
            end if
            if (pbuf(m)%z .lt. dzl) then
                nud = nud - 9
            else if (pbuf(m)%z .ge. dzu) then
                nud = nud + 9
            end if
            if (nud .eq. 14) then
                xlocala = pbuf(m)%x - dxl
                ylocala = pbuf(m)%y - dyl
                zlocala = pbuf(m)%z - dzl

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
            else
                rid = nborps(nud, is, ps)
                pbuf(m)%nid = rid
                nphgram(sdid(ps) + 1, is, ps) = nphgram(sdid(ps) + 1, is, ps) - 1

                if (pbuf(m)%x .lt. 0.0d0) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(1, is) .eq. 0) then
                        xlocala = pbuf(m)%x
                        pbuf(m)%x = pbuf(m)%x + slx
                        ia = -1
                        xa = -1.0d0
                        i1a = 0
                        xr = 0.0d0
                    else if (npbnd(1, is) .eq. 1) then
                        pbuf(m)%x = -pbuf(m)%x
                        pbuf(m)%vx = -pbuf(m)%vx
                        xlocala = pbuf(m)%x
                        ia = 0
                        xa = 0.0d0
                        i1a = 1
                        xr = 0.0d0
                    else
                        xlocala = -1.0d0
                        ia = -1
                        xa = -1.0d0
                        i1a = 0
                        xr = 0.0d0
                    end if
                else if (pbuf(m)%x .ge. slx) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(1, is) .eq. 0) then
                        xlocala = pbuf(m)%x - dxl
                        pbuf(m)%x = pbuf(m)%x - slx
                        ia = ngx
                        xa = dngx
                        i1a = ngx + 1
                        xr = dngx
                    else if (npbnd(1, is) .eq. 1) then
                        pbuf(m)%x = tslx - pbuf(m)%x
                        pbuf(m)%vx = -pbuf(m)%vx
                        xlocala = pbuf(m)%x - dxl
                        ia = ngx - 1
                        xa = dngx - 1.0d0
                        i1a = ngx
                        xr = dngx
                    else
                        xlocala = dngx + 1.0d0
                        ia = ngx
                        xa = dngx
                        i1a = ngx + 1
                        xr = dngx
                    end if
                else
                    xlocala = pbuf(m)%x - dxl
                    ia = floor(xlocala)
                    xa = ia
                    i1a = ia + 1
                    if (ib .eq. ia) then
                        xr = (xlocalb + xlocala)*0.5d0
                    else
                        xr = max(xb, xa)
                    end if
                end if

                if (pbuf(m)%y .lt. 0.0d0) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(2, is) .eq. 0) then
                        ylocala = pbuf(m)%y
                        pbuf(m)%y = pbuf(m)%y + sly
                        ja = -1
                        ya = -1.0d0
                        j1a = 0
                        yr = 0.0d0
                    else if (npbnd(2, is) .eq. 1) then
                        pbuf(m)%y = -pbuf(m)%y
                        pbuf(m)%vy = -pbuf(m)%vy
                        ylocala = pbuf(m)%y
                        ja = 0
                        ya = 0.0d0
                        j1a = 1
                        yr = 0.0d0
                    else
                        ylocala = -1.0d0
                        ja = -1
                        ya = -1.0d0
                        j1a = 0
                        yr = 0.0d0
                    end if
                else if (pbuf(m)%y .ge. sly) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(2, is) .eq. 0) then
                        ylocala = pbuf(m)%y - dyl
                        pbuf(m)%y = pbuf(m)%y - sly
                        ja = ngy
                        ya = dngy
                        j1a = ngy + 1
                        yr = dngy
                    else if (npbnd(2, is) .eq. 1) then
                        pbuf(m)%y = tsly - pbuf(m)%y
                        pbuf(m)%vy = -pbuf(m)%vy
                        ylocala = pbuf(m)%y - dyl
                        ja = ngy - 1
                        ya = dngy - 1.0d0
                        j1a = ngy
                        yr = dngy
                    else
                        ylocala = dngy + 1.0d0
                        ja = ngy
                        ya = dngy
                        j1a = ngy + 1
                        yr = dngy
                    end if
                else
                    ylocala = pbuf(m)%y - dyl
                    ja = floor(ylocala)
                    ya = ja
                    j1a = ja + 1
                    if (jb .eq. ja) then
                        yr = (ylocalb + ylocala)*0.5d0
                    else
                        yr = max(yb, ya)
                    end if
                end if

                if (pbuf(m)%z .lt. 0.0d0) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(3, is) .eq. 0) then
                        zlocala = pbuf(m)%z
                        pbuf(m)%z = pbuf(m)%z + slz
                        ka = -1
                        za = -1.0d0
                        k1a = 0
                        zr = 0.0d0
                    else if (npbnd(3, is) .eq. 1) then
                        pbuf(m)%z = -pbuf(m)%z
                        pbuf(m)%vz = -pbuf(m)%vz
                        zlocala = pbuf(m)%z
                        ka = 0
                        za = 0.0d0
                        k1a = 1
                        zr = 0.0d0
                    else
                        zlocala = -1.0d0
                        ka = -1
                        za = -1.0d0
                        k1a = 0
                        zr = 0.0d0
                    end if
                else if (pbuf(m)%z .ge. slz) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(3, is) .eq. 0) then
                        zlocala = pbuf(m)%z - dzl
                        pbuf(m)%z = pbuf(m)%z - slz
                        ka = ngz
                        za = dngz
                        k1a = ngz + 1
                        zr = dngz
                    else if (npbnd(3, is) .eq. 1) then
                        pbuf(m)%z = tslz - pbuf(m)%z
                        pbuf(m)%vz = -pbuf(m)%vz
                        zlocala = pbuf(m)%z - dzl
                        ka = ngz - 1
                        za = dngz - 1.0d0
                        k1a = ngz
                        zr = dngz
                    else
                        zlocala = dngz + 1.0d0
                        ka = ngz
                        za = dngz
                        k1a = ngz + 1
                        zr = dngz
                    end if
                else
                    zlocala = pbuf(m)%z - dzl
                    ka = floor(zlocala)
                    za = ka
                    k1a = ka + 1
                    if (kb .eq. ka) then
                        zr = (zlocalb + zlocala)*0.5d0
                    else
                        zr = max(zb, za)
                    end if
                end if

                if (pbuf(m)%nid .ne. -1) then
                    nphgram(pbuf(m)%nid + 1, is, ps) = nphgram(pbuf(m)%nid + 1, is, ps) + 1
                else
                    gcount(1)%nesc(is) = gcount(1)%nesc(is) + 1
                end if
            end if

            ! If some particle boundary condition are applied and the position of the particle has moved,
            ! perform the collision detection with the internal boundary again.
            if ((pbuf(m)%preside /= OH_PCL_TO_BE_ACCUMULATED) .and. is_boundary_condition_applied) then
                block
                    double precision :: p1(3), p2(3), v(3)
                    integer :: iboundary
                    type(t_CollisionRecord) :: tmp_record

                    v(:) = [pbuf(m)%vx, pbuf(m)%vy, pbuf(m)%vz]
                    p2(:) = [pbuf(m)%x, pbuf(m)%y, pbuf(m)%z]
                    p1(:) = p2(:) - v(:)*gustep
                    
                    record%is_collided = .false.
                    record%t = 100.0d0  ! "0 < t < 1" should be established when collided.
                    do iboundary = 1, boundaries%nboundaries
                        tmp_record = boundaries%boundaries(iboundary)%check_collision(p1, p2)
                        if (tmp_record%is_collided .and. tmp_record%t < record%t) then
                            record = tmp_record
                        end if
                    end do

                    ! record = boundaries%check_collision(p1(:), p2(:))
                    if (record%is_collided) then
                        pbuf(m)%x = record%position(1)
                        pbuf(m)%y = record%position(2)
                        pbuf(m)%z = record%position(3)
                        pbuf(m)%preside = OH_PCL_TO_BE_ACCUMULATED
                    end if
                end block
            end if

            ! internal boundary treatment 1
            if (line_mode .ne. 1) then
                do ipc = 1, npc
                    if (pbuf(m)%nid .eq. -1) cycle
                    if (geotype(ipc) .eq. 0 .or. geotype(ipc) .eq. 1) then
                        if (xlocalb .ge. xlc(ipc) .and. xlocalb .le. xuc(ipc) .and. &
                       &   ylocalb .ge. ylc(ipc) .and. ylocalb .le. yuc(ipc) .and. &
                       &   zlocalb .ge. zlc(ipc) .and. zlocalb .le. zuc(ipc)) then
                            xlocalb = xr
                            ylocalb = yr
                            zlocalb = zr
                        end if
                        if (xlocala .ge. xlc(ipc) .and. xlocala .le. xuc(ipc) .and. &
                       &   ylocala .ge. ylc(ipc) .and. ylocala .le. yuc(ipc) .and. &
                       &   zlocala .ge. zlc(ipc) .and. zlocala .le. zuc(ipc)) then
                            nphgram(pbuf(m)%nid + 1, is, ps) = nphgram(pbuf(m)%nid + 1, is, ps) - 1
                            pbuf(m)%nid = -1
                            gcount(1)%influx(ipc, is) = gcount(1)%influx(ipc, is) + 1
                            gcount(1)%infhist(pbuf(m)%pid, ipc, is) = &
                           &  gcount(1)%infhist(pbuf(m)%pid, ipc, is) + 1
                            xlocala = xr
                            ylocala = yr
                            zlocala = zr
                        end if
                    else if (geotype(ipc) .eq. 2) then
                        if (cylinder(ipc)%align .eq. 1) then
                            disp1 = pbuf(m)%y - cylinder(ipc)%axis(1)
                            disp2 = pbuf(m)%z - cylinder(ipc)%axis(2)
                            if (pbuf(m)%x .ge. cylinder(ipc)%edge(1) .and. pbuf(m)%x .le. cylinder(ipc)%edge(2) .and. &
                           &   disp1*disp1 + disp2*disp2 .lt. radsq(ipc)) then
                                nphgram(pbuf(m)%nid + 1, is, ps) = &
                               &  nphgram(pbuf(m)%nid + 1, is, ps) - 1
                                pbuf(m)%nid = -1
                                gcount(1)%influx(ipc, is) = gcount(1)%influx(ipc, is) + 1
                                gcount(1)%infhist(pbuf(m)%pid, ipc, is) = &
                               &  gcount(1)%infhist(pbuf(m)%pid, ipc, is) + 1
                            end if
                        else if (cylinder(ipc)%align .eq. 2) then
                            disp1 = pbuf(m)%z - cylinder(ipc)%axis(1)
                            disp2 = pbuf(m)%x - cylinder(ipc)%axis(2)
                            if (pbuf(m)%y .ge. cylinder(ipc)%edge(1) .and. pbuf(m)%y .le. cylinder(ipc)%edge(2) .and. &
                           &   disp1*disp1 + disp2*disp2 .lt. radsq(ipc)) then
                                nphgram(pbuf(m)%nid + 1, is, ps) = &
                              &  nphgram(pbuf(m)%nid + 1, is, ps) - 1
                                pbuf(m)%nid = -1
                                gcount(1)%influx(ipc, is) = gcount(1)%influx(ipc, is) + 1
                                gcount(1)%infhist(pbuf(m)%pid, ipc, is) = &
                               &  gcount(1)%infhist(pbuf(m)%pid, ipc, is) + 1
                            end if
                        else if (cylinder(ipc)%align .eq. 3) then
                            disp1 = pbuf(m)%x - cylinder(ipc)%axis(1)
                            disp2 = pbuf(m)%y - cylinder(ipc)%axis(2)
                            if (pbuf(m)%z .ge. cylinder(ipc)%edge(1) .and. pbuf(m)%z .le. cylinder(ipc)%edge(2) .and. &
                           &   disp1*disp1 + disp2*disp2 .lt. radsq(ipc)) then
                                nphgram(pbuf(m)%nid + 1, is, ps) = &
                               &  nphgram(pbuf(m)%nid + 1, is, ps) - 1
                                pbuf(m)%nid = -1
                                gcount(1)%influx(ipc, is) = gcount(1)%influx(ipc, is) + 1
                                gcount(1)%infhist(pbuf(m)%pid, ipc, is) = &
                               &  gcount(1)%infhist(pbuf(m)%pid, ipc, is) + 1
                            end if
                        end if
                    else if (geotype(ipc) .eq. 3) then
                        disp1 = pbuf(m)%x - sphere(ipc)%center(1)
                        disp2 = pbuf(m)%y - sphere(ipc)%center(2)
                        disp3 = pbuf(m)%z - sphere(ipc)%center(3)
                        if (disp1*disp1 + disp2*disp2 + disp3*disp3 .lt. radsq(ipc)) then
                            nphgram(pbuf(m)%nid + 1, is, ps) = &
                           &  nphgram(pbuf(m)%nid + 1, is, ps) - 1
                            pbuf(m)%nid = -1
                            gcount(1)%influx(ipc, is) = gcount(1)%influx(ipc, is) + 1
                            gcount(1)%infhist(pbuf(m)%pid, ipc, is) = &
                           &  gcount(1)%infhist(pbuf(m)%pid, ipc, is) + 1
                        end if
                    end if
                end do
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
        end do

        ! store current value from work to ajx
        qs = q(is)
        do k = -1, zu - zl + 1
        do j = -1, yu - yl + 1
        do i = -1, xu - xl + 1
            aj(JX, i, j, k, ps) = aj(JX, i, j, k, ps) + wrk(TJX, i, j, k)*qs
            aj(JY, i, j, k, ps) = aj(JY, i, j, k, ps) + wrk(TJY, i, j, k)*qs
            aj(JZ, i, j, k, ps) = aj(JZ, i, j, k, ps) + wrk(TJZ, i, j, k)*qs
        end do
        end do
        end do
        if (func .eq. 1) then
            iss = (is - 1)*3
            do k = -1, zu - zl + 1
            do j = -1, yu - yl + 1
            do i = -1, xu - xl + 1
                ajdg(iss + JX, i, j, k, ps) = ajdg(iss + JX, i, j, k, ps) &
               &                      + wrk(TJX, i, j, k)*qs
                ajdg(iss + JY, i, j, k, ps) = ajdg(iss + JY, i, j, k, ps) &
               &                      + wrk(TJY, i, j, k)*qs
                ajdg(iss + JZ, i, j, k, ps) = ajdg(iss + JZ, i, j, k, ps) &
               &                      + wrk(TJZ, i, j, k)*qs
            end do
            end do
            end do
        end if

    end do

    do is = 1, nspec
        do ipc = 1, npc
            if (ncond(ipc) .gt. 0) then
                gcount(1)%chgacm(:, ipc) = gcount(1)%chgacm(:, ipc) &
               &  + gcount(1)%influx(ipc, is)*q(is)*sqdscaled(ipc)
            end if
        end do
    end do

    call boundaries%destroy

    return

contains

    function sey_model_1(delta, efrac, costh)
        real(kind=8) :: sey_model_1
        real(kind=8), intent(in) :: delta, efrac, costh

        sey_model_1 = delta*efrac &
                      *exp(2.0d0 - 2.0d0*sqrt(efrac)) &
                      *exp(2.0d0*(1.0d0 - costh))

        return
    end function sey_model_1

    function sey_model_2(delta, efrac, costh)
        real(kind=8) :: sey_model_2
        real(kind=8), intent(in) :: delta, efrac, costh

        sey_model_2 = 1.114d0*delta/costh*efrac**(-0.35d0) &
                      *(1.0d0 - exp(-2.28d0*costh*efrac**1.35d0))

        return
    end function sey_model_2

end subroutine

subroutine add_boundary_current(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!              A D D _ B O U N D A R Y _ C U R R E N T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a charge distribution to      .
!   .      grid points from particle locations by the          .
!   .      "area-sharing" or "area-weighting" scheme.          .
!   ............................................................
!
!
!-------------------- parameter and common block
    use oh_type
    use paramt
    use allcom
    implicit none

    integer(kind=4) :: i, j, k
    integer(kind=4) :: xu, yu, zu
    integer(kind=4) :: sl, su, dl, du, nl, nu
    integer(kind=4) :: ps

!--------------------
    xu = sdoms(2, 1, sdid(ps) + 1) - sdoms(1, 1, sdid(ps) + 1)
    yu = sdoms(2, 2, sdid(ps) + 1) - sdoms(1, 2, sdid(ps) + 1)
    zu = sdoms(2, 3, sdid(ps) + 1) - sdoms(1, 3, sdid(ps) + 1)

!--------------------
    sl = ctypes(2, 2, 1, CAJ)        !=-4
    su = ctypes(2, 1, 1, CAJ)        !=+2

    nl = ctypes(3, 2, 1, CAJ)        != 3
    nu = ctypes(3, 1, 1, CAJ)        != 3

    dl = sl + nl                !=-1
    du = su - nu                !=-1

!--------------------
    if (bounds(1, 3, sdid(ps) + 1) .eq. 1) then
        do k = 0, nl - 1                !(i.e., do k=0,2)
        do j = 0, yu + (su + nu - sl) - 1        !(i.e., do j=0,yu+8)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
            aj(JX, sl + i, sl + j, dl + k, ps) = aj(JX, sl + i, sl + j, dl + k, ps) &
           &                         + aj(JX, sl + i, sl + j, sl + k, ps)
            aj(JY, sl + i, sl + j, dl + k, ps) = aj(JY, sl + i, sl + j, dl + k, ps) &
           &                         + aj(JY, sl + i, sl + j, sl + k, ps)
            aj(JZ, sl + i, sl + j, dl + k, ps) = aj(JZ, sl + i, sl + j, dl + k, ps) &
           &                         + aj(JZ, sl + i, sl + j, sl + k, ps)
        end do
        end do
        end do
    else if (nfbnd(3) .eq. 1) then
        do k = 0, 0
        do j = 0, yu + (su + nu - sl) - 1        !(i.e., do j=0,yu+8)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
            aj(JX, sl + i, sl + j, k, ps) = aj(JX, sl + i, sl + j, k, ps)*2.0d0
            aj(JY, sl + i, sl + j, k, ps) = aj(JY, sl + i, sl + j, k, ps)*2.0d0
            aj(JZ, sl + i, sl + j, k - 1, ps) = -aj(JZ, sl + i, sl + j, -k, ps)
        end do
        end do
        end do
    end if

!--------------------
    if (bounds(2, 3, sdid(ps) + 1) .eq. 1) then
        do k = 0, nu - 1                !(i.e., do k=0,2)
        do j = 0, yu + (su + nu - sl) - 1        !(i.e., do j=0,yu+8)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
            aj(JX, sl + i, sl + j, zu + du + k, ps) = aj(JX, sl + i, sl + j, zu + du + k, ps) &
           &                            + aj(JX, sl + i, sl + j, zu + su + k, ps)
            aj(JY, sl + i, sl + j, zu + du + k, ps) = aj(JY, sl + i, sl + j, zu + du + k, ps) &
           &                            + aj(JY, sl + i, sl + j, zu + su + k, ps)
            aj(JZ, sl + i, sl + j, zu + du + k, ps) = aj(JZ, sl + i, sl + j, zu + du + k, ps) &
           &                            + aj(JZ, sl + i, sl + j, zu + su + k, ps)
        end do
        end do
        end do
    else if (nfbnd(3) .eq. 1) then
        do k = 0, 0
        do j = 0, yu + (su + nu - sl) - 1        !(i.e., do j=0,yu+8)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
            aj(JX, sl + i, sl + j, zu + k, ps) = aj(JX, sl + i, sl + j, zu + k, ps)*2.0d0
            aj(JY, sl + i, sl + j, zu + k, ps) = aj(JY, sl + i, sl + j, zu + k, ps)*2.0d0
            aj(JZ, sl + i, sl + j, zu + k, ps) = -aj(JZ, sl + i, sl + j, zu - k - 1, ps)
        end do
        end do
        end do
        do k = 1, 1
        do j = 0, yu + (su + nu - sl) - 1        !(i.e., do j=0,yu+8)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
            aj(JX, sl + i, sl + j, zu + k, ps) = aj(JX, sl + i, sl + j, zu - k, ps)
            aj(JY, sl + i, sl + j, zu + k, ps) = aj(JY, sl + i, sl + j, zu - k, ps)
            aj(JZ, sl + i, sl + j, zu + k, ps) = -aj(JZ, sl + i, sl + j, zu - k - 1, ps)
        end do
        end do
        end do
    end if

!--------------------
    if (bounds(1, 2, sdid(ps) + 1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, nl - 1                !(i.e., do j=0,2)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
            aj(JX, sl + i, dl + j, dl + k, ps) = aj(JX, sl + i, dl + j, dl + k, ps) &
           &                         + aj(JX, sl + i, sl + j, dl + k, ps)
            aj(JY, sl + i, dl + j, dl + k, ps) = aj(JY, sl + i, dl + j, dl + k, ps) &
           &                         + aj(JY, sl + i, sl + j, dl + k, ps)
            aj(JZ, sl + i, dl + j, dl + k, ps) = aj(JZ, sl + i, dl + j, dl + k, ps) &
           &                         + aj(JZ, sl + i, sl + j, dl + k, ps)
        end do
        end do
        end do
    else if (nfbnd(2) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, 0
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
            aj(JX, sl + i, j, dl + k, ps) = aj(JX, sl + i, j, dl + k, ps)*2.0d0
            aj(JY, sl + i, j - 1, dl + k, ps) = -aj(JY, sl + i, -j, dl + k, ps)
            aj(JZ, sl + i, j, dl + k, ps) = aj(JZ, sl + i, j, dl + k, ps)*2.0d0
        end do
        end do
        end do
    end if

!--------------------
    if (bounds(2, 2, sdid(ps) + 1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, nu - 1                !(i.e., do j=0,2)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
            aj(JX, sl + i, yu + du + j, dl + k, ps) = aj(JX, sl + i, yu + du + j, dl + k, ps) &
           &                            + aj(JX, sl + i, yu + su + j, dl + k, ps)
            aj(JY, sl + i, yu + du + j, dl + k, ps) = aj(JY, sl + i, yu + du + j, dl + k, ps) &
           &                            + aj(JY, sl + i, yu + su + j, dl + k, ps)
            aj(JZ, sl + i, yu + du + j, dl + k, ps) = aj(JZ, sl + i, yu + du + j, dl + k, ps) &
           &                            + aj(JZ, sl + i, yu + su + j, dl + k, ps)
        end do
        end do
        end do
    else if (nfbnd(2) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, 0
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
            aj(JX, sl + i, yu + j, dl + k, ps) = aj(JX, sl + i, yu + j, dl + k, ps)*2.0d0
            aj(JY, sl + i, yu + j, dl + k, ps) = -aj(JY, sl + i, yu - j - 1, dl + k, ps)
            aj(JZ, sl + i, yu + j, dl + k, ps) = aj(JZ, sl + i, yu + j, dl + k, ps)*2.0d0
        end do
        end do
        do j = 1, 1
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
            aj(JX, sl + i, yu + j, dl + k, ps) = aj(JX, sl + i, yu - j, dl + k, ps)
            aj(JY, sl + i, yu + j, dl + k, ps) = -aj(JY, sl + i, yu - j - 1, dl + k, ps)
            aj(JZ, sl + i, yu + j, dl + k, ps) = aj(JZ, sl + i, yu - j, dl + k, ps)
        end do
        end do
        end do
    end if

!--------------------
    if (bounds(1, 1, sdid(ps) + 1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, yu + (du + nu - dl) - 1        !(i.e., do j=0,yu+2)
        do i = 0, nl - 1                !(i.e., do i=0,2)
            aj(JX, dl + i, dl + j, dl + k, ps) = aj(JX, dl + i, dl + j, dl + k, ps) &
           &                         + aj(JX, sl + i, dl + j, dl + k, ps)
            aj(JY, dl + i, dl + j, dl + k, ps) = aj(JY, dl + i, dl + j, dl + k, ps) &
           &                         + aj(JY, sl + i, dl + j, dl + k, ps)
            aj(JZ, dl + i, dl + j, dl + k, ps) = aj(JZ, dl + i, dl + j, dl + k, ps) &
           &                         + aj(JZ, sl + i, dl + j, dl + k, ps)
        end do
        end do
        end do
    else if (nfbnd(1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, yu + (du + nu - dl) - 1        !(i.e., do j=0,yu+2)
        do i = 0, 0
            aj(JX, i - 1, dl + j, dl + k, ps) = -aj(JX, -i, dl + j, dl + k, ps)
            aj(JY, i, dl + j, dl + k, ps) = aj(JY, i, dl + j, dl + k, ps)*2.0d0
            aj(JZ, i, dl + j, dl + k, ps) = aj(JZ, i, dl + j, dl + k, ps)*2.0d0
        end do
        end do
        end do
    end if

!--------------------
    if (bounds(2, 1, sdid(ps) + 1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, yu + (du + nu - dl) - 1        !(i.e., do j=0,yu+2)
        do i = 0, nu - 1                !(i.e., do i=0,2)
            aj(JX, xu + du + i, dl + j, dl + k, ps) = aj(JX, xu + du + i, dl + j, dl + k, ps) &
           &                            + aj(JX, xu + su + i, dl + j, dl + k, ps)
            aj(JY, xu + du + i, dl + j, dl + k, ps) = aj(JY, xu + du + i, dl + j, dl + k, ps) &
           &                            + aj(JY, xu + su + i, dl + j, dl + k, ps)
            aj(JZ, xu + du + i, dl + j, dl + k, ps) = aj(JZ, xu + du + i, dl + j, dl + k, ps) &
           &                            + aj(JZ, xu + su + i, dl + j, dl + k, ps)
        end do
        end do
        end do
    else if (nfbnd(1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, yu + (du + nu - dl) - 1        !(i.e., do j=0,yu+2)
        do i = 0, 0
            aj(JX, xu + i, dl + j, dl + k, ps) = -aj(JX, xu - i - 1, dl + j, dl + k, ps)
            aj(JY, xu + i, dl + j, dl + k, ps) = aj(JY, xu + i, dl + j, dl + k, ps)*2.0d0
            aj(JZ, xu + i, dl + j, dl + k, ps) = aj(JZ, xu + i, dl + j, dl + k, ps)*2.0d0
        end do
        do i = 1, 1
            aj(JX, xu + i, dl + j, dl + k, ps) = -aj(JX, xu - i - 1, dl + j, dl + k, ps)
            aj(JY, xu + i, dl + j, dl + k, ps) = aj(JY, xu - i, dl + j, dl + k, ps)
            aj(JZ, xu + i, dl + j, dl + k, ps) = aj(JZ, xu - i, dl + j, dl + k, ps)
        end do
        end do
        end do
    end if

    return
end subroutine

subroutine add_boundary_current2(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!              A D D _ B O U N D A R Y _ C U R R E N T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a charge distribution to      .
!   .      grid points from particle locations by the          .
!   .      "area-sharing" or "area-weighting" scheme.          .
!   ............................................................
!
!
!-------------------- parameter and common block
    use oh_type
    use paramt
    use allcom
    implicit none

    integer(kind=4) :: i, j, k
    integer(kind=4) :: xu, yu, zu
    integer(kind=4) :: sl, su, dl, du, nl, nu
    integer(kind=4) :: is
    integer(kind=4) :: ps

!--------------------
    xu = sdoms(2, 1, sdid(ps) + 1) - sdoms(1, 1, sdid(ps) + 1)
    yu = sdoms(2, 2, sdid(ps) + 1) - sdoms(1, 2, sdid(ps) + 1)
    zu = sdoms(2, 3, sdid(ps) + 1) - sdoms(1, 3, sdid(ps) + 1)

!--------------------
    sl = ctypes(2, 2, 1, CJD)        !=-4
    su = ctypes(2, 1, 1, CJD)        !=+2

    nl = ctypes(3, 2, 1, CJD)        != 3
    nu = ctypes(3, 1, 1, CJD)        != 3

    dl = sl + nl                !=-1
    du = su - nu                !=-1

!--------------------
    if (bounds(1, 3, sdid(ps) + 1) .eq. 1) then
        do k = 0, nl - 1                !(i.e., do k=0,2)
        do j = 0, yu + (su + nu - sl) - 1        !(i.e., do j=0,yu+8)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, sl + i, sl + j, dl + k, ps) = ajdg(is + JX, sl + i, sl + j, dl + k, ps) &
           &                              + ajdg(is + JX, sl + i, sl + j, sl + k, ps)
            ajdg(is + JY, sl + i, sl + j, dl + k, ps) = ajdg(is + JY, sl + i, sl + j, dl + k, ps) &
           &                              + ajdg(is + JY, sl + i, sl + j, sl + k, ps)
            ajdg(is + JZ, sl + i, sl + j, dl + k, ps) = ajdg(is + JZ, sl + i, sl + j, dl + k, ps) &
           &                              + ajdg(is + JZ, sl + i, sl + j, sl + k, ps)
        end do
        end do
        end do
        end do
    else if (nfbnd(3) .eq. 1) then
        do k = 0, 0
        do j = 0, yu + (su + nu - sl) - 1        !(i.e., do j=0,yu+8)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, sl + i, sl + j, k, ps) = ajdg(is + JX, sl + i, sl + j, k, ps)*2.0d0
            ajdg(is + JY, sl + i, sl + j, k, ps) = ajdg(is + JY, sl + i, sl + j, k, ps)*2.0d0
            ajdg(is + JZ, sl + i, sl + j, k - 1, ps) = -ajdg(is + JZ, sl + i, sl + j, -k, ps)
        end do
        end do
        end do
        end do
    end if

    if (bounds(2, 3, sdid(ps) + 1) .eq. 1) then
        do k = 0, nu - 1                !(i.e., do k=0,2)
        do j = 0, yu + (su + nu - sl) - 1        !(i.e., do j=0,yu+8)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, sl + i, sl + j, zu + du + k, ps) = ajdg(is + JX, sl + i, sl + j, zu + du + k, ps) &
           &                                 + ajdg(is + JX, sl + i, sl + j, zu + su + k, ps)
            ajdg(is + JY, sl + i, sl + j, zu + du + k, ps) = ajdg(is + JY, sl + i, sl + j, zu + du + k, ps) &
           &                                 + ajdg(is + JY, sl + i, sl + j, zu + su + k, ps)
            ajdg(is + JZ, sl + i, sl + j, zu + du + k, ps) = ajdg(is + JZ, sl + i, sl + j, zu + du + k, ps) &
           &                                 + ajdg(is + JZ, sl + i, sl + j, zu + su + k, ps)
        end do
        end do
        end do
        end do
    else if (nfbnd(3) .eq. 1) then
        do k = 0, 0
        do j = 0, yu + (su + nu - sl) - 1        !(i.e., do j=0,yu+8)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, sl + i, sl + j, zu + k, ps) = ajdg(is + JX, sl + i, sl + j, zu + k, ps)*2.0d0
            ajdg(is + JY, sl + i, sl + j, zu + k, ps) = ajdg(is + JY, sl + i, sl + j, zu + k, ps)*2.0d0
            ajdg(is + JZ, sl + i, sl + j, zu + k, ps) = -ajdg(is + JZ, sl + i, sl + j, zu - k - 1, ps)
        end do
        end do
        end do
        end do
        do k = 1, 1
        do j = 0, yu + (su + nu - sl) - 1        !(i.e., do j=0,yu+8)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, sl + i, sl + j, zu + k, ps) = ajdg(is + JX, sl + i, sl + j, zu - k, ps)
            ajdg(is + JY, sl + i, sl + j, zu + k, ps) = ajdg(is + JY, sl + i, sl + j, zu - k, ps)
            ajdg(is + JZ, sl + i, sl + j, zu + k, ps) = -ajdg(is + JZ, sl + i, sl + j, zu - k - 1, ps)
        end do
        end do
        end do
        end do
    end if

!--------------------
    if (bounds(1, 2, sdid(ps) + 1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, nl - 1                !(i.e., do j=0,2)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, sl + i, dl + j, dl + k, ps) = ajdg(is + JX, sl + i, dl + j, dl + k, ps) &
           &                              + ajdg(is + JX, sl + i, sl + j, dl + k, ps)
            ajdg(is + JY, sl + i, dl + j, dl + k, ps) = ajdg(is + JY, sl + i, dl + j, dl + k, ps) &
           &                              + ajdg(is + JY, sl + i, sl + j, dl + k, ps)
            ajdg(is + JZ, sl + i, dl + j, dl + k, ps) = ajdg(is + JZ, sl + i, dl + j, dl + k, ps) &
           &                              + ajdg(is + JZ, sl + i, sl + j, dl + k, ps)
        end do
        end do
        end do
        end do
    else if (nfbnd(2) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, 0
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, sl + i, j, dl + k, ps) = ajdg(is + JX, sl + i, j, dl + k, ps)*2.0d0
            ajdg(is + JY, sl + i, j - 1, dl + k, ps) = -ajdg(is + JY, sl + i, -j, dl + k, ps)
            ajdg(is + JZ, sl + i, j, dl + k, ps) = ajdg(is + JZ, sl + i, j, dl + k, ps)*2.0d0
        end do
        end do
        end do
        end do
    end if

!--------------------
    if (bounds(2, 2, sdid(ps) + 1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, nu - 1                !(i.e., do j=0,2)
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, sl + i, yu + du + j, dl + k, ps) = ajdg(is + JX, sl + i, yu + du + j, dl + k, ps) &
           &                                 + ajdg(is + JX, sl + i, yu + su + j, dl + k, ps)
            ajdg(is + JY, sl + i, yu + du + j, dl + k, ps) = ajdg(is + JY, sl + i, yu + du + j, dl + k, ps) &
           &                                 + ajdg(is + JY, sl + i, yu + su + j, dl + k, ps)
            ajdg(is + JZ, sl + i, yu + du + j, dl + k, ps) = ajdg(is + JZ, sl + i, yu + du + j, dl + k, ps) &
           &                                 + ajdg(is + JZ, sl + i, yu + su + j, dl + k, ps)
        end do
        end do
        end do
        end do
    else if (nfbnd(2) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, 0
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, sl + i, yu + j, dl + k, ps) = ajdg(is + JX, sl + i, yu + j, dl + k, ps)*2.0d0
            ajdg(is + JY, sl + i, yu + j, dl + k, ps) = -ajdg(is + JY, sl + i, yu - j - 1, dl + k, ps)
            ajdg(is + JZ, sl + i, yu + j, dl + k, ps) = ajdg(is + JZ, sl + i, yu + j, dl + k, ps)*2.0d0
        end do
        end do
        end do
        do j = 1, 1
        do i = 0, xu + (su + nu - sl) - 1        !(i.e., do i=0,xu+8)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, sl + i, yu + j, dl + k, ps) = ajdg(is + JX, sl + i, yu - j, dl + k, ps)
            ajdg(is + JY, sl + i, yu + j, dl + k, ps) = -ajdg(is + JY, sl + i, yu - j - 1, dl + k, ps)
            ajdg(is + JZ, sl + i, yu + j, dl + k, ps) = ajdg(is + JZ, sl + i, yu - j, dl + k, ps)
        end do
        end do
        end do
        end do
    end if

!--------------------
    if (bounds(1, 1, sdid(ps) + 1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, yu + (du + nu - dl) - 1        !(i.e., do j=0,yu+2)
        do i = 0, nl - 1                !(i.e., do i=0,2)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, dl + i, dl + j, dl + k, ps) = ajdg(is + JX, dl + i, dl + j, dl + k, ps) &
           &                              + ajdg(is + JX, sl + i, dl + j, dl + k, ps)
            ajdg(is + JY, dl + i, dl + j, dl + k, ps) = ajdg(is + JY, dl + i, dl + j, dl + k, ps) &
           &                              + ajdg(is + JY, sl + i, dl + j, dl + k, ps)
            ajdg(is + JZ, dl + i, dl + j, dl + k, ps) = ajdg(is + JZ, dl + i, dl + j, dl + k, ps) &
           &                              + ajdg(is + JZ, sl + i, dl + j, dl + k, ps)
        end do
        end do
        end do
        end do
    else if (nfbnd(1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, yu + (du + nu - dl) - 1        !(i.e., do j=0,yu+2)
        do i = 0, 0
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, i - 1, dl + j, dl + k, ps) = -ajdg(is + JX, -i, dl + j, dl + k, ps)
            ajdg(is + JY, i, dl + j, dl + k, ps) = ajdg(is + JY, i, dl + j, dl + k, ps)*2.0d0
            ajdg(is + JZ, i, dl + j, dl + k, ps) = ajdg(is + JZ, i, dl + j, dl + k, ps)*2.0d0
        end do
        end do
        end do
        end do
    end if

!--------------------
    if (bounds(2, 1, sdid(ps) + 1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, yu + (du + nu - dl) - 1        !(i.e., do j=0,yu+2)
        do i = 0, nu - 1                !(i.e., do i=0,2)
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, xu + du + i, dl + j, dl + k, ps) = ajdg(is + JX, xu + du + i, dl + j, dl + k, ps) &
           &                                 + ajdg(is + JX, xu + su + i, dl + j, dl + k, ps)
            ajdg(is + JY, xu + du + i, dl + j, dl + k, ps) = ajdg(is + JY, xu + du + i, dl + j, dl + k, ps) &
           &                                 + ajdg(is + JY, xu + su + i, dl + j, dl + k, ps)
            ajdg(is + JZ, xu + du + i, dl + j, dl + k, ps) = ajdg(is + JZ, xu + du + i, dl + j, dl + k, ps) &
           &                                 + ajdg(is + JZ, xu + su + i, dl + j, dl + k, ps)
        end do
        end do
        end do
        end do
    else if (nfbnd(1) .eq. 1) then
        do k = 0, zu + (du + nu - dl) - 1        !(i.e., do k=0,zu+2)
        do j = 0, yu + (du + nu - dl) - 1        !(i.e., do j=0,yu+2)
        do i = 0, 0
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, xu + i, dl + j, dl + k, ps) = -ajdg(is + JX, xu - i - 1, dl + j, dl + k, ps)
            ajdg(is + JY, xu + i, dl + j, dl + k, ps) = ajdg(is + JY, xu + i, dl + j, dl + k, ps)*2.0d0
            ajdg(is + JZ, xu + i, dl + j, dl + k, ps) = ajdg(is + JZ, xu + i, dl + j, dl + k, ps)*2.0d0
        end do
        end do
        do i = 1, 1
        do is = 0, nspec*3 - 1, 3
            ajdg(is + JX, xu + i, dl + j, dl + k, ps) = -ajdg(is + JX, xu - i - 1, dl + j, dl + k, ps)
            ajdg(is + JY, xu + i, dl + j, dl + k, ps) = ajdg(is + JY, xu - i, dl + j, dl + k, ps)
            ajdg(is + JZ, xu + i, dl + j, dl + k, ps) = ajdg(is + JZ, xu - i, dl + j, dl + k, ps)
        end do
        end do
        end do
        end do
    end if

    return
end subroutine

subroutine add_source_current(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!              A D D _ S O U R C E _ C U R R E N T
!   ____________________________________________________________
!
!   ............................................................
!   ............................................................
!
!
!-------------------- parameter and common block
    use oh_type
    use paramt
    use allcom
    implicit none
!
    integer(kind=4) :: i, j, k, ijs, jcomp
    integer(kind=4) :: xl, yl, zl, xu, yu, zu
    integer(kind=4) :: ps

!--------------------
    xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
    yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
    zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)

!--------------------
    do ijs = 1, njs
        i = rjs(1, ijs); j = rjs(2, ijs); k = rjs(3, ijs)
        if (i .ge. xl .and. i .lt. xu .and. &
       &   j .ge. yl .and. j .lt. yu .and. &
       &   k .ge. zl .and. k .lt. zu) then
            aj(JX:JZ, i - xl, j - yl, k - zl, ps) &
           &  = aj(JX:JZ, i - xl, j - yl, k - zl, ps) &
           &  + ajs(JX:JZ, ijs) &
           &  *sin(wjs(ijs)*(2.0d0*istep + 1.0d0) + th0js(ijs))
        end if
    end do

    return
end subroutine
