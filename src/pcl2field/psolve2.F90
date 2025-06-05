#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"

subroutine psolve2(ps, ustep)
!
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
    implicit none

    integer(kind=8) :: m, mm, ns, ne, nee
    integer(kind=8) :: nprs, npre
    integer(kind=8) :: isfluxt(inpc, ispec)
    integer(kind=8) :: influxsum
    integer(kind=4) :: i, j, k
    integer(kind=4) :: i1, j1, k1
    integer(kind=4) :: is, ustep, ibdy, ipc, itch
    integer(kind=4) :: nb(3)
    integer(kind=4) :: ps, rid, nud
    integer(kind=4) :: nyield
    integer(kind=4) :: icon, iran
    integer(kind=4) :: addr0, addr1, addr2
!  integer(kind=4) :: oh3_map_region_to_node
    real(kind=8) :: xlocal, ylocal, zlocal
    real(kind=8) :: x, y, z
    real(kind=8) :: x1, y1, z1, z2, xy1, xz1, yz1, xz2, yz2
    real(kind=8) :: tslx, tsly, tslz, gustep, rustep, qmp
    real(kind=8) :: v1, v2, v3, v4, v5, v6, v7, v8
    real(kind=8) :: eex, eey, eez, bbx, bby, bbz, boris, vxt, vyt, vzt, vxyz
    real(kind=8) :: disp1, disp2, disp3
    real(kind=8) :: radsq(inpc), rbwlsq, rdomsq
    real(kind=8) :: vxx, vyy, vzz, vxy, vyx, vyz, vzy, vzx, vxz
    real(kind=8) :: vxymyx, vyzmzy, vzxmxz
    real(kind=8) :: discrim, v2norm, v2normi, paramt1, paramt2, llmt, ulmt
    real(kind=8) :: perg, pergfrac, pemaxinvh(ispec), costhi
    real(kind=8) :: xpast0, ypast0, zpast0, xpast1, ypast1, zpast1
    real(kind=8) :: yield, weightr, delta

    real(kind=8) :: displx, disply, displz, displ3i

    real(kind=8) :: xmove, ymove, zmove, xsepa, ysepa, zsepa, termA, termB
    real(kind=8) :: tfrac, tfrac1, tfrac2
    real(kind=8) :: xpast1a, ypast1a, xpast1b, ypast1b
    real(kind=8) :: xpast2a, ypast2a, xpast2b, ypast2b
    real(kind=8) :: cosbetav, sinbetav

    real(kind=8) :: nsigmadt, colfp, uppb
    real(kind=8) :: ewx, ewy, tew
    integer(kind=4) :: chkflg

    real(kind=8) :: wrk_e(EX:EZ, 8)

    type(t_CollisionRecord) :: record
    type(t_BoundaryList) :: boundaries

    logical :: is_boundary_condition_applied

    call ohhinfo_update(sdoms(:, :, sdid(ps) + 1))

    boundaries = create_simple_collision_boundaries(sdoms(1:2, 1:3, sdid(ps) + 1))

    pemaxinvh(1:nspec) = 0.5d0/pemax(1:nspec)

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

!============================== species loop
    nee = pbase(ps)
    npre = 0
    do is = 1, isse
        nprs = npre + 1
        npre = npre + npr(is)
    end do
    rustep = 1.0d0/ustep
    tslx = 2.0d0*slx
    tsly = 2.0d0*sly
    tslz = 2.0d0*slz
    ISL1: do is = 1, nspec
        ns = nee + 1
        ne = nee + totalp(is, ps)
        nee = ne
!       -------------
        weightr = q(is)/q(isse)
        delta = weightr*deltaemax(is)
        if (ewmodel .eq. 1 .or. ewmodel .eq. 2) then
            tew = t - dt*nretard
            ewx = Ew(1, is, 1)*cos(omegaw(1)*tew)
            ewy = Ew(2, is, 1)*sin(omegaw(1)*tew)
        else
            ewx = 0.0d0
            ewy = 0.0d0
        end if
!       -------------
        if (mfpath(is) .gt. 1.0d0) then
            nsigmadt = ustep/mfpath(is)
        else if (mfpath(is) .lt. -1.0d0) then
            nsigmadt = ustep/mfpath(is)
        else
            nsigmadt = 0.0d0
        end if
        uppb = 1.0d0 - exp(-vxyzmax(is)*abs(nsigmadt))
!       ------------- inner loop
        do m = ns, ne
            if (pbuf(m)%preside .lt. 0) then
                nphgram(pbuf(m)%nid + 1, is, ps) = nphgram(pbuf(m)%nid + 1, is, ps) - 1
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

            if (pbuf(m)%nid .eq. -1) cycle

            xlocal = pbuf(m)%x - dxl
            ylocal = pbuf(m)%y - dyl
            zlocal = pbuf(m)%z - dzl

            i = int(xlocal)
            j = int(ylocal)
            k = int(zlocal)
!
            i1 = i + 1
            j1 = j + 1
            k1 = k + 1
!
            x1 = xlocal - i
            y1 = ylocal - j
            z1 = zlocal - k
!
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

            wrk_e(EX:EZ, 1) = wrk(EX:EZ, i, j, k1)
            wrk_e(EX:EZ, 2) = wrk(EX:EZ, i1, j, k1)
            wrk_e(EX:EZ, 3) = wrk(EX:EZ, i1, j1, k1)
            wrk_e(EX:EZ, 4) = wrk(EX:EZ, i, j1, k1)
            wrk_e(EX:EZ, 5) = wrk(EX:EZ, i, j, k)
            wrk_e(EX:EZ, 6) = wrk(EX:EZ, i1, j, k)
            wrk_e(EX:EZ, 7) = wrk(EX:EZ, i1, j1, k)
            wrk_e(EX:EZ, 8) = wrk(EX:EZ, i, j1, k)
!
!         ----------- eb-fields interpolation
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
            bbx = wrk(BX, i, j, k1)*v1 + wrk(BX, i1, j, k1)*v2 &
           &    + wrk(BX, i1, j1, k1)*v3 + wrk(BX, i, j1, k1)*v4 &
           &    + wrk(BX, i, j, k)*v5 + wrk(BX, i1, j, k)*v6 &
           &    + wrk(BX, i1, j1, k)*v7 + wrk(BX, i, j1, k)*v8
            bby = wrk(BY, i, j, k1)*v1 + wrk(BY, i1, j, k1)*v2 &
           &    + wrk(BY, i1, j1, k1)*v3 + wrk(BY, i, j1, k1)*v4 &
           &    + wrk(BY, i, j, k)*v5 + wrk(BY, i1, j, k)*v6 &
           &    + wrk(BY, i1, j1, k)*v7 + wrk(BY, i, j1, k)*v8
            bbz = wrk(BZ, i, j, k1)*v1 + wrk(BZ, i1, j, k1)*v2 &
           &    + wrk(BZ, i1, j1, k1)*v3 + wrk(BZ, i, j1, k1)*v4 &
           &    + wrk(BZ, i, j, k)*v5 + wrk(BZ, i1, j, k)*v6 &
           &    + wrk(BZ, i1, j1, k)*v7 + wrk(BZ, i, j1, k)*v8
!
!         -----------
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
!
!         ----------- charge-to-mass ratio is taken into consideration
            eex = (eex + ewx)*qmp
            eey = (eey + ewy)*qmp
            eez = eez*qmp
            bbx = bbx*qmp
            bby = bby*qmp
            bbz = bbz*qmp
!
!         ----------- update particle velocities (Buneman-Boris method)
            boris = 2.0d0/(1.0d0 + bbx*bbx + bby*bby + bbz*bbz)
!
            pbuf(m)%vx = pbuf(m)%vx + eex
            pbuf(m)%vy = pbuf(m)%vy + eey
            pbuf(m)%vz = pbuf(m)%vz + eez
!
            vxt = pbuf(m)%vx + pbuf(m)%vy*bbz - pbuf(m)%vz*bby
            vyt = pbuf(m)%vy + pbuf(m)%vz*bbx - pbuf(m)%vx*bbz
            vzt = pbuf(m)%vz + pbuf(m)%vx*bby - pbuf(m)%vy*bbx
!
            pbuf(m)%vx = pbuf(m)%vx + boris*(vyt*bbz - vzt*bby)
            pbuf(m)%vy = pbuf(m)%vy + boris*(vzt*bbx - vxt*bbz)
            pbuf(m)%vz = pbuf(m)%vz + boris*(vxt*bby - vyt*bbx)
!
            pbuf(m)%vx = pbuf(m)%vx + eex
            pbuf(m)%vy = pbuf(m)%vy + eey
            pbuf(m)%vz = pbuf(m)%vz + eez

!         ----------- velocity modifications associated with collisions
            colfp = colf(1, i, j, k1, ps)*v1 + colf(1, i1, j, k1, ps)*v2 &
           &      + colf(1, i1, j1, k1, ps)*v3 + colf(1, i, j1, k1, ps)*v4 &
           &      + colf(1, i, j, k, ps)*v5 + colf(1, i1, j, k, ps)*v6 &
           &      + colf(1, i1, j1, k, ps)*v7 + colf(1, i, j1, k, ps)*v8
            if ((nsigmadt .gt. 0.0d0) .and. &
           &   (pbuf(m)%x .le. xlcol(1) .or. pbuf(m)%x .ge. xucol(1) .or. &
           &    pbuf(m)%y .le. ylcol(1) .or. pbuf(m)%y .ge. yucol(1) .or. &
           &    pbuf(m)%z .le. zlcol(1) .or. pbuf(m)%z .ge. zucol(1))) then
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
           &   (pbuf(m)%x .le. xlcol(1) .or. pbuf(m)%x .ge. xucol(1) .or. &
           &    pbuf(m)%y .le. ylcol(1) .or. pbuf(m)%y .ge. yucol(1) .or. &
           &    pbuf(m)%z .le. zlcol(1) .or. pbuf(m)%z .ge. zucol(1))) then
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

!         ----------- update particle positions
            xmove = pbuf(m)%vx*gustep
            ymove = pbuf(m)%vy*gustep
            zmove = pbuf(m)%vz*gustep
            pbuf(m)%x = pbuf(m)%x + xmove
            pbuf(m)%y = pbuf(m)%y + ymove
            pbuf(m)%z = pbuf(m)%z + zmove

!         ----------- internal boundary treatment 0
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
            end block
            if (record%is_collided) then
                pbuf(m)%x = record%position(1)
                pbuf(m)%y = record%position(2)
                pbuf(m)%z = record%position(3)
                pbuf(m)%preside = OH_PCL_TO_BE_ACCUMULATED
            end if
! #include "defsurf.fnc"

           is_boundary_condition_applied = .false.

!         ----------- external boundary treatment
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
            if (nud .ne. 14) then
                rid = nborps(nud, is, ps)
                pbuf(m)%nid = rid
                nphgram(sdid(ps) + 1, is, ps) = nphgram(sdid(ps) + 1, is, ps) - 1
!
                if (pbuf(m)%x .lt. 0.0d0) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(1, is) .eq. 0) then
                        pbuf(m)%x = pbuf(m)%x + slx
                    else if (npbnd(1, is) .eq. 1) then
                        pbuf(m)%x = -pbuf(m)%x
                        pbuf(m)%vx = -pbuf(m)%vx
                    end if
                else if (pbuf(m)%x .ge. slx) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(1, is) .eq. 0) then
                        pbuf(m)%x = pbuf(m)%x - slx
                    else if (npbnd(1, is) .eq. 1) then
                        pbuf(m)%x = tslx - pbuf(m)%x
                        pbuf(m)%vx = -pbuf(m)%vx
                    end if
                end if
!
                if (pbuf(m)%y .lt. 0.0d0) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(2, is) .eq. 0) then
                        pbuf(m)%y = pbuf(m)%y + sly
                    else if (npbnd(2, is) .eq. 1) then
                        pbuf(m)%y = -pbuf(m)%y
                        pbuf(m)%vy = -pbuf(m)%vy
                    end if
                else if (pbuf(m)%y .ge. sly) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(2, is) .eq. 0) then
                        pbuf(m)%y = pbuf(m)%y - sly
                    else if (npbnd(2, is) .eq. 1) then
                        pbuf(m)%y = tsly - pbuf(m)%y
                        pbuf(m)%vy = -pbuf(m)%vy
                    end if
                end if
!
                if (pbuf(m)%z .lt. 0.0d0) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(3, is) .eq. 0) then
                        pbuf(m)%z = pbuf(m)%z + slz
                    else if (npbnd(3, is) .eq. 1) then
                        pbuf(m)%z = -pbuf(m)%z
                        pbuf(m)%vz = -pbuf(m)%vz
                    end if
                else if (pbuf(m)%z .ge. slz) then
                    is_boundary_condition_applied = .true.
                    if (npbnd(3, is) .eq. 0) then
                        pbuf(m)%z = pbuf(m)%z - slz
                    else if (npbnd(3, is) .eq. 1) then
                        pbuf(m)%z = tslz - pbuf(m)%z
                        pbuf(m)%vz = -pbuf(m)%vz
                    end if
                end if
!
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

!         ----------- internal boundary treatment 1
!          if(line_mode.ne.1.and.pbuf(m)%nid.ne.-1) then
            if (line_mode .ne. 1) then
                do ipc = 1, npc
                    if (pbuf(m)%nid .eq. -1) cycle
                    if ((geotype(ipc) .eq. 0 .or. geotype(ipc) .eq. 1) .and. &
                   &   (pbuf(m)%x .ge. xlpc(ipc) .and. pbuf(m)%x .le. xupc(ipc) .and. &
                   &    pbuf(m)%y .ge. ylpc(ipc) .and. pbuf(m)%y .le. yupc(ipc) .and. &
                   &    pbuf(m)%z .ge. zlpc(ipc) .and. pbuf(m)%z .le. zupc(ipc))) then
                        nphgram(pbuf(m)%nid + 1, is, ps) = nphgram(pbuf(m)%nid + 1, is, ps) - 1
                        pbuf(m)%nid = -1
                        gcount(1)%influx(ipc, is) = gcount(1)%influx(ipc, is) + 1
                        gcount(1)%infhist(pbuf(m)%pid, ipc, is) = &
                       &  gcount(1)%infhist(pbuf(m)%pid, ipc, is) + 1
!
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
        end do

    end do ISL1
!
    ISL2: do is = 1, nspec
        do ipc = 1, npc
            if (ncond(ipc) .gt. 0) then
                gcount(1)%chgacm(:, ipc) = gcount(1)%chgacm(:, ipc) &
               &  + gcount(1)%influx(ipc, is)*q(is)*sqdscaled(ipc)
            end if
        end do
    end do ISL2

    call boundaries%destroy

end subroutine
