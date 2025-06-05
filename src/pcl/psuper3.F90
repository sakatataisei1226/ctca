#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
subroutine psuper3(ustep)
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
    use finbound
    use m_raycast
    use particle_collision
#define MCW local_comm
#define MSS MPI_STATUS_SIZE
    implicit none

    integer, intent(in) :: ustep

    integer :: icon

    integer :: xl, xu, yl, yu, zl, zu
    double precision :: align, radius

    integer(kind=8) :: nprs, npre
    integer :: nns, nne
    integer :: nemit(inepl), omni

    integer :: m
    integer :: ie, ipcplane, nemitw
    integer :: iepl, ipc, is

    xl = sdoms(1, 1, sdid(1) + 1); xu = sdoms(2, 1, sdid(1) + 1)
    yl = sdoms(1, 2, sdid(1) + 1); yu = sdoms(2, 2, sdid(1) + 1)
    zl = sdoms(1, 3, sdid(1) + 1); zu = sdoms(2, 3, sdid(1) + 1)

    npre = 0
    nne = 0

    do is = 1, nspec
        nprs = npre + 1
        npre = npre + npr(is)
        nns = nne + 1
        nne = nne + nepl(is)

        if (pe_ray_cast .and. nflag_emit(is) == 2) then
            call emit_with_raycasting
            cycle
        end if

        if (nflag_emit(is) == 1 .or. nflag_emit(is) == 2) then
            if (injct(is) == 0) cycle

            ie = mod(istep, injct(is))
            if (ie /= 0 .or. inpf(is) == 0) cycle

            ! Emit from prescribed surfaces.
            do iepl = nns, nne
                nemit(iepl) = (int((istep + lstep)*fluxf(iepl), 8) &
                               - int((istep + lstep - injct(is))*fluxf(iepl), 8))*ustep
                ipcplane = ipcpl(iepl)

                nemitw = nemit(iepl)
                omni = omniemit(iepl)

                ! Generate random numbers.
                if (nemitw > 0) then
                    call RANU0(dranu, nemitw*6, icon)
                    if (icon /= 0) print *, "[psuper] icon=", icon
                end if

                ! pos. & vel. assignment
                if (geotype(ipcplane) == 0 .or. geotype(ipcplane) == 1) then
                    call emit_from_rectangle
                else if (geotype(ipcplane) == 2) then
                    call emit_from_cylinder
                end if

                ! Calculate accumulated charge.
                if (ncond(ipc) > 0) then
                    gcount(1)%chgacm(:, ipcplane) = gcount(1)%chgacm(:, ipcplane) &
                                                    - q(is)*nemitw*sqdscaled(ipcplane)
                end if
                gcount(1)%outflux(ipcplane, is) = gcount(1)%outflux(ipcplane, is) &
                                                  + nemitw
            end do

        end if
    end do

contains

    subroutine emit_with_raycasting
        double precision :: rnrays

        integer :: nrays
        integer :: npgen

        type(t_ParallelCamera) :: camera
        type(t_BoundaryList) :: boundaries

        block
            integer :: doms(2, 3)
            doms = reshape([[0, nx], [0, ny], [0, nz]], [2, 3])
            boundaries = create_simple_collision_boundaries(doms, cover_all=.true.)
        end block

        camera = new_ParallelCamera_optimized(vdthz(1)/180d0*pi, (vdthxy(1) - 180d0)/180d0*pi, &
                                              reshape([[0, nx], [0, ny], [0, nz]], [2, 3]))
        ! reshape([[0, nx], [0, ny], [0, int(zssurf)]], [2, 3]))
        ! camera = new_ParallelCamera(vdthz(1)/180d0*pi, (vdthxy(1) - 180d0)/180d0*pi, nx, ny, nz)

        ! nrays is multiplyed by 2 times for ustep = 2.
        rnrays = abs(curf(is)*renj*camera%S/q(is)*ustep)
        nrays = int(rnrays)
        if (nrays == 0) then
            return
        end if

        call RANU0(dranu, 1, icon)
        if (dranu(1) < (rnrays - nrays)) then
            nrays = nrays + 1
        end if

        call generate_particles_with_raycast(camera, boundaries, &
                                             is, nrays, &
                                             npgen)

        call boundaries%destroy
    end subroutine

    subroutine emit_from_rectangle
        integer :: rid

        do m = 1, nemit(iepl)
            if (nemd(iepl) == 1) then
                call emit_from_rectangular_surf &
                    (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, +1, &
                     peject(iepl)%grd, peject(iepl)%yl, peject(iepl)%yu, &
                     peject(iepl)%zl, peject(iepl)%zu, m*6 - 5)
            else if (nemd(iepl) == -1) then
                call emit_from_rectangular_surf &
                    (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, -1, &
                     peject(iepl)%grd, peject(iepl)%yl, peject(iepl)%yu, &
                     peject(iepl)%zl, peject(iepl)%zu, m*6 - 5)
            else if (nemd(iepl) == 2) then
                call emit_from_rectangular_surf &
                    (pinj(1)%vy, pinj(1)%vz, pinj(1)%vx, pinj(1)%y, pinj(1)%z, pinj(1)%x, +1, &
                     peject(iepl)%grd, peject(iepl)%zl, peject(iepl)%zu, &
                     peject(iepl)%xl, peject(iepl)%xu, m*6 - 5)
            else if (nemd(iepl) == -2) then
                call emit_from_rectangular_surf &
                    (pinj(1)%vy, pinj(1)%vz, pinj(1)%vx, pinj(1)%y, pinj(1)%z, pinj(1)%x, -1, &
                     peject(iepl)%grd, peject(iepl)%zl, peject(iepl)%zu, &
                     peject(iepl)%xl, peject(iepl)%xu, m*6 - 5)
            else if (nemd(iepl) == 3) then
                call emit_from_rectangular_surf &
                    (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, +1, &
                     peject(iepl)%grd, peject(iepl)%xl, peject(iepl)%xu, &
                     peject(iepl)%yl, peject(iepl)%yu, m*6 - 5)
            else if (nemd(iepl) == -3) then
                call emit_from_rectangular_surf &
                    (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, -1, &
                     peject(iepl)%grd, peject(iepl)%xl, peject(iepl)%xu, &
                     peject(iepl)%yl, peject(iepl)%yu, m*6 - 5)
            end if

            if (pinj(1)%x < xl .or. pinj(1)%x >= xu .or. &
                pinj(1)%y < yl .or. pinj(1)%y >= yu .or. &
                pinj(1)%z < zl .or. pinj(1)%z >= zu) then
                rid = oh3_map_particle_to_subdomain(pinj(1)%x, pinj(1)%y, pinj(1)%z)
                pinj(1)%nid = rid
            else
                pinj(1)%nid = sdid(1)
            end if
            pinj(1)%spec = is
            pinj(1)%pid = 0

            pinj(2) = pinj(1)

            if (nemd(iepl) == 1) then
                pinj(2)%x = pinj(2)%x + 1e-5
            else if (nemd(iepl) == -1) then
                pinj(2)%x = pinj(2)%x - 1e-5
            else if (nemd(iepl) == 2) then
                pinj(2)%y = pinj(2)%y + 1e-5
            else if (nemd(iepl) == -2) then
                pinj(2)%y = pinj(2)%y - 1e-5
            else if (nemd(iepl) == 3) then
                pinj(2)%z = pinj(2)%z + 1e-5
            else if (nemd(iepl) == -3) then
                pinj(2)%z = pinj(2)%z - 1e-5
            end if

            pinj(1)%preside = -2
            call oh2_inject_particle(pinj(1))
            pinj(2)%preside = +1
            call oh2_inject_particle(pinj(2))
        end do
    end subroutine

    subroutine emit_from_rectangular_surf(vnorm, vtang1, vtang2, &
                                          rnorm, rtang1, rtang2, &
                                          ud, surfp, &
                                          hl1, hu1, hl2, hu2, &
                                          ir)
        double precision, intent(out) :: vnorm, vtang1, vtang2
        double precision, intent(out) :: rnorm, rtang1, rtang2
        integer, intent(in) :: ud
        double precision, intent(in) :: surfp
        double precision, intent(in) :: hl1, hu1, hl2, hu2
        integer, intent(in) :: ir

        block ! Pick up from velocity distribution
            integer :: addr0, addr1, addr2

            addr0 = nprs + npr(is)*dranu(ir)
            addr1 = nprs + npr(is)*dranu(ir + 1)
            addr2 = nprs + npr(is)*dranu(ir + 2)

            vnorm = ud*vnormf(addr0)
            vtang1 = vtangf(addr1)
            vtang2 = vtangf(addr2)
        end block

        rnorm = surfp
        rtang1 = hl1 + (hu1 - hl1)*dranu(ir + 4)
        rtang2 = hl2 + (hu2 - hl2)*dranu(ir + 5)
    end subroutine emit_from_rectangular_surf

    subroutine emit_from_cylinder
        integer :: rid

        align = cylinder(ipcplane)%align
        radius = cylinder(ipcplane)%radius

        if (abs(nemd(iepl)) == align) then
            do m = 1, nemit(iepl)
                if (nemd(iepl) == 1) then
                    call emit_from_circular_surf &
                        (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, +1, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%grd, m*6 - 5)
                else if (nemd(iepl) == -1) then
                    call emit_from_circular_surf &
                        (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, -1, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%grd, m*6 - 5)
                else if (nemd(iepl) == 2) then
                    call emit_from_circular_surf &
                        (pinj(1)%vy, pinj(1)%vz, pinj(1)%vx, pinj(1)%y, pinj(1)%z, pinj(1)%x, +1, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%grd, m*6 - 5)
                else if (nemd(iepl) == -2) then
                    call emit_from_circular_surf &
                        (pinj(1)%vy, pinj(1)%vz, pinj(1)%vx, pinj(1)%y, pinj(1)%z, pinj(1)%x, -1, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%grd, m*6 - 5)
                else if (nemd(iepl) == 3) then
                    call emit_from_circular_surf &
                        (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, +1, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%grd, m*6 - 5)
                else if (nemd(iepl) == -3) then
                    call emit_from_circular_surf &
                        (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, -1, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%grd, m*6 - 5)
                end if

                if (pinj(1)%x < xl .or. pinj(1)%x >= xu .or. &
                    pinj(1)%y < yl .or. pinj(1)%y >= yu .or. &
                    pinj(1)%z < zl .or. pinj(1)%z >= zu) then
                    rid = oh3_map_particle_to_subdomain(pinj(1)%x, pinj(1)%y, pinj(1)%z)
                    pinj(1)%nid = rid
                else
                    pinj(1)%nid = sdid(1)
                end if
                pinj(1)%spec = is
                pinj(1)%pid = 0

                pinj(2) = pinj(1)

                pinj(1)%preside = +1
                call oh2_inject_particle(pinj(1))
                if (ncond(ipc) <= 0) then
                    pinj(2)%preside = -2
                    call oh2_inject_particle(pinj(2))
                end if
            end do
        else
            do m = 1, nemit(iepl)
                if (align == 1 .and. omni == 0 .and. radius >= 1.0d0) then
                    call emit_from_cylindrical_surf1 &
                        (pinj(1)%vz, pinj(1)%vy, pinj(1)%vx, pinj(1)%z, pinj(1)%y, pinj(1)%x, psizy, &
                         cylinder(ipcplane)%axis(2), cylinder(ipcplane)%axis(1), &
                         peject(iepl)%xl, peject(iepl)%xu, m*6 - 5)
                else if (align == 1 .and. omni == 0 .and. radius < 1.0d0) then
                    call emit_from_cylindrical_surf2 &
                        (pinj(1)%vz, pinj(1)%vy, pinj(1)%vx, pinj(1)%z, pinj(1)%y, pinj(1)%x, psizy, &
                         cylinder(ipcplane)%axis(2), cylinder(ipcplane)%axis(1), &
                         peject(iepl)%xl, peject(iepl)%xu, m*5 - 4)
                else if (align == 1 .and. omni == 1) then
                    call emit_from_cylindrical_surf3 &
                        (pinj(1)%vz, pinj(1)%vy, pinj(1)%vx, pinj(1)%z, pinj(1)%y, pinj(1)%x, psizy, &
                         cylinder(ipcplane)%axis(2), cylinder(ipcplane)%axis(1), &
                         peject(iepl)%xl, peject(iepl)%xu, m*6 - 5)
                else if (align == 2 .and. omni == 0 .and. radius >= 1.0d0) then
                    call emit_from_cylindrical_surf1 &
                        (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, psizx, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%yl, peject(iepl)%yu, m*6 - 5)
                else if (align == 2 .and. omni == 0 .and. radius < 1.0d0) then
                    call emit_from_cylindrical_surf2 &
                        (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, psizx, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%yl, peject(iepl)%yu, m*5 - 4)
                else if (align == 2 .and. omni == 1) then
                    call emit_from_cylindrical_surf3 &
                        (pinj(1)%vz, pinj(1)%vx, pinj(1)%vy, pinj(1)%z, pinj(1)%x, pinj(1)%y, psizx, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%yl, peject(iepl)%yu, m*6 - 5)
                else if (align == 3 .and. omni == 0 .and. radius >= 1.0d0) then
                    call emit_from_cylindrical_surf1 &
                        (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, psixy, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%zl, peject(iepl)%zu, m*6 - 5)
                else if (align == 3 .and. omni == 0 .and. radius < 1.0d0) then
                    call emit_from_cylindrical_surf2 &
                        (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, psixy, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%zl, peject(iepl)%zu, m*5 - 4)
                else if (align == 3 .and. omni == 1) then
                    call emit_from_cylindrical_surf3 &
                        (pinj(1)%vx, pinj(1)%vy, pinj(1)%vz, pinj(1)%x, pinj(1)%y, pinj(1)%z, psixy, &
                         cylinder(ipcplane)%axis(1), cylinder(ipcplane)%axis(2), &
                         peject(iepl)%zl, peject(iepl)%zu, m*6 - 5)
                end if

                if (pinj(1)%x < xl .or. pinj(1)%x >= xu .or. &
                    pinj(1)%y < yl .or. pinj(1)%y >= yu .or. &
                    pinj(1)%z < zl .or. pinj(1)%z >= zu) then
                    rid = oh3_map_particle_to_subdomain(pinj(1)%x, pinj(1)%y, pinj(1)%z)
                    pinj(1)%nid = rid
                else
                    pinj(1)%nid = sdid(1)
                end if
                pinj(1)%spec = is
                pinj(1)%pid = 0

                pinj(2) = pinj(1)

                pinj(1)%preside = +1
                call oh2_inject_particle(pinj(1))
                if (ncond(ipc) <= 0) then
                    pinj(2)%preside = -2
                    call oh2_inject_particle(pinj(2))
                end if
            end do
        end if
    end subroutine

    subroutine emit_from_circular_surf(vnorm, vtang1, vtang2, &
                                       rnorm, rtang1, rtang2, &
                                       ud, &
                                       axis1, axis2, &
                                       surfp, &
                                       ir)
        double precision, intent(out)   :: vnorm, vtang1, vtang2
        double precision, intent(out)   :: rnorm, rtang1, rtang2
        integer, intent(in) :: ud
        double precision, intent(in)    :: axis1, axis2
        double precision, intent(in)    :: surfp
        integer, intent(in) :: ir

        block ! Pick up from velocity distribution
            integer :: addr0, addr1, addr2
            addr0 = nprs + npr(is)*dranu(ir)
            addr1 = nprs + npr(is)*dranu(ir + 1)
            addr2 = nprs + npr(is)*dranu(ir + 2)

            vnorm = ud*vnormf(addr0)
            vtang1 = vtangf(addr1)
            vtang2 = vtangf(addr2)
        end block

        rnorm = surfp

        block
            integer :: j3
            double precision :: radial
            double precision :: angular

            j3 = nprs + npr(is)*dranu(ir + 3)
            radial = radius*runiform(j3)
            angular = pi2*dranu(ir + 4)
            rtang1 = axis1 + radial*cos(angular)
            rtang2 = axis2 + radial*sin(angular)
        end block
    end subroutine emit_from_circular_surf

    subroutine emit_from_cylindrical_surf1(vnorm, vtang1, vtang2, &
                                           rnorm, rtang1, rtang2, &
                                           psi, &
                                           axis1, axis2, &
                                           hl, hu, &
                                           ir)
        double precision, intent(out)   :: vnorm, vtang1, vtang2
        double precision, intent(out)   :: rnorm, rtang1, rtang2
        double precision, intent(in)    :: psi
        double precision, intent(in)    :: axis1, axis2
        double precision, intent(in)    :: hl, hu
        integer, intent(in) :: ir

        block ! Pick up from velocity distribution
            integer :: addr0, addr1, addr2

            addr0 = nprs + npr(is)*dranu(ir)
            addr1 = nprs + npr(is)*dranu(ir + 1)
            addr2 = nprs + npr(is)*dranu(ir + 2)

            vnorm = vnormc(addr0)*dcos(psicos(addr2) + psi) &
                    - vtangc(addr1)*dsin(psicos(addr2) + psi)
            vtang1 = vnormc(addr0)*dsin(psicos(addr2) + psi) &
                     + vtangc(addr1)*dcos(psicos(addr2) + psi)

            rnorm = axis1 &
                    + radius*dcos(psicos(addr2) + psi)
            rtang1 = axis2 &
                     + radius*dsin(psicos(addr2) + psi)
        end block

        block
            integer :: j3
            j3 = nprs + npr(is)*dranu(ir + 3)
            vtang2 = vtangf(j3)
        end block

        rtang2 = hl + (hu - hl)*dranu(ir + 5)
    end subroutine emit_from_cylindrical_surf1

    subroutine emit_from_cylindrical_surf2(vnorm, vtang1, vtang2, &
                                           rnorm, rtang1, rtang2, &
                                           psi, &
                                           axis1, axis2, &
                                           hl, hu, &
                                           ir)
        double precision, intent(out)   :: vnorm, vtang1, vtang2
        double precision, intent(out)   :: rnorm, rtang1, rtang2
        double precision, intent(in)    :: psi
        double precision, intent(in)    :: axis1, axis2
        double precision, intent(in)    :: hl, hu
        integer, intent(in) :: ir

        block ! Pick up from velocity distribution
            integer :: addr0, addr1, addr2

            addr0 = nprs + npr(is)*dranu(ir)
            addr1 = nprs + npr(is)*dranu(ir + 1)
            addr2 = nprs + npr(is)*dranu(ir + 2)

            vnorm = vnormc(addr0)*dcos(psi) &
                    - vtangc(addr1)*dsin(psi)
            vtang1 = vnormc(addr0)*dsin(psi) &
                     + vtangc(addr1)*dcos(psi)
            vtang2 = vtangf(addr2)
        end block

        block
            double precision :: v2normi

            v2normi = radius/dsqrt(vnorm*vnorm + vtang1*vtang1)
            rnorm = axis1 + vnorm*v2normi
            rtang1 = axis2 + vtang1*v2normi
        end block

        rtang2 = hl + (hu - hl)*dranu(ir + 4)
    end subroutine emit_from_cylindrical_surf2

    subroutine emit_from_cylindrical_surf3(vnorm, vtang1, vtang2, &
                                           rnorm, rtang1, rtang2, &
                                           psi, &
                                           axis1, axis2, &
                                           hl, hu, &
                                           ir)
        double precision, intent(out)   :: vnorm, vtang1, vtang2
        double precision, intent(out)   :: rnorm, rtang1, rtang2
        double precision, intent(in)    :: psi
        double precision, intent(in)    :: axis1, axis2
        double precision, intent(in)    :: hl, hu
        integer, intent(in) :: ir

        double precision :: eazimuth

        eazimuth = pi2*dranu(ir + 2)

        block ! Pick up from velocity distribution
            integer :: addr0, addr1, addr2

            addr0 = nprs + npr(is)*dranu(ir)
            addr1 = nprs + npr(is)*dranu(ir + 1)

            vnorm = vnormc(addr0)*dcos(eazimuth) &
                    - vtangc(addr1)*dsin(eazimuth)
            vtang1 = vnormc(addr0)*dsin(eazimuth) &
                     + vtangc(addr1)*dcos(eazimuth)
        end block

        block
            integer :: j3
            j3 = nprs + npr(is)*dranu(ir + 3)
            vtang2 = vtangf(j3)
        end block

        rnorm = axis1 + radius*dcos(eazimuth)
        rtang1 = axis2 + radius*dsin(eazimuth)
        rtang2 = hl + (hu - hl)*dranu(ir + 5)
    end subroutine emit_from_cylindrical_surf3

end subroutine psuper3
