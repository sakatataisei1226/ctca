module m_raycast
    use oh_type
    use paramt
    use allcom
    use finbound
    use m_vector, only: cross, dot, normalize
    use m_grad_ema
    implicit none

    private
    public generate_particles_with_raycast

contains

    subroutine generate_particles_with_raycast(camera, boundaries, is, nrays, npgen)
        type(t_ParallelCamera), intent(inout) :: camera
        type(t_BoundaryList), intent(in) :: boundaries
        integer, intent(in) :: is !> Species number (usually 1: electron, 2: ion, 3: pe).
        integer, intent(in) :: nrays
        integer, intent(out) :: npgen

        integer :: icon
        integer :: iray
        logical :: is_generated

        integer :: nprs
        double precision :: xl, xu, yl, yu, zl, zu
        xl = sdoms(1, 1, sdid(1) + 1); xu = sdoms(2, 1, sdid(1) + 1)
        yl = sdoms(1, 2, sdid(1) + 1); yu = sdoms(2, 2, sdid(1) + 1)
        zl = sdoms(1, 3, sdid(1) + 1); zu = sdoms(2, 3, sdid(1) + 1)

        if (is == 1) then
            nprs = 0
        else
            nprs = sum(npr(1:is - 1)) + 1
        end if

        npgen = 0
        do iray = 1, nrays
            call RANU0(dranu, 5, icon)
            call generate_particle_with_raycast(dranu(1:5), is_generated)
            if (is_generated) then
                npgen = npgen + 1
            end if
        end do

    contains

        subroutine generate_particle_with_raycast(rands, generated)
            double precision, intent(in) :: rands(5)
            logical, intent(out), optional :: generated

            type(t_Ray) :: ray
            type(t_HitRecord) :: hit_record

            double precision :: nc1(3), nc2(3) !> THe vectors perpendicular to hit_record%n

            type(oh_particle) :: pgen(2)
            double precision :: vgen(3)
            integer :: iboundary
            double precision :: v_maxwell(3)

            integer :: rid

            if (present(generated)) then
                generated = .false.
            end if

            ! Execute raycasting.
            ray = camera%generate_randray(rands(4:5))

            block
                type(t_HitRecord) :: tmp_hit_record

                hit_record%is_hit = .false.
                hit_record%t = 1d5
                do iboundary = 1, boundaries%nboundaries
                    tmp_hit_record = boundaries%boundaries(iboundary)%hit(ray)
                    if (tmp_hit_record%is_hit .and. &
                        (.not. hit_record%is_hit .or. tmp_hit_record%t < hit_record%t)) then
                        hit_record = tmp_hit_record
                    end if
                end do

            end block

            ! hit_record = boundaries%hit(ray)

            ! Return if the ray does not hit.
            if (.not. hit_record%is_hit) then
                return
            end if

            ! Return if the surface does not emit photoelectrons.
            if (hit_record%material%tag == -1) then
                return
            end if

            ! if (hit_record%position(1) < 0.0d0 .or. nx < hit_record%position(1) &
            !     .or. hit_record%position(2) < 0.0d0 .or. ny < hit_record%position(2) &
            !     .or. hit_record%position(3) < 0.0d0 .or. nz < hit_record%position(3)) then
            !     return
            ! end if

            if (hit_record%position(1) < xl .or. xu <= hit_record%position(1) &
                .or. hit_record%position(2) < yl .or. yu <= hit_record%position(2) &
                .or. hit_record%position(3) < zl .or. zu <= hit_record%position(3)) then
                return
            end if

            ! Find two vectors perpendicular to the normal vector.
            block
                double precision :: vtmp(3)
                vtmp = [1.0d0, 0.0d0, 0.0d0]
                if (abs(dot(vtmp, hit_record%n(:)) - 1.0d0) < 1d-3) then
                    vtmp = [0.0d0, 1.0d0, 0.0d0]
                end if

                nc1(:) = cross(vtmp(:), hit_record%n(:))
                call normalize(nc1(:))
                nc2(:) = cross(hit_record%n(:), nc1(:))
            end block

            ! Get particle position.
            pgen(1)%x = hit_record%position(1)
            pgen(1)%y = hit_record%position(2)
            pgen(1)%z = hit_record%position(3)

            ! Get particle velocity from half-maxwellian destribution.
            block
                integer(kind=4) :: addr0, addr1, addr2
                addr0 = nprs + npr(is)*rands(1)
                addr1 = nprs + npr(is)*rands(2)
                addr2 = nprs + npr(is)*rands(3)

                vgen(:) = vnormf(addr0)*hit_record%n(:) &
                          + vtangf(addr1)*nc1(:) &
                          + vtangf(addr2)*nc2(:)
            end block

            pgen(1)%vx = vgen(1)
            pgen(1)%vy = vgen(2)
            pgen(1)%vz = vgen(3)

            ! Set particle infomation.
            pgen(1)%spec = is
            pgen(1)%pid = 0

            ! Copy particle infomation for the patricle to be accumulated as anti-particle.
            pgen(2) = pgen(1)

            ! Shift particle position for convenience of inner boundary collision detection.
            pgen(2)%x = pgen(2)%x + hit_record%n(1)*1d-5
            pgen(2)%y = pgen(2)%y + hit_record%n(2)*1d-5
            pgen(2)%z = pgen(2)%z + hit_record%n(3)*1d-5

            ! Get the subdomain ID where the particle should exist.
            rid = oh3_map_particle_to_subdomain(pgen(1)%x, pgen(1)%y, pgen(1)%z)
            pgen(1)%nid = rid

            rid = oh3_map_particle_to_subdomain(pgen(2)%x, pgen(2)%y, pgen(2)%z)
            pgen(2)%nid = rid

            ! Inject paticles
            pgen(1)%preside = OH_PCL_TO_BE_ACCUMULATED_AS_ANTIPCL
            call oh2_inject_particle(pgen(1))

            pgen(2)%preside = OH_PCL_INJECTED
            call oh2_inject_particle(pgen(2))

            if (present(generated)) then
                generated = .true.
            end if
        end subroutine

    end subroutine

end module
