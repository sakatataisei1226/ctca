module m_probability_grid
    use m_vector
    use m_probability_base
    implicit none

    !> Phase probability defined on the grid.
    type, extends(t_PhaseProbability) :: t_GridPhaseProbability
        !> Domain of definition of phase probability.
        !> [[xs, xe, nx], [ys, ye, ny], ..., [vzs, vze, nvz]]
        double precision :: domains(3, 6)

        !> Grid width.
        !> [nx, ny, nz, nvx, nvy, nvz]
        double precision :: gridwidth(6)

        !> Phase probability defined on the grid.
        double precision, allocatable :: probs(:, :, :, :, :, :)
        !> Max phase probability.
        double precision :: max_prob
        double precision :: mean_prob

        !> Emittion direction
        double precision :: direction(3)

    contains
        procedure :: sample => GridPhaseProbability_sample
        procedure :: isvalid_in_local => GridPhaseProbability_isvalid_in_local
        procedure :: hitrate => GridPhaseProbability_hitrate
        procedure :: subdomain_rate => GridPhaseProbability_subdomain_rate
        procedure :: interpolate => GridPhaseProbability_interpolate
        procedure :: load => GridPhaseProbability_load
        procedure :: get_x => GridPhaseProbability_get_z
        procedure :: get_y => GridPhaseProbability_get_y
        procedure :: get_z => GridPhaseProbability_get_z
        procedure :: get_vx => GridPhaseProbability_get_vx
        procedure :: get_vy => GridPhaseProbability_get_vy
        procedure :: get_vz => GridPhaseProbability_get_vz
    end type

    private
    public t_GridPhaseProbability
    public new_GridPhaseProbability

contains

    function new_GridPhaseProbability(domains, subdomains, direction) result(obj)
        double precision, intent(in) :: domains(3, 6)
        double precision, intent(in) :: subdomains(2, 3)
        double precision, optional, intent(in) :: direction(3)
        type(t_GridPhaseProbability) :: obj

        integer :: i
        double precision :: ps(6), pe(6)
        integer :: np(6)

        obj%domains = domains
        obj%subdomains = subdomains
        if (present(direction)) then
            obj%direction = direction
        else
            obj%direction = 0.0d0
        end if

        ps = domains(1, :)  ! start of domain
        pe = domains(2, :)  ! end of domain
        np = int(domains(3, :))  ! number of domain grids

        do i = 1, 6
            if (np(i) == 1) then
                obj%gridwidth(i) = 0d0; 
            else
                obj%gridwidth(i) = (pe(i) - ps(i))/max(1, np(i) - 1)
            end if
        end do

        allocate (obj%probs(0:np(1) - 1, &
                            0:np(2) - 1, &
                            0:np(3) - 1, &
                            0:np(4) - 1, &
                            0:np(5) - 1, &
                            0:np(6) - 1))
    end function

    subroutine GridPhaseProbability_load(self, filename)
        class(t_GridPhaseProbability) :: self
        character(len=*), intent(in) :: filename

        integer idf
        integer :: ivx, ivy, ivz
        double precision :: d
        double precision :: v(3)
        integer :: np(6)
        double precision :: probs_velocity_sum(3)

        ! Load datafiles
        open (newunit=idf, file=filename, form='unformatted', access='stream')
        read (idf) self%probs
        close (idf)

        self%probs = abs(self%probs)

        np(1:6) = int(self%domains(3, 1:6))

        ! Correct according to the emission direction
        if (norm2(self%direction) /= 0) then
            self%max_prob = 0.0d0
            self%mean_prob = 0.0d0
            do ivx = 0, np(4) - 1
                v(1) = self%get_vx(ivx)
                do ivy = 0, np(5) - 1
                    v(2) = self%get_vy(ivy)
                    do ivz = 0, np(6) - 1
                        v(3) = self%get_vz(ivz)

                        d = abs(dot(self%direction, v))

                        self%max_prob = max(self%max_prob, &
                                            maxval(self%probs(:, :, :, ivx, ivy, ivz))*d)
                        self%mean_prob = self%mean_prob &
                                         + sum(self%probs(:, :, :, ivx, ivy, ivz))*d
                        ! self%probs(:, :, :, ivx, ivy, ivz) = &
                        !     self%probs(:, :, :, ivx, ivy, ivz)*d
                    end do
                end do
            end do
            self%mean_prob = self%mean_prob/size(self%probs)
        else
            self%max_prob = maxval(self%probs)
            self%mean_prob = sum(self%probs)/size(self%probs)
        end if
    end subroutine

    subroutine GridPhaseProbability_sample(self, rands, success)
        use wrapper
        class(t_GridPhaseProbability) :: self
        double precision, intent(out) :: rands(6)
        logical, intent(out) :: success

        double precision :: ds(6), de(6)

        double precision :: rand
        integer :: i

        double precision :: r
        integer(kind=4) :: icon

        ds(1:3) = max(self%domains(1, 1:3), self%subdomains(1, 1:3))
        de(1:3) = min(self%domains(2, 1:3), self%subdomains(2, 1:3))

        ds(4:6) = self%domains(1, 4:6)
        de(4:6) = self%domains(2, 4:6)

        call ranu0(rands, 6, icon)
        rands = rands*(de - ds) + ds

        call ranu0(rand, 1, icon)
        rand = rand*self%max_prob

        r = self%interpolate(rands)
        if (rand < r) then
            success = .true.
        else
            success = .false.
        end if
    end subroutine

    function GridPhaseProbability_isvalid_in_local(self) result(ret)
        class(t_GridPhaseProbability) :: self
        logical :: ret

        integer :: i

        !> Domain range.
        double precision :: ds, de
        !> Subdomain range.
        double precision :: sds, sde

        do i = 1, 3
            ds = self%domains(1, i)
            de = self%domains(2, i)
            sds = self%subdomains(1, i)
            sde = self%subdomains(2, i)

            if (de < sds .or. sde < ds) then
                ret = .false.
                return
            end if
        end do

        ret = .true.
    end function

    function GridPhaseProbability_hitrate(self) result(ret)
        class(t_GridPhaseProbability) :: self
        double precision :: ret

        ret = self%mean_prob/self%max_prob
    end function

    function GridPhaseProbability_subdomain_rate(self) result(ret)
        class(t_GridPhaseProbability) :: self
        double precision :: ret

        double precision :: lx, ly, lz
        double precision :: slx, sly, slz

        lx = self%domains(2, 1) - self%domains(1, 1)
        ly = self%domains(2, 2) - self%domains(1, 2)
        lz = self%domains(2, 3) - self%domains(1, 3)

        slx = min(self%subdomains(2, 1), self%domains(2, 1)) &
              - max(self%subdomains(1, 1), self%domains(1, 1))
        sly = min(self%subdomains(2, 2), self%domains(2, 2)) &
              - max(self%subdomains(1, 2), self%domains(1, 2))
        slz = min(self%subdomains(2, 3), self%domains(2, 3)) &
              - max(self%subdomains(1, 3), self%domains(1, 3))

        if (lx == 0) then
            lx = 1
            slx = 1
        end if

        if (ly == 0) then
            ly = 1
            sly = 1
        end if

        if (lz == 0) then
            lz = 1
            slz = 1
        end if

        ret = (slx*sly*slz)/(lx*ly*lz)
    end function

    function GridPhaseProbability_interpolate(self, p) result(val)
        class(t_GridPhaseProbability) :: self
        double precision, intent(in) :: p(6)
        double precision :: val

        double precision :: d

        double precision :: ps(6)
        double precision :: dp(6)

        integer :: ip0(6)
        integer :: ip1(6)
        integer :: ipd(6)
        double precision :: lp(6)
        double precision :: rp0(6)
        double precision :: rp1(6)

        double precision :: probs(0:1, 0:1, 0:1, 0:1, 0:1, 0:1)

        integer :: i

        d = dot(self%direction, p(4:6))

        ps = self%domains(1, :)
        dp = self%gridwidth(:)

        do i = 1, 6
            if (abs(dp(i)) < 1d-10) then
                ip0(i) = 0
                ip1(i) = 0

                ! If dp=0, use-range of probs will be halved,
                ! so both of rp0 and rp1 be set to 1.
                rp0(i) = 1
                rp1(i) = 1
            else
                lp(i) = (p(i) - ps(i))/dp(i)

                ip0(i) = int(lp(i))
                ip1(i) = ip0(i) + 1

                rp0(i) = lp(i) - ip0(i)
                rp1(i) = 1 - rp0(i)
            end if
        end do

        ipd = ip1 - ip0
        probs(0:ipd(1), &
              0:ipd(2), &
              0:ipd(3), &
              0:ipd(4), &
              0:ipd(5), &
              0:ipd(6)) = self%probs(ip0(1):ip1(1), &
                                     ip0(2):ip1(2), &
                                     ip0(3):ip1(3), &
                                     ip0(4):ip1(4), &
                                     ip0(5):ip1(5), &
                                     ip0(6):ip1(6))

        probs(0, :, :, :, :, :) = probs(0, :, :, :, :, :)*rp1(1)
        probs(:, 0, :, :, :, :) = probs(:, 0, :, :, :, :)*rp1(2)
        probs(:, :, 0, :, :, :) = probs(:, :, 0, :, :, :)*rp1(3)
        probs(:, :, :, 0, :, :) = probs(:, :, :, 0, :, :)*rp1(4)
        probs(:, :, :, :, 0, :) = probs(:, :, :, :, 0, :)*rp1(5)
        probs(:, :, :, :, :, 0) = probs(:, :, :, :, :, 0)*rp1(6)
        probs(ipd(1), :, :, :, :, :) = probs(ipd(1), :, :, :, :, :)*rp0(1)
        probs(:, ipd(2), :, :, :, :) = probs(:, ipd(2), :, :, :, :)*rp0(2)
        probs(:, :, ipd(3), :, :, :) = probs(:, :, ipd(3), :, :, :)*rp0(3)
        probs(:, :, :, ipd(4), :, :) = probs(:, :, :, ipd(4), :, :)*rp0(4)
        probs(:, :, :, :, ipd(5), :) = probs(:, :, :, :, ipd(5), :)*rp0(5)
        probs(:, :, :, :, :, ipd(6)) = probs(:, :, :, :, :, ipd(6))*rp0(6)

        val = sum(probs(0:ipd(1), &
                        0:ipd(2), &
                        0:ipd(3), &
                        0:ipd(4), &
                        0:ipd(5), &
                        0:ipd(6)))*d
    end function

    function GridPhaseProbability_get_x(self, index) result(ret)
        class(t_GridPhaseProbability) :: self
        integer, intent(in) :: index
        double precision :: ret

        ret = self%domains(1, 1) + self%gridwidth(1)*index
    end function

    function GridPhaseProbability_get_y(self, index) result(ret)
        class(t_GridPhaseProbability) :: self
        integer, intent(in) :: index
        double precision :: ret

        ret = self%domains(1, 2) + self%gridwidth(2)*index
    end function

    function GridPhaseProbability_get_z(self, index) result(ret)
        class(t_GridPhaseProbability) :: self
        integer, intent(in) :: index
        double precision :: ret

        ret = self%domains(1, 3) + self%gridwidth(3)*index
    end function

    function GridPhaseProbability_get_vx(self, index) result(ret)
        class(t_GridPhaseProbability) :: self
        integer, intent(in) :: index
        double precision :: ret

        ret = self%domains(1, 4) + self%gridwidth(4)*index
    end function

    function GridPhaseProbability_get_vy(self, index) result(ret)
        class(t_GridPhaseProbability) :: self
        integer, intent(in) :: index
        double precision :: ret

        ret = self%domains(1, 5) + self%gridwidth(5)*index
    end function

    function GridPhaseProbability_get_vz(self, index) result(ret)
        class(t_GridPhaseProbability) :: self
        integer, intent(in) :: index
        double precision :: ret

        ret = self%domains(1, 6) + self%gridwidth(6)*index
    end function
end module
