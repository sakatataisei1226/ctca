module m_probability_maxwell
    use m_vector
    use m_probability_base
    implicit none

    type, extends(t_PhaseProbability) :: t_MaxwellWithVelocityProbability
        double precision :: domains(2, 3)
        double precision :: vel_mean
        double precision :: vel_variance
        double precision :: vel_direction(3)
        double precision :: vel_perp1(3)
        double precision :: vel_perp2(3)

    contains
        procedure :: sample => MVProbability_sample
        procedure :: isvalid_in_local => MVProbability_isvalid_in_local
        procedure :: hitrate => MVProbability_hitrate
        procedure :: subdomain_rate => MVProbability_subdomain_rate
    end type

contains

    function new_MVProbability(domains, subdomains, &
                               vel_variance, vel_direction) result(obj)
        double precision, intent(in) :: domains(2, 3)
        double precision, intent(in) :: subdomains(2, 3)
        double precision, intent(in) :: vel_variance
        double precision, intent(in) :: vel_direction(3)
        type(t_MaxwellWithVelocityProbability) :: obj

        obj%domains = domains
        obj%subdomains = subdomains
        obj%vel_variance = vel_variance
        obj%vel_direction = vel_direction

        obj%vel_perp1 = perp_vec(vel_direction)
        obj%vel_perp2 = cross(obj%vel_direction, obj%vel_perp1)
    end function

    subroutine MVProbability_sample(self, rands, success)
        use wrapper
        class(t_MaxwellWithVelocityProbability) :: self
        double precision, intent(out) :: rands(6)
        logical, intent(out) :: success

        ! Upper and lower of self%domains.
        double precision :: dl(3), du(3)
        double precision :: v(3)

        double precision :: sigma

        double precision :: PI

        integer :: icon

        PI = acos(-1.0d0)

        call ranu0(rands, 6, icon)

        ! Assign position.
        dl(1:3) = max(self%domains(1, 1:3), self%subdomains(1, 1:3))
        du(1:3) = min(self%domains(2, 1:3), self%subdomains(2, 1:3))
        rands(1:3) = rands(1:3)*(du(1:3) - dl(1:3)) + dl(1:3)

        ! Assign maxwellian velocity.
        sigma = self%vel_variance

        ! v in (-inf, +inf)
        ! Probability density function:
        !   f(v) = 1/sqrt(2*pi*sigma^2)*exp(-v^2/(2*sigma^2))
        v(1) = 2.0d0*sigma*sqrt(-2.0d0*log(rands(4)))*cos(2.0d0*PI*rands(5))
        v(2) = 2.0d0*sigma*sqrt(-2.0d0*log(rands(4)))*sin(2.0d0*PI*rands(5))

        ! v in (0, +inf)
        ! Probability density function:
        !   f(v) = 1/sigma^2*v*exp(-v^2/(2*sigma^2))
        !   F(v) = 1 - exp(-v^2/(2*sigma^2))
        ! v(3) = sqrt(-2*sigma*sigma*log(rands(6)))
        v(3) = sqrt(-2*sigma*sigma*log(rands(6)))

        rands(4:6) = self%vel_direction*v(3) &
                     + self%vel_perp1*v(1) &
                     + self%vel_perp2*v(2)

        success = .true.
    end subroutine

    function MVProbability_isvalid_in_local(self) result(ret)
        class(t_MaxwellWithVelocityProbability) :: self
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

            if (ds == de .and. ds == sds) then
                ret = .false.
                return
            end if

            if (ds == de .and. ds == sde) then
                cycle
            end if

            if (de <= sds .or. sde <= ds) then
                ret = .false.
                return
            end if
        end do

        ret = .true.
    end function

    function MVProbability_hitrate(self) result(ret)
        class(t_MaxwellWithVelocityProbability) :: self
        double precision :: ret

        ret = 1.0d0
    end function

    function MVProbability_subdomain_rate(self) result(ret)
        class(t_MaxwellWithVelocityProbability) :: self
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

end module
