module m_probability_base
    implicit none

    !> Phase probability.
    type, abstract :: t_PhaseProbability
        !> Range of local subdomain.
        !> [[xl, xu], [yl, yu], [zl, zu]]
        double precision :: subdomains(2, 3)

    contains
        !> Sampling a phase that follows a probability distribution.
        procedure(PhaseProbability_sample), deferred :: sample

        !> Sampling a phase that follows a probability distribution (get a phase absolutely).
        procedure :: fsample => PhaseProbability_fsample

        !> Return if it needs to be considered in this subdomain.
        procedure(PhaseProbability_isvalid_in_local), deferred :: isvalid_in_local

        !> Return the hit rate of sampling.
        procedure(PhaseProbability_hitrate), deferred :: hitrate

        procedure(PhaseProbability_subdomain_rate), deferred :: subdomain_rate
    end type

    abstract interface
        subroutine PhaseProbability_sample(self, rands, success)
            import t_PhaseProbability
            class(t_PhaseProbability) :: self
            double precision, intent(out) :: rands(6)
            logical, intent(out) :: success
        end subroutine

        function PhaseProbability_isvalid_in_local(self) result(ret)
            import t_PhaseProbability
            class(t_PhaseProbability) :: self
            logical :: ret
        end function

        function PhaseProbability_hitrate(self) result(ret)
            import t_PhaseProbability
            class(t_PhaseProbability) :: self
            double precision :: ret
        end function

        function PhaseProbability_subdomain_rate(self) result(ret)
            import t_PhaseProbability
            class(t_PhaseProbability) :: self
            double precision :: ret
        end function
    end interface

    private
    public t_PhaseProbability

contains

    subroutine PhaseProbability_fsample(self, rands)
        class(t_PhaseProbability) :: self
        double precision, intent(out) :: rands(6)

        logical :: success

        success = .false.

        do while (.not. success)
            call self%sample(rands, success)
        end do
    end subroutine

end module
