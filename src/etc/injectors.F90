#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"

module m_particle_injectors
    use oh_type
    use m_phase_probabilities
    use m_vector

    implicit none

    double precision, parameter :: EPSILON = 1d-5
    integer(kind=4), external :: oh3_map_particle_to_subdomain

    type :: t_ParticleInjector

        class(t_PhaseProbability), pointer :: probability
        double precision :: ninj_per_step
        double precision :: ntries_per_step_in_subdomain
        double precision :: inject_direction(3)
        double precision :: initial_acceleration
        integer :: is
        integer :: sdid
        logical :: remain_antiparticles

    contains

        procedure :: inject => InjectionFace_inject

    end type

    interface t_ParticleInjector
        module procedure :: new_ParticleInjector
    end interface

    private
    public t_ParticleInjector
    public new_ParticleInjector

contains

    function new_ParticleInjector(probability, &
                                  ninj_per_step, &
                                  inject_direction, &
                                  initial_acceleration, &
                                  is, sdid, &
                                  remain_antiparticles) result(obj)
        class(t_PhaseProbability), pointer, intent(in) :: probability
        double precision, intent(in) :: ninj_per_step
        double precision, intent(in) :: inject_direction(3)
        double precision, intent(in) :: initial_acceleration
        integer, intent(in) :: is
        integer, intent(in) :: sdid
        logical, intent(in) :: remain_antiparticles
        type(t_ParticleInjector) :: obj

        double precision :: hitrate, sdrate

        obj%probability => probability
        obj%ninj_per_step = ninj_per_step

        hitrate = probability%hitrate()
        sdrate = probability%subdomain_rate()
        obj%ntries_per_step_in_subdomain = ninj_per_step/hitrate*sdrate

        obj%inject_direction = normalized(inject_direction)
        obj%initial_acceleration = initial_acceleration

        obj%is = is
        obj%sdid = sdid

        obj%remain_antiparticles = remain_antiparticles
    end function

    !> 粒子を生成する
    subroutine InjectionFace_inject(self, ustep)
        use wrapper

        class(t_ParticleInjector) :: self
        double precision, intent(in) :: ustep

        type(oh_particle) :: pinj(2)
        double precision :: phase_sampled(6)
        integer :: i

        integer :: ntries
        double precision :: rand
        logical :: success
        integer :: icon
        double precision :: ustep_valid

        double precision :: velocity(3)
        double precision :: initial_move(3)
        double precision :: d

        double precision :: generate_offsets(3)

        ntries = int(self%ntries_per_step_in_subdomain*ustep)

        ! Map the remaining minority number to 0 or 1 randomly.
        call ranu0(rand, 1, icon)
        if (rand < self%ntries_per_step_in_subdomain*ustep - ntries) then
            ntries = ntries + 1
        end if

        ! When judging a collision with an internal boundary,
        ! if the position of the particle one step before is above the boundary,
        ! the collision is not judged, so the particle is generated with a slight shift.
        generate_offsets = self%inject_direction*EPSILON

        do i = 1, ntries
            call self%probability%sample(phase_sampled, success)

            if (.not. success) then
                cycle
            end if

            call ranu0(rand, 1, icon)
            ustep_valid = ustep*rand

            pinj(1)%x = phase_sampled(1)
            pinj(1)%y = phase_sampled(2)
            pinj(1)%z = phase_sampled(3)
            pinj(1)%vx = phase_sampled(4) + self%inject_direction(1)*self%initial_acceleration*ustep_valid
            pinj(1)%vy = phase_sampled(5) + self%inject_direction(2)*self%initial_acceleration*ustep_valid
            pinj(1)%vz = phase_sampled(6) + self%inject_direction(3)*self%initial_acceleration*ustep_valid

            pinj(1)%nid = self%sdid
            pinj(1)%spec = self%is
            pinj(1)%pid = 0

            if (self%remain_antiparticles) then
                pinj(2) = pinj(1)

                pinj(2)%preside = OH_PCL_TO_BE_ACCUMULATED_AS_ANTIPCL
                call oh2_inject_particle(pinj(2))
            end if

            ! Particles are generated in continuous time,
            ! so they are moved by generation time (random).
            velocity(:) = (/pinj(1)%vx, pinj(1)%vy, pinj(1)%vz/)
            d = dot(velocity, self%inject_direction)
            initial_move = d*self%inject_direction*ustep_valid

            pinj(1)%x = (pinj(1)%x + generate_offsets(1)) + initial_move(1)
            pinj(1)%y = (pinj(1)%y + generate_offsets(2)) + initial_move(2)
            pinj(1)%z = (pinj(1)%z + generate_offsets(3)) + initial_move(3)

            pinj(1)%preside = OH_PCL_ALIVE
            ! pinj(1)%nid = oh3_map_particle_to_subdomain(pinj(1)%x,pinj(1)%y,pinj(1)%z)
            call oh2_inject_particle(pinj(1))
        end do
    end subroutine
end module
