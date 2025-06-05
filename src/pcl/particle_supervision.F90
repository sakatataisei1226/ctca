#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"

module m_particle_supervision
    use m_phase_probabilities
    use m_particle_injectors
    use allcom
    use paramt
    use m_ohhinfo
    use m_str
    use m_logging
    use m_vector
    implicit none

    ! ==== Namelist parameters 'plasma.inp' ====

    !> True if use this module's particle supervision.
    logical :: use_pinj = .false.

    !> True when the particle generation leaves an antiparticle
    !> (particle with opposite charge) at the generation position.
    integer :: remain_antipcl(inepl) = 0

    !> Direction of particle emission.
    !> Valid only if nemd(iepl) = 0
    double precision :: emit_directions(inepl, 3)

    !> Acceleration applied during particle creation
    !> This is applied for the step time (ustep*rand) at which the generated particle becomes active.
    !>
    !> pcl%velocity += acc * ustep*rand
    double precision :: initial_accelerations(inepl) = 0

    !> True if use t_GridPhaseProbability.
    character(len=200) :: interp_param_files(inepl)

    namelist /pclinj/ &
        use_pinj, &
        remain_antipcl, emit_directions, &
        initial_accelerations, &
        interp_param_files
    ! ==========================================

    ! ==== Namelist parameters '<interp>.inp' ====
    character(len=200) :: interp_filename
    double precision :: interp_domain(3, 6)

    namelist /interp/ interp_filename, interp_domain
    ! ============================================

    !> Parameter filename of t_GridPhaseProbability.
    character(len=200) :: interp_filenames(inepl)

    !> Domain of t_GridPhaseProbability.
    double precision :: interp_domains(inepl, 3, 6)

    integer, parameter :: MAX_NINJECTORS = 30
    integer :: ninjectors = 0

    type(t_ParticleInjector), target :: injectors(inepl)

    private
    public use_pinj
    public psuper_read_namelist
    public psuper_initialize
    public psuper_inject_particles
    public psuper_show_settings

contains

    subroutine psuper_read_namelist(namelist_filename)
        character(len=*), intent(in) :: namelist_filename

        integer :: idf

        integer :: iepl
        integer :: iosta

        if (myid == 0) print *, 'call psuper_read_namelist...'

        interp_param_files = ''
        interp_filenames = ''

        open (newunit=idf, file=namelist_filename)
        read (idf, nml=pclinj, IOSTAT=iosta)
        close (idf)

        if (myid .eq. 0 .and. iosta .eq. -1) then
            print *, "Warning.Input: nml=pclinj not found"
            return
        end if

        if (myid == 0) print *, '  use_pinj:', use_pinj

        if (.not. use_pinj) then
            if (myid == 0) print *, 'end psuper_read_namelist'
            return
        end if

        do iepl = 1, inepl
            if (interp_param_files(iepl) /= '') then
                if (myid == 0) print *, '  read interp_param_file: '//trim(interp_param_files(iepl))
                open (newunit=idf, file=interp_param_files(iepl))
                read (idf, nml=interp)
                close (idf)

                interp_filenames(iepl) = interp_filename
                interp_domains(iepl, :, :) = interp_domain(:, :)
                if (myid == 0) print *, '  end read interp_param_file: '//trim(interp_param_files(iepl))
            end if
        end do

        if (myid == 0) print *, 'end psuper_read_namelist'
    end subroutine

    subroutine psuper_initialize
        integer :: is
        integer :: iepl, iepl_start, iepl_end

        if (myid == 0) print *, 'start psuper_initialize...'
        call ohhinfo_update(sdoms(:, :, sdid(1) + 1))

        iepl_end = 0
        do is = 1, nspec
            iepl_start = iepl_end + 1
            iepl_end = iepl_end + nepl(is)
            do iepl = iepl_start, iepl_end
                if (myid == 0) print *, '  start add_emission:', iepl
                if (interp_param_files(iepl) == '') then
                    call add_inner_emmision(iepl, is)
                else
                    call add_interp_emission(iepl, is)
                end if
                if (myid == 0) print *, '  end add_emission:', iepl
            end do
        end do

        if (myid == 0) print *, 'end psuper_initialize'
    end subroutine

    subroutine add_inner_emmision(iepl, is)
        integer, intent(in) :: iepl
        integer, intent(in) :: is

        double precision :: domains(2, 3)
        type(t_MaxwellWithVelocityProbability), pointer :: mprob
        class(t_PhaseProbability), pointer :: pprob
        type(t_ParticleInjector) :: injector

        double precision :: vth
        double precision :: direction(3)
        double precision :: ninj_per_step
        double precision :: area
        double precision :: width
        double precision :: subdomains(2, 3)

        integer :: i

        domains = reshape( &
                  (/xmine(iepl), xmaxe(iepl), &
                    ymine(iepl), ymaxe(iepl), &
                    zmine(iepl), zmaxe(iepl)/), &
                  (/2, 3/))
        vth = path(is)
        direction = to_direction(iepl)
        subdomains = sdoms(:, :, sdid(1) + 1)

        area = 1.0d0
        do i = 1, 3
            width = domains(2, i) - domains(1, i)
            if (width == 0) then
                width = 1.0d0
            end if
            area = area*width
        end do
        area = area*dr*dr

        allocate (mprob)
        mprob = new_MVProbability( &
                domains, &
                subdomains, &
                vth, &
                direction)

        if (mprob%isvalid_in_local()) then
            pprob => mprob

            ninj_per_step = abs( &
                            (curfs(iepl)*renj) &
                            /abs(q(is)) &
                            *area &
                            )

            injector = new_ParticleInjector( &
                       pprob, &
                       ninj_per_step, &
                       direction, &
                       initial_accelerations(iepl)*rena, &
                       is, &
                       sdid(1), &
                       remain_antipcl(iepl) == 1 &
                       )
            ninjectors = ninjectors + 1
            injectors(ninjectors) = injector
        end if
    end subroutine

    subroutine add_interp_emission(iepl, is)
        integer, intent(in) :: iepl
        integer, intent(in) :: is

        double precision :: domains(3, 6)
        double precision :: subdomains(2, 3)
        double precision :: direction(3)

        type(t_GridPhaseProbability), pointer :: gprob
        class(t_PhaseProbability), pointer :: pprob
        type(t_ParticleInjector) :: injector

        double precision :: ninj_per_step
        double precision :: area
        double precision :: width

        integer :: i

        domains = interp_domains(iepl, :, :)
        domains(1:2, 1:3) = domains(1:2, 1:3)*renr
        domains(1:2, 4:6) = domains(1:2, 4:6)*renv
        direction = to_direction(iepl)
        subdomains = sdoms(:, :, sdid(1) + 1)

        area = 1.0d0
        do i = 1, 3
            width = domains(2, i) - domains(1, i)
            if (width == 0) then
                width = 1.0d0
            end if
            area = area*width
        end do

        allocate (gprob)
        gprob = new_GridPhaseProbability(domains, subdomains, direction)
        if (gprob%isvalid_in_local()) then
            call gprob%load(interp_filenames(iepl))
            pprob => gprob

            ninj_per_step = abs( &
                            (curfs(iepl)*renj) &
                            /abs(q(is)) &
                            *area &
                            )

            ninjectors = ninjectors + 1
            injector = new_ParticleInjector( &
                       pprob, &
                       ninj_per_step, &
                       direction, &
                       initial_accelerations(iepl)*rena, &
                       is, &
                       sdid(1), &
                       remain_antipcl(iepl) == 1)
            injectors(ninjectors) = injector
        end if
    end subroutine

    function to_direction(iepl) result(ret)
        integer, intent(in) :: iepl
        double precision :: ret(3)

        ret = 0.0d0

        if (nemd(iepl) == 0) then
            ret(1:3) = emit_directions(iepl, 1:3)
        else if (nemd(iepl) == 1) then
            ret(1) = 1.0d0
        else if (nemd(iepl) == -1) then
            ret(1) = -1.0d0
        else if (nemd(iepl) == 2) then
            ret(2) = 1.0d0
        else if (nemd(iepl) == -2) then
            ret(2) = -1.0d0
        else if (nemd(iepl) == 3) then
            ret(3) = 1.0d0
        else if (nemd(iepl) == -3) then
            ret(3) = -1.0d0
        end if
    end function

    subroutine psuper_inject_particles(ustep)
        double precision, intent(in) :: ustep

        integer :: i

        do i = 1, ninjectors
            call injectors(i)%inject(ustep)
        end do
    end subroutine

    subroutine psuper_show_settings
        double precision :: ninj_per_step
        double precision :: ntries_per_step_in_subdomain
        integer :: i

        ntries_per_step_in_subdomain = 0
        ninj_per_step = 0
        do i = 1, ninjectors
            ninj_per_step = ninj_per_step + injectors(i)%ninj_per_step
            ntries_per_step_in_subdomain = ntries_per_step_in_subdomain + injectors(i)%ntries_per_step_in_subdomain
            call logging_debuglog('ninj/step ('//str(i)//'): '//str(injectors(i)%ninj_per_step))
            call logging_debuglog('hitrate ('//str(i)//'): '//str(injectors(i)%probability%hitrate()))
            call logging_debuglog('sdrate ('//str(i)//'): '//str(injectors(i)%probability%subdomain_rate()))
        end do

        print *, 'sum of ninj_per_step: '//str(ninj_per_step)
        print *, 'ninjectors: '//str(ninjectors)
        print *, 'ntries_per_step_in_subdomain: '//str(ntries_per_step_in_subdomain)

        if (myid == 0) print *, 'particle_injection activated'
    end subroutine

end module
