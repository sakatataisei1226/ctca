!***************************************************************
! ElectroMagnetic Spacecraft Environment Simulator named EMSES
! ~      ~        ~          ~           ~
!      load-balanced by One-handed Help (OhHelp) Algorithm
!
!           controle program for emses
!
!***************************************************************
module m_esses
    !
    ! **** "esses" is an electrostatic version of emses ****
    !
    !---      Three Dimensional ElectroMagnetic Simulation Code     --
    !   This code is based on
    !        Kyoto university ElectroMagnetic Particle cOde (KEMPO)
    !
    !========================= History of kempo =========================
    !
    !       Programed by H.Matsumoto          July, 1980
    !                 by H.tahara             February, 1981
    !                 by I.Ayukawa            February, 1982
    !                 by M.Ohashi             April, 1982
    !                 by Y.Omura              September, 1982
    !                 by Y.Omura, K.Fukuchi
    !                          and T.Yamada   September, 1983
    !           Optimized for vector processor
    !                 by Y.Omura              February, 1984
    !                 by Y.Omura, T.Kimura    June, 1984
    !           Revised as version 7 & 8
    !                 by Y.Omura and N.Komori November, 1985
    !           Revised as version 9
    !                 by Y.Omura & T.Tanaka   April, 1986
    !           Data output format is revised as version 10
    !                 by Y.Omura & K.Inagaki  April,1986
    !           Velocity distribution diagnostics is added as version 11
    !                 by Y.Omura              July, 1986
    !
    !           Expanded to three dimensional system
    !                 by H.Yokoyama           January, 1992
    !           Free boundary is added as version 2 (3D)
    !                 by M.Yamane             May, 1993
    !
    !========================= History of emses =========================
    !
    !       Programed by Y.Miyake             October, 2008
    !           OhHelp Library (version 0.9) is developed
    !                 by H.Nakashima          August, 2009
    !           OhHelp is applied to ESSES
    !                 by Y.Miyake             July, 2010
    !
    !=====================================================================
    !
    !            Radio Atomospheric Science Center (-2004)
    !        Research Institute for Sustainable Humanosphere (2004-)
    !          Academic Center for Computing and Media Studies
    !              Kyoto University, Kyoto, Japan.
    !
    !
    !-------------------- parameter and common blocks
    use oh_type
    use paramt
    use allcom
    use m_grad_ema
    use m_ohhinfo
    use m_particle_supervision
    use m_logging
    use m_str
    use mpi
    use m_conductivity

#include "oh_stats.h"
#define MCW local_comm

    implicit none

    private
    public esses_run
    public esses_initialize
    public esses_mainstep
    public esses_finalize

contains

    !> Run electric static simulation.
    subroutine esses_run
        integer(kind=4) :: ustep
        integer(kind=4) :: mpierr
        integer(kind=8) :: mnmtotal = 0
        integer(kind=8) :: mnmlocal = 0
        real(kind=8) :: smltime, emltime
        real(kind=8) :: performance

        call MPI_Init(mpierr)
        local_comm=MPI_COMM_WORLD
        call MPI_Comm_size(MCW, nnode, mpierr)
        call MPI_Comm_rank(MCW, myid, mpierr)
        call oh1_fam_comm(MPI_COMM_WORLD)

        call hdfinit()

        emflag = 0
        call input

        call esses_initialize

        call MPI_Barrier(MCW, mpierr)
        smltime = MPI_Wtime()

        ! Main loop
        do istep = 1, nstep
            ! Update simulation time
            t = t + dt
            itime = istep

            if (istep /= nstep) then
                ustep = 2
            else
                ustep = 1
            end if

            if (myid == 0) then
                write (6, *) '**** step --------- ', itime
            end if

            call esses_mainstep(ustep)

            mnmlocal = mnmlocal + sum(totalp(:, :))
        end do

        call MPI_Barrier(MCW, mpierr)
        emltime = MPI_Wtime()

        call esses_finalize

        call MPI_Reduce(mnmlocal, mnmtotal, 1, &
                        MPI_INTEGER8, MPI_SUM, 0, &
                        MCW, mpierr)
        if (myid == 0) then
            print *, "time main loop", (emltime - smltime)
            print *, "total N. of processed part.", mnmtotal

            performance = dble(mnmtotal)/(emltime - smltime)*1.0d-6
            print *, "performance   ", performance, "Mpart/s"
        end if

        call hdffinalize()
        call MPI_Finalize(mpierr)
    end subroutine

    !> Initialize simulation.
    subroutine esses_initialize
        integer(kind=4) :: mpierr

        call logging_initialize

        ! Initialization
        call logging_debuglog("start ohinit")
        call ohinit
        call logging_debuglog("end ohinit")
        call logging_debuglog("start inital")
        call inital
        call logging_debuglog("end inital")

        ! Particle initialization
        call logging_debuglog("start inipcl")
        call inipcl
        call logging_debuglog("end inipcl")
        call logging_debuglog("start inipin")
        call inipin
        call logging_debuglog("end inipin")

        ! Field initialization
        call logging_debuglog("start inifld")
        call inifld
        call logging_debuglog("end inifld")

        call logging_debuglog("start grad_ema_init")
        call grad_ema_init
        call logging_debuglog("end grad_ema_init")

        ! Loadbalancing
        currmode = esses_loadbalance(0, 0, .true.)

        ! Variable unit conversion
        call rescal

        ! Calculate capacity matrix
        call getccl
        call capmtx

        ! Charge neutrality
        call logging_debuglog("start chgntr")
        call chgntr
        call logging_debuglog("end chgntr")

        ! r(0.5) to r(0.5+m/2)
        call logging_debuglog("start positn")
        if (jobnum(1) == 1) then
            call positn(1, 1)
            if (sdid(2) >= 0) then
                call positn(2, 1)
            end if

            if (.not. use_pinj) then
                rhobk(:, :, :, :, 1:2) = 0d0
                rhobksp(:, :, :, :, 1:2, :) = 0d0

                call psuper2(1)
            end if
        end if
        call logging_debuglog("end positn")

        ! Particle supervisor 3
        call logging_debuglog("Read psuper's namelist")
        call psuper_read_namelist('plasma.inp')
        if (use_pinj) then
            ! Initialize particle injection
            call logging_debuglog('Psuper initialize')
            call psuper_initialize
            call psuper_show_settings
            call logging_debuglog('Start psuper inject partciles')
            call psuper_inject_particles(2.0d0)
            call logging_debuglog('End psuper inject partciles')
        else
            call logging_debuglog('Call super3')
            call psuper3(2)
        end if

        ! Loadbalancing
        call logging_debuglog('Start loadbalance')
        currmode = esses_loadbalance(currmode, 0, .false.)
        call logging_debuglog('End loadbalance')

        ! Gather particle charge
        call logging_debuglog('Start charge')
        call esses_charge  ! rho(0.5+m/2)

        call logging_debuglog('Update rhobk')
        call grad_ema_update_rhobk  ! rhobk(:,:,:,:,3) = rhobk(:,:,:,:,3) + rhobk(:,:,:,:,1)

        call logging_debuglog('rho += rhobk')
        call add_background_charge(1)  ! rho += rhobk

        call logging_debuglog('Fsmask')
        ! Charge density masking
        if (imask(5) == 1) then
            call fsmask(1, 4)
        end if

        ! e(0.5+m/2):ES
        if (npc >= 1) call surchg

        call logging_debuglog('Update electric static field (esfld1)')
        call esfld1(0)
        if (currmode /= OHH_PRIMARY_MODE) then
            call oh3_bcast_field(phi(1, 0, 0, 0, 1), phi(1, 0, 0, 0, 2), FPH)
        end if

        call logging_debuglog('Update electric static field (esfld2)')
        call esfld2(1)
        if (sdid(2) >= 0) then
            call esfld2(2)
        end if

        call logging_debuglog('Poichk')
        call poichk(1)

        ! Check of parameters
        call logging_debuglog('Check parameters')
        call chkprm

        ! Initialize diagnostics
        call logging_debuglog('Initialize diagnostics')
        call energy

        ! Run particle sort
        if (isort(1) /= 0) then
            call logging_debuglog('Start partcile sort')
            call spsort(1)
            if (sdid(2) >= 0) then
                call spsort(2)
            end if
            call logging_debuglog('End particle sort')
        end if

        ! Initial diagnostics
        call MPI_Barrier(MCW, mpierr)
        call frmttd(1)
        call digset

        call logging_debuglog('start Hdfdig')
        call logging_debuglog('Hdfdig(0, 0)')
        call hdfdig(0, 0)
        call logging_debuglog('Hdfdig(1, 0)')
        call hdfdig(1, 0)
        call logging_debuglog('end Hdfdig')
    end subroutine

    !> Run main step.
    subroutine esses_mainstep(ustep)
        integer, intent(in) :: ustep

        !> Set 1 if start data output, otherwise 0.
        integer :: fajdg

        call logging_debuglog('Start mainstep')

        if (ijdiag == 0 .or. istep < hdfdigstart) then
            fajdg = 0
        else
            fajdg = 1
        end if

        call logging_debuglog('Start particle supervisor')

        ! particle supervisor 1
        if (use_pinj) then
            rhobk(:,:,:,:,1:2) = 0.0d0
            if (fajdg == 1) then
                aj(:, :, :, :, 1:2) = 0.0d0
                ajdg(:, :, :, :, 1:2) = 0.0d0
            end if

            ! Loadbalancing
            ! call logging_debuglog('Loadbalance')
            ! currmode = esses_loadbalance(currmode, 1, .false.)
            ! call logging_debuglog('Totalp: '//str(sum(totalp(:, :))))
        else
            rhobk(:, :, :, :, 1:2) = 0d0
            rhobksp(:, :, :, :, 1:2, :) = 0d0

            call logging_debuglog('  Start psuper12')
            if (fajdg == 0) then
                call psuper2(ustep)
            else
                call psuper1(ustep, fajdg)
            end if
            call logging_debuglog('  End psuper12')
        end if

        call logging_debuglog('Start primary psolve')

        ! relocation of field pe and pb
        call grdcrt(1)

        ! v(n+0.5-m) to v(n+0.5)
        ! r(n+0.5-m/2) to r(n+0.5+m/2)
        ! rho(n+0.5+m/2)_main
        if (fajdg == 0) then
            call logging_debuglog('  Psolve2')
            call psolve2(1, ustep)
        else
            call logging_debuglog('  Psolve1')
            call psolve1(1, ustep, fajdg)
        end if
        call logging_debuglog('End primary psolve')

        if (sdid(2) >= 0) then
            call logging_debuglog('Start secondary psolve')
            call grdcrt(2)

            if (fajdg == 0) then
                call psolve2(2, ustep)
            else
                call psolve1(2, ustep, fajdg)
            end if
            call logging_debuglog('End secondary psolve')
        end if

        ! j(n+0.5)_diag
        if (fajdg == 1) then
            call logging_debuglog('Exchange borders')
            if (currmode /= OHH_PRIMARY_MODE) then
                call logging_debuglog('  Start reduce field')
                call oh3_allreduce_field(aj(1, 0, 0, 0, 1), aj(1, 0, 0, 0, 2), FAJ)
                call logging_debuglog('  End reduce field')
            end if
            call logging_debuglog('  Exchange')
            call oh3_exchange_borders(aj(1,0,0,0,1),aj(1,0,0,0,2),CAJ,0)
            call logging_debuglog('  Add boundary current')
            call add_boundary_current(1)

            if (currmode /= OHH_PRIMARY_MODE) then
                call logging_debuglog('  Start reduce field')
                call oh3_allreduce_field(ajdg(1,0,0,0,1),ajdg(1,0,0,0,2),FJD)
                call logging_debuglog('  End reduce field')
            end if
            call logging_debuglog('  Exchange')
            call oh3_exchange_borders(ajdg(1, 0, 0, 0, 1), ajdg(1, 0, 0, 0, 2), CJD, 0)
            call logging_debuglog('  Add boundary current')
            call add_boundary_current2(1)
        end if

        ! Particle supervisor 3 [INSSERT generation]
        if (use_pinj) then
            call logging_debuglog('Start particle injection by psuper module')
            call psuper_inject_particles(2.0d0)
            call logging_debuglog('End particle injection by psuper module')
        else
            call logging_debuglog('Start particle injection by psuper3')
            call psuper3(2)
            call logging_debuglog('End particle injection by psuper3')
        end if

        ! Loadbalancing
        call logging_debuglog('Loadbalance')
        currmode = esses_loadbalance(currmode, 1, .false.)
        call logging_debuglog('Totalp: '//str(sum(totalp(:, :))))

        ! particle sort
        if (isort(1) /= 0 &
            .and. mod(istep, isort(1)) == 0 &
            .and. istep /= nstep) then
            call logging_debuglog('Particle sort')
            call spsort(1)

            if (sdid(2) >= 0) then
                call spsort(2)
            end if
        end if

        call logging_debuglog('Start charge')

        ! Gather particle charge
        call esses_charge  ! rho(n+0.5+m/2)_sub

        call logging_debuglog('Update rhobk')
        
        call update_conductive_rhobk(2d0)

        call grad_ema_update_rhobk  ! rhobk(:,:,:,:,3) = rhobk(:,:,:,:,3) + rhobk(:,:,:,:,1)

        call logging_debuglog('Add background charge')
        call add_background_charge(1)

        call logging_debuglog('End charge')

        ! Charge density masking
        if (imask(5) == 1) call fsmask(1, 4)

        ! Correction of e
        if (npc >= 1) call surchg
        call esfld1(0)
        if (currmode /= OHH_PRIMARY_MODE) then
            call oh3_bcast_field(phi(1, 0, 0, 0, 1), phi(1, 0, 0, 0, 2), FPH)
        end if

        call esfld2(1)
        if (sdid(2) >= 0) then
            call esfld2(2)
        end if
        call poichk(1)

        call logging_debuglog('Start fsmasking')
        ! e-field masking
        if (imask(2) == 1) then
            call fsmask(1, 1)
            if (sdid(2) >= 0) then
                call fsmask(2, 1)
            end if
        end if

        if (imask(3) == 1) then
            call fsmask(1, 2)
            if (sdid(2) >= 0) then
                call fsmask(2, 2)
            end if
        end if
        call logging_debuglog('End fsmasking')

        ! diagnostics
        call logging_debuglog('Diagnostics')
        call energy

        call frmttd(ustep)

        call hdfdig(1, istep)
    end subroutine

    !> Finalize simulation.
    subroutine esses_finalize
        integer :: mpierr

        currmode = oh3_transbound(currmode, 1)
        call hdfdig(2, istep)

        ! Save data for next job
        if (jobnum(2) == 1 .and. myid == 0) then
            call system('mkdir SNAPSHOT1')
        end if

        call MPI_Barrier(MCW, mpierr)
        call save_esdat
    end subroutine

    !> Run loadbalancing.
    function esses_loadbalance(mode, stats, need_primary_nborps) result(newmode)
        integer(kind=4), intent(in) :: mode
        integer(kind=4), intent(in) :: stats
        logical, intent(in) :: need_primary_nborps
        integer :: newmode

        integer(kind=4) :: mpierr

        newmode = oh3_transbound(mode, stats)

        nphgram(sdid(1) + 1, 1:nspec, 1) = totalp(1:nspec, 1)
        if (sdid(2) >= 0) then
            nphgram(sdid(2) + 1, 1:nspec, 2) = totalp(1:nspec, 2)
        end if

        if (need_primary_nborps) then
            call create_nborps(1)
        end if

        if (newmode == OHH_REQUIRES_BCAST) then
            call create_nborps(2)
            call oh3_bcast_field(eb(1, 0, 0, 0, 1), eb(1, 0, 0, 0, 2), FEB)
            call oh3_bcast_field(colf(1, 0, 0, 0, 1), colf(1, 0, 0, 0, 2), FPH)
            newmode = OHH_SECONDARY_MODE
        end if

        gcount(1)%globalp(1:nspec, :) = totalp(1:nspec, :)
        call MPI_Allreduce(gcount(1), gcount(2), 1, &
                           mpi_type_count, mpi_sum_count, &
                           MCW, mpierr)
    end function

    !> Gather particle charges.
    subroutine esses_charge()
        integer(kind=4) :: is
        integer(kind=4) :: ps
        logical :: require_output_data
        integer :: ierr

        require_output_data = ifdiag /= 0 .and. (mod(istep, ifdiag) == 0 .or. daverg >= 1)
        if (require_output_data) then
            call charge(1, 2)
            if (sdid(2) >= 0) then
                call charge(2, 2)
            end if
        else
            call charge(1, 0)
            if (sdid(2) >= 0) then
                call charge(2, 0)
            end if
        end if

        call logging_debuglog('Reduce rho')

        if (currmode /= OHH_PRIMARY_MODE) then
            call oh3_reduce_field(rho(1, 0, 0, 0, 1), rho(1, 0, 0, 0, 2), FRH)
        end if
        call logging_debuglog('Exchange rho')
        call oh3_exchange_borders(rho(1, 0, 0, 0, 1), rho(1, 0, 0, 0, 2), CRH, 0)
        call logging_debuglog('Add boundary rho')
        call add_boundary_charge(rho(:, :, :, :, 1), &
                                 1, 1, &
                                 sdoms(:, :, sdid(1) + 1), &
                                 bounds(:, :, sdid(1) + 1), &
                                 ctypes(:, :, 1, CRH), &
                                 fsizes(:, :, FRH), myid)

        call logging_debuglog('Reduce rhobk')
        if (currmode /= OHH_PRIMARY_MODE) then
            call oh3_reduce_field(rhobk(1, 0, 0, 0, 1), rhobk(1, 0, 0, 0, 2), FRH)
        end if
        call logging_debuglog('Exchange rhobk')
        call oh3_exchange_borders(rhobk(1, 0, 0, 0, 1), rhobk(1, 0, 0, 0, 2), CRH, 0)
        call logging_debuglog('Add boundary rhobk')
        call add_boundary_charge(rhobk(:, :, :, :, 1), &
                                 1, 1, &
                                 sdoms(:, :, sdid(1) + 1), &
                                 bounds(:, :, sdid(1) + 1), &
                                 ctypes(:, :, 1, CRH), &
                                 fsizes(:, :, FRH), myid)

        do is = 1, nspec
            if (currmode /= OHH_PRIMARY_MODE) then
                call oh3_reduce_field(rhobksp(1, 0, 0, 0, 1, is), rhobksp(1, 0, 0, 0, 2, is), FRH)
            end if
            call logging_debuglog('Exchange rhobksp '//str(is))
            call oh3_exchange_borders(rhobksp(1, 0, 0, 0, 1, is), rhobksp(1, 0, 0, 0, 2, is), CRH, 0)
            call logging_debuglog('Add boundary rhobksp'//str(is))
            call add_boundary_charge(rhobksp(:, :, :, :, 1, is), &
                                     1, 1, &
                                     sdoms(:, :, sdid(1) + 1), &
                                     bounds(:, :, sdid(1) + 1), &
                                     ctypes(:, :, 1, CRH), &
                                     fsizes(:, :, FRH), myid)
        end do

        if (require_output_data) then
            call logging_debuglog('Reduce rhodg')
            if (currmode /= OHH_PRIMARY_MODE) then
                call oh3_reduce_field(rhodg(1, 0, 0, 0, 1), rhodg(1, 0, 0, 0, 2), FRD)
            end if
            call logging_debuglog('Exchange rhodg')
            call oh3_exchange_borders(rhodg(1, 0, 0, 0, 1), rhodg(1, 0, 0, 0, 2), CRD, 0)
            call logging_debuglog('Add boundary rhodg')
            call add_boundary_charge(rhodg(:, :, :, :, 1), &
                                     minspec*2 + 1, nspec*2, &
                                     sdoms(:, :, sdid(1) + 1), &
                                     bounds(:, :, sdid(1) + 1), &
                                     ctypes(:, :, 1, CRD), &
                                     fsizes(:, :, FRD), myid)
        end if
    end subroutine
end module
