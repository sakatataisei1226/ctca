module m_esses_worker
    use allcom
    use m_ctcamain
    use m_esses

#include "oh_stats.h"
#define MCW local_comm

    implicit none

    private
    public esses_worker_run

contains

    !> Run electric static simulation as worker.
    subroutine esses_worker_run
        integer(kind=4) :: ustep
        integer(kind=4) :: mpierr
        integer(kind=4) :: nprocs,myrank
        integer(kind=8) :: mnmtotal = 0
        integer(kind=8) :: mnmlocal = 0
        real(kind=8) :: smltime, emltime
        real(kind=8) :: performance

        call CTCAW_init(0,1)
        local_comm=CTCA_subcomm
        call MPI_Comm_size(MPI_COMM_WORLD,nprocs,mpierr)
        call MPI_Comm_rank(MPI_COMM_WORLD,myrank,mpierr)
        call MPI_Comm_size(MCW, nnode, mpierr)
        call MPI_Comm_rank(MCW, myid, mpierr)
        call oh1_fam_comm(CTCA_subcomm)

        call hdfinit()

        emflag = 0
        call input

        call esses_initialize

        print *, "es_wrk check 1"

        call cotocoa_init

        print *, "es_wrk check 2"
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

            print *, "es_wrk check 1"

            call esses_mainstep(ustep)
            print *, "es_wrk check 2"
            call cotocoa_mainstep
            print *, "es_wrk check 3"

            mnmlocal = mnmlocal + sum(totalp(:, :))
        end do

        call MPI_Barrier(MCW, mpierr)
        emltime = MPI_Wtime()

        call esses_finalize
        call cotocoa_finalize

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
        call CTCAW_finalize()
    end subroutine
end module