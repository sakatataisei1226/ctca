module m_logging
    use allcom
    use mpi
    implicit none

    logical :: logging_debug = .false.
    logical :: logging_barrier = .false.
    character(:), allocatable :: log_filename

    private
    public logging_initialize
    public logging_debuglog

contains

    subroutine get_env(name, value, success)
        character(:), intent(in), allocatable :: name
        character(:), intent(out), allocatable :: value
        logical, intent(out) :: success
        integer :: length
        integer :: status

        success = .false.

        call get_environment_variable(name, status=status, length=length)

        if (status == 1) then
            if (myid == 0) print *, 'Environment variable "', name, '" does not exist.'
        else if (status /= 0) then
            if (myid == 0) print *, 'Error', status, 'for environment variable "', name, '"'
        else
            allocate (character(length) :: value)
            call get_environment_variable(name, value=value)
            success = .true.
        end if
    end subroutine

    subroutine logging_initialize
        character(:), allocatable :: name
        character(:), allocatable :: debug
        character(:), allocatable :: isbarrier
        character(:), allocatable :: log_filename_
        logical :: success

        integer :: idf

        log_filename = '__none__'

        name = 'EMSES_DEBUG'
        call get_env(name, debug, success)

        if (.not. success) then
            return
        end if

        logging_debug = (debug == 'yes')

        name = 'EMSES_DEBUG_BARRIER'
        call get_env(name, isbarrier, success)

        if (success) then
            logging_barrier = (isbarrier == 'yes')
        end if

        name = 'EMSES_LOGFILE'
        call get_env(name, log_filename_, success)

        if (.not. success) then
            return
        end if

        log_filename = log_filename_
        open (newunit=idf, file=log_filename, action='write', &
              status='replace')
        close (idf)
    end subroutine

    subroutine logging_debuglog(message, barrier)
        character(len=*), intent(in) :: message
        logical, intent(in), optional :: barrier

        integer :: ierr
        integer :: idf

        if (logging_debug) then
            if ((present(barrier) .and. barrier) .or. logging_barrier) then
                call MPI_Barrier(MPI_COMM_WORLD, ierr)
            end if

            if (log_filename == '__none__') then
                print *, '[myid=', myid, ']', message
            else
                open (newunit=idf, file=log_filename, action='write', &
                      status='old', position='append')
                write (idf, *) '[myid=', myid, ']', message
                close (idf)
            end if
        end if
    end subroutine

end module
