module m_steplog
    use m_str
    use mpi
    implicit none

    type t_StepLogger
        character, allocatable :: filename(:)
        double precision :: interval_sec
        double precision, private :: prev_sec
    contains
        procedure update => stepLogger_update
    end type

    private
    public t_StepLogger, new_StepLogger

contains

    function new_StepLogger(filename, interval_sec) result(obj)
        character(len=*), intent(in), optional :: filename
        double precision, intent(in), optional :: interval_sec
        type(t_StepLogger) :: obj

        if (present(filename)) then
            obj%filename = filename
        else
            obj%filename = 'current.step'
        end if

        if (present(interval_sec)) then
            obj%interval_sec = interval_sec
        else
            obj%interval_sec = 10
        end if
        obj%prev_sec = -1e10
    end function

    subroutine stepLogger_update(self, current_step, max_step)
        class(t_StepLogger) :: self
        integer, intent(in) :: current_step
        integer, intent(in) :: max_step

        !> output identifier
        integer :: idf
        double precision :: cur_sec, elapsed_sec

        cur_sec = MPI_Wtime()
        elapsed_sec = cur_sec - self%prev_sec
        if (elapsed_sec <= self%interval_sec) then
            return
        end if

        self%prev_sec = cur_sec

        open (newunit=idf, file=self%filename, status='replace')

        write (idf, *) str(current_step)//' / '//str(max_step)

        close (idf)
    end subroutine

end module
