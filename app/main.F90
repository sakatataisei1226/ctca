#ifndef __VERSION__ 
#define __VERSION__ '0.0.0'
#endif
program main
#ifdef COTOCOA_DISABLED
    use m_esses, only: esses_run
#endif
#ifdef COTOCOA_REQUESTER
    use m_esses_requester
#endif
#ifdef COTOCOA_WORKER
    use m_esses_worker
#endif

    implicit none

    character(len=20) :: inputfilename
    integer(kind=4) :: emflag = 1

    namelist /esorem/ emflag

    call get_command_argument(1, inputfilename)

    if (inputfilename == '--version') then
        print *, __VERSION__
        stop
    endif

    open (1, file=inputfilename)
    read (1, esorem)
    close (1)

    if (emflag .eq. 0) then
#ifdef COTOCOA_DISABLED
        call esses_run
#endif
#ifdef COTOCOA_REQUESTER
        call esses_requester_run
#endif
#ifdef COTOCOA_WORKER
        call esses_worker_run
#endif
    elseif (emflag .eq. 1) then
        call emses
    else if (emflag .eq. 2) then
        call imses
    end if

end program
