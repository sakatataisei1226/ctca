!> Module for surface charging accelaration
!!
!! Reduce the charging time by applying a factor to the EMA of the current flowing into the surface.
!!
!! Usage:
!!   use m_grad_ema
!!
!!   call grad_ema_init
!!
!!   do mainloop
!!     ! rhobk(:, 3) = rhobk(:, 3) + rhobk(:, 1)
!!     call grad_ema_update_rhobk
!!   end do
!!
module m_grad_ema
    use oh_type
    use paramt
    use allcom
    implicit none

    double precision, allocatable    :: rhobk_ema(:, :, :, :)
    double precision, allocatable    :: gcs(:, :, :, :)
    integer, parameter :: max_coefs = 10
    character(len=30) :: grad_type = "constant" ! "linear", "bezier"
    double precision :: coef_start(max_coefs) = 1.0d0
    double precision :: smooth_coef(max_coefs) = 1.0d0
    double precision :: grad_coef(max_coefs) = 1.0d0
    double precision :: grad_coef2(max_coefs) = 1.0d0
    double precision :: control_grad_coef(2, max_coefs) = 1.0d0

    double precision, allocatable :: moments(:, :, :, :)
    double precision, allocatable :: variances(:, :, :, :)

    namelist /gradema/ coef_start, smooth_coef, grad_coef, grad_coef2, grad_type, control_grad_coef

    private
    public grad_ema_init
    public grad_ema_update_rhobk

contains

    !> Initialize this module.
    subroutine grad_ema_init
        character(len=20) :: inputfilename
        integer :: iosta

        call get_command_argument(1, inputfilename)
        open (1, file=inputfilename)

        rewind (1); read (1, nml=gradema, IOSTAT=iosta)
        if (myid .eq. 0 .and. iosta .eq. -1) then
            print *, "Warning.Input: nml=gradema not found"
        end if

        close (1)

        allocate (rhobk_ema(1, &
                            fsizes(1, 1, FRH):fsizes(2, 1, FRH), &
                            fsizes(1, 2, FRH):fsizes(2, 2, FRH), &
                            fsizes(1, 3, FRH):fsizes(2, 3, FRH)))
        rhobk_ema(:, :, :, :) = 0.0d0

        allocate (gcs, mold=rhobk_ema)

        if (myid == 0) then
            print *, "...initialized grad_ema", grad_coef, smooth_coef
        end if
    end subroutine

    !> Add rhobk by gradient accelaration and EMA.
    !!
    !! rhobk(t) = rhobk(t-1) + b * rhobk_ema(t)
    !! rhobk_ema(t) = (1 - a) * rhobk_ema(t-1) + a * drhobk(t)
    !! a = smooth_coef
    !! b = grad_coef
    subroutine grad_ema_update_rhobk
        double precision :: sc
        double precision :: gc
        double precision :: gc2
        double precision :: control_height
        double precision :: control_gc

        integer :: is

        sc = 1.0d0
        gc = 1.0d0
        gc2 = 1.0d0
        control_height = zssurf
        control_gc = 1.0d0

        block
            integer :: i
            double precision :: step_rate

            step_rate = (1.0d0*istep)/nstep
            do i = 1, max_coefs
                if (step_rate < coef_start(i)) then
                    sc = smooth_coef(i)
                    gc = grad_coef(i)
                    gc2 = grad_coef2(i)
                    control_height = control_grad_coef(1, i)
                    control_gc = control_grad_coef(2, i)
                    exit
                end if
            end do
        end block

        rhobk_ema(:, :, :, :) = (1 - sc)*rhobk_ema(:, :, :, :) &
                                + sc*rhobk(:, :, :, :, 1)

        if (grad_type == "constant") then

            rhobk(:, :, :, :, 3) = rhobk(:, :, :, :, 3) + gc*rhobk_ema(:, :, :, :)
            rhobksp(:, :, :, :, 3, :) = rhobksp(:, :, :, :, 3, :) + gc*rhobksp(:, :, :, :, 1, :)

        else if (grad_type == "linear") then

            block
                integer :: z
                double precision :: a
                integer :: zl, zu

                zl = sdoms(1, 3, sdid(1) + 1)
                zu = sdoms(2, 3, sdid(1) + 1)

                a = (gc - gc2)/(zssurf - zlrechole(2))
                do z = zl - 1, zu + 1
                    gcs(1, :, :, z - zl) = a*(z - zlrechole(2)) + gc2
                end do
            end block

            rhobk(:, :, :, :, 3) = rhobk(:, :, :, :, 3) + gcs(:, :, :, :)*rhobk_ema(:, :, :, :)

        else if (grad_type == "bezier") then

            block
                double precision :: t
                double precision :: a, b, c, d
                integer :: zl, zu
                integer :: z

                zl = sdoms(1, 3, sdid(1) + 1)
                zu = sdoms(2, 3, sdid(1) + 1)

                a = zlrechole(2) - 2.0d0*control_height + zssurf
                b = control_height - zlrechole(2)
                c = zlrechole(2)
                do z = zl - 1, zu + 1
                    if (z < zlrechole(2) .or. zssurf < z) then
                        gcs(1, :, :, z - zl) = 1.0d0
                        cycle
                    end if

                    if (abs(a) < 1d-8) then
                        t = (z - c) / (2.0d0 * b)
                    else
                        d = b*b - a*(c - z)
                        t = (-b + dsqrt(d)) / a
                    end if

                    gcs(1, :, :, z-zl) = (gc2 - 2.0d0*control_gc + gc)*t*t + 2.0d0*(control_gc - gc2)*t + gc2
                end do
            end block

            rhobk(:, :, :, :, 3) = rhobk(:, :, :, :, 3) + gcs(:, :, :, :)*rhobk_ema(:, :, :, :)

        end if
    end subroutine

end module
