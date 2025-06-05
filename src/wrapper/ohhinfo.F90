module m_ohhinfo
    implicit none

    !> Helpand-helper configuration in (re)build and in secondary mode (= -1).
    integer, parameter :: OHH_REQUIRES_BCAST = -1
    !> In primary mode (= 0).
    integer, parameter :: OHH_PRIMARY_MODE = 0
    !> In secondary mode (= 1).
    integer, parameter :: OHH_SECONDARY_MODE = 1

    
    ! Region of subdomain (integer ver.)
    integer :: xl, xu, yl, yu, zl, zu
    ! Region of subdomain (double precision ver.)
    double precision :: dxl, dxu, dyl, dyu, dzl, dzu

    !> Subdomain's width (integer ver.)
    integer :: ngx, ngy, ngz
    !> Subdomain's width (double precision ver.)
    double precision :: dngx, dngy, dngz

contains

    subroutine ohhinfo_update(sdom)
        !> Subdomain range. (i.e. sdoms(:, :, sdid(ps) + 1))
        integer, intent(in) :: sdom(:, :)

        xl = sdom(1, 1)
        xu = sdom(2, 1)
        yl = sdom(1, 2)
        yu = sdom(2, 2)
        zl = sdom(1, 3)
        zu = sdom(2, 3)
    
        dxl = xl
        dxu = xu
        dyl = yl
        dyu = yu
        dzl = zl
        dzu = zu

        ngx = xu - xl
        ngy = yu - yl
        ngz = zu - zl
        dngx = ngx
        dngy = ngy
        dngz = ngz
    end subroutine

end module