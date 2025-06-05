module m_conductivity
    use allcom
    use m_ohhinfo
    use particle_collision
    use finbound
    implicit none
    private

    public create_conductivity_array
    public calculate_conductive_current
    public update_conductive_rhobk

    logical :: has_initialized = .false.
    logical, allocatable :: is_conductive(:, :, :)

    double precision, allocatable :: condx(:, :, :)
    double precision, allocatable :: condy(:, :, :)
    double precision, allocatable :: condz(:, :, :)

contains

    subroutine update_conductive_rhobk(dt)
        double precision, intent(in) :: dt

        integer :: lnx
        integer :: lny
        integer :: lnz

        double precision, allocatable :: current_x(:, :, :)
        double precision, allocatable :: current_y(:, :, :)
        double precision, allocatable :: current_z(:, :, :)

        double precision, allocatable :: current(:, :, :)
        integer :: i, j, k

        if (conductivity == 0d0) then
            return
        end if

        if (.not. has_initialized) then
            call initialize_conductivities
        end if

        call ohhinfo_update(sdoms(:, :, sdid(1) + 1))

        lnx = xu - xl
        lny = yu - yl
        lnz = zu - zl

        allocate (current(0:lnx, 0:lny, 0:lnz))

        allocate (current_x(-1:lnx, 0:lny, 0:lnz))
        allocate (current_y(0:lnx, -1:lny, 0:lnz))
        allocate (current_z(0:lnx, 0:lny, -1:lnz))

        current_x(-1:lnx, 0:lny, 0:lnz) = condx(-1:lnx, 0:lny, 0:lnz)*eb(EX, -1:lnx, 0:lny, 0:lnz, 1)
        current_y(0:lnx, -1:lny, 0:lnz) = condy(0:lnx, -1:lny, 0:lnz)*eb(EY, 0:lnx, -1:lny, 0:lnz, 1)
        current_z(0:lnx, 0:lny, -1:lnz) = condz(0:lnx, 0:lny, -1:lnz)*eb(EZ, 0:lnx, 0:lny, -1:lnz, 1)

        current(0:lnx, 0:lny, 0:lnz) = 0d0

        current(0:lnx, 0:lny, 0:lnz) = current(0:lnx, 0:lny, 0:lnz) + current_x(-1:lnx - 1, :, :)
        current(0:lnx, 0:lny, 0:lnz) = current(0:lnx, 0:lny, 0:lnz) - current_x(0:lnx, :, :)

        current(0:lnx, 0:lny, 0:lnz) = current(0:lnx, 0:lny, 0:lnz) + current_y(:, -1:lny - 1, :)
        current(0:lnx, 0:lny, 0:lnz) = current(0:lnx, 0:lny, 0:lnz) - current_y(:, 0:lny, :)

        current(0:lnx, 0:lny, 0:lnz) = current(0:lnx, 0:lny, 0:lnz) + current_z(:, :, -1:lnz - 1)
        current(0:lnx, 0:lny, 0:lnz) = current(0:lnx, 0:lny, 0:lnz) - current_z(:, :, 0:lnz)

        rhobk(1, 0:lnx, 0:lny, 0:lnz, 3) = rhobk(1, 0:lnx, 0:lny, 0:lnz, 3) + dt*current(:, :, :)
    end subroutine

    subroutine initialize_conductivities
        integer :: i, j, k
        integer :: lnx, lny, lnz

        type(t_BoundaryList) :: boundaries

        call ohhinfo_update(sdoms(:, :, sdid(1) + 1))

        lnx = xu - xl
        lny = yu - yl
        lnz = zu - zl

        allocate (is_conductive(-1:lnx + 1, -1:lny + 1, -1:lnz + 1))
        allocate (condx(-1:lnx, 0:lny, 0:lnz))
        allocate (condy(0:lnx, -1:lny, 0:lnz))
        allocate (condz(0:lnx, 0:lny, -1:lnz))

        boundaries = create_simple_collision_boundaries(sdoms(1:2, 1:3, sdid(1) + 1))
        do concurrent(i=-1:lnx + 1, j=-1:lny + 1, k=-1:lnz + 1)
            ! block
            !     double precision :: grid_box(2, 3)
            !     double precision :: extent(2, 3)

            !     extent(:, :) = reshape([[0d0, 0d0], [0d0, 0d0], [0d0, 0d0]], [2, 3])

            !     grid_box(:, :) = reshape([[xl+i-0.5d0, xl+i+0.5d0], [yl+j-0.5d0, yl+j+0.5d0], [zl+k-0.5d0, zl+k+0.5d0]], shape=[2, 3])
            !     is_conductive(i, j, k) = boundaries%is_overlap(grid_box, extent=extent)
            ! end block

            ! Determine if a boundary exists within the grid box
            ! by whether the diagonal of each grid box intersects the internal boundary. 
            block
                double precision :: grid_box(2, 3)
                double precision :: p1(3), p2(3)
                type(t_CollisionRecord) :: record
                grid_box(:, :) = reshape([[xl + i - 0.5d0, xl + i + 0.5d0], &
                                          [yl + j - 0.5d0, yl + j + 0.5d0], &
                                          [zl + k - 0.5d0, zl + k + 0.5d0]], &
                                         shape=[2, 3])

                p1(:) = [grid_box(1, 1), grid_box(1, 2), grid_box(1, 3)]
                p2(:) = [grid_box(2, 1), grid_box(2, 2), grid_box(2, 3)]
                record = boundaries%check_collision(p1, p2)
                if (record%is_collided) then
                    is_conductive(i, j, k) = .true.
                    cycle
                end if

                p1(:) = [grid_box(2, 1), grid_box(1, 2), grid_box(1, 3)]
                p2(:) = [grid_box(1, 1), grid_box(2, 2), grid_box(2, 3)]
                record = boundaries%check_collision(p1, p2)
                if (record%is_collided) then
                    is_conductive(i, j, k) = .true.
                    cycle
                end if

                p1(:) = [grid_box(1, 1), grid_box(2, 2), grid_box(1, 3)]
                p2(:) = [grid_box(2, 1), grid_box(1, 2), grid_box(2, 3)]
                record = boundaries%check_collision(p1, p2)
                if (record%is_collided) then
                    is_conductive(i, j, k) = .true.
                    cycle
                end if

                p1(:) = [grid_box(2, 1), grid_box(2, 2), grid_box(1, 3)]
                p2(:) = [grid_box(1, 1), grid_box(1, 2), grid_box(2, 3)]
                record = boundaries%check_collision(p1, p2)
                if (record%is_collided) then
                    is_conductive(i, j, k) = .true.
                    cycle
                end if

                is_conductive(i, j, k) = .false.
            end block
        end do

        do concurrent(i=-1:lnx, j=0:lny, k=0:lnz)
            if (is_conductive(i, j, k) .and. is_conductive(i + 1, j, k)) then
                condx(i, j, k) = conductivity
            else
                condx(i, j, k) = 0d0
            end if
        end do

        do concurrent(i=0:lnx, j=-1:lny, k=0:lnz)
            if (is_conductive(i, j, k) .and. is_conductive(i, j + 1, k)) then
                condy(i, j, k) = conductivity
            else
                condy(i, j, k) = 0d0
            end if
        end do

        do concurrent(i=0:lnx, j=0:lny, k=-1:lnz)
            if (is_conductive(i, j, k) .and. is_conductive(i, j, k + 1)) then
                condz(i, j, k) = conductivity
            else
                condz(i, j, k) = 0d0
            end if
        end do

        has_initialized = .true.
    end subroutine

    subroutine create_conductivity_array( &
        nx, ny, nz, is_conductive, condx, condy, condz, val)
        integer, intent(in) :: nx
        integer, intent(in) :: ny
        integer, intent(in) :: nz
        logical, intent(in) :: is_conductive(0:nx, 0:ny, 0:nz)
        double precision, intent(out) :: condx(0:nx - 1, 0:ny, 0:nz)
        double precision, intent(out) :: condy(0:nx, 0:ny - 1, 0:nz)
        double precision, intent(out) :: condz(0:nx, 0:ny, 0:nz - 1)
        double precision, intent(in) :: val

        condx(:, :, :) = 0d0
        condy(:, :, :) = 0d0
        condz(:, :, :) = 0d0

        where (is_conductive(0:nx - 1, :, :) .and. is_conductive(1:nx, :, :))
            condx = val
        end where

        where (is_conductive(:, 0:ny - 1, :) .and. is_conductive(:, 1:ny, :))
            condy = val
        end where

        where (is_conductive(:, :, 0:nz - 1) .and. is_conductive(:, :, 1:nz))
            condz = val
        end where
    end subroutine

    subroutine calculate_conductive_current(nx, ny, nz, current, ex, ey, ez, condx, condy, condz)
        integer, intent(in) :: nx
        integer, intent(in) :: ny
        integer, intent(in) :: nz
        double precision, intent(out) :: current(0:nx, 0:ny, 0:nz)
        double precision, intent(in) :: ex(0:nx - 1, 0:ny, 0:nz)
        double precision, intent(in) :: ey(0:nx, 0:ny - 1, 0:nz)
        double precision, intent(in) :: ez(0:nx, 0:ny, 0:nz - 1)
        double precision, intent(in) :: condx(0:nx - 1, 0:ny, 0:nz)
        double precision, intent(in) :: condy(0:nx, 0:ny - 1, 0:nz)
        double precision, intent(in) :: condz(0:nx, 0:ny, 0:nz - 1)

        double precision, allocatable :: current_x(:, :, :)
        double precision, allocatable :: current_y(:, :, :)
        double precision, allocatable :: current_z(:, :, :)

        allocate (current_x(0:nx - 1, 0:ny, 0:nz))
        allocate (current_y(0:nx, 0:ny - 1, 0:nz))
        allocate (current_z(0:nx, 0:ny, 0:nz - 1))

        current_x(:, :, :) = condx(:, :, :)*ex(:, :, :)
        current_y(:, :, :) = condy(:, :, :)*ey(:, :, :)
        current_z(:, :, :) = condz(:, :, :)*ez(:, :, :)

        current(:, :, :) = 0d0

        current(1:nx, 0:ny, 0:nz) = current(1:nx, 0:ny, 0:nz) + current_x(:, :, :)
        current(0:nx - 1, 0:ny, 0:nz) = current(0:nx - 1, 0:ny, 0:nz) - current_x(:, :, :)

        current(0:nx, 1:ny, 0:nz) = current(0:nx, 1:ny, 0:nz) + current_y(:, :, :)
        current(0:nx, 0:ny - 1, 0:nz) = current(0:nx, 0:ny - 1, 0:nz) - current_y(:, :, :)

        current(0:nx, 0:ny, 1:nz) = current(0:nx, 0:ny, 1:nz) + current_z(:, :, :)
        current(0:nx, 0:ny, 0:nz - 1) = current(0:nx, 0:ny, 0:nz - 1) - current_z(:, :, :)

        current(:, :, :) = current(:, :, :)
    end

end module
