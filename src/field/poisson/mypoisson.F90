module m_mypoisson
    use m_poisson_solver, only: t_PoissonSolver3d, new_PoissonSolver3d
    use mpifft
    use allcom, only: sdoms, sdid, myid, nnode, &
                      nx, ny, nz, &
                      mtd_vbnd, &
                      boundary_conditions, &
                      rho, poi, phi, &
                      CPH
    use mpi
    implicit none

    type(t_MPIFFTW3_Factory) :: fftw3_factory
    type(t_MPI_FFTExecutor3d), target :: fft3d
    type(t_PoissonSolver3d) :: poisson_solver
    logical :: initialized = .false.

    private
    public solve_poisson

contains

    subroutine init_poisson
        type(t_Block) :: local_block
        type(t_Block) :: global_block

        integer :: local_start(3), local_end(3)
        integer :: global_start(3), global_end(3)
        integer :: boundary_types(3)
        double precision :: boundary_values(2, 3)

        class(t_MPI_FFTExecutor3d), pointer :: pfft3d

        integer :: xl, yl, zl, xu, yu, zu
        integer :: axis
        integer :: ngrids(3)
        integer :: ps

        ps = 1
        xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
        yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
        zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)

        local_block = new_Block([xl, yl, zl], [xu, yu, zu])

        ngrids(:) = [nx, ny, nz]
        do axis = 1, 3
            select case (mtd_vbnd(axis))
            case (0)
                global_start(axis) = 0
                global_end(axis) = ngrids(axis) - 1
                boundary_types(axis) = BoundaryType_Periodic

            case (1)
                global_start(axis) = 1
                global_end(axis) = ngrids(axis) - 1
                boundary_types(axis) = BoundaryType_Dirichlet

            case (2)
                global_start(axis) = 0
                global_end(axis) = ngrids(axis)
                boundary_types(axis) = BoundaryType_Neumann
            end select
        end do
        global_block = new_Block(global_start, global_end)

        boundary_values(:, :) = boundary_conditions(:, :)
        ! boundary_types(:) = [BoundaryType_Dirichlet, BoundaryType_Dirichlet, BoundaryType_Dirichlet]

        fftw3_factory = new_MPIFFTW3_Factory(local_block, global_block, myid, nnode, MPI_COMM_WORLD)
        fft3d = fftw3_factory%create(boundary_types)

        pfft3d => fft3d
        poisson_solver = new_PoissonSolver3d(local_block, global_block, pfft3d, 1.0d0, boundary_values=boundary_values)
    end subroutine

    subroutine solve_poisson
        integer :: ps
        integer :: xl, yl, zl, xu, yu, zu

        if (.not. initialized) then
            call init_poisson
            initialized = .true.
        end if

        ps = 1
        xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
        yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
        zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)

        ! print *, "start poisson solve"
        call poisson_solver%solve(rho(1, 0:xu - xl, 0:yu - yl, 0:zu - zl, 1), &
                                  phi(1, 0:xu - xl, 0:yu - yl, 0:zu - zl, 1))
        ! print *, "end poisson solve"

        ! print *, "start exchange borders in poisson"
        call oh3_exchange_borders(phi(1, 0, 0, 0, 1), phi(1, 0, 0, 0, 2), CPH, 0)
        ! print *, "end exchange borders in poisson"
    end subroutine

end module
