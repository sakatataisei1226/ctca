module m_output
    use mpi
    use m_h5fortran
    implicit none

contains

    subroutine output(myid, nodes, sdom)
        integer, intent(in) :: myid
        integer, intent(in) :: nodes(3)
        integer, intent(in) :: sdom(2, 3)

        integer(kind=HSIZE_T) :: local_data_shape(3)
        integer(kind=HSIZE_T) :: offset(3)
        integer(kind=HSIZE_T) :: all_data_shape(3)
        double precision :: local_data(1, 1, 2)

        integer(kind=HSIZE_T) lb(3)
        integer(kind=HSIZE_T) ub(3)

        block
            integer :: coord(3)
            integer :: i

            integer :: nxsd, nysd, nzsd
            integer :: dsize(3)

            nxsd = sdom(2, 1) - sdom(1, 1)
            nysd = sdom(2, 2) - sdom(1, 2)
            nzsd = sdom(2, 3) - sdom(1, 3)
            dsize(1) = nxsd + 1
            dsize(2) = nysd + 1
            dsize(3) = nzsd + 1

            coord(1) = mod(myid, nodes(1))
            coord(2) = mod(myid/nodes(1), nodes(2))
            coord(3) = myid/(nodes(1)*nodes(2))

            do i = 1, 3
                if (coord(i) .eq. 0 .and. coord(i) .eq. nodes(i) - 1) then
                    lb(i) = 0
                    ub(i) = dsize(i) - 1
                    offset(i) = 0
                else if (coord(i) .eq. 0) then
                    lb(i) = 0
                    ub(i) = dsize(i) - 2
                    offset(i) = 0
                else if (coord(i) .eq. nodes(i) - 1) then
                    lb(i) = 0
                    ub(i) = dsize(i) - 1
                    offset(i) = (dsize(i) - 1)*coord(i)
                else
                    lb(i) = 0
                    ub(i) = dsize(i) - 2
                    offset(i) = (dsize(i) - 1)*coord(i)
                end if
                local_data_shape(i) = ub(i) - lb(i) + 1
            end do
        end block

        call write_field3d('ex00_0000.h5', 'ex', '0000', &
                           local_data_shape, local_data, &
                           offset, all_data_shape)
    end subroutine

    subroutine write_field3d(filename, group_name, data_name, &
                             local_data_shape, local_data, &
                             offset, all_data_shape)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: group_name
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: local_data_shape(3)
        double precision, intent(in) :: local_data(local_data_shape(1), &
                                                   local_data_shape(2), &
                                                   local_data_shape(3))
        integer(kind=HSIZE_T), intent(in) :: offset(3)
        integer(kind=HSIZE_T), intent(in) :: all_data_shape(3)

        type(t_h5file) :: h5f
        type(t_h5group) :: group
        integer :: ierr

        call h5f%open(filename, 'w', MPI_COMM_WORLD)
        group = h5f%get_group(group_name, MPI_COMM_WORLD)

        call group%create_dataset_double_3d_parallel(data_name, &
                                                     local_data_shape, local_data, &
                                                     offset, all_data_shape, &
                                                     ierr)

        call group%close()
        call h5f%close()
    end subroutine

end module
