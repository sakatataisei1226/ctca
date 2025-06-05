!> Wrapper module for HDF5 library.
!>
!> example:
!> program example
!>     use mpi
!>     use m_h5fortran
!>     implicit none
!>
!>     type(t_h5file) :: h5f
!>     type(t_h5group) :: group
!>
!>     integer :: size, id, ierr
!>     INTEGER(kind=HSIZE_T) :: data_shape(1)
!>
!>     integer :: d(3), d2(3)
!>
!>     call MPI_Init(ierr)
!>     call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
!>     call MPI_Comm_rank(MPI_COMM_WORLD, id, ierr)
!>
!>     call h5fortran_init
!>
!>     call h5f%open('example.h5', 'w', MPI_COMM_WORLD)
!>
!>     group = h5f%create_group('a', MPI_COMM_WORLD)
!>     data_shape = [3]
!>     d = [0, id, id*3]
!>     call group%create_dataset_integer_1d('sample', data_shape, d, ierr)
!>     call group%close()
!>
!>     call h5f%close()
!>
!>     call h5f%open('example.h5', 'r', MPI_COMM_WORLD)
!>     group = h5f%get_group('a', MPI_COMM_WORLD)
!>     data_shape = [3]
!>     call group%get_dataset_integer_1d('sample', data_shape, d2, ierr)
!>     print *, d2
!>
!>     call group%close()
!>     call h5f%close()
!>
!>     call h5fortran_finalize
!>
!>     call MPI_Finalize(ierr)
!> end program
module m_h5fortran
    use HDF5
    use mpi
    implicit none

    integer, parameter :: MAX_NAME_LENGTH = 100
    integer, parameter :: MAX_GROUPS = 10

    type :: t_h5node
        integer(kind=HID_T) :: id
    contains
        procedure :: create_group => h5node_create_group
        procedure :: get_group => h5node_get_group

        procedure :: create_dataset_integer_1d => h5node_create_dataset_integer_1d
        procedure :: create_dataset_integer_2d => h5node_create_dataset_integer_2d
        procedure :: create_dataset_integer_3d => h5node_create_dataset_integer_3d
        procedure :: create_dataset_real_1d => h5node_create_dataset_real_1d
        procedure :: create_dataset_real_2d => h5node_create_dataset_real_2d
        procedure :: create_dataset_real_3d => h5node_create_dataset_real_3d
        procedure :: create_dataset_double_1d => h5node_create_dataset_double_1d
        procedure :: create_dataset_double_2d => h5node_create_dataset_double_2d
        procedure :: create_dataset_double_3d => h5node_create_dataset_double_3d
        procedure :: create_dataset_double_3d_parallel => h5node_create_dataset_double_3d_parallel

        generic :: create_dataset => &
            create_dataset_integer_1d, &
            create_dataset_integer_2d, &
            create_dataset_integer_3d, &
            create_dataset_real_1d, &
            create_dataset_real_2d, &
            create_dataset_real_3d, &
            create_dataset_double_1d, &
            create_dataset_double_2d, &
            create_dataset_double_3d, &
            create_dataset_double_3d_parallel

        procedure :: get_dataset_integer_1d => h5node_get_dataset_integer_1d
        procedure :: get_dataset_integer_2d => h5node_get_dataset_integer_2d
        procedure :: get_dataset_integer_3d => h5node_get_dataset_integer_3d
        procedure :: get_dataset_real_1d => h5node_get_dataset_real_1d
        procedure :: get_dataset_real_2d => h5node_get_dataset_real_2d
        procedure :: get_dataset_real_3d => h5node_get_dataset_real_3d
        procedure :: get_dataset_double_1d => h5node_get_dataset_double_1d
        procedure :: get_dataset_double_2d => h5node_get_dataset_double_2d
        procedure :: get_dataset_double_3d => h5node_get_dataset_double_3d

        generic :: get_dataset => &
            get_dataset_integer_1d, &
            get_dataset_integer_2d, &
            get_dataset_integer_3d, &
            get_dataset_real_1d, &
            get_dataset_real_2d, &
            get_dataset_real_3d, &
            get_dataset_double_1d, &
            get_dataset_double_2d, &
            get_dataset_double_3d
    end type

    type, extends(t_h5node) :: t_h5file
        character :: filename(MAX_NAME_LENGTH)
    contains
        procedure :: open => h5file_open
        procedure :: close => h5file_close
    end type

    type, extends(t_h5node) :: t_h5group
        character :: name(MAX_NAME_LENGTH)
    contains
        procedure :: close => h5group_close
    end type

contains

    subroutine h5fortran_init(force)
        logical, intent(in), optional :: force
        integer :: hdferr
        logical :: check_error

        check_error = .not. (present(force) .and. force)

        call h5open_f(hdferr)
        if (check_error .and. hdferr /= 0) then
            write (0, *) '[Error@h5fortran_init] Cannot initialize HDF5.'
            error stop
        end if

        call h5eset_auto_f(0, hdferr)
        if (check_error .and. hdferr /= 0) then
            write (0, *) '[Error@h5fortran_init] Cannot set auto_f=0.'
            error stop 1
        end if
    end subroutine

    subroutine h5fortran_finalize(force)
        logical, intent(in), optional :: force
        integer :: hdferr
        logical :: check_error

        check_error = .not. (present(force) .and. force)

        call h5close_f(hdferr)
        if (check_error .and. hdferr /= 0) then
            write (0, *) '[Error@h5fortran_finalize] Cannot finalize HDF5.'
            error stop 1
        end if
    end subroutine

    function h5node_create_group(self, group_name, comm) result(group)
        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: group_name
        integer, intent(in) :: comm
        type(t_h5group) :: group

        !> Property list identifier.
        integer(kind=HID_T):: prp_id
        integer :: hdferr

        group%name = group_name

        call h5pcreate_f(H5P_FILE_ACCESS_F, prp_id, hdferr)
        call h5pset_fapl_mpio_f(prp_id, comm, MPI_INFO_NULL, hdferr)

        call h5gcreate_f(self%id, group_name, group%id, hdferr)

        call h5pclose_f(prp_id, hdferr)
    end function

    function h5node_get_group(self, group_name, comm) result(group)
        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: group_name
        integer, intent(in) :: comm
        type(t_h5group) :: group

        !> Property list identifier.
        integer(kind=HID_T):: prp_id
        integer :: hdferr

        group%name = group_name

        call h5pcreate_f(H5P_FILE_ACCESS_F, prp_id, hdferr)
        call h5pset_fapl_mpio_f(prp_id, comm, MPI_INFO_NULL, hdferr)

        call h5gopen_f(self%id, group_name, group%id, hdferr)

        call h5pclose_f(prp_id, hdferr)
    end function

    subroutine h5file_open(self, filename, action, comm)
        class(t_h5file), intent(inout) :: self
        character(len=*), intent(in) :: filename
        !> Action char (Read: 'r', Write: 'w', Read-Write: 'rw', create: 'c')
        character(len=*), intent(in) :: action
        !> MPI communicator.
        integer, intent(in) :: comm

        !> Property list identifier.
        integer(kind=HID_T):: prp_id
        integer :: hdferr

        self%filename = filename

        call h5pcreate_f(H5P_FILE_ACCESS_F, prp_id, hdferr)
        call h5pset_fapl_mpio_f(prp_id, comm, MPI_INFO_NULL, hdferr)

        if (action == 'r') then
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, self%id, hdferr, access_prp=prp_id)
        else if (action == 'c') then
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, self%id, hdferr, access_prp=prp_id)
        else if (action == 'w' .or. action == 'rw' .or. action == 'wr') then
            call h5fopen_f(filename, H5F_ACC_RDWR_F, self%id, hdferr, access_prp=prp_id)
            if (hdferr /= 0) then
                call h5fcreate_f(filename, H5F_ACC_TRUNC_F, self%id, hdferr, access_prp=prp_id)
            end if
        else
            write (0, *) '[Error@h5file_open] Action "'//action//'" is invalid.'
            error stop 1
        end if

        call h5pclose_f(prp_id, hdferr)
    end subroutine

    subroutine h5file_close(self)
        class(t_h5file), intent(inout) :: self
        integer :: ierr

        call h5fclose_f(self%id, ierr)
    end subroutine

    subroutine h5group_close(self)
        class(t_h5group), intent(inout) :: self
        integer :: ierr

        call h5gclose_f(self%id, ierr)
    end subroutine

    subroutine h5node_create_dataset_integer_1d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 1

        class(t_h5node), intent(inout) ::self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        integer, intent(in) :: data(data_shape(1))
        integer, intent(out) :: ierr

        integer(kind=HID_T) data_id, space_id

        call h5screate_simple_f(rank, data_shape, space_id, ierr)
        if (ierr /= 0) return
        call h5dcreate_f(self%id, data_name, H5T_NATIVE_INTEGER, space_id, data_id, ierr)
        if (ierr /= 0) return
        call h5dwrite_f(data_id, H5T_NATIVE_INTEGER, data, data_shape, ierr)
        if (ierr /= 0) return
        call h5dclose_f(data_id, ierr)
        if (ierr /= 0) return
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_create_dataset_integer_2d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 2

        class(t_h5node), intent(inout) ::self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        integer, intent(in) :: data(data_shape(1), data_shape(2))
        integer, intent(out) :: ierr

        integer(kind=HID_T) data_id, space_id

        call h5screate_simple_f(rank, data_shape, space_id, ierr)
        if (ierr /= 0) return
        call h5dcreate_f(self%id, data_name, H5T_NATIVE_INTEGER, space_id, data_id, ierr)
        if (ierr /= 0) return
        call h5dwrite_f(data_id, H5T_NATIVE_INTEGER, data, data_shape, ierr)
        if (ierr /= 0) return
        call h5dclose_f(data_id, ierr)
        if (ierr /= 0) return
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_create_dataset_integer_3d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 3

        class(t_h5node), intent(inout) ::self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        integer, intent(in) :: data(data_shape(1), data_shape(2), data_shape(3))
        integer, intent(out) :: ierr

        integer(kind=HID_T) data_id, space_id

        call h5screate_simple_f(rank, data_shape, space_id, ierr)
        if (ierr /= 0) return
        call h5dcreate_f(self%id, data_name, H5T_NATIVE_INTEGER, space_id, data_id, ierr)
        if (ierr /= 0) return
        call h5dwrite_f(data_id, H5T_NATIVE_INTEGER, data, data_shape, ierr)
        if (ierr /= 0) return
        call h5dclose_f(data_id, ierr)
        if (ierr /= 0) return
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_create_dataset_real_1d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 1

        class(t_h5node), intent(inout) ::self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        real, intent(in) :: data(data_shape(1))
        integer, intent(out) :: ierr

        integer(kind=HID_T) data_id, space_id

        call h5screate_simple_f(rank, data_shape, space_id, ierr)
        if (ierr /= 0) return
        call h5dcreate_f(self%id, data_name, H5T_NATIVE_REAL, space_id, data_id, ierr)
        if (ierr /= 0) return
        call h5dwrite_f(data_id, H5T_NATIVE_REAL, data, data_shape, ierr)
        if (ierr /= 0) return
        call h5dclose_f(data_id, ierr)
        if (ierr /= 0) return
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_create_dataset_real_2d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 2

        class(t_h5node), intent(inout) ::self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        real, intent(in) :: data(data_shape(1), data_shape(2))
        integer, intent(out) :: ierr

        integer(kind=HID_T) data_id, space_id

        call h5screate_simple_f(rank, data_shape, space_id, ierr)
        if (ierr /= 0) return
        call h5dcreate_f(self%id, data_name, H5T_NATIVE_REAL, space_id, data_id, ierr)
        if (ierr /= 0) return
        call h5dwrite_f(data_id, H5T_NATIVE_REAL, data, data_shape, ierr)
        if (ierr /= 0) return
        call h5dclose_f(data_id, ierr)
        if (ierr /= 0) return
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_create_dataset_real_3d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 3

        class(t_h5node), intent(inout) ::self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        real, intent(in) :: data(data_shape(1), data_shape(2), data_shape(3))
        integer, intent(out) :: ierr

        integer(kind=HID_T) data_id, space_id

        call h5screate_simple_f(rank, data_shape, space_id, ierr)
        if (ierr /= 0) return
        call h5dcreate_f(self%id, data_name, H5T_NATIVE_REAL, space_id, data_id, ierr)
        if (ierr /= 0) return
        call h5dwrite_f(data_id, H5T_NATIVE_REAL, data, data_shape, ierr)
        if (ierr /= 0) return
        call h5dclose_f(data_id, ierr)
        if (ierr /= 0) return
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_create_dataset_double_1d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 1

        class(t_h5node), intent(inout) ::self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        double precision, intent(in) :: data(data_shape(1))
        integer, intent(out) :: ierr

        integer(kind=HID_T) data_id, space_id

        call h5screate_simple_f(rank, data_shape, space_id, ierr)
        if (ierr /= 0) return
        call h5dcreate_f(self%id, data_name, H5T_NATIVE_DOUBLE, space_id, data_id, ierr)
        if (ierr /= 0) return
        call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, data, data_shape, ierr)
        if (ierr /= 0) return
        call h5dclose_f(data_id, ierr)
        if (ierr /= 0) return
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_create_dataset_double_2d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 2

        class(t_h5node), intent(inout) ::self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        double precision, intent(in) :: data(data_shape(1), data_shape(2))
        integer, intent(out) :: ierr

        integer(kind=HID_T) data_id, space_id

        call h5screate_simple_f(rank, data_shape, space_id, ierr)
        if (ierr /= 0) return
        call h5dcreate_f(self%id, data_name, H5T_NATIVE_DOUBLE, space_id, data_id, ierr)
        if (ierr /= 0) return
        call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, data, data_shape, ierr)
        if (ierr /= 0) return
        call h5dclose_f(data_id, ierr)
        if (ierr /= 0) return
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_create_dataset_double_3d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 3

        class(t_h5node), intent(inout) ::self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        double precision, intent(in) :: data(data_shape(1), data_shape(2), data_shape(3))
        integer, intent(out) :: ierr

        integer(kind=HID_T) data_id, space_id

        call h5screate_simple_f(rank, data_shape, space_id, ierr)
        if (ierr /= 0) return
        call h5dcreate_f(self%id, data_name, H5T_NATIVE_DOUBLE, space_id, data_id, ierr)
        if (ierr /= 0) return
        call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, data, data_shape, ierr)
        if (ierr /= 0) return
        call h5dclose_f(data_id, ierr)
        if (ierr /= 0) return
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_get_dataset_integer_1d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 1

        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        integer, intent(out) :: data(data_shape(1))
        integer, intent(out) :: ierr

        integer(kind=HID_T) :: data_id, space_id
        integer :: space_rank
        integer(kind=HSIZE_T) :: space_dim(rank), max_space_dim(rank)

        call h5dopen_f(self%id, data_name, data_id, ierr)

        ! Pre-check the data to be read.
        call h5dget_space_f(data_id, space_id, ierr)

        ! Data rank check
        call h5sget_simple_extent_ndims_f(space_id, space_rank, ierr)
        if (space_rank .ne. rank) then
            write (0, *) "[ERROR@read1i] RANK of dataset is invalid:", space_rank
            write (0, *) "[ERROR@read1i] no data are outputted"
            return
        end if

        ! Data shape check
        call h5sget_simple_extent_dims_f(space_id, space_dim, max_space_dim, ierr)
        if (space_dim(1) .gt. data_shape(1)) then
            write (0, *) "[ERROR@read1i] DIM of dataset is too large:", space_dim(1)
            write (0, *) "[ERROR@read1i] no data are outputted"
            return
        end if

        call h5dread_f(data_id, H5T_NATIVE_INTEGER, data, data_shape, ierr)

        call h5dclose_f(data_id, ierr)
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_get_dataset_integer_2d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 2

        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        integer, intent(out) :: data(data_shape(1), data_shape(2))
        integer, intent(out) :: ierr

        integer(kind=HID_T) :: data_id, space_id
        integer :: space_rank
        integer(kind=HSIZE_T) :: space_dim(rank), max_space_dim(rank)

        call h5dopen_f(self%id, data_name, data_id, ierr)

        ! Pre-check the data to be read.
        call h5dget_space_f(data_id, space_id, ierr)

        ! Data rank check.
        call h5sget_simple_extent_ndims_f(space_id, space_rank, ierr)
        if (space_rank .ne. rank) then
            write (0, *) "[ERROR@read2i] RANK of dataset is invalid:", space_rank
            write (0, *) "[ERROR@read2i] no data are outputted"
            return
        end if

        ! Data shape check.
        call h5sget_simple_extent_dims_f(space_id, space_dim, max_space_dim, ierr)
        if (space_dim(1) > data_shape(1) &
            .or. space_dim(2) > data_shape(2)) then
            write (0, *) "[ERROR@read2i] DIM of dataset is too large:", space_dim(1)
            write (0, *) "[ERROR@read2i] no data are outputted"
            return
        end if

        call h5dread_f(data_id, H5T_NATIVE_INTEGER, data, data_shape, ierr)

        call h5dclose_f(data_id, ierr)
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_get_dataset_integer_3d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 3

        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        integer, intent(out) :: data(data_shape(1), data_shape(2), data_shape(3))
        integer, intent(out) :: ierr

        integer(kind=HID_T) :: data_id, space_id
        integer :: space_rank
        integer(kind=HSIZE_T) :: space_dim(rank), max_space_dim(rank)

        call h5dopen_f(self%id, data_name, data_id, ierr)

        ! Pre-check the data to be read.
        call h5dget_space_f(data_id, space_id, ierr)

        ! Data rank check.
        call h5sget_simple_extent_ndims_f(space_id, space_rank, ierr)
        if (space_rank .ne. rank) then
            write (0, *) "[ERROR@read3i] RANK of dataset is invalid:", space_rank
            write (0, *) "[ERROR@read3i] no data are outputted"
            return
        end if

        ! Data shape check.
        call h5sget_simple_extent_dims_f(space_id, space_dim, max_space_dim, ierr)
        if (space_dim(1) > data_shape(1) &
            .or. space_dim(2) > data_shape(2) &
            .or. space_dim(3) > data_shape(3)) then
            write (0, *) "[ERROR@read3i] DIM of dataset is too large:", space_dim(1)
            write (0, *) "[ERROR@read3i] no data are outputted"
            return
        end if

        call h5dread_f(data_id, H5T_NATIVE_INTEGER, data, data_shape, ierr)

        call h5dclose_f(data_id, ierr)
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_get_dataset_real_1d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 1

        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        real, intent(out) :: data(data_shape(1))
        integer, intent(out) :: ierr

        integer(kind=HID_T) :: data_id, space_id
        integer :: space_rank
        integer(kind=HSIZE_T) :: space_dim(rank), max_space_dim(rank)

        call h5dopen_f(self%id, data_name, data_id, ierr)

        ! Pre-check the data to be read.
        call h5dget_space_f(data_id, space_id, ierr)

        ! Data rank check
        call h5sget_simple_extent_ndims_f(space_id, space_rank, ierr)
        if (space_rank .ne. rank) then
            write (0, *) "[ERROR@read1r] RANK of dataset is invalid:", space_rank
            write (0, *) "[ERROR@read1r] no data are outputted"
            return
        end if

        ! Data shape check
        call h5sget_simple_extent_dims_f(space_id, space_dim, max_space_dim, ierr)
        if (space_dim(1) .gt. data_shape(1)) then
            write (0, *) "[ERROR@read1r] DIM of dataset is too large:", space_dim(1)
            write (0, *) "[ERROR@read1r] no data are outputted"
            return
        end if

        call h5dread_f(data_id, H5T_NATIVE_REAL, data, data_shape, ierr)

        call h5dclose_f(data_id, ierr)
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_get_dataset_real_2d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 2

        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        real, intent(out) :: data(data_shape(1), data_shape(2))
        integer, intent(out) :: ierr

        integer(kind=HID_T) :: data_id, space_id
        integer :: space_rank
        integer(kind=HSIZE_T) :: space_dim(rank), max_space_dim(rank)

        call h5dopen_f(self%id, data_name, data_id, ierr)

        ! Pre-check the data to be read.
        call h5dget_space_f(data_id, space_id, ierr)

        ! Data rank check.
        call h5sget_simple_extent_ndims_f(space_id, space_rank, ierr)
        if (space_rank .ne. rank) then
            write (0, *) "[ERROR@read2r] RANK of dataset is invalid:", space_rank
            write (0, *) "[ERROR@read2r] no data are outputted"
            return
        end if

        ! Data shape check.
        call h5sget_simple_extent_dims_f(space_id, space_dim, max_space_dim, ierr)
        if (space_dim(1) > data_shape(1) &
            .or. space_dim(2) > data_shape(2)) then
            write (0, *) "[ERROR@read2r] DIM of dataset is too large:", space_dim(1)
            write (0, *) "[ERROR@read2r] no data are outputted"
            return
        end if

        call h5dread_f(data_id, H5T_NATIVE_REAL, data, data_shape, ierr)

        call h5dclose_f(data_id, ierr)
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_get_dataset_real_3d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 3

        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        real, intent(out) :: data(data_shape(1), data_shape(2), data_shape(3))
        integer, intent(out) :: ierr

        integer(kind=HID_T) :: data_id, space_id
        integer :: space_rank
        integer(kind=HSIZE_T) :: space_dim(rank), max_space_dim(rank)

        call h5dopen_f(self%id, data_name, data_id, ierr)

        ! Pre-check the data to be read.
        call h5dget_space_f(data_id, space_id, ierr)

        ! Data rank check.
        call h5sget_simple_extent_ndims_f(space_id, space_rank, ierr)
        if (space_rank .ne. rank) then
            write (0, *) "[ERROR@read3r] RANK of dataset is invalid:", space_rank
            write (0, *) "[ERROR@read3r] no data are outputted"
            return
        end if

        ! Data shape check.
        call h5sget_simple_extent_dims_f(space_id, space_dim, max_space_dim, ierr)
        if (space_dim(1) > data_shape(1) &
            .or. space_dim(2) > data_shape(2) &
            .or. space_dim(3) > data_shape(3)) then
            write (0, *) "[ERROR@read3r] DIM of dataset is too large:", space_dim(1)
            write (0, *) "[ERROR@read3r] no data are outputted"
            return
        end if

        call h5dread_f(data_id, H5T_NATIVE_REAL, data, data_shape, ierr)

        call h5dclose_f(data_id, ierr)
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_get_dataset_double_1d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 1

        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        double precision, intent(out) :: data(data_shape(1))
        integer, intent(out) :: ierr

        integer(kind=HID_T) :: data_id, space_id
        integer :: space_rank
        integer(kind=HSIZE_T) :: space_dim(rank), max_space_dim(rank)

        call h5dopen_f(self%id, data_name, data_id, ierr)

        ! Pre-check the data to be read.
        call h5dget_space_f(data_id, space_id, ierr)

        ! Data rank check
        call h5sget_simple_extent_ndims_f(space_id, space_rank, ierr)
        if (space_rank .ne. rank) then
            write (0, *) "[ERROR@read1d] RANK of dataset is invalid:", space_rank
            write (0, *) "[ERROR@read1d] no data are outputted"
            return
        end if

        ! Data shape check
        call h5sget_simple_extent_dims_f(space_id, space_dim, max_space_dim, ierr)
        if (space_dim(1) .gt. data_shape(1)) then
            write (0, *) "[ERROR@read1d] DIM of dataset is too large:", space_dim(1)
            write (0, *) "[ERROR@read1d] no data are outputted"
            return
        end if

        call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, data_shape, ierr)

        call h5dclose_f(data_id, ierr)
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_get_dataset_double_2d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 2

        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        double precision, intent(out) :: data(data_shape(1), data_shape(2))
        integer, intent(out) :: ierr

        integer(kind=HID_T) :: data_id, space_id
        integer :: space_rank
        integer(kind=HSIZE_T) :: space_dim(rank), max_space_dim(rank)

        call h5dopen_f(self%id, data_name, data_id, ierr)

        ! Pre-check the data to be read.
        call h5dget_space_f(data_id, space_id, ierr)

        ! Data rank check.
        call h5sget_simple_extent_ndims_f(space_id, space_rank, ierr)
        if (space_rank .ne. rank) then
            write (0, *) "[ERROR@read2d] RANK of dataset is invalid:", space_rank
            write (0, *) "[ERROR@read2d] no data are outputted"
            return
        end if

        ! Data shape check.
        call h5sget_simple_extent_dims_f(space_id, space_dim, max_space_dim, ierr)
        if (space_dim(1) > data_shape(1) &
            .or. space_dim(2) > data_shape(2)) then
            write (0, *) "[ERROR@read2d] DIM of dataset is too large:", space_dim(1)
            write (0, *) "[ERROR@read2d] no data are outputted"
            return
        end if

        call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, data_shape, ierr)

        call h5dclose_f(data_id, ierr)
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_get_dataset_double_3d(self, data_name, data_shape, data, ierr)
        integer, parameter :: rank = 3

        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: data_shape(rank)
        double precision, intent(out) :: data(data_shape(1), data_shape(2), data_shape(3))
        integer, intent(out) :: ierr

        integer(kind=HID_T) :: data_id, space_id
        integer :: space_rank
        integer(kind=HSIZE_T) :: space_dim(rank), max_space_dim(rank)

        call h5dopen_f(self%id, data_name, data_id, ierr)

        ! Pre-check the data to be read.
        call h5dget_space_f(data_id, space_id, ierr)

        ! Data rank check.
        call h5sget_simple_extent_ndims_f(space_id, space_rank, ierr)
        if (space_rank .ne. rank) then
            write (0, *) "[ERROR@read3d] RANK of dataset is invalid:", space_rank
            write (0, *) "[ERROR@read3d] no data are outputted"
            return
        end if

        ! Data shape check.
        call h5sget_simple_extent_dims_f(space_id, space_dim, max_space_dim, ierr)
        if (space_dim(1) > data_shape(1) &
            .or. space_dim(2) > data_shape(2) &
            .or. space_dim(3) > data_shape(3)) then
            write (0, *) "[ERROR@read3d] DIM of dataset is too large:", space_dim(1)
            write (0, *) "[ERROR@read3d] no data are outputted"
            return
        end if

        call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, data_shape, ierr)

        call h5dclose_f(data_id, ierr)
        call h5sclose_f(space_id, ierr)
    end subroutine

    subroutine h5node_create_dataset_double_3d_parallel(self, data_name, &
                                                        local_data_shape, local_data, &
                                                        offset, all_data_shape, &
                                                        ierr)
        integer, parameter :: rank = 3

        class(t_h5node), intent(inout) :: self
        character(len=*), intent(in) :: data_name
        integer(kind=HSIZE_T), intent(in) :: local_data_shape(rank)
        double precision, intent(in) :: local_data(local_data_shape(1), local_data_shape(2), local_data_shape(3))
        integer(kind=HSIZE_T), intent(in) :: offset(rank)
        integer(kind=HSIZE_T), intent(in) :: all_data_shape(rank)
        integer, intent(out) :: ierr

        integer(kind=HID_T) :: data_id
        integer(kind=HID_T) :: file_space_id, mem_space_id, xfer_prp

        ! Create dataspace (file_space) and dataset for all data.
        call h5screate_simple_f(rank, all_data_shape, file_space_id, ierr)
        call h5dcreate_f(self%id, data_name, H5T_NATIVE_DOUBLE, file_space_id, data_id, ierr)
        call h5sclose_f(file_space_id, ierr)

        ! Create dataspace for local data (mem_space) + get dataset (data_id) + set offset for local index.
        call h5screate_simple_f(rank, local_data_shape, mem_space_id, ierr)
        call h5dget_space_f(data_id, file_space_id, ierr)
        call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset, local_data_shape, ierr)

        ! Generate transfer property list (xfer_prp).
        call h5pcreate_f(H5P_DATASET_XFER_F, xfer_prp, ierr)
        call h5pset_dxpl_mpio_f(xfer_prp, H5FD_MPIO_COLLECTIVE_F, ierr)

        ! Write local data in each process.
        call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, local_data, all_data_shape, ierr, &
                        file_space_id=file_space_id, mem_space_id=mem_space_id, xfer_prp=xfer_prp)

        ! Close dataset and dataspace.
        call h5sclose_f(file_space_id, ierr)
        call h5sclose_f(mem_space_id, ierr)
        call h5dclose_f(data_id, ierr)
        call h5pclose_f(xfer_prp, ierr)
    end subroutine

end module
