  module hdf5io
!
    implicit none
!
  contains

  subroutine hdfinit()
    use HDF5
    implicit none
!
    integer(kind=4) :: hdferr
!
      call h5open_f(hdferr)
      call h5eset_auto_f(0, hdferr)
!
    return
  end subroutine hdfinit


  subroutine hdffinalize()
    use HDF5
    implicit none
!
    integer(kind=4) :: hdferr
!
      call h5close_f(hdferr)
!
    return
  end subroutine hdffinalize


  subroutine hdfopen(filename,id,acc)
    use HDF5
    implicit none
!
    character(len=*),intent(in) :: filename
    integer(kind=HID_T),intent(out) :: id
    integer(kind=4),intent(in) :: acc
!
!    integer(kind=4) :: n
    integer(kind=4) :: hdferr
!    integer(kind=4),external :: lnblnk
!     
!      n  = lnblnk(filename)
      if(acc .eq. 1) then
         call h5fopen_f(filename, H5F_ACC_RDONLY_F, id, hdferr)
      else if(acc .eq. 4) then
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, id, hdferr)
      else
         call h5fopen_f(filename, H5F_ACC_RDWR_F, id, hdferr)
      end if
!     
    return
  end subroutine hdfopen


  subroutine hdfopen_p(filename, id, acc, comm)
    use HDF5
    include 'mpif.h'
    character(len=*),intent(in):: filename
    integer(kind=HID_T),intent(out):: id
    integer(kind=4),intent(in) :: acc
    integer(kind=4),intent(in) :: comm
!
    integer(kind=4) :: n
    integer(kind=4) :: hdferr
    integer(kind=HID_T):: plist_id
    integer(kind=4),external :: lnblnk
!
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
    call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, hdferr)

    n  = lnblnk(filename)
    if(acc .eq. 1) then
       call h5fopen_f(filename(1:n), H5F_ACC_RDONLY_F, id, hdferr, access_prp = plist_id)
    else if(acc .eq. 4) then
       call h5fcreate_f(filename(1:n), H5F_ACC_TRUNC_F, id, hdferr, access_prp = plist_id)
    else
       call h5fopen_f(filename(1:n), H5F_ACC_RDWR_F, id, hdferr, access_prp = plist_id)
    endif

    call h5pclose_f(plist_id, hdferr)
!
    return
  end subroutine hdfopen_p


  subroutine hdfopen_pg(filename, groupname, f_id, g_id, acc, comm)
    use HDF5
    include 'mpif.h'
    character(len=*),intent(in):: filename, groupname
    integer(kind=HID_T),intent(out):: f_id, g_id
    integer(kind=4),intent(in) :: acc
    integer(kind=4),intent(in) :: comm
!
    integer(kind=4) :: hdferr
    integer(kind=HID_T):: plist_id
!
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
    call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, hdferr)

    if(acc .eq. 1) then
       call h5fopen_f(filename, H5F_ACC_RDONLY_F, f_id, hdferr, access_prp = plist_id)
       call h5gopen_f(f_id, groupname, g_id, hdferr)
    else if(acc .eq. 4) then
       call h5fcreate_f(filename, H5F_ACC_TRUNC_F, f_id, hdferr, access_prp = plist_id)
       call h5gcreate_f(f_id, groupname, g_id, hdferr)
    else
       call h5fopen_f(filename, H5F_ACC_RDWR_F, f_id, hdferr, access_prp = plist_id)
       call h5gopen_f(f_id, groupname, g_id, hdferr)
    endif

    call h5pclose_f(plist_id, hdferr)
!
    return
  end subroutine hdfopen_pg


  subroutine hdfclose(id,istate)
    use HDF5
    implicit none
!
    integer(kind=HID_T),intent(in) :: id
    integer(kind=4),intent(out) :: istate
!
      call h5fclose_f(id, istate)
!     
    return
  end subroutine hdfclose


  subroutine hdfclose_g(f_id,g_id,istate1,istate2)
    use HDF5
    implicit none
!
    integer(kind=HID_T),intent(in) :: f_id, g_id
    integer(kind=4),intent(out) :: istate1, istate2
!
      call h5gclose_f(g_id, istate2)
      call h5fclose_f(f_id, istate1)
!     
    return
  end subroutine hdfclose_g


  subroutine read1i(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 1
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    integer(kind=4),intent(out) :: data(dim(1))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) d_id, s_id
    integer(kind=4) datatype
    integer(kind=4) frank
    integer(kind=HSIZE_T) :: fdim(rank), mdim(rank)
!
      call h5dopen_f(f_id, dsname, d_id, status)
      if(status .ne. 0) then
         print*, "[ERROR@read1i] dataset name is invalid:", dsname
         return
      end if
!
      call h5dget_space_f(d_id, s_id, status)
      ! data rank check
      call h5sget_simple_extent_ndims_f(s_id, frank, status)
      if(frank .ne. rank) then
         print*, "[ERROR@read1i] RANK of dataset is invalid:", frank
         print*, "[ERROR@read1i] no data are outputted"
         return
      end if
      ! data shape check
      call h5sget_simple_extent_dims_f(s_id, fdim, mdim, status)
      if(fdim(1).gt.dim(1)) then
         print*, "[ERROR@read1i] DIM of dataset is too large:", fdim
         print*, "[ERROR@read1i] no data are outputted"
         return
      end if
!
      call h5dread_f(d_id, H5T_NATIVE_INTEGER, data, dim, istat)
      call h5sclose_f(s_id, status)
      call h5dclose_f(d_id, status)
!
    return
  end subroutine read1i


  subroutine read2i(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 2
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    integer(kind=4),intent(out) :: data(dim(1),dim(2))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) d_id, s_id
    integer(kind=4) datatype
    integer(kind=4) frank
    integer(kind=HSIZE_T) :: fdim(rank), mdim(rank)
!
      call h5dopen_f(f_id, dsname, d_id, status)
      if(status .ne. 0) then
         print*, "[ERROR@read2i] dataset name is invalid:", dsname
         return
      end if
!
      call h5dget_space_f(d_id, s_id, status)
      ! data rank check
      call h5sget_simple_extent_ndims_f(s_id, frank, status)
      if(frank .ne. rank) then
         print*, "[ERROR@read2i] RANK of dataset is invalid:", frank
         print*, "[ERROR@read2i] no data are outputted"
         return
      end if
      ! data shape check
      call h5sget_simple_extent_dims_f(s_id, fdim, mdim, status)
      if(fdim(1).gt.dim(1).or.fdim(2).gt.dim(2)) then
         print*, "[ERROR@read2i] DIM of dataset is too large:", fdim
         print*, "[ERROR@read2i] no data are outputted"
         return
      end if
!
      call h5dread_f(d_id, H5T_NATIVE_INTEGER, data, dim, istat)
      call h5sclose_f(s_id, status)
      call h5dclose_f(d_id, status)
!
    return
  end subroutine read2i


  subroutine read3i(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 3
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    integer(kind=4),intent(out) :: data(dim(1),dim(2),dim(3))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) d_id, s_id
    integer(kind=4) datatype
    integer(kind=4) frank
    integer(kind=HSIZE_T) :: fdim(rank), mdim(rank)
!
      call h5dopen_f(f_id, dsname, d_id, status)
      if(status .ne. 0) then
         print*, "[ERROR@read3i] dataset name is invalid:", dsname
         return
      end if
!
      call h5dget_space_f(d_id, s_id, status)
      ! data rank check
      call h5sget_simple_extent_ndims_f(s_id, frank, status)
      if(frank .ne. rank) then
         print*, "[ERROR@read3i] RANK of dataset is invalid:", frank
         print*, "[ERROR@read3i] no data are outputted"
         return
      end if
      ! data shape check
      call h5sget_simple_extent_dims_f(s_id, fdim, mdim, status)
      if(fdim(1).gt.dim(1).or.fdim(2).gt.dim(2).or.fdim(3).gt.dim(3)) then
         print*, "[ERROR@read3i] DIM of dataset is too large:", fdim
         print*, "[ERROR@read3i] no data are outputted"
         return
      end if
!
      call h5dread_f(d_id, H5T_NATIVE_INTEGER, data, dim, istat)
      call h5sclose_f(s_id, status)
      call h5dclose_f(d_id, status)
!
    return
  end subroutine read3i


  subroutine read1r(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 1
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=4),intent(out) :: data(dim(1))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) d_id, s_id
    integer(kind=4) datatype
    integer(kind=4) frank
    integer(kind=HSIZE_T) :: fdim(rank), mdim(rank)
!
      call h5dopen_f(f_id, dsname, d_id, status)
      if(status .ne. 0) then
         print*, "[ERROR@read1r] dataset name is invalid:", dsname
         return
      end if
      call h5dget_space_f(d_id, s_id, status)
      ! data rank check
      call h5sget_simple_extent_ndims_f(s_id, frank, status)
      if(frank .ne. rank) then
         print*, "[ERROR@read1r] RANK of dataset is invalid:", frank
         print*, "[ERROR@read1r] no data are outputted"
         return
      end if
      ! data shape check
      call h5sget_simple_extent_dims_f(s_id, fdim, mdim, status)
      if(fdim(1).gt.dim(1)) then
         print*, "[ERROR@read1r] DIM of dataset is too large:", fdim
         print*, "[ERROR@read1r] no data are outputted"
         return
      end if
!
      call h5dread_f(d_id, H5T_NATIVE_REAL, data, dim, istat)
      call h5sclose_f(s_id, status)
      call h5dclose_f(d_id, status)
!
    return
  end subroutine read1r


  subroutine read2r(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 2
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=4),intent(out) :: data(dim(1),dim(2))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) d_id, s_id
    integer(kind=4) datatype
    integer(kind=4) frank
    integer(kind=HSIZE_T) :: fdim(rank), mdim(rank)
!
      call h5dopen_f(f_id, dsname, d_id, status)
      if(status .ne. 0) then
         print*, "[ERROR@read2r] dataset name is invalid:", dsname
         return
      end if
!
      call h5dget_space_f(d_id, s_id, status)
      ! data rank check
      call h5sget_simple_extent_ndims_f(s_id, frank, status)
      if(frank .ne. rank) then
         print*, "[ERROR@read2r] RANK of dataset is invalid:", frank
         print*, "[ERROR@read2r] no data are outputted"
         return
      end if
      ! data shape check
      call h5sget_simple_extent_dims_f(s_id, fdim, mdim, status)
      if(fdim(1).gt.dim(1).or.fdim(2).gt.dim(2)) then
         print*, "[ERROR@read2r] DIM of dataset is too large:", fdim
         print*, "[ERROR@read2r] no data are outputted"
         return
      end if
!
      call h5dread_f(d_id, H5T_NATIVE_REAL, data, dim, istat)
      call h5sclose_f(s_id, status)
      call h5dclose_f(d_id, status)
!
    return
  end subroutine read2r


  subroutine read3r(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 3
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=4),intent(out) :: data(dim(1),dim(2),dim(3))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) d_id, s_id
    integer(kind=4) datatype
    integer(kind=4) frank
    integer(kind=HSIZE_T) :: fdim(rank), mdim(rank)
!
      call h5dopen_f(f_id, dsname, d_id, status)
      if(status .ne. 0) then
         print*, "[ERROR@read3r] dataset name is invalid:", dsname
         return
      end if
!
      call h5dget_space_f(d_id, s_id, status)
      ! data rank check
      call h5sget_simple_extent_ndims_f(s_id, frank, status)
      if(frank .ne. rank) then
         print*, "[ERROR@read3r] RANK of dataset is invalid:", frank
         print*, "[ERROR@read3r] no data are outputted"
         return
      end if
      ! data shape check
      call h5sget_simple_extent_dims_f(s_id, fdim, mdim, status)
      if(fdim(1).gt.dim(1).or.fdim(2).gt.dim(2).or.fdim(3).gt.dim(3)) then
         print*, "[ERROR@read3r] DIM of dataset is too large:", fdim
         print*, "[ERROR@read3r] no data are outputted"
         return
      end if
!
      call h5dread_f(d_id, H5T_NATIVE_REAL, data, dim, istat)
      call h5sclose_f(s_id, status)
      call h5dclose_f(d_id, status)
!
    return
  end subroutine read3r


  subroutine read1d(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 1
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=8),intent(out) :: data(dim(1))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) d_id, s_id, t_id
    integer(kind=HID_T) datatype
    integer(kind=4) frank
    integer(kind=HSIZE_T) :: fdim(rank), mdim(rank)
      !
      call h5dopen_f(f_id, dsname, d_id, status)
      if(status .ne. 0) then
         print*, "[ERROR@read1d] dataset name is invalid:", dsname
         return
      end if
!
      call h5dget_space_f(d_id, s_id, status)
      ! data rank check
      call h5sget_simple_extent_ndims_f(s_id, frank, status)
      if(frank .ne. rank) then
         print*, "[ERROR@read1d] RANK of dataset is invalid:", frank
         print*, "[ERROR@read1d] no data are outputted"
         return
      end if
      ! data shape check
      call h5sget_simple_extent_dims_f(s_id, fdim, mdim, status)
      if(fdim(1).gt.dim(1)) then
         print*, "[ERROR@read1d] DIM of dataset is too large:", fdim
         print*, "[ERROR@read1d] no data are outputted"
         return
      end if
!
      call h5dread_f(d_id, H5T_NATIVE_DOUBLE, data, dim, istat)
      call h5sclose_f(s_id, status)
      call h5dclose_f(d_id, status)
!
    return
  end subroutine read1d


  subroutine read2d(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 2
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=8),intent(out) :: data(dim(1),dim(2))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) d_id, s_id, t_id
    integer(kind=HID_T) datatype
    integer(kind=4) frank
    integer(kind=HSIZE_T) :: fdim(rank), mdim(rank)
!
      call h5dopen_f(f_id, dsname, d_id, status)
      if(status .ne. 0) then
         print*, "[ERROR@read2d] dataset name is invalid:", dsname
         return
      end if
!
      call h5dget_space_f(d_id, s_id, status)
      ! data rank check
      call h5sget_simple_extent_ndims_f(s_id, frank, status)
      if(frank .ne. rank) then
         print*, "[ERROR@read2d] RANK of dataset is invalid:", frank
         print*, "[ERROR@read2d] no data are outputted"
         return
      end if
      ! data shape check
      call h5sget_simple_extent_dims_f(s_id, fdim, mdim, status)
      if(fdim(1).gt.dim(1).or.fdim(2).gt.dim(2)) then
         print*, "[ERROR@read2d] DIM of dataset is too large:", fdim
         print*, "[ERROR@read2d] no data are outputted"
         return
      end if
!
      call h5dread_f(d_id, H5T_NATIVE_DOUBLE, data, dim, istat)
      call h5sclose_f(s_id, status)
      call h5dclose_f(d_id, status)
!
    return
  end subroutine read2d


  subroutine read3d(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 3
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=8),intent(out) :: data(dim(1),dim(2),dim(3))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) d_id, s_id, t_id
    integer(kind=HID_T) datatype
    integer(kind=4) frank
    integer(kind=HSIZE_T) :: fdim(rank), mdim(rank)
!
      call h5dopen_f(f_id, dsname, d_id, status)
      if(status .ne. 0) then
         print*, "[ERROR@read3d] dataset name is invalid:", dsname
         return
      end if
!
      call h5dget_space_f(d_id, s_id, status)
      ! data rank check
      call h5sget_simple_extent_ndims_f(s_id, frank, status)
      if(frank .ne. rank) then
         print*, "[ERROR@read3d] RANK of dataset is invalid:", frank
         print*, "[ERROR@read3d] no data are outputted"
         return
      end if
      ! data shape check
      call h5sget_simple_extent_dims_f(s_id, fdim, mdim, status)
      if(fdim(1).gt.dim(1).or.fdim(2).gt.dim(2).or.fdim(3).gt.dim(3)) then
         print*, "[ERROR@read3d] DIM of dataset is too large:", fdim
         print*, "[ERROR@read3d] no data are outputted"
         return
      end if
!
      call h5dread_f(d_id, H5T_NATIVE_DOUBLE, data, dim, istat)
      call h5sclose_f(s_id, status)
      call h5dclose_f(d_id, status)
!
    return
  end subroutine read3d


  subroutine wrt1i(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
    integer(kind=4),parameter:: rank=1
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    integer(kind=4),intent(in) :: data(dim(1))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) s_id, d_id
!
      istat = -1
      call h5screate_simple_f(rank, dim, s_id, status)
      if(status .ne. 0) return
      call h5dcreate_f(f_id, dsname, H5T_NATIVE_INTEGER, s_id, d_id, status)
      if(status .ne. 0) return
      call h5dwrite_f(d_id, H5T_NATIVE_INTEGER, data, dim, istat)
      if(istat .ne. 0) return
      call h5dclose_f(d_id, status)
      if(status .ne. 0) return
      call h5sclose_f(s_id, status)
!
    return
  end subroutine wrt1i


  subroutine wrt2i(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
    integer(kind=4),parameter:: rank=2
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    integer(kind=4),intent(in) :: data(dim(1),dim(2))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) s_id, d_id
!
      istat = -1
      call h5screate_simple_f(rank, dim, s_id, status)
      if(status .ne. 0) return
      call h5dcreate_f(f_id, dsname, H5T_NATIVE_INTEGER, s_id, d_id, status)
      if(status .ne. 0) return
      call h5dwrite_f(d_id, H5T_NATIVE_INTEGER, data, dim, istat)
      if(istat .ne. 0) return
      call h5dclose_f(d_id, status)
      if(status .ne. 0) return
      call h5sclose_f(s_id, status)
!
    return
  end subroutine wrt2i


  subroutine wrt3i(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
    integer(kind=4),parameter:: rank=3
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    integer(kind=4),intent(in) :: data(dim(1),dim(2),dim(3))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) s_id, d_id
!
      istat = -1
      call h5screate_simple_f(rank, dim, s_id, status)
      if(status .ne. 0) return
      call h5dcreate_f(f_id, dsname, H5T_NATIVE_INTEGER, s_id, d_id, status)
      if(status .ne. 0) return
      call h5dwrite_f(d_id, H5T_NATIVE_INTEGER, data, dim, istat)
      if(istat .ne. 0) return
      call h5dclose_f(d_id, status)
      if(status .ne. 0) return
      call h5sclose_f(s_id, status)
!
    return
  end subroutine wrt3i


  subroutine wrt1r(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
    integer(kind=4),parameter:: rank=1
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=4),intent(in) :: data(dim(1))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) s_id, d_id
!
      istat = -1
      call h5screate_simple_f(rank, dim, s_id, status)
      if(status .ne. 0) return
      call h5dcreate_f(f_id, dsname, H5T_NATIVE_REAL, s_id, d_id, status)
      if(status .ne. 0) return
      call h5dwrite_f(d_id, H5T_NATIVE_REAL, data, dim, istat)
      if(istat .ne. 0) return
      call h5dclose_f(d_id, status)
      if(status .ne. 0) return
      call h5sclose_f(s_id, status)
!
    return
  end subroutine wrt1r


  subroutine wrt2r(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
    integer(kind=4),parameter:: rank=2
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=4),intent(in) :: data(dim(1),dim(2))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) s_id, d_id
!
      istat = -1
      call h5screate_simple_f(rank, dim, s_id, status)
      if(status .ne. 0) return
      call h5dcreate_f(f_id, dsname, H5T_NATIVE_REAL, s_id, d_id, status)
      if(status .ne. 0) return
      call h5dwrite_f(d_id, H5T_NATIVE_REAL, data, dim, istat)
      if(istat .ne. 0) return
      call h5dclose_f(d_id, status)
      if(status .ne. 0) return
      call h5sclose_f(s_id, status)
!
    return
  end subroutine wrt2r


  subroutine wrt3r(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
    integer(kind=4),parameter:: rank=3
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=4),intent(in) :: data(dim(1),dim(2),dim(3))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) s_id, d_id
!
      istat = -1
      call h5screate_simple_f(rank, dim, s_id, status)
      if(status .ne. 0) return
      call h5dcreate_f(f_id, dsname, H5T_NATIVE_REAL, s_id, d_id, status)
      if(status .ne. 0) return
      call h5dwrite_f(d_id, H5T_NATIVE_REAL, data, dim, istat)
      if(istat .ne. 0) return
      call h5dclose_f(d_id, status)
      if(status .ne. 0) return
      call h5sclose_f(s_id, status)
!
    return
  end subroutine wrt3r


  subroutine wrt1d(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
    integer(kind=4),parameter:: rank=1
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=8),intent(in) :: data(dim(1))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) s_id, d_id
!
      istat = -1
      call h5screate_simple_f(rank, dim, s_id, status)
      if(status .ne. 0) return
      call h5dcreate_f(f_id, dsname, H5T_NATIVE_DOUBLE, s_id, d_id, status)
      if(status .ne. 0) return
      call h5dwrite_f(d_id, H5T_NATIVE_DOUBLE, data, dim, istat)
      if(istat .ne. 0) return
      call h5dclose_f(d_id, status)
      if(status .ne. 0) return
      call h5sclose_f(s_id, status)
!
    return
  end subroutine wrt1d


  subroutine wrt2d(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
    integer(kind=4),parameter:: rank=2
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=8),intent(in) :: data(dim(1),dim(2))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) s_id, d_id
!
      istat = -1
      call h5screate_simple_f(rank, dim, s_id, status)
      if(status .ne. 0) return
      call h5dcreate_f(f_id, dsname, H5T_NATIVE_DOUBLE, s_id, d_id, status)
      if(status .ne. 0) return
      call h5dwrite_f(d_id, H5T_NATIVE_DOUBLE, data, dim, istat)
      if(istat .ne. 0) return
      call h5dclose_f(d_id, status)
      if(status .ne. 0) return
      call h5sclose_f(s_id, status)
!
    return
  end subroutine wrt2d


  subroutine wrt3d(f_id,dsname,dim,data,istat,status)
    use HDF5
    implicit none
    integer(kind=4),parameter:: rank=3
!
    integer(kind=HID_T),intent(in) :: f_id
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, status
    real(kind=8),intent(in) :: data(dim(1),dim(2),dim(3))
    character(len=*),intent(in) :: dsname
!
    integer(kind=HID_T) s_id, d_id
!
      istat = -1
      call h5screate_simple_f(rank, dim, s_id, status)
      if(status .ne. 0) return
      call h5dcreate_f(f_id, dsname, H5T_NATIVE_DOUBLE, s_id, d_id, status)
      if(status .ne. 0) return
      call h5dwrite_f(d_id, H5T_NATIVE_DOUBLE, data, dim, istat)
      if(istat .ne. 0) return
      call h5dclose_f(d_id, status)
      if(status .ne. 0) return
      call h5sclose_f(s_id, status)
!
    return
  end subroutine wrt3d


  subroutine fwrt1i(filename,sdsname,dim,sdsdata,istat,jstat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 1
!
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, jstat, status
    integer(kind=4),intent(in) :: sdsdata(dim(1))
    character(len=*),intent(in) :: filename, sdsname
!
!    integer(kind=4) :: n
    integer(kind=HID_T) :: sd_id
    integer(kind=4) :: hdferr
!    integer(kind=4),external :: lnblnk
!     
!      n  = lnblnk(filename)
      istat = -1
      jstat = -1
      call h5fopen_f(filename, H5F_ACC_RDWR_F, sd_id, hdferr)
      if(hdferr.ne.0) then
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, sd_id, status)
      end if
      if(status .ne. 0) return
!     
      call wrt1i(sd_id,sdsname,dim,sdsdata,istat,jstat)
!
      call h5fclose_f(sd_id, status)
!
    return
  end subroutine fwrt1i


  subroutine fwrt2i(filename,sdsname,dim,sdsdata,istat,jstat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 2
!
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, jstat, status
    integer(kind=4),intent(in) :: sdsdata(dim(1),dim(2))
    character(len=*),intent(in) :: filename, sdsname
!
!    integer(kind=4) :: n
    integer(kind=HID_T) :: sd_id
    integer(kind=4) :: hdferr
!    integer(kind=4),external :: lnblnk
!     
!      n  = lnblnk(filename)
      istat = -1
      jstat = -1
      call h5fopen_f(filename, H5F_ACC_RDWR_F, sd_id, hdferr)
      if(hdferr.ne.0) then
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, sd_id, status)
      end if
      if(status .ne. 0) return
!     
      call wrt2i(sd_id,sdsname,dim,sdsdata,istat,jstat)
!
      call h5fclose_f(sd_id, status)
!
    return
  end subroutine fwrt2i


  subroutine fwrt3i(filename,sdsname,dim,sdsdata,istat,jstat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 3
!
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, jstat, status
    integer(kind=4),intent(in) :: sdsdata(dim(1),dim(2),dim(3))
    character(len=*),intent(in) :: filename, sdsname
!
!    integer(kind=4) :: n
    integer(kind=HID_T) :: sd_id
    integer(kind=4) :: hdferr
!    integer(kind=4),external :: lnblnk
!     
!      n  = lnblnk(filename)
      istat = -1
      jstat = -1
      call h5fopen_f(filename, H5F_ACC_RDWR_F, sd_id, hdferr)
      if(hdferr.ne.0) then
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, sd_id, status)
      end if
      if(status .ne. 0) return
!     
      call wrt3i(sd_id,sdsname,dim,sdsdata,istat,jstat)
!
      call h5fclose_f(sd_id, status)
!
    return
  end subroutine fwrt3i


  subroutine fwrt1r(filename,sdsname,dim,sdsdata,istat,jstat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 1
!
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, jstat, status
    real(kind=4),intent(in) :: sdsdata(dim(1))
    character(len=*),intent(in) :: filename, sdsname
!
!    integer(kind=4) :: n
    integer(kind=HID_T) :: sd_id
    integer(kind=4) :: hdferr
!    integer(kind=4),external :: lnblnk
!     
!      n  = lnblnk(filename)
      istat = -1
      jstat = -1
      call h5fopen_f(filename, H5F_ACC_RDWR_F, sd_id, hdferr)
      if(hdferr.ne.0) then
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, sd_id, status)
      end if
      if(status .ne. 0) return
!     
      call wrt1r(sd_id,sdsname,dim,sdsdata,istat,jstat)
!
      call h5fclose_f(sd_id, status)
!
    return
  end subroutine fwrt1r


  subroutine fwrt2r(filename,sdsname,dim,sdsdata,istat,jstat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 2
!
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, jstat, status
    real(kind=4),intent(in) :: sdsdata(dim(1),dim(2))
    character(len=*),intent(in) :: filename, sdsname
!
!    integer(kind=4) :: n
    integer(kind=HID_T) :: sd_id
    integer(kind=4) :: hdferr
!    integer(kind=4),external :: lnblnk
!     
!      n  = lnblnk(filename)
      istat = -1
      jstat = -1
      call h5fopen_f(filename, H5F_ACC_RDWR_F, sd_id, hdferr)
      if(hdferr.ne.0) then
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, sd_id, status)
      end if
      if(status .ne. 0) return
!     
      call wrt2r(sd_id,sdsname,dim,sdsdata,istat,jstat)
!
      call h5fclose_f(sd_id, status)
!
    return
  end subroutine fwrt2r


  subroutine fwrt3r(filename,sdsname,dim,sdsdata,istat,jstat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 3
!
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, jstat, status
    real(kind=4),intent(in) :: sdsdata(dim(1),dim(2),dim(3))
    character(len=*),intent(in) :: filename, sdsname
!
!    integer(kind=4) :: n
    integer(kind=HID_T) :: sd_id
    integer(kind=4) :: hdferr
!    integer(kind=4),external :: lnblnk
!     
!      n  = lnblnk(filename)
      istat = -1
      jstat = -1
      call h5fopen_f(filename, H5F_ACC_RDWR_F, sd_id, hdferr)
      if(hdferr.ne.0) then
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, sd_id, status)
      end if
      if(status .ne. 0) return
!     
      call wrt3r(sd_id,sdsname,dim,sdsdata,istat,jstat)
!
      call h5fclose_f(sd_id, status)
!
    return
  end subroutine fwrt3r


  subroutine fwrt1d(filename,sdsname,dim,sdsdata,istat,jstat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 1
!
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, jstat, status
    real(kind=8),intent(in) :: sdsdata(dim(1))
    character(len=*),intent(in) :: filename, sdsname
!
!    integer(kind=4) :: n
    integer(kind=HID_T) :: sd_id
    integer(kind=4) :: hdferr
!    integer(kind=4),external :: lnblnk
!     
!      n  = lnblnk(filename)
      istat = -1
      jstat = -1
      call h5fopen_f(filename, H5F_ACC_RDWR_F, sd_id, hdferr)
      if(hdferr.ne.0) then
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, sd_id, status)
      end if
      if(status .ne. 0) return
!     
      call wrt1d(sd_id,sdsname,dim,sdsdata,istat,jstat)
!
      call h5fclose_f(sd_id, status)
!
    return
  end subroutine fwrt1d


  subroutine fwrt2d(filename,sdsname,dim,sdsdata,istat,jstat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 2
!
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, jstat, status
    real(kind=8),intent(in) :: sdsdata(dim(1),dim(2))
    character(len=*),intent(in) :: filename, sdsname
!
!    integer(kind=4) :: n
    integer(kind=HID_T) :: sd_id
    integer(kind=4) :: hdferr
!    integer(kind=4),external :: lnblnk
!     
!      n  = lnblnk(filename)
      istat = -1
      jstat = -1
      call h5fopen_f(filename, H5F_ACC_RDWR_F, sd_id, hdferr)
      if(hdferr.ne.0) then
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, sd_id, status)
      end if
      if(status .ne. 0) return
!     
      call wrt2d(sd_id,sdsname,dim,sdsdata,istat,jstat)
!
      call h5fclose_f(sd_id, status)
!
    return
  end subroutine fwrt2d


  subroutine fwrt3d(filename,sdsname,dim,sdsdata,istat,jstat,status)
    use HDF5
    implicit none
!
    integer(kind=4),parameter :: rank = 3
!
    integer(kind=HSIZE_T),intent(in) :: dim(rank)
    integer(kind=4),intent(out) :: istat, jstat, status
    real(kind=8),intent(in) :: sdsdata(dim(1),dim(2),dim(3))
    character(len=*),intent(in) :: filename, sdsname
!
!    integer(kind=4) :: n
    integer(kind=HID_T) :: sd_id
    integer(kind=4) :: hdferr
!    integer(kind=4),external :: lnblnk
!     
!      n  = lnblnk(filename)
      istat = -1
      jstat = -1
      call h5fopen_f(filename, H5F_ACC_RDWR_F, sd_id, hdferr)
      if(hdferr.ne.0) then
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, sd_id, status)
      end if
      if(status .ne. 0) return
!     
      call wrt3d(sd_id,sdsname,dim,sdsdata,istat,jstat)
!
      call h5fclose_f(sd_id, status)
!
    return
  end subroutine fwrt3d


  subroutine wrt3d_p(f_id,dsname,myid,nodes,dsize,tab,data,istat,status)
    use HDF5
    implicit none
    integer(kind=4),parameter:: rank=3
!
    integer(kind=HID_T),intent(in) :: f_id
    character(len=*),intent(in) :: dsname
    integer(kind=4),intent(in) :: myid
    integer(kind=4),intent(in) :: nodes(rank)
    integer(kind=HSIZE_T),intent(in) :: dsize(rank)
    integer(kind=HSIZE_T),intent(in) :: tab
    real(kind=8),intent(in) :: data(-tab:dsize(1)-1+tab,-tab:dsize(2)-1+tab,-tab:dsize(3)-1+tab)
    integer(kind=4),intent(out) :: istat, status
!
    integer(kind=HID_T) s_id, d_id, filespace, memspace, p_id
    integer(kind=HSIZE_T) all_data_size(rank)
    integer(kind=HSIZE_T) local_data_size(rank)
    integer(kind=HSIZE_T) lb(rank)
    integer(kind=HSIZE_T) ub(rank)
    integer coord(rank)
    integer(kind=HSIZE_T) offset(rank)
    integer i
!
    istat = -1

    coord(1) = mod(myid, nodes(1))
    coord(2) = mod(myid/nodes(1), nodes(2))
    coord(3) = myid/(nodes(1)*nodes(2))

    do i=1, rank
       all_data_size(i) = (dsize(i)-1)*nodes(i)+2*tab+1
!       all_data_size(i) = dsize(i)*nodes(i)
    enddo

    call h5screate_simple_f(rank, all_data_size, filespace, status)
    call h5dcreate_f(f_id, dsname, H5T_NATIVE_DOUBLE, filespace, d_id, status)
    call h5sclose_f(filespace, status)

    do i=1, rank
       if(coord(i) .eq. 0.and.coord(i) .eq. nodes(i)-1) then
          lb(i) = -tab
          ub(i) = dsize(i)+tab-1
          offset(i) = 0
       else if(coord(i) .eq. 0) then
          lb(i) = -tab
          ub(i) = dsize(i)-2
          offset(i) = 0
       else if(coord(i) .eq. nodes(i)-1) then
          lb(i) = 0
          ub(i) = dsize(i)+tab-1
          offset(i) = (dsize(i)-1)*coord(i)+tab
       else
          lb(i) = 0
          ub(i) = dsize(i)-2
          offset(i) = (dsize(i)-1)*coord(i)+tab
       endif
       local_data_size(i) = ub(i)-lb(i)+1
    enddo

    call h5screate_simple_f(rank, local_data_size, memspace, status)
    call h5dget_space_f(d_id, filespace, status)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, local_data_size, status)

    call h5pcreate_f(H5P_DATASET_XFER_F, p_id, status)
    call h5pset_dxpl_mpio_f(p_id, H5FD_MPIO_COLLECTIVE_F, status)

    call h5dwrite_f(d_id, H5T_NATIVE_DOUBLE, data(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), all_data_size, istat, &
         &file_space_id=filespace, mem_space_id=memspace, xfer_prp=p_id)

    call h5sclose_f(filespace, status)
    call h5sclose_f(memspace, status)
    call h5dclose_f(d_id, status)
    call h5pclose_f(p_id, status)
!
    return
  end subroutine wrt3d_p


  subroutine wrtfield3d(f_id,dsname,myid,nodes,dsize,tab,data,istat,status)
    use HDF5
    implicit none
    integer(kind=4),parameter:: rank=3
!
    integer(kind=HID_T),intent(in) :: f_id
    character(len=*),intent(in) :: dsname
    integer(kind=4),intent(in) :: myid
    integer(kind=4),intent(in) :: nodes(rank)
    integer(kind=HSIZE_T),intent(in) :: dsize(rank)
    integer(kind=HSIZE_T),intent(in) :: tab
    real(kind=8),intent(in) :: data(-tab:dsize(1)-1+tab,-tab:dsize(2)-1+tab,-tab:dsize(3)-1+tab)
    integer(kind=4),intent(out) :: istat, status
!
    integer(kind=HID_T) s_id, d_id, filespace, memspace, p_id
    integer(kind=HSIZE_T) all_data_size(rank)
    integer(kind=HSIZE_T) local_data_size(rank)
    integer(kind=HSIZE_T) lb(rank)
    integer(kind=HSIZE_T) ub(rank)
    integer coord(rank)
    integer(kind=HSIZE_T) offset(rank)
    integer i
!
    istat = -1

    coord(1) = mod(myid, nodes(1))
    coord(2) = mod(myid/nodes(1), nodes(2))
    coord(3) = myid/(nodes(1)*nodes(2))

    do i=1, rank
       all_data_size(i) = (dsize(i)-1)*nodes(i)+2*tab+1
!       all_data_size(i) = dsize(i)*nodes(i)
    enddo

    call h5screate_simple_f(rank, all_data_size, filespace, status)
    call h5dcreate_f(f_id, dsname, H5T_NATIVE_DOUBLE, filespace, d_id, status)
    call h5sclose_f(filespace, status)

    do i=1, rank
       if(coord(i) .eq. 0.and.coord(i) .eq. nodes(i)-1) then
          lb(i) = -tab
          ub(i) = dsize(i)+tab-1
          offset(i) = 0
       else if(coord(i) .eq. 0) then
          lb(i) = -tab
          ub(i) = dsize(i)-2
          offset(i) = 0
       else if(coord(i) .eq. nodes(i)-1) then
          lb(i) = 0
          ub(i) = dsize(i)+tab-1
          offset(i) = (dsize(i)-1)*coord(i)+tab
       else
          lb(i) = 0
          ub(i) = dsize(i)-2
          offset(i) = (dsize(i)-1)*coord(i)+tab
       endif
       local_data_size(i) = ub(i)-lb(i)+1
    enddo

    call h5screate_simple_f(rank, local_data_size, memspace, status)
    call h5dget_space_f(d_id, filespace, status)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, local_data_size, status)

    call h5pcreate_f(H5P_DATASET_XFER_F, p_id, status)
    call h5pset_dxpl_mpio_f(p_id, H5FD_MPIO_COLLECTIVE_F, status)

    call h5dwrite_f(d_id, H5T_NATIVE_DOUBLE, data(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), all_data_size, istat, &
         &file_space_id=filespace, mem_space_id=memspace, xfer_prp=p_id)

    call h5sclose_f(filespace, status)
    call h5sclose_f(memspace, status)
    call h5dclose_f(d_id, status)
    call h5pclose_f(p_id, status)
!
    return
  end subroutine wrtfield3d


  end module hdf5io
