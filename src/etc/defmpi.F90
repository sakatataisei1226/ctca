#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine defmpi
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   D E F M P I
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  use hdf
#define MCW local_comm
#define MSS MPI_STATUS_SIZE
#define MAK MPI_ADDRESS_KIND
  implicit none
  integer(kind=4) :: i,j, m
  integer(kind=4) :: ncomp(2)=(/3, 1/), ebyte
  integer(kind=4) :: mpierr
  interface
    subroutine sumcount(in, inout, len, dtype)
      use paramt
      integer :: len
      integer :: dtype
      type(counter) :: in(len), inout(len)
    end subroutine sumcount
  end interface


!-------------------- 
      do m=1,2
        ebyte = 8*ncomp(m)
!
        do j=0,1
        do i=0,1
          call MPI_Type_Create_Hvector((nysd+j), (nxsd+i)*ncomp(m), &
         &                             int(lwslc*ebyte,MAK), &
         &                             MPI_REAL8, mptype_rs(i,j,m,1), mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int((/0,0,lwslc*ldslc/)*ebyte,MAK), &
         &                            (/MPI_LB,mptype_rs(i,j,m,1),MPI_UB/), &
         &                            mptype_rs(i,j,m,1), mpierr)
          call MPI_Type_Commit(mptype_rs(i,j,m,1), mpierr)
        end do
        end do
!
        do j=0,1
        do i=0,1
          call MPI_Type_Create_Hvector((nysd+j), (nxsd+i)*ncomp(m), &
         &                             int(lxsd(FDB)*ebyte,MAK), &
         &                             MPI_REAL8, mptype_sc(i,j,m,1), mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int((/0,0,lxsd(FDB)*lysd(FDB)/)*ebyte,MAK), &
         &                            (/MPI_LB,mptype_sc(i,j,m,1),MPI_UB/), &
         &                            mptype_sc(i,j,m,1), mpierr)
          call MPI_Type_Commit(mptype_sc(i,j,m,1), mpierr)
        end do
        end do
!
        do i=-1,0
          call MPI_Type_Create_Hvector((nysd+3), (nxsd+3)*ncomp(m), &
         &                             int(lxsd(FDB)*ebyte,MAK), &
         &                             MPI_REAL8, mptype_rc(i,m,1), mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int( &
         &                            (-1-lxsd(FDB)+lxsd(FDB)*lysd(FDB)*(/i,i,i+1/))*ebyte, &
         &                            MAK), &
         &                            (/MPI_LB,mptype_rc(i,m,1),MPI_UB/), &
         &                            mptype_rc(i,m,1), mpierr)
          call MPI_Type_Commit(mptype_rc(i,m,1), mpierr)
        end do
!
        do i=-1,0
          call MPI_Type_Create_Hvector((nysd+3), (nxsd+3)*ncomp(m), &
         &                             int(lwslc*ebyte,MAK), &
         &                             MPI_REAL8, mptype_ss(i,m,1), mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int((-1-lwslc+lwslc*ldslc*(/i,i,i+1/))*ebyte, &
         &                                MAK), &
         &                            (/MPI_LB,mptype_ss(i,m,1),MPI_UB/), &
         &                            mptype_ss(i,m,1), mpierr)
          call MPI_Type_Commit(mptype_ss(i,m,1), mpierr)
        end do
!
        do i=0,1
          call MPI_Type_Create_Hvector((ny+1), (nxslc+1)*ncomp(m), &
         &                             int(lwslc*ebyte,MAK), &
         &                             MPI_REAL8, mptype_aa(i,0,m,1), mpierr)
          call MPI_Type_Create_Hvector((nzslc+1), 1, &
         &                             int(lwslc*ldslc*ebyte,MAK), &
         &                             mptype_aa(i,0,m,1), mptype_aa(i,0,m,1), &
         &                             mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int((/0,0,(nxslc+i)/)*ebyte,MAK), &
         &                            (/MPI_LB,mptype_aa(i,0,m,1),MPI_UB/), &
         &                            mptype_aa(i,0,m,1), mpierr)
          call MPI_Type_Commit(mptype_aa(i,0,m,1), mpierr)
!
          call MPI_Type_Create_Hvector((ny+1), (nzslc+1)*ncomp(m), &
         &                             int(lwslc*ebyte,MAK), &
         &                             MPI_REAL8, mptype_aa(i,1,m,1), mpierr)
          call MPI_Type_Create_Hvector((nxslc+1), 1, &
         &                             int(lwslc*ldslc*ebyte,MAK), &
         &                             mptype_aa(i,1,m,1), mptype_aa(i,1,m,1), &
         &                             mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int((/0,0,(nzslc+i)/)*ebyte,MAK), &
         &                            (/MPI_LB,mptype_aa(i,1,m,1),MPI_UB/), &
         &                            mptype_aa(i,1,m,1), mpierr)
          call MPI_Type_Commit(mptype_aa(i,1,m,1), mpierr)
        end do
      end do


!-------------------- 
      do m=1,2
        ebyte = 8*ncomp(m)
!
        do j=0,1
        do i=0,1
          call MPI_Type_Create_Hvector((nzsd+j), (nxsd+i)*ncomp(m), &
         &                             int(lwslc*ldslc*ebyte,MAK), &
         &                             MPI_REAL8, mptype_rs(i,j,m,2), mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int((/0,0,lwslc/)*ebyte,MAK), &
         &                            (/MPI_LB,mptype_rs(i,j,m,2),MPI_UB/), &
         &                            mptype_rs(i,j,m,2), mpierr)
          call MPI_Type_Commit(mptype_rs(i,j,m,2), mpierr)
        end do
        end do
!
        do j=0,1
        do i=0,1
          call MPI_Type_Create_Hvector((nzsd+j), (nxsd+i)*ncomp(m), &
         &                             int(lxsd(FDB)*lysd(FDB)*ebyte,MAK), &
         &                             MPI_REAL8, mptype_sc(i,j,m,2), mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int((/0,0,lxsd(FDB)/)*ebyte,MAK), &
         &                            (/MPI_LB,mptype_sc(i,j,m,2),MPI_UB/), &
         &                            mptype_sc(i,j,m,2), mpierr)
          call MPI_Type_Commit(mptype_sc(i,j,m,2), mpierr)
        end do
        end do
!
        do i=-1,0
          call MPI_Type_Create_Hvector((nzsd+3), (nxsd+3)*ncomp(m), &
         &                             int(lxsd(FDB)*lysd(FDB)*ebyte,MAK), &
         &                             MPI_REAL8, mptype_rc(i,m,2), mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int( &
         &                            (-1-lxsd(FDB)*lysd(FDB)+lxsd(FDB)*(/i,i,i+1/))*ebyte, &
         &                            MAK), &
         &                            (/MPI_LB,mptype_rc(i,m,2),MPI_UB/), &
         &                            mptype_rc(i,m,2), mpierr)
          call MPI_Type_Commit(mptype_rc(i,m,2), mpierr)
        end do
!
        do i=-1,0
          call MPI_Type_Create_Hvector((nzsd+3), (nxsd+3)*ncomp(m), &
         &                             int(lwslc*ldslc*ebyte,MAK), &
         &                             MPI_REAL8, mptype_ss(i,m,2), mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int((-1-lwslc*ldslc+lwslc*(/i,i,i+1/))*ebyte, &
         &                                MAK), &
         &                            (/MPI_LB,mptype_ss(i,m,2),MPI_UB/), &
         &                            mptype_ss(i,m,2), mpierr)
          call MPI_Type_Commit(mptype_ss(i,m,2), mpierr)
        end do
!
        do i=0,1
          call MPI_Type_Create_Hvector((nyslc+1), (nxslc+1)*ncomp(m), &
         &                             int(lwslc*ebyte,MAK), &
         &                             MPI_REAL8, mptype_aa(i,0,m,2), mpierr)
          call MPI_Type_Create_Hvector((nz+1), 1, &
         &                             int(lwslc*ldslc*ebyte,MAK), &
         &                             mptype_aa(i,0,m,2), mptype_aa(i,0,m,2), &
         &                             mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int((/0,0,(nxslc+i)/)*ebyte,MAK), &
         &                            (/MPI_LB,mptype_aa(i,0,m,2),MPI_UB/), &
         &                            mptype_aa(i,0,m,2), mpierr)
          call MPI_Type_Commit(mptype_aa(i,0,m,2), mpierr)
!
          call MPI_Type_Create_Hvector((nxslc+1), (nyslc+1)*ncomp(m), &
         &                             int(lwslc*ebyte,MAK), &
         &                             MPI_REAL8, mptype_aa(i,1,m,2), mpierr)
          call MPI_Type_Create_Hvector((nz+1), 1, &
         &                             int(lwslc*ldslc*ebyte,MAK), &
         &                             mptype_aa(i,1,m,2), mptype_aa(i,1,m,2), &
         &                             mpierr)
          call MPI_Type_Create_Struct(3, (/1,1,1/), &
         &                            int((/0,0,(nyslc+i)/)*ebyte,MAK), &
         &                            (/MPI_LB,mptype_aa(i,1,m,2),MPI_UB/), &
         &                            mptype_aa(i,1,m,2), mpierr)
          call MPI_Type_Commit(mptype_aa(i,1,m,2), mpierr)
        end do
      end do


!-------------------- 
      do m=1,2
        ebyte = 8*ncomp(m)
        call MPI_Type_Contiguous(lwslc*ldslc*ncomp(m), &
       &                         MPI_REAL8, mptype_xy(m), mpierr)
        call MPI_Type_Commit(mptype_xy(m), mpierr)
      end do


!-------------------- 
      do m=1,2
        ebyte = 8*ncomp(m)
        call MPI_Type_Vector(ldslc, ncomp(m), lwslc*ncomp(m), &
       &                     MPI_REAL8, mptype_yz(m), mpierr)
        call MPI_Type_Create_Hvector(lhslc, 1, int(lwslc*ldslc*ebyte,MAK), &
       &                     mptype_yz(m), mptype_yz(m), mpierr)
        call MPI_Type_Create_Struct(2, (/1,1/), int((/0,1/)*ebyte,MAK), &
       &                            (/mptype_yz(m),MPI_UB/), &
       &                            mptype_yz(m), mpierr)
        call MPI_Type_Commit(mptype_yz(m), mpierr)
      end do


!-------------------- 
      do m=1,2
        ebyte = 8*ncomp(m)
        call MPI_Type_Vector(lhslc, lwslc*ncomp(m), lwslc*ldslc*ncomp(m), &
       &                     MPI_REAL8, mptype_xz(m), mpierr)
        call MPI_Type_Create_Struct(2, (/1,1/), int((/0,lwslc/)*ebyte,MAK), &
       &                            (/mptype_xz(m),MPI_UB/), &
       &                            mptype_xz(m), mpierr)
        call MPI_Type_Commit(mptype_xz(m), mpierr)
      end do


!-------------------- 
      call MPI_Comm_Split(MCW, int(myid/snode), mod(myid,snode), &
     &                    subcomm, mpierr)


!-------------------- MPI datatype & operation definitions for counter
      call MPI_Type_Create_Struct(6, &
     &  (/ ispec*2, inpc*ispec, inpc*ispec, &
     &     (inpc+1)*inpc*ispec, (ispec+1)/2*2, 2*inpc /), &
     &  int((/ 0, ispec*2*8, ispec*2*8 + inpc*ispec*4, &
     &         ispec*2*8 + 2*inpc*ispec*4, &
     &         ispec*2*8 + (inpc + 3)*inpc*ispec*4, &
     &         ispec*2*8 + ((inpc + 3)*inpc*ispec + (ispec+1)/2*2)*4 /), &
     &      MAK), &
     &  (/ MPI_INTEGER8, MPI_INTEGER, MPI_INTEGER, &
     &     MPI_INTEGER, MPI_INTEGER, MPI_REAL8 /), &
     &  mpi_type_count, mpierr)
      call MPI_Type_Commit(mpi_type_count, mpierr)
!
      call MPI_Op_Create(sumcount, .true., mpi_sum_count, mpierr)


  return
  end subroutine defmpi


  subroutine sumcount(in, inout, len, dtype)
    use paramt
    implicit none
    integer :: len
    integer :: dtype
    type(counter) :: in(len), inout(len)
    integer :: i

!------------------------------ 
      do i=1,len
        inout(i)%globalp(:,:) = inout(i)%globalp(:,:) + in(i)%globalp(:,:)
        inout(i)%influx(:,:) = inout(i)%influx(:,:) + in(i)%influx(:,:)
        inout(i)%outflux(:,:) = inout(i)%outflux(:,:) + in(i)%outflux(:,:)
        inout(i)%infhist(:,:,:) = inout(i)%infhist(:,:,:) + in(i)%infhist(:,:,:)
        inout(i)%nesc(:) = inout(i)%nesc(:) + in(i)%nesc(:)
        inout(i)%chgacm(:,:) = inout(i)%chgacm(:,:) + in(i)%chgacm(:,:)
      end do


  return
  end subroutine sumcount
