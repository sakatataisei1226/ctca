#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  module hdf
!
!   ____________________________________________________________
!
!                    M O D U L E   H   D   F
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .         common declaration of variables for HDF          .
!   ............................................................
!
!-------------------- common declaration of variables for HDF
  use paramt
  implicit none

    integer(kind=4) :: DFACC_READ, DFACC_WRITE, DFACC_RW, DFACC_CREATE
    integer(kind=4) :: DF_INT4, DF_REAL8
    integer(kind=4) :: irank

    parameter(DFACC_READ=1, DFACC_WRITE=2, DFACC_RW=3, DFACC_CREATE=4)
    parameter(DF_INT4=24, DF_REAL8=6)
    parameter(irank=OH_DIMENSION)

!    integer(kind=4) :: sd_id, sds_id, status
!    integer(kind=4),external :: sfstart, sfcreate, sfwdata, sfendacc, sfend
!    integer(kind=4) :: datatype, attributus
!    integer(kind=4) :: dim(irank),start(irank),edges(irank),stride(irank)
!    character(len=7) :: sdname

    integer(kind=HID_T) :: id_sd(2,400)
    integer(kind=4) :: istats(2,400)
    integer(kind=8) :: idimfd(irank,FNR), idimft(irank)
!    integer(kind=4) :: istafd(irank,FNR), istaft(irank)
!    integer(kind=4) :: istrfd(irank,FNR), istrft(irank)
!    integer(kind=4) :: iendfd(irank,FNR), iendft(irank)

    common /idnum/  id_sd
    common /hdfcom/ idimfd,idimft
!    common /hdfcom/ idimfd,idimft, &
!   &                istafd,istaft, &
!   &                istrfd,istrft, &
!   &                iendfd,iendft


  end module

