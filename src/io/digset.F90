#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine digset
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   D I G S E T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   . this subroutine sets parameters required for HDF output. .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  use hdf

  implicit none


!-------------------- digset: idim
      idimfd(1,FEB) = nxsd + 1
      idimfd(2,FEB) = nysd + 1
      idimfd(3,FEB) = nzsd + 1
!
      idimfd(1,FAJ) = nxsd + 1
      idimfd(2,FAJ) = nysd + 1
      idimfd(3,FAJ) = nzsd + 1
!
      idimfd(1,FRH) = nxsd + 1
      idimfd(2,FRH) = nysd + 1
      idimfd(3,FRH) = nzsd + 1
!
      idimfd(1,FPH) = nxsd + 1
      idimfd(2,FPH) = nysd + 1
      idimfd(3,FPH) = nzsd + 1
!
      idimfd(1,FRD) = nxsd + 1
      idimfd(2,FRD) = nysd + 1
      idimfd(3,FRD) = nzsd + 1
!
      idimft(1) = lwfft
      idimft(2) = ldfft
      idimft(3) = lhfft


!-------------------- digset: istart
!      istafd(1:3,FEB:FRD) = 0
!      istaft(1:3) = 0


!-------------------- digset: istride
!      istrfd(1:3,FEB:FRD) = 1
!      istrft(1:3) = 1


!-------------------- digset: iend
!      iendfd(:,:) = idimfd(:,:)
!      iendft(:)  = idimft(:)


!-------------------- 
      h5fcount = 0
      h5jcount = 0


  return
  end subroutine
