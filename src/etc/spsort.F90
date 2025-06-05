#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine spsort(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   S P S O R T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .     this subroutine updates particle position.           .
!   .     the leap-frog scheme is used.                        .
!   ............................................................
!

!-------------------- parameters and variables
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=8) :: ns,ne
  integer(kind=4) :: is
  integer(kind=4) :: ps


!============================== species loop
      ne = pbase(ps)
 ISL: do is=1,nspec
        ns = ne + 1
        ne = ne + totalp(is,ps)
        call particle_sort(pbuf(ns:ne),totalp(is,ps))
      end do ISL


  return
  end subroutine

!
  subroutine specsort
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   S P S O R T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .     this subroutine updates particle position.           .
!   .     the leap-frog scheme is used.                        .
!   ............................................................
!

!-------------------- parameters and variables
  use oh_type
  use paramt
  use allcom
  implicit none


!============================== species loop
      call particle_spec(pbuf(1:sum(totalp(1:nspec,1))), &
     &  sum(totalp(1:nspec,1)))


  return
  end subroutine
