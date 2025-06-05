#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine create_nborps(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   N B O R P S
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
  implicit none
!
  integer(kind=4) :: is, ps


!--------------------
      if(sdid(ps).eq.-1) return


!--------------------
      do is=1,nspec
      nborps(14,is,ps) = sdid(ps)
!! 
      if(mod(nborps(14,is,ps),nodes(1)).ne.0) then
        nborps(13,is,ps) = nborps(14,is,ps) - 1
      else if(npbnd(1,is).eq.0) then
        nborps(13,is,ps) = nborps(14,is,ps) + nodes(1) - 1
      else if(npbnd(1,is).eq.1) then
        nborps(13,is,ps) = nborps(14,is,ps)
      else
        nborps(13,is,ps) = -1
      end if
      if(mod(nborps(14,is,ps),nodes(1)).ne.nodes(1)-1) then
        nborps(15,is,ps) = nborps(14,is,ps) + 1
      else if(npbnd(1,is).eq.0) then
        nborps(15,is,ps) = nborps(14,is,ps) - nodes(1) + 1
      else if(npbnd(1,is).eq.1) then
        nborps(15,is,ps) = nborps(14,is,ps)
      else
        nborps(15,is,ps) = -1
      end if
!
      if(mod(nborps(14,is,ps),nodes(1)*nodes(2)).ge.nodes(1)) then
        nborps(11,is,ps) = nborps(14,is,ps) - nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(11,is,ps) = nborps(14,is,ps) + nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(11,is,ps) = nborps(14,is,ps)
      else
        nborps(11,is,ps) = -1
      end if
      if(mod(nborps(14,is,ps),nodes(1)*nodes(2)).lt.nodes(1)*(nodes(2)-1)) then
        nborps(17,is,ps) = nborps(14,is,ps) + nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(17,is,ps) = nborps(14,is,ps) - nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(17,is,ps) = nborps(14,is,ps)
      else
        nborps(17,is,ps) = -1
      end if
!
      if(nborps(14,is,ps)       .ge.nodes(1)*nodes(2)) then
        nborps( 5,is,ps) = nborps(14,is,ps) - nodes(1)*nodes(2)
      else if(npbnd(3,is).eq.0) then
        nborps( 5,is,ps) = nborps(14,is,ps) + nodes(1)*nodes(2)*(nodes(3) - 1)
      else if(npbnd(3,is).eq.1) then
        nborps( 5,is,ps) = nborps(14,is,ps)
      else
        nborps( 5,is,ps) = -1
      end if
      if(nborps(14,is,ps)       .lt.nodes(1)*nodes(2)*(nodes(3)-1)) then
        nborps(23,is,ps) = nborps(14,is,ps) + nodes(1)*nodes(2)
      else if(npbnd(3,is).eq.0) then
        nborps(23,is,ps) = nborps(14,is,ps) - nodes(1)*nodes(2)*(nodes(3) - 1)
      else if(npbnd(3,is).eq.1) then
        nborps(23,is,ps) = nborps(14,is,ps)
      else
        nborps(23,is,ps) = -1
      end if
!!
      if(nborps(13,is,ps).eq.-1) then
        nborps(10,is,ps) = -1
      else if(mod(nborps(13,is,ps),nodes(1)*nodes(2)).ge.nodes(1)) then
        nborps(10,is,ps) = nborps(13,is,ps) - nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(10,is,ps) = nborps(13,is,ps) + nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(10,is,ps) = nborps(13,is,ps)
      else
        nborps(10,is,ps) = -1
      end if
      if(nborps(13,is,ps).eq.-1) then
        nborps(16,is,ps) = -1
      else if(mod(nborps(13,is,ps),nodes(1)*nodes(2)).lt.nodes(1)*(nodes(2)-1)) then
        nborps(16,is,ps) = nborps(13,is,ps) + nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(16,is,ps) = nborps(13,is,ps) - nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(16,is,ps) = nborps(13,is,ps)
      else
        nborps(16,is,ps) = -1
      end if
!
      if(nborps(15,is,ps).eq.-1) then
        nborps(12,is,ps) = -1
      else if(mod(nborps(15,is,ps),nodes(1)*nodes(2)).ge.nodes(1)) then
        nborps(12,is,ps) = nborps(15,is,ps) - nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(12,is,ps) = nborps(15,is,ps) + nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(12,is,ps) = nborps(15,is,ps)
      else
        nborps(12,is,ps) = -1
      end if
      if(nborps(15,is,ps).eq.-1) then
        nborps(18,is,ps) = -1
      else if(mod(nborps(15,is,ps),nodes(1)*nodes(2)).lt.nodes(1)*(nodes(2)-1)) then
        nborps(18,is,ps) = nborps(15,is,ps) + nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(18,is,ps) = nborps(15,is,ps) - nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(18,is,ps) = nborps(15,is,ps)
      else
        nborps(18,is,ps) = -1
      end if
!
      if(nborps( 5,is,ps).eq.-1) then
        nborps( 4,is,ps) = -1
      else if(mod(nborps( 5,is,ps),nodes(1)).ne.0) then
        nborps( 4,is,ps) = nborps( 5,is,ps) - 1
      else if(npbnd(1,is).eq.0) then
        nborps( 4,is,ps) = nborps( 5,is,ps) + nodes(1) - 1
      else if(npbnd(1,is).eq.1) then
        nborps( 4,is,ps) = nborps( 5,is,ps)
      else
        nborps( 4,is,ps) = -1
      end if
      if(nborps( 5,is,ps).eq.-1) then
        nborps( 6,is,ps) = -1
      else if(mod(nborps( 5,is,ps),nodes(1)).ne.nodes(1)-1) then
        nborps( 6,is,ps) = nborps( 5,is,ps) + 1
      else if(npbnd(1,is).eq.0) then
        nborps( 6,is,ps) = nborps( 5,is,ps) - nodes(1) + 1
      else if(npbnd(1,is).eq.1) then
        nborps( 6,is,ps) = nborps( 5,is,ps)
      else
        nborps( 6,is,ps) = -1
      end if
!
      if(nborps( 5,is,ps).eq.-1) then
        nborps( 2,is,ps) = -1
      else if(mod(nborps( 5,is,ps),nodes(1)*nodes(2)).ge.nodes(1)) then
        nborps( 2,is,ps) = nborps( 5,is,ps) - nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps( 2,is,ps) = nborps( 5,is,ps) + nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps( 2,is,ps) = nborps( 5,is,ps)
      else
        nborps( 2,is,ps) = -1
      end if
      if(nborps( 5,is,ps).eq.-1) then
        nborps( 8,is,ps) = -1
      else if(mod(nborps( 5,is,ps),nodes(1)*nodes(2)).lt.nodes(1)*(nodes(2)-1)) then
        nborps( 8,is,ps) = nborps( 5,is,ps) + nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps( 8,is,ps) = nborps( 5,is,ps) - nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps( 8,is,ps) = nborps( 5,is,ps)
      else
        nborps( 8,is,ps) = -1
      end if
!
      if(nborps(23,is,ps).eq.-1) then
        nborps(22,is,ps) = -1
      else if(mod(nborps(23,is,ps),nodes(1)).ne.0) then
        nborps(22,is,ps) = nborps(23,is,ps) - 1
      else if(npbnd(1,is).eq.0) then
        nborps(22,is,ps) = nborps(23,is,ps) + nodes(1) - 1
      else if(npbnd(1,is).eq.1) then
        nborps(22,is,ps) = nborps(23,is,ps)
      else
        nborps(22,is,ps) = -1
      end if
      if(nborps(23,is,ps).eq.-1) then
        nborps(24,is,ps) = -1
      else if(mod(nborps(23,is,ps),nodes(1)).ne.nodes(1)-1) then
        nborps(24,is,ps) = nborps(23,is,ps) + 1
      else if(npbnd(1,is).eq.0) then
        nborps(24,is,ps) = nborps(23,is,ps) - nodes(1) + 1
      else if(npbnd(1,is).eq.1) then
        nborps(24,is,ps) = nborps(23,is,ps)
      else
        nborps(24,is,ps) = -1
      end if
!
      if(nborps(23,is,ps).eq.-1) then
        nborps(20,is,ps) = -1
      else if(mod(nborps(23,is,ps),nodes(1)*nodes(2)).ge.nodes(1)) then
        nborps(20,is,ps) = nborps(23,is,ps) - nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(20,is,ps) = nborps(23,is,ps) + nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(20,is,ps) = nborps(23,is,ps)
      else
        nborps(20,is,ps) = -1
      end if
      if(nborps(23,is,ps).eq.-1) then
        nborps(26,is,ps) = -1
      else if(mod(nborps(23,is,ps),nodes(1)*nodes(2)).lt.nodes(1)*(nodes(2)-1)) then
        nborps(26,is,ps) = nborps(23,is,ps) + nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(26,is,ps) = nborps(23,is,ps) - nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(26,is,ps) = nborps(23,is,ps)
      else
        nborps(26,is,ps) = -1
      end if
!!
      if(nborps( 4,is,ps).eq.-1) then
        nborps( 1,is,ps) = -1
      else if(mod(nborps( 4,is,ps),nodes(1)*nodes(2)).ge.nodes(1)) then
        nborps( 1,is,ps) = nborps( 4,is,ps) - nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps( 1,is,ps) = nborps( 4,is,ps) + nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps( 1,is,ps) = nborps( 4,is,ps)
      else
        nborps( 1,is,ps) = -1
      end if
      if(nborps( 4,is,ps).eq.-1) then
        nborps( 7,is,ps) = -1
      else if(mod(nborps( 4,is,ps),nodes(1)*nodes(2)).lt.nodes(1)*(nodes(2)-1)) then
        nborps( 7,is,ps) = nborps( 4,is,ps) + nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps( 7,is,ps) = nborps( 4,is,ps) - nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps( 7,is,ps) = nborps( 4,is,ps)
      else
        nborps( 7,is,ps) = -1
      end if
!
      if(nborps( 6,is,ps).eq.-1) then
        nborps( 3,is,ps) = -1
      else if(mod(nborps( 6,is,ps),nodes(1)*nodes(2)).ge.nodes(1)) then
        nborps( 3,is,ps) = nborps( 6,is,ps) - nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps( 3,is,ps) = nborps( 6,is,ps) + nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps( 3,is,ps) = nborps( 6,is,ps)
      else
        nborps( 3,is,ps) = -1
      end if
      if(nborps( 6,is,ps).eq.-1) then
        nborps( 9,is,ps) = -1
      else if(mod(nborps( 6,is,ps),nodes(1)*nodes(2)).lt.nodes(1)*(nodes(2)-1)) then
        nborps( 9,is,ps) = nborps( 6,is,ps) + nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps( 9,is,ps) = nborps( 6,is,ps) - nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps( 9,is,ps) = nborps( 6,is,ps)
      else
        nborps( 9,is,ps) = -1
      end if
!
      if(nborps(22,is,ps).eq.-1) then
        nborps(19,is,ps) = -1
      else if(mod(nborps(22,is,ps),nodes(1)*nodes(2)).ge.nodes(1)) then
        nborps(19,is,ps) = nborps(22,is,ps) - nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(19,is,ps) = nborps(22,is,ps) + nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(19,is,ps) = nborps(22,is,ps)
      else
        nborps(19,is,ps) = -1
      end if
      if(nborps(22,is,ps).eq.-1) then
        nborps(25,is,ps) = -1
      else if(mod(nborps(22,is,ps),nodes(1)*nodes(2)).lt.nodes(1)*(nodes(2)-1)) then
        nborps(25,is,ps) = nborps(22,is,ps) + nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(25,is,ps) = nborps(22,is,ps) - nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(25,is,ps) = nborps(22,is,ps)
      else
        nborps(25,is,ps) = -1
      end if
!
      if(nborps(24,is,ps).eq.-1) then
        nborps(21,is,ps) = -1
      else if(mod(nborps(24,is,ps),nodes(1)*nodes(2)).ge.nodes(1)) then
        nborps(21,is,ps) = nborps(24,is,ps) - nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(21,is,ps) = nborps(24,is,ps) + nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(21,is,ps) = nborps(24,is,ps)
      else
        nborps(21,is,ps) = -1
      end if
      if(nborps(24,is,ps).eq.-1) then
        nborps(27,is,ps) = -1
      else if(mod(nborps(24,is,ps),nodes(1)*nodes(2)).lt.nodes(1)*(nodes(2)-1)) then
        nborps(27,is,ps) = nborps(24,is,ps) + nodes(1)
      else if(npbnd(2,is).eq.0) then
        nborps(27,is,ps) = nborps(24,is,ps) - nodes(1)*(nodes(2) - 1)
      else if(npbnd(2,is).eq.1) then
        nborps(27,is,ps) = nborps(24,is,ps)
      else
        nborps(27,is,ps) = -1
      end if
      end do


  return
  end subroutine
