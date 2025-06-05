  module interf
!
    interface
!      subroutine arange(ii,aa,bb,nam,id)
!        implicit none
!        integer(kind=4) :: ii, id
!        integer(kind=8) :: nam
!        real(kind=8) :: aa(nam),bb(nam)
!      end subroutine arange
!
!      subroutine arymod(aa,bb,in,nam,is,id)
!        implicit none
!        integer(kind=8) :: in, nam
!        integer(kind=4) :: is, id
!        real(kind=8) :: aa(in),bb(in)
!      end subroutine arymod
!
!      subroutine mtsatn(m,dt,nstep,w,attn,damp)
!        implicit none
!        integer(kind=4) :: m, nstep, i
!        real(kind=8) :: dt, w, attn, damp
!      end subroutine mtsatn
!
!      subroutine qsort(x,n)
!        implicit none
!        integer(kind=8) :: n
!        real(kind=8) :: x(n)
!      end subroutine qsort
!
!      subroutine rexch(ii,n,i1,i2,i3,i4,r1,r2,r3,r4, &
!     &                      d1,d2,d3,d4,d5,d6,d7,d8)
!        implicit none
!        integer(kind=4),intent(inout) :: ii
!        integer(kind=8),intent(in) :: n
!        integer(kind=4),optional,intent(inout) :: i1(n),i2(n),i3(n),i4(n)
!        real(kind=4),optional,intent(inout) :: r1(n),r2(n),r3(n),r4(n)
!        real(kind=8),optional,intent(inout) :: d1(n),d2(n),d3(n),d4(n)
!        real(kind=8),optional,intent(inout) :: d5(n),d6(n),d7(n),d8(n)
!      end subroutine rexch
!
!      subroutine ringds(ii,ar,na,vt,bbb)
!        implicit none
!        integer(kind=4) :: ii
!        integer(kind=8) :: na
!        real(kind=8) :: ar(na), vt, bbb
!      end subroutine ringds
!
!      subroutine vflux(ii,ar,na,v0,vt,isl,icon)
!        implicit none
!        integer(kind=4) :: ii, isl, icon
!        integer(kind=8) :: na
!        real(kind=8) :: ar(na), v0, vt
!      end subroutine vflux
!
!      subroutine vinit(ii,ar,na,v0,vt,xlos)
!        implicit none
!        integer(kind=4) :: ii
!        integer(kind=8) :: na
!        real(kind=8) :: ar(na), v0, vt, xlos
!      end subroutine vinit
    end interface

  end module interf
