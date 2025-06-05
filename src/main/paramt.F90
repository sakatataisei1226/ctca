module paramt
!
!   ____________________________________________________________
!
!                    M O D U L E   P A R A M T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   this file defines the scale of each array.             .
!   .   first of all, you should determin each parameter       .
!   .   following this.                                        .
!   ............................................................
!
!
!       ix >= nx + 4
!       iy >= ny + 4
!       iz >= nz + 4
!       ixy = max(ix,iy,iz)
!       in   >= np(1) + np(2) + np(3) + ...
!       inr  >= npr(1)+ npr(2)+ npr(3)+ ...
!       iw : buffer size for data output
!       ispec : number of particle species
!       ilbsum >= ilb(1) + ilb(2) + ilb(3) + ...
!           ilb(1): stuck size for particle 1
!       ivx >= nvx
!       ivy >= nvy
!       ivz >= nvz
!       ixw >= ncpmx=(nxw2-nxw1+2)*(nyw2-nyw1+2)*(nzw2-nzw1+2)
!
!--------------------
    use HDF5
    implicit none
!
    integer(kind=4) :: ix, iy, iz, ixy, ispec, ivx, ivy, ivz
    integer(kind=4) :: ixw, inpc, ingap, isj, injs, inepl, inwave, intch
    integer(kind=4) :: nblimit, lpara
    integer(kind=8) :: inr, i_work1_num
    integer(kind=4) :: OH_IPBUF_SIZE
    integer(kind=4) :: ncomreq

    integer, parameter :: MAXFRAC = 10
    integer, parameter :: FNR = 7, FEB = 1, FDB = 2, FAJ = 3, FRH = 4, FPH = 5, FRD = 6, FJD = 7
    integer, parameter :: BNR = 2
    integer, parameter :: CNR = 6, CEB = 1, CAJ = 2, CRH = 3, CRD = 4, CJD = 5, CPH = 6
    integer, parameter :: EX = 1, EY = 2, EZ = 3, BX = 4, BY = 5, BZ = 6
    integer, parameter :: CX = 1, CY = 2, CZ = 3, CP = 4
    integer, parameter :: JX = 1, JY = 2, JZ = 3, JS = 4, TRH = 7, TJX = 7, TJY = 8, TJZ = 9

    parameter(ix=4)
    parameter(iy=4)
    parameter(iz=4)
    parameter(ixy=4)
    parameter(inr=2097152)
    parameter(ispec=10)
    parameter(ivx=32)
    parameter(ivy=32)
    parameter(ivz=32)
!    parameter( ixw    = 11*11*11+5 )
    parameter(ixw=5000)
    parameter(inpc=15)
    parameter(injs=4)
    parameter(inwave=4)
    parameter(intch=4)
    parameter(ingap=4)
    parameter(isj=100*6*3)
!    parameter( ipoint = 1 )
    parameter(inepl=16)
!    parameter( lvec   = 1 )
    parameter(nblimit=16)
    parameter(i_work1_num=10_8)
    parameter(lpara=1)
    parameter(OH_IPBUF_SIZE=16*1024)
    parameter(ncomreq=10000)

    type counter
        sequence
        integer(kind=8) :: globalp(ispec, 2)
        integer(kind=4) :: influx(inpc, ispec), outflux(inpc, ispec)
        integer(kind=4) :: infhist(0:inpc, inpc, ispec)
        integer(kind=4) :: nesc((ispec + 1)/2*2)
        real(kind=8) :: chgacm(2, inpc)
    end type

    type connect
        sequence
        integer(kind=4) :: from, to
        real(kind=8) :: val
    end type

    type cylindrical
        sequence
        integer(kind=4) :: align, pad
        real(kind=8) :: radius, edge(2), axis(2)
    end type

    type spherical
        sequence
        real(kind=8) :: radius, center(3)
    end type

    type wire
        sequence
        integer(kind=4) :: align, pad
        real(kind=8) :: hlength, rradius, eradius, origin(3)
    end type

    type particle3
        sequence
        real(kind=8) :: x, y, z
    end type

end module
