#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
subroutine charge(ps, func)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   C H A R G E
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a charge distribution to      .
!   .      grid points from particle locations by the          .
!   .      "area-sharing" or "area-weighting" scheme.          .
!   ............................................................
!
!-------------------- parameter and common block
    use oh_type
    use paramt
    use allcom
    use m_ohhinfo
    implicit none

    type t_InterpolationInfo
        double precision :: v1, v2, v3, v4, v5, v6, v7, v8
        integer :: i, j, k, i1, j1, k1
    end type

!
    integer(kind=4) :: func ! 2 if require data for output otherwise 0
    integer(kind=4) :: ps

    integer(kind=8) :: m
    integer(kind=8) :: ns, ne, nee
    integer :: is, nis

    type(t_InterpolationInfo) :: iw

    call ohhinfo_update(sdoms(:, :, sdid(ps) + 1))

!-------------------- for three dimensional system
!     --------------- zero clear rho
    if (func .eq. 0 .or. func .eq. 2) then
        rho(:, :, :, :, ps) = 0.0d0
    end if
!     --------------- zero clear rhodg
    if (func .eq. 1 .or. func .eq. 2) then
        rhodg(:, :, :, :, ps) = 0.0d0
    end if

!-------------------- sum up charge component
    nee = pbase(ps)
    ISL: do is = 1, nspec
        ns = nee + 1
        ne = nee + totalp(is, ps)
        nee = ne

        if (func .eq. 0) then
            do m = ns, ne
                if (pbuf(m)%nid .eq. -1) cycle

                iw = get_interpolation_weight(pbuf(m))

                if (pbuf(m)%preside >= 0) then
                    call distribute_charge(rho(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), q(is), iw)

                else if (pbuf(m)%preside == OH_PCL_TO_BE_ACCUMULATED) then
                    call distribute_charge(rhobk(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), q(is), iw)
                    call distribute_charge(rhobksp(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps, is), q(is), iw)
                    nphgram(pbuf(m)%nid + 1, is, ps) = nphgram(pbuf(m)%nid + 1, is, ps) - 1
                    pbuf(m)%nid = -1
                    pbuf(m)%preside = 0

                else if (pbuf(m)%preside == OH_PCL_TO_BE_ACCUMULATED_AS_ANTIPCL) then
                    call distribute_charge(rhobk(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), -q(is), iw)
                    call distribute_charge(rhobksp(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps, is), -q(is), iw)
                    nphgram(pbuf(m)%nid + 1, is, ps) = nphgram(pbuf(m)%nid + 1, is, ps) - 1
                    pbuf(m)%nid = -1
                    pbuf(m)%preside = 0

                end if
            end do
        else if (func .eq. 1) then
            nis = nspec + is
            do m = ns, ne
                if (pbuf(m)%nid .eq. -1) cycle

                iw = get_interpolation_weight(pbuf(m))

                if (pbuf(m)%preside >= 0) then
                    call distribute_charge(rhodg(is, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), q(is), iw)

                else if (pbuf(m)%preside == OH_PCL_TO_BE_ACCUMULATED) then
                    call distribute_charge(rhobk(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), q(is), iw)
                    call distribute_charge(rhobksp(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps, is), q(is), iw)
                    call distribute_charge(rhodg(nis, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), q(is), iw)
                    nphgram(pbuf(m)%nid + 1, is, ps) = nphgram(pbuf(m)%nid + 1, is, ps) - 1
                    pbuf(m)%nid = -1
                    pbuf(m)%preside = 0

                else if (pbuf(m)%preside == OH_PCL_TO_BE_ACCUMULATED_AS_ANTIPCL) then
                    call distribute_charge(rhobk(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), -q(is), iw)
                    call distribute_charge(rhobksp(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps, is), -q(is), iw)
                    call distribute_charge(rhodg(nis, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), -q(is), iw)
                    nphgram(pbuf(m)%nid + 1, is, ps) = nphgram(pbuf(m)%nid + 1, is, ps) - 1
                    pbuf(m)%nid = -1
                    pbuf(m)%preside = 0

                end if
            end do
        else if (func .eq. 2) then
            nis = nspec + is
            do m = ns, ne
                if (pbuf(m)%nid .eq. -1) cycle

                iw = get_interpolation_weight(pbuf(m))

                if (pbuf(m)%preside >= 0) then
                    call distribute_charge(rho(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), q(is), iw)
                    call distribute_charge(rhodg(is, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), q(is), iw)

                else if (pbuf(m)%preside == OH_PCL_TO_BE_ACCUMULATED) then
                    call distribute_charge(rhobk(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), q(is), iw)
                    call distribute_charge(rhobksp(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps, is), q(is), iw)
                    call distribute_charge(rhodg(nis, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), q(is), iw)
                    nphgram(pbuf(m)%nid + 1, is, ps) = nphgram(pbuf(m)%nid + 1, is, ps) - 1
                    pbuf(m)%nid = -1
                    pbuf(m)%preside = 0

                else if (pbuf(m)%preside == OH_PCL_TO_BE_ACCUMULATED_AS_ANTIPCL) then
                    call distribute_charge(rhobk(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), -q(is), iw)
                    call distribute_charge(rhobksp(1, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps, is), -q(is), iw)
                    call distribute_charge(rhodg(is, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), -q(is), iw)
                    call distribute_charge(rhodg(nis, iw%i:iw%i1, iw%j:iw%j1, iw%k:iw%k1, ps), -q(is), iw)
                    nphgram(pbuf(m)%nid + 1, is, ps) = nphgram(pbuf(m)%nid + 1, is, ps) - 1
                    pbuf(m)%nid = -1
                    pbuf(m)%preside = 0

                end if
            end do
        end if
!
    end do ISL
!============================== species loop end

!-------------------- masking of charge
!      if(iimsk.eq.1) call fsmask(5)
contains

    !> Get interpolation weight with particle position.
    pure function get_interpolation_weight(pcl) result(iweight)
        type(oh_particle), intent(in) :: pcl
        type(t_InterpolationInfo) :: iweight

        real(kind=8) :: xlocal, ylocal, zlocal

        xlocal = pcl%x - xl
        ylocal = pcl%y - yl
        zlocal = pcl%z - zl

        iweight%i = int(xlocal)
        iweight%j = int(ylocal)
        iweight%k = int(zlocal)

        iweight%i1 = iweight%i + 1
        iweight%j1 = iweight%j + 1
        iweight%k1 = iweight%k + 1

        block
            double precision :: x1, y1, z1, z2, xy1, xz1, yz1, xz2, yz2

            x1 = xlocal - iweight%i
            y1 = ylocal - iweight%j
            z1 = zlocal - iweight%k

            xy1 = x1*y1
            xz1 = x1*z1
            yz1 = y1*z1
            z2 = 1.0d0 - z1
            xz2 = x1*z2
            yz2 = y1*z2

            iweight%v3 = xy1*z1
            iweight%v2 = xz1 - iweight%v3
            iweight%v4 = yz1 - iweight%v3
            iweight%v1 = z1 - xz1 - iweight%v4

            iweight%v7 = xy1*z2
            iweight%v6 = xz2 - iweight%v7
            iweight%v8 = yz2 - iweight%v7
            iweight%v5 = z2 - xz2 - iweight%v8
        end block
    end function

    !> Distribute charge with interpolation weight.
    pure subroutine distribute_charge(rh, q_charge, iweight)
        double precision, intent(inout) :: rh(1:2, 1:2, 1:2)
        double precision, intent(in) :: q_charge
        type(t_InterpolationInfo), intent(in) :: iweight

        rh(1, 1, 2) = rh(1, 1, 2) + iweight%v1*q_charge
        rh(2, 1, 2) = rh(2, 1, 2) + iweight%v2*q_charge
        rh(2, 2, 2) = rh(2, 2, 2) + iweight%v3*q_charge
        rh(1, 2, 2) = rh(1, 2, 2) + iweight%v4*q_charge
        rh(1, 1, 1) = rh(1, 1, 1) + iweight%v5*q_charge
        rh(2, 1, 1) = rh(2, 1, 1) + iweight%v6*q_charge
        rh(2, 2, 1) = rh(2, 2, 1) + iweight%v7*q_charge
        rh(1, 2, 1) = rh(1, 2, 1) + iweight%v8*q_charge
    end subroutine

end subroutine

subroutine add_boundary_charge(rh, nelem, uelem, sdom, bound, ctype, fs, myid)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!              A D D _ B O U N D A R Y _ C H A R G E
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a charge distribution to      .
!   .      grid points from particle locations by the          .
!   .      "area-sharing" or "area-weighting" scheme.          .
!   ............................................................
!
!
!-------------------- parameter and common block
    implicit none
    integer(kind=4), intent(in) :: nelem
    integer(kind=4), intent(in) :: fs(2, 3)
    double precision, intent(inout) :: rh(1:nelem, fs(1, 1):fs(2, 1), fs(1, 2):fs(2, 2), fs(1, 3):fs(2, 3))
    integer(kind=4), intent(in) :: uelem
    integer(kind=4), intent(in) :: sdom(2, OH_DIMENSION)
    integer(kind=4), intent(in) :: bound(2, 3)
    integer(kind=4), intent(in) :: ctype(3, 2)
    integer(kind=4), intent(in) :: myid

    integer(kind=4) :: xl, yl, zl, xu, yu, zu
    integer(kind=4) :: lx, ly, lz
    integer(kind=4) :: sl, su, dl, du, nl, nu

    integer(kind=4) :: i, j, k
    integer(kind=4) :: elem

    xl = sdom(1, 1); xu = sdom(2, 1)
    yl = sdom(1, 2); yu = sdom(2, 2)
    zl = sdom(1, 3); zu = sdom(2, 3)

    lx = xu - xl
    ly = yu - yl
    lz = zu - zl

    sl = ctype(2, 2)        !=-1
    su = ctype(2, 1)        !=+1

    nl = ctype(3, 2)        != 1
    nu = ctype(3, 1)        != 1

    dl = sl + nl        != 0
    du = su - nu        != 0

    if (bound(1, 3) .eq. 1) then
        do k = 0, nl - 1                !(i.e., do k=0,0)
        do j = 0, ly + (su + nu - sl) - 1        !(i.e., do j=0,ly+2)
        do i = 0, lx + (su + nu - sl) - 1        !(i.e., do i=0,lx+2)
        do elem = 1, uelem
            rh(elem, sl + i, sl + j, dl + k) = &
                rh(elem, sl + i, sl + j, dl + k) &
                + rh(elem, sl + i, sl + j, sl + k)
        end do
        end do
        end do
        end do
    else
        do k = 0, 0
        do j = 0, ly + (su + nu - sl) - 1        !(i.e., do j=0,ly+2)
        do i = 0, lx + (su + nu - sl) - 1        !(i.e., do i=0,lx+2)
        do elem = 1, uelem
            rh(elem, sl + i, sl + j, dl + k) = rh(elem, sl + i, sl + j, dl + k)*2.0d0
        end do
        end do
        end do
        end do
    end if

!--------------------
    if (bound(2, 3) .eq. 1) then
        do k = 0, nu - 1                !(i.e., do k=0,0)
        do j = 0, ly + (su + nu - sl) - 1        !(i.e., do j=0,ly+2)
        do i = 0, lx + (su + nu - sl) - 1        !(i.e., do i=0,lx+2)
        do elem = 1, uelem
            rh(elem, sl + i, sl + j, lz + du + k) = &
                rh(elem, sl + i, sl + j, lz + du + k) &
                + rh(elem, sl + i, sl + j, lz + su + k)
        end do
        end do
        end do
        end do
    else
        do k = 0, 0
        do j = 0, ly + (su + nu - sl) - 1        !(i.e., do j=0,ly+2)
        do i = 0, lx + (su + nu - sl) - 1        !(i.e., do i=0,lx+2)
        do elem = 1, uelem
            rh(elem, sl + i, sl + j, lz + du + k) = &
                rh(elem, sl + i, sl + j, lz + du + k)*2.0d0
        end do
        end do
        end do
        end do
    end if

    if (bound(1, 2) .eq. 1) then
        do k = 0, lz + (du + nu - dl) - 1        !(i.e., do k=0,lz)
        do j = 0, nl - 1                !(i.e., do j=0,0)
        do i = 0, lx + (su + nu - sl) - 1        !(i.e., do i=0,lx+2)
        do elem = 1, uelem
            rh(elem, sl + i, dl + j, dl + k) = &
                rh(elem, sl + i, dl + j, dl + k) &
                + rh(elem, sl + i, sl + j, dl + k)
        end do
        end do
        end do
        end do
    else
        do k = 0, lz + (du + nu - dl) - 1        !(i.e., do k=0,lz)
        do j = 0, 0
        do i = 0, lx + (su + nu - sl) - 1        !(i.e., do i=0,lx+2)
        do elem = 1, uelem
            rh(elem, sl + i, dl + j, dl + k) = &
                rh(elem, sl + i, dl + j, dl + k)*2.0d0
        end do
        end do
        end do
        end do
    end if

    if (bound(2, 2) .eq. 1) then
        do k = 0, lz + (du + nu - dl) - 1        !(i.e., do k=0,lz)
        do j = 0, nu - 1                !(i.e., do j=0,0)
        do i = 0, lx + (su + nu - sl) - 1        !(i.e., do i=0,lx+2)
        do elem = 1, uelem
            rh(elem, sl + i, ly + du + j, dl + k) = &
                rh(elem, sl + i, ly + du + j, dl + k) &
                + rh(elem, sl + i, ly + su + j, dl + k)
        end do
        end do
        end do
        end do
    else
        do k = 0, lz + (du + nu - dl) - 1        !(i.e., do k=0,lz)
        do j = 0, 0
        do i = 0, lx + (su + nu - sl) - 1        !(i.e., do i=0,lx+2)
        do elem = 1, uelem
            rh(elem, sl + i, ly + du + j, dl + k) = &
                rh(elem, sl + i, ly + du + j, dl + k)*2.0d0
        end do
        end do
        end do
        end do
    end if

    if (bound(1, 1) .eq. 1) then
        do k = 0, lz + (du + nu - dl) - 1        !(i.e., do k=0,lz)
        do j = 0, ly + (du + nu - dl) - 1        !(i.e., do j=0,ly)
        do i = 0, nl - 1                !(i.e., do i=0,0)
        do elem = 1, uelem
            rh(elem, dl + i, dl + j, dl + k) = rh(elem, dl + i, dl + j, dl + k) &
         &                          + rh(elem, sl + i, dl + j, dl + k)
        end do
        end do
        end do
        end do
    else
        do k = 0, lz + (du + nu - dl) - 1        !(i.e., do k=0,lz)
        do j = 0, ly + (du + nu - dl) - 1        !(i.e., do j=0,ly)
        do i = 0, 0
        do elem = 1, uelem
            rh(elem, dl + i, dl + j, dl + k) = rh(elem, dl + i, dl + j, dl + k)*2.0d0
        end do
        end do
        end do
        end do
    end if

    if (bound(2, 1) .eq. 1) then
        do k = 0, lz + (du + nu - dl) - 1        !(i.e., do k=0,lz)
        do j = 0, ly + (du + nu - dl) - 1        !(i.e., do j=0,ly)
        do i = 0, nu - 1                        !(i.e., do i=0,0)
        do elem = 1, uelem
            rh(elem, lx + du + i, dl + j, dl + k) = rh(elem, lx + du + i, dl + j, dl + k) &
         &                             + rh(elem, lx + su + i, dl + j, dl + k)
        end do
        end do
        end do
        end do
    else
        do k = 0, lz + (du + nu - dl) - 1        !(i.e., do k=0,lz)
        do j = 0, ly + (du + nu - dl) - 1        !(i.e., do j=0,ly)
        do i = 0, 0
        do elem = 1, uelem
            rh(elem, lx + du + i, dl + j, dl + k) = rh(elem, lx + du + i, dl + j, dl + k)*2.0d0
        end do
        end do
        end do
        end do
    end if

    return
end subroutine

subroutine add_background_charge(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!            A D D _ B A C K G R O U N D _ C H A R G E
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a charge distribution to      .
!   .      grid points from particle locations by the          .
!   .      "area-sharing" or "area-weighting" scheme.          .
!   ............................................................
!
!
!-------------------- parameter and common block
    use oh_type
    use paramt
    use allcom
    use m_ohhinfo
    implicit none

    integer(kind=4), intent(in) :: ps

    call ohhinfo_update(sdoms(:, :, sdid(ps) + 1))

    if (nflag_testp .ne. 1) then
        rho(1, 0:ngx, 0:ngy, 0:ngz, ps) = &
            rho(1, 0:ngx, 0:ngy, 0:ngz, ps) &
            + rhobk(1, 0:ngx, 0:ngy, 0:ngz, 3)
    else
        rho(1, 0:ngx, 0:ngy, 0:ngz, ps) = &
            rhobk(1, 0:ngx, 0:ngy, 0:ngz, 3)
    end if
end subroutine
