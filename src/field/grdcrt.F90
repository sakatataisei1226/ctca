#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"

subroutine interpolate_zssurf3_electric_field(ps)
    use oh_type
    use paramt
    use allcom
    implicit none

    integer, intent(in) :: ps
    integer :: xl, yl, zl, xu, yu, zu
    integer :: i, j, k

    integer :: ixl, ixu, iyl, iyu, izl, izu

    integer :: isurf
    integer :: zsurfs(2)  ! z-coodinates of XY-surfaces

    logical :: exists_xl_wall, exists_xu_wall
    logical :: exists_yl_wall, exists_yu_wall
    integer :: xl_min, xl_max, yl_min, yl_max, zl_min, zl_max  ! local range of the hole ([xl_min, xl_max] or [yl_min, yl_max])

    ! Calculate local space
    xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
    yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
    zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)

    ! Case: zssurf3 not in local space
    if ((zu < zlrechole(2)) .or. (zssurf < zl)) then
        return
    end if

    ixl = int(xlrechole(1)); ixu = int(xurechole(1))
    iyl = int(ylrechole(1)); iyu = int(yurechole(1))
    izl = int(zlrechole(2)); izu = int(zurechole(1))

    zsurfs = (/izu, izl/)

    ! Interpolate electric field at 1 grid above XY-surface
    do isurf = 1, size(zsurfs)
        k = zsurfs(isurf) - zl
        if ((k < 0) .or. (zu - zl < k)) cycle
        do j = 0, yu - yl
            do i = 0, xu - xl
                if ((ixl < i + xl) .and. (i + xl < ixu) .and. &
                    (iyl < j + yl) .and. (j + yl < iyu)) then
                    cycle
                end if

                wrk(EZ, i, j, k) = 2*eb(EZ, i, j, k, ps) + e0z
            end do
        end do
    end do

    do isurf = 1, size(zsurfs)
        k = zsurfs(isurf) - zl + 1
        if ((k < 0) .or. (zu - zl < k)) cycle
        do j = 0, yu - yl
            do i = 0, xu - xl
                if ((ixl < i + xl) .and. (i + xl < ixu) .and. &
                    (iyl < j + yl) .and. (j + yl < iyu)) then
                    cycle
                end if

                wrk(EZ, i, j, k) = 0.5*eb(EZ, i, j, k, ps) + e0z
            end do
        end do
    end do

    ! Check if a wall of the hole exists
    exists_xl_wall = (xl <= ixl) .and. (ixl <= xu)
    exists_xu_wall = (xl <= ixu) .and. (ixu <= xu)
    exists_yl_wall = (yl <= iyl) .and. (iyl <= yu)
    exists_yu_wall = (yl <= iyu) .and. (iyu <= yu)

    ! Calculate local range of the hole
    xl_min = max(ixl - xl, 0)
    xl_max = min(ixu - xl, xu - xl)
    yl_min = max(iyl - yl, 0)
    yl_max = min(iyu - yl, yu - yl)
    zl_min = max(izl - zl, 0)
    zl_max = min(izu - zl, zu - zl)

    ! Interpolate electric field at 1 grid above XZ-surface and YZ-surface
    do k = zl_min, zl_max
        ! Case: left XZ-surface
        j = iyl - yl
        if (exists_yl_wall) then
            do i = xl_min, xl_max
                wrk(EY, i, j, k) = 2*eb(EY, i, j, k, ps) + e0y
                wrk(EY, i, j + 1, k) = 0.5*eb(EY, i, j + 1, k, ps) + e0y
            end do
        end if

        ! Case: right XZ-surface
        j = iyu - yl
        if (exists_yu_wall) then
            do i = xl_min, xl_max
                wrk(EY, i, j, k) = 2*eb(EY, i, j - 1, k, ps) + e0y
                wrk(EY, i, j - 1, k) = 0.5*eb(EY, i, j - 2, k, ps) + e0y
            end do
        end if

        ! Case: left YZ-surface
        i = int(xlrechole(1)) - xl
        if (exists_xl_wall) then
            do j = yl_min, yl_max
                wrk(EX, i, j, k) = 2*eb(EX, i, j, k, ps) + e0x
                wrk(EX, i + 1, j, k) = 0.5*eb(EX, i + 1, j, k, ps) + e0x
            end do
        end if

        ! Case: right YZ-surface
        i = int(xurechole(1)) - xl
        if (exists_xu_wall) then
            do j = yl_min, yl_max
                wrk(EX, i, j, k) = 2*eb(EX, i - 1, j, k, ps) + e0x
                wrk(EX, i - 1, j, k) = 0.5*eb(EX, i - 2, j, k, ps) + e0x
            end do
        end if
    end do
end subroutine

!
subroutine grdcrt(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   G R D C R T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine relocates electric and magnetic fields  .
!   .  defined at grids for interpolation to particle          .
!   .  positions.  (self-force of particle is eliminated by    .
!   .  the relocation)                                         .
!   ............................................................
!
!-------------------- parameter and common blocks
    use oh_type
    use paramt
    use allcom
    use m_logging
    implicit none

    integer(kind=4) :: i, j, k, ipc, icap
    integer(kind=4) :: xl, yl, zl, xu, yu, zu
    integer(kind=4) :: ps
    real(kind=8) :: eserg

    xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
    yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
    zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)

!-------------------- relocation of pe-fields
    call logging_debuglog('Relocate pe-fields')
    do k = 0, zu - zl
    do j = 0, yu - yl
    do i = 0, xu - xl
        wrk(EX, i, j, k) = 0.5d0*(eb(EX, i, j, k, ps) + eb(EX, i - 1, j, k, ps)) &
                            + e0x
        wrk(EY, i, j, k) = 0.5d0*(eb(EY, i, j, k, ps) + eb(EY, i, j - 1, k, ps)) &
                            + e0y
        wrk(EZ, i, j, k) = 0.5d0*(eb(EZ, i, j, k, ps) + eb(EZ, i, j, k - 1, ps)) &
                            + e0z
    end do
    end do
    end do

    if (zurechole(1) >= 0 .and. boundary_type /= 'flat-surface' .and. boundary_type /= 'hyperboloid-hole') then
        call logging_debuglog('Start interpolate zssurf3')
        call interpolate_zssurf3_electric_field(ps)
        call logging_debuglog('End interpolate zssurf3')
    end if

!-------------------- body surface electric field
    if (sfecrrct .ge. 1) then
        eserg = path(1)*path(1)
        do ipc = 1, npc
            do icap = nscpmx(ipc), nmxcpmx(ipc) - 1
                i = bdygrid(1, icap) - xl
                j = bdygrid(2, icap) - yl
                k = bdygrid(3, icap) - zl
                if (i .ge. 0 .and. i .le. xu - xl .and. &
               &   j .ge. 0 .and. j .le. yu - yl .and. &
               &   k .ge. 0 .and. k .le. zu - zl) then
                    if (abs(eb(EX, i, j, k, ps)) .ge. 1.0d-8*eserg .and. &
                   &   abs(eb(EX, i - 1, j, k, ps)) .lt. 1.0d-8*eserg) then
                        wrk(EX, i, j, k) = eb(EX, i, j, k, ps) + e0x
                    else if (abs(eb(EX, i - 1, j, k, ps)) .ge. 1.0d-8*eserg .and. &
                   &        abs(eb(EX, i, j, k, ps)) .lt. 1.0d-8*eserg) then
                        wrk(EX, i, j, k) = eb(EX, i - 1, j, k, ps) + e0x
                    end if
                    if (abs(eb(EY, i, j, k, ps)) .ge. 1.0d-8*eserg .and. &
                   &   abs(eb(EY, i, j - 1, k, ps)) .lt. 1.0d-8*eserg) then
                        wrk(EY, i, j, k) = eb(EY, i, j, k, ps) + e0y
                    else if (abs(eb(EY, i, j - 1, k, ps)) .ge. 1.0d-8*eserg .and. &
                   &        abs(eb(EY, i, j, k, ps)) .lt. 1.0d-8*eserg) then
                        wrk(EY, i, j, k) = eb(EY, i, j - 1, k, ps) + e0y
                    end if
                    if (abs(eb(EZ, i, j, k, ps)) .ge. 1.0d-8*eserg .and. &
                   &   abs(eb(EZ, i, j, k - 1, ps)) .lt. 1.0d-8*eserg) then
                        wrk(EZ, i, j, k) = eb(EZ, i, j, k, ps) + e0z
                    else if (abs(eb(EZ, i, j, k - 1, ps)) .ge. 1.0d-8*eserg .and. &
                   &        abs(eb(EZ, i, j, k, ps)) .lt. 1.0d-8*eserg) then
                        wrk(EZ, i, j, k) = eb(EZ, i, j, k - 1, ps) + e0z
                    end if
                end if
            end do
        end do
    end if

!-------------------- relocation of pb-fields
    call logging_debuglog('Relocation of pb-fields')
    do k = 0, zu - zl
    do j = 0, yu - yl
    do i = 0, xu - xl
        wrk(BX, i, j, k) = 0.25d0*(eb(BX, i, j, k, ps) + eb(BX, i, j - 1, k, ps) &
       &                      + eb(BX, i, j, k - 1, ps) + eb(BX, i, j - 1, k - 1, ps)) &
       &              + b0x
        wrk(BY, i, j, k) = 0.25d0*(eb(BY, i, j, k, ps) + eb(BY, i - 1, j, k, ps) &
       &                      + eb(BY, i, j, k - 1, ps) + eb(BY, i - 1, j, k - 1, ps)) &
       &              + b0y
        wrk(BZ, i, j, k) = 0.25d0*(eb(BZ, i, j, k, ps) + eb(BZ, i - 1, j, k, ps) &
       &                      + eb(BZ, i, j - 1, k, ps) + eb(BZ, i - 1, j - 1, k, ps)) &
       &              + b0z
    end do
    end do
    end do

!--------------------- masking of longitudinal e-field
!      call fsmask(6)
!-------------------- boundary treatment
!      call fbound(4)
!      call fbound(3)

!-------------------- filtering pex, pey and pez in y and x respectively
!      call filtr3(pex,ix,iy,iz,nxm,nym,nzm,ipexfl)
!      call filtr3(pey,ix,iy,iz,nxm,nym,nzm,ipeyfl)
!      call filtr3(pez,ix,iy,iz,nxm,nym,nzm,ipezfl)
!      call fbound(3)

    return
end subroutine
