#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
subroutine poisson(func)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   P O I S S O N
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   subroutine for solving linear equations for poisson eq..
!   ............................................................
!
!-------------------- parameter and common blocks
    use oh_type
    use paramt
    use allcom
    use m_mypoisson
#define MCW local_comm
#define MIP MPI_IN_PLACE
    implicit none
!
    integer(kind=4), intent(in) :: func
    integer(kind=4) :: i, j, k, ii, jj, kk, iii, ifloor, islice
    integer(kind=4) :: xl, xu, yl, yu, zl, zu
    integer(kind=4) :: iu, ju, ku
    integer(kind=4) :: siu, sju, sku, stiu
    integer(kind=4) :: ps
    integer(kind=4) :: sd, ifm, from, to
    real(kind=8) :: capc1
    real(kind=8) :: divisor

    integer(kind=4) :: lfloor, ufloor
    integer(kind=4) :: nxsbar, nysbar, nzsbar
    integer(kind=4) :: count, icomreq
    integer(kind=4) :: src, dst
    integer(kind=4) :: stag = 0, rtag = 0
    integer(kind=4) :: mpierr

    ! call solve_poisson
    ! return

!--------------------
    xl = sdoms(1, 1, sdid(1) + 1); xu = sdoms(2, 1, sdid(1) + 1)
    yl = sdoms(1, 2, sdid(1) + 1); yu = sdoms(2, 2, sdid(1) + 1)
    zl = sdoms(1, 3, sdid(1) + 1); zu = sdoms(2, 3, sdid(1) + 1)
    if (idxsd .ne. nodes(1) - 1) then
        iu = xu - xl - 1
    else
        iu = xu - xl
    end if
    if (idysd .ne. nodes(2) - 1) then
        ju = yu - yl - 1
    else
        ju = yu - yl
    end if
    if (idzsd .ne. nodes(3) - 1) then
        ku = zu - zl - 1
    else
        ku = zu - zl
    end if
    if (pftmode .eq. 3 .and. myid .ne. snode - 1) then
        siu = sxu - sxl - 1
    else
        siu = sxu - sxl
    end if
    if ((pftmode .eq. 1 .or. pftmode .eq. 2) .and. myid .ne. snode - 1) then
        stiu = stxu - stxl - 1
    else
        stiu = stxu - stxl
    end if
    if (pftmode .eq. 2 .and. myid .ne. snode - 1) then
        sju = syu - syl - 1
    else
        sju = syu - syl
    end if
    if ((pftmode .eq. 1 .or. pftmode .eq. 4) .and. myid .ne. snode - 1) then
        sku = szu - szl - 1
    else
        sku = szu - szl
    end if
    rcnts(1:nnode - 1) = slcs(1:nnode - 1, 3)
    rcnts(nnode) = slcs(nnode, 3) + 1

!--------------------
    ! Create poisson working array.

    ! Caluclate rho + boundary conditions (Neumman's condition)
    rho_tilde(1, :, :, :, 1) = rho(1, :, :, :, 1)

    ! X-axis
    if (mtd_vbnd(1) == 2) then
        if (xl <= 0) then
            rho_tilde(1, 0, :, :, 1) = rho_tilde(1, 0, :, :, 1) &
                                       - boundary_conditions(1, 1)*rene*2.0d0
        end if
        if (xu >= nx) then
            rho_tilde(1, xu - xl, :, :, 1) = rho_tilde(1, xu - xl, :, :, 1) &
                                             + boundary_conditions(2, 1)*rene*2.0d0
        end if
    end if

    ! Y-axis
    if (mtd_vbnd(2) == 2) then
        if (yl <= 0) then
            rho_tilde(1, :, 0, :, 1) = rho_tilde(1, :, 0, :, 1) &
                                       - boundary_conditions(1, 2)*rene*2.0d0
        end if
        if (yu >= ny) then
            rho_tilde(1, :, yu - yl, :, 1) = rho_tilde(1, :, yu - yl, :, 1) &
                                             + boundary_conditions(2, 2)*rene*2.0d0
        end if
    end if

    ! Z-axis
    if (mtd_vbnd(3) == 2) then
        if (zl <= 0) then
            rho_tilde(1, :, :, 0, 1) = rho_tilde(1, :, :, 0, 1) &
                                       - boundary_conditions(1, 3)*rene*2.0d0
        end if
        if (zu >= nz) then
            rho_tilde(1, :, :, zu - zl, 1) = rho_tilde(1, :, :, zu - zl, 1) &
                                             + boundary_conditions(2, 3)*rene*2.0d0
        end if
    end if

    if (mod(func, 2) .eq. 1) then
    if (pftmode .eq. 1) then
        icomreq = 0
        if (myid .lt. snode) then
!          lfloor = floor(dble(szl/nzsd)); ufloor = floor(dble((szu-1)/nzsd))
            lfloor = int(szl/nzsd); ufloor = int((szu - 1)/nzsd)
            do ifloor = lfloor, ufloor
                count = min(szu, (ifloor + 1)*nzsd) - max(szl, ifloor*nzsd)
                if (myid .eq. snode - 1 .and. ifloor .eq. ufloor) &
               &   count = count + 1
                do j = 0, nodes(2) - 1
                do i = 0, nodes(1) - 1
                    src = i + j*nodes(1) + ifloor*nodes(1)*nodes(2)
                    icomreq = icomreq + 1
                    call MPI_Irecv(poi(1, i*nxsd, j*nysd, max(0, ifloor*nzsd - szl), 1), &
                   &               count, mptype_rs(xtype(i), ytype(j), 2, 1), &
                   &               src, 0, MCW, ireqs(icomreq), mpierr)
                end do
                end do
            end do
        end if
!
!        do dst=floor(dble(zl/nzslc)),floor(dble((zu-1)/nzslc))
        do dst = int(zl/nzslc), int((zu - 1)/nzslc)
            count = 0
            do islice = dst*nzslc, (dst + 1)*nzslc - 1
                if (zl .le. islice .and. islice .lt. zu) &
               &  count = count + 1
            end do
            if (dst .eq. snode - 1 .and. (dst + 1)*nzslc .eq. zu) count = count + 1
            icomreq = icomreq + 1
            call MPI_Isend(rho_tilde(1, 0, 0, max(dst*nzslc - zl, 0), 1), &
           &               count, mptype_sc(xtype(idxsd), ytype(idysd), 2, 1), &
           &               dst, 0, MCW, ireqs(icomreq), mpierr)
        end do
!
        call MPI_Waitall(icomreq, ireqs, istatus, mpierr)
!
!     ---------------
    else if (pftmode .eq. 2) then
        icomreq = 0
        if (myid .lt. snode) then
!          lfloor = floor(dble(syl/nysd)); ufloor = floor(dble((syu-1)/nysd))
            lfloor = int(syl/nysd); ufloor = int((syu - 1)/nysd)
            do ifloor = lfloor, ufloor
                count = min(syu, (ifloor + 1)*nysd) - max(syl, ifloor*nysd)
                if (myid .eq. snode - 1 .and. ifloor .eq. ufloor) &
               &   count = count + 1
                do k = 0, nodes(3) - 1
                do i = 0, nodes(1) - 1
                    src = i + ifloor*nodes(1) + k*nodes(1)*nodes(2)
                    icomreq = icomreq + 1
                    call MPI_Irecv(poi(1, i*nxsd, max(0, ifloor*nysd - syl), k*nzsd, 1), &
                   &               count, mptype_rs(xtype(i), ztype(k), 2, 2), &
                   &               src, 0, MCW, ireqs(icomreq), mpierr)
                end do
                end do
            end do
        end if
!
!        do dst=floor(dble(yl/nyslc)),floor(dble((yu-1)/nyslc))
        do dst = int(yl/nyslc), int((yu - 1)/nyslc)
            count = 0
            do islice = dst*nyslc, (dst + 1)*nyslc - 1
                if (yl .le. islice .and. islice .lt. yu) &
               &  count = count + 1
            end do
            if (dst .eq. snode - 1) count = count + 1
            icomreq = icomreq + 1
            call MPI_Isend(rho_tilde(1, 0, max(dst*nyslc - yl, 0), 0, 1), &
           &               count, mptype_sc(xtype(idxsd), ztype(idzsd), 2, 2), &
           &               dst, 0, MCW, ireqs(icomreq), mpierr)
        end do
!
        call MPI_Waitall(icomreq, ireqs, istatus, mpierr)
!
!     ---------------
    else if (pftmode .eq. 3) then
        poi(1, 0:iu, 0:ny, 0:nz, 1) = rho_tilde(1, 0:iu, 0:ny, 0:nz, 1)
!
!     ---------------
    else if (pftmode .eq. 4) then
        poi(1, 0:nx, 0:ny, 0:ku, 1) = rho_tilde(1, 0:nx, 0:ny, 0:ku, 1)
    end if
    end if

!--------------------
    ! Apply FFT to poisson working array.
    if (myid .lt. snode) then
    if (pftmode .eq. 1) then
!       ------------- FFT-xy
        do k = 0 - szl, nz - szl
        if (k .ge. 0 .and. k .le. sku) then
            call dfftw_execute_r2r(fftplan(CP, 1, 1), &
           &  poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1), &
           &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1))
            poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1) = &
           &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1)
        end if
        end do
!       ------------- transpose x-z
        call MPI_Alltoall(poi(1, 0, 0, 0, 1), 1, mptype_aa(0, 0, 2, 1), &
       &                  poi(1, 0, 0, 0, 2), 1, mptype_aa(1, 0, 2, 1), &
       &                  subcomm, mpierr)
        iii = 0
        do ii = 0, (nxslc + 1)*(snode - 1), nxslc + 1
            nxsbar = stiu
            if (ii .eq. (nxslc + 1)*(snode - 1)) then
                nzsbar = nzslc
            else
                nzsbar = nzslc - 1
            end if
            do k = 0, nxsbar
                do i = 0, nzsbar
                    poi(1, iii + i, 0:ny, k, 1) = poi(1, ii + k, 0:ny, i, 2)
                end do
            end do
            iii = iii + nzslc
        end do
!       ------------- FFT-z
        do k = 0 - stxl, nx - stxl
        do j = 0, ny
            if (k .ge. 0 .and. k .le. stiu) then
                call dfftw_execute_r2r(fftplan(CP, 2, 1), &
               &  poi(1, lzfft(CP):uzfft(CP), j, k, 1), &
               &  warray1d(lzfft(CP):uzfft(CP), 1))
                poi(1, lzfft(CP):uzfft(CP), j, k, 1) = &
               &  warray1d(lzfft(CP):uzfft(CP), 1)
            end if
        end do
        end do
!
!     ---------------
    else if (pftmode .eq. 2) then
!       ------------- FFT-xz
        do j = 0 - syl, ny - syl
        if (j .ge. 0 .and. j .le. sju) then
            call dfftw_execute_r2r(fftplan(CP, 1, 1), &
           &  poi(1, lxfft(CP):uxfft(CP), j, lzfft(CP):uzfft(CP), 1), &
           &  warray2d(lxfft(CP):uxfft(CP), lzfft(CP):uzfft(CP), 1))
            poi(1, lxfft(CP):uxfft(CP), j, lzfft(CP):uzfft(CP), 1) = &
           &  warray2d(lxfft(CP):uxfft(CP), lzfft(CP):uzfft(CP), 1)
        end if
        end do
!       ------------- transpose x-y
        call MPI_Alltoall(poi(1, 0, 0, 0, 1), 1, mptype_aa(0, 0, 2, 2), &
       &                  poi(1, 0, 0, 0, 2), 1, mptype_aa(1, 0, 2, 2), &
       &                  subcomm, mpierr)
        iii = 0
        do ii = 0, (nxslc + 1)*(snode - 1), nxslc + 1
            nxsbar = stiu
            if (ii .eq. (nxslc + 1)*(snode - 1)) then
                nysbar = nyslc
            else
                nysbar = nyslc - 1
            end if
            do j = 0, nxsbar
                do i = 0, nysbar
                    poi(1, iii + i, j, 0:nz, 1) = poi(1, ii + j, i, 0:nz, 2)
                end do
            end do
            iii = iii + nyslc
        end do
!       ------------- FFT-y
        do k = 0, nz
        do j = 0 - stxl, nx - stxl
            if (j .ge. 0 .and. j .le. stiu) then
                call dfftw_execute_r2r(fftplan(CP, 2, 1), &
               &  poi(1, lyfft(CP):uyfft(CP), j, k, 1), &
               &  warray1d(lyfft(CP):uyfft(CP), 1))
                poi(1, lyfft(CP):uyfft(CP), j, k, 1) = &
               &  warray1d(lyfft(CP):uyfft(CP), 1)
            end if
        end do
        end do
!
!     ---------------
    else if (pftmode .eq. 3) then
!       ------------- FFT-yz
        do i = 0 - sxl, nx - sxl
        if (i .ge. 0 .and. i .le. siu) then
            call dfftw_execute_r2r(fftplan(CP, 1, 1), &
           &  poi(1, i, lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1), &
           &  warray2d(lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1))
            poi(1, i, lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1) = &
           &  warray2d(lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1)
        end if
        end do
!       ------------- gather in x-dir
        call MPI_Allgatherv(poi(1, 0, -1, -1, 1), siu + 1, mptype_yz(2), &
       &                    poi(1, 0, -1, -1, 2), rcnts, slcs(:, 1), mptype_yz(2), &
       &                    subcomm, mpierr)
!       ------------- FFT-x
        do k = 0, nz
        do j = 0, ny
            call dfftw_execute_r2r(fftplan(CP, 2, 1), &
           &  poi(1, lxfft(CP):uxfft(CP), j, k, 2), &
           &  warray1d(lxfft(CP):uxfft(CP), 1))
            poi(1, lxfft(CP):uxfft(CP), j, k, 1) = &
           &  warray1d(lxfft(CP):uxfft(CP), 1)
        end do
        end do
!
!     ---------------
    else if (pftmode .eq. 4) then
!       ------------- FFT-xy
        do k = 0 - szl, nz - szl
        if (k .ge. 0 .and. k .le. sku) then
            call dfftw_execute_r2r(fftplan(CP, 1, 1), &
           &  poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1), &
           &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1))
            poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1) = &
           &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1)
        end if
        end do
!       ------------- gather in z-dir
        call MPI_Allgatherv(poi(1, -1, -1, 0, 1), sku + 1, mptype_xy(2), &
       &                    poi(1, -1, -1, 0, 2), rcnts, slcs(:, 1), mptype_xy(2), &
       &                    subcomm, mpierr)
!       ------------- FFT-z
        do j = 0, ny
        do i = 0, nx
            call dfftw_execute_r2r(fftplan(CP, 2, 1), &
           &  poi(1, i, j, lzfft(CP):uzfft(CP), 2), &
           &  warray1d(lzfft(CP):uzfft(CP), 1))
            poi(1, i, j, lzfft(CP):uzfft(CP), 1) = &
           &  warray1d(lzfft(CP):uzfft(CP), 1)
        end do
        end do
    end if
    end if

    ! poi

!--------------------
    ! Solve equation in frequency space
    if (myid .lt. snode) then
    if (pftmode .eq. 1) then
        do k = 0 - stxl, nx - stxl
        do j = 0, ny
        do i = 0, nz
            if (k .ge. 0 .and. k .le. stiu) then
                kk = k + stxl
                divisor = (kmod(1, kk, 2) + kmod(2, j, 2) + kmod(3, i, 2))*mfactor(2)
                if (divisor .ne. 0.0d0) then
                    poi(1, i, j, k, 1) = poi(1, i, j, k, 1)/divisor
                else
                    poi(1, i, j, k, 1) = 0.0d0
                end if
            end if
        end do
        end do
        end do
    else if (pftmode .eq. 2) then
        do k = 0, nz
        do j = 0 - stxl, nx - stxl
        do i = 0, ny
            if (j .ge. 0 .and. j .le. stiu) then
                jj = j + stxl
                divisor = (kmod(1, jj, 2) + kmod(2, i, 2) + kmod(3, k, 2))*mfactor(2)
                if (divisor .ne. 0.0d0) then
                    poi(1, i, j, k, 1) = poi(1, i, j, k, 1)/divisor
                else
                    poi(1, i, j, k, 1) = 0.0d0
                end if
            end if
        end do
        end do
        end do
    else if (pftmode .eq. 3 .or. pftmode .eq. 4) then
        do k = 0, nz
        do j = 0, ny
        do i = 0, nx
            divisor = (kmod(1, i, 2) + kmod(2, j, 2) + kmod(3, k, 2))*mfactor(2)
            if (divisor .ne. 0.0d0) then
                poi(1, i, j, k, 1) = poi(1, i, j, k, 1)/divisor
            else
                poi(1, i, j, k, 1) = 0.0d0
            end if
        end do
        end do
        end do
    end if
    end if

!--------------------
    ! Apply IFFT to poisson working array.
    if (myid .lt. snode) then
    if (pftmode .eq. 1) then
!       ------------- IFFT-z
        do k = 0 - stxl, nx - stxl
        do j = 0, ny
            if (k .ge. 0 .and. k .le. stiu) then
                call dfftw_execute_r2r(fftplan(CP, 2, 2), &
               &  poi(1, lzfft(CP):uzfft(CP), j, k, 1), &
               &  warray1d(lzfft(CP):uzfft(CP), 1))
                poi(1, lzfft(CP):uzfft(CP), j, k, 1) = &
               &  warray1d(lzfft(CP):uzfft(CP), 1)
            end if
        end do
        end do
!       ------------- transpose z-x
        call MPI_Alltoall(poi(1, 0, 0, 0, 1), 1, mptype_aa(0, 1, 2, 1), &
       &                  poi(1, 0, 0, 0, 2), 1, mptype_aa(1, 1, 2, 1), &
       &                  subcomm, mpierr)
        iii = 0
        do ii = 0, (nzslc + 1)*(snode - 1), nzslc + 1
            if (ii .eq. (nzslc + 1)*(snode - 1)) then
                nxsbar = nxslc
            else
                nxsbar = nxslc - 1
            end if
            nzsbar = sku
            do k = 0, nzsbar
                do i = 0, nxsbar
                    poi(1, iii + i, 0:ny, k, 1) = poi(1, ii + k, 0:ny, i, 2)
                end do
            end do
            iii = iii + nxslc
        end do
!       ------------- IFFT-xy
        do k = 0 - szl, nz - szl
            if (k .ge. 0 .and. k .le. sku) then
                call dfftw_execute_r2r(fftplan(CP, 1, 2), &
               &  poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1), &
               &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1))
                poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1) = &
               &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1)
            end if
        end do
!
!     ---------------
    else if (pftmode .eq. 2) then
!       ------------- IFFT-y
        do k = 0, nz
        do j = 0 - stxl, nx - stxl
            if (j .ge. 0 .and. j .le. stiu) then
                call dfftw_execute_r2r(fftplan(CP, 2, 2), &
               &  poi(1, lyfft(CP):uyfft(CP), j, k, 1), &
               &  warray1d(lyfft(CP):uyfft(CP), 1))
                poi(1, lyfft(CP):uyfft(CP), j, k, 1) = &
               &  warray1d(lyfft(CP):uyfft(CP), 1)
            end if
        end do
        end do
!       ------------- transpose y-x
        call MPI_Alltoall(poi(1, 0, 0, 0, 1), 1, mptype_aa(0, 1, 2, 2), &
       &                  poi(1, 0, 0, 0, 2), 1, mptype_aa(1, 1, 2, 2), &
       &                  subcomm, mpierr)
        iii = 0
        do ii = 0, (nyslc + 1)*(snode - 1), nyslc + 1
            if (ii .eq. (nyslc + 1)*(snode - 1)) then
                nxsbar = nxslc
            else
                nxsbar = nxslc - 1
            end if
            nysbar = sju
            do j = 0, nysbar
                do i = 0, nxsbar
                    poi(1, iii + i, j, 0:nz, 1) = poi(1, ii + j, i, 0:nz, 2)
                end do
            end do
            iii = iii + nxslc
        end do
!       ------------- IFFT-xz
        do j = 0 - syl, ny - syl
            if (j .ge. 0 .and. j .le. sju) then
                call dfftw_execute_r2r(fftplan(CP, 1, 2), &
               &  poi(1, lxfft(CP):uxfft(CP), j, lzfft(CP):uzfft(CP), 1), &
               &  warray2d(lxfft(CP):uxfft(CP), lzfft(CP):uzfft(CP), 1))
                poi(1, lxfft(CP):uxfft(CP), j, lzfft(CP):uzfft(CP), 1) = &
               &  warray2d(lxfft(CP):uxfft(CP), lzfft(CP):uzfft(CP), 1)
            end if
        end do
!
!     ---------------
    else if (pftmode .eq. 3) then
!       ------------- IFFT-x
        do k = 0, nz
        do j = 0, ny
            call dfftw_execute_r2r(fftplan(CP, 2, 2), &
           &  poi(1, lxfft(CP):uxfft(CP), j, k, 1), &
           &  warray1d(lxfft(CP):uxfft(CP), 1))
            poi(1, lxfft(CP):uxfft(CP), j, k, 1) = &
           &  warray1d(lxfft(CP):uxfft(CP), 1)
        end do
        end do
!       ------------- IFFT-yz
        if (mtd_vbnd(1) .eq. 0 .and. myid .eq. snode - 1) then
        do i = 0, 1
            call dfftw_execute_r2r(fftplan(CP, 1, 2), &
           &  poi(1, i, lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1), &
           &  warray2d(lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1))
            poi(1, i, lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1) = &
           &  warray2d(lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1)
        end do
        end if
        do i = max(sxl - 1, lxfft(CP)), min(sxu + 1, uxfft(CP))
            call dfftw_execute_r2r(fftplan(CP, 1, 2), &
           &  poi(1, i, lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1), &
           &  warray2d(lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1))
            poi(1, i, lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1) = &
           &  warray2d(lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1)
        end do
        if (mtd_vbnd(1) .eq. 0 .and. myid .eq. 0) then
            i = nx - 1
            call dfftw_execute_r2r(fftplan(CP, 1, 2), &
           &  poi(1, i, lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1), &
           &  warray2d(lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1))
            poi(1, i, lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1) = &
           &  warray2d(lyfft(CP):uyfft(CP), lzfft(CP):uzfft(CP), 1)
        end if
!
!     ---------------
    else if (pftmode .eq. 4) then
!       ------------- IFFT-z
        do j = 0, ny
        do i = 0, nx
            call dfftw_execute_r2r(fftplan(CP, 2, 2), &
           &  poi(1, i, j, lzfft(CP):uzfft(CP), 1), &
           &  warray1d(lzfft(CP):uzfft(CP), 1))
            poi(1, i, j, lzfft(CP):uzfft(CP), 1) = &
           &  warray1d(lzfft(CP):uzfft(CP), 1)
        end do
        end do
!       ------------- IFFT-xy
        if (mtd_vbnd(3) .eq. 0 .and. myid .eq. snode - 1) then
        do k = 0, 1
            call dfftw_execute_r2r(fftplan(CP, 1, 2), &
           &  poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1), &
           &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1))
            poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1) = &
           &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1)
        end do
        end if
        do k = max(szl - 1, lzfft(CP)), min(szu + 1, uzfft(CP))
            call dfftw_execute_r2r(fftplan(CP, 1, 2), &
           &  poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1), &
           &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1))
            poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1) = &
           &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1)
        end do
        if (mtd_vbnd(3) .eq. 0 .and. myid .eq. 0) then
            k = nz - 1
            call dfftw_execute_r2r(fftplan(CP, 1, 2), &
           &  poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1), &
           &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1))
            poi(1, lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), k, 1) = &
           &  warray2d(lxfft(CP):uxfft(CP), lyfft(CP):uyfft(CP), 1)
        end if
    end if
    end if

!--------------------
    if (myid .lt. snode) then
    if (pftmode .eq. 1) then
        poi(1, -1:lxfft(CP) - 1, :, :, 1) = 0.0d0
        poi(1, uxfft(CP) + 1:nx + 1, :, :, 1) = 0.0d0
        poi(1, :, -1:lyfft(CP) - 1, :, 1) = 0.0d0
        poi(1, :, uyfft(CP) + 1:ny + 1, :, 1) = 0.0d0
!
        if (myid .eq. 0) &
       &  poi(1, :, :, -1:lzfft(CP) - 1, 1) = 0.0d0
        if (myid .eq. snode - 1) &
       &  poi(1, :, :, uzfft(CP) + 1 - szl:nhslc + 1, 1) = 0.0d0
!
        if (mtd_vbnd(1) .eq. 0) then
            poi(1, -1:-1, :, :, 1) = poi(1, nx - 1:nx - 1, :, :, 1)
            poi(1, nx:nx, :, :, 1) = poi(1, 0:0, :, :, 1)
            poi(1, nx + 1:nx + 1, :, :, 1) = poi(1, +1:+1, :, :, 1)
        end if
!
        if (mtd_vbnd(2) .eq. 0) then
            poi(1, :, -1:-1, :, 1) = poi(1, :, ny - 1:ny - 1, :, 1)
            poi(1, :, ny:ny, :, 1) = poi(1, :, 0:0, :, 1)
            poi(1, :, ny + 1:ny + 1, :, 1) = poi(1, :, +1:+1, :, 1)
        end if
!
        call MPI_sendrecv(poi(1, -1, -1, szu - szl - 1, 1), 1, &
       &                  mptype_xy(2), uslice(2), stag, &
       &                  poi(1, -1, -1, -1, 1), 1, &
       &                  mptype_xy(2), lslice(2), rtag, &
       &                  subcomm, MPI_STATUS_IGNORE, mpierr)
        if (nzslc .eq. 1) then
            call MPI_sendrecv(poi(1, -1, -1, 0, 1), 1, &
           &                  mptype_xy(2), lslice(2), stag, &
           &                  poi(1, -1, -1, szu - szl, 1), 1, &
           &                  mptype_xy(2), uslice(2), rtag, &
           &                  subcomm, MPI_STATUS_IGNORE, mpierr)
            call MPI_sendrecv(poi(1, -1, -1, 1, 1), 1, &
           &                  mptype_xy(2), lslice(2), stag, &
           &                  poi(1, -1, -1, szu - szl + 1, 1), 1, &
           &                  mptype_xy(2), uslice(2), rtag, &
           &                  subcomm, MPI_STATUS_IGNORE, mpierr)
        else
            call MPI_sendrecv(poi(1, -1, -1, 0, 1), 2, &
           &                  mptype_xy(2), lslice(2), stag, &
           &                  poi(1, -1, -1, szu - szl, 1), 2, &
           &                  mptype_xy(2), uslice(2), rtag, &
           &                  subcomm, MPI_STATUS_IGNORE, mpierr)
        end if
!
!     ---------------
    else if (pftmode .eq. 2) then
        poi(1, -1:lxfft(CP) - 1, :, :, 1) = 0.0d0
        poi(1, uxfft(CP) + 1:nx + 1, :, :, 1) = 0.0d0
        poi(1, :, :, -1:lzfft(CP) - 1, 1) = 0.0d0
        poi(1, :, :, uzfft(CP) + 1:nz + 1, 1) = 0.0d0
!
        if (myid .eq. 0) &
       &  poi(1, :, -1:lyfft(CP) - 1, :, 1) = 0.0d0
        if (myid .eq. snode - 1) &
       &  poi(1, :, uyfft(CP) + 1 - syl:ndslc + 1, :, 1) = 0.0d0
!
        if (mtd_vbnd(1) .eq. 0) then
            poi(1, -1:-1, :, :, 1) = poi(1, nx - 1:nx - 1, :, :, 1)
            poi(1, nx:nx, :, :, 1) = poi(1, 0:0, :, :, 1)
            poi(1, nx + 1:nx + 1, :, :, 1) = poi(1, +1:+1, :, :, 1)
        end if
!
        if (mtd_vbnd(3) .eq. 0) then
            poi(1, :, :, -1:-1, 1) = poi(1, :, :, nz - 1:nz - 1, 1)
            poi(1, :, :, nz:nz, 1) = poi(1, :, :, 0:0, 1)
            poi(1, :, :, nz + 1:nz + 1, 1) = poi(1, :, :, +1:+1, 1)
        end if
!
        call MPI_sendrecv(poi(1, -1, syu - syl - 1, -1, 1), 1, &
       &                  mptype_xz(2), uslice(2), stag, &
       &                  poi(1, -1, -1, -1, 1), 1, &
       &                  mptype_xz(2), lslice(2), rtag, &
       &                  subcomm, MPI_STATUS_IGNORE, mpierr)
        if (nyslc .eq. 1) then
            call MPI_sendrecv(poi(1, -1, 0, -1, 1), 1, &
           &                  mptype_xz(2), lslice(2), stag, &
           &                  poi(1, -1, syu - syl, -1, 1), 1, &
           &                  mptype_xz(2), uslice(2), rtag, &
           &                  subcomm, MPI_STATUS_IGNORE, mpierr)
            call MPI_sendrecv(poi(1, -1, 1, -1, 1), 1, &
           &                  mptype_xz(2), lslice(2), stag, &
           &                  poi(1, -1, syu - syl + 1, -1, 1), 1, &
           &                  mptype_xz(2), uslice(2), rtag, &
           &                  subcomm, MPI_STATUS_IGNORE, mpierr)
        else
            call MPI_sendrecv(poi(1, -1, 0, -1, 1), 2, &
           &                  mptype_xz(2), lslice(2), stag, &
           &                  poi(1, -1, syu - syl, -1, 1), 2, &
           &                  mptype_xz(2), uslice(2), rtag, &
           &                  subcomm, MPI_STATUS_IGNORE, mpierr)
        end if
!
!     ---------------
    else if (pftmode .eq. 3) then
        poi(1, :, -1:lyfft(CP) - 1, :, 1) = 0.0d0
        poi(1, :, uyfft(CP) + 1:ny + 1, :, 1) = 0.0d0
        poi(1, :, :, -1:lzfft(CP) - 1, 1) = 0.0d0
        poi(1, :, :, uzfft(CP) + 1:nz + 1, 1) = 0.0d0
!
        if (myid .eq. 0) &
       &  poi(1, -1:lxfft(CP) - 1, :, :, 1) = 0.0d0
        if (myid .eq. snode - 1) &
       &  poi(1, uxfft(CP) + 1:nwslc + 1, :, :, 1) = 0.0d0
!
        if (mtd_vbnd(1) .eq. 0) then
            if (myid .eq. 0) &
           &  poi(1, -1:-1, :, :, 1) = poi(1, nx - 1:nx - 1, :, :, 1)
            if (myid .eq. snode - 1) then
                poi(1, nx:nx, :, :, 1) = poi(1, 0:0, :, :, 1)
                poi(1, nx + 1:nx + 1, :, :, 1) = poi(1, +1:+1, :, :, 1)
            end if
        end if
!
        if (mtd_vbnd(2) .eq. 0) then
            poi(1, :, -1:-1, :, 1) = poi(1, :, ny - 1:ny - 1, :, 1)
            poi(1, :, ny:ny, :, 1) = poi(1, :, 0:0, :, 1)
            poi(1, :, ny + 1:ny + 1, :, 1) = poi(1, :, +1:+1, :, 1)
        end if
!
        if (mtd_vbnd(3) .eq. 0) then
            poi(1, :, :, -1:-1, 1) = poi(1, :, :, nz - 1:nz - 1, 1)
            poi(1, :, :, nz:nz, 1) = poi(1, :, :, 0:0, 1)
            poi(1, :, :, nz + 1:nz + 1, 1) = poi(1, :, :, +1:+1, 1)
        end if
!
!     ---------------
    else if (pftmode .eq. 4) then
        poi(1, -1:lxfft(CP) - 1, :, :, 1) = 0.0d0
        poi(1, uxfft(CP) + 1:nx + 1, :, :, 1) = 0.0d0
        poi(1, :, -1:lyfft(CP) - 1, :, 1) = 0.0d0
        poi(1, :, uyfft(CP) + 1:ny + 1, :, 1) = 0.0d0
!
        if (myid .eq. 0) &
       &  poi(1, :, :, -1:lzfft(CP) - 1, 1) = 0.0d0
        if (myid .eq. snode - 1) &
       &  poi(1, :, :, uzfft(CP) + 1:nhslc + 1, 1) = 0.0d0
!
        if (mtd_vbnd(1) .eq. 0) then
            poi(1, -1:-1, :, :, 1) = poi(1, nx - 1:nx - 1, :, :, 1)
            poi(1, nx:nx, :, :, 1) = poi(1, 0:0, :, :, 1)
            poi(1, nx + 1:nx + 1, :, :, 1) = poi(1, +1:+1, :, :, 1)
        end if
!
        if (mtd_vbnd(2) .eq. 0) then
            poi(1, :, -1:-1, :, 1) = poi(1, :, ny - 1:ny - 1, :, 1)
            poi(1, :, ny:ny, :, 1) = poi(1, :, 0:0, :, 1)
            poi(1, :, ny + 1:ny + 1, :, 1) = poi(1, :, +1:+1, :, 1)
        end if
!
        if (mtd_vbnd(3) .eq. 0) then
            if (myid .eq. 0) &
           &  poi(1, :, :, -1:-1, 1) = poi(1, :, :, nz - 1:nz - 1, 1)
            if (myid .eq. snode - 1) then
                poi(1, :, :, nz:nz, 1) = poi(1, :, :, 0:0, 1)
                poi(1, :, :, nz + 1:nz + 1, 1) = poi(1, :, :, +1:+1, 1)
            end if
        end if
    end if
    end if

!--------------------
    ! Send from poisson working array to phi.
    if (int(func/2) .eq. 1) then
    if (pftmode .eq. 1) then
        icomreq = 0
        do ps = 1, 1
            if (sdid(ps) .eq. -1) cycle
            xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
            yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
            zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)
            do src = int(zl/nzslc), int((zu - 1)/nzslc)
                count = 0
                do islice = src*nzslc, (src + 1)*nzslc - 1
                    if (zl .le. islice .and. islice .lt. zu) then
                        count = count + 1
                    end if
                end do
!            if(count.ne.nzslc) print*, "nzslc1", count, nzslc
                ztype(1) = 0
                if (src .eq. int(zl/nzslc)) then
                    count = count + 1
                    ztype(1) = ztype(1) - 1
                end if
                if (src .eq. int((zu - 1)/nzslc)) then
                    count = count + 2
                end if
                icomreq = icomreq + 1
                call MPI_Irecv(phi(1, 0, 0, max(src*nzslc - zl, 0), ps), &
               &               count, mptype_rc(ztype(1), 2, 1), &
               &               src, 0, MCW, ireqs(icomreq), mpierr)
            end do
        end do
!
        if (myid .lt. snode) then
            lfloor = int(szl/nzsd); ufloor = int((szu - 1)/nzsd)
            do ifloor = lfloor, ufloor
                count = min(szu, (ifloor + 1)*nzsd) - max(szl, ifloor*nzsd)
                ! if (count .ne. nzslc) print *, "nzslc2", count, nzslc
                ztype(1) = 0
                if (mod(max(szl, ifloor*nzsd), nzsd) .eq. 0) then
                    count = count + 1
                    ztype(1) = ztype(1) - 1
                end if
                if (mod(min(szu, (ifloor + 1)*nzsd), nzsd) .eq. 0) then
                    count = count + 2
                end if
                do j = 0, nodes(2) - 1
                do i = 0, nodes(1) - 1
                    sd = i + j*nodes(1) + ifloor*nodes(1)*nodes(2)
!              from = famind(sd+1) + 1
!              to = famind(sd+2)
!              do ifm=from,to
!                dst = fammbr(ifm)
                    dst = sd
                    icomreq = icomreq + 1
                    call MPI_Isend(poi(1, i*nxsd, j*nysd, max(ifloor*nzsd - szl, 0), 1), &
                   &               count, mptype_ss(ztype(1), 2, 1), &
                   &               dst, 0, MCW, ireqs(icomreq), mpierr)
!              end do
                end do
                end do
            end do
        end if
!
        call MPI_Waitall(icomreq, ireqs, istatus, mpierr)
!
!     ---------------
    else if (pftmode .eq. 2) then
        icomreq = 0
        do ps = 1, 1
            if (sdid(ps) .eq. -1) cycle
            xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
            yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
            zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)
            do src = int(yl/nyslc), int((yu - 1)/nyslc)
                count = 0
                do islice = src*nyslc, (src + 1)*nyslc - 1
                    if (yl .le. islice .and. islice .lt. yu) then
                        count = count + 1
                    end if
                end do
!            if(count.ne.nyslc) print*, "nyslc1", count, nyslc
                ytype(1) = 0
                if (src .eq. int(yl/nyslc)) then
                    count = count + 1
                    ytype(1) = ytype(1) - 1
                end if
                if (src .eq. int((yu - 1)/nyslc)) then
                    count = count + 2
                end if
                icomreq = icomreq + 1
                call MPI_Irecv(phi(1, 0, max(src*nyslc - yl, 0), 0, ps), &
               &               count, mptype_rc(ytype(1), 2, 2), &
               &               src, 0, MCW, ireqs(icomreq), mpierr)
            end do
        end do
!
        if (myid .lt. snode) then
            lfloor = int(syl/nysd); ufloor = int((syu - 1)/nysd)
            do ifloor = lfloor, ufloor
                count = min(syu, (ifloor + 1)*nysd) - max(syl, ifloor*nysd)
                ! if (count .ne. nyslc) print *, "nyslc2", count, nyslc
                ytype(1) = 0
                if (mod(max(syl, ifloor*nysd), nysd) .eq. 0) then
                    count = count + 1
                    ytype(1) = ytype(1) - 1
                end if
                if (mod(min(syu, (ifloor + 1)*nysd), nysd) .eq. 0) then
                    count = count + 2
                end if
                do k = 0, nodes(3) - 1
                do i = 0, nodes(1) - 1
                    sd = i + ifloor*nodes(1) + k*nodes(1)*nodes(2)
!              from = famind(sd+1) + 1
!              to = famind(sd+2)
!              do ifm=from,to
!                dst = fammbr(ifm)
                    dst = sd
                    icomreq = icomreq + 1
                    call MPI_Isend(poi(1, i*nxsd, max(ifloor*nysd - syl, 0), k*nzsd, 1), &
                   &               count, mptype_ss(ytype(1), 2, 2), &
                   &               dst, 0, MCW, ireqs(icomreq), mpierr)
!              end do
                end do
                end do
            end do
        end if
!
        call MPI_Waitall(icomreq, ireqs, istatus, mpierr)
!
!     ---------------
    else if (pftmode .eq. 3) then
        do ps = 1, 1
            if (sdid(ps) .eq. -1) cycle
            xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
            yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
            zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)
            phi(1, -1:xu - xl + 1, -1:ny + 1, -1:nz + 1, ps) = &
           &  poi(1, xl - 1:xu + 1, -1:ny + 1, -1:nz + 1, 1)
        end do
!
!     ---------------
    else if (pftmode .eq. 4) then
        do ps = 1, 1
            if (sdid(ps) .eq. -1) cycle
            xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
            yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
            zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)
            phi(1, -1:nx + 1, -1:ny + 1, -1:zu - zl + 1, ps) = &
           &  poi(1, -1:nx + 1, -1:ny + 1, zl - 1:zu + 1, 1)
        end do
    end if
    end if

    return
end subroutine poisson

!
subroutine poichk(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   P O I C H K
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
    integer(kind=4) :: i, j, k, ii, jj, kk, ncntl, ncnt
    integer(kind=4) :: xl, xu, yl, yu, zl, zu
    integer(kind=4) :: ierr
    integer(kind=4) :: ps
    real(kind=8) :: diffrhol, diffrho, drho

!--------------------
    xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
    yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
    zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)

!--------------------
    ncntl = 0
    diffrhol = 0.0d0
    do k = max(zl, lzfft(CP)), min(zu - 1, uzfft(CP))
    do j = max(yl, lyfft(CP)), min(yu - 1, uyfft(CP))
    do i = max(xl, lxfft(CP)), min(xu - 1, uxfft(CP))
        ii = i - xl; jj = j - yl; kk = k - zl
        drho = rho(1, ii, jj, kk, ps) + (phi(1, ii - 1, jj, kk, ps) + phi(1, ii + 1, jj, kk, ps) &
       &                           + phi(1, ii, jj - 1, kk, ps) + phi(1, ii, jj + 1, kk, ps) &
       &                           + phi(1, ii, jj, kk - 1, ps) + phi(1, ii, jj, kk + 1, ps)) &
       &                          - 6.0d0*phi(1, ii, jj, kk, ps)
        ncntl = ncntl + 1
        diffrhol = diffrhol + drho*drho
    end do
    end do
    end do

!--------------------
    call MPI_Reduce(ncntl, ncnt, 1, MPI_INTEGER, MPI_SUM, 0, MCW, ierr)
    call MPI_Reduce(diffrhol, diffrho, 1, MPI_REAL8, MPI_SUM, 0, MCW, ierr)

!--------------------
    if (myid .eq. 0) then
        diffrho = sqrt(diffrho)
        print *, "ncnt,diffrho =", ncnt, diffrho/ncnt
    end if

    return
end subroutine poichk
