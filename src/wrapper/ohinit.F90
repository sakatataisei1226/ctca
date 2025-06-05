#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"

subroutine ohinit
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   O H I N I T
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
#define MSS MPI_STATUS_SIZE
    implicit none

    integer(kind=8) :: npmax, nranmax
    integer(kind=4) :: n
    integer(kind=4) :: is

    integer(kind=4) :: inttmp(2, 3)
    integer(kind=HID_T) :: fileid
    integer(kind=4) :: stats0, stats1
    integer(kind=8) :: dims(2)
    character(len=30) :: filename, dsname

    integer :: i

    pcoord(:) = nodes(:)
    n = pcoord(1)*pcoord(2)*pcoord(3)
    minspec = max(nspec, 1)

    allocate (nphgram(n, minspec, 2))
    allocate (nphgram2(n, minspec, 2))
    allocate (sdoms(2, OH_DIMENSION, n))
    allocate (bounds(2, OH_DIMENSION, n))
    allocate (bared(2, OH_DIMENSION + 1, n))
    allocate (medges(2, OH_DIMENSION, n))
    allocate (famind(n + 1), fammbr(2*n))
    allocate (slcs(n, 3))
    allocate (scnts(n), rcnts(n))
    allocate (ncpmxs(n), scpmxs(n))

    allocate (totalp(minspec, 2))

    npmax = 0
    do is = 1, nspec
        npmax = npmax + np(is)
        if (myid .eq. 0) print *, "is,np =", is, np(is)
    end do
    npmax = max(npmax, 1000)
    if (myid .eq. 0) print *, "npmax =", npmax

    maxlocalp = oh_max_local_particles(npmax, MAXFRAC, OH_IPBUF_SIZE)
    if (myid .eq. 0) print *, "maxlocalp =", maxlocalp

    allocate (pbuf(maxlocalp))

    nranmax = 0
    do is = 1, nspec
        if (nranmax .lt. np(is)/nnode + 1) nranmax = np(is)/nnode + 1
    end do
    nranmax = max(nranmax*(100 + MAXFRAC)/100 + 1, 1000)
    allocate (dranu(nranmax))

    nbor(1, 1, 1) = -1
    sdoms(1, 1, 1) = 0; sdoms(2, 1, 1) = -1

    idxsd = mod(myid, nodes(1))
    idysd = mod(int(myid/nodes(1)), nodes(2))
    idzsd = int(myid/(nodes(1)*nodes(2)))

    scoord(1, 1) = 0; scoord(2, 1) = nx
    scoord(1, 2) = 0; scoord(2, 2) = ny
    scoord(1, 3) = 0; scoord(2, 3) = nz
    if (myid .eq. 0) print *, "scoord", scoord

    nbound = BNR + 1

    do i = 1, 3
        if (nfbnd(i) == 0) then
            bcond(:, i) = 1
        else
            bcond(:, i) = 2
        end if
    end do

    ! ftypes(:, *) = [ ndim, el, eu, elb, eub, elr, eur ]
    ! e: extention number
    ! e*b: for broadcast
    ! e*r: for reduce or allreduce
    ftypes(:, FEB) = (/6, -1, 2, -1, 2, 0, 0/)  ! for eb()
    ftypes(:, FDB) = (/3, -1, 2, -1, 2, 0, 0/)  ! for db()
    ftypes(:, FAJ) = (/3, 0, 0, 0, 0, -1, 2/)  ! for aj()
    ! ftypes(:, FRH) = (/1, 0, 0, 0, 0, 0, 1/)  ! for rho()
    ftypes(:, FRH) = (/1, 0, 1, 0, 0, 0, 1/)  ! for rho()
    ftypes(:, FPH) = (/1, 0, 0, -1, 2, 0, 0/)  ! for phi()
    ! ftypes(:, FRD) = (/minspec*2 + 1, 0, 0, 0, 0, 0, 1/)  ! for rhodg()
    ftypes(:, FRD) = (/minspec*2 + 1, 0, 1, 0, 0, 0, 1/)  ! for rhodg()
    ftypes(:, FJD) = (/minspec*3, 0, 0, 0, 0, -1, 2/)  ! for ajdg()
    ftypes(1, FNR + 1) = 0  ! terminator

    cfields(:) = (/FEB, FAJ, FRH, FRD, FJD, CPH, 0/)
    ! ctypes(:,:,*,***) = reshape((/downward,  upward/), (/3,2)))  ! sample
    ctypes(:, :, 1, CEB) = reshape((/0, 0, 2, -1, -1, 1/), (/3, 2/))  ! for eb()
    ctypes(:, :, 2, CEB) = reshape((/0, 0, 2, -1, -1, 1/), (/3, 2/))  ! for eb()
    ctypes(:, :, 1, CAJ) = reshape((/-1, 2, 3, -1, -4, 3/), (/3, 2/))  ! for aj()
    ctypes(:, :, 2, CAJ) = reshape((/-1, 2, 3, -1, -4, 3/), (/3, 2/))  ! for aj()
    ctypes(:, :, 1, CRH) = reshape([[0, 1, 1], [0, -1, 1]], (/3, 2/))  ! for rho()
    ctypes(:, :, 2, CRH) = reshape([[0, 1, 1], [0, -1, 1]], (/3, 2/))  ! for rho()
    ctypes(:, :, 1, CRD) = reshape([[0, 1, 1], [0, -1, 1]], (/3, 2/))  ! for rhodg()
    ctypes(:, :, 2, CRD) = reshape([[0, 1, 1], [0, -1, 1]], (/3, 2/))  ! for rhodg()
    ctypes(:, :, 1, CJD) = reshape((/-1, 2, 3, -1, -4, 3/), (/3, 2/))  ! for ajdg()
    ctypes(:, :, 2, CJD) = reshape((/-1, 2, 3, -1, -4, 3/), (/3, 2/))  ! for ajdg()
    ctypes(:, :, 1, CPH) = reshape([[0, 0, 1], [-1, -1, 1]], [3, 2])  ! for phi()
    ctypes(:, :, 2, CPH) = reshape([[0, 0, 1], [-1, -1, 1]], [3, 2])  ! for phi()

    call oh3_init(sdid(:), minspec, MAXFRAC, nphgram(:, :, :), totalp(:, :), &
                  pbuf(:), pbase(:), maxlocalp, mycomm, nbor(:, :, :), &
                  pcoord(:), sdoms(:, :, :), scoord(:, :), nbound, bcond(:, :), &
                  bounds(:, :, :), ftypes(:, :), cfields(:), ctypes(:, :, :, :), &
                  fsizes(:, :, :), OHHELP_stats, OHHELP_repiter, &
                  OHHELP_verbose)

    allocate (eb(6, &
                 fsizes(1, 1, FEB):fsizes(2, 1, FEB), &
                 fsizes(1, 2, FEB):fsizes(2, 2, FEB), &
                 fsizes(1, 3, FEB):fsizes(2, 3, FEB), &
                 4))
    allocate (ebav(6, &
                   fsizes(1, 1, FEB):fsizes(2, 1, FEB), &
                   fsizes(1, 2, FEB):fsizes(2, 2, FEB), &
                   fsizes(1, 3, FEB):fsizes(2, 3, FEB)))
    allocate (mp(6, &
                 fsizes(1, 1, FEB):fsizes(2, 1, FEB), &
                 fsizes(1, 2, FEB):fsizes(2, 2, FEB), &
                 fsizes(1, 3, FEB):fsizes(2, 3, FEB), &
                 4))
    allocate (db(3, &
                 fsizes(1, 1, FDB):fsizes(2, 1, FDB), &
                 fsizes(1, 2, FDB):fsizes(2, 2, FDB), &
                 fsizes(1, 3, FDB):fsizes(2, 3, FDB), &
                 2))
    allocate (aj(3, &
                 fsizes(1, 1, FAJ):fsizes(2, 1, FAJ), &
                 fsizes(1, 2, FAJ):fsizes(2, 2, FAJ), &
                 fsizes(1, 3, FAJ):fsizes(2, 3, FAJ), &
                 2))
    allocate (rho(1, &
                  fsizes(1, 1, FRH):fsizes(2, 1, FRH), &
                  fsizes(1, 2, FRH):fsizes(2, 2, FRH), &
                  fsizes(1, 3, FRH):fsizes(2, 3, FRH), &
                  2))
    allocate (rho_tilde(1, &
                        fsizes(1, 1, FRH):fsizes(2, 1, FRH), &
                        fsizes(1, 2, FRH):fsizes(2, 2, FRH), &
                        fsizes(1, 3, FRH):fsizes(2, 3, FRH), &
                        1))
    allocate (phi(1, &
                  fsizes(1, 1, FPH):fsizes(2, 1, FPH), &
                  fsizes(1, 2, FPH):fsizes(2, 2, FPH), &
                  fsizes(1, 3, FPH):fsizes(2, 3, FPH), &
                  2))
    allocate (phiav(1, &
                    fsizes(1, 1, FPH):fsizes(2, 1, FPH), &
                    fsizes(1, 2, FPH):fsizes(2, 2, FPH), &
                    fsizes(1, 3, FPH):fsizes(2, 3, FPH)))
    allocate (wrk(9, &
                  fsizes(1, 1, FAJ):fsizes(2, 1, FAJ), &
                  fsizes(1, 2, FAJ):fsizes(2, 2, FAJ), &
                  fsizes(1, 3, FAJ):fsizes(2, 3, FAJ)))
    allocate (rhobk(1, &
                    fsizes(1, 1, FRH):fsizes(2, 1, FRH), &
                    fsizes(1, 2, FRH):fsizes(2, 2, FRH), &
                    fsizes(1, 3, FRH):fsizes(2, 3, FRH), &
                    3))
    allocate (rhobksp(1, &
                      fsizes(1, 1, FRH):fsizes(2, 1, FRH), &
                      fsizes(1, 2, FRH):fsizes(2, 2, FRH), &
                      fsizes(1, 3, FRH):fsizes(2, 3, FRH), &
                      3, &
                      minspec))
    allocate (rhodg(minspec*2 + 1, &
                    fsizes(1, 1, FRD):fsizes(2, 1, FRD), &
                    fsizes(1, 2, FRD):fsizes(2, 2, FRD), &
                    fsizes(1, 3, FRD):fsizes(2, 3, FRD), &
                    2))
    allocate (rhoav(minspec*2 + 3, &
                    fsizes(1, 1, FRD):fsizes(2, 1, FRD), &
                    fsizes(1, 2, FRD):fsizes(2, 2, FRD), &
                    fsizes(1, 3, FRD):fsizes(2, 3, FRD)))
    allocate (ajdg(minspec*3, &
                   fsizes(1, 1, FJD):fsizes(2, 1, FJD), &
                   fsizes(1, 2, FJD):fsizes(2, 2, FJD), &
                   fsizes(1, 3, FJD):fsizes(2, 3, FJD), &
                   2))
    allocate (ajav(minspec*3 + 3, &
                   fsizes(1, 1, FJD):fsizes(2, 1, FJD), &
                   fsizes(1, 2, FJD):fsizes(2, 2, FJD), &
                   fsizes(1, 3, FJD):fsizes(2, 3, FJD)))
    allocate (colf(1, &
                   fsizes(1, 1, FPH):fsizes(2, 1, FPH), &
                   fsizes(1, 2, FPH):fsizes(2, 2, FPH), &
                   fsizes(1, 3, FPH):fsizes(2, 3, FPH), &
                   2))

    ! sub-domain/array dimensions
    nxsd = sdoms(2, 1, sdid(1) + 1) - sdoms(1, 1, sdid(1) + 1)
    nysd = sdoms(2, 2, sdid(1) + 1) - sdoms(1, 2, sdid(1) + 1)
    nzsd = sdoms(2, 3, sdid(1) + 1) - sdoms(1, 3, sdid(1) + 1)

    do n = 1, FNR
        lxsd(n) = fsizes(2, 1, n) - fsizes(1, 1, n) + 1
        lysd(n) = fsizes(2, 2, n) - fsizes(1, 2, n) + 1
        lzsd(n) = fsizes(2, 3, n) - fsizes(1, 3, n) + 1
    end do

    ! Load previous simulation settings to continue.
    if (jobnum(1) .gt. 0) then
        if (jobnum(1) .eq. 1) then
            write (filename, '(a,i4.4,a)') './SNAPSHOT0/esdat', myid, '.h5'
        else
            write (filename, '(a,i4.4,a)') './SNAPSHOT0/emdat', myid, '.h5'
        end if
        call hdfopen(filename, fileid, DFACC_READ)

        dsname = 'sdoms'
        dims(1) = 2; dims(2) = 3
        call read2i(fileid, dsname, dims(1:2), inttmp(1:2, 1:3), stats0, stats1)
        if (inttmp(1, 1) .ne. sdoms(1, 1, 1) .or. inttmp(2, 1) .ne. sdoms(2, 1, 1) .or. &
            inttmp(1, 2) .ne. sdoms(1, 2, 1) .or. inttmp(2, 2) .ne. sdoms(2, 2, 1) .or. &
            inttmp(1, 3) .ne. sdoms(1, 3, 1) .or. inttmp(2, 3) .ne. sdoms(2, 3, 1)) then
            if (myid .eq. 0) print *, "sdoms is inconsistent with continued job data: STOP"
            stop
        end if

        call hdfclose(fileid, stats0)
    end if

    if (myid .eq. 0) print *, 'maxlocalp=', maxlocalp

    if (myid .eq. 0) print *, 'EB:'
    if (myid .eq. 0) print *, '  ', fsizes(1, 1, FEB), fsizes(2, 1, FEB), lxsd(FEB)
    if (myid .eq. 0) print *, '  ', fsizes(1, 2, FEB), fsizes(2, 2, FEB), lysd(FEB)
    if (myid .eq. 0) print *, '  ', fsizes(1, 3, FEB), fsizes(2, 3, FEB), lzsd(FEB)
    if (myid .eq. 0) print *, 'DB:'
    if (myid .eq. 0) print *, '  ', fsizes(1, 1, FDB), fsizes(2, 1, FDB), lxsd(FDB)
    if (myid .eq. 0) print *, '  ', fsizes(1, 2, FDB), fsizes(2, 2, FDB), lysd(FDB)
    if (myid .eq. 0) print *, '  ', fsizes(1, 3, FDB), fsizes(2, 3, FDB), lzsd(FDB)
    if (myid .eq. 0) print *, 'AJ:'
    if (myid .eq. 0) print *, '  ', fsizes(1, 1, FAJ), fsizes(2, 1, FAJ), lxsd(FAJ)
    if (myid .eq. 0) print *, '  ', fsizes(1, 2, FAJ), fsizes(2, 2, FAJ), lysd(FAJ)
    if (myid .eq. 0) print *, '  ', fsizes(1, 3, FAJ), fsizes(2, 3, FAJ), lzsd(FAJ)
    if (myid .eq. 0) print *, 'RH:'
    if (myid .eq. 0) print *, '  ', fsizes(1, 1, FRH), fsizes(2, 1, FRH), lxsd(FRH)
    if (myid .eq. 0) print *, '  ', fsizes(1, 2, FRH), fsizes(2, 2, FRH), lysd(FRH)
    if (myid .eq. 0) print *, '  ', fsizes(1, 3, FRH), fsizes(2, 3, FRH), lzsd(FRH)
    if (myid .eq. 0) print *, 'PH:'
    if (myid .eq. 0) print *, '  ', fsizes(1, 1, FPH), fsizes(2, 1, FPH), lxsd(FPH)
    if (myid .eq. 0) print *, '  ', fsizes(1, 2, FPH), fsizes(2, 2, FPH), lysd(FPH)
    if (myid .eq. 0) print *, '  ', fsizes(1, 3, FPH), fsizes(2, 3, FPH), lzsd(FPH)
    if (myid .eq. 0) print *, 'RD:'
    if (myid .eq. 0) print *, '  ', fsizes(1, 1, FRD), fsizes(2, 1, FRD), lxsd(FRD)
    if (myid .eq. 0) print *, '  ', fsizes(1, 2, FRD), fsizes(2, 2, FRD), lysd(FRD)
    if (myid .eq. 0) print *, '  ', fsizes(1, 3, FRD), fsizes(2, 3, FRD), lzsd(FRD)
    if (myid .eq. 0) print *, 'JD:'
    if (myid .eq. 0) print *, '  ', fsizes(1, 1, FJD), fsizes(2, 1, FJD), lxsd(FJD)
    if (myid .eq. 0) print *, '  ', fsizes(1, 2, FJD), fsizes(2, 2, FJD), lysd(FJD)
    if (myid .eq. 0) print *, '  ', fsizes(1, 3, FJD), fsizes(2, 3, FJD), lzsd(FJD)

    ! if(myid.eq.0) print*,'bounds:'
    ! if(myid.eq.0) print*,'  ',bounds
end subroutine
