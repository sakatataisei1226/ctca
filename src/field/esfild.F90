#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
#include "oh_stats.h"

  subroutine esfld1(func)
!   ____________________________________________________________
!
!               S U B R O U T I N E   E S F I L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   this subroutine gives a solution of electrostatic      .
!   .   by solving the Poisson's equation.                     .
!   ............................................................
!
!-------------------- parameter and common blocks
      use oh_type
      use paramt
      use allcom
      use m_logging
#define MCW local_comm
      implicit none
!
      integer(kind=4) :: i, j, k, i1, j1, k1, ii, jj, kk, ipc, icap
      integer(kind=4) :: xl, xu, yl, yu, zl, zu
      integer(kind=4) :: siu, sju, sku
      integer(kind=4) :: xlower, xupper, ylower, yupper, zlower, zupper
      integer(kind=4) :: func
      integer(kind=4) :: ierr
      real(kind=8) :: xlocal, ylocal, zlocal
      real(kind=8) :: x1, y1, z1, z2, xy1, xz1, yz1, xz2, yz2
      real(kind=8) :: v1, v2, v3, v4, v5, v6, v7, v8

!--------------------
      xl = sdoms(1, 1, sdid(1) + 1); xu = sdoms(2, 1, sdid(1) + 1)
      yl = sdoms(1, 2, sdid(1) + 1); yu = sdoms(2, 2, sdid(1) + 1)
      zl = sdoms(1, 3, sdid(1) + 1); zu = sdoms(2, 3, sdid(1) + 1)
      if (pftmode .eq. 3 .and. myid .ne. snode - 1) then
          siu = sxu - sxl - 1
      else
          siu = sxu - sxl
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

!-------------------- zero clear
      phi(:, :, :, :, :) = 0.0d0
      poi(:, :, :, :, :) = 0.0d0

!-------------------- Poisson Solver
      if (func .eq. 0) then
          call logging_debuglog('Start call poisson(3)')
          call poisson(3)
          call logging_debuglog('End call poisson(3)')
      else if (func .eq. 1) then
          do icap = 1, ncpmx
              xlocal = bdygrid(1, icap) - sxl
              ylocal = bdygrid(2, icap) - syl
              zlocal = bdygrid(3, icap) - szl

              i = floor(xlocal)
              j = floor(ylocal)
              k = floor(zlocal)

              i1 = i + 1
              j1 = j + 1
              k1 = k + 1

              x1 = xlocal - i
              y1 = ylocal - j
              z1 = zlocal - k

              xy1 = x1*y1
              xz1 = x1*z1
              yz1 = y1*z1
              z2 = 1.0d0 - z1
              xz2 = x1*z2
              yz2 = y1*z2

              v3 = xy1*z1
              v2 = xz1 - v3
              v4 = yz1 - v3
              v1 = z1 - xz1 - v4

              v7 = xy1*z2
              v6 = xz2 - v7
              v8 = yz2 - v7
              v5 = z2 - xz2 - v8
!         ----------- distribute charge with weight
              if (i .ge. 0 .and. i .le. siu .and. &
             &   j .ge. 0 .and. j .le. sju .and. &
             &   k1 .ge. 0 .and. k1 .le. sku) then
                  poi(1, i, j, k1, 1) = poi(1, i, j, k1, 1) + v1*sfrho(icap)
              end if
              if (i1 .ge. 0 .and. i1 .le. siu .and. &
             &   j .ge. 0 .and. j .le. sju .and. &
             &   k1 .ge. 0 .and. k1 .le. sku) then
                  poi(1, i1, j, k1, 1) = poi(1, i1, j, k1, 1) + v2*sfrho(icap)
              end if
              if (i1 .ge. 0 .and. i1 .le. siu .and. &
             &   j1 .ge. 0 .and. j1 .le. sju .and. &
             &   k1 .ge. 0 .and. k1 .le. sku) then
                  poi(1, i1, j1, k1, 1) = poi(1, i1, j1, k1, 1) + v3*sfrho(icap)
              end if
              if (i .ge. 0 .and. i .le. siu .and. &
             &   j1 .ge. 0 .and. j1 .le. sju .and. &
             &   k1 .ge. 0 .and. k1 .le. sku) then
                  poi(1, i, j1, k1, 1) = poi(1, i, j1, k1, 1) + v4*sfrho(icap)
              end if
              if (i .ge. 0 .and. i .le. siu .and. &
             &   j .ge. 0 .and. j .le. sju .and. &
             &   k .ge. 0 .and. k .le. sku) then
                  poi(1, i, j, k, 1) = poi(1, i, j, k, 1) + v5*sfrho(icap)
              end if
              if (i1 .ge. 0 .and. i1 .le. siu .and. &
             &   j .ge. 0 .and. j .le. sju .and. &
             &   k .ge. 0 .and. k .le. sku) then
                  poi(1, i1, j, k, 1) = poi(1, i1, j, k, 1) + v6*sfrho(icap)
              end if
              if (i1 .ge. 0 .and. i1 .le. siu .and. &
             &   j1 .ge. 0 .and. j1 .le. sju .and. &
             &   k .ge. 0 .and. k .le. sku) then
                  poi(1, i1, j1, k, 1) = poi(1, i1, j1, k, 1) + v7*sfrho(icap)
              end if
              if (i .ge. 0 .and. i .le. siu .and. &
             &   j1 .ge. 0 .and. j1 .le. sju .and. &
             &   k .ge. 0 .and. k .le. sku) then
                  poi(1, i, j1, k, 1) = poi(1, i, j1, k, 1) + v8*sfrho(icap)
              end if
          end do
          call logging_debuglog('Start call poisson(2)')
          call poisson(2)
          call logging_debuglog('End call poisson(2)')
      end if

!--------------------
      ngref(1) = 0
      phiref(1) = 0.0d0
      if (mtd_vbnd(1) .ne. 1 .or. mtd_vbnd(2) .ne. 1 .or. mtd_vbnd(3) .ne. 1) then
          xlower = 0
          if (bared(2, 1, myid + 1) .eq. 1) then
              xupper = xu - xl
          else
              xupper = xu - xl - 1
          end if

          ylower = 0
          if (bared(2, 2, myid + 1) .eq. 1) then
              yupper = yu - yl
          else
              yupper = yu - yl - 1
          end if

          zlower = 0
          if (bared(2, 3, myid + 1) .eq. 1) then
              zupper = zu - zl
          else
              zupper = zu - zl - 1
          end if

          if (iphiref(1, 1) .eq. 1 .and. bared(1, 1, myid + 1) .eq. 1) then
              do k = zlower, zupper
              do j = ylower, yupper
                  ngref(1) = ngref(1) + 1
                  phiref(1) = phiref(1) + phi(1, 0, j, k, 1)
              end do
              end do
              xlower = 1
          end if

          if (iphiref(2, 1) .eq. 1 .and. bared(2, 1, myid + 1) .eq. 1) then
              do k = zlower, zupper
              do j = ylower, yupper
                  ngref(1) = ngref(1) + 1
                  phiref(1) = phiref(1) + phi(1, xu - xl, j, k, 1)
              end do
              end do
              xupper = xu - xl - 1
          end if

          if (iphiref(1, 2) .eq. 1 .and. bared(1, 2, myid + 1) .eq. 1) then
              do k = zlower, zupper
              do i = xlower, xupper
                  ngref(1) = ngref(1) + 1
                  phiref(1) = phiref(1) + phi(1, i, 0, k, 1)
              end do
              end do
              ylower = 1
          end if

          if (iphiref(2, 2) .eq. 1 .and. bared(2, 2, myid + 1) .eq. 1) then
              do k = zlower, zupper
              do i = xlower, xupper
                  ngref(1) = ngref(1) + 1
                  phiref(1) = phiref(1) + phi(1, i, yu - yl, k, 1)
              end do
              end do
              yupper = yu - yl - 1
          end if

          if (iphiref(1, 3) .eq. 1 .and. bared(1, 3, myid + 1) .eq. 1) then
              do j = ylower, yupper
              do i = xlower, xupper
                  ngref(1) = ngref(1) + 1
                  phiref(1) = phiref(1) + phi(1, i, j, 0, 1)
              end do
              end do
          end if

          if (iphiref(2, 3) .eq. 1 .and. bared(2, 3, myid + 1) .eq. 1) then
              do j = ylower, yupper
              do i = xlower, xupper
                  ngref(1) = ngref(1) + 1
                  phiref(1) = phiref(1) + phi(1, i, j, zu - zl, 1)
              end do
              end do
          end if

          call logging_debuglog('Start call MPI_Reduce in esfld1')
          call MPI_Reduce(ngref(1), ngref(2), 1, MPI_INTEGER, MPI_SUM, 0, MCW, ierr)
          call MPI_Reduce(phiref(1), phiref(2), 1, MPI_REAL8, MPI_SUM, 0, MCW, ierr)
          call logging_debuglog('End call MPI_Reduce in esfld1')

          if (myid .eq. 0 .and. ngref(2) .ne. 0) then
              phiref(2) = phiref(2)/ngref(2)
          else
              phiref(2) = 0.0d0
          end if
      end if

      return
  end subroutine

  subroutine esfld2(ps)

!   ____________________________________________________________
!
!               S U B R O U T I N E   E S F I L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   this subroutine gives a solution of electrostatic      .
!   .   by solving the Poisson's equation.                     .
!   ............................................................
!
!-------------------- parameter and common blocks
      use oh_type
      use paramt
      use allcom
      implicit none
!
      integer(kind=4) :: i, j, k, i1, j1, k1
      integer(kind=4) :: ipc, icap
      integer(kind=4) :: xl, yl, zl, xu, yu, zu, lx, ly, lz
      integer(kind=4) :: ps
      real(kind=8) :: xlocal, ylocal, zlocal
      real(kind=8) :: x1, y1, z1, x2, y2, z2, etmp
      real(kind=8) :: vfactor, disp, r1, r2

      xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
      yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
      zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)
      lx = xu - xl
      ly = yu - yl
      lz = zu - zl

!-------------------- smoothing potentials inside bodies
      if (sfecrrct .eq. 2) then
          do ipc = 1, npc
              if (boom(ipc)%align .eq. 1) then
                  vfactor = 0.5d0/log(2.0d0*boom(ipc)%hlength/boom(ipc)%rradius)
                  do icap = nscpmx(ipc), nmxcpmx(ipc) - 1
                      disp = bdygrid(1, icap) - boom(ipc)%origin(1)
                      r1 = sqrt((disp - boom(ipc)%hlength)*(disp - boom(ipc)%hlength) &
                     &          + boom(ipc)%eradius*boom(ipc)%eradius)
                      r2 = sqrt((disp + boom(ipc)%hlength)*(disp + boom(ipc)%hlength) &
                     &          + boom(ipc)%eradius*boom(ipc)%eradius)
                      i = bdygrid(1, icap) - xl
                      j = bdygrid(2, icap) - yl
                      k = bdygrid(3, icap) - zl
                      if (i .ge. -1 .and. i .le. lx + 1 .and. &
                     &   j .ge. -1 .and. j .le. ly + 1 .and. &
                     &   k .ge. -1 .and. k .le. lz + 1) then
                          phi(1, i, j, k, ps) = selfp(ipc)*vfactor &
                       &            *log((+boom(ipc)%hlength - disp + r1) &
                       &                /(-boom(ipc)%hlength - disp + r2))
                      end if
                  end do
              else if (boom(ipc)%align .eq. 2) then
                  vfactor = 0.5d0/log(2.0d0*boom(ipc)%hlength/boom(ipc)%rradius)
                  do icap = nscpmx(ipc), nmxcpmx(ipc) - 1
                      disp = bdygrid(2, icap) - boom(ipc)%origin(2)
                      r1 = sqrt((disp - boom(ipc)%hlength)*(disp - boom(ipc)%hlength) &
                     &          + boom(ipc)%eradius*boom(ipc)%eradius)
                      r2 = sqrt((disp + boom(ipc)%hlength)*(disp + boom(ipc)%hlength) &
                     &          + boom(ipc)%eradius*boom(ipc)%eradius)
                      i = bdygrid(1, icap) - xl
                      j = bdygrid(2, icap) - yl
                      k = bdygrid(3, icap) - zl
                      if (i .ge. -1 .and. i .le. lx + 1 .and. &
                     &   j .ge. -1 .and. j .le. ly + 1 .and. &
                     &   k .ge. -1 .and. k .le. lz + 1) then
                          phi(1, i, j, k, ps) = selfp(ipc)*vfactor &
                       &            *log((+boom(ipc)%hlength - disp + r1) &
                       &                /(-boom(ipc)%hlength - disp + r2))
                      end if
                  end do
              else if (boom(ipc)%align .eq. 3) then
                  vfactor = 0.5d0/log(2.0d0*boom(ipc)%hlength/boom(ipc)%rradius)
                  do icap = nscpmx(ipc), nmxcpmx(ipc) - 1
                      disp = bdygrid(3, icap) - boom(ipc)%origin(3)
                      r1 = sqrt((disp - boom(ipc)%hlength)*(disp - boom(ipc)%hlength) &
                     &          + boom(ipc)%eradius*boom(ipc)%eradius)
                      r2 = sqrt((disp + boom(ipc)%hlength)*(disp + boom(ipc)%hlength) &
                     &          + boom(ipc)%eradius*boom(ipc)%eradius)
                      i = bdygrid(1, icap) - xl
                      j = bdygrid(2, icap) - yl
                      k = bdygrid(3, icap) - zl
                      if (i .ge. -1 .and. i .le. lx + 1 .and. &
                     &   j .ge. -1 .and. j .le. ly + 1 .and. &
                     &   k .ge. -1 .and. k .le. lz + 1) then
                          phi(1, i, j, k, ps) = selfp(ipc)*vfactor &
                       &            *log((+boom(ipc)%hlength - disp + r1) &
                       &                /(-boom(ipc)%hlength - disp + r2))
                      end if
                  end do
              else
                  do icap = nscpmx(ipc), nmxcpmx(ipc) - 1
                      i = bdygrid(1, icap) - xl
                      j = bdygrid(2, icap) - yl
                      k = bdygrid(3, icap) - zl
                      if (i .ge. -1 .and. i .le. lx + 1 .and. &
                     &   j .ge. -1 .and. j .le. ly + 1 .and. &
                     &   k .ge. -1 .and. k .le. lz + 1) then
                          phi(1, i, j, k, ps) = selfp(ipc)
                      end if
                  end do
              end if
          end do
      end if

!-------------------- zero clear
      eb(EX:EZ, :, :, :, ps) = 0.0d0

!-------------------- electrostatic field
      eb(EX, -1:xu - xl, 0:yu - yl, 0:zu - zl, ps) = &
          (phi(1, -1:xu - xl, 0:yu - yl, 0:zu - zl, ps) &
           - phi(1, 0:xu - xl + 1, 0:yu - yl, 0:zu - zl, ps))*mltstp

      eb(EY, 0:xu - xl, -1:yu - yl, 0:zu - zl, ps) = &
          (phi(1, 0:xu - xl, -1:yu - yl, 0:zu - zl, ps) &
           - phi(1, 0:xu - xl, 0:yu - yl + 1, 0:zu - zl, ps))*mltstp

      eb(EZ, 0:xu - xl, 0:yu - yl, -1:zu - zl, ps) = &
          (phi(1, 0:xu - xl, 0:yu - yl, -1:zu - zl, ps) &
           - phi(1, 0:xu - xl, 0:yu - yl, 0:zu - zl + 1, ps))*mltstp

      if (xl == 0) then
          i = -1
          if (mtd_vbnd(1) == 0) then
              ! do nothing
          else if (mtd_vbnd(1) == 1) then
              eb(EX, i, 0:yu - yl, 0:zu - zl, ps) = eb(EX, 0, 0:yu - yl, 0:zu - zl, ps)
          else
              eb(EX, i, 0:yu - yl, 0:zu - zl, ps) = -eb(EX, 0, 0:yu - yl, 0:zu - zl, ps)
          end if
      end if

      if (xu == nx) then
          i = xu - xl
          if (mtd_vbnd(1) == 0) then
              ! do nothing
          else if (mtd_vbnd(1) == 1) then
              eb(EX, i, 0:yu - yl, 0:zu - zl, ps) = eb(EX, i - 1, 0:yu - yl, 0:zu - zl, ps)
          else
              eb(EX, i, 0:yu - yl, 0:zu - zl, ps) = -eb(EX, i - 1, 0:yu - yl, 0:zu - zl, ps)
          end if
      end if

      if (yl == 0) then
          j = -1
          if (mtd_vbnd(2) == 0) then
              ! do nothing
          else if (mtd_vbnd(2) == 1) then
              eb(EY, 0:xu - xl, j, 0:zu - zl, ps) = eb(EY, 0:xu - xl, 0, 0:zu - zl, ps)
          else
              eb(EY, 0:xu - xl, j, 0:zu - zl, ps) = -eb(EY, 0:xu - xl, 0, 0:zu - zl, ps)
          end if
      end if

      if (yu == ny) then
          j = yu - yl
          if (mtd_vbnd(2) == 0) then
              ! do nothing
          else if (mtd_vbnd(2) == 1) then
              eb(EY, 0:xu - xl, j, 0:zu - zl, ps) = eb(EY, 0:xu - xl, j - 1, 0:zu - zl, ps)
          else
              eb(EY, 0:xu - xl, j, 0:zu - zl, ps) = -eb(EY, 0:xu - xl, j - 1, 0:zu - zl, ps)
          end if
      end if

      if (zl == 0) then
          k = -1
          if (mtd_vbnd(3) == 0) then
              ! do nothing
          else if (mtd_vbnd(3) == 1) then
              eb(EZ, 0:xu - xl, 0:yu - yl, k, ps) = eb(EZ, 0:xu - xl, 0:yu - yl, 0, ps)
          else
              eb(EZ, 0:xu - xl, 0:yu - yl, k, ps) = -eb(EZ, 0:xu - xl, 0:yu - yl, 0, ps)
          end if
      end if

      if (zu == nz) then
          k = zu - zl
          if (mtd_vbnd(3) == 0) then
              ! do nothing
          else if (mtd_vbnd(3) == 1) then
              eb(EZ, 0:xu - xl, 0:yu - yl, k, ps) = eb(EZ, 0:xu - xl, 0:yu - yl, k - 1, ps)
          else
              eb(EZ, 0:xu - xl, 0:yu - yl, k, ps) = -eb(EZ, 0:xu - xl, 0:yu - yl, k - 1, ps)
          end if
      end if

!-------------------- smoothing potentials inside bodies
      if (sfecrrct .eq. 2) then
          do ipc = 1, npc
              do icap = nmxcpmx(ipc), nmycpmx(ipc) - 1
                  xlocal = bdygrid(1, icap) - xl
                  ylocal = bdygrid(2, icap) - yl
                  zlocal = bdygrid(3, icap) - zl

                  if (xlocal .ge. -1 .and. xlocal .le. lx .and. &
                 &   ylocal .ge. 0 .and. ylocal .le. ly .and. &
                 &   zlocal .ge. 0 .and. zlocal .le. lz) then
                      i = floor(xlocal)
                      j = floor(ylocal)
                      k = floor(zlocal)

                      i1 = i + 1
                      x1 = xlocal - i
                      x2 = 1.0d0 - x1

                      if (phi(1, i, j, k, ps) .eq. selfp(ipc)) then
                          eb(EX, i, j, k, ps) = eb(EX, i, j, k, ps)/x2
                      end if

                      if (phi(1, i1, j, k, ps) .eq. selfp(ipc)) then
                          eb(EX, i, j, k, ps) = eb(EX, i, j, k, ps)/x1
                      end if
                  end if
              end do

              do icap = nmycpmx(ipc), nmzcpmx(ipc) - 1
                  xlocal = bdygrid(1, icap) - xl
                  ylocal = bdygrid(2, icap) - yl
                  zlocal = bdygrid(3, icap) - zl

                  if (xlocal .ge. 0 .and. xlocal .le. lx .and. &
                 &   ylocal .ge. -1 .and. ylocal .le. ly .and. &
                 &   zlocal .ge. 0 .and. zlocal .le. lz) then
                      i = floor(xlocal)
                      j = floor(ylocal)
                      k = floor(zlocal)

                      j1 = j + 1
                      y1 = ylocal - j
                      y2 = 1.0d0 - y1

                      if (phi(1, i, j, k, ps) .eq. selfp(ipc)) then
                          eb(EY, i, j, k, ps) = eb(EY, i, j, k, ps)/y2
                      end if

                      if (phi(1, i, j1, k, ps) .eq. selfp(ipc)) then
                          eb(EY, i, j, k, ps) = eb(EY, i, j, k, ps)/y1
                      end if
                  end if
              end do

              do icap = nmzcpmx(ipc), nmxycpmx(ipc) - 1
                  xlocal = bdygrid(1, icap) - xl
                  ylocal = bdygrid(2, icap) - yl
                  zlocal = bdygrid(3, icap) - zl

                  if (xlocal .ge. 0 .and. xlocal .le. lx .and. &
                 &   ylocal .ge. 0 .and. ylocal .le. ly .and. &
                 &   zlocal .ge. -1 .and. zlocal .le. lz) then
                      i = floor(xlocal)
                      j = floor(ylocal)
                      k = floor(zlocal)

                      k1 = k + 1
                      z1 = zlocal - k
                      z2 = 1.0d0 - z1

                      if (phi(1, i, j, k, ps) .eq. selfp(ipc)) then
                          eb(EZ, i, j, k, ps) = eb(EZ, i, j, k, ps)/z2
                      end if

                      if (phi(1, i, j, k1, ps) .eq. selfp(ipc)) then
                          eb(EZ, i, j, k, ps) = eb(EZ, i, j, k, ps)/z1
                      end if
                  end if
              end do
          end do
      end if

      return
  end subroutine

!
  subroutine esfld3(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   E S F I L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   this subroutine gives a solution of electrostatic      .
!   .   by solving the Poisson's equation.                     .
!   ............................................................
!
!-------------------- parameter and common blocks
      use oh_type
      use paramt
      use allcom
      implicit none
!
      integer(kind=4) :: i, j, k
      integer(kind=4) :: xl, yl, zl, xu, yu, zu
      integer(kind=4) :: ps

!--------------------
      xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
      yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
      zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)

!-------------------- electrostatic field
      do k = -1, zu - zl
      do j = -1, yu - yl
      do i = -1, xu - xl
          if (j /= -1 .and. k /= -1) then
              eb(EX, i, j, k, ps) = eb(EX, i, j, k, ps) &
                                    + (phi(1, i, j, k, ps) - phi(1, i + 1, j, k, ps))*mltstp
          end if
          if (k /= -1 .and. i /= -1) then
              eb(EY, i, j, k, ps) = eb(EY, i, j, k, ps) &
                                    + (phi(1, i, j, k, ps) - phi(1, i, j + 1, k, ps))*mltstp
          end if
          if (i /= -1 .and. j /= -1) then
              eb(EZ, i, j, k, ps) = eb(EZ, i, j, k, ps) &
                                    + (phi(1, i, j, k, ps) - phi(1, i, j, k + 1, ps))*mltstp
          end if
      end do
      end do
      end do
  end subroutine
