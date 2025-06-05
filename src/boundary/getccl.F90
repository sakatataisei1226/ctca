#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine getccl
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   G E T C C L
!   ____________________________________________________________
!
!-------------------- parameter and common blocks
      use oh_type
      use paramt
      use allcom
      implicit none
!
      interface rfold
          integer(kind=4) function rfold_i(r, nr, dir)
              integer(kind=4), intent(in) :: r, nr
              integer(kind=4), intent(in) :: dir
          end function rfold_i
!
          real(kind=4) function rfold_r(r, nr, dir)
              real(kind=4), intent(in) :: r, nr
              integer(kind=4), intent(in) :: dir
          end function rfold_r
!
          real(kind=8) function rfold_d(r, nr, dir)
              real(kind=8), intent(in) :: r, nr
              integer(kind=4), intent(in) :: dir
          end function rfold_d
      end interface rfold
!
      integer(kind=4) :: i, j, k, i0, j0, k0
      integer(kind=4) :: icap, jcap, iicap, ibdy
      integer(kind=4) :: ipc, islice
      integer(kind=4) :: ibedge, jbedge, kbedge
      integer(kind=4) :: maxbedge
      real(kind=8) :: x, y, z, x1, y1, z1
      real(kind=8) :: cylaln, cylrad, cylaxs1, cylaxs2, cyledg1, cyledg2
      real(kind=8) :: sphrad, sphcnt1, sphcnt2, sphcnt3
      real(kind=8) :: radsq, radsqin
      real(kind=8) :: bdisp

!-------------------- getccl
!
!       ncpcnd(n) = ncpmx1(n) - ncpmx1(n-1)
!       (ncpcnd: number of grids for each body)
!
! 1    ncpmx1(1)       ncpmx1(2)  ... ncpmx1(n_body)
! |------|---------------|----------------|
!
!-------------------- get n{xyz}body in the global coordinate
      icap = 0
      do ipc = 1, npc
          nscpmx(ipc) = icap + 1
!
          cylaln = cylinder(ipc)%align
          cylrad = cylinder(ipc)%radius
          cylaxs1 = cylinder(ipc)%axis(1)
          cylaxs2 = cylinder(ipc)%axis(2)
          cyledg1 = cylinder(ipc)%edge(1)
          cyledg2 = cylinder(ipc)%edge(2)
!
          sphrad = sphere(ipc)%radius
          sphcnt1 = sphere(ipc)%center(1)
          sphcnt2 = sphere(ipc)%center(2)
          sphcnt3 = sphere(ipc)%center(3)
!
!GEOTYPE0
          if (geotype(ipc) .eq. 0) then
              do k = 0, nz; do j = 0, ny; do i = 0, nx
                      if (k .ge. zlpc(ipc) .and. k .le. zupc(ipc) .and. &
                     &   j .ge. ylpc(ipc) .and. j .le. yupc(ipc) .and. &
                     &   i .ge. xlpc(ipc) .and. i .le. xupc(ipc)) then
                          icap = icap + 1
                          bdygrid(1, icap) = rfold(i, nx, 1)
                          bdygrid(2, icap) = rfold(j, ny, 2)
                          bdygrid(3, icap) = rfold(k, nz, 3)
                          if (chk_dupli(icap)) icap = icap - 1
                      end if
                  end do; end do; end do

              nmxcpmx(ipc) = icap + 1

              if (nflag_subcell(1) .eq. 1) then
                  do k = 0, nz; do j = 0, ny; do i = 0, nx
                          if (k .ge. zlpc(ipc) .and. k .le. zupc(ipc) .and. &
                         &   j .ge. ylpc(ipc) .and. j .le. yupc(ipc)) then
                              if (xlpc(ipc) .ne. ceiling(xlpc(ipc))) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(xlpc(ipc), slx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
!
                              if (xupc(ipc) .ne. floor(xupc(ipc))) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(xupc(ipc), slx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end if
                      end do; end do; end do
              end if

              nmycpmx(ipc) = icap + 1

              if (nflag_subcell(1) .eq. 1) then
                  do k = 0, nz; do j = 0, ny; do i = 0, nx
                          if (k .ge. zlpc(ipc) .and. k .le. zupc(ipc) .and. &
                         &   i .ge. xlpc(ipc) .and. i .le. xupc(ipc)) then
                              if (ylpc(ipc) .ne. ceiling(ylpc(ipc))) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(ylpc(ipc), sly, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
!
                              if (yupc(ipc) .ne. floor(yupc(ipc))) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(yupc(ipc), sly, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end if
                      end do; end do; end do
              end if

              nmzcpmx(ipc) = icap + 1

              if (nflag_subcell(1) .eq. 1) then
                  do k = 0, nz; do j = 0, ny; do i = 0, nx
                          if (j .ge. ylpc(ipc) .and. j .le. yupc(ipc) .and. &
                         &   i .ge. xlpc(ipc) .and. i .le. xupc(ipc)) then
                              if (zlpc(ipc) .ne. ceiling(zlpc(ipc))) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(zlpc(ipc), slz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
!
                              if (zupc(ipc) .ne. floor(zupc(ipc))) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(zupc(ipc), slz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end if
                      end do; end do; end do
              end if

              nmxycpmx(ipc) = icap + 1
              nmyzcpmx(ipc) = icap + 1
              nmzxcpmx(ipc) = icap + 1

!GEOTYPE1
          else if (geotype(ipc) .eq. 1) then
              do k = 0, nz; do j = 0, ny; do i = 0, nx
                      if ((k .eq. ceiling(zlpc(ipc)) .or. k .eq. floor(zupc(ipc))) .and. &
                     &         j .ge. ylpc(ipc) .and. j .le. yupc(ipc) .and. &
                     &         i .ge. xlpc(ipc) .and. i .le. xupc(ipc)) then
                          icap = icap + 1
                          bdygrid(1, icap) = rfold(i, nx, 1)
                          bdygrid(2, icap) = rfold(j, ny, 2)
                          bdygrid(3, icap) = rfold(k, nz, 3)
                          if (chk_dupli(icap)) icap = icap - 1
                      else if (k .ge. zlpc(ipc) .and. k .le. zupc(ipc) .and. &
                     &        (j .eq. ceiling(ylpc(ipc)) .or. j .eq. floor(yupc(ipc))) .and. &
                     &         i .ge. xlpc(ipc) .and. i .le. xupc(ipc)) then
                          icap = icap + 1
                          bdygrid(1, icap) = rfold(i, nx, 1)
                          bdygrid(2, icap) = rfold(j, ny, 2)
                          bdygrid(3, icap) = rfold(k, nz, 3)
                          if (chk_dupli(icap)) icap = icap - 1
                      else if (k .ge. zlpc(ipc) .and. k .le. zupc(ipc) .and. &
                     &         j .ge. ylpc(ipc) .and. j .le. yupc(ipc) .and. &
                     &        (i .eq. ceiling(xlpc(ipc)) .or. i .eq. floor(xupc(ipc)))) then
                          icap = icap + 1
                          bdygrid(1, icap) = rfold(i, nx, 1)
                          bdygrid(2, icap) = rfold(j, ny, 2)
                          bdygrid(3, icap) = rfold(k, nz, 3)
                          if (chk_dupli(icap)) icap = icap - 1
                      end if
                  end do; end do; end do

              if (nflag_subcell(1) .eq. 1) then
                  if (myid .eq. 0) print *, &
                 &  "[getccl] nflag_subcell(1)==1 unavailable, => cancelled."
              end if

              nmxcpmx(ipc) = icap + 1
              nmycpmx(ipc) = icap + 1
              nmzcpmx(ipc) = icap + 1
              nmxycpmx(ipc) = icap + 1
              nmyzcpmx(ipc) = icap + 1
              nmzxcpmx(ipc) = icap + 1

!GEOTYPE2
          else if (geotype(ipc) .eq. 2) then
              radsq = cylrad**2
              radsqin = (cylrad - 1.0d0)**2
!GEOTYPE2-1
              if (cylaln .eq. 1) then
                  do k = 0, nz
                      do j = 0, ny
                          if ((j - cylaxs1)**2 + (k - cylaxs2)**2 .le. radsq) then
                              do i = ceiling(cyledg1), ceiling(cyledg1)
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end do
                          end if
!
                          if ((j - cylaxs1)**2 + (k - cylaxs2)**2 .le. radsq .and. &
                         &   (j - cylaxs1)**2 + (k - cylaxs2)**2 .ge. radsqin) then
                              do i = ceiling(cyledg1) + 1, floor(cyledg2) - 1
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end do
                          end if
!
                          if ((j - cylaxs1)**2 + (k - cylaxs2)**2 .le. radsq) then
                              do i = floor(cyledg2), floor(cyledg2)
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end do
                          end if
                      end do
                  end do

                  nmxcpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      do k = 0, nz
                          do j = 0, ny
                              if ((j - cylaxs1)**2 + (k - cylaxs2)**2 .le. radsq) then
                                  if (cyledg1 .ne. ceiling(cyledg1)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(cyledg1, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  if (cyledg2 .ne. floor(cyledg2)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(cyledg2, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end do
                      end do
                  end if

                  nmycpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      do k = 0, nz
                          if (abs(k - cylaxs2) .lt. cylrad) then
                              bdisp = sqrt(radsq - (k - cylaxs2)**2)
!
                              y = cylaxs1 - bdisp
                              if (abs(y - ceiling(y)) .gt. mingap) then
                                  do i = ceiling(cyledg1), floor(cyledg2)
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end do
                              end if
!
                              y = cylaxs1 + bdisp
                              if (abs(y - floor(y)) .gt. mingap) then
                                  do i = ceiling(cyledg1), floor(cyledg2)
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end do
                              end if
                          end if
                      end do
                  end if

                  nmzcpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      do j = 0, ny
                          if (abs(j - cylaxs1) .lt. cylrad) then
                              bdisp = sqrt(radsq - (j - cylaxs1)**2)
!
                              z = cylaxs2 - bdisp
                              if (abs(z - ceiling(z)) .gt. mingap) then
                                  do i = ceiling(cyledg1), floor(cyledg2)
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end do
                              end if
!
                              z = cylaxs2 + bdisp
                              if (abs(z - floor(z)) .gt. mingap) then
                                  do i = ceiling(cyledg1), floor(cyledg2)
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end do
                              end if
                          end if
                      end do
                  end if

                  nmxycpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      do k = 0, nz
                          if (abs(k - cylaxs2) .lt. cylrad) then
                              bdisp = sqrt(radsq - (k - cylaxs2)**2)
!
                              y = cylaxs1 - bdisp
                              if (abs(y - ceiling(y)) .gt. mingap) then
                                  if (cyledg1 .ne. ceiling(cyledg1)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(cyledg1, slx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  if (cyledg2 .ne. floor(cyledg2)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(cyledg2, slx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
!
                              y = cylaxs1 + bdisp
                              if (abs(y - floor(y)) .gt. mingap) then
                                  if (cyledg1 .ne. ceiling(cyledg1)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(cyledg1, slx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  if (cyledg2 .ne. floor(cyledg2)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(cyledg2, slx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end if
                      end do
                  end if

                  nmyzcpmx(ipc) = icap + 1
                  nmzxcpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      do j = 0, ny
                          if (abs(j - cylaxs1) .lt. cylrad) then
                              bdisp = sqrt(radsq - (j - cylaxs1)**2)
!
                              z = cylaxs2 - bdisp
                              if (abs(z - ceiling(z)) .gt. mingap) then
                                  if (cyledg1 .ne. ceiling(cyledg1)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(cyledg1, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                                  if (cyledg2 .ne. floor(cyledg2)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(cyledg2, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
!
                              z = cylaxs2 + bdisp
                              if (abs(z - floor(z)) .gt. mingap) then
                                  if (cyledg1 .ne. ceiling(cyledg1)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(cyledg1, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                                  if (cyledg2 .ne. floor(cyledg2)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(cyledg2, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end if
                      end do
                  end if

!GEOTYPE2-2
              else if (cylaln .eq. 2) then
                  do k = 0, nz
                      do j = ceiling(cyledg1), ceiling(cyledg1)
                          do i = 0, nx
                              if ((k - cylaxs1)**2 + (i - cylaxs2)**2 .le. radsq) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end do
                      end do
!
                      do j = ceiling(cyledg1) + 1, floor(cyledg2) - 1
                          do i = 0, nx
                              if ((k - cylaxs1)**2 + (i - cylaxs2)**2 .le. radsq .and. &
                             &   (k - cylaxs1)**2 + (i - cylaxs2)**2 .ge. radsqin) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end do
                      end do
!
                      do j = floor(cyledg2), floor(cyledg2)
                          do i = 0, nx
                              if ((k - cylaxs1)**2 + (i - cylaxs2)**2 .le. radsq) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end do
                      end do
                  end do

                  nmxcpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      do k = 0, nz
                          do j = ceiling(cyledg1), floor(cyledg2)
                              if (abs(k - cylaxs1) .lt. cylrad) then
                                  bdisp = sqrt(radsq - (k - cylaxs1)**2)
!
                                  x = cylaxs2 - bdisp
                                  if (abs(x - ceiling(x)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  x = cylaxs2 + bdisp
                                  if (abs(x - floor(x)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end do
                      end do
                  end if

                  nmycpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      do k = 0, nz
                          if (cyledg1 .ne. ceiling(cyledg1)) then
                              do i = 0, nx
                                  if ((k - cylaxs1)**2 + (i - cylaxs2)**2 .le. radsq) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(cyledg1, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end do
                          end if
!
                          if (cyledg2 .ne. floor(cyledg2)) then
                              do i = 0, nx
                                  if ((k - cylaxs1)**2 + (i - cylaxs2)**2 .le. radsq) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(cyledg2, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end do
                          end if
                      end do
                  end if

                  nmzcpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      do j = ceiling(cyledg1), floor(cyledg2)
                          do i = 0, nx
                              if (abs(i - cylaxs2) .lt. cylrad) then
                                  bdisp = sqrt(radsq - (i - cylaxs2)**2)
!
                                  z = cylaxs1 - bdisp
                                  if (abs(z - ceiling(z)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  z = cylaxs1 + bdisp
                                  if (abs(z - floor(z)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end do
                      end do
                  end if

                  nmxycpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      do k = 0, nz
                          if (abs(k - cylaxs1) .lt. cylrad) then
                              bdisp = sqrt(radsq - (k - cylaxs1)**2)
!
                              x = cylaxs2 - bdisp
                              if (abs(x - ceiling(x)) .gt. mingap) then
                                  if (cyledg1 .ne. ceiling(cyledg1)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(cyledg1, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                                  if (cyledg2 .ne. floor(cyledg2)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(cyledg2, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
!
                              x = cylaxs2 + bdisp
                              if (abs(x - floor(x)) .gt. mingap) then
                                  if (cyledg1 .ne. ceiling(cyledg1)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(cyledg1, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                                  if (cyledg2 .ne. floor(cyledg2)) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(cyledg2, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end if
                      end do
                  end if

                  nmyzcpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      if (cyledg1 .ne. ceiling(cyledg1)) then
                          do i = 0, nx
                              if (abs(i - cylaxs2) .lt. cylrad) then
                                  bdisp = sqrt(radsq - (i - cylaxs2)**2)
!
                                  z = cylaxs1 - bdisp
                                  if (abs(z - ceiling(z)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(cyledg1, sly, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  z = cylaxs1 + bdisp
                                  if (abs(z - floor(z)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(cyledg1, sly, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end do
                      end if
!
                      if (cyledg2 .ne. floor(cyledg2)) then
                          do i = 0, nx
                              if (abs(i - cylaxs2) .lt. cylrad) then
                                  bdisp = sqrt(radsq - (i - cylaxs2)**2)
!
                                  z = cylaxs1 - bdisp
                                  if (abs(z - ceiling(z)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(cyledg2, sly, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  z = cylaxs1 + bdisp
                                  if (abs(z - floor(z)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(cyledg2, sly, 2)
                                      bdygrid(3, icap) = rfold(z, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end do
                      end if
                  end if

                  nmzxcpmx(ipc) = icap + 1

!GEOTYPE2-3
              else if (cylaln .eq. 3) then
                  do k = ceiling(cyledg1), ceiling(cyledg1)
                      do j = 0, ny
                          do i = 0, nx
                              if ((i - cylaxs1)**2 + (j - cylaxs2)**2 .le. radsq) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end do
                      end do
                  end do
!
                  do k = ceiling(cyledg1) + 1, floor(cyledg2) - 1
                      do j = 0, ny
                          do i = 0, nx
                              if ((i - cylaxs1)**2 + (j - cylaxs2)**2 .le. radsq .and. &
                             &   (i - cylaxs1)**2 + (j - cylaxs2)**2 .ge. radsqin) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end do
                      end do
                  end do
!
                  do k = floor(cyledg2), floor(cyledg2)
                      do j = 0, ny
                          do i = 0, nx
                              if ((i - cylaxs1)**2 + (j - cylaxs2)**2 .le. radsq) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end do
                      end do
                  end do

                  nmxcpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      do k = ceiling(cyledg1), floor(cyledg2)
                          do j = 0, ny
                              if (abs(j - cylaxs2) .lt. cylrad) then
                                  bdisp = sqrt(radsq - (j - cylaxs2)**2)
!
                                  x = cylaxs1 - bdisp
                                  if (abs(x - ceiling(x)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  x = cylaxs1 + bdisp
                                  if (abs(x - floor(x)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end do
                      end do
                  end if

                  nmycpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      do k = ceiling(cyledg1), floor(cyledg2)
                          do i = 0, nx
                              if (abs(i - cylaxs1) .lt. cylrad) then
                                  bdisp = dsqrt(radsq - (i - cylaxs1)**2)
!
                                  y = cylaxs2 - bdisp
                                  if (abs(y - ceiling(y)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  y = cylaxs2 + bdisp
                                  if (abs(y - floor(y)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(k, nz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end do
                      end do
                  end if

                  nmzcpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      if (cyledg1 .ne. ceiling(cyledg1)) then
                          do j = 0, ny
                              do i = 0, nx
                                  if ((i - cylaxs1)**2 + (j - cylaxs2)**2 .le. radsq) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(cyledg1, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end do
                          end do
                      end if
!
                      if (cyledg2 .ne. floor(cyledg2)) then
                          do j = 0, ny
                              do i = 0, nx
                                  if ((i - cylaxs1)**2 + (j - cylaxs2)**2 .le. radsq) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(cyledg2, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end do
                          end do
                      end if
                  end if

                  nmxycpmx(ipc) = icap + 1
                  nmyzcpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      if (cyledg1 .ne. ceiling(cyledg1)) then
                          do i = 0, nx
                              if (abs(i - cylaxs1) .lt. cylrad) then
                                  bdisp = sqrt(radsq - (i - cylaxs1)**2)
!
                                  y = cylaxs2 - bdisp
                                  if (abs(y - ceiling(y)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(cyledg1, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  y = cylaxs2 + bdisp
                                  if (abs(y - floor(y)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(cyledg1, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end do
                      end if
!
                      if (cyledg2 .ne. floor(cyledg2)) then
                          do i = 0, nx
                              if (abs(i - cylaxs1) .lt. cylrad) then
                                  bdisp = sqrt(radsq - (i - cylaxs1)**2)
!
                                  y = cylaxs2 - bdisp
                                  if (abs(y - ceiling(y)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(cyledg2, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  y = cylaxs2 + bdisp
                                  if (abs(y - floor(y)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(i, nx, 1)
                                      bdygrid(2, icap) = rfold(y, sly, 2)
                                      bdygrid(3, icap) = rfold(cyledg2, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end do
                      end if
                  end if

                  nmzxcpmx(ipc) = icap + 1

                  if (nflag_subcell(1) .eq. 1) then
                      if (cyledg1 .ne. ceiling(cyledg1)) then
                          do j = 0, ny
                              if (abs(j - cylaxs2) .lt. cylrad) then
                                  bdisp = sqrt(radsq - (j - cylaxs2)**2)
!
                                  x = cylaxs1 - bdisp
                                  if (abs(x - ceiling(x)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(cyledg1, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  x = cylaxs1 + bdisp
                                  if (abs(x - floor(x)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(cyledg1, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end do
                      end if

                      if (cyledg2 .ne. floor(cyledg2)) then
                          do j = 0, ny
                              if (abs(j - cylaxs2) .lt. cylrad) then
                                  bdisp = sqrt(radsq - (j - cylaxs2)**2)
!
                                  x = cylaxs1 - bdisp
                                  if (abs(x - ceiling(x)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(cyledg2, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
!
                                  x = cylaxs1 + bdisp
                                  if (abs(x - floor(x)) .gt. mingap) then
                                      icap = icap + 1
                                      bdygrid(1, icap) = rfold(x, slx, 1)
                                      bdygrid(2, icap) = rfold(j, ny, 2)
                                      bdygrid(3, icap) = rfold(cyledg2, slz, 3)
                                      if (chk_dupli(icap)) icap = icap - 1
                                  end if
                              end if
                          end do
                      end if
                  end if
              end if

!GEOTYPE3
          else if (geotype(ipc) .eq. 3) then
              radsq = sphrad**2
              radsqin = (sphrad - 1.0d0)**2
              do k = 0, nz; do j = 0, ny; do i = 0, nx
                      if ((i - sphcnt1)**2 + (j - sphcnt2)**2 + (k - sphcnt3)**2 .le. radsq .and. &
                     &   (i - sphcnt1)**2 + (j - sphcnt2)**2 + (k - sphcnt3)**2 .ge. radsqin) then
                          icap = icap + 1
                          bdygrid(1, icap) = rfold(i, nx, 1)
                          bdygrid(2, icap) = rfold(j, ny, 2)
                          bdygrid(3, icap) = rfold(k, nz, 3)
                          if (chk_dupli(icap)) icap = icap - 1
                      end if
                  end do; end do; end do

              nmxcpmx(ipc) = icap + 1

              if (nflag_subcell(1) .eq. 1) then
                  do k = 0, nz
                      do j = 0, ny
                          if ((j - sphcnt2)**2 + (k - sphcnt3)**2 .lt. radsq) then
                              bdisp = sqrt(radsq - (j - sphcnt2)**2 - (k - sphcnt3)**2)
!
                              x = sphcnt1 - bdisp
                              if (abs(x - ceiling(x)) .gt. mingap) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(x, slx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
!
                              x = sphcnt1 + bdisp
                              if (abs(x - floor(x)) .gt. mingap) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(x, slx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end if
                      end do
                  end do
              end if

              nmycpmx(ipc) = icap + 1

              if (nflag_subcell(1) .eq. 1) then
                  do k = 0, nz
                      do i = 0, nx
                          if ((i - sphcnt1)**2 + (k - sphcnt3)**2 .lt. radsq) then
                              bdisp = sqrt(radsq - (i - sphcnt1)**2 - (k - sphcnt3)**2)
!
                              y = sphcnt2 - bdisp
                              if (abs(y - ceiling(y)) .gt. mingap) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(y, sly, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
!
                              y = sphcnt2 + bdisp
                              if (abs(y - floor(y)) .gt. mingap) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(y, sly, 2)
                                  bdygrid(3, icap) = rfold(k, nz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end if
                      end do
                  end do
              end if

              nmzcpmx(ipc) = icap + 1

              if (nflag_subcell(1) .eq. 1) then
                  do j = 0, ny
                      do i = 0, nx
                          if ((i - sphcnt1)**2 + (j - sphcnt2)**2 .lt. radsq) then
                              bdisp = sqrt(radsq - (i - sphcnt1)**2 - (j - sphcnt2)**2)
!
                              z = sphcnt3 - bdisp
                              if (abs(z - ceiling(z)) .gt. mingap) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(z, slz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
!
                              z = sphcnt3 + bdisp
                              if (abs(z - floor(z)) .gt. mingap) then
                                  icap = icap + 1
                                  bdygrid(1, icap) = rfold(i, nx, 1)
                                  bdygrid(2, icap) = rfold(j, ny, 2)
                                  bdygrid(3, icap) = rfold(z, slz, 3)
                                  if (chk_dupli(icap)) icap = icap - 1
                              end if
                          end if
                      end do
                  end do
              end if

              nmxycpmx(ipc) = icap + 1
              nmyzcpmx(ipc) = icap + 1
              nmzxcpmx(ipc) = icap + 1

          end if
!
          necpmx(ipc) = icap
          ncpcnt(ipc) = necpmx(ipc) - nscpmx(ipc) + 1
!
          if (myid .eq. 0) then
              print *, 'n{s,e}cpmx(', ipc, '), ncpcnt(', ipc, ') = ', &
             &       nscpmx(ipc), necpmx(ipc), ncpcnt(ipc)
              print *, 'myid=', myid, ': nbdsf1(', ipc, ')=', nbdsf1(ipc)
          end if
      end do

!--------------------
!      nbdglocal(:) = 0
!      do ipc=1,npc
!      do icap=nscpmx(ipc),necpmx(ipc)
!        i = bdygrid(1,icap)
!        j = bdygrid(2,icap)
!        k = bdygrid(3,icap)
!  slil: do islice=0,snode-1
!          if(k.ge.sinfo(2,islice).and.k.le.sinfo(3,islice)) then
!            nbdglocal(islice+1) = nbdglocal(islice+1) + 1
!            exit slil
!          end if
!        end do slil
!      end do
!      end do
!      do islice=0,snode-1
!        dispbdg(islice+1) = sum(nbdglocal(1:islice))
!      end do

!--------------------
      ibdy = 0
      do ipc = 1, npc
          if (geotype(ipc) .lt. 0) then
              do k = 0, nz
              do j = 0, ny
              do i = 0, nx
!           -----------
                  if ((k .eq. nint(zlpc(ipc)) .and. k .le. (xupc(ipc) - 1)) .and. &
                 &        (j .ge. ylpc(ipc) .and. j .le. (yupc(ipc) - 1)) .and. &
                 &        (i .ge. xlpc(ipc) .and. i .le. (xupc(ipc) - 1))) then
                      ibdy = ibdy + 1
                      bdyvoxel(1, ibdy) = i
                      bdyvoxel(2, ibdy) = j
                      bdyvoxel(3, ibdy) = k
                  else if ((k .ge. zlpc(ipc) .and. k .eq. nint((xupc(ipc) - 1))) .and. &
                 &        (j .ge. ylpc(ipc) .and. j .le. (yupc(ipc) - 1)) .and. &
                 &        (i .ge. xlpc(ipc) .and. i .le. (xupc(ipc) - 1))) then
                      ibdy = ibdy + 1
                      bdyvoxel(1, ibdy) = i
                      bdyvoxel(2, ibdy) = j
                      bdyvoxel(3, ibdy) = k
                  else if ((k .ge. zlpc(ipc) .and. k .le. (xupc(ipc) - 1)) .and. &
                 &        (j .eq. nint(ylpc(ipc)) .and. j .le. (yupc(ipc) - 1)) .and. &
                 &        (i .ge. xlpc(ipc) .and. i .le. (xupc(ipc) - 1))) then
                      ibdy = ibdy + 1
                      bdyvoxel(1, ibdy) = i
                      bdyvoxel(2, ibdy) = j
                      bdyvoxel(3, ibdy) = k
                  else if ((k .ge. zlpc(ipc) .and. k .le. (xupc(ipc) - 1)) .and. &
                 &        (j .ge. ylpc(ipc) .and. j .eq. nint((yupc(ipc) - 1))) .and. &
                 &        (i .ge. xlpc(ipc) .and. i .le. (xupc(ipc) - 1))) then
                      ibdy = ibdy + 1
                      bdyvoxel(1, ibdy) = i
                      bdyvoxel(2, ibdy) = j
                      bdyvoxel(3, ibdy) = k
                  else if ((k .ge. zlpc(ipc) .and. k .le. (xupc(ipc) - 1)) .and. &
                 &        (j .ge. ylpc(ipc) .and. j .le. (yupc(ipc) - 1)) .and. &
                 &        (i .eq. nint(xlpc(ipc)) .and. i .le. (xupc(ipc) - 1))) then
                      ibdy = ibdy + 1
                      bdyvoxel(1, ibdy) = i
                      bdyvoxel(2, ibdy) = j
                      bdyvoxel(3, ibdy) = k
                  else if ((k .ge. zlpc(ipc) .and. k .le. (xupc(ipc) - 1)) .and. &
                 &        (j .ge. ylpc(ipc) .and. j .le. (yupc(ipc) - 1)) .and. &
                 &        (i .ge. xlpc(ipc) .and. i .eq. nint((xupc(ipc) - 1)))) then
                      ibdy = ibdy + 1
                      bdyvoxel(1, ibdy) = i
                      bdyvoxel(2, ibdy) = j
                      bdyvoxel(3, ibdy) = k
                  end if
              end do
              end do
              end do
          end if
          nbdsf1(ipc) = ibdy
      end do

!-------------------- get i{xyz}body contained in each domain
      maxbedge = 0
      nbedge(:, :) = 0
      do ipc = 1, npc
          ibedge = 0
          jbedge = 0
          kbedge = 0
          JCAPL1: do jcap = nscpmx(ipc), necpmx(ipc)
              i = bdygrid(1, jcap)
              j = bdygrid(2, jcap)
              k = bdygrid(3, jcap)
              ICAPL1: do icap = nscpmx(ipc), necpmx(ipc)
                  if (i .eq. bdygrid(1, icap) - 1 .and. &
                 &   j .eq. bdygrid(2, icap) .and. &
                 &   k .eq. bdygrid(3, icap)) then
                      ibedge = ibedge + 1
                  end if
                  if (i .eq. bdygrid(1, icap) .and. &
                 &   j .eq. bdygrid(2, icap) - 1 .and. &
                 &   k .eq. bdygrid(3, icap)) then
                      jbedge = jbedge + 1
                  end if
                  if (i .eq. bdygrid(1, icap) .and. &
                 &   j .eq. bdygrid(2, icap) .and. &
                 &   k .eq. bdygrid(3, icap) - 1) then
                      kbedge = kbedge + 1
                  end if
              end do ICAPL1
          end do JCAPL1
          nbedge(1, ipc) = ibedge
          nbedge(2, ipc) = jbedge
          nbedge(3, ipc) = kbedge
          if (ibedge .gt. maxbedge) maxbedge = ibedge
          if (jbedge .gt. maxbedge) maxbedge = jbedge
          if (kbedge .gt. maxbedge) maxbedge = kbedge
      end do
!
      if (myid .eq. 0) print *, "maxbedge=", maxbedge
      allocate (bdyedges(3, maxbedge, 3, npc))
!
      do ipc = 1, npc
          ibedge = 0
          jbedge = 0
          kbedge = 0
          JCAPL2: do jcap = nscpmx(ipc), necpmx(ipc)
              i = bdygrid(1, jcap)
              j = bdygrid(2, jcap)
              k = bdygrid(3, jcap)
              ICAPL2: do icap = nscpmx(ipc), necpmx(ipc)
                  if (i .eq. bdygrid(1, icap) - 1 .and. &
                 &   j .eq. bdygrid(2, icap) .and. &
                 &   k .eq. bdygrid(3, icap)) then
                      ibedge = ibedge + 1
                      bdyedges(1, ibedge, 1, ipc) = i
                      bdyedges(2, ibedge, 1, ipc) = j
                      bdyedges(3, ibedge, 1, ipc) = k
                  end if
                  if (i .eq. bdygrid(1, icap) .and. &
                 &   j .eq. bdygrid(2, icap) - 1 .and. &
                 &   k .eq. bdygrid(3, icap)) then
                      jbedge = jbedge + 1
                      bdyedges(1, jbedge, 2, ipc) = i
                      bdyedges(2, jbedge, 2, ipc) = j
                      bdyedges(3, jbedge, 2, ipc) = k
                  end if
                  if (i .eq. bdygrid(1, icap) .and. &
                 &   j .eq. bdygrid(2, icap) .and. &
                 &   k .eq. bdygrid(3, icap) - 1) then
                      kbedge = kbedge + 1
                      bdyedges(1, kbedge, 3, ipc) = i
                      bdyedges(2, kbedge, 3, ipc) = j
                      bdyedges(3, kbedge, 3, ipc) = k
                  end if
              end do ICAPL2
          end do JCAPL2
      end do

      return

  contains

      function chk_dupli(i3cap)
          implicit none
          logical :: chk_dupli
          integer(kind=4), intent(in) :: i3cap
          integer(kind=4) :: i4cap

          chk_dupli = .false.
          do i4cap = 1, i3cap - 1
              if (bdygrid(1, i3cap) .eq. bdygrid(1, i4cap) .and. &
             &   bdygrid(2, i3cap) .eq. bdygrid(2, i4cap) .and. &
             &   bdygrid(3, i3cap) .eq. bdygrid(3, i4cap)) then
                  chk_dupli = .true.
                  exit
              end if
          end do

          return
      end function chk_dupli

  end subroutine getccl

  function rfold_i(r, nr, dir)
      use allcom
      implicit none
      integer(kind=4) :: rfold_i
      integer(kind=4), intent(in) :: r, nr
      integer(kind=4), intent(in) :: dir

      if (mtd_vbnd(dir) .ne. 0) then
          rfold_i = r
      else
          rfold_i = mod(r, nr)
      end if

      return
  end function rfold_i

  function rfold_r(r, nr, dir)
      use allcom
      implicit none
      real(kind=4) :: rfold_r
      real(kind=4), intent(in) :: r, nr
      integer(kind=4), intent(in) :: dir

      if (mtd_vbnd(dir) .ne. 0) then
          rfold_r = r
      else
          rfold_r = mod(r, nr)
      end if

      return
  end function rfold_r

  function rfold_d(r, nr, dir)
      use allcom
      implicit none
      real(kind=8) :: rfold_d
      real(kind=8), intent(in) :: r, nr
      integer(kind=4), intent(in) :: dir

      if (mtd_vbnd(dir) .ne. 0) then
          rfold_d = r
      else
          rfold_d = mod(r, nr)
      end if

      return
  end function rfold_d
