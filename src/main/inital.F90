#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine inital
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   I N I T A L
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   this subroutine gives initialization of the run,       .
!   .   including constant calculation, field setting,  and    .
!   .   and particle setting.                                  .
!   .   velocities, positions and b(magnetic field) should be  .
!   .   at t=0.5dt, and e(electric field) should be given at   .
!   .   t=dt .                                                 .
!   ............................................................
!
!-------------------- parameter and common blocks
      use oh_type
      use paramt
      use allcom
      use hdf
#define MCW local_comm
#define MSS MPI_STATUS_SIZE
      implicit none

      interface
          subroutine sumcount(in, inout, len, dtype)
              use paramt
              integer :: len
              integer :: dtype
              type(counter) :: in(len), inout(len)
          end subroutine sumcount
      end interface
      integer(kind=4) :: i, j, k, m, iran, ipc, ipcg, ig, icon, is, cnt, itch
      integer(kind=4) :: nmaxfft
      integer(kind=4) :: nranu
      integer(kind=4) :: ierr
      integer(kind=4) :: icomm, ncomm
      integer(kind=4) :: slid, layer
      integer(kind=4) :: sd
      real(kind=8) :: ce0x, ce0y, ce0z, ce0
      real(kind=8) :: totalc, rhoval
      real(kind=8) :: csz, snz, csxy, snxy, qmrss, prqmrs, tmppc
      real(kind=8) :: dsckx, dscky, dsckz
      real(kind=8) :: dran
!
      integer(kind=HID_T) :: fileid
      integer(kind=4) :: stats0, stats1
      integer(kind=8) :: dims(1)
      character(len=30) :: filename, dsname

!-------------------- basic constant values
!     --------------- value of pi
      pi = acos(-1.0d0)
      pi2 = pi*2.0d0
      pih = pi*0.5d0
!     --------------- inverse of grid span dr
      dri = 1.0d0/dr
      si = dri**2.0d0
!     --------------- global system length
      slx = nx*dr
      sly = ny*dr
      slz = nz*dr
!     --------------- number of grids for FFT
      nxfft = nx
      nyfft = ny
      nzfft = nz
!     --------------- spatial resolution in fourier spectrum
      dkx = pi2/(nx*dr)
      dky = pi2/(ny*dr)
      dkz = pi2/(nz*dr)
!     --------------- squared/doubled light speed
      cs = cv**2
      tcs = cs*2.0d0
!     --------------- for implicit solver
      cfactor(1:3) = gfactor*cv*dt*dri

!-------------------- initialiazation of time step count
      istep = 0

!-------------------- static magnetic field setting
!
!
!           z                inital direction of B : +z
!             |
!         b0z |.            ex.) phiz=0, phixy=0 --> B = (0, 0, B0)
!             | .
!             |  .
!             |   . |B|=B0
!             |   /.
!            phiz/ .
!             | /  .
!             |/   .   B0y
!             +----.---.------- y
!            / \   .  .
!          phixy\  . .
!          /     \ ..
!     B0x /.......\.
!        /
!       x
!
!
!  (0,0,B0) ----> T(phiz,phixy) ---> (B0x,B0y,B0z)
!
!                        +-           -+
!                        | t11 t12 t13 |
!   where T(phiz,phixy)= | t21 t22 t23 |
!                        | t31 t32 t33 |
!                        +-           -+
!
!     ---------------
      b0 = wc/qm(1)
!
      csz = dcos(pi*phiz/180.0d0)
      snz = dsin(pi*phiz/180.0d0)
      csxy = dcos(pi*phixy/180.0d0)
      snxy = dsin(pi*phixy/180.0d0)
!
      t11 = csz*csxy
      t12 = -snxy
      t13 = snz*csxy
      t21 = csz*snxy
      t22 = csxy
      t23 = snz*snxy
      t31 = -snz
      t32 = 0.0d0
      t33 = csz
!
      b0x = b0*t13
      b0y = b0*t23
      b0z = b0*t33
!
      if (b0 .ne. 0.0d0) then
          if (dabs(b0x/b0) .lt. 1.0d-10) b0x = 0.0d0
          if (dabs(b0y/b0) .lt. 1.0d-10) b0y = 0.0d0
          if (dabs(b0z/b0) .lt. 1.0d-10) b0z = 0.0d0
      end if
!
      if (myid .eq. 0) write (*, *) 'B0 =', b0
      if (myid .eq. 0) write (*, *) 'B0x =', b0x
      if (myid .eq. 0) write (*, *) 'B0y =', b0y
      if (myid .eq. 0) write (*, *) 'B0z =', b0z
!
! * transfer matrix for vpara direction
!   If we initialize (Vx,Vy,Vz) to (Vperp1, Vperp2, Vpara),
!   we should transfer it by T(phiz,phixy).
!   And if we want to know value of B-para element,
!   x-,y-,z-components should be transferd by T**(-1).
!
!                   +-           -+
!        -1   t     | t11 t21 t31 |
!       T  ==  T =  | t12 t22 t32 |
!                   | t13 t23 t33 | .
!                   +-           -+

!-------------------- static electric field setting
      do is = 1, nspec
          vdtx(is) = vdx(is) + spe(is)*cos(speth(is)/180.0d0*pi)*t11 &
           &       + spe(is)*sin(speth(is)/180.0d0*pi)*t12 + spa(is)*t13 &
           &       + vdri(is)*sin(vdthz(is)/180.0d0*pi)*cos(vdthxy(is)/180.0d0*pi)
          vdty(is) = vdy(is) + spe(is)*cos(speth(is)/180.0d0*pi)*t21 &
           &       + spe(is)*sin(speth(is)/180.0d0*pi)*t22 + spa(is)*t23 &
           &       + vdri(is)*sin(vdthz(is)/180.0d0*pi)*sin(vdthxy(is)/180.0d0*pi)
          vdtz(is) = vdz(is) + spe(is)*cos(speth(is)/180.0d0*pi)*t31 &
           &       + spe(is)*sin(speth(is)/180.0d0*pi)*t32 + spa(is)*t33 &
           &       + vdri(is)*cos(vdthz(is)/180.0d0*pi)
      end do

      if (myid .eq. 0) write (*, *) 'Vd =', &
     &  dsqrt(vdtx(1)*vdtx(1) + vdty(1)*vdty(1) + vdtz(1)*vdtz(1))
      if (myid .eq. 0) write (*, *) 'Vdx =', vdtx(1)
      if (myid .eq. 0) write (*, *) 'Vdy =', vdty(1)
      if (myid .eq. 0) write (*, *) 'Vdz =', vdtz(1)

      ce0x = -(vdty(1)*b0z - vdtz(1)*b0y)
      ce0y = -(vdtz(1)*b0x - vdtx(1)*b0z)
      ce0z = -(vdtx(1)*b0y - vdty(1)*b0x)
      ce0 = dsqrt(ce0x*ce0x + ce0y*ce0y + ce0z*ce0z)

      e0x = e0x + ce0x
      e0y = e0y + ce0y
      e0z = e0z + ce0z
      e0 = dsqrt(e0x*e0x + e0y*e0y + e0z*e0z)

      if (myid .eq. 0 .and. e0 .gt. 0.0d0) write (*, *) 'E0conv =', ce0
      if (myid .eq. 0 .and. e0 .gt. 0.0d0) write (*, *) 'E0xconv =', ce0x
      if (myid .eq. 0 .and. e0 .gt. 0.0d0) write (*, *) 'E0yconv =', ce0y
      if (myid .eq. 0 .and. e0 .gt. 0.0d0) write (*, *) 'E0zconv =', ce0z
      if (myid .eq. 0 .and. e0 .gt. 0.0d0) write (*, *) 'E0 =', e0
      if (myid .eq. 0 .and. e0 .gt. 0.0d0) write (*, *) 'E0x =', e0x
      if (myid .eq. 0 .and. e0 .gt. 0.0d0) write (*, *) 'E0y =', e0y
      if (myid .eq. 0 .and. e0 .gt. 0.0d0) write (*, *) 'E0z =', e0z

      ! Cyclotron frequency & Larmor radius display
      do is = 1, nspec
          if (myid .eq. 0) print *, "Cyclotron freq. for is =", is, &
         &                      "is ", abs(wc*qm(is)/qm(1))
          if (myid .eq. 0 .and. wc .ne. 0.0d0) then
              print *, "Larmor radius   for is =", is, &
             &        "is ", abs(peth(is)/wc/qm(is)*qm(1))
          else if (myid .eq. 0 .and. wc .eq. 0.0d0) then
              print *, "Larmor radius   for is =", is, &
             &        "is INFINITE"
          end if
      end do

      ! plasma frequency redefinition
      totalc = 0.0d0
      do is = 1, nspec
          if (wp(is) .lt. 0.0d0 .and. npin(is) .gt. 0) totalc = totalc + abs(wp(is))
      end do
      rhoval = wp(1)**2/qm(1)
      do is = 1, nspec
          if (wp(is) .lt. 0.0d0) then
              wp(is) = sqrt(abs(rhoval*wp(is)/totalc*qm(is)))
              if (myid .eq. 0) print *, "wp(", is, ") is redefined as", wp(is)
          end if
      end do

      ! macro-particle definition
      qmrss = 1.0d0
      npsum = 0
      do is = 1, ispec
          if (is .le. nspec) then
              ! particle charge for each species
              if (npin(is) .gt. 0) then
                  q(is) = (wp(is)**2)*slx*sly*slz/npin(is)/qm(is)
                  if (myid .eq. 0) write (6, *) 'species', is, &
                      'q(is) was determined by npin(is)', npin(is)
              else if (q(is) .ne. 0.0d0) then
                  if (myid .eq. 0) write (6, *) 'species', is, &
                      'q(is) was determined by qp(is)', q(is)
              else if (is .gt. 1 .and. qpr(is) .ne. 0.0d0) then
                  q(is) = q(1)*qpr(is)
                  if (myid .eq. 0) write (6, *) 'species', is, &
                      'q(is) was determined by qpr(is)', qpr(is)
              else if (dnsf(is) .gt. 0.0d0) then
                  q(is) = (wp(is)**2)*dr*dr*dr/dnsf(is)/qm(is)
                  if (myid .eq. 0) write (6, *) 'species', is, &
                      'q(is) was determined by dnsf(is)', dnsf(is)
              else if (dnsb(is) .gt. 0.0d0) then
                  q(is) = (wp(is)**2)*dr*dr*dr/dnsb(is)/qm(is)
                  if (myid .eq. 0) write (6, *) 'species', is, &
                      'q(is) was determined by dnsb(is)', dnsb(is)
              else
                  q(is) = 0.0d0
                  if (myid .eq. 0) write (6, *) 'species', is, 'q(is) = 0.0d0'
              end if
              if (myid .eq. 0) print *, 'q(', is, ')', q(is)

              ! charge density ref.
              if (rho0 .le. 0.0d0 .and. is .eq. int(abs(rho0))) then
                  if (wp(is) .ne. 0.0d0) then
                      rho0 = abs((wp(is)**2)/qm(is))
                      if (myid .eq. 0) print *, 'rho0 was determined by wp(', is, ')', rho0
                  else if (qp(is) .ne. 0.0d0) then
                      rho0 = abs(q(is))
                      if (myid .eq. 0) print *, 'rho0 was determined by q(', is, ')', rho0
                  else
                      rho0 = 1.0d0
                      if (myid .eq. 0) print *, 'rho0 = 1.0d0'
                  end if
              end if

              ! particle mass for each species
              ! qmrs is relative q-m ratio: i.e., qmrs(is)=qmr(is)/qmr(is-1) [qmr(0)=1.0]
              rm(is) = q(is)/qm(is)
              qmr(is) = qm(is)/qm(1)
              prqmrs = qmrss
              qmrss = qmr(is)
              qmrs(is) = qmr(is)/prqmrs
              ! number of total particles
              npsum = npsum + np(is)
          else
              q(is) = 0.0d0
              rm(is) = 0.0d0
              qmr(is) = 0.0d0
              qmrs(is) = 0.0d0
          end if
      end do

      ! absorbing region setting
      if (nfbnd(1) .eq. 2) then
          if (nxl .lt. 0) nxl = nx/4
          if (nxr .lt. 0) nxr = nx/4
          nxl0 = nxl
          nxr0 = nx - nxr
      else
          nxl = 0
          nxr = 0
          nxl0 = -2
          nxr0 = nx + 1
      end if
      if (nfbnd(2) .eq. 2) then
          if (nyl .lt. 0) nyl = ny/4
          if (nyr .lt. 0) nyr = ny/4
          nyl0 = nyl
          nyr0 = ny - nyr
      else
          nyl = 0
          nyr = 0
          nyl0 = -2
          nyr0 = ny + 1
      end if
      if (nfbnd(3) .eq. 2) then
          if (nzl .lt. 0) nyl = nz/4
          if (nzr .lt. 0) nyr = nz/4
          nzl0 = nzl
          nzr0 = nz - nzr
      else
          nzl = 0
          nzr = 0
          nzl0 = -2
          nzr0 = nz + 1
      end if

      do sd = 0, nnode - 1
          if (nxl0 .le. sdoms(2, 1, sd + 1)) then
              medges(1, 1, sd + 1) = nxl0
          else
              medges(1, 1, sd + 1) = sdoms(2, 1, sd + 1)
          end if
          if (nxr0 .ge. sdoms(1, 1, sd + 1)) then
              medges(2, 1, sd + 1) = nxr0
          else
              medges(2, 1, sd + 1) = sdoms(1, 1, sd + 1)
          end if

          if (nyl0 .le. sdoms(2, 2, sd + 1)) then
              medges(1, 2, sd + 1) = nyl0
          else
              medges(1, 2, sd + 1) = sdoms(2, 2, sd + 1)
          end if
          if (nxr0 .ge. sdoms(1, 2, sd + 1)) then
              medges(2, 2, sd + 1) = nyr0
          else
              medges(2, 2, sd + 1) = sdoms(1, 2, sd + 1)
          end if

          if (nzl0 .le. sdoms(2, 3, sd + 1)) then
              medges(1, 3, sd + 1) = nzl0
          else
              medges(1, 3, sd + 1) = sdoms(2, 3, sd + 1)
          end if
          if (nzr0 .ge. sdoms(1, 3, sd + 1)) then
              medges(2, 3, sd + 1) = nzr0
          else
              medges(2, 3, sd + 1) = sdoms(1, 3, sd + 1)
          end if
      end do

      ! initialization for potential calculation
      if ((pftmode .eq. 4) .or. &
          (pftmode .eq. 0 .and. nodes(1)*nodes(2) .eq. 1)) then
          snode = nnode
          pftmode = 4
      else if ((pftmode .eq. 3) .or. &
               (pftmode .eq. 0 .and. nodes(2)*nodes(3) .eq. 1)) then
          snode = nnode
          pftmode = 3
      else if (pftmode .eq. 2 .or. (pftmode .eq. 0 .and. ny .gt. nz)) then
          snode = getgcd(nx, ny)
          snode = lediv(snode, nnode)
          pftmode = 2
      else
          snode = getgcd(nx, nz)
          snode = lediv(snode, nnode)
          pftmode = 1
      end if
      if (myid .eq. 0) print *, "pftmode,snode =", pftmode, snode, nx, ny, nz

      nxslc = nx/snode
      nyslc = ny/snode
      nzslc = nz/snode

      if (myid .lt. snode) then
          if (pftmode .eq. 1) then
              sxl = 0; sxu = sxl + nx
              stxl = nxslc*myid; stxu = stxl + nxslc
              syl = 0; syu = syl + ny
              szl = nzslc*myid; szu = szl + nzslc
          else if (pftmode .eq. 2) then
              sxl = 0; sxu = sxl + nx
              stxl = nxslc*myid; stxu = stxl + nxslc
              syl = nyslc*myid; syu = syl + nyslc
              szl = 0; szu = szl + nz
          else if (pftmode .eq. 3) then
              sxl = nxslc*myid; sxu = sxl + nxslc
              stxl = 0; stxu = stxl + nx
              syl = 0; syu = syl + ny
              szl = 0; szu = szl + nz
          else if (pftmode .eq. 4) then
              sxl = 0; sxu = sxl + nx
              stxl = 0; stxu = stxl + nx
              syl = 0; syu = syl + ny
              szl = nzslc*myid; szu = szl + nzslc
          else
              sxl = nxslc*myid; sxu = sxl + nxslc
              stxl = nxslc*myid; stxu = stxl + nxslc
              syl = nyslc*myid; syu = syl + nyslc
              szl = nzslc*myid; szu = szl + nzslc
          end if
      else
          sxl = nx; sxu = sxl
          stxl = nx; stxu = stxl
          syl = ny; syu = syl
          szl = nz; szu = szl
      end if

      if (pftmode .eq. 1) then
          do sd = 1, snode
              slcs(sd, 1) = nzslc*(sd - 1)
              slcs(sd, 2) = nzslc*sd
              slcs(sd, 3) = slcs(sd, 2) - slcs(sd, 1)
          end do
      else if (pftmode .eq. 2) then
          do sd = 1, snode
              slcs(sd, 1) = nyslc*(sd - 1)
              slcs(sd, 2) = nyslc*sd
              slcs(sd, 3) = slcs(sd, 2) - slcs(sd, 1)
          end do
      else if (pftmode .eq. 3) then
          do sd = 1, snode
              slcs(sd, 1) = nxslc*(sd - 1)
              slcs(sd, 2) = nxslc*sd
              slcs(sd, 3) = slcs(sd, 2) - slcs(sd, 1)
          end do
      else if (pftmode .eq. 4) then
          do sd = 1, snode
              slcs(sd, 1) = nzslc*(sd - 1)
              slcs(sd, 2) = nzslc*sd
              slcs(sd, 3) = slcs(sd, 2) - slcs(sd, 1)
          end do
      else
          if (myid .eq. 0) print *, "Error@inital: something wrong..."
      end if

      do sd = snode + 1, nnode
          slcs(sd, 1) = slcs(snode, 2)
          slcs(sd, 2) = slcs(snode, 2)
          slcs(sd, 3) = slcs(sd, 2) - slcs(sd, 1)
      end do

      if (pftmode .eq. 1) then
          nwslc = max(nx, nz) + snode
          ndslc = ny
          nhslc = max(nxslc, nzslc) + 1
      else if (pftmode .eq. 2) then
          nwslc = max(nx, ny) + snode
          ndslc = max(nxslc, nyslc) + 1
          nhslc = nz
      else if (pftmode .eq. 3 .or. pftmode .eq. 4) then
          nwslc = nx
          ndslc = ny
          nhslc = nz
      else
          if (myid .eq. 0) print *, "Error@inital: something wrong..."
      end if
      if (myid .eq. 0) print *, "n{w,d,h}slc =", nx, ny, nz, snode, nwslc, ndslc, nhslc

      allocate (dbs(3, -1:nwslc + 1, -1:ndslc + 1, -1:nhslc + 1, 2))
      allocate (poi(1, -1:nwslc + 1, -1:ndslc + 1, -1:nhslc + 1, 2))
      lwslc = size(dbs, 2)
      ldslc = size(dbs, 3)
      lhslc = size(dbs, 4)
      if (myid .eq. 0) print *, "l{w,d,h}slc =", lwslc, ldslc, lhslc
      if ((pftmode .eq. 1 .and. nfbnd(3) .eq. 0) .or. &
          (pftmode .eq. 2 .and. nfbnd(2) .eq. 0) .or. &
          (pftmode .eq. 3 .and. nfbnd(1) .eq. 0) .or. &
          (pftmode .eq. 4 .and. nfbnd(3) .eq. 0)) then
          lslice(1) = mod(myid - 1 + snode, snode)
          uslice(1) = mod(myid + 1 + snode, snode)
      else
          lslice(1) = myid - 1
          uslice(1) = myid + 1
          if (myid .eq. 0) lslice(1) = MPI_PROC_NULL
          if (myid .eq. snode - 1) uslice(1) = MPI_PROC_NULL
      end if
      if ((pftmode .eq. 1 .and. mtd_vbnd(3) .eq. 0) .or. &
          (pftmode .eq. 2 .and. mtd_vbnd(2) .eq. 0) .or. &
          (pftmode .eq. 3 .and. mtd_vbnd(1) .eq. 0) .or. &
          (pftmode .eq. 4 .and. mtd_vbnd(3) .eq. 0)) then
          lslice(2) = mod(myid - 1 + snode, snode)
          uslice(2) = mod(myid + 1 + snode, snode)
      else
          lslice(2) = myid - 1
          uslice(2) = myid + 1
          if (myid .eq. 0) lslice(2) = MPI_PROC_NULL
          if (myid .eq. snode - 1) uslice(2) = MPI_PROC_NULL
      end if

      if (pftmode .eq. 1) then
          allocate (xtype(0:nodes(1) - 1), ytype(0:nodes(2) - 1), ztype(1))
          xtype(0:nodes(1) - 2) = 0; xtype(nodes(1) - 1) = 1
          ytype(0:nodes(2) - 2) = 0; ytype(nodes(2) - 1) = 1
          ztype(1) = 0
          allocate (istatus(MSS, NCOMREQ), ireqs(NCOMREQ))
      else if (pftmode .eq. 2) then
          allocate (xtype(0:nodes(1) - 1), ytype(1), ztype(0:nodes(3) - 1))
          xtype(0:nodes(1) - 2) = 0; xtype(nodes(1) - 1) = 1
          ytype(1) = 0
          ztype(0:nodes(3) - 2) = 0; ztype(nodes(3) - 1) = 1
          allocate (istatus(MSS, NCOMREQ), ireqs(NCOMREQ))
      else
          allocate (xtype(1), ytype(1), ztype(1))
          xtype(1) = 0
          ytype(1) = 0
          ztype(1) = 0
          allocate (istatus(MSS, NCOMREQ), ireqs(NCOMREQ))
      end if

      if (pftmode .eq. 1) then
          allocate (warray1d(-1:nz + 1, 2))
          allocate (warray2d(-1:nx + 1, -1:ny + 1, 2))
      else if (pftmode .eq. 2) then
          allocate (warray1d(-1:ny + 1, 2))
          allocate (warray2d(-1:nx + 1, -1:nz + 1, 2))
      else if (pftmode .eq. 3) then
          allocate (warray1d(-1:nx + 1, 2))
          allocate (warray2d(-1:ny + 1, -1:nz + 1, 2))
      else if (pftmode .eq. 4) then
          allocate (warray1d(-1:nz + 1, 2))
          allocate (warray2d(-1:nx + 1, -1:ny + 1, 2))
      else
          if (myid .eq. 0) print *, "Error@inital: something wrong..."
      end if
!
      allocate (kmod(3, 0:max(nx, ny, nz), 2))
      warray1d(:, :) = 0.0d0
      warray2d(:, :, :) = 0.0d0
      kmod(:, :, :) = 0.0d0
      rnx = 1.0d0/nx
      rny = 1.0d0/ny
      rnz = 1.0d0/nz

      call defmpi
      call inifft

      ! initialization for conducting body charging
      prefcrd = (/nx/2, ny/2, nz/2/)

      do i = 1, ixw
          sfrho(i) = 0.0d0
          sfrhoh(i) = 0.0d0
      end do

      if (jobnum(1) .eq. 0) then
          gcount(1)%chgacm(1:2, 1:npc) = 0.0d0
      else
          if (jobnum(1) .eq. 1) then
              write (filename, '(a,i4.4,a)') './SNAPSHOT0/esdat', myid, '.h5'
          else
              write (filename, '(a,i4.4,a)') './SNAPSHOT0/emdat', myid, '.h5'
          end if
          call hdfopen(filename, fileid, DFACC_READ)

          dsname = 'chgacm'
          dims(1) = npc
          if (npc .gt. 0) then
              call read1d(fileid, dsname, dims(1:1), &
             &  gcount(1)%chgacm(1, 1:npc), stats0, stats1)
          else
              gcount(1)%chgacm(1, 1:npc) = 0.0d0
          end if
          gcount(1)%chgacm(2, 1:npc) = 0.0d0

          call hdfclose(fileid, stats0)
      end if

      allocate (ndcpmx(nnode), ndcpmx2(nnode))
      allocate (nsdcpmx(nnode), nsdcpmx2(nnode), nedcpmx(nnode))
      allocate (nbdglocal(nnode), dispbdg(nnode))

      ! initialization for conducting body charging
      if (npc .gt. 0) then
          if (npcg .le. 0 .or. npcg .gt. npc) then
              if (myid .eq. 0) print *, "npcg is forcely changed... ", npcg, "->", npc
              npcg = npc
              do ipc = 1, npc + 1
                  pcgs(ipc) = ipc - 1
                  ccgs(ipc) = ipc - 1
              end do
          else
              if (sum(pcgs(1:npcg)) .ne. npc .or. sum(ccgs(1:npcg)) .ne. npc) then
                  print *, "sum of pcgs & ccgs should be equal to npc => stop"
                  stop
              end if

              cnt = 0
              do ipcg = 1, npcg
                  pcgs(npcg + 1) = cnt
                  cnt = cnt + pcgs(ipcg)
                  pcgs(ipcg) = pcgs(npcg + 1)
              end do
              pcgs(npcg + 1) = cnt

              cnt = 0
              do ipcg = 1, npcg
                  ccgs(npcg + 1) = cnt
                  cnt = cnt + ccgs(ipcg)
                  ccgs(ipcg) = ccgs(npcg + 1)
              end do
              ccgs(npcg + 1) = cnt
          end if
      end if

      do ipcg = 1, npcg
          if (mtd_vchg(ipcg) .eq. 1) then
              i = 1
              do while (biasc(i)%to .ne. 0)
                  if (abs(biasc(i)%to) .ge. ccgs(ipcg) + 1 .and. &
                 &   abs(biasc(i)%to) .le. ccgs(ipcg + 1)) then
                      if (myid .eq. 0) &
                     &  print *, "Bias current magnitude for", i, &
                     &         "is treated as a variable."
                      biasc(i)%to = -abs(biasc(i)%to)
                      biasc(i)%val = 0.0d0
                  end if
                  i = i + 1
              end do
          end if
          groupph(ipcg) = 0.0d0
      end do

      ! conducting body dimensions
      do ipc = 1, npc
          if (nxpc1(ipc) .gt. nxpc2(ipc)) then
              tmppc = nxpc1(ipc)
              nxpc1(ipc) = nxpc2(ipc)
              nxpc2(ipc) = tmppc
          end if
          if (nypc1(ipc) .gt. nypc2(ipc)) then
              tmppc = nypc1(ipc)
              nypc1(ipc) = nypc2(ipc)
              nypc2(ipc) = tmppc
          end if
          if (nzpc1(ipc) .gt. nzpc2(ipc)) then
              tmppc = nzpc1(ipc)
              nzpc1(ipc) = nzpc2(ipc)
              nzpc2(ipc) = tmppc
          end if
          if (nxpc1(ipc) .lt. 0) nxpc1(ipc) = 0
          if (nypc1(ipc) .lt. 0) nypc1(ipc) = 0
          if (nzpc1(ipc) .lt. 0) nzpc1(ipc) = 0
          if (nxpc2(ipc) .gt. nx) nxpc2(ipc) = nx
          if (nypc2(ipc) .gt. ny) nypc2(ipc) = ny
          if (nzpc2(ipc) .gt. nz) nzpc2(ipc) = nz
!
          if (xlpc(ipc) .eq. -1.0d0) then
              xlpc(ipc) = nxpc1(ipc) - oradius(1, ipc)
              if (myid .eq. 0) print *, "ipc,xlpc =", ipc, xlpc(ipc)
          end if
          if (ylpc(ipc) .eq. -1.0d0) then
              ylpc(ipc) = nypc1(ipc) - oradius(2, ipc)
              if (myid .eq. 0) print *, "ipc,ylpc =", ipc, ylpc(ipc)
          end if
          if (zlpc(ipc) .eq. -1.0d0) then
              zlpc(ipc) = nzpc1(ipc) - oradius(3, ipc)
              if (myid .eq. 0) print *, "ipc,zlpc =", ipc, zlpc(ipc)
          end if
          if (xupc(ipc) .eq. -1.0d0) then
              xupc(ipc) = nxpc2(ipc) + oradius(1, ipc)
              if (myid .eq. 0) print *, "ipc,xupc =", ipc, xupc(ipc)
          end if
          if (yupc(ipc) .eq. -1.0d0) then
              yupc(ipc) = nypc2(ipc) + oradius(2, ipc)
              if (myid .eq. 0) print *, "ipc,yupc =", ipc, yupc(ipc)
          end if
          if (zupc(ipc) .eq. -1.0d0) then
              zupc(ipc) = nzpc2(ipc) + oradius(3, ipc)
              if (myid .eq. 0) print *, "ipc,zupc =", ipc, zupc(ipc)
          end if

          sqdscaled(ipc) = dscaled(ipc)*dscaled(ipc)
      end do

      ! initialization for testch
      do itch = 1, ntch
          if (qtch(itch) .ne. 0.0d0) then
              e1tch(itch) = qtch(itch)/4.0d0/pi
              if (myid .eq. 0) print *, "e1tch", itch, "=", e1tch(itch), "based on qtch"
          else if (p1tch(itch) .ne. 0.0d0) then
              e1tch(itch) = p1tch(itch)
              if (myid .eq. 0) print *, "e1tch", itch, "=", e1tch(itch), "based on p1tch"
          else
              if (myid .eq. 0) print *, "e1tch", itch, "=", e1tch(itch), "based on e1tch"
          end if
          r2cutoff(itch) = rcutoff(itch)*rcutoff(itch)
      end do

!-------------------- initialization for antenna feed treatment
      do ig = 1, ngap
          offez(ig) = 0.0d0
          annrms(ig) = 0.0d0
      end do

!-------------------- initialization of time counters
      t = 0.0d0
      elatime = 0.0d0

      ! random number generator test
      if (myid .eq. 0) print *, 'random number test...'
      iseed = (myid + 1)*10
      call randinit0(iseed)
      nranu = 1
      call RANU0(dran, nranu, icon)
      if (icon .ne. 0) then
          print *, "Warning(RANU0): myid,icon=", myid, icon
          stop
      end if
      call MPI_Gather(dran, nranu, MPI_REAL8, dranu, nranu, MPI_REAL8, 0, MCW, ierr)
      if (myid .eq. 0) then
          m = 0
          do iran = 1, 100
              m = m + max(1, nranu*nnode/100)
              if (m .gt. nranu*nnode) exit
              print *, "  RAU: dranu(',iran,') =", dranu(m)
          end do
      end if
      call RANN0(0.0d0, 1.0d0, dran, nranu, icon)
      if (icon .ne. 0) then
          print *, "Warning(RANN0): myid,icon=", myid, icon
          stop
      end if
      call MPI_Gather(dran, nranu, MPI_REAL8, dranu, nranu, MPI_REAL8, 0, MCW, ierr)
      if (myid .eq. 0) then
          do iran = 1, 100
              m = m + max(1, nranu*nnode/100)
              if (m .gt. nranu*nnode) exit
              print *, "  RAN: dranu(',iran,') =", dranu(m)
          end do
      end if

      lfdiag = -1
      ljdiag = -1
      ladiag = -1
      lpdiag = -1

      ! labeling sub-domains located at global-domain edges
      do sd = 0, nnode - 1
          if (sdoms(1, 1, sd + 1) .eq. 0) then
              bared(1, 1, sd + 1) = 1
          else
              bared(1, 1, sd + 1) = 0
          end if
          if (sdoms(2, 1, sd + 1) .eq. nx) then
              bared(2, 1, sd + 1) = 1
          else
              bared(2, 1, sd + 1) = 0
          end if
          if (sdoms(1, 2, sd + 1) .eq. 0) then
              bared(1, 2, sd + 1) = 1
          else
              bared(1, 2, sd + 1) = 0
          end if
          if (sdoms(2, 2, sd + 1) .eq. ny) then
              bared(2, 2, sd + 1) = 1
          else
              bared(2, 2, sd + 1) = 0
          end if
          if (sdoms(1, 3, sd + 1) .eq. 0) then
              bared(1, 3, sd + 1) = 1
          else
              bared(1, 3, sd + 1) = 0
          end if
          if (sdoms(2, 3, sd + 1) .eq. nz) then
              bared(2, 3, sd + 1) = 1
          else
              bared(2, 3, sd + 1) = 0
          end if
          if (sd .eq. 0) then
              bared(1, 4, sd + 1) = 1
          else
              bared(1, 4, sd + 1) = 0
          end if
          if (sd .eq. snode - 1) then
              bared(2, 4, sd + 1) = 1
          else
              bared(2, 4, sd + 1) = 0
          end if
      end do

      ! memory allocation for non-blocking comm.
      ncomm = max(ixw*4, &
                  nodes(1)*nodes(2)*4 + nnode/nodes(3) + 4, &
                  nodes(1)*nodes(2)*4 + nz/nodes(3) + 3)
      if (myid .eq. 0) print *, "Max times of non-blocking comm.:", ncomm

      ! creating output Files
      if (intfoc .ne. 0 .and. myid .eq. 0) then
          open (50, file='oltime', status='replace')
          close (50)
          open (79, file='energy', status='replace')
          close (79)
          open (80, file='currnt', status='replace')
          close (80)
          open (81, file='volt', status='replace')
          close (81)
          if (mode_dipole .eq. 2) then
              open (82, file='jsource', status='replace')
              close (82)
          end if
          if (mode_dipole .eq. 3) then
              open (83, file='curnte', status='replace')
              close (83)
          end if
          open (84, file='pbodyr', status='replace')
          close (84)
          open (85, file='pbodyd', status='replace')
          close (85)
          open (86, file='pbody', status='replace')
          close (86)
          open (89, file='influx', status='replace')
          close (89)
          open (90, file='noflux', status='replace')
          close (90)
          open (91, file='isflux', status='replace')
          close (91)
          open (92, file='nesc', status='replace')
          close (92)
          open (93, file='chgacm1', status='replace')
          close (93)
          open (94, file='chgacm2', status='replace')
          close (94)
          open (95, file='chgmov', status='replace')
          close (95)
          open (96, file='seyield', status='replace')
          close (96)
          open (97, file='icur', status='replace')
          close (97)
          open (98, file='ocur', status='replace')
          close (98)
          open (99, file='ewave', status='replace')
          close (99)
          open (3100, file='energy1', status='replace')
          close (3100)
          open (3101, file='energy2', status='replace')
          close (3101)
      end if

      return

  contains

      integer function getgcd(x, y)
          implicit none
          integer, intent(in) :: x, y
          integer :: a, b, r, tmp

          a = x
          b = y

          if (a < b) then
              tmp = a
              a = b
              b = tmp
          end if

          r = mod(a, b)
          do while (r .ne. 0)
              a = b
              b = r
              r = mod(a, b)
          end do
          getgcd = b

          return
      end function getgcd

      !> Find a specific value.
      !    where
      !      value <= y
      !      mod(x, value) == 0
      !
      !   example:
      !     (x, y) =  1,  2  => 1
      !     (x, y) =  4,  2  => 2
      !     (x, y) = 12, 10  => 6
      integer function lediv(x, y)
          implicit none
          integer, intent(in) :: x, y
          integer :: a, b, i

          a = x
          b = y

          do i = a, 1, -1
              if (mod(a, i) .eq. 0 .and. i .le. b) then
                  lediv = i
                  return
              end if
          end do
          print *, "[lediv]Warning", a, b, i
          lediv = i

          return
      end function lediv

  end subroutine inital
