! - gap voltage -----
! ------
!      81: voltage at the gap (should be minus of the efield)
!           ez is the unit of amjz

      if(ps.eq.1.and. &
     &   ip(1).gt.0.and.ip(1).le.xu-xl.and. &
     &   jp(1).gt.0.and.jp(1).le.yu-yl.and. &
     &   kp(1).gt.0.and.kp(1).le.zu-zl) then
        if(intfoc.ne.0) then
        open(80,file='currnt',position='append')
        open(81,file='volt',position='append')
        if(mode_dipole.eq.2) open(82,file='jsource',position='append')
        if(mode_dipole.eq.3) open(83,file='curnte',position='append')
        write(81,*) t, (-eb(EZ,ip(ig),jp(ig),kp(ig),ps)/rene, &
       &                ig=1,ngap)
        if(mode_dipole.eq.2) write(82,*) t, (gap_jsrc(ig)/renj, ig=1,ngap)
        if(mode_dipole.eq.3) then
          do ig=1,ngap
            currnt(ig) = (eb(EZ,ip(ig),jp(ig),kp(ig),ps) + ezgap(ig)) &
         &                *0.5d0/resc(ig)
          end do
          write(83,*) t, (currnt(ig)/renb*(cv*cv/renv/renv), ig=1,ngap)
        end if
        end if


! - current  -----
!--------
!      80: currnt from ampere's law
!
        do ig=1,ngap
          currnt(ig) = eb(BX,ip(ig)  ,jp(ig)-1,kp(ig),ps) &
         &           + eb(BY,ip(ig)  ,jp(ig)  ,kp(ig),ps) &
         &           - eb(BX,ip(ig)  ,jp(ig)  ,kp(ig),ps) &
         &           - eb(BY,ip(ig)-1,jp(ig)  ,kp(ig),ps)
!
!          currnt(ig) = 0.d0
! Bottom
!          do ii=ip(ig), ip(ig)
!            currnt(ig) = currnt(ig) + bx(ii, jp(ig)-1, kp(ig))
!          enddo
! Right
!          do jj=jp(ig), jp(ig)
!            currnt(ig) = currnt(ig) + by(ip(ig),   jj, kp(ig))
!          enddo
! Left
!          do ii=ip(ig), ip(ig)
!            currnt(ig) = currnt(ig) - bx(ii, jp(ig),   kp(ig))
!          enddo
! Up
!          do jj=jp(ig), jp(ig)
!            currnt(ig) = currnt(ig) - by(ip(ig)-1, jj, kp(ig))
!          enddo
        end do

        write(80,*) t, (currnt(ig)/renb*(cv*cv/renv/renv),ig=1,ngap)
!


!----------  calculate surface current  -----------------c
!
!        ncont = 0
!        do ipc = 1, npc 
!        if(ipc.eq.1) then
!          do zzp = nzpc1(ipc)+nz0-4, nzpc2(ipc)+nz0 
!
!            ncont = ncont+1
!            scurnt = 0.d0
! Bottom
!            do ii=nxpc1(ipc)+nx0, nxpc2(ipc)+nx0
!              scurnt = scurnt + bx(ii, nypc1(ipc)+ny0m1, zzp)
!            enddo
! Right
!            do jj=nypc1(ipc)+ny0, nypc2(ipc)+ny0
!              scurnt = scurnt + by(nxpc2(ipc)+nx0,   jj, zzp)
!            enddo
! Left
!            do ii=nxpc1(ipc)+nx0, nxpc2(ipc)+nx0
!              scurnt = scurnt - bx(ii, nypc2(ipc)+ny0,   zzp)
!            enddo
! Up
!            do jj=nypc1(ipc)+ny0, nypc2(ipc)+ny0
!              scurnt = scurnt - by(nxpc1(ipc)+nx0m1, jj, zzp)
!            enddo
!
!-------- rotB = mu J --> J = rotB / mu = rotB * cv^2 (ep=1)
!            surfj(ncont) = scurnt/renb*(cv*cv/renv/renv)
!          enddo
!
!        else if(ipc.gt.1.and.ipc.lt.npc) then
!          do zzp = nzpc1(ipc)+nz0, nzpc2(ipc)+nz0 
!
!            ncont = ncont+1
!            scurnt = 0.d0
! Bottom
!            do ii=nxpc1(ipc)+nx0, nxpc2(ipc)+nx0
!              scurnt = scurnt + bx(ii, nypc1(ipc)+ny0m1, zzp)
!            enddo
! Right
!            do jj=nypc1(ipc)+ny0, nypc2(ipc)+ny0
!              scurnt = scurnt + by(nxpc2(ipc)+nx0,   jj, zzp)
!            enddo
! Left
!            do ii=nxpc1(ipc)+nx0, nxpc2(ipc)+nx0
!              scurnt = scurnt - bx(ii, nypc2(ipc)+ny0,   zzp)
!            enddo
! Up
!            do jj=nypc1(ipc)+ny0, nypc2(ipc)+ny0
!              scurnt = scurnt - by(nxpc1(ipc)+nx0m1, jj, zzp)
!            enddo
!
!-------- rotB = mu J --> J = rotB / mu = rotB * cv^2 (ep=1)
!            surfj(ncont) = scurnt/renb*(cv*cv/renv/renv)
!          enddo
!
!        else if(ipc.eq.npc) then
!          do zzp = nzpc1(ipc)+nz0, nzpc2(ipc)+nz0m1+4
!
!            ncont = ncont+1
!            scurnt = 0.d0
! Bottom
!            do ii=nxpc1(ipc)+nx0, nxpc2(ipc)+nx0
!              scurnt = scurnt + bx(ii, nypc1(ipc)+ny0m1, zzp)
!            enddo
! Right
!            do jj=nypc1(ipc)+ny0, nypc2(ipc)+ny0
!              scurnt = scurnt + by(nxpc2(ipc)+nx0,   jj, zzp)
!            enddo
! Left
!            do ii=nxpc1(ipc)+nx0, nxpc2(ipc)+nx0
!              scurnt = scurnt - bx(ii, nypc2(ipc)+ny0,   zzp)
!            enddo
! Up
!            do jj=nypc1(ipc)+ny0, nypc2(ipc)+ny0
!              scurnt = scurnt - by(nxpc1(ipc)+nx0m1, jj, zzp)
!            enddo
!
!-------- rotB = mu J --> J = rotB / mu = rotB * cv^2 (ep=1)
!            surfj(ncont) = scurnt/renb*(cv*cv/renv/renv)
!          enddo
!
!        else
!          print*,"[vical].WARN01: ipc=",ipc
!
!        end if
!        end do
!        write(6,*) 'surfj:number of data',ncont

        close(80)
        close(81)
        if(mode_dipole.eq.2) close(82)
        if(mode_dipole.eq.3) close(83)
      end if


! - line electric field and potential -----
!--------
!
!

!        do iline = 1,8
!        do igrid = 1,ixy
!          surfe(iline,igrid) = 0.0d0
!        end do
!        end do

!        do i = nx0,nx1p1
!          surfe(1,i) = ez(i,nypc2(1)+ny0,int(nzpc1(1)+nzpc2(npc))/2+nz0)
!        end do

!        do i = nx0,nx1p1
!          surfe(2,i) = ez(i,nypc2(1)+ny0,nzpc2(1)+nz0)
!          surfe(3,i) = ez(i,nypc2(2)+ny0,nzpc2(2)+nz0)
!        end do

!        do k = nz0,nz1p1
!          surfe(6,k) = ex(nx0+nx1p1-k,nypc2(1)+ny0,k)
!        end do

!        do k = nz0,nz1p1
!          surfe(7,k) = ez(nx0+nx1p1-k,nypc2(1)+ny0,k)
!        end do

!        do k = nz0,nz1p1
!          surfe(8,k) = ez(nxpc2(1)+nx0,nypc2(1)+ny0,k)
!        end do
