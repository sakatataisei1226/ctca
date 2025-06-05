#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine bdyfld(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   B D Y F L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   subroutine to update e-field by current contribution   .
!   ............................................................
!
!-------------------- parameters and variables
  use oh_type
  use paramt
  use allcom
#define MCW local_comm
#define MSS MPI_STATUS_SIZE
  implicit none
!
  integer(kind=4) :: i,j,k, ii,jj 
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ipc, iline, igrid, ncont
  integer(kind=4) :: ibedge
  integer(kind=4) :: ig, ip(ingap),jp(ingap),kp(ingap)
  integer(kind=4) :: ps
  integer(kind=4) :: inode, ierr
  real(kind=8) :: egap, jgap
  real(kind=8) :: tlarge, t_over_tlarge, scurnt, zzp
  real(kind=8) :: currnt(ingap)
  real(kind=8) :: gap_amp_tmp(ingap)


!-------------------- test particle simulation
      if(juncan.ge.1000) return


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- 
      do k=0,zu-zl
      do j=0,yu-yl
      do i=0,xu-xl
        eb(EX,i,j,k,ps) = mp(EX,i,j,k,ps)*eb(EX,i,j,k,ps) - e0x
        eb(EY,i,j,k,ps) = mp(EY,i,j,k,ps)*eb(EY,i,j,k,ps) - e0y
        eb(EZ,i,j,k,ps) = mp(EZ,i,j,k,ps)*eb(EZ,i,j,k,ps) - e0z
      end do
      end do
      end do


!-------------------- artificial handling of e-field
!  One can change the values of e-field in this section if necessary.
!  For example in the case,
!    *he wants to place the currnt source,
!    *he wants to feed the gap voltage at a certain grid,
!    *and so on.............
!


!  *******************************
!  *antenna gap voltage treatment*
!  *******************************
!start--->

!-------------------- gap volt start
      if(mode_dipole.eq.1) then
!        print*,'w_c  ====???',  w_c
!        do ig=1,ngap
!          kp(ig)=k_gap(ig)-nzoffs+nz0
!          print*,'ig,i_gap,j_gap,k_gap',ig,i_gap(ig),j_gap(ig),k_gap(ig)
!          print*,'gap_amp',gap_amp(ig)
!        end do
!        print*,'nxpc1(1), nypc1(1) ', nxpc1(1), nypc1(1)
!        print*,'nxpc2(1), nypc2(1) ', nxpc2(1), nypc2(1)

!       ------------- oscillation with a specific frequency
!        if(w_c.lt.0.0)  then
!          if(istep.ge.ifeedst-nstp_oe.and.istep.lt.ifeedst) then
!            if(istep.eq.ifeedst-nstp_oe) tfeedst = t-dt
!            do ig=1,ngap
!              gap_amp_tmp(ig) = 0.5d0*gap_amp(ig) &
!       &                      *(1-dcos(pi*(t-tfeedst)/dt/nstp_oe))
!            end do
!          else if(istep.ge.ifeedst) then
!            do ig=1,ngap
!              gap_amp_tmp(ig) = gap_amp(ig)
!            end do
!          else
!            do ig=1,ngap
!              gap_amp_tmp(ig) = 0.0d0
!            end do
!          end if

!          do ig=1,ngap
!          do ip = nxpc1(1)-nxoffs+nx0, nxpc2(1)-nxoffs+nx0
!          do jp = nypc1(1)-nyoffs+ny0, nypc2(1)-nyoffs+ny0
!            if(ip    .ge.nx0m1.and.ip    .le.nx1p1.and. &
!         &     jp    .ge.ny0m1.and.jp    .le.ny1p1.and. &
!         &     kp(ig).ge.nz0m1.and.kp(ig).le.nz1p1) then
!              pez(ip,jp,kp(ig))=pez(ip,jp,kp(ig))-ez(ip,jp,kp(ig))
!              ez(ip,jp,kp(ig))=gap_amp_tmp(ig)*dsin(-w_c*t)
!              pez(ip,jp,kp(ig))=pez(ip,jp,kp(ig))+ez(ip,jp,kp(ig))
!              print*,'myid,ig,ez(ip,jp,kp)', &
!         &           myid,ig,ez(ip+nxoffs,jp+nyoffs,kp(ig)+nzoffs)
!            end if
!          end do
!          end do
!          end do
!#include "vical.fnc"
!        end if

!       ------------- gaussian pulse
!           center frequency : f_c = w_pulse * 2.*pi
!        if(w_c.gt.0.0) then
!          if(istep.ge.ifeedst-nstp_oe+1.and.istep.le.ifeedst) then
!            do ig=1,ngap
!              inode=int(nxpc1(1)/nx)+int(nypc1(1)/ny)*nodes(1) &
!           &       +int(k_gap(ig)/nz)*nodes(2)
!              if(myid.eq.inode) then
!                egap=ez(nxpc1(1)-nxoffs+nx0,nypc1(1)-nyoffs+ny0,kp(ig))
!              end if
!              call MPI_Bcast(egap, 1, MPI_REAL8, inode, MCW, ierr)
!              annrms(ig) = annrms(ig) + egap**2.0d0
!              offez(ig) = offez(ig) + egap
!            end do
!          end if
!          if(istep.eq.ifeedst) then
!            tfeedst = t
!            do ig=1,ngap
!              annrms(ig) = annrms(ig)/nstp_oe !rms**2
!              offez(ig) = offez(ig)/nstp_oe !am
!              annrms(ig) = dsqrt(annrms(ig)-offez(ig)**2) !sdv
!              if(pgamp(ig).gt.0) then
!                gap_amp(ig) = annrms(ig)*pgamp(ig)
!                if(myid.eq.0) &
!           &      print*,'modified gap_amp: ',gap_amp(ig),' for ig=',ig
!              end if
!            end do
!            if(myid.eq.0) print*,'tfeedst: ',tfeedst
!            if(myid.eq.0) print*,'gap voltage offez: ',(offez(ig),ig=1,ngap)
!            if(myid.eq.0) &
!         &    print*,'RMS of v-noise annrms: ',(annrms(ig),ig=1,ngap)
!            if(myid.eq.0) print*,'gap voltage feeding: start'
!          end if
!
!          if(istep.ge.ifeedst) then
!            f_c = w_c/(2.0*pi)
!            tlarge = 0.152 / f_c
!            t_over_tlarge =  (t-tfeedst)/tlarge
!
!            if(t_over_tlarge.ge.0.0d0) then
!              do ig=1,ngap
!              do ip = nxpc1(1)-nxoffs+nx0, nxpc2(1)-nxoffs+nx0
!              do jp = nypc1(1)-nyoffs+ny0, nypc2(1)-nyoffs+ny0
!                if(ip    .ge.nx0m1.and.ip    .le.nx1p1.and. &
!             &     jp    .ge.ny0m1.and.jp    .le.ny1p1.and. &
!             &     kp(ig).ge.nz0m1.and.kp(ig).le.nz1p1) then
!                  pez(ip,jp,kp(ig))=pez(ip,jp,kp(ig))-ez(ip,jp,kp(ig))
!                  ez(ip,jp,kp(ig))=offez(ig) + gap_amp(ig) &
!             &      *(t_over_tlarge**4.d0*(-1.0/tlarge)*dexp(-t_over_tlarge) &
!             &      +4.d0*(t_over_tlarge**3.d0)/tlarge*dexp(-t_over_tlarge))
!                  pez(ip,jp,kp(ig))=pez(ip,jp,kp(ig))+ez(ip,jp,kp(ig))
!                end if
!              end do
!              end do
!              end do
!            end if
!          end if
!
!#include "vical.fnc"
!        end if 
!
!        if(w_c.eq.0.0) then
!#include "vical.fnc"
!        end if


!-------------------- gap volt start
      else if(mode_dipole.eq.2) then
!        print*,'w_c  ====???',  w_c
        do ig=1,ngap
          ip(ig) = i_gap(ig) - xl
          jp(ig) = j_gap(ig) - yl
          kp(ig) = k_gap(ig) - zl
!          print*,'ig,i_gap,j_gap,k_gap',ig,i_gap(ig),j_gap(ig),k_gap(ig)
!          print*,'gap_amp',gap_amp(ig)
        end do
!        print*,'nxpc1(1), nypc1(1) ', nxpc1(1), nypc1(1)
!        print*,'nxpc2(1), nypc2(1) ', nxpc2(1), nypc2(1)

!       ------------- oscillation with a specific frequency
!        if(w_c.lt.0.0)  then
!          do ig=1,ngap
!          do ip = nxpc1(1)-nxoffs+nx0, nxpc2(1)-nxoffs+nx0
!          do jp = nypc1(1)-nyoffs+ny0, nypc2(1)-nyoffs+ny0
!            if(ip    .ge.nx0m1.and.ip    .le.nx1p1.and. &
!         &     jp    .ge.ny0m1.and.jp    .le.ny1p1.and. &
!         &     kp(ig).ge.nz0m1.and.kp(ig).le.nz1p1) then
!            ez(ip,jp,kp(ig))=ez(ip,jp,kp(ig))-2.0d0*ajs_amp(ig)*dsin(-w_c*t)
!            pez(ip,jp,kp(ig))=pez(ip,jp,kp(ig))-2.0d0*ajs_amp(ig)*dsin(-w_c*t)
!            print*,'myid,ig,ez(ip,jp,kp)', &
!         &         myid,ig,ez(ip+nxoffs,jp+nyoffs,kp(ig)+nzoffs)
!            end if
!          end do
!          end do
!          end do
!#include "vical.fnc"
!        end if

!       ------------- gaussian pulse
!           center frequency : f_c = w_pulse * 2.*pi
        if(w_c.gt.0.0) then
!          if(istep.ge.ifeedst-nstp_oe+1.and.istep.le.ifeedst) then
!            do ig=1,ngap
!              inode=int(nxpc1(1)/nx)+int(nypc1(1)/ny)*nodes(1) &
!           &       +int(k_gap(ig)/nz)*nodes(2)
!              if(myid.eq.inode) then
!                jgap=ajz(nxpc1(1)-xl,nypc1(1)-yl,kp(ig))
!              end if
!              call MPI_Bcast(jgap, 1, MPI_REAL8, inode, MCW, ierr)
!              annrms(ig) = annrms(ig) + jgap**2.0d0
!            end do
!          end if
          if(istep.eq.ifeedst) then
            tfeedst = t
!            do ig=1,ngap
!              annrms(ig) = dsqrt(annrms(ig)/nstp_oe)
!              if(pgamp(ig).gt.0) then
!                ajs_amp(ig) = annrms(ig)*pgamp(ig)
!                if(myid.eq.0) &
!           &      print*,'modified ajs_amp: ',ajs_amp(ig),' for ig=',ig
!              end if
!            end do
            if(myid.eq.0) print*,'tfeedst: ',tfeedst
!            if(myid.eq.0) &
!           &  print*,'RMS of j-noise annrms: ',(annrms(ig),ig=1,ngap)
            if(myid.eq.0) print*,'gap voltage feeding: start'
            do ig=1,ngap
              if(ps.eq.1.and. &
             &   ip(ig).gt.0.and.ip(ig).le.xu-xl.and. &
             &   jp(ig).gt.0.and.jp(ig).le.yu-yl.and. &
             &   kp(ig).gt.0.and.kp(ig).le.zu-zl) then
!              if(ps.eq.1) then
                print*,'myid,ig,{i,j,k}p',myid,ig,ip(ig),jp(ig),kp(ig)
              end if
            end do
          end if

          if(istep.ge.ifeedst) then
            f_c = w_c/(2.0*pi)
            tlarge = 0.152/f_c
            t_over_tlarge = (t - tfeedst)/tlarge

            if(t_over_tlarge.ge.0.0d0) then
              do ig=1,ngap
                gap_jsrc(ig) = ajs_amp(ig) &
               &              *(t_over_tlarge**4.d0*(-1.0/tlarge) &
               &              *dexp(-t_over_tlarge) &
               &             + 4.d0*(t_over_tlarge**3.d0)/tlarge &
               &              *dexp(-t_over_tlarge))
!                gap_jsrc(ig) = ajs_amp(ig) &
!               &              *(t_over_tlarge**4.d0/(tlarge**2.0d0) &
!               &              *dexp(-t_over_tlarge) &
!               &             - 8.0d0*(t_over_tlarge**3.d0)/(tlarge**2.0d0) &
!               &              *dexp(-t_over_tlarge) &
!               &             + 12.0d0*(t_over_tlarge**2.d0)/(tlarge**2.0d0) &
!               &              *dexp(-t_over_tlarge))
                if(ip(ig).ge.0.and.ip(ig).le.xu-xl.and. &
               &   jp(ig).ge.0.and.jp(ig).le.yu-yl.and. &
               &   kp(ig).ge.0.and.kp(ig).le.zu-zl) then
                  eb(EZ,ip(ig),jp(ig),kp(ig),ps) = &
               &  eb(EZ,ip(ig),jp(ig),kp(ig),ps) - 2.0d0*gap_jsrc(ig)
                end if
              end do
            end if
          end if

#include "vical.fnc"
        end if 
!
!        if(w_c.eq.0.0) then
!#include "vical.fnc"
!        end if
      end if


!-------------------- gap volt start
!      if(mode_dipole.eq.3) then
!        print*,'w_c  ====???',  w_c
!        do ig=1,ngap
!          kp(ig)=k_gap(ig)-nzoffs+nz0
!          print*,'ig,i_gap,j_gap,k_gap',ig,i_gap(ig),j_gap(ig),k_gap(ig)
!          print*,'gap_amp',gap_amp(ig)
!        end do
!        print*,'nxpc1(1), nypc1(1) ', nxpc1(1), nypc1(1)
!        print*,'nxpc2(1), nypc2(1) ', nxpc2(1), nypc2(1)
!
!#include "vical.fnc"
!      end if


!<---end "antenna gap voltage treatment"


  return
  end subroutine
