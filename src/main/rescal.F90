#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine rescal
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   R E S C A L
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   in this subroutine, all the parameters except dt,dr,   .
!   .   wp(n),qm(n), are re-scaled  by dt/2, dr, and qm(1)     .
!   .   for further calculation                                .
!   ............................................................

!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  use hdf
  implicit none
!
  integer(kind=8) :: m, ns,ne
  integer(kind=4) :: i,j,k, ig, ipc, ipcg, is, iepl, ijs, itch
  integer(kind=4) :: neplsum
  integer(kind=4) :: ps
  real(kind=8) :: qm1, dth, epc1, amu1
!
  integer(kind=HID_T) :: fileid
  integer(kind=4) :: stats0,stats1
  integer(kind=8) :: dims(1)
  real(kind=8) :: rens(4),irens(4)
  character(len=30) :: filename,dsname


!-------------------- S-unit system
!   In S-unit system, four physical quatities are fixed as follows;
!       permittivity ---> 1.0d0,
!       qm(1) ---> 1.0d0,
!       dr    ---> 1.0d0,
!       dt    ---> 2.0d0.
!
        qm1     = qm(1)
        dri     = 1.0d0/dr
        dth     = dt/2.0d0


!-------------------- rescaling coefficients
!       ------------- distance r
        renr    = dri
!       ------------- time t
        rent    = 2.0d0/dt
!       ------------- velocity v
        renv    = dth*dri
!       ------------- charge q
        renq    = qm1*((dth*dri)**2)*dri
!       ------------- electric field e
        rene    = qm1*(dth**2)*dri
!       ------------- magnetic field b
        renb    = qm1*dth
!       ------------- charge density rho
        renrho  = qm1*(dth**2)
!       ------------- current density j
        renj    = qm1*(dth**3)*dri
!       ------------- mass m
        renm    = (qm1*dth*dri)**2*dri
!       ------------- electric potential phi
        renphi  = qm1*(dth**2)*(dri**2)
!       ------------- energy density
        rened   = (qm1**2)*(dth**4)*(dri**2)
!       ------------- acceleration
        rena    = renv / rent

        reni = renj * renr*renr ! Current [A]
        renreg = renphi / reni ! Registance [Ω]
        reng = 1 / renreg ! Conductance [S]
        renec = reng / renr ! Electric conductivity [S/m]

        conductivity = conductivity * renec



      if (myid == 0) then
        print *, 'renr =', renr
        print *, 'rent =', rent
        print *, 'renv =', renv
        print *, 'renq =', renq
        print *, 'rene =', rene
        print *, 'renb =', renb
        print *, 'renrho =', renrho
        print *, 'renj =', renj
        print *, 'renm =', renm
        print *, 'renphi =', renphi
        print *, 'rened =', rened
        print *, 'rena =', rena
      end if


!-------------------- 
      if(jobnum(1).gt.0) then
        if(jobnum(1).eq.1) then
          write(filename,'(a,i4.4,a)') './SNAPSHOT0/esdat', myid, '.h5'
        else
          write(filename,'(a,i4.4,a)') './SNAPSHOT0/emdat', myid, '.h5'
        end if
        call hdfopen(filename,fileid,DFACC_READ)
!
        dsname = 'rens'
        dims(1) = 4
        call read1d(fileid,dsname,dims(1:1),rens(1:4),stats0,stats1)
        irens(1:4) = 1.0d0/rens(1:4)
!        print*, "rens, irens",myid,rens,irens
!
        call hdfclose(fileid,stats0)
      else
        rens(1:4) = 1.0d0; irens(1:4) = 1.0d0
      end if



!----------------------------------!
!            re-scaling            !
!----------------------------------!
        slx     = slx*renr
        sly     = sly*renr
        slz     = slz*renr
        fftofs  = fftofs*renr
        cv      = cv*renv
        cs      = cs*renv*renv
        tcs     = tcs*renv*renv
        b0      = b0*renb
        b0x     = b0x*renb
        b0y     = b0y*renb
        b0z     = b0z*renb
        e0      = e0*rene
        e0x     = e0x*rene
        e0y     = e0y*rene
        e0z     = e0z*rene
!       ------------- gap voltage amplitude
        do ig=1,ngap
          gap_amp(ig) = gap_amp(ig)*rene
          ajs_amp(ig) = ajs_amp(ig)*renj
        end do
        Ew(:,:,:) = Ew(:,:,:)*rene

!-------------------- species loop
      do is=1,nspec
        qm(is)  = qm(is)/qm1
        q(is)   = q(is)*renq
        rm(is)  = rm(is)*renm
        path(is)= path(is)*renv
        spa(is) = spa(is)*renv
        spe(is) = spe(is)*renv
        peth(is)= peth(is)*renv
        vdx(is) = vdx(is)*renv
        vdy(is) = vdy(is)*renv
        vdz(is) = vdz(is)*renv
        vdri(is) = vdri(is)*renv
        vdtx(is)= vdtx(is)*renv
        vdty(is)= vdty(is)*renv
        vdtz(is)= vdtz(is)*renv
        vpa(is) = vpa(is)*renv
        vpb(is) = vpb(is)*renv
        vpe(is) = vpe(is)*renv
        vxyzmax(is) = vxyzmax(is)*renv
        if(myid.eq.0) print*,'is,q(is)=',is,q(is)
      end do

!-------------------- inverse position
      si      = si /renr/renr
      dkx     = dkx/renr
      dky     = dky/renr
      dkz     = dkz/renr

      kmod(:,:,2) = kmod(:,:,2)/renr/renr

      npsum = 0.0d0
      do ps=1,2
      do is=1,nspec
        npsum = npsum + totalp(is,ps)
      end do
      end do

      do m=1,pbase(3)
        pbuf(m)%x  = pbuf(m)%x*irens(1)*renr
        pbuf(m)%y  = pbuf(m)%y*irens(1)*renr
        pbuf(m)%z  = pbuf(m)%z*irens(1)*renr
        pbuf(m)%vx = pbuf(m)%vx*irens(2)*renv
        pbuf(m)%vy = pbuf(m)%vy*irens(2)*renv
        pbuf(m)%vz = pbuf(m)%vz*irens(2)*renv
      end do

      ne = 0
      do is=1,nspec
        ns = ne + 1
        ne = ne + npr(is)
        if(nflag_emit(is).eq.0) then
          do i=ns,ne
            vxf(i)  = vxf(i)*renv
            vxb(i)  = vxb(i)*renv
            vyf(i)  = vyf(i)*renv
            vyb(i)  = vyb(i)*renv
            vzf(i)  = vzf(i)*renv
            vzb(i)  = vzb(i)*renv
            vxr(i)  = vxr(i)*renv
            vyr(i)  = vyr(i)*renv
            vzr(i)  = vzr(i)*renv
          end do
        else
          do i=ns,ne
            vnormf(i)  = vnormf(i)*renv
            vtangf(i)  = vtangf(i)*renv
            vnormc(i)  = vnormc(i)*renv
            vtangc(i)  = vtangc(i)*renv
            vnorms(i)  = vnorms(i)*renv
            vtangs1(i) = vtangs1(i)*renv
            vtangs2(i) = vtangs2(i)*renv
          end do
        end if
      end do

!-------------------- normalization for fields
      eb(EX,:,:,:,:) = eb(EX,:,:,:,:)*rene
      eb(EY,:,:,:,:) = eb(EY,:,:,:,:)*rene
      eb(EZ,:,:,:,:) = eb(EZ,:,:,:,:)*rene
      eb(BX,:,:,:,:) = eb(BX,:,:,:,:)*renb
      eb(BY,:,:,:,:) = eb(BY,:,:,:,:)*renb
      eb(BZ,:,:,:,:) = eb(BZ,:,:,:,:)*renb
      aj(:,:,:,:,:) = aj(:,:,:,:,:)*renj
      rho(:,:,:,:,:) = rho(:,:,:,:,:)*renrho
      phi(:,:,:,:,:) = phi(:,:,:,:,:)*renphi
      rhobk(:,:,:,:,3) = rhobk(:,:,:,:,3)*irens(3)*renrho

!      do k=1,nzm
!      do j=1,nym
!      do i=1,nxm
!         ex(i,j,k) =  ex(i,j,k)*rene
!         ey(i,j,k) =  ey(i,j,k)*rene
!         ez(i,j,k) =  ez(i,j,k)*rene
!        pex(i,j,k) = pex(i,j,k)*rene
!        pey(i,j,k) = pey(i,j,k)*rene
!        pez(i,j,k) = pez(i,j,k)*rene
!         bx(i,j,k) =  bx(i,j,k)*renb
!         by(i,j,k) =  by(i,j,k)*renb
!         bz(i,j,k) =  bz(i,j,k)*renb
!        pbx(i,j,k) = pbx(i,j,k)*renb
!        pby(i,j,k) = pby(i,j,k)*renb
!        pbz(i,j,k) = pbz(i,j,k)*renb
!        pbx0(i,j,k) = pbx0(i,j,k)*renb
!        pby0(i,j,k) = pby0(i,j,k)*renb
!        pbz0(i,j,k) = pbz0(i,j,k)*renb
!      end do
!      end do
!      end do
!      if(jobnum.lt.0) then
!      do k=1,nzm
!      do j=1,nym
!      do i=1,nxm
!        rhoion(i,j,k) = rhoion(i,j,k)*renrho
!      end do
!      end do
!      end do
!      end if

!-------------------- normalization of accumulated charge
      do ipc=1,npc
        gcount(1)%chgacm(1,ipc) = gcount(1)%chgacm(1,ipc)*irens(4)*renq
        gcount(1)%chgacm(2,ipc) = gcount(1)%chgacm(2,ipc)*renq
        do is=1,nspec
          fluxv(is,ipc) = fluxv(is,ipc)*renq
        end do
        biasp(ipc) = biasp(ipc)*renphi
        pswspn(ipc) = pswspn(ipc)*renphi
      end do
      i = 1
      do while(biasc(i)%to.ne.0)
        biasc(i)%val = biasc(i)%val*renq/rent
        i = i + 1
      end do
      do ipcg=1,npcg
        pfixed(ipcg) = pfixed(ipcg)*renphi
      end do

!-------------------- 
      neplsum = 0
      do is=1,nspec
        neplsum = neplsum + nepl(is)
      end do
      do iepl=1,neplsum
        xwdte(iepl) = xwdte(iepl)*renr
        ywdte(iepl) = ywdte(iepl)*renr
        zwdte(iepl) = zwdte(iepl)*renr
      end do

!-------------------- 
      do ijs=1,njs
        ajs(1:3,ijs) = ajs(1:3,ijs)*renj
        wjs(ijs)     = wjs(ijs)/rent
        th0js(ijs)   = th0js(ijs)/180.0d0*pi
      end do

!-------------------- 
      do itch=1,ntch
        qtch(itch) = qtch(itch)*renq
        e1tch(itch) = e1tch(itch)*rene
        eMtch(itch) = eMtch(itch)*rene
        p1tch(itch) = p1tch(itch)*renphi
      end do

!-------------------- 
      do itch=1,nqclst
        rhopeak(itch) = rhopeak(itch)*renrho
      end do

!-------------------- 
      do is=1,nspec
        pemax(is) = pemax(is)*renv*renv
      end do

!-------------------- for FDTD
        if(myid.eq.0) print*,'dt,dth',dt,dth
        epc1    = 1.0
        amu1    = 1./(cv*cv*epc1)
        amu2    = amu1*amu2
        epc2    = epc1*epc2
        sigma   = sigma*dth*dr
      if(mode_dipole.eq.3) then
        write(6,*) "resc_org:",(resc(ig),ig=1,ngap)
        write(6,*) "capc_org:",(capc(ig),ig=1,ngap)
        do ig=1,ngap
          if(resc(ig).lt.0) then
            write(6,*) "resc < 0 --> Infinite resistance"
          else
            resc(ig) = resc(ig)/120.0d0/pi/cv
          end if
          capc(ig) = capc(ig)*120.0d0*pi*cv*rent
        end do
        write(6,*) "resc-->",(resc(ig),ig=1,ngap)
        write(6,*) "capc-->",(capc(ig),ig=1,ngap)
      end if


  return
  end subroutine
