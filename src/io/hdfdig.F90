#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine hdfdig(func,idiag)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   H D F D I G
!   ____________________________________________________________
!
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  use hdf
  use m_str
#define MCW local_comm
  implicit none
!
  integer(kind=8) :: m, nsP,nsS,neP,neS, neeP,neeS, nns,nne, nout
  integer(kind=4) :: func, idiag
  integer(kind=4) :: is, h5ind
  integer(kind=4) :: istat1,istat2
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xl,yl,zl, xu,yu,zu
  integer(kind=HSIZE_T) :: h5tab=0
  real(kind=8) :: cfbj
  real(kind=8) :: pout(1)
  real(kind=8) :: elemchg
  real(kind=8) :: massele
  real(kind=8) :: massion
  real(kind=8) :: epsiln0
  real(kind=8) :: boltfct
  real(kind=8) :: coe_v
  real(kind=8) :: coe_x
  real(kind=8) :: coe_t
  real(kind=8) :: coe_w
  real(kind=8) :: coe_ep
  real(kind=8) :: coe_qm
  real(kind=8) :: coe_j
  real(kind=8) :: coe_ph
  real(kind=8) :: coe_e
  real(kind=8) :: coe_m
  character(len=30) :: plhdfn
  character(len=10) :: grname
  character(len=8) :: dsname


!--------------------
      if(func.eq.1) then
        if(idiag.lt.hdfdigstart.or.idiag.gt.hdfdigend) return
        if(ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) lfdiag = lfdiag + 1
        if(ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) ljdiag = ljdiag + 1
        elemchg = 1.602177d-19
        massele = 9.109382d-31
        massion = massele*massratio
        epsiln0 = 8.854188d-12
        boltfct = 1.380651E-23
        if(deltax.gt.0.0d0) then
          coe_v = cvsim/2.99792458d8
          coe_x = 1.0d0/deltax/1.0d-2
          coe_t = coe_x/coe_v
          coe_w = 1.0d0/coe_t
          coe_ep = 1.0d0/epsiln0
          coe_qm = massele/elemchg
          coe_j = coe_ep/coe_qm*coe_v*coe_v*coe_v/coe_x/coe_x
          coe_ph = coe_v*coe_v/coe_qm
          coe_e = coe_v*coe_v/coe_qm/coe_x
          coe_m = coe_v/coe_qm/coe_x
        else
          coe_v = 1.0d0
          coe_x = 1.0d0
          coe_t = 1.0d0
          coe_w = 1.0d0
          coe_ep = 1.0d0
          coe_qm = 1.0d0
          coe_j = 1.0d0
          coe_ph = 1.0d0
          coe_e = 1.0d0
          coe_m = 1.0d0
        end if
      end if


!-------------------- 
      if(func.eq.1) then
        if((ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0).or.daverg.ge.1) then
          xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
          yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
          zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)
!
          phi(:,:,:,:,:) = 0.0d0
          poi(:,:,:,:,:) = 0.0d0
!
          call poisson(3)
!
          do k=zl-1,zu+1
          do j=yl-1,yu+1
          do i=xl-1,xu+1
            phi(1,i-xl,j-yl,k-zl,1) = phi(1,i-xl,j-yl,k-zl,1) &
           &                        - (i - prefcrd(1))*e0x &
           &                        - (j - prefcrd(2))*e0y &
           &                        - (k - prefcrd(3))*e0z
          end do
          end do
          end do
        end if
!
        if(ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          if(daverg.ge.1) then
            ebav(EX:BZ,:,:,:) = ebav(EX:BZ,:,:,:) &
           &                  + eb(EX:BZ,:,:,:,1)
            rhoav(1:nspec,:,:,:) = rhoav(1:nspec,:,:,:) &
           &                     + abs(rhodg(1:nspec,:,:,:,1))
            rhoav(nspec+1:nspec*2,:,:,:) = rhoav(nspec+1:nspec*2,:,:,:) &
           &                             + abs(rhodg(nspec+1:nspec*2,:,:,:,1))
            rhoav(nspec*2+1,:,:,:) = rhoav(nspec*2+1,:,:,:) &
           &                       + rho(1,:,:,:,1)
            rhoav(nspec*2+2,:,:,:) = rhoav(nspec*2+2,:,:,:) &
           &                       + rhobk(1,:,:,:,3)
            phiav(1,:,:,:) = phiav(1,:,:,:) &
           &               + phi(1,:,:,:,1)
          else if(daverg.eq.0) then
            ebav(EX:BZ,:,:,:) = eb(EX:BZ,:,:,:,1)
            rhoav(1:nspec,:,:,:) = abs(rhodg(1:nspec,:,:,:,1))
            rhoav(nspec+1:nspec*2,:,:,:) = abs(rhodg(nspec+1:nspec*2,:,:,:,1))
            rhoav(nspec*2+1,:,:,:) = rho(1,:,:,:,1)
            rhoav(nspec*2+2,:,:,:) = rhobk(1,:,:,:,3)
            phiav(1,:,:,:) = phi(1,:,:,:,1)
          end if
          h5fcount = h5fcount + 1
          ebav(:,:,:,:) = ebav(:,:,:,:)/h5fcount
          rhoav(:,:,:,:) = rhoav(:,:,:,:)/h5fcount
          phiav(:,:,:,:) = phiav(:,:,:,:)/h5fcount
        else
          if(daverg.ge.1) then
            ebav(EX:BZ,:,:,:) = ebav(EX:BZ,:,:,:) &
           &                  + eb(EX:BZ,:,:,:,1)
            rhoav(1:nspec,:,:,:) = rhoav(1:nspec,:,:,:) &
           &                     + abs(rhodg(1:nspec,:,:,:,1))
            rhoav(nspec+1:nspec*2,:,:,:) = rhoav(nspec+1:nspec*2,:,:,:) &
           &                             + abs(rhodg(nspec+1:nspec*2,:,:,:,1))
            rhoav(nspec*2+1,:,:,:) = rhoav(nspec*2+1,:,:,:) &
           &                       + rho(1,:,:,:,1)
            rhoav(nspec*2+2,:,:,:) = rhoav(nspec*2+2,:,:,:) &
           &                       + rhobk(1,:,:,:,3)
            phiav(1,:,:,:) = phiav(1,:,:,:) &
           &               + phi(1,:,:,:,1)
            h5fcount = h5fcount + 1
          end if
        end if
      end if
!
      if(func.eq.1.and.(ijxyz(JX).ge.1.or.ijxyz(JY).ge.1.or.ijxyz(JZ).ge.1)) then
        if(ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          if(daverg.eq.1) then
            ajav(1:nspec*3,:,:,:) = ajav(1:nspec*3,:,:,:) &
           &                      + ajdg(1:nspec*3,:,:,:,1)
            ajav(nspec*3+1:nspec*3+3,:,:,:) = ajav(nspec*3+1:nspec*3+3,:,:,:) &
           &                                + aj(1:3,:,:,:,1)
          else if(daverg.eq.0) then
            ajav(1:nspec*3,:,:,:) = ajdg(1:nspec*3,:,:,:,1)
            ajav(nspec*3+1:nspec*3+3,:,:,:) = aj(1:3,:,:,:,1)
          end if
          h5jcount = h5jcount + 1
          ajav(:,:,:,:) = ajav(:,:,:,:)/h5jcount
        else
          if(daverg.eq.1) then
            ajav(1:nspec*3,:,:,:) = ajav(1:nspec*3,:,:,:) &
           &                      + ajdg(1:nspec*3,:,:,:,1)
            ajav(nspec*3+1:nspec*3+3,:,:,:) = ajav(nspec*3+1:nspec*3+3,:,:,:) &
           &                                + aj(1:3,:,:,:,1)
            h5jcount = h5jcount + 1
          end if
        end if
      end if


      h5ind = 0
!-------------------- fields
!     --------------- electromagnetic field:
      if(func.eq.1) write(dsname,'(i4.4)') lfdiag
      if(ifxyz(EX).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'ex'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : ex create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FEB),h5tab, &
       &               ebav(EX,0:nxsd,0:nysd,0:nzsd)/rene/coe_e,istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : ex close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(EY).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'ey'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : ey create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FEB),h5tab, &
       &               ebav(EY,0:nxsd,0:nysd,0:nzsd)/rene/coe_e,istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : ey close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(EZ).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'ez'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : ez create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FEB),h5tab, &
       &               ebav(EZ,0:nxsd,0:nysd,0:nzsd)/rene/coe_e,istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : ez close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(BX).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'bx'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : bx create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FEB),h5tab, &
       &               ebav(BX,0:nxsd,0:nysd,0:nzsd)/renb/coe_m,istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : bx close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(BY).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'by'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : by create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FEB),h5tab, &
       &               ebav(BX,0:nxsd,0:nysd,0:nzsd)/renb/coe_m,istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : by close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(BZ).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'bz'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : bz create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FEB),h5tab, &
       &               ebav(BZ,0:nxsd,0:nysd,0:nzsd)/renb/coe_m,istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : bz close h5 file'
          end if
        end if
      end if

!     --------------- charge density + potential:
      if(func.eq.1) write(dsname,'(i4.4)') lfdiag
      if(ifxyz(7).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'rho'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : rho create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FRH),h5tab, &
         &             rhoav(nspec*2+1,0:nxsd,0:nysd,0:nzsd)/renrho, &
         &             istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : rho close h5 file'
          end if
        end if
      end if
!
      do is=1,nspec
        if(irhsp(is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'nd',is,'p'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : ndp', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
            call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FRD),h5tab, &
           &             abs(rhoav(nspec*0+is,0:nxsd,0:nysd,0:nzsd)/renrho/rho0), &
           &             istat1,istat2)
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : ndp', is, 'close h5 file'
            end if
          end if
        end if
      end do
!
      if(ifxyz(7).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'rhobk'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : rhobk create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FRH),h5tab, &
         &             rhoav(nspec*2+2,0:nxsd,0:nysd,0:nzsd)/renrho, &
         &             istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : rhobk close h5 file'
          end if
        end if
      end if

      if(ifxyz(7).ge.1) then
        do is = 1, nspec
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a)') 'rhobksp'//str(is)
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : rhobksp create h5 file'
            end if
            print *, 'rhobksp h5 created', grname, plhdfn 
          else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
            call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FRH),h5tab, &
                        rhobksp(1,0:nxsd,0:nysd,0:nzsd,3,is)/renrho, &
                        istat1,istat2)
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : rhobksp close h5 file'
            end if
          end if
        end do
      end if
!
      if(ifxyz(7).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'phisp'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : phisp create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FPH),h5tab, &
         &             phiav(1,0:nxsd,0:nysd,0:nzsd)/renphi/coe_ph, &
         &             istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : phisp close h5 file'
          end if
        end if
      end if


!     --------------- current density:
      if(func.eq.1) write(dsname,'(i4.4)') ljdiag
      if(ijxyz(JX).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'jx'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jx create h5 file'
          end if
        else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FAJ),h5tab, &
         &             aj(JX,0:nxsd,0:nysd,0:nzsd,1)/renj/coe_j,istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : jx close h5 file'
          end if
        end if
!
        do is=1,nspec
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'j',is,'x'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jx', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
            call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FAJ),h5tab, &
           &             ajav((is-1)*3+JX,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,istat1,istat2)
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : jx', is, 'close h5 file'
            end if
          end if
        end do
      end if
!
      if(ijxyz(JY).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'jy'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jy create h5 file'
          end if
        else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FAJ),h5tab, &
         &             aj(JY,0:nxsd,0:nysd,0:nzsd,1)/renj/coe_j,istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : jy close h5 file'
          end if
        end if
!
        do is=1,nspec
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'j',is,'y'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jy', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
            call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FAJ),h5tab, &
           &             ajav((is-1)*3+JY,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,istat1,istat2)
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : jy', is, 'close h5 file'
            end if
          end if
        end do
      end if
!
      if(ijxyz(JZ).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'jz'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jz create h5 file'
          end if
        else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FAJ),h5tab, &
         &             aj(JZ,0:nxsd,0:nysd,0:nzsd,1)/renj/coe_j,istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : jz close h5 file'
          end if
        end if
!
        do is=1,nspec
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'j',is,'z'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jz', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
            call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FAJ),h5tab, &
           &             ajav((is-1)*3+JZ,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,istat1,istat2)
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : jz', is, 'close h5 file'
            end if
          end if
        end do
      end if
!
      if(ijxyz(JS).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'js'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : js create h5 file'
          end if
        else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FRD),h5tab, &
         &             sum(rhoav(nspec*1+1:nspec*1+minspec,0:nxsd,0:nysd,0:nzsd),dim=1) &
         &                *0.5d0/renj/coe_j,istat1,istat2)
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : js close h5 file'
          end if
        end if
!
        do is=1,nspec
          if(irhsp(is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'j',is,'s'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : js', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
            call wrt3d_p(id_sd(2,h5ind),dsname,myid,nodes,idimfd(:,FRD),h5tab, &
           &             abs(rhoav(nspec*1+is,0:nxsd,0:nysd,0:nzsd)*0.5d0/renj/coe_j), &
           &             istat1,istat2)
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : js', is, 'close h5 file'
            end if
          end if
          end if
        end do
      end if


!-------------------- 
!     --------------- clear rho/phiave
      if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
        if(myid.eq.0) print*,"h5fcount =",h5fcount
        h5fcount = 0
        ebav(:,:,:,:) = 0.0d0
        rhoav(:,:,:,:) = 0.0d0
        phiav(:,:,:,:) = 0.0d0
      end if
!
!     --------------- clear ajav
      if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
        if(myid.eq.0) print*,"h5jcount =",h5jcount
        h5jcount = 0
        ajav(:,:,:,:) = 0.0d0
      end if


!-------------------- other components
!     --------------- antenna surface current/e-field:
      if(iadiag.ge.1) then
!        write(plhdfn,'(a,i2.2,a)') 'surfj',istep/i2hdf,'.h5'
!        call hdfopen_pg(plhdfn,grname,id_sd(24),DFACC_CREATE,MCW)
!        if(id_sd(24).le.0) then
!          write(*,*) id_sd(24), plhdfn, ' : surfj create h5 file'
!        end if
!
!        write(plhdfn,'(a,i2.2,a)') 'surfe',istep/i2hdf,'.h5'
!        call hdfopen_pg(plhdfn,grname,id_sd(25),DFACC_CREATE,MCW)
!        if(id_sd(25).le.0) then
!          write(*,*) id_sd(25), plhdfn, ' : surfe create h5 file'
!        end if
!
!        write(plhdfn,'(a,i2.2,a)') 'phi1d',istep/i2hdf,'.h5'
!        call hdfopen_pg(plhdfn,grname,id_sd(26),DFACC_CREATE,MCW)
!        if(id_sd(26).le.0) then
!          write(*,*) id_sd(26), plhdfn, ' : phi1d create h5 file'
!        end if
!
!        write(plhdfn,'(a,i2.2,a)') 'surfeU',istep/i2hdf,'.h5'
!        call hdfopen_pg(plhdfn,grname,id_sd(27),DFACC_CREATE,MCW)
!        if(id_sd(27).le.0) then
!          write(*,*) id_sd(27), plhdfn, ' : surfeU create h5 file'
!        end if
!
!        write(plhdfn,'(a,i2.2,a)') 'surfeT',istep/i2hdf,'.h5'
!        call hdfopen_pg(plhdfn,grname,id_sd(28),DFACC_CREATE,MCW)
!        if(id_sd(28).le.0) then
!          write(*,*) id_sd(28), plhdfn, ' : surfeT create h5 file'
!        end if
      end if


!-------------------- particles
!     --------------- particle positions + velocities:
      neeP = pbase(1)
      neeS = pbase(2)
      nne = 0
      do is=1,nspec
        if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
          nsP = neeP + 1
          neP = neeP + totalp(is,1)
          neeP = neeP + totalp(is,1)
          nsS = neeS + 1
          neS = neeS + totalp(is,2)
          neeS = neeS + totalp(is,2)
          nns = nne + 1
          nne = nne + npr(is)
          lpdiag(is) = lpdiag(is) + 1
        end if
!
        if(ipaxyz(1,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a,i4.4,a)') 'p',is,'xe',istep/i3hdf,'_',nfsnap,'_',myid,'.h5'
            call hdfopen(plhdfn,id_sd(1,h5ind),DFACC_CREATE)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : px', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            nout = int((neP - nsP + ipadig(is))/ipadig(is)) &
           &     + int((neS - nsS + ipadig(is))/ipadig(is))
            if(nout.gt.0) then
              write(dsname,'(a,i4.4)') 'xe',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/), &
             &      (/pbuf(nsP:neP:ipadig(is))%x/renr, &
             &        pbuf(nsS:neS:ipadig(is))%x/renr/), &
             &      istat1,istat2)
            else if(nflag_emit(is).eq.0.and.idiag.eq.0) then
              write(dsname,'(a,i4.4)') 'vxr',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/),vxr(nns:nne:ipadig(is))/renv, &
             &      istat1,istat2)
            else
              nout = 1
              pout(1) = dmiss
              write(dsname,'(a,i4.4)') 'xe',lpdiag(is)
              call wrt1d(id_sd(1,h5ind),dsname,(/nout/),pout(1),istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose(id_sd(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0) then
              write(*,*) h5ind, istats(1,h5ind), ' : px', is, 'close h5 file'
            end if
          end if
        end if
!
        if(ipaxyz(2,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a,i4.4,a)') 'p',is,'ye',istep/i3hdf,'_',nfsnap,'_',myid,'.h5'
            call hdfopen(plhdfn,id_sd(1,h5ind),DFACC_CREATE)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : py', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            nout = int((neP - nsP + ipadig(is))/ipadig(is)) &
           &     + int((neS - nsS + ipadig(is))/ipadig(is))
            if(nout.gt.0) then
              write(dsname,'(a,i4.4)') 'ye',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/), &
             &      (/pbuf(nsP:neP:ipadig(is))%y/renr, &
             &        pbuf(nsS:neS:ipadig(is))%y/renr/), &
             &      istat1,istat2)
            else if(nflag_emit(is).eq.0.and.idiag.eq.0) then
              write(dsname,'(a,i4.4)') 'vyr',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/),vyr(nns:nne:ipadig(is))/renv, &
             &      istat1,istat2)
            else
              nout = 1
              pout(1) = dmiss
              write(dsname,'(a,i4.4)') 'ye',lpdiag(is)
              call wrt1d(id_sd(1,h5ind),dsname,(/nout/),pout(1),istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose(id_sd(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0) then
              write(*,*) h5ind, istats(1,h5ind), ' : py', is, 'close h5 file'
            end if
          end if
        end if
!
        if(ipaxyz(3,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a,i4.4,a)') 'p',is,'ze',istep/i3hdf,'_',nfsnap,'_',myid,'.h5'
            call hdfopen(plhdfn,id_sd(1,h5ind),DFACC_CREATE)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : pz', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            nout = int((neP - nsP + ipadig(is))/ipadig(is)) &
           &     + int((neS - nsS + ipadig(is))/ipadig(is))
            if(nout.gt.0) then
              write(dsname,'(a,i4.4)') 'ze',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/), &
             &      (/pbuf(nsP:neP:ipadig(is))%z/renr, &
             &        pbuf(nsS:neS:ipadig(is))%z/renr/), &
             &      istat1,istat2)
            else if(nflag_emit(is).eq.0.and.idiag.eq.0) then
              write(dsname,'(a,i4.4)') 'vzr',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/),vzr(nns:nne:ipadig(is))/renv, &
             &      istat1,istat2)
            else
              nout = 1
              pout(1) = dmiss
              write(dsname,'(a,i4.4)') 'ze',lpdiag(is)
              call wrt1d(id_sd(1,h5ind),dsname,(/nout/),pout(1),istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose(id_sd(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0) then
              write(*,*) h5ind, istats(1,h5ind), ' : pz', is, 'close h5 file'
            end if
          end if
        end if
!
        if(ipaxyz(4,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a,i4.4,a)') 'p',is,'vxe',istep/i3hdf,'_',nfsnap,'_',myid,'.h5'
            call hdfopen(plhdfn,id_sd(1,h5ind),DFACC_CREATE)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : pvx', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            nout = int((neP - nsP + ipadig(is))/ipadig(is)) &
           &     + int((neS - nsS + ipadig(is))/ipadig(is))
            if(nout.gt.0) then
              write(dsname,'(a,i4.4)') 'vxe',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/), &
             &      (/pbuf(nsP:neP:ipadig(is))%vx/renv, &
             &        pbuf(nsS:neS:ipadig(is))%vx/renv/), &
             &      istat1,istat2)
            else if(nflag_emit(is).eq.0.and.idiag.eq.0) then
              write(dsname,'(a,i4.4)') 'vxf',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/),vxf(nns:nne:ipadig(is))/renv, &
             &      istat1,istat2)
            else
              nout = 1
              pout(1) = dmiss
              write(dsname,'(a,i4.4)') 'vxe',lpdiag(is)
              call wrt1d(id_sd(1,h5ind),dsname,(/nout/),pout(1),istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose(id_sd(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0) then
              write(*,*) h5ind, istats(1,h5ind), ' : pvx', is, 'close h5 file'
            end if
          end if
        end if
!
        if(ipaxyz(5,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a,i4.4,a)') 'p',is,'vye',istep/i3hdf,'_',nfsnap,'_',myid,'.h5'
            call hdfopen(plhdfn,id_sd(1,h5ind),DFACC_CREATE)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : pvy', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            nout = int((neP - nsP + ipadig(is))/ipadig(is)) &
           &     + int((neS - nsS + ipadig(is))/ipadig(is))
            if(nout.gt.0) then
              write(dsname,'(a,i4.4)') 'vye',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/), &
             &      (/pbuf(nsP:neP:ipadig(is))%vy/renv, &
             &        pbuf(nsS:neS:ipadig(is))%vy/renv/), &
             &      istat1,istat2)
            else if(nflag_emit(is).eq.0.and.idiag.eq.0) then
              write(dsname,'(a,i4.4)') 'vyf',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/),vyf(nns:nne:ipadig(is))/renv, &
             &      istat1,istat2)
            else
              nout = 1
              pout(1) = dmiss
              write(dsname,'(a,i4.4)') 'vye',lpdiag(is)
              call wrt1d(id_sd(1,h5ind),dsname,(/nout/),pout(1),istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose(id_sd(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0) then
              write(*,*) h5ind, istats(1,h5ind), ' : pvy', is, 'close h5 file'
            end if
          end if
        end if
!
        if(ipaxyz(6,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a,i4.4,a)') 'p',is,'vze',istep/i3hdf,'_',nfsnap,'_',myid,'.h5'
            call hdfopen(plhdfn,id_sd(1,h5ind),DFACC_CREATE)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : pvz', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            nout = int((neP - nsP + ipadig(is))/ipadig(is)) &
           &     + int((neS - nsS + ipadig(is))/ipadig(is))
            if(nout.gt.0) then
              write(dsname,'(a,i4.4)') 'vze',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/), &
             &      (/pbuf(nsP:neP:ipadig(is))%vz/renv, &
             &        pbuf(nsS:neS:ipadig(is))%vz/renv/), &
             &      istat1,istat2)
            else if(nflag_emit(is).eq.0.and.idiag.eq.0) then
              write(dsname,'(a,i4.4)') 'vzf',lpdiag(is)
              call wrt1d &
             &     (id_sd(1,h5ind),dsname,(/nout/),vzf(nns:nne:ipadig(is))/renv, &
             &      istat1,istat2)
            else
              nout = 1
              pout(1) = dmiss
              write(dsname,'(a,i4.4)') 'vze',lpdiag(is)
              call wrt1d(id_sd(1,h5ind),dsname,(/nout/),pout(1),istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose(id_sd(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0) then
              write(*,*) h5ind, istats(1,h5ind), ' : pvz', is, 'close h5 file'
            end if
          end if
        end if
      end do


  return
  end subroutine hdfdig
