#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine frmttd(ustep)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   F R M T T D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine supervises particle injection from      .
!   .  the outer/inner boundaries of the simulation box.       .
!   ............................................................
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
#define MCW local_comm
#define MSS MPI_STATUS_SIZE
  implicit none
!
  integer(kind=4) :: nintgrtdstep = 4000
  integer(kind=4) :: ustep, ndnstep
  integer(kind=4) :: ipc,is
  integer(kind=4) :: minpc
  real(kind=8) :: icur(3,inpc,ispec), ocur(3,inpc,ispec)
  real(kind=8) :: expema
  real(kind=8) :: tew
  character(len=40) :: sformat


!-------------------- 
      if(nstep.gt.0) then
        ndnstep = int(log10(real(nstep*2))) + 1
      else
        ndnstep = 6
      end if


!-------------------- 
      if(istep.gt.0.and.emarlxt.gt.0) then
        expema = exp(-1.0d0/emarlxt)
      else
        expema = 0.0d0
      end if


!-------------------- 
      minpc = max(npc,1)


!-------------------- diagnostics
      if(myid.eq.0) then
        if(istep.eq.0.and.nintgrtdstep.gt.nstep) then
          if(nstep.gt.0) then
            nintgrtdstep = nstep
          else
            nintgrtdstep = 1
          end if
        end if
        if(istep.eq.0) then
          icur(:,:,:) = 0.0d0
          ocur(:,:,:) = 0.0d0
        end if
        if(istep.eq.nstep-nintgrtdstep) then
          selfpintgrtd(1:npc+1) = 0.0d0
          outfluxintgrtd(1:npc,1:nspec) = 0
          infhistintgrtd(0:npc,1:npc,1:nspec) = 0
          rhoindintgrtd(1:npc) = 0.0d0
          vidata(1:10) = 0.0d0
        end if
        if(istep.gt.nstep-nintgrtdstep) then
          selfpintgrtd(1:npc) = selfpintgrtd(1:npc) + selfp(1:npc)
          selfpintgrtd(npc+1) = selfpintgrtd(npc+1) + phiref(2)
          outfluxintgrtd(1:npc,1:nspec) = &
         &  outfluxintgrtd(1:npc,1:nspec) &
         &  + gcount(2)%outflux(1:npc,1:nspec)
          infhistintgrtd(0:npc,1:npc,1:nspec) = &
         &  infhistintgrtd(0:npc,1:npc,1:nspec) &
         &  + gcount(2)%infhist(0:npc,1:npc,1:nspec)
          rhoindintgrtd(1:npc) = &
         &  rhoindintgrtd(1:npc) &
         &  + rhoind(1:npc)/wrelax - rhoindh(1:npc)/wrelax - gcount(2)%chgacm(2,1:npc)
!          vidata( 1) = vidata( 1) + (selfp(1)+selfp(npc))*0.5d0
!          vidata( 2) = vidata( 2) + (biasc(1)%val+biasc(2)%val)*0.5d0
!          vidata( 3) = vidata( 3) + (sum(gcount(2)%outflux(1:2,nspec)) &
!         &                          +sum(gcount(2)%outflux(npc-1:npc,nspec)) &
!         &                          -sum(gcount(2)%infhist(1,1:npc,nspec)) &
!         &                          -sum(gcount(2)%infhist(2,1:npc,nspec)) &
!         &                          -sum(gcount(2)%infhist(npc-1,1:npc,nspec)) &
!         &                          -sum(gcount(2)%infhist(npc,1:npc,nspec)))*q(3)*0.5d0
!          vidata( 4) = vidata( 4) - (sum(gcount(2)%infhist(0,1:2,1)) &
!         &                          +sum(gcount(2)%infhist(0,npc-1:npc,1)))*q(1)*0.5d0
!          vidata( 5) = vidata( 5) - (sum(gcount(2)%infhist(0,1:2,2)) &
!         &                          +sum(gcount(2)%infhist(0,npc-1:npc,2)))*q(2)*0.5d0
!          vidata( 6) = vidata( 6) + (sum(gcount(2)%infhist(1,3:npc,nspec)) &
!         &                          +sum(gcount(2)%infhist(2,3:npc,nspec)) &
!         &                          +sum(gcount(2)%infhist(npc-1,1:npc-2,nspec)) &
!         &                          +sum(gcount(2)%infhist(npc,1:npc-2,nspec)))*q(3)*0.5d0
!          vidata( 7) = vidata( 7) - (sum(gcount(2)%infhist(3:npc,1,nspec)) &
!         &                          +sum(gcount(2)%infhist(3:npc,2,nspec)) &
!         &                          +sum(gcount(2)%infhist(1:npc-2,npc-1,nspec)) &
!         &                          +sum(gcount(2)%infhist(1:npc-2,npc,nspec)))*q(3)*0.5d0
!          vidata( 8) = vidata( 8) - (   (gcount(2)%infhist(3,1,nspec)) &
!         &                          +   (gcount(2)%infhist(3,2,nspec)) &
!         &                          +   (gcount(2)%infhist(npc-2,npc-1,nspec)) &
!         &                          +   (gcount(2)%infhist(npc-2,npc,nspec)))*q(3)*0.5d0
!          vidata( 9) = vidata( 9) - (   (gcount(2)%infhist(4,1,nspec)) &
!         &                          +   (gcount(2)%infhist(4,2,nspec)) &
!         &                          +   (gcount(2)%infhist(npc-3,npc-1,nspec)) &
!         &                          +   (gcount(2)%infhist(npc-3,npc,nspec)))*q(3)*0.5d0
!          vidata(10) = vidata(10) - (sum(gcount(2)%infhist(5:npc-4,1,nspec)) &
!         &                          +sum(gcount(2)%infhist(5:npc-4,2,nspec)) &
!         &                          +sum(gcount(2)%infhist(5:npc-4,npc-1,nspec)) &
!         &                          +sum(gcount(2)%infhist(5:npc-4,npc,nspec)))*q(3)*0.5d0
        end if
        if(istep.eq.nstep) then
          selfpintgrtd(1:npc+1) = selfpintgrtd(1:npc+1)/nintgrtdstep/renphi
          outfluxintgrtd(1:npc,1:nspec) = &
         &  outfluxintgrtd(1:npc,1:nspec)/nintgrtdstep
          infhistintgrtd(0:npc,1:npc,1:nspec) = &
         &  infhistintgrtd(0:npc,1:npc,1:nspec)/nintgrtdstep
          rhoindintgrtd(1:npc) = &
         &  rhoindintgrtd(1:npc)/nintgrtdstep/renq*rent
!          vidata(1) = vidata(1)/nintgrtdstep/renphi
!          vidata(2:10) = vidata(2:10)/nintgrtdstep/renq
        end if
!
        if(intfoc.ne.0.and.mod(istep,intfoc).eq.0) then
          print*,"  globalp =",gcount(2)%globalp(1:nspec,:)
        end if
        if(intfoc.gt.0) then
          if(istep.eq.0.or.mod(istep-1,intfoc).eq.0) then
            open(79,file='energy',position='append')
            open(86,file='pbody',position='append')
            open(89,file='influx',position='append')
            open(90,file='noflux',position='append')
            open(91,file='isflux',position='append')
            open(92,file='nesc',position='append')
            open(93,file='chgacm1',position='append')
            open(94,file='chgacm2',position='append')
            open(95,file='chgmov',position='append')
            open(96,file='seyield',position='append')
            open(97,file='icur',position='append')
            open(98,file='ocur',position='append')
            open(99,file='ewave',position='append')
          end if
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') 2+3*nspec
          write(79,sformat) istep, engebg(1:2)/rened, engkg(1:3,1:nspec)/rened
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') npc+1
          write(86,sformat) istep, selfp(1:npc)/renphi, phiref(2)/renphi
!          if(istep.eq.nstep) &
!         &  write(86,sformat) istep*2, selfpintgrtd(1:npc+1)
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(I,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') minpc*minspec
          write(89,sformat) istep, int(gcount(2)%influx(1:minpc,1:minspec)*2/ustep)
!          if(istep.eq.nstep) &
!         &  write(89,sformat) istep*2, influxintgrtd(1:npc,1:nspec)
!         -----------
          sformat = '(xxxx(Ix.x,xxxx(1X,ES9.4e1),2X))'
          write(sformat(2:5),'(I4.4)') minspec
          write(sformat(8:10),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(12:15),'(I4.4)') minpc*2
          do is=1,nspec
          do ipc=1,npc
            icur(1,ipc,is) = abs(gcount(2)%influx(ipc,is)*2/ustep*q(is)/renq/dt)
            icur(1,ipc,is) = icur(1,ipc,is)*sqdscaled(ipc)
            icur(3,ipc,is) = icur(3,ipc,is)*expema &
           &               + icur(1,ipc,is)*(1.0d0-expema)
            if(abs(icur(3,ipc,is)).ge.1.0d-9) then
              icur(2,ipc,is) = icur(3,ipc,is)
            else
              icur(2,ipc,is) = 0.0d0
            end if
          end do
          end do
          write(97,sformat) (istep, icur(1:2,1:minpc,is), is=1,minspec)
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(I,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') minpc*minspec
          write(90,sformat) istep, int(gcount(2)%outflux(1:minpc,1:minspec)*2/ustep)
!          if(istep.eq.nstep) &
!         &  write(90,sformat) istep*2, outfluxintgrtd(1:npc,1:nspec)
!         ----------- 
          sformat = '(xxxx(Ix.x,xxxx(1X,ES9.4e1),2X))'
          write(sformat(2:5),'(I4.4)') minspec
          write(sformat(8:10),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(12:15),'(I4.4)') minpc*2
          do is=1,nspec
          do ipc=1,npc
            ocur(1,ipc,is) = abs(gcount(2)%outflux(ipc,is)*2/ustep*q(is)/renq/dt)
            ocur(1,ipc,is) = ocur(1,ipc,is)*sqdscaled(ipc)
            ocur(3,ipc,is) = ocur(3,ipc,is)*expema &
           &               + ocur(1,ipc,is)*(1.0d0-expema)
            if(abs(ocur(3,ipc,is)).ge.1.0d-9) then
              ocur(2,ipc,is) = ocur(3,ipc,is)
            else
              ocur(2,ipc,is) = 0.0d0
            end if
          end do
          end do
          write(98,sformat) (istep, ocur(1:2,1:minpc,is), is=1,minspec)
!         -----------
          sformat = '(Ix.x,1X,xxxx(I,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') (1+npc)*minpc*minspec
          write(91,sformat) istep, int(gcount(2)%infhist(0:npc,1:minpc,1:minspec))
!          if(istep.eq.nstep) &
!         &  write(91,sformat) istep*2, infhistintgrtd(0:npc,1:npc,1:nspec)
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(I,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') nspec
          write(92,sformat) istep, int(gcount(2)%nesc(1:nspec)*2/ustep)
!          if(istep.eq.nstep) &
!         &  write(92,sformat) istep*2, nescintgrtd(0:npc,1:npc,1:nspec)
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') minpc
          write(93,sformat) istep, gcount(2)%chgacm(1,1:npc)/renq
!          if(istep.eq.nstep) &
!         &  write(93,sformat) istep*2, chgacmintgrtd(1,1:npc)
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') minpc
          write(94,sformat) istep, gcount(2)%chgacm(2,1:npc)/renq
!          if(istep.eq.nstep) &
!         &  write(94,sformat) istep*2, chgacmintgrtd(2.1:npc)
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') minpc
          write(95,sformat) istep, (rhoind(1:npc)/wrelax-rhoindh(1:npc)/wrelax-gcount(2)%chgacm(2,1:npc))/renq*rent*2.0d0/ustep
!          if(istep.eq.nstep) &
!         &  write(95,sformat) istep*2, rhoindintgrtd(1:npc)
!         ----------- 
          write(96,'(f7.2,1X,36(f9.2,1X))') &
         &  t, seygl(1:11,1),seygl(12,1)*4.0d0,seygl(1:11,3), &
         &     seygl(12,3)*4.0d0,seygl(1:11,4),seygl(12,4)*4.0d0
!         ----------- 
          sformat = '(Ix.x,xxxx(1X,xxxx(1X,ES9.4e1)))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(7:10),'(I4.4)') nspec + 1
          write(sformat(15:18),'(I4.4)') max(2,1)
          tew = t - dt*nretard
          write(99,sformat) istep,(Ew(1,is,1)*cos(omegaw(1)*tew)/rene, &
         &                         Ew(2,is,1)*sin(omegaw(1)*tew)/rene,is=0,nspec)
!         ----------- 
!          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
!          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
!          write(sformat(10:13),'(I4.4)') 10
!          write(99,sformat) istep, (selfp(1)+selfp(npc))*0.5d0/renphi, &

!         &               (biasc(1)%val+biasc(2)%val)*0.5d0/renq, &

!         &              +(sum(gcount(2)%outflux(1:2,nspec)) &
!         &               +sum(gcount(2)%outflux(npc-1:npc,nspec)) &
!         &               -sum(gcount(2)%infhist(1,1:npc,nspec)) &
!         &               -sum(gcount(2)%infhist(2,1:npc,nspec)) &
!         &               -sum(gcount(2)%infhist(npc-1,1:npc,nspec)) &
!         &               -sum(gcount(2)%infhist(npc,1:npc,nspec)))*q(3)*0.5d0/renq, &

!         &              -(sum(gcount(2)%infhist(0,1:2,1)) &
!         &               +sum(gcount(2)%infhist(0,npc-1:npc,1)))*q(1)*0.5d0/renq, &

!         &              -(sum(gcount(2)%infhist(0,1:2,2)) &
!         &               +sum(gcount(2)%infhist(0,npc-1:npc,2)))*q(2)*0.5d0/renq, &

!         &              +(sum(gcount(2)%infhist(1,3:npc,nspec)) &
!         &               +sum(gcount(2)%infhist(2,3:npc,nspec)) &
!         &               +sum(gcount(2)%infhist(npc-1,1:npc-2,nspec)) &
!         &               +sum(gcount(2)%infhist(npc,1:npc-2,nspec)))*q(3)*0.5d0/renq, &

!         &              -(sum(gcount(2)%infhist(3:npc,1,nspec)) &
!         &               +sum(gcount(2)%infhist(3:npc,2,nspec)) &
!         &               +sum(gcount(2)%infhist(1:npc-2,npc-1,nspec)) &
!         &               +sum(gcount(2)%infhist(1:npc-2,npc,nspec)))*q(3)*0.5d0/renq, &

!         &              -(   (gcount(2)%infhist(3,1,nspec)) &
!         &               +   (gcount(2)%infhist(3,2,nspec)) &
!         &               +   (gcount(2)%infhist(npc-2,npc-1,nspec)) &
!         &               +   (gcount(2)%infhist(npc-2,npc,nspec)))*q(3)*0.5d0/renq, &

!         &              -(   (gcount(2)%infhist(4,1,nspec)) &
!         &               +   (gcount(2)%infhist(4,2,nspec)) &
!         &               +   (gcount(2)%infhist(npc-3,npc-1,nspec)) &
!         &               +   (gcount(2)%infhist(npc-3,npc,nspec)))*q(3)*0.5d0/renq, &

!         &              -(sum(gcount(2)%infhist(5:npc-4,1,nspec)) &
!         &               +sum(gcount(2)%infhist(5:npc-4,2,nspec)) &
!         &               +sum(gcount(2)%infhist(5:npc-4,npc-1,nspec)) &
!         &               +sum(gcount(2)%infhist(5:npc-4,npc,nspec)))*q(3)*0.5d0/renq
!          if(istep.eq.nstep) &
!         &  write(99,sformat) istep*2, vidata(1:10)
          if(mod(istep,intfoc).eq.0.or.istep.eq.nstep) then
            close(79)
            close(86)
            close(89)
            close(90)
            close(91)
            close(92)
            close(93)
            close(94)
            close(95)
            close(96)
            close(97)
            close(98)
            close(99)
          end if
        end if
      end if
!
!-------------------- clear particle counter
      gcount(1)%influx = 0
      gcount(1)%outflux = 0
      gcount(1)%infhist = 0
      gcount(1)%nesc = 0
      gcount(1)%chgacm(2,:) = 0.0d0


  return
  end subroutine
