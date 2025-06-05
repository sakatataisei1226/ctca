#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine load_essnap
!
!   ____________________________________________________________
!
!            S U B R O U T I N E   L O A D _ E S S N A P
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine load ES-snap data of the previous job.  .
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
  integer(kind=4) :: m, mm, ns(ispec),ne(ispec), nee
  integer(kind=4) :: is, ierr, icon
  integer(kind=4) :: ipc
  integer(kind=4) :: rid
  integer(kind=4) :: prevnspec,prevxu,prevyu,prevzu,prevnpc
  integer(kind=4),allocatable :: prevnp(:)
  real(kind=8) :: vxtemp, vytemp, vztemp
  character(len=20) :: finpfn


!-------------------- 

!-------------------- read the ES-snap data from the file
      if(myid.eq.0) &
     &  write(6,*) 'load ES-snap data from ../FINAL/essnapsht'
      write(finpfn,'(a,i4.4)') '../FINAL/essnapsht',myid
      open(60,file=finpfn,status='old',form='unformatted')
      rewind 60
      read(60) prevnspec
!
!     --------------- 
      allocate(prevnp(prevnspec))
!
!     --------------- 
      read(60) prevnp(1:prevnspec),prevxu,prevyu,prevzu,prevnpc
!
!     --------------- 
      if(prevnspec.gt.nspec) then
        if(myid.eq.0) print*, &
       &  "Error: nspec should be equal to or greater than the previous one"
        if(myid.eq.0) &
       &  print*,"Error: stop..."
        stop
      end if
!
      if(prevxu.ne.sdoms(2,1,sdid(1)+1)-sdoms(1,1,sdid(1)+1).or. &
         prevyu.ne.sdoms(2,2,sdid(1)+1)-sdoms(1,2,sdid(1)+1).or. &
         prevzu.ne.sdoms(2,3,sdid(1)+1)-sdoms(1,3,sdid(1)+1)) then
        if(myid.eq.0) &
       &  print*,"Error: subdomain coords. do not match the loaded data"
        if(myid.eq.0) &
       &  print*,"Error: stop..."
        stop
      end if
!
      if(prevnpc.ne.npc) then
        if(myid.eq.0) &
       &  print*,"Error: npc does not match the loaded data"
        if(myid.eq.0) &
       &  print*,"Error: stop..."
        stop
      end if
!
!     --------------- 
      totalp(1:prevnspec,1) = prevnp(1:prevnspec)
      nee = 0
      do is=1,prevnspec
        ns(is) = nee + 1
        ne(is) = nee + prevnp(is)
        nee = ne(is)
      end do
!
!     --------------- 
      read(60) (pbuf(ns(is):ne(is)),is=1,prevnspec)
      read(60) rhobk(1,0:prevxu,0:prevyu,0:prevzu,3)
      read(60) gcount(1)%chgacm(1,1:prevnpc)
      close(60, status='keep')
      call MPI_Barrier(MCW,ierr)
      if(myid.eq.0) &
     &  write(6,*) 'ES-snap data loaded'


!-------------------- 
      do is=1,prevnspec
        do m=ns(is),ne(is)
          if(pbuf(m)%nid.ge.0) then
            nphgram(pbuf(m)%nid+1,is,1) = nphgram(pbuf(m)%nid+1,is,1) + 1
          end if
        end do
      end do


!-------------------- 
      nee = sum(prevnp(1:prevnspec))
ISL1: do is=prevnspec+1,nspec
        totalp(is,1) = npin(is)/nnode
        if(myid.lt.mod(npin(is),nnode)) then
          totalp(is,1) = totalp(is,1) + 1
        end if
        ns(is) = nee + 1
        ne(is) = nee + totalp(is,1)
!
        call RANU0(dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(DVRAU4): myid,icon=",myid,icon
        do m=ns(is),ne(is)
          pbuf(m)%vx = dranu(m-ns(is)+1)
        end do
        call RANU0(dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(DVRAU4): myid,icon=",myid,icon
        do m=ns(is),ne(is)
          pbuf(m)%vy = dranu(m-ns(is)+1)
        end do
        call RANU0(dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(DVRAU4): myid,icon=",myid,icon
        do m=ns(is),ne(is)
          pbuf(m)%vz = dranu(m-ns(is)+1)
        end do
!
        m = ns(is)
        do mm=0,totalp(is,1)-1
          pbuf(m)%x = slx*pbuf(ns(is)+mm)%vx
          pbuf(m)%y = sly*pbuf(ns(is)+mm)%vy
          pbuf(m)%z = slz*pbuf(ns(is)+mm)%vz
          rid = oh3_map_particle_to_subdomain &
         &        (pbuf(m)%x/dr,pbuf(m)%y/dr,pbuf(m)%z/dr)
          pbuf(m)%nid = rid
          do ipc=1,npc
            if(pbuf(m)%x.ge.xlpc(ipc)*dr.and.pbuf(m)%x.le.xupc(ipc)*dr.and. &
           &   pbuf(m)%y.ge.ylpc(ipc)*dr.and.pbuf(m)%y.le.yupc(ipc)*dr.and. &
           &   pbuf(m)%z.ge.zlpc(ipc)*dr.and.pbuf(m)%z.le.zupc(ipc)*dr) then
              pbuf(m)%x = 0.0d0
              pbuf(m)%y = 0.0d0
              pbuf(m)%z = 0.0d0
              pbuf(m)%nid = -1
            end if
          end do
          if(pbuf(m)%nid.ge.0) then
            nphgram(pbuf(m)%nid+1,is,1) = nphgram(pbuf(m)%nid+1,is,1) + 1
            m = m + 1
          end if
        end do
        totalp(is,1) = m - ns(is)
        ne(is) = m - 1
        nee = ne(is)
!
        call RANN0(spe(is)*dcos(speth(is)/180.0d0*pi),peth(is), &
       &            dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
        do m=ns(is),ne(is)
          pbuf(m)%vx = dranu(m-ns(is)+1)
        end do
        call RANN0(spe(is)*dsin(speth(is)/180.0d0*pi),peth(is), &
       &            dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
        do m=ns(is),ne(is)
          pbuf(m)%vy = dranu(m-ns(is)+1)
        end do
        call RANN0(spa(is),path(is),dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
        do m=ns(is),ne(is)
          pbuf(m)%vz = dranu(m-ns(is)+1)
        end do
!
        do m=ns(is),ne(is)
          vxtemp = pbuf(m)%vx*t11 + pbuf(m)%vy*t12 + pbuf(m)%vz*t13
          vytemp = pbuf(m)%vx*t21 + pbuf(m)%vy*t22 + pbuf(m)%vz*t23
          vztemp = pbuf(m)%vx*t31 + pbuf(m)%vy*t32 + pbuf(m)%vz*t33
          pbuf(m)%vx = vxtemp
          pbuf(m)%vy = vytemp
          pbuf(m)%vz = vztemp
        end do
!
        do m=ns(is),ne(is)
          pbuf(m)%spec = is
        end do
      end do ISL1


!-------------------- 
      deallocate(prevnp)


  return
  end subroutine



!
  subroutine save_esdat
!
!   ____________________________________________________________
!
!            S U B R O U T I N E   S A V E _ E S D A T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine save ES data for continuous job.        .
!   ............................................................
!

!-------------------- parameters and variables
  use oh_type
  use paramt
  use allcom
  use hdf
#define MCW local_comm
#define MSS MPI_STATUS_SIZE
  implicit none
!
  integer(kind=4) :: m, ns(ispec,2),ne(ispec,2), nee
  integer(kind=4) :: ps, is, ierr
  integer(kind=4) :: xu,yu,zu
!
  integer(kind=HID_T) :: fileid
  integer(kind=4) :: stats0,stats1
  integer(kind=8) :: dims(3)
  character(len=30) :: filename,dsname


!-------------------- 
      xu = sdoms(2,1,sdid(1)+1) - sdoms(1,1,sdid(1)+1)
      yu = sdoms(2,2,sdid(1)+1) - sdoms(1,2,sdid(1)+1)
      zu = sdoms(2,3,sdid(1)+1) - sdoms(1,3,sdid(1)+1)


!-------------------- 
      nee = 0
      do ps=1,2
      do is=1,nspec
        ns(is,ps) = nee + 1
        ne(is,ps) = nee + totalp(is,ps)
        nee = ne(is,ps)
      end do
      end do


!-------------------- rescale
!      do m=1,pbase(3)
!        pbuf(m)%x  = pbuf(m)%x/renr
!        pbuf(m)%y  = pbuf(m)%y/renr
!        pbuf(m)%z  = pbuf(m)%z/renr
!        pbuf(m)%vx = pbuf(m)%vx/renv
!        pbuf(m)%vy = pbuf(m)%vy/renv
!        pbuf(m)%vz = pbuf(m)%vz/renv
!      end do
!      rhobk(1,0:xu,0:yu,0:zu,3) = rhobk(1,0:xu,0:yu,0:zu,3)/renrho
!      gcount(1)%chgacm(1,1:npc) = gcount(1)%chgacm(1,1:npc)/renq


!-------------------- write the ES-snap data into a file
      if(myid.eq.0) &
     &  write(6,*) 'write ES-snapshot data into ./SNAPSHOT1'
      write(filename,'(a,i4.4,a)') './SNAPSHOT1/esdat', myid, '.h5'
      call hdfopen(filename,fileid,DFACC_CREATE)
!
      if(myid.eq.0) then
        dsname = 'nnode'
        dims(1) = 1
        call wrt1i(fileid,dsname,dims(1:1),(/nnode/),stats0,stats1)
      end if
!
      dsname = 'nspec'
      dims(1) = 1
      call wrt1i(fileid,dsname,dims(1:1),(/nspec/),stats0,stats1)
!
      dsname = 'np'
      dims(1) = nspec
      call wrt1i(fileid,dsname,dims(1:1), &
     &  (totalp(1:nspec,1)+totalp(1:nspec,2)),stats0,stats1)
!
      dsname = 'npc'
      dims(1) = 1
      call wrt1i(fileid,dsname,dims(1:1),(/npc/),stats0,stats1)
!
      dsname = 'sdoms'
      dims(1) = 2; dims(2) = 3
      call wrt2i(fileid,dsname,dims(1:2),sdoms(:,:,1),stats0,stats1)
!
      dsname = 'rens'
      dims(1) = 4
      call wrt1d(fileid,dsname,dims(1:1), &
     &  (/renr, renv, renrho, renq/),stats0,stats1)
!
      dsname = 'px'
      dims(1) = pbase(3)
      call wrt1d(fileid,dsname,dims(1:1), &
     &  pbuf(1:pbase(3))%x,stats0,stats1)
!
      dsname = 'py'
      dims(1) = pbase(3)
      call wrt1d(fileid,dsname,dims(1:1), &
     &  pbuf(1:pbase(3))%y,stats0,stats1)
!
      dsname = 'pz'
      dims(1) = pbase(3)
      call wrt1d(fileid,dsname,dims(1:1), &
     &  pbuf(1:pbase(3))%z,stats0,stats1)
!
      dsname = 'pvx'
      dims(1) = pbase(3)
      call wrt1d(fileid,dsname,dims(1:1), &
     &  pbuf(1:pbase(3))%vx,stats0,stats1)
!
      dsname = 'pvy'
      dims(1) = pbase(3)
      call wrt1d(fileid,dsname,dims(1:1), &
     &  pbuf(1:pbase(3))%vy,stats0,stats1)
!
      dsname = 'pvz'
      dims(1) = pbase(3)
      call wrt1d(fileid,dsname,dims(1:1), &
     &  pbuf(1:pbase(3))%vz,stats0,stats1)
!
      dsname = 'pres'
      dims(1) = pbase(3)
      call wrt1i(fileid,dsname,dims(1:1), &
     &  pbuf(1:pbase(3))%preside,stats0,stats1)
!
      dsname = 'pis'
      dims(1) = pbase(3)
      call wrt1i(fileid,dsname,dims(1:1), &
     &  pbuf(1:pbase(3))%spec,stats0,stats1)
!
      dsname = 'rhobk'
      dims(1) = xu + 1; dims(2) = yu + 1; dims(3) = zu + 1
      call wrt3d(fileid,dsname,dims(1:3), &
     &  rhobk(1,0:xu,0:yu,0:zu,3),stats0,stats1)
!
      dsname = 'chgacm'
      dims(1) = npc
      call wrt1d(fileid,dsname,dims(1:1), &
     &  gcount(1)%chgacm(1,1:npc),stats0,stats1)
!
      call hdfclose(fileid,stats0)
      call MPI_Barrier(MCW,ierr)
      if(myid.eq.0) &
     &  write(6,*) 'ESDAT written'


  return
  end subroutine



!
  subroutine load_emsnap
!
!   ____________________________________________________________
!
!            S U B R O U T I N E   L O A D _ E M S N A P
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine load EM-snap data of the previous job.  .
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
  integer(kind=4) :: m, mm, ns(ispec),ne(ispec), nee
  integer(kind=4) :: is, ierr, icon
  integer(kind=4) :: ipc
  integer(kind=4) :: rid
  integer(kind=4) :: prevnspec,prevxu,prevyu,prevzu,prevnpc
  integer(kind=4),allocatable :: prevnp(:)
  real(kind=8) :: vxtemp, vytemp, vztemp
  character(len=20) :: finpfn


!-------------------- 

!-------------------- read the ES-snap data from the file
      if(myid.eq.0) &
     &  write(6,*) 'load EM-snap data from ../FINAL/emsnapsht'
      write(finpfn,'(a,i4.4)') '../FINAL/emsnapsht',myid
      open(60,file=finpfn,status='old',form='unformatted')
      rewind 60
      read(60) prevnspec
!
!     --------------- 
      allocate(prevnp(prevnspec))
!
!     --------------- 
      read(60) prevnp(1:prevnspec),prevxu,prevyu,prevzu,prevnpc
!
!     --------------- 
      if(prevnspec.gt.nspec) then
        if(myid.eq.0) print*, &
       &  "Error: nspec should be equal to or greater than the previous one"
        if(myid.eq.0) &
       &  print*,"Error: stop..."
        stop
      end if
!
      if(prevxu.ne.sdoms(2,1,sdid(1)+1)-sdoms(1,1,sdid(1)+1).or. &
         prevyu.ne.sdoms(2,2,sdid(1)+1)-sdoms(1,2,sdid(1)+1).or. &
         prevzu.ne.sdoms(2,3,sdid(1)+1)-sdoms(1,3,sdid(1)+1)) then
        if(myid.eq.0) &
       &  print*,"Error: subdomain coords. do not match the loaded data"
        if(myid.eq.0) &
       &  print*,"Error: stop..."
        stop
      end if
!
      if(prevnpc.ne.npc) then
        if(myid.eq.0) &
       &  print*,"Error: npc does not match the loaded data"
        if(myid.eq.0) &
       &  print*,"Error: stop..."
        stop
      end if
!
!     --------------- 
      totalp(1:prevnspec,1) = prevnp(1:prevnspec)
      nee = 0
      do is=1,prevnspec
        ns(is) = nee + 1
        ne(is) = nee + prevnp(is)
        nee = ne(is)
      end do
!
!     --------------- 
      read(60) (pbuf(ns(is):ne(is)),is=1,prevnspec)
      read(60) eb(EX:BZ,-1:prevxu+1,-1:prevyu+1,-1:prevzu+1,1)
      read(60) rho(1,0:prevxu,0:prevyu,0:prevzu,1)
      read(60) rhobk(1,0:prevxu,0:prevyu,0:prevzu,3)
      read(60) gcount(1)%chgacm(1,1:prevnpc)
      close(60, status='keep')
      call MPI_Barrier(MCW,ierr)
      if(myid.eq.0) &
     &  write(6,*) 'EM-snap data loaded'


!-------------------- 
      do is=1,prevnspec
        do m=ns(is),ne(is)
          if(pbuf(m)%nid.ge.0) then
            nphgram(pbuf(m)%nid+1,is,1) = nphgram(pbuf(m)%nid+1,is,1) + 1
          end if
        end do
      end do


!-------------------- 
      nee = sum(prevnp(1:prevnspec))
ISL1: do is=prevnspec+1,nspec
        totalp(is,1) = npin(is)/nnode
        if(myid.lt.mod(npin(is),nnode)) then
          totalp(is,1) = totalp(is,1) + 1
        end if
        ns(is) = nee + 1
        ne(is) = nee + totalp(is,1)
!
        call RANU0(dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
        do m=ns(is),ne(is)
          pbuf(m)%vx = dranu(m-ns(is)+1)
        end do
        call RANU0(dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
        do m=ns(is),ne(is)
          pbuf(m)%vy = dranu(m-ns(is)+1)
        end do
        call RANU0(dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
        do m=ns(is),ne(is)
          pbuf(m)%vz = dranu(m-ns(is)+1)
        end do
!
        m = ns(is)
        do mm=0,totalp(is,1)-1
          pbuf(m)%x = slx*pbuf(ns(is)+mm)%vx
          pbuf(m)%y = sly*pbuf(ns(is)+mm)%vy
          pbuf(m)%z = slz*pbuf(ns(is)+mm)%vz
          rid = oh3_map_particle_to_subdomain &
         &        (pbuf(m)%x/dr,pbuf(m)%y/dr,pbuf(m)%z/dr)
          pbuf(m)%nid = rid
          do ipc=1,npc
            if(pbuf(m)%x.ge.xlpc(ipc)*dr.and.pbuf(m)%x.le.xupc(ipc)*dr.and. &
           &   pbuf(m)%y.ge.ylpc(ipc)*dr.and.pbuf(m)%y.le.yupc(ipc)*dr.and. &
           &   pbuf(m)%z.ge.zlpc(ipc)*dr.and.pbuf(m)%z.le.zupc(ipc)*dr) then
              pbuf(m)%x = 0.0d0
              pbuf(m)%y = 0.0d0
              pbuf(m)%z = 0.0d0
              pbuf(m)%nid = -1
            end if
          end do
          if(pbuf(m)%nid.ge.0) then
            nphgram(pbuf(m)%nid+1,is,1) = nphgram(pbuf(m)%nid+1,is,1) + 1
            m = m + 1
          end if
        end do
        totalp(is,1) = m - ns(is)
        ne(is) = m - 1
        nee = ne(is)
!
        call RANN0(spe(is)*dcos(speth(is)/180.0d0*pi),peth(is), &
       &            dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
        do m=ns(is),ne(is)
          pbuf(m)%vx = dranu(m-ns(is)+1)
        end do
        call RANN0(spe(is)*dsin(speth(is)/180.0d0*pi),peth(is), &
       &            dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
        do m=ns(is),ne(is)
          pbuf(m)%vy = dranu(m-ns(is)+1)
        end do
        call RANN0(spa(is),path(is),dranu,totalp(is,1),icon)

        if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
        do m=ns(is),ne(is)
          pbuf(m)%vz = dranu(m-ns(is)+1)
        end do
!
        do m=ns(is),ne(is)
          vxtemp = pbuf(m)%vx*t11 + pbuf(m)%vy*t12 + pbuf(m)%vz*t13
          vytemp = pbuf(m)%vx*t21 + pbuf(m)%vy*t22 + pbuf(m)%vz*t23
          vztemp = pbuf(m)%vx*t31 + pbuf(m)%vy*t32 + pbuf(m)%vz*t33
          pbuf(m)%vx = vxtemp
          pbuf(m)%vy = vytemp
          pbuf(m)%vz = vztemp
        end do
!
        do m=ns(is),ne(is)
          pbuf(m)%spec = is
        end do
      end do ISL1


!-------------------- 
      deallocate(prevnp)


  return
  end subroutine



!
  subroutine save_emsnap
!
!   ____________________________________________________________
!
!            S U B R O U T I N E   S A V E _ E M S N A P
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine save EM-snap data for continuous job.   .
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
  integer(kind=4) :: m, ns(ispec,2),ne(ispec,2), nee
  integer(kind=4) :: ps, is, ierr
  integer(kind=4) :: xu,yu,zu
  character(len=20) :: finpfn


!-------------------- 
      xu = sdoms(2,1,sdid(1)+1) - sdoms(1,1,sdid(1)+1)
      yu = sdoms(2,2,sdid(1)+1) - sdoms(1,2,sdid(1)+1)
      zu = sdoms(2,3,sdid(1)+1) - sdoms(1,3,sdid(1)+1)


!-------------------- 
      nee = 0
      do ps=1,2
      do is=1,nspec
        ns(is,ps) = nee + 1
        ne(is,ps) = nee + totalp(is,ps)
        nee = ne(is,ps)
      end do
      end do


!-------------------- rescale
      do m=1,pbase(3)
        pbuf(m)%x  = pbuf(m)%x/renr
        pbuf(m)%y  = pbuf(m)%y/renr
        pbuf(m)%z  = pbuf(m)%z/renr
        pbuf(m)%vx = pbuf(m)%vx/renv
        pbuf(m)%vy = pbuf(m)%vy/renv
        pbuf(m)%vz = pbuf(m)%vz/renv
      end do
      eb(EX:EZ,-1:xu+1,-1:yu+1,-1:zu+1,1) = &
     &  eb(EX:EZ,-1:xu+1,-1:yu+1,-1:zu+1,1)/rene
      eb(BX:BZ,-1:xu+1,-1:yu+1,-1:zu+1,1) = &
     &  eb(BX:BZ,-1:xu+1,-1:yu+1,-1:zu+1,1)/renb
      rho(1,0:xu,0:yu,0:zu,1) = rho(1,0:xu,0:yu,0:zu,1)/renrho
      rhobk(1,0:xu,0:yu,0:zu,3) = rhobk(1,0:xu,0:yu,0:zu,3)/renrho
      gcount(1)%chgacm(1,1:npc) = gcount(1)%chgacm(1,1:npc)/renq


!-------------------- write the ES-snap data into a file
      if(myid.eq.0) &
     &  write(6,*) 'write EM-snap data into FINAL/emsnapsht'
      write(finpfn,'(a,i4.4)') './FINAL/emsnapsht',myid
      open(61,file=finpfn,status='unknown',form='unformatted')
      rewind 61
      write(61) nspec
      write(61) totalp(1:nspec,1)+totalp(1:nspec,2),xu,yu,zu,npc
      write(61) ((pbuf(ns(is,ps):ne(is,ps)),ps=1,2),is=1,nspec)
      write(61) eb(EX:BZ,-1:xu+1,-1:yu+1,-1:zu+1,1)
      write(61) rho(1,0:xu,0:yu,0:zu,1)
      write(61) rhobk(1,0:xu,0:yu,0:zu,3)
      write(61) gcount(1)%chgacm(1,1:npc)
      close(61, status='keep')
      call MPI_Barrier(MCW,ierr)
      if(myid.eq.0) &
     &  write(6,*) 'EM-snap data written'


  return
  end subroutine
