#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  module allcom
!
!   ____________________________________________________________
!
!                    M O D U L E   A L L C O M
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   module for common declaration of variables and arrays  .
!   ............................................................
!
!-------------------- declaration of global variables and arrays
  use mpi
  use oh_type
  use paramt
  use lemses
  use ctca
#define MCW local_comm
  implicit none
!
!   ----------------- /esorem/
    real*8                :: deltax,deltat,massratio,density,etemp,itemp,phtemp
    real*8                :: flowvx,flowvy,flowvz,magb,jph,cvsim,weightph

!   ----------------- /esorem/
    integer               :: emflag
    integer               :: nflag_testp

!   ----------------- /ohhelp3/
    integer               :: sdid(2)
    integer,allocatable   :: nphgram(:,:,:)
    integer,allocatable   :: nphgram2(:,:,:)
    integer,allocatable   :: totalp(:,:)
    integer,allocatable   :: famind(:)
    integer,allocatable   :: fammbr(:)
    type(oh_particle),allocatable &
   &                      :: pbuf(:)
    type(oh_particle)     :: pinj(2)
    integer               :: pbase(3)
    integer               :: maxlocalp
    type(oh_mycomm)       :: mycomm
    integer               :: nbor(3,3,3)
    integer               :: nborps(27,ispec,2)
    integer               :: pcoord(OH_DIMENSION)
    integer,allocatable   :: sdoms(:,:,:)
    integer,allocatable   :: slcs(:,:)
    integer,allocatable   :: scnts(:),rcnts(:)
    integer               :: scoord(2,OH_DIMENSION)
    integer               :: nbound
    integer               :: bcond(2,OH_DIMENSION)
    integer,allocatable   :: bounds(:,:,:)
    integer               :: ftypes(7,FNR+1)
    integer               :: cfields(CNR+1)
    integer               :: ctypes(3,2,BNR+1,CNR)
    integer               :: fsizes(2,OH_DIMENSION,FNR)
    real*8,allocatable    :: eb(:,:,:,:,:)
    real*8,allocatable    :: ebav(:,:,:,:)
    real*8,allocatable    :: mp(:,:,:,:,:)
    real*8,allocatable    :: db(:,:,:,:,:)
    real*8,allocatable    :: dbs(:,:,:,:,:)
    real*8,allocatable    :: aj(:,:,:,:,:)
    real*8,allocatable    :: ajdg(:,:,:,:,:)
    real*8,allocatable    :: ajav(:,:,:,:)
    real*8,allocatable    :: rho(:,:,:,:,:)
    real*8,allocatable    :: rho_tilde(:,:,:,:,:)  ! rho + boundary conditions
    real*8,allocatable    :: phi(:,:,:,:,:)
    real*8,allocatable    :: phiav(:,:,:,:)
    real*8,allocatable    :: wrk(:,:,:,:)
    real*8,allocatable    :: rhobk(:,:,:,:,:)
    real*8,allocatable    :: rhobksp(:,:,:,:,:,:)
    real*8,allocatable    :: rhodg(:,:,:,:,:)
    real*8,allocatable    :: rhoav(:,:,:,:)
    real*8,allocatable    :: colf(:,:,:,:,:)
    integer               :: currmode
    integer               :: istep
    integer               :: nxsd,nysd,nzsd
    integer               :: nxrho,nyrho,nzrho
    integer               :: nxphi,nyphi,nzphi
    integer               :: idxsd,idysd,idzsd
    integer               :: lxsd(FNR),lysd(FNR),lzsd(FNR)

!   ----------------- /ohhelpf/
    integer(kind=4),external :: oh2_max_local_particles
    integer(kind=4),external :: oh3_transbound
    integer(kind=4),external :: oh3_map_particle_to_neighbor
    integer(kind=4),external :: oh3_map_particle_to_subdomain

    double precision :: conductivity

!   ----------------- /constc/
    real(kind=8)    :: pi,pi2,pih, cv,dr,dt,dtf,dri,si,tcs,cs,phixy,phiz
    real(kind=8)    :: anx,any,anz, dkx,dky,dkz, slx,sly,slz
    real(kind=8)    :: b0,b0x,b0y,b0z, wc, alpha
    real(kind=8)    :: e0,e0x,e0y,e0z
    real(kind=8)    :: rho0
    real(kind=8)    :: gfactor, cfactor(3)

!   ----------------- /consar/
    real(kind=8)    :: wp(ispec), bor(ispec)
    real(kind=8)    :: rm(ispec), qm(ispec), q(ispec), qmr(ispec), qp(ispec)
    real(kind=8)    :: qpr(ispec),qmrs(ispec), wpmax(ispec)
    real(kind=8)    :: denmod(ispec),denk(ispec),omegalw(ispec)

!   ----------------- /consit/
    integer(kind=8) :: npsum
    integer(kind=4) :: nx,ny,nz
    integer(kind=4) :: nxl,nxr,nyl,nyr,nzl,nzr
    integer(kind=4) :: nxl0,nxr0,nyl0,nyr0,nzl0,nzr0
    real(kind=8)    :: xlcol(2),xucol(2),ylcol(2),yucol(2),zlcol(2),zucol(2)
    real(kind=8)    :: mfpath(ispec)
    real(kind=8)    :: vxyzmax(ispec)

!   ----------------- /consia/
    integer(kind=8) :: np(ispec)
    integer(kind=4) :: nfbnd(3), imask(7), npbnd(3,ispec)

!   ----------------- /cpmxcm/
    real(kind=8)    :: cpsum2,sfsum,phic,dpsum

!   ----------------- /cpmxar/
!    real(kind=8)    :: cgacm(inpc), cgimpg(inpc)
!    real(kind=8)    :: cgacmt(inpc), cgimpg(inpc)
    real(kind=8)    :: selfp(inpc)
    real(kind=8)    :: groupp(inpc), groupph(inpc), grouppl(inpc)
    real(kind=8)    :: phiref(2)
    real(kind=8)    :: phipr(20)
    real(kind=8)    :: biasp(inpc)
    real(kind=8)    :: dscaled(inpc), sqdscaled(inpc)
    type(connect)   :: biasc(inpc*(inpc-1)/2)
    integer(kind=4) :: bcfrom(inpc),bcto(inpc)
    real(kind=8)    :: bcval(inpc)
    real(kind=8)    :: apmx(inpc,inpc), bpmx(inpc,inpc)
    real(kind=8)    :: bcmx(inpc,inpc)
    real(kind=8)    :: diff(ixw)
    real(kind=8)    :: sfrho(ixw), sfrhoh(ixw)
    real(kind=8)    :: rhoind(inpc),rhoindh(inpc)
    real(kind=8),allocatable :: cpmx(:,:)
    integer(kind=4) :: npcg, prefcrd(3)
    integer(kind=4) :: pcgs(inpc+1), ccgs(inpc+1)
    integer(kind=4) :: ngref(2)

!   ----------------- /cpmxit/
    integer(kind=4) :: bdyvoxel(3,ixw)
    integer(kind=4),allocatable :: bdyedges(:,:,:,:)
    integer(kind=4) :: nbedge(3,inpc)
    real(kind=8)    :: bdygrid(3,ixw)

!   ----------------- /ecorcm/
    integer(kind=4) :: pftmode
    integer(kind=4) :: nxfft, nyfft, nzfft
    integer(kind=4) :: lslice(2),uslice(2), l2slice(2),u2slice(2)
    integer(kind=4) :: nwfft, ndfft, nhfft, lwfft, ldfft, lhfft
    integer(kind=4) :: nwdbs, nddbs, nhdbs, lwdbs, lddbs, lhdbs
    integer(kind=4) :: nxslc,nyslc,nzslc
    integer(kind=4) :: nwslc,ndslc,nhslc, lwslc,ldslc,lhslc
    integer(kind=4) :: mptype_xy(2), mptype_yz(2), mptype_xz(2)
    integer(kind=4) :: mptype_rs(0:1,0:1,2,4),mptype_sc(0:1,0:1,2,4)
    integer(kind=4) :: mptype_rc(-1:0,2,4),mptype_ss(-1:0,2,4)
    integer(kind=4) :: mptype_aa(0:1,0:1,2,4)
    integer(kind=4),allocatable :: xtype(:),ytype(:),ztype(:)
    integer(kind=kind(MPI_COMM_WORLD)) :: subcomm
    integer(kind=4) :: sxl,sxu, syl,syu, szl,szu, stxl,stxu
    real(kind=8),allocatable :: poi(:,:,:,:,:)
    real(kind=8),allocatable :: warray1d(:,:),warray2d(:,:,:)
    real(kind=8)    :: rnx,rny,rnz
    real(kind=8)    :: mfactor(2)
    integer(kind=8) :: fftplan(4,2,2)
    integer(kind=kind(FFTW_R2HC)) :: fftw_type(4,3,2)
    integer(kind=4) :: lxfft(4),uxfft(4),lyfft(4),uyfft(4),lzfft(4),uzfft(4)

!   ----------------- /mediac/
    integer(kind=4) :: idex(ix,iy,iz),idey(ix,iy,iz),idez(ix,iy,iz)
    integer(kind=4) :: idbx(ix,iy,iz),idby(ix,iy,iz),idbz(ix,iy,iz)

!   ----------------- /initp/
    integer(kind=4) :: ndst(ispec),nphi(ispec),ioptd(ispec)
    real(kind=8)    :: path(ispec), spa(ispec),spe(ispec),speth(ispec)
    real(kind=8)    :: peth(ispec), f(ispec)
    real(kind=8)    :: vdx(ispec),vdy(ispec),vdz(ispec)
    real(kind=8)    :: vdtx(ispec),vdty(ispec),vdtz(ispec)
!    real(kind=8)    :: xe0(ispec),ye0(ispec),ze0(ispec)
!    real(kind=8)    :: xed(ispec),yed(ispec),zed(ispec)
    real(kind=8)    :: vpa(ispec),vpb(ispec),vpe(ispec)
    real(kind=8)    :: vdri(ispec),vdthz(ispec),vdthxy(ispec)
    real(kind=8)    :: lcgamma(ispec),lcbeta(ispec)

!   ----------------- /inpcom/
    integer(kind=4) :: i1hdf,i2hdf,i3hdf
    integer(kind=4) :: isort(ispec),idim_sp
    integer(kind=4) :: intfoc
    integer(kind=4) :: hdfdigstart, hdfdigend
    integer(kind=4) :: ifdiag,isdiag,nstep,jobnum(2),nspec,minspec
    integer(kind=4) :: iediag,ijdiag,iddiag,imdig(ispec),iadiag
    integer(kind=4) :: ipadig(ispec),ipahdf(ispec),ipaxyz(6,ispec)
    integer(kind=4) :: daverg
    integer(kind=4) :: itchck, jxfltr,jyfltr,jzfltr,ipexfl,ipeyfl,ipezfl
    integer(kind=4) :: mltstp,mltstpf,juncan,maxt,ifxyz(7),ijxyz(4),irhsp(ispec)
    integer(kind=4) :: ildig(ispec),nfltrx,nfltry,nfltrz
    integer(kind=4) :: ikdiag,ikxmax,ikymax,ikzmax,nkx,nky,nkz
    integer(kind=4) :: ivdiag,nvx,nvy,nvz,ionchg
    integer(kind=4) :: lfdiag,ljdiag,ladiag,lpdiag(ispec)
    integer(kind=4) :: nfsnap,njsnap,npsnap(ispec)

!   ----------------- /potcon/
    real(kind=8),allocatable    :: dkxm(:),dkyn(:),dkzl(:)
    real(kind=8),allocatable    :: kmod(:,:,:)

!   ----------------- /rescom/
    real(kind=8),target :: vxf(inr), vyf(inr), vzf(inr)
    real(kind=8),target :: vxb(inr), vyb(inr), vzb(inr)
    real(kind=8),target :: vxr(inr), vyr(inr), vzr(inr)
    real(kind=8),pointer :: vnormf(:), vtangf(:), runiform(:)
    real(kind=8),pointer :: vnormc(:), vtangc(:), psicos(:)
    real(kind=8),pointer :: vnorms(:), vtangs1(:), vtangs2(:)

!   ----------------- /rsclcm/
    real(kind=8)    :: renr,rent,renv,renq,rene,renb
    real(kind=8)    :: renrho,renj,renm,renphi,rened
    real(kind=8)    :: rena
    real(kind=8)    :: reni, renreg, reng, renec

!   ----------------- /svpcom/
    integer(kind=8) :: npin(ispec)
    integer(kind=8) :: npr(ispec), nprsum
    integer(kind=4) :: inpf(ispec),inpb(ispec)
    integer(kind=4) :: injct(ispec)

!   ----------------- /timecm/
    integer(kind=4) :: lstep,itime
    real(kind=8)    :: t
    real(kind=8)    :: elatime

!   ----------------- /trnscm/
    real(kind=8)    :: t11,t12,t13,t21,t22,t23,t31,t32,t33

!   ----------------- /ptcond/
    integer(kind=4) :: npc,nxpc1(inpc),nypc1(inpc),nzpc1(inpc)
    integer(kind=4) :: nxpc2(inpc),nypc2(inpc),nzpc2(inpc)
    integer(kind=4) :: ncond(inpc),mtd_vchg(inpc),reducecm(inpc)
    integer(kind=4) :: nscpmx(inpc),necpmx(inpc)
    integer(kind=4) :: nmxcpmx(inpc),nmycpmx(inpc),nmzcpmx(inpc)
    integer(kind=4) :: nmxycpmx(inpc),nmyzcpmx(inpc),nmzxcpmx(inpc)
    integer(kind=4) :: nbdsf1(inpc)
    integer(kind=4) :: ncpmx,nbdsf,modeww,ncpcnt(inpc)
    integer(kind=4) :: geotype(inpc)
    integer(kind=4) :: nflag_subcell(2)
    integer(kind=4) :: pswper,pswstr
    integer(kind=4),allocatable :: ncpmxs(:), scpmxs(:)
    integer(kind=4),allocatable :: ndcpmx(:), ndcpmx2(:)
    integer(kind=4),allocatable :: nsdcpmx(:), nsdcpmx2(:), nedcpmx(:)
    integer(kind=4),allocatable :: nbdglocal(:),dispbdg(:)
    real(kind=8)    :: xlpc(inpc),ylpc(inpc),zlpc(inpc)
    real(kind=8)    :: xupc(inpc),yupc(inpc),zupc(inpc)
    real(kind=8)    :: pfixed(inpc),v_max(inpc),v_omega(inpc),sfecrrct,wrelax
    real(kind=8)    :: pswspn(inpc),pswini
    real(kind=8)    :: epc2,amu2,sigma,oradius(3,inpc)
    real(kind=8)    :: xc,yc,zc,rc
    real(kind=8)    :: mingap
    type(cylindrical) :: cylinder(inpc)
    type(spherical) :: sphere(inpc)
    type(wire)      :: boom(inpc)
    integer(kind=4) :: bdyalign(inpc)
    real(kind=8)    :: bdyradius(inpc), bdyedge(2,inpc), bdycoord(3,inpc)
    integer(kind=4) :: wirealign(inpc)
    real(kind=8) :: wirehlength(inpc), wirerradius(inpc), wireeradius(inpc)
    real(kind=8) :: wireorigin(3,inpc)

!   ----------------- /solidsurf/
    real(kind=8)    :: zssurf, zsbuf
    real(kind=8)    :: xbowlc,ybowlc,zbowlc,rbowl
    real(kind=8)    :: xdomec,ydomec,zdomec,rdome
    real(kind=8)    :: dray(3), laray, dcdome,dcbowl
    real(kind=8)    :: xholec,yholec,rhole,rholesq,rholeinv,dhole,dlhole
    real(kind=8)    :: lbhole,ubhole,flhole,llbhole,lubhole
    real(kind=8)    :: xlrechole(2),ylrechole(2),zlrechole(2)
    real(kind=8)    :: xurechole(2),yurechole(2),zurechole(2)

    integer, parameter :: nboundary_types = 10
    ! Boundary type
    ! "flat-surface"
    ! "rectangle-hole"
    ! "cylinder-hole"
    ! "hyperboloid-hole"
    ! "rectangle"
    ! "circle"
    ! "cuboid"
    ! "disk"
    character(len=30) :: boundary_type = "flat-surface"
    character(len=30) :: boundary_types(nboundary_types) = "none"
    double precision :: cylinder_origin(nboundary_types, 3) = 0.0d0
    double precision :: cylinder_radius(nboundary_types) = 0.0d0
    double precision :: cylinder_height(nboundary_types) = 0.0d0
    !> Minimum radius relative to the entrance radius of the hyperbolic hole (r_min == rcurv*r_max).
    !> If rcurv == 1.0d0, then r_min == r_max.
    !> Else if rcurv == 0.5d0 then r_min == 0.5d0 * r_max
    double precision :: rcurv = 1.0d0
    double precision :: rectangle_shape(nboundary_types, 6) = 0.0d0
    double precision :: sphere_origin(nboundary_types, 3) = 0.0d0
    double precision :: sphere_radius(nboundary_types) = 0.0d0
    double precision :: circle_origin(nboundary_types, 3) = 0.0d0
    double precision :: circle_radius(nboundary_types) = 0.0d0
    double precision :: cuboid_shape(nboundary_types, 6) = 0.0d0
    double precision :: disk_origin(nboundary_types, 3) = 0.0d0
    double precision :: disk_height(nboundary_types) = 0.0d0
    double precision :: disk_radius(nboundary_types) = 0.0d0
    double precision :: disk_inner_radius(nboundary_types) = 0.0d0
    double precision :: plane_with_circle_hole_zlower(nboundary_types) = 0.0d0
    double precision :: plane_with_circle_hole_height(nboundary_types) = 0.0d0
    double precision :: plane_with_circle_hole_radius(nboundary_types) = 0.0d0

    logical :: pe_ray_cast = .false.  ! PE emission is determined using raycasting.
    

!   ----------------- /ecrctc/
    integer(kind=4) :: nflag_ecrct, mtd_vbnd(3), iphiref(2,3)
    double precision :: boundary_conditions(2, 3)

!   ----------------- /energy/
    real(kind=8)    :: rkt1(ispec),rkt2(ispec),rkt3(ispec)
    real(kind=8)    :: rkd1(ispec),rkd2(ispec),rkd3(ispec)
    real(kind=8)    :: eneary(10)
    real(kind=8)    :: engk(3,ispec), engeb(2), engkg(3,ispec), engebg(2)
    integer(kind=8) :: nactv(ispec)

!   ----------------- /workcm/
    real(kind=8),allocatable :: rfttab(:,:),csttab(:,:),snttab(:,:)

!   ----------------- /vicalf/
    real(kind=8)    :: surfj(ixy),surfe(8,ixy),phi1d(8,ixy)

!   ----------------- /jsrc/
    integer(kind=4) :: njs
    integer(kind=4) :: rjs(3,injs)
    real(kind=8)    :: ajs(3,injs),wjs(injs),th0js(injs)

!   ----------------- /wave/
    integer(kind=4) :: nwave, twave(inwave)
    integer(kind=4) :: imode
    real(kind=8)    :: ebamp(6,inwave),ldwv(inwave)
    real(kind=8)    :: orgwv(3,inwave),angwv(2,inwave)
    real(kind=8)    :: lambdatw, omegatw, amptw, thetatw

!   ----------------- /testch/
    integer(kind=4) :: ntch, nqclst
    integer(kind=4) :: dimtch(3,intch),dimclst(3,intch)
    real(kind=8)    :: rtch(3,intch)
    real(kind=8)    :: qtch(intch)
    real(kind=8)    :: e1tch(intch), eMtch(intch)
    real(kind=8)    :: p1tch(intch)
    real(kind=8)    :: rcutoff(intch), r2cutoff(intch)
    real(kind=8)    :: rclst(5,intch)
    real(kind=8)    :: rhopeak(intch)

!   ----------------- /dpolec/
    integer(kind=4) :: mode_dipole,line_mode,n_shth_rgn,ifeedst,ifeeded
    integer(kind=4) :: ngap,i_gap(ingap),j_gap(ingap),k_gap(ingap)
    integer(kind=4) :: ewmodel,nretard
    real(kind=8)    :: gap_amp(ingap),ajs_amp(ingap),f_c,w_c
    real(kind=8)    :: annrms(ingap), pgamp(ingap)
    real(kind=8)    :: tfeedst, offez(ingap), gap_jsrc(ingap), nstp_oe
    real(kind=8)    :: resc(ingap),capc(ingap),ezgap(ingap)
    real(kind=8)    :: xbpcc1,xbpcc2,ybpcc1,ybpcc2
    real(kind=8)    :: zbpcc1,zbpcc2,zbpcc3,zbpcc4
    real(kind=8)    :: zbpcc5,zbpcc6,zbpcc7,zbpcc8
    real(kind=8)    :: xline,yline
    real(kind=8)    :: Ew0,Ew(2,0:ispec,10),omegaw(10)

!   ----------------- /emissn/
    integer(kind=4) :: nflag_emit(ispec),nepl(ispec),ipcpl(inepl),nemd(inepl)
    integer(kind=4) :: omniemit(inepl)
!    integer(kind=4) :: nh1s(ispec),nh1e(ispec),nh2s(ispec),nh2e(ispec)
    integer(kind=8) :: narrxf(ispec),narryf(ispec),narrzf(ispec)
    integer(kind=8) :: narrxb(ispec),narryb(ispec),narrzb(ispec)
    real(kind=8)    :: arrxf(ispec),arryf(ispec),arrzf(ispec)
    real(kind=8)    :: arrxb(ispec),arryb(ispec),arrzb(ispec)
    real(kind=8)    :: flpf(ispec),flpb(ispec),dnsf(ispec),dnsb(ispec)
    real(kind=8)    :: curf(ispec),curb(ispec)
    real(kind=8)    :: flpfs(inepl),flpbs(inepl),dnsfs(inepl),dnsbs(inepl)
    real(kind=8)    :: curfs(inepl),curbs(inepl)
    real(kind=8)    :: abvdem(ispec),fluxf(inepl),fluxb(inepl)
    real(kind=8)    :: fluxub, fluxlb
    real(kind=8)    :: fluxv(ispec,inpc),imarea(inpc)
    real(kind=8)    :: remf(inepl),remb(inepl),xmaxe(inepl),xmine(inepl)
    real(kind=8)    :: ymaxe(inepl),ymine(inepl),zmaxe(inepl),zmine(inepl)
    real(kind=8)    :: xwdte(inepl),ywdte(inepl),zwdte(inepl)
    real(kind=8)    :: thetaz, thetaxy, thetax, thetay
    real(kind=8)    :: countertxy, ltxy,utxy
    real(kind=8)    :: psixy, psizx, psizy
    real(kind=8)    :: plreloc
    real(kind=8)    :: mainprop
    type pejection
      real(kind=8) :: grd, xl,yl,zl, xu,yu,zu, xlu,ylu,zlu, area, tarea
    end type
    type(pejection) :: peject(inepl)
    type flanks
      real(kind=8) :: lb, ub, flux, sinlb, sinlub
    end type
    type(flanks)    :: beta(inepl)

!   ----------------- /see/
    integer(kind=4) :: isse
    real(kind=8)    :: pemax(ispec), deltaemax(ispec)
    real(kind=8)    :: sey(12,ispec), seygl(12,ispec)

!   ----------------- /count/
    type(counter)   :: gcount(2)
    integer(kind=4) :: mpi_type_count
    integer(kind=4) :: mpi_sum_count
!    integer(kind=8) :: numnpc(inpc),ionnumber(inpc),influx(inpc)
!    integer(kind=8) :: noflux(inpc),nfront(ispec),nend(ispec)
!    integer(kind=8) :: isflux(inpc,ispec)
!    integer(kind=8) :: nesc(ispec)
    integer(kind=4) :: infhistintgrtd(0:inpc,inpc,ispec)
    integer(kind=4) :: outfluxintgrtd(inpc,ispec)
    integer(kind=4) :: h5fcount, h5jcount
    integer(kind=4) :: emarlxt
    real(kind=8)    :: selfpintgrtd(inpc)
    real(kind=8)    :: rhoindintgrtd(inpc)
    real(kind=8)    :: vidata(10)
    real(kind=8)    :: phiprintgrtd(20)

!   ----------------- /randcm/
    integer(kind=4) :: iseed
    integer(kind=4) :: nrawork1=100000,nrawork2=100000
    real(kind=8),allocatable :: dranu(:)
    real(kind=8)    :: drawork1(100000),drawork2(100000)

!   ----------------- /lblval/
    real(kind=8)    :: dmiss=-9999.d0

!   ----------------- /mpi/
    integer(kind=4) :: nnode,myid,nodes(3)
    integer(kind=4) :: snode, ssnode
    integer(kind=4) :: nfftofs
    integer(kind=4),allocatable :: bared(:,:,:)
    integer(kind=4),allocatable :: medges(:,:,:)
!    integer(kind=4) :: wsend,snend,tbend
!    integer(kind=4) :: fftend
    integer(kind=4),allocatable :: istatus(:,:), ireqs(:)
    real(kind=8)    :: fftofs

!   ----------------- /verbose/
    integer(kind=4) :: EMSES_verbose, OHHELP_verbose
    integer(kind=4) :: OHHELP_stats, OHHELP_repiter
    
!   ----------------- /CoToCoA/
    integer :: local_comm


  end module
