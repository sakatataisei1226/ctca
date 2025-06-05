  subroutine input
!
!   ____________________________________________________________
!
!                S U B R O U T I N E   I N P U T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   This subroutine reads input information and print the  .
!   .   read-in data.                                          .
!   ............................................................

!-------------------- parameter, common, namelist blocks
      use oh_type
      use paramt
      use allcom
      use namels
      use hdf
      implicit none
!
      integer(kind=4) :: i, is, ipc, iepl, ig, ijs, iw, itch
      character(len=20) :: inputfilename
      integer(kind=4) :: iosta = 0
!
      integer(kind=4) :: inttmp(ispec)
      integer(kind=HID_T) :: fileid
      integer(kind=4) :: stats0, stats1
      integer(kind=8) :: dims(1)
!
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
      character(len=30) :: filename, dsname

!---------- unit system for em3d code -----------------------------
!
!   Maxwell's equations and equations of motion are written
!   in the form of MKS unit system.  However, the electric
!   permittivity (epsi) and the magnetic permeability (rmu)
!   are chosen as
!          epsi = 1  ,  rmu = 1/(light speed)**2   .
!   Basic quantities for scaling are
!          plasma angular frequency          wp(is)
!          system length                     (slx, sly, slz)
!          charge-to-mass ratio              qm(is)
!          number of particles in the system np(is).
!
!---------- content of input parameters and default values----------
!
! /esorem/
      nflag_testp = 0
! /real/ & /realcv/
      deltax = -1.0d0
      deltat = -1.0d0
      density = -1.0d0
      etemp = -1.0d0
      itemp = -1.0d0
      phtemp = -1.0d0
      flowvx = -1.0d0
      flowvy = -1.0d0
      flowvz = -1.0d0
      magb = -1.0d0
      jph = -1.0d0
      cvsim = -1.0d0
      weightph = -1.0d0
! /jobcon/
      nstep = 128
!                   number of total time steps in one job
      jobnum(1:2) = 0
!                   sequence number of the job.  i.e.,
!                     "0" or "1" means a new job, while "2" means
!                     a continuous job for the job "1".
!                     if "0" is specified, the data file(ft56f001)
!                     for a continuous job is not created.
      lstep = 0
      maxt = 100
!                   maximum cpu time available (in minutes)
!                   the job is automatically terminated within
!                   maxt (minutes), even before completing the
!                   time steps specified by nstep.
!                   cpu time is monitored at an interval specified
!                   by itchck.
!
! /plasma/
      wp(1) = 1.0d0
      do is = 2, ispec
          wp(is) = 0.0d0
      end do
!                   plasma frequencies of each species
      wc = -1.0d0
!                   cyclotron frequency of species 1
      cv = 100.d0
!                   speed of light
      gfactor = 0.5d0
!
      phixy = 0.0d0
!                   direction of static magnetic field
!                   in x-y plane from x-axis.
      phiz = 0.0d0
!                   direction of static magnetic field from z-axis;
!                   phixy, phiz are specified in degree.
      do is = 1, ispec
          wpmax(is) = 0.0d0
      end do
!                   if not 0, the max value of wp will be equal to wpmax
!                     at the inital stage.
      e0x = 0.0d0
      e0y = 0.0d0
      e0z = 0.0d0

      rho0 = -1.0d0

      denmod(1) = 0.0d0
      denk(1) = 0.0d0
      omegalw(1) = 0.0d0
      do is = 2, ispec
          denmod(is) = 0.0d0
          denk(is) = 0.0d0
          omegalw(is) = 0.0d0
      end do

! /tmgrid/
      dt = 0.005d0
      dtf = 1.0d0
!                   time step
!                   in the same unit as 2*pi/wc and 2*pi/wp1.
      dr = 1.0d0
!                   grid spacing (delta-x = delta-y = dr)
      nx = 8
      ny = 8
      nz = 8
!                   grid numbers in x- and y-directions
!                     which should be 2**(integer)  (integer=0,1,2,...)
!                     nx >= ny is recomended for better efficiency
!                     vector processor.
! /system/
      nspec = 1
!                   number of species; select 1,2,3
      mltstp = 1
      mltstpf = 1
!                   number of field updating loops between
!                     particle pushig loops.(multiple time step)
      juncan = 0
!                   option for cancellation of uniform components of
!                   the current.  if(juncan.eq.1), uniform componetns
!                   are subtracted from the current.
!                   if(juncan.eq.2) electrostatic component is not
!                   solved as well as the uniform componet of the
!                   current.  if(juncan.eq.-2) only electrostatic
!                   component is not solved, i.e., the uniform
!                   component of the current is retained.
!
!                   if (juncan.ge.1000) kempo is in test particle
!                   simulation mode.  the wave field e and b must be
!                   given externally by a new routine called from
!
      jxfltr = 0
      jyfltr = 0
      jzfltr = 0
!                   option for filtering the current ajx, ajy, ajz
!                     -1 : no filtering
!                     mod(j?fltr/1,2)=1  : filter for x component
!                     mod(j?fltr/2,2)=1  : filter for y component
!                     mod(j?fltr/4,2)=1  : filter for z component
      ipexfl = 0
      ipeyfl = 0
      ipezfl = 0
!                   option for filtering the particle-pushig fields
!                   pex, pey and pez.  option specification is the
!                   same as for j?fltr.
!                   the subroutine pebfld.
      nfltrx = 0
!                   option  for filtering electrostatic filed ex.
!                   if nfltrx > 0, ex field is filtered by fltrf3
!                   routine.
!                   modes of ex with wave-mode number larger than
!                   nfltrx is filtered out.
      nfltry = 0
!                   option  for filtering electrostatic filed ey.
!                   if nfltry > 0, ey field is filtered by fltrf3
!                   routine.
!                   modes of ey with wave-mode number larger than
!                   nfltry is filtered out.
      nfltrz = 0
!                   option  for filtering electrostatic filed ez.
!                   if nfltrz > 0, ez field is filtered by fltrf3
!                   routine.
!                   modes of ez with wave-mode number larger than
!                   nfltrz is filtered out.
      ionchg = 0
!                   option for background fixed ions
!                   = 0 : uniform neutralization in the system
!                   = 1 : local newtralization at each grid point
      alpha = 0.0d0
!                  increment of static bx-field at each dt.
      nfbnd(1) = 0
      nfbnd(2) = 0
      nfbnd(3) = 0
!                   field boundary treatment
!                     (1): x-axis  (2): y-axis  (3): z-axis
!                     = 0:periodic   = 1:free    = 2:damping
      nxl = 2
      nxr = 2
      nyl = 2
      nyr = 2
      nzl = 2
      nzr = 2
!                   damping region
!                     default: nxl=nx/4, nxr=nx/4, ...
      xlcol = -1.0d6
      xucol = 1.0d6
      ylcol = -1.0d6
      yucol = 1.0d6
      zlcol = -1.0d6
      zucol = 1.0d6
      do is = 1, ispec
          mfpath(is) = 0.0d0
      end do
!
      do i = 1, 5
          imask(1) = 1
          imask(2) = 1
          imask(3) = 1
          imask(4) = 0
          imask(5) = 0
      end do
!                   select component to mask
!                     if i = 1 then b-field
!                     if i = 2 then e-field (transverse)
!                     if i = 3 then e-field (longitudinal)
!                     if i = 4 then current
!                     if i = 5 then charge
!
      npbnd(1, 1) = -1
      npbnd(2, 1) = -1
      npbnd(3, 1) = -1
      do is = 2, nspec
          npbnd(1, is) = -1
          npbnd(2, is) = -1
          npbnd(3, is) = -1
      end do
!                   perticle boundary treatment
!                   (1): x-axis  (2): y-axis  (3): z-axis
!                   = 0:periodic   = 1:reflect  = 2:free   = -1:=nfbnd
!
!
!  option for ecrrct  and the method of poisson's solver
!
      nflag_ecrct = 1
!                   option for ecrrct routine
!                   = 1 : use ecrrct ( correction of electrostatic field)
!                   = 0 : neglect ecrrct
!
      pftmode = 0
!
      mtd_vbnd(1) = 0
      mtd_vbnd(2) = 0
      mtd_vbnd(3) = 0
!                   boundary treatment for the Poisson solver
!                   = 0 : Periodic
!                   = 1 : Fixed (phi=0)
!                   e.g. mtd_vbnd = 0,0,0 : all periodic boundary condition
!                   This option is valid only for  nflag_ecrct = 1
      iphiref(:, :) = 1
      boundary_conditions(:, :) = 0.0d0
!
! /digcon/
      i1hdf = -1
      i2hdf = -1
      i3hdf = -1
      intfoc = 128
      hdfdigstart = 0
      hdfdigend = -2
      iddiag = 0
      iediag = 0
      ifdiag = 1
      ijdiag = 0
      iadiag = 0
      isort(1) = 500
      idim_sp = 3
      imdig(1) = 0
      ipadig(1) = 0
      ildig(1) = 0
      ipahdf(1) = 0
      isdiag = 0
      ikdiag = 0
      daverg = 0
      ikxmax = 16
      ikymax = 16
      ikzmax = 16
      ivdiag = 0
      emarlxt = 100
      nvx = nx
      nvy = ny
      nvz = nz
      do is = 2, ispec
          isort(is) = 500
!        imdig(is) = 0
          ipadig(is) = 0
!        ildig(is) = 0
          ipahdf(is) = 0
      end do

!                   number of main loop repetitions between diagnostics
!                     ( should be integer multiples of 'mltstp')
!                     iddiag : particle distribution
!                     iediag : energy (if iediag<0, ignore damping region)
!                     ifdiag : field (component is specified by ifxyz)
!                     ijdiag : current (component is specified by ijxyz)
!                     ipadig(is) : particle diagnostics
!                                        for particle species is
!                     ipahdf(is) : data skip of species is
!                     ipaxyz(i,is) :
!                     imdig(is) : moment diagnostics
!                              for particle species is
!                     isdiag : sub-diagnostics
!                               if isdiag<0, no print-out on lp)
!                     ikdiag : k spectrum analysis
!                     ikxmax : maximun kx mode number to be analyzed.
!                     ikymax : maximun ky mode number to be analyzed.
!                     ikzmax : maximun kz mode number to be analyzed.
!                     ildig(is) : local energy diagnostics for species is.
!                     ivdiag : vpara-vperp distribution function
!                     nvx,nvy,nvz: number of divisions in vpara & vperp
!                               distribution function
      itchck = 64
!                   number of main loop repetitions between
!                     cpu time checks (see  maxt)
      ifxyz(1) = 1
      ifxyz(2) = 1
      ifxyz(3) = 1
      ifxyz(4) = 1
      ifxyz(5) = 1
      ifxyz(6) = 1
      ifxyz(7) = 1
!                   option for field data storage
!                     if(ifxyz(1).eq.1)  ex is to be stored.
!                     if(ifxyz(2).eq.1)  ey is to be stored.
!                     if(ifxyz(3).eq.1)  ez is to be stored.
!                     if(ifxyz(4).eq.1)  bx is to be stored.
!                     if(ifxyz(5).eq.1)  by is to be stored.
!                     if(ifxyz(6).eq.1)  bz is to be stored.
!                     if(ifxyz(7).eq.1) rho is to be stored.
      ijxyz(1) = 0
      ijxyz(2) = 0
      ijxyz(3) = 0
      ijxyz(4) = 0
!                   option for current data storage
!                     if(ijxyz(1).eq.1) ajx is to be stored.
!                     if(ijxyz(2).eq.1) ajy is to be stored.
!                     if(ijxyz(3).eq.1) ajz is to be stored.
      do is = 1, ispec
          irhsp(is) = 1
      enddo
!
! /intp/
      np(1) = 0
!                   maximun number of particles of each species
      qm(1) = -1.0d0
!          qm(2)  = 1.00d0
!          qm(3)  = 10.0d0
!                   charge to mass ratio of each species
      path(1) = 2.0d0
!                   parallel thermal velocities
      spa(1) = 0.0d0
!                   parallel drift velocities
      spe(1) = 0.0d0
!                   perpendicular drift velocities
      speth(1) = 0.0d0
!                   direction of perpendicular drift velocities
      peth(1) = 2.0d0
!                   perpendicular thermal velocity / sqrt((2))
!                     ( v-perp = sqrt(2.0) * peth )
      f(1) = 1.0d0
!                   loss cone factors
!           ***** distribution functions *****
!               para  =  exp -( (v - spa)**2/(2*path**2) )
!               perp  =  v *(f-1)*  exp -( v**2/perp**2 )
!               for isotropic distribution, f=1.0 and perp=path.
!            if  0 < f < 1.0 then
!               loscone distribution is realized by the
!               subtracted maxwellian distribution
!
!               perp  = 1/(1-f) * ( exp(-v**2/2*perp**2)
!                                     - exp(-v**2/(f*2*perp**2)) )
!
      vdx(1) = 0.0d0
      vdy(1) = 0.0d0
      vdz(1) = 0.0d0
!
      vdri(1) = 0.0d0
      vdthz(1) = 0.0d0
      vdthxy(1) = 0.0d0
!
      ioptd(1) = 0
!                   if 1, use optional distribution function
!           xe0(1)   = 0.0d0
!           ye0(1)   = 0.0d0
!           ze0(1)   = 0.0d0
!                   center position of spatial gaussian distribution
!           xed(1)   = 0.0d0
!           yed(1)   = 0.0d0
!           zed(1)   = 0.0d0
!                   standard deviation of spatial gaussian distribution
!                     if =0., uniform distribution is realized.
      nphi(1) = 1
!                  number of phase for single vperp value
      ndst(1) = 1
!                  number of distribution in a system
!
      vpa(1) = -16.d0
      vpb(1) = 16.d0
      vpe(1) = 32.d0
!
!             vpa - vpb  : minimum and maximum of vpara distribution
!             vpe        : maximum of vperp distribution
!
!          npsum   = number of total particles in the simulation
!                     (np=np(1)+np(2)+np(3).....)
!
      lcgamma(1) = 1.0d0
      lcbeta(1) = 1.0d0
!
      do is = 2, ispec
          np(is) = 0
          qm(is) = 10.0d0
          path(is) = 1.0d0
          spa(is) = 0.0d0
          spe(is) = 0.0d0
          speth(is) = 0.0d0
          peth(is) = 1.0d0
          f(is) = 1.0d0
          vdx(is) = 0.0d0
          vdy(is) = 0.0d0
          vdz(is) = 0.0d0
          vdri(is) = 0.0d0
          vdthz(is) = 0.0d0
          vdthxy(is) = 0.0d0
          ioptd(is) = 0
!           xe0(is) = 0.0d0
!           ye0(is) = 0.0d0
!           ze0(is) = 0.0d0
!           xed(is) = 0.0d0
!           yed(is) = 0.0d0
!           zed(is) = 0.0d0
          nphi(is) = 1
          ndst(is) = 1
          vpa(is) = -16.d0
          vpb(is) = 16.d0
          vpe(is) = 32.d0
          lcgamma(is) = 1.0d0
          lcbeta(is) = 1.0d0
      end do
!
!  /inp/
!
      npin(1) = 0
!                   number of particles at t=0
!                   if(npin(is).lt.0) you can inject particles only
!                   as beam so particles are not distributed at t=0
      inpf(1) = 0
      inpb(1) = 0
!                   inp&f is index for injection at x=0
!                   inp&b is index for injection at x=slx
!                   if '0' is specified, particles are not injected
!                   if '(1)' is specified, particles are injected
      injct(1) = 0
!                   interval time step of injection
      npr(1) = 0
!                   number of standard distribution array for
!                     injecting particles
!!!           iapdg(1) = 0
!                   number of main loop repetitions between diagnostics
!                     for active particles
!!!           idens(1) = 0
!                   index of density gradient in x-direction
!                   if '0' is specified, uniform density
!                   if '(1)' is specified, density has linear gradient
!!!           adens(1) = 0
!                   coeffecient of density gradient
!                   density(x)=d0*((1)+adens*x)
!!!           nxld(1)  = 0
!                   the begining of plasma gradient from the starting area
!!!           nxrd(1)  = 0
!                   the end of plasma gradient from the end area
!!!           nxmd1(1) = 0
!                   the corner point of plasma disposition from the starting
!                   area
!                   size of the region where density has gradient
!!!           nxmd2(1) = 0
!                    the next corner point op plasma disposition from the
!                    nxmd1
      do is = 2, ispec
          npin(is) = 0
          inpf(is) = 0
          inpb(is) = 0
          injct(is) = 0
          npr(is) = 0
!!!          iapdg(is) = 0
!!!          idens(is) = 0
!!!          adens(is) = 0
!!!           nxld(is) = 0
!!!           nxrd(is) = 0
!!!          nxmd1(is) = 0
!!!          nxmd2(is) = 0
      end do
!
! /drft/
!     vd = 0.0
!
!!! /lcbm/
!!!         beamd(1) = 0.0d0
!!!                   Local beam density for super particles ( /grid )
!!!                   if "0" is specified beam is not localized
!!!                   Note these are not real density
!!!                   Note these are not "0"
!!!                        then npin((1)-3) must be negative value "<0"
!!!                             and inp((1)-3)f or inp(1-3)b must not be "0"
!!!         nbmyt(1) = 0
!!!                   Number of top grid
!!!                   of local beam ( in y-direction )
!!!         nbmyb(1) = 0
!!!                   Number of bottom grid
!!!                   of local beam ( in y-direction )
!!!      do is=2,ispec
!!!        beamd(is) = 0.0d0
!!!        nbmyt(is) = 0
!!!        nbmyb(is) = 0
!!!      end do
!
! /ptcond/
! npc is the number of conductors
! epc2 is the relative dielectric constant
! amu2 is the relative permeabilty constant
! if ncond(ipc).eq.1------> perfect conductor
! if ncond(ipc).eq.2------> general conductor
      npc = 0
      epc2 = 0.d0
      amu2 = 0.d0
      sigma = 0.d0
      nflag_subcell(1:2) = 0
      do ipc = 1, inpc
          ncond(ipc) = 1
          geotype(ipc) = 0
      end do
      do ipc = 1, inpc
          nxpc1(ipc) = 0
          nypc1(ipc) = 0
          nzpc1(ipc) = 0
          nxpc2(ipc) = 0
          nypc2(ipc) = 0
          nzpc2(ipc) = 0
          bdyalign(ipc) = 1
          bdyradius(ipc) = 0.0d0
          bdyedge(1, ipc) = dmiss
          bdyedge(2, ipc) = dmiss
          bdycoord(1, ipc) = dmiss
          bdycoord(2, ipc) = dmiss
          bdycoord(3, ipc) = dmiss
          xlpc(ipc) = -1.0d0
          ylpc(ipc) = -1.0d0
          zlpc(ipc) = -1.0d0
          xupc(ipc) = -1.0d0
          yupc(ipc) = -1.0d0
          zupc(ipc) = -1.0d0
          wirealign(ipc) = 0
          wirehlength(ipc) = 0.0d0
          wirerradius(ipc) = 0.0d0
          wireeradius(ipc) = 0.0d0
          wireorigin(1:3, ipc) = dmiss
      end do
      npcg = npc
      do ipc = 1, npc + 1
          pcgs(ipc) = ipc - 1
          ccgs(ipc) = -1
      end do
      do ipc = 1, inpc
          biasp(ipc) = 0.0d0
          dscaled(ipc) = 1.0d0
      end do
      biasc(1)%to = 0
      bcto(1) = 0
      xc = -9999.0d0
      yc = -9999.0d0
      zc = -9999.0d0
      rc = 0.0d0
!
      zsbuf = 0d0
      zssurf = -9999.0d0
      xbowlc = -9999.0d0
      ybowlc = -9999.0d0
      zbowlc = -9999.0d0
      rbowl = -9999.0d0
      xdomec = -9999.0d0
      ydomec = -9999.0d0
      zdomec = -9999.0d0
      rdome = -9999.0d0
      xholec = -9999.0d0
      yholec = -9999.0d0
      rhole = 0.0d0
      lbhole = -9999.0d0
      ubhole = -9999.0d0
      flhole = -9999.0d0
      xlrechole(1:2) = -9999.0d0
      ylrechole(1:2) = -9999.0d0
      zlrechole(1:2) = -9999.0d0
      xurechole(1:2) = -9999.0d0
      yurechole(1:2) = -9999.0d0
      zurechole(1:2) = -9999.0d0
!
!                   This is test version of perfect conductor wall.
!                   npc : number of perfect conductor wall
!                   (nxpc1(i),nypc1(i),nzpc1(i))-(nxpc2(i),nypc2(i),nzpc2(i))
!
      isse = 1
      do is = 1, ispec
          pemax(is) = 1.0d10
          deltaemax(is) = 0.0d0
      end do
!
!  -- mtd_vchg    : Voltage control for bodies
!           0: Selfconsistent (conserving all charges on bodies)
!           1: force voltage oscillation ( v_max with v_omega frequency)
!
!
      pswper = 0
      pswstr = 0
      pswini = 0.0d0
      do ipc = 1, inpc
          mtd_vchg(ipc) = 0
          reducecm(ipc) = 0
          pfixed(ipc) = 0.0d0
          pswspn(ipc) = 0.0d0
          v_max(ipc) = 0.0d0
          v_omega(ipc) = 0.0d0
          oradius(1:3, ipc) = 0.0d0
      end do
!
      sfecrrct = 0
      wrelax = 1.0d0
      mingap = 1.0d-4
      modeww = -1
!            if(abs(modeww).eq.1) capacity matrix method is not used
!            if(abs(modeww).eq.2) capacity matrix method is used
!            if(modeww.gt.0) charge accumulation is not considered
!            if(modeww.lt.0) charge accumulation is considered

! /jsrc/
      njs = 0
      do ijs = 1, injs
          rjs(1:3, ijs) = -9999
          ajs(1:3, ijs) = (/0.0d0, 0.0d0, 0.0d0/)
          wjs(ijs) = 1.0d0
          th0js(ijs) = 0.0d0
      end do

! /wave/
      nwave = 0
      do iw = 1, inwave
          twave(iw) = 0
          ebamp(1:6, iw) = 0.0d0
          ldwv(iw) = 10.0d0
          orgwv(1:3, iw) = 0.0d0
          angwv(1:2, iw) = 0.0d0
      end do

! /scrnt/
      imode = 0
      omegatw = 0.0d0
      lambdatw = 0
      amptw = 0.0d0
      thetatw = 0.0d0

!
!!! /neut/
!!!              neutr      = 0
!!!              neuele(1)  = 1
!!!              colfq(1)   = 0.0
!!!      do is=2,ispec
!!!              neuele(is) = 0
!!!          colfq(is)  = 0.0
!!!      end do
!!!          if(neutr.eq.1) collision effect on
!!!          if(neuele(is).eq.1) electron scattering
!!!          Ion scattering can't be supported
!
! /emissn/
      do is = 1, ispec
          nflag_emit(is) = 0
          nepl(is) = 0
          qp(is) = 0.0d0
          qpr(is) = 0.0d0
          abvdem(is) = 0.0d0
          flpf(is) = 0.0d0
          flpb(is) = 0.0d0
          curf(is) = 0.0d0
          curb(is) = 0.0d0
          dnsf(is) = 0.0d0
          dnsb(is) = 0.0d0
      enddo
!
      do iepl = 1, inepl
          ipcpl(iepl) = 1
          nemd(iepl) = 1
          omniemit(iepl) = 0
          flpfs(iepl) = 0.0d0
          flpbs(iepl) = 0.0d0
          curfs(iepl) = 0.0d0
          curbs(iepl) = 0.0d0
          dnsfs(iepl) = 0.0d0
          dnsbs(iepl) = 0.0d0
          remf(iepl) = dmiss
          remb(iepl) = dmiss
          xmaxe(iepl) = dmiss
          xmine(iepl) = dmiss
          ymaxe(iepl) = dmiss
          ymine(iepl) = dmiss
          zmaxe(iepl) = dmiss
          zmine(iepl) = dmiss
          peject(iepl)%grd = dmiss
          peject(iepl)%xl = dmiss
          peject(iepl)%yl = dmiss
          peject(iepl)%zl = dmiss
          peject(iepl)%xu = dmiss
          peject(iepl)%yu = dmiss
          peject(iepl)%zu = dmiss
          peject(iepl)%xlu = 0.0d0
          peject(iepl)%ylu = 0.0d0
          peject(iepl)%zlu = 0.0d0
          peject(iepl)%area = 0.0d0
!              ephiz(iepl)=0.0d0
!              ephixy(iepl)=0.0d0
!              emitcx(iepl)=dmiss
!              emitcy(iepl)=dmiss
!              emitcz(iepl)=dmiss
      enddo
!
      do ipc = 1, inpc
          imarea(ipc) = 0.0d0
      enddo
!
      thetaz = 0.0d0
      thetaxy = 0.0d0
      plreloc = 0.0d0
      mainprop = 1.0d0
!
! /testch/
      ntch = 0
      do itch = 1, intch
          rtch(1:3, itch) = -9999.0d0
          qtch(itch) = 0.0d0
          e1tch(itch) = 0.0d0
          p1tch(itch) = 0.0d0
          rcutoff(itch) = epsilon(1.0d0)
      end do
      nqclst = 0
      do itch = 1, intch
          dimclst(1:3, itch) = 1
          rclst(1:3, itch) = -9999.0d0
          rclst(4, itch) = 1.0d0
          rclst(5, itch) = 10.0d0
          rhopeak(itch) = 0.0d0
      end do
!
! /dpolec/
      mode_dipole = 0
      line_mode = 0
      n_shth_rgn = 0
      ifeedst = 0
      ifeeded = 10000000
      nstp_oe = 1
      ngap = 0
      xbpcc1 = 0.0
      xbpcc2 = 0.0
      ybpcc1 = 0.0
      ybpcc2 = 0.0
      zbpcc1 = 0.0
      zbpcc2 = 0.0
      zbpcc3 = 0.0
      zbpcc4 = 0.0
      zbpcc5 = 0.0
      zbpcc6 = 0.0
      zbpcc7 = 0.0
      zbpcc8 = 0.0
      xline = 0.0
      yline = 0.0
      do ig = 1, ingap
          i_gap(ig) = 0
          j_gap(ig) = 0
          k_gap(ig) = 0
          gap_amp(ig) = 0.0
          ajs_amp(ig) = 0.0
          pgamp(ig) = 0.0
          resc(ig) = -1.0d0
          capc(ig) = 0.0d0
      end do
      w_c = 0.0d0
      Ew0 = 0.0d0
      Ew = 0.0d0
      omegaw = 0.0d0
      ewmodel = 0
      nretard = 0
!
!  mode_dipole  = 1 : Gap Kyodenn in 'bdyfld.f'
!  line_mode    = 1 : Virtual antenna (no particle accumulate)
!  n_shth_rgn   = 0 : Virtual thickness of sheath around the body
!  w_c          : Center frequency for pulse Ez kyuden
!                  negative value: monochro frequency
!
! /mpi/
      nodes(1) = nnode
      nodes(2) = 1
      nodes(3) = 1

! /verbose/
      EMSES_verbose = 0
      OHHELP_verbose = 0
      OHHELP_stats = 0
      OHHELP_repiter = 0

!  time
      t = 0.0d0
!
!    notice : other parameters are explained in the subroutine "inital"
!
!        never forget to set real parameters as double precision!
!                       good luck!!!
!
!-------------------------------------------------------------------
!
!--- read data ---
      call get_command_argument(1, inputfilename)
!
      open (1, file=inputfilename)
!
      rewind (1); read (1, nml=real, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=real not found"
!
      rewind (1); read (1, nml=realcv, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=realcv not found"
!
      rewind (1); read (1, nml=esorem, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=esorem not found"
!
      rewind (1); read (1, nml=jobcon, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=jobcon not found"
!
      rewind (1); read (1, nml=digcon, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=digcon not found"
!
      rewind (1); read (1, nml=plasma, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=plasma not found"
!
      rewind (1); read (1, nml=tmgrid, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=tmgrid not found"
!
      rewind (1); read (1, nml=system, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=system not found"
!
      rewind (1); read (1, nml=intp, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=intp not found"
!
      rewind (1); read (1, nml=inp, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=inp not found"
!
!           rewind(1);read(1,drft,IOSTAT=iosta)
!           if(myid.eq.0.and.iosta.eq.-1) &
!          &  print*, "Warning.Input: nml=drft not found"
!
      rewind (1); read (1, nml=ptcond, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=ptcond not found"
!
      rewind (1); read (1, nml=jsrc, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=jsrc not found"
!
      rewind (1); read (1, nml=wave, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=wave not found"
!
      rewind (1); read (1, nml=scrnt, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=scrnt not found"
!
      rewind (1); read (1, nml=emissn, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=emissn not found"
!
      rewind (1); read (1, nml=testch, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=testch not found"
!
      rewind (1); read (1, nml=dipole, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=dipole not found"
!
      rewind (1); read (1, nml=mpi, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=mpi not found"
!
      rewind (1); read (1, nml=verbose, IOSTAT=iosta)
      if (myid .eq. 0 .and. iosta .eq. -1) &
     &  print *, "Warning.Input: nml=verbose not found"
!
      close (1)

      if (myid .eq. 0) print *, "nx,ny,nz =", nx, ny, nz
!-------------------- job type
      if (jobnum(1) .lt. 0 .or. jobnum(1) .gt. 2 .or. &
     &   jobnum(2) .lt. 0 .or. jobnum(2) .gt. 1) then
          if (myid .eq. 0) print *, "jobnum(1) should be in a range [0, 2]"
          if (myid .eq. 0) print *, "jobnum(2) should be in a range [0, 1]"
          if (myid .eq. 0) print *, "jobnum is set to 0"
          jobnum(1:2) = 0
      end if
      if (emflag .eq. 0 .and. jobnum(1) .eq. 2) then
          if (myid .eq. 0) print *, "ES sim. can only load ES-snapshot data"
          if (myid .eq. 0) print *, "jobnum(1) is decremented by one"
          jobnum(1) = jobnum(1) - 1
      end if

!------------------------------
      if (hdfdigend .eq. -2) hdfdigend = nstep
      if (hdfdigstart .lt. 0) then
          hdfdigstart = hdfdigend + hdfdigstart + 1
      end if
      nfsnap = 0
      njsnap = 0
      npsnap(1:nspec) = 0
      do istep = 0, nstep
          if (istep .ge. hdfdigstart .and. istep .le. hdfdigend) then
              if (ifdiag .gt. 0 .and. mod(istep, ifdiag) .eq. 0) then
!            nfsnap = nfsnap + 1
              end if
              if (ijdiag .gt. 0 .and. mod(istep, ijdiag) .eq. 0) then
!            njsnap = njsnap + 1
              end if
              do is = 1, nspec
                  if (ipahdf(is) .gt. 0 .and. mod(istep, ipahdf(is)) .eq. 0) then
!              npsnap(is) = npsnap(is) + 1
                  end if
              end do
          end if
      end do

!------------------------------
      dt = dt*dtf

!------------------------------
      do is = 1, nspec
          if (np(is) .eq. 0) then
              np(is) = npin(is)
          end if
      end do

!------------------------------
      if (abs(xlcol(2)) .eq. 1.0d6) xlcol(2) = xlcol(1)
      if (abs(xucol(2)) .eq. 1.0d6) xucol(2) = xucol(1)
      if (abs(ylcol(2)) .eq. 1.0d6) ylcol(2) = ylcol(1)
      if (abs(yucol(2)) .eq. 1.0d6) yucol(2) = yucol(1)
      if (abs(zlcol(2)) .eq. 1.0d6) zlcol(2) = zlcol(1)
      if (abs(zucol(2)) .eq. 1.0d6) zucol(2) = zucol(1)

!------------------------------
      Ew = Ew(:, :, :)*Ew0

!------------------------------
      if (jobnum(1) .gt. 0) then
          if (jobnum(1) .eq. 1) then
              write (filename, '(a,i4.4,a)') './SNAPSHOT0/esdat', myid, '.h5'
          else
              write (filename, '(a,i4.4,a)') './SNAPSHOT0/emdat', myid, '.h5'
          end if
          call hdfopen(filename, fileid, DFACC_READ)
!
          dsname = 'nspec'
          dims(1) = 1
          call read1i(fileid, dsname, dims(1:1), inttmp(1:1), stats0, stats1)
          if (inttmp(1) .ne. nspec) then
              if (myid .eq. 0) print *, "nspec is inconsistent with continued job data: STOP"
              stop
          end if
!
          dsname = 'np'
          dims(1) = nspec
          call read1i(fileid, dsname, dims(1:1), inttmp(1:nspec), stats0, stats1)
!        do is=1,nspec
!          if(inttmp(is).gt.int((100.0d0+MAXFRAC)*0.01d0*np(is))) then
!            if(myid.eq.0) print*, "np(", is, ") is changed to... ", &
!           &                      int(inttmp(is)*100.0d0/(100.0d0+MAXFRAC)) + 1
!            np(is) = int(inttmp(is)*100.0d0/(100.0d0+MAXFRAC)) + 1
!          end if
!        end do
!
          dsname = 'npc'
          dims(1) = 1
          call read1i(fileid, dsname, dims(1:1), inttmp(1:1), stats0, stats1)
          if (inttmp(1) .ne. npc) then
              if (myid .eq. 0) print *, "npc is inconsistent with continued job data: STOP"
              stop
          end if
!
          call hdfclose(fileid, stats0)
      end if

!------------------------------
      do ipc = 1, npc
          if (geotype(ipc) .eq. 2) then
              cylinder(ipc)%align = bdyalign(ipc)
              cylinder(ipc)%radius = bdyradius(ipc)
              cylinder(ipc)%edge(1) = bdyedge(1, ipc)
              cylinder(ipc)%edge(2) = bdyedge(2, ipc)
              cylinder(ipc)%axis(1) = bdycoord(1, ipc)
              cylinder(ipc)%axis(2) = bdycoord(2, ipc)
              if (cylinder(ipc)%align .eq. 1) then
                  xlpc(ipc) = cylinder(ipc)%edge(1)
                  xupc(ipc) = cylinder(ipc)%edge(2)
              else if (cylinder(ipc)%align .eq. 2) then
                  ylpc(ipc) = cylinder(ipc)%edge(1)
                  yupc(ipc) = cylinder(ipc)%edge(2)
              else
                  zlpc(ipc) = cylinder(ipc)%edge(1)
                  zupc(ipc) = cylinder(ipc)%edge(2)
              end if
          else if (geotype(ipc) .eq. 3) then
              sphere(ipc)%radius = bdyradius(ipc)
              sphere(ipc)%center(1) = bdycoord(1, ipc)
              sphere(ipc)%center(2) = bdycoord(2, ipc)
              sphere(ipc)%center(3) = bdycoord(3, ipc)
          end if
!
          if (abs(wirealign(ipc)) .ge. 1 .and. abs(wirealign(ipc)) .le. 3) then
              boom(ipc)%align = wirealign(ipc)
              boom(ipc)%hlength = wirehlength(ipc)
              if (wirerradius(ipc) .ne. 0.0d0) then
                  boom(ipc)%rradius = wirerradius(ipc)
                  if (myid .eq. 0) print *, "ipc, boomrrad =", ipc, boom(ipc)%rradius
              else
                  boom(ipc)%rradius = cylinder(ipc)%radius
                  if (myid .eq. 0) print *, "ipc, boomrrad =", ipc, boom(ipc)%rradius
              end if
              if (wireeradius(ipc) .ne. 0.0d0) then
                  boom(ipc)%eradius = wireeradius(ipc)
                  if (myid .eq. 0) print *, "ipc, boomerad =", ipc, boom(ipc)%eradius
              else
                  boom(ipc)%eradius = 1.0d0/exp(0.5d0*dacos(-1.0d0))
                  if (myid .eq. 0) print *, "ipc, boomerad =", ipc, boom(ipc)%eradius
              end if
              boom(ipc)%origin(1:3) = wireorigin(1:3, ipc)
          else
              boom(ipc)%align = 0
              boom(ipc)%hlength = 0.0d0
              boom(ipc)%rradius = 0.0d0
              boom(ipc)%eradius = 0.0d0
              boom(ipc)%origin(1:3) = dmiss
          end if
      end do
!
      if (ccgs(1) .lt. 0) then
          ccgs(:) = pcgs(:)
      end if
!
      biasc(1:3)%from = bcfrom(1:3)
      biasc(1:3)%to = bcto(1:3)
      biasc(1:3)%val = bcval(1:3)
      if (myid .eq. 0) then
          print *, "biasc%from =", biasc(1:3)%from
          print *, "biasc%to   =", biasc(1:3)%to
          print *, "biasc%val  =", biasc(1:3)%val
      end if

      open (2, file="plasma.out", status='replace')
      write (2, '(A)') "&param"
      write (2, '(A,I,A,I,A,I,A)') "  ndim      = ", nx, ",", ny, ",", nz, ","
      write (2, '(A,I,A,I,A,I,A)') "  nodes     = ", nodes(1), ",", nodes(2), ",", nodes(3), ","
      if (ifdiag .gt. 0) then
          write (2, '(A,I,A,I,A)') "  seindices = ", 0, ",", int((hdfdigend - hdfdigstart)/ifdiag), ","
      else
          write (2, '(A,I,A,I,A)') "  seindices = ", 0, ",", 0, ","
      end if
      write (2, '(A)') "/"
      close (2)
!--- nstep must be a multiple of mltstp --
      nstep = int(DBLE(nstep)/DBLE(mltstp))*mltstp
!----------------------------------------
      do is = 1, nspec
          q(is) = qp(is)
      enddo
      do is = 1, nspec
          if (npbnd(1, is) .eq. -1) npbnd(1, is) = nfbnd(1)
          if (npbnd(2, is) .eq. -1) npbnd(2, is) = nfbnd(2)
          if (npbnd(3, is) .eq. -1) npbnd(3, is) = nfbnd(3)
      end do
      if (xline .eq. 0.0 .or. yline .eq. 0.0) then
          xline = dble(i_gap(1))
          yline = dble(j_gap(2))
      end if
!----------------------------------------
      elemchg = 1.602177d-19
      massele = 9.109382d-31
      massion = massele*massratio
      epsiln0 = 8.854188d-12
      boltfct = 1.380651E-23
      coe_v = cvsim/2.99792458d8
      coe_x = 1.0d0/deltax/1.0d-2
      coe_t = coe_x/coe_v
      coe_w = 1.0d0/coe_t
      coe_ep = 1.0d0/epsiln0
      coe_qm = massele/elemchg
      coe_j = coe_ep/coe_qm*coe_v*coe_v*coe_v/coe_x/coe_x
      if (deltax .gt. 0.0d0) then
          if (myid .eq. 0) print *, "Accept Real Parameters"
          cv = cvsim
          dt = deltat*coe_t
          qm(2) = -qm(1)/massratio
          wp(1) = sqrt(density*1d6*elemchg*elemchg &
         &             /massele/epsiln0)*coe_w
          wp(2) = wp(1)/sqrt(massratio)
          wp(3) = wp(1)
          wc = -magb*elemchg/massele*coe_w
          path(1) = sqrt(etemp*boltfct/massele)*coe_v
          peth(1) = path(1)
          path(2) = sqrt(itemp*boltfct/massion)*coe_v
          peth(2) = path(2)
          path(3) = sqrt(phtemp*boltfct/massele)*coe_v
          peth(3) = path(3)
          vdx(1) = flowvx*coe_v
          vdx(2) = vdx(1)
          vdy(1) = flowvy*coe_v
          vdy(2) = vdy(1)
          vdz(1) = flowvz*coe_v
          vdz(2) = vdz(1)
          spa(1) = 0.0d0; spe(1) = 0.0d0
          spa(2) = 0.0d0; spe(2) = 0.0d0
          qp(1) = wp(1)*wp(1)*nx*ny*nz/npin(1)/qm(1)
          qp(3) = qp(1)*weightph
          q(3) = qp(3)
          if (jph .gt. 0.0d0) then
              nflag_emit(3) = 1
          else
              nflag_emit(3) = 0
              nspec = 2
          end if
          flpf(3) = abs(jph*coe_j/q(3))
          if (myid .eq. 0) print *, "  dt =", dt
          if (myid .eq. 0) print *, "  qm(1:3) =", qm(1:3)
          if (myid .eq. 0) print *, "  wp(1:3) =", wp(1:3)
          if (myid .eq. 0) print *, "  wc =", wc
          if (myid .eq. 0) print *, "  path(1:3) =", path(1:3)
          if (myid .eq. 0) print *, "  peth(1:3) =", peth(1:3)
          if (myid .eq. 0) print *, "  vdx(1:3) =", vdx(1:3)
          if (myid .eq. 0) print *, "  vdy(1:3) =", vdy(1:3)
          if (myid .eq. 0) print *, "  vdz(1:3) =", vdz(1:3)
          if (myid .eq. 0) print *, "  q(3) =", q(3)
          if (myid .eq. 0) print *, "  flpf(3) =", flpf(3)
      end if
!--- print out of the input data ---
!
!          write(6,100)
!  100     format(1h ///1h ,10x,'input data list'///)
!
!           write(6,jobcon)
!           write(6,plasma)
!           write(6,tmgrid)
!           write(6,system)
!           write(6,digcon)
!           write(6,intp)
!
!         if(iabs(nwave).ge.1) write(6,wave1)
!         if(iabs(nwave).ge.2) write(6,wave2)
!
!ccYasugi
!           write(6,crnt)
!           write(6,inp)
!           write(6,wall)
!           write(6,lcbm)
!           write(6,ptcond)
!           write(6,scrnt)
!           write(6,neut)
!ccYSGI
!           write(6,beamij)
!           write(6,dipole)
!c
!ccc   endif
      if (myid .eq. 0) print *, "dt =", dt
      if (myid .eq. 0) print *, "Ew =", Ew(:, :, 1)

      return
  end subroutine

