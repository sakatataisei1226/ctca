## esorem
    emflag: emflag = 1  Full-electromagnetic treatment
    emflag = 0  Electrostatic approximation
## jobcon
    jobnum specifying if this is the continuous job or not
    nstep number of time steps

## digcon definition of some diagnostic parameters

## plasma
    wp(1:2) electron/ion plasma frequencies (EMSES-U)
    wc electron cyclotron frequency (EMSES-U)
    cv speed of light (EMSES-U)

## tmgrid
    dt time step width (EMSES-U)
    nx, ny, nz number of grid points in x, y, and z directions

## system
    nspec number of plasma species
    nfbnd, npbnd boundary treatments of field & particles (0: periodic, 2:free)
    mtd_vbnd boundary cond. for potential (0: periodic, 1: Dirichlet, 2: Neumann)

## intp
    qm(1:2) charge-to-mass ratios (EMSES-U)
    npin(1:2) numbers of initial macro-(numerical-)particles
    path(1:2) thermal velocities along static B-field (EMSES-U)
    peth(1:2) thermal velocities perpendicular to B- field (EMSES-U)
    vdri(1:2) magnitudes of flow velocities (EMSES-U)
    vdthz(1:2), vdthxy(1:2) definition of flow direction (in degrees)

## ptcond
    npc number of solid bodies
    npcg number of conducting bodies
    pcgs, ccgs definition for grouping bodies to form one conductor
    mtd_vchg treatment of potential (0: floating, -1: fixed)
    pfixed fixed body potential (if mtd_vchg=-1) (EMSES-U)
    geotype body shape (0&1: rectangular, 2: cylinder, 3: sphere)
    bdyalign, bdyradius, bdyedge, bdycoord definition of cylinder/sphere geometry
    {x,y,z}{l,u}pc definition of rectangular geometry

## mpi
    nodes(1:3) number of mpi processes in x, y, and z directions