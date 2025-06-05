!
module namels
!
!   ____________________________________________________________
!
!                    M O D U L E   N A M E L S
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .        module for common declaration of namelist         .
!   ............................................................
!
!-------------------- common declaration of namelist
    use paramt
    use allcom
    namelist /real/ deltax, deltat, massratio, density, etemp, itemp, phtemp, &
        flowvx, flowvy, flowvz, magb, phiz, phixy, jph, thetaz, thetaxy

    namelist /realcv/ cvsim, weightph

    namelist /esorem/ emflag, nflag_testp

    namelist /jobcon/ jobnum, nstep, lstep, maxt

    namelist /plasma/ wp, wc, cv, gfactor, phixy, phiz, b0x, b0y, b0z, b0, wpmax, &
        e0x, e0y, e0z, rho0, denmod, denk, omegalw

    namelist /tmgrid/ dt, dtf, dr, nx, ny, nz

    namelist /system/ nspec, mltstp, mltstpf, juncan, alpha, ionchg, &
        jxfltr, jyfltr, jzfltr, ipexfl, ipeyfl, ipezfl, &
        nfltrx, nfltry, nfltrz, nfbnd, &
        nxl, nxr, nyl, nyr, nzl, nzr, imask, npbnd, &
        xlcol, xucol, ylcol, yucol, zlcol, zucol, mfpath, &
        nflag_ecrct, pftmode, mtd_vbnd, iphiref, &
        boundary_conditions

    namelist /digcon/ i1hdf, i2hdf, i3hdf, isort, idim_sp, ifxyz, ijxyz, irhsp, &
        itchck, intfoc, hdfdigstart, hdfdigend, &
        iddiag, iediag, ifdiag, ijdiag, isdiag, iadiag, &
        ipadig, ipahdf, ipaxyz, ildig, ivdiag, imdig, ikdiag, &
        daverg, &
        nvx, nvy, nvz, ikxmax, ikymax, ikzmax, emarlxt

    namelist /intp/ np, npin, qm, path, spa, spe, speth, peth, f, vdx, vdy, vdz, &
        vpa, vpb, vpe, nphi, ndst, ioptd, vdri, vdthz, vdthxy, &
        lcgamma, lcbeta
!                   xe0,ye0,ze0, xed,yed,zed, nphi,ndst,ioptd

    namelist /time/ t, dt

    namelist /grid/ dr, nx, ny, nz, dri, si, &
        anx, any, anz, &
        slx, sly, slz, dkx, dky, dkz

    namelist /out/ npsum, qm, q, rm, qmr

    namelist /others/ tcs, cs

    namelist /inp/ inpf, inpb, injct, npr

    namelist /ptcond/ npc, epc2, amu2, sigma, ncond, geotype, &
        nxpc1, nypc1, nzpc1, nxpc2, nypc2, nzpc2, &
        bdyalign, bdyradius, bdyedge, bdycoord, &
        wirealign, wirehlength, &
        wirerradius, wireeradius, wireorigin, &
        xlpc, ylpc, zlpc, xupc, yupc, zupc, &
        nflag_subcell, &
        npcg, pcgs, ccgs, &
        biasp, biasc, dscaled, &
        bcval, bcfrom, bcto, &
        xc, yc, zc, rc, &
        isse, pemax, deltaemax, &
        mtd_vchg, pfixed, mingap, &
        pswper, pswstr, pswspn, pswini, &
        v_omega, v_max, sfecrrct, wrelax, modeww, oradius, &
        zssurf, zsbuf, xbowlc, ybowlc, zbowlc, rbowl, &
        xdomec, ydomec, zdomec, rdome, &
        xholec, yholec, rhole, lbhole, ubhole, flhole, &
        xlrechole, ylrechole, zlrechole, xurechole, yurechole, zurechole, &
        boundary_type, boundary_types, &
        cylinder_origin, cylinder_radius, cylinder_height, rcurv, &
        rectangle_shape, &
        sphere_origin, sphere_radius, &
        circle_origin, circle_radius, &
        cuboid_shape, &
        disk_origin, disk_height, disk_radius, disk_inner_radius, &
        conductivity, &
        plane_with_circle_hole_zlower, &
        plane_with_circle_hole_height, &
        plane_with_circle_hole_radius

    namelist /wave/ nwave, twave, ebamp, ldwv, orgwv, angwv

    namelist /scrnt/ imode, omegatw, lambdatw, amptw, thetatw

    namelist /jsrc/ njs, rjs, ajs, wjs, th0js

    namelist /testch/ ntch, dimtch, rtch, qtch, e1tch, p1tch, rcutoff, &
        nqclst, dimclst, rclst, rhopeak

    namelist /dipole/ mode_dipole, gap_amp, f_c, line_mode, n_shth_rgn, &
        ajs_amp, pgamp, ifeedst, ifeeded, nstp_oe, &
        w_c, ngap, i_gap, j_gap, k_gap, &
        resc, capc, &
        xbpcc1, xbpcc2, ybpcc1, ybpcc2, &
        zbpcc1, zbpcc2, zbpcc3, zbpcc4, &
        zbpcc5, zbpcc6, zbpcc7, zbpcc8, &
        xline, yline, Ew0, Ew, omegaw, ewmodel, nretard

    namelist /emissn/ nflag_emit, nepl, flpf, flpb, curf, curb, qp, qpr, &
        abvdem, dnsf, dnsb, imarea, ipcpl, nemd, omniemit, &
        flpfs, flpbs, curfs, curbs, dnsfs, dnsbs, &
        remf, remb, xmaxe, xmine, ymaxe, ymine, zmaxe, zmine, peject, &
        thetaz, thetaxy, plreloc, mainprop, pe_ray_cast

    namelist /mpi/ nodes

    namelist /verbose/ EMSES_verbose, OHHELP_verbose, &
        OHHELP_stats, OHHELP_repiter

end module
