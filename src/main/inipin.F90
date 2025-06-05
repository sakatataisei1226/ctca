#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
subroutine inipin
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   I N I P I N
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .    this routine gives an initial setting of injected     .
!   .    particles.                                            .
!   .    available for injection at ex/internal boundary.      .
!   ............................................................
!
!
!-------------------- parameter common block
    use oh_type
    use paramt
    use allcom
#define MCW local_comm
    implicit none
!
    integer(kind=8) :: m, ns, ne, nee
    integer(kind=8) :: i, j, k
    integer(kind=8) :: nns, nne, npq(ispec)
    integer(kind=8) :: nprfrac1, nprfrac2
    integer(kind=4) :: l, is, ipc, iepl
    integer(kind=4) :: iran
    integer(kind=4) :: icon, iconf, iconb
    integer(kind=4) :: ierr
    integer(kind=4) :: xl, xu, yl, yu, zl, zu
    integer(kind=4) :: ipcplane, align
    real(kind=8) :: psi1, psi2
    real(kind=8) :: v0, rad, rad0, vs, pdf, pdb, aa, bb
    real(kind=8) :: radius, emarea(inepl)
    real(kind=8) :: betaul(2, 2, 2)
    real(kind=8) :: betalist(10), betatmp
    real(kind=8) :: tanthetaz, costhetaxy, sinthetaxy
    real(kind=8) :: tantzcostxy, tantzsintxy

!--------------------
    vnormf => vxf
    vtangf => vyf
    runiform => vzf
    vnormc => vxb
    vtangc => vyb
    psicos => vzb
    vnorms => vxr
    vtangs1 => vyr
    vtangs2 => vzr

!--------------------
    xl = sdoms(1, 1, sdid(1) + 1)
    xu = sdoms(2, 1, sdid(1) + 1)
    yl = sdoms(1, 2, sdid(1) + 1)
    yu = sdoms(2, 2, sdid(1) + 1)
    zl = sdoms(1, 3, sdid(1) + 1)
    zu = sdoms(2, 3, sdid(1) + 1)

!-------------------- definition of sunlight incident angle
!
!
!           z                inital direction of sunlight : +z
!             |
!             |.
!             | .
!             |  .
!             |   . incident angle
!             |   /.
!          thetaz/ .
!             | /  .
!             |/   .
!             +----.---.------- y
!            / \   .  .
!        thetaxy\  . .
!          /     \ ..
!         /.......\.
!        /
!       x
!
!     --------------- 0 <= thetaz <= 180.0d0, 0 <= thetaxy < 360.0d0
    thetaz = 180.0d0 - abs(mod(abs(thetaz), 360.0d0) - 180.0d0)
    thetaxy = mod(thetaxy + 360.0d0, 360.0d0)
!
!     --------------- angle formed with x-, and y-axes
    thetax = acos(sin(thetaz*pi/180.0d0)*cos(thetaxy*pi/180.0d0))
    thetay = acos(sin(thetaz*pi/180.0d0)*sin(thetaxy*pi/180.0d0))
!
!     --------------- projection to xy-, xz-, and yz-planes
    psixy = thetaxy*pi/180.0d0
    if (thetaz .ne. 90.0d0) then
        if (thetaxy .le. 90.0d0 .or. thetaxy .gt. 270.0d0) then
            psizx = atan(tan(thetaz*pi/180.0d0)*cos(thetaxy*pi/180.0d0))
            psizx = mod(psizx + pi, pi)
        else if (thetaxy .gt. 90.0d0 .or. thetaxy .le. 270.0d0) then
            psizx = atan(tan(thetaz*pi/180.0d0)*cos(thetaxy*pi/180.0d0))
            psizx = mod(psizx + pi, pi) + pi
        end if
        if (thetaxy .lt. 180.0d0) then
            psizy = atan(tan(thetaz*pi/180.0d0)*sin(thetaxy*pi/180.0d0))
            psizy = mod(psizy + pi, pi)
        else if (thetaxy .ge. 180.0d0) then
            psizy = atan(tan(thetaz*pi/180.0d0)*sin(thetaxy*pi/180.0d0))
            psizy = mod(psizy + pi, pi) + pi
        end if
    else
        if (thetaxy .lt. 90.0d0 .or. thetaxy .gt. 270.0d0) then
            psizx = 0.5d0*pi
        else if (thetaxy .gt. 90.0d0 .or. thetaxy .lt. 270.0d0) then
            psizx = 1.5d0*pi
        else
            psizx = 0.0d0
        end if
        if (thetaxy .lt. 180.0d0) then
            psizy = 0.5d0*pi
        else if (thetaxy .gt. 180.0d0) then
            psizy = 1.5d0*pi
        else
            psizy = 0.0d0
        end if
    end if
!
!     --------------- unit conversion from degree to radian
    thetaz = thetaz*pi/180.0d0
    thetaxy = thetaxy*pi/180.0d0

!     ---------------
    countertxy = thetaxy + pi
    ltxy = thetaxy + pi*0.5d0
    utxy = thetaxy + pi*1.5d0

!     ---------------
    tanthetaz = tan(thetaz)
    costhetaxy = cos(thetaxy)
    sinthetaxy = sin(thetaxy)
    tantzcostxy = tanthetaz*costhetaxy
    tantzsintxy = tanthetaz*sinthetaxy

!--------------------
    dray(1:3) = -(/sin(thetaz)*cos(thetaxy), &
                   sin(thetaz)*sin(thetaxy), &
                   cos(thetaz)/)
    laray = dray(1)*dray(1) + dray(2)*dray(2) + dray(3)*dray(3)
    dcdome = dray(1)*xdomec + dray(2)*ydomec + dray(3)*zdomec
    dcbowl = dray(1)*xbowlc + dray(2)*ybowlc + dray(3)*zbowlc

!-------------------- effective emitting plane for line antenna
    nne = 0
    do is = 1, nspec
        nns = nne + 1
        nne = nne + nepl(is)
        do iepl = nns, nne
            ipcplane = ipcpl(iepl)
            if (peject(iepl)%grd .eq. dmiss) then
                if (nemd(iepl) .eq. 1 .and. xmaxe(iepl) .eq. dmiss) then
                    peject(iepl)%grd = xupc(ipcplane)
                    peject(iepl)%xl = dmiss
                    peject(iepl)%xu = dmiss
                    peject(iepl)%yl = ylpc(ipcplane)
                    peject(iepl)%yu = yupc(ipcplane)
                    peject(iepl)%zl = zlpc(ipcplane)
                    peject(iepl)%zu = zupc(ipcplane)
                    peject(iepl)%area = (peject(iepl)%yu - peject(iepl)%yl)*dr &
                                        &                     *(peject(iepl)%zu - peject(iepl)%zl)*dr
                else if (nemd(iepl) .eq. 1 .and. xmaxe(iepl) .ne. dmiss) then
                    peject(iepl)%grd = xmaxe(iepl)
                    peject(iepl)%xl = dmiss
                    peject(iepl)%xu = dmiss
                    peject(iepl)%yl = ymine(iepl)
                    peject(iepl)%yu = ymaxe(iepl)
                    peject(iepl)%zl = zmine(iepl)
                    peject(iepl)%zu = zmaxe(iepl)
                    peject(iepl)%area = (peject(iepl)%yu - peject(iepl)%yl)*dr &
                                        &                     *(peject(iepl)%zu - peject(iepl)%zl)*dr
                else if (nemd(iepl) .eq. -1 .and. xmine(iepl) .eq. dmiss) then
                    peject(iepl)%grd = xlpc(ipcplane)
                    peject(iepl)%xl = dmiss
                    peject(iepl)%xu = dmiss
                    peject(iepl)%yl = ylpc(ipcplane)
                    peject(iepl)%yu = yupc(ipcplane)
                    peject(iepl)%zl = zlpc(ipcplane)
                    peject(iepl)%zu = zupc(ipcplane)
                    peject(iepl)%area = (peject(iepl)%yu - peject(iepl)%yl)*dr &
                                        &                     *(peject(iepl)%zu - peject(iepl)%zl)*dr
                else if (nemd(iepl) .eq. -1 .and. xmine(iepl) .ne. dmiss) then
                    peject(iepl)%grd = xmine(iepl)
                    peject(iepl)%xl = dmiss
                    peject(iepl)%xu = dmiss
                    peject(iepl)%yl = ymine(iepl)
                    peject(iepl)%yu = ymaxe(iepl)
                    peject(iepl)%zl = zmine(iepl)
                    peject(iepl)%zu = zmaxe(iepl)
                    peject(iepl)%area = (peject(iepl)%yu - peject(iepl)%yl)*dr &
                                        &                     *(peject(iepl)%zu - peject(iepl)%zl)*dr
                else if (nemd(iepl) .eq. 2 .and. ymaxe(iepl) .eq. dmiss) then
                    peject(iepl)%grd = yupc(ipcplane)
                    peject(iepl)%xl = xlpc(ipcplane)
                    peject(iepl)%xu = xupc(ipcplane)
                    peject(iepl)%yl = dmiss
                    peject(iepl)%yu = dmiss
                    peject(iepl)%zl = zlpc(ipcplane)
                    peject(iepl)%zu = zupc(ipcplane)
                    peject(iepl)%area = (peject(iepl)%xu - peject(iepl)%xl)*dr &
                                        &                     *(peject(iepl)%zu - peject(iepl)%zl)*dr
                else if (nemd(iepl) .eq. 2 .and. ymaxe(iepl) .ne. dmiss) then
                    peject(iepl)%grd = ymaxe(iepl)
                    peject(iepl)%xl = xmine(iepl)
                    peject(iepl)%xu = xmaxe(iepl)
                    peject(iepl)%yl = dmiss
                    peject(iepl)%yu = dmiss
                    peject(iepl)%zl = zmine(iepl)
                    peject(iepl)%zu = zmaxe(iepl)
                    peject(iepl)%area = (peject(iepl)%xu - peject(iepl)%xl)*dr &
                                        &                     *(peject(iepl)%zu - peject(iepl)%zl)*dr
                else if (nemd(iepl) .eq. -2 .and. ymine(iepl) .eq. dmiss) then
                    peject(iepl)%grd = ylpc(ipcplane)
                    peject(iepl)%xl = xlpc(ipcplane)
                    peject(iepl)%xu = xupc(ipcplane)
                    peject(iepl)%yl = dmiss
                    peject(iepl)%yu = dmiss
                    peject(iepl)%zl = zlpc(ipcplane)
                    peject(iepl)%zu = zupc(ipcplane)
                    peject(iepl)%area = (peject(iepl)%xu - peject(iepl)%xl)*dr &
                                        &                     *(peject(iepl)%zu - peject(iepl)%zl)*dr
                else if (nemd(iepl) .eq. -2 .and. ymine(iepl) .ne. dmiss) then
                    peject(iepl)%grd = ymine(iepl)
                    peject(iepl)%xl = xmine(iepl)
                    peject(iepl)%xu = xmaxe(iepl)
                    peject(iepl)%yl = dmiss
                    peject(iepl)%yu = dmiss
                    peject(iepl)%zl = zmine(iepl)
                    peject(iepl)%zu = zmaxe(iepl)
                    peject(iepl)%area = (peject(iepl)%xu - peject(iepl)%xl)*dr &
                                        &                     *(peject(iepl)%zu - peject(iepl)%zl)*dr
                else if (nemd(iepl) .eq. 3 .and. zmaxe(iepl) .eq. dmiss) then
                    peject(iepl)%grd = zupc(ipcplane)
                    peject(iepl)%xl = xlpc(ipcplane)
                    peject(iepl)%xu = xupc(ipcplane)
                    peject(iepl)%yl = ylpc(ipcplane)
                    peject(iepl)%yu = yupc(ipcplane)
                    peject(iepl)%zl = dmiss
                    peject(iepl)%zu = dmiss
                    peject(iepl)%area = (peject(iepl)%xu - peject(iepl)%xl)*dr &
                                        &                     *(peject(iepl)%yu - peject(iepl)%yl)*dr
                else if (nemd(iepl) .eq. 3 .and. zmaxe(iepl) .ne. dmiss) then
                    peject(iepl)%grd = zmaxe(iepl)
                    peject(iepl)%xl = xmine(iepl)
                    peject(iepl)%xu = xmaxe(iepl)
                    peject(iepl)%yl = ymine(iepl)
                    peject(iepl)%yu = ymaxe(iepl)
                    peject(iepl)%zl = dmiss
                    peject(iepl)%zu = dmiss
                    peject(iepl)%area = (peject(iepl)%xu - peject(iepl)%xl)*dr &
                                        &                     *(peject(iepl)%yu - peject(iepl)%yl)*dr
                else if (nemd(iepl) .eq. -3 .and. zmine(iepl) .eq. dmiss) then
                    peject(iepl)%grd = zlpc(ipcplane)
                    peject(iepl)%xl = xlpc(ipcplane)
                    peject(iepl)%xu = xupc(ipcplane)
                    peject(iepl)%yl = ylpc(ipcplane)
                    peject(iepl)%yu = yupc(ipcplane)
                    peject(iepl)%zl = dmiss
                    peject(iepl)%zu = dmiss
                    peject(iepl)%area = (peject(iepl)%xu - peject(iepl)%xl)*dr &
                                        &                     *(peject(iepl)%yu - peject(iepl)%yl)*dr
                else if (nemd(iepl) .eq. -3 .and. zmine(iepl) .ne. dmiss) then
                    peject(iepl)%grd = zmine(iepl)
                    peject(iepl)%xl = xmine(iepl)
                    peject(iepl)%xu = xmaxe(iepl)
                    peject(iepl)%yl = ymine(iepl)
                    peject(iepl)%yu = ymaxe(iepl)
                    peject(iepl)%zl = dmiss
                    peject(iepl)%zu = dmiss
                    peject(iepl)%area = (peject(iepl)%xu - peject(iepl)%xl)*dr &
                                        &                     *(peject(iepl)%yu - peject(iepl)%yl)*dr
                end if
            end if
!
            if (nemd(iepl) .eq. 1) then
                if (peject(iepl)%grd .lt. xl .or. peject(iepl)%grd .ge. xu) then
                    peject(iepl)%grd = dmiss
                end if
            else if (nemd(iepl) .eq. -1) then
                if (peject(iepl)%grd .le. xl .or. peject(iepl)%grd .gt. xu) then
                    peject(iepl)%grd = dmiss
                end if
            else if (nemd(iepl) .eq. 2) then
                if (peject(iepl)%grd .lt. yl .or. peject(iepl)%grd .ge. yu) then
                    peject(iepl)%grd = dmiss
                end if
            else if (nemd(iepl) .eq. -2) then
                if (peject(iepl)%grd .le. yl .or. peject(iepl)%grd .gt. yu) then
                    peject(iepl)%grd = dmiss
                end if
            else if (nemd(iepl) .eq. 3) then
                if (peject(iepl)%grd .lt. zl .or. peject(iepl)%grd .ge. zu) then
                    peject(iepl)%grd = dmiss
                end if
            else if (nemd(iepl) .eq. -3) then
                if (peject(iepl)%grd .le. zl .or. peject(iepl)%grd .gt. zu) then
                    peject(iepl)%grd = dmiss
                end if
            end if
!
            if (peject(iepl)%xl .ge. xu .or. peject(iepl)%xu .le. xl) then
                peject(iepl)%xl = dmiss
                peject(iepl)%xu = dmiss
            else
                if (peject(iepl)%xl .lt. xl) then
                    peject(iepl)%xl = xl
                end if
                if (peject(iepl)%xu .gt. xu) then
                    peject(iepl)%xu = xu
                end if
            end if
            if (peject(iepl)%yl .ge. yu .or. peject(iepl)%yu .le. yl) then
                peject(iepl)%yl = dmiss
                peject(iepl)%yu = dmiss
            else
                if (peject(iepl)%yl .lt. yl) then
                    peject(iepl)%yl = yl
                end if
                if (peject(iepl)%yu .gt. yu) then
                    peject(iepl)%yu = yu
                end if
            end if
            if (peject(iepl)%zl .ge. zu .or. peject(iepl)%zu .le. zl) then
                peject(iepl)%zl = dmiss
                peject(iepl)%zu = dmiss
            else
                if (peject(iepl)%zl .lt. zl) then
                    peject(iepl)%zl = zl
                end if
                if (peject(iepl)%zu .gt. zu) then
                    peject(iepl)%zu = zu
                end if
            end if
        end do
    end do

!-------------------- calculate area of particle emission
    nne = 0
    do is = 1, nspec
        nns = nne + 1
        nne = nne + nepl(is)
        do iepl = nns, nne
            ipcplane = ipcpl(iepl)
            if (geotype(ipcplane) .eq. 0 .or. &
           &   geotype(ipcplane) .eq. 1) then
                peject(iepl)%xlu = (peject(iepl)%xu - peject(iepl)%xl)*dr
                peject(iepl)%ylu = (peject(iepl)%yu - peject(iepl)%yl)*dr
                peject(iepl)%zlu = (peject(iepl)%zu - peject(iepl)%zl)*dr
                if (abs(nemd(iepl)) .eq. 1 .and. peject(iepl)%grd .ne. dmiss) then
                    peject(iepl)%area = peject(iepl)%ylu*peject(iepl)%zlu
                else if (abs(nemd(iepl)) .eq. 2 .and. peject(iepl)%grd .ne. dmiss) then
                    peject(iepl)%area = peject(iepl)%zlu*peject(iepl)%xlu
                else if (abs(nemd(iepl)) .eq. 3 .and. peject(iepl)%grd .ne. dmiss) then
                    peject(iepl)%area = peject(iepl)%xlu*peject(iepl)%ylu
                else
                    peject(iepl)%area = 0.0d0
                end if
            else if (geotype(ipcplane) .eq. 2) then
                if (abs(nemd(iepl)) .eq. cylinder(ipcplane)%align) then
                    if (abs(nemd(iepl)) .eq. 1) then
                        if (peject(iepl)%grd .ne. dmiss .and. &
                       &   cylinder(ipcplane)%axis(1) .ge. yl .and. &
                       &   cylinder(ipcplane)%axis(1) .lt. yu .and. &
                       &   cylinder(ipcplane)%axis(2) .ge. zl .and. &
                       &   cylinder(ipcplane)%axis(2) .lt. zu) then
                            radius = cylinder(ipcplane)%radius
                            peject(iepl)%area = pi*radius*radius
                        else
                            peject(iepl)%area = 0.0d0
                        end if
                    else if (abs(nemd(iepl)) .eq. 2) then
                        if (peject(iepl)%grd .ne. dmiss .and. &
                       &   cylinder(ipcplane)%axis(1) .ge. zl .and. &
                       &   cylinder(ipcplane)%axis(1) .lt. zu .and. &
                       &   cylinder(ipcplane)%axis(2) .ge. xl .and. &
                       &   cylinder(ipcplane)%axis(2) .lt. xu) then
                            radius = cylinder(ipcplane)%radius
                            peject(iepl)%area = pi*radius*radius
                        else
                            peject(iepl)%area = 0.0d0
                        end if
                    else if (abs(nemd(iepl)) .eq. 3) then
                        if (peject(iepl)%grd .ne. dmiss .and. &
                       &   cylinder(ipcplane)%axis(1) .ge. xl .and. &
                       &   cylinder(ipcplane)%axis(1) .lt. xu .and. &
                       &   cylinder(ipcplane)%axis(2) .ge. yl .and. &
                       &   cylinder(ipcplane)%axis(2) .lt. yu) then
                            radius = cylinder(ipcplane)%radius
                            peject(iepl)%area = pi*radius*radius
                        else
                            peject(iepl)%area = 0.0d0
                        end if
                    end if
                else
                    if (cylinder(ipcplane)%align .eq. 1) then
                        peject(iepl)%xlu = (peject(iepl)%xu - peject(iepl)%xl)*dr
                    else
                        peject(iepl)%xlu = 2.0d0*cylinder(ipcplane)%radius*dr
                    end if
                    if (cylinder(ipcplane)%align .eq. 2) then
                        peject(iepl)%ylu = (peject(iepl)%yu - peject(iepl)%yl)*dr
                    else
                        peject(iepl)%ylu = 2.0d0*cylinder(ipcplane)%radius*dr
                    end if
                    if (cylinder(ipcplane)%align .eq. 3) then
                        peject(iepl)%zlu = (peject(iepl)%zu - peject(iepl)%zl)*dr
                    else
                        peject(iepl)%zlu = 2.0d0*cylinder(ipcplane)%radius*dr
                    end if
                    if (cylinder(ipcplane)%align .eq. 1 .and. &
                   &   cylinder(ipcplane)%axis(1) .ge. yl .and. &
                   &   cylinder(ipcplane)%axis(1) .lt. yu .and. &
                   &   cylinder(ipcplane)%axis(2) .ge. zl .and. &
                   &   cylinder(ipcplane)%axis(2) .lt. zu) then
                        peject(iepl)%area = &
                       &  peject(iepl)%xlu*2.0d0*cylinder(ipcplane)%radius*dr
                        if (omniemit(iepl) .eq. 1) &
                       &  peject(iepl)%area = peject(iepl)%area*pi
                    else if (cylinder(ipcplane)%align .eq. 2 .and. &
                   &   cylinder(ipcplane)%axis(1) .ge. zl .and. &
                   &   cylinder(ipcplane)%axis(1) .lt. zu .and. &
                   &   cylinder(ipcplane)%axis(2) .ge. xl .and. &
                   &   cylinder(ipcplane)%axis(2) .lt. xu) then
                        peject(iepl)%area = &
                       &  peject(iepl)%ylu*2.0d0*cylinder(ipcplane)%radius*dr
                        if (omniemit(iepl) .eq. 1) &
                       &  peject(iepl)%area = peject(iepl)%area*pi
                    else if (cylinder(ipcplane)%align .eq. 3 .and. &
                   &   cylinder(ipcplane)%axis(1) .ge. xl .and. &
                   &   cylinder(ipcplane)%axis(1) .lt. xu .and. &
                   &   cylinder(ipcplane)%axis(2) .ge. yl .and. &
                   &   cylinder(ipcplane)%axis(2) .lt. yu) then
                        peject(iepl)%area = &
                       &  peject(iepl)%zlu*2.0d0*cylinder(ipcplane)%radius*dr
                        if (omniemit(iepl) .eq. 1) &
                       &  peject(iepl)%area = peject(iepl)%area*pi
                    else
                        peject(iepl)%area = 0.0d0
                    end if
                end if
            else if (geotype(ipcplane) .eq. 3) then
                if (sphere(ipcplane)%center(1) .ge. xl .and. &
               &   sphere(ipcplane)%center(1) .lt. xu .and. &
               &   sphere(ipcplane)%center(2) .ge. yl .and. &
               &   sphere(ipcplane)%center(2) .lt. yu .and. &
               &   sphere(ipcplane)%center(3) .ge. zl .and. &
               &   sphere(ipcplane)%center(3) .lt. zu) then
                    radius = sphere(ipcplane)%radius
                    peject(iepl)%area = pi*radius*radius
                else
                    peject(iepl)%area = 0.0d0
                end if
            else
                peject(iepl)%area = 0.0d0
            end if
            if (peject(iepl)%area .gt. 0.0d0) then
                write (6, *) " area of emission: iepl,myid,nemd,area", &
               &           iepl, myid, nemd(iepl), peject(iepl)%area
            end if
        end do
    end do

!--------------------
    nprsum = 0
    ISL1: do is = 1, nspec
        nprsum = nprsum + npr(is)
    end do ISL1
    vxf(:) = 0.0d0
    vyf(:) = 0.0d0
    vzf(:) = 0.0d0
    vxb(:) = 0.0d0
    vyb(:) = 0.0d0
    vzb(:) = 0.0d0
    vxr(:) = 0.0d0
    vyr(:) = 0.0d0
    vzr(:) = 0.0d0

!**************************************** top of species loop 2
    nee = 0
    nne = 0
    ISL2: do is = 1, nspec
        ns = nee + 1
        npq(is) = npr(is)/nphi(is)
        ne = nee + npq(is)
        nee = nee + npr(is)
        nns = nne + 1
        nne = nne + nepl(is)
        dnsf(is) = dnsf(is)/dr/dr/dr
        dnsb(is) = dnsb(is)/dr/dr/dr
!
!
!-------------------- injection from outer boundary
        if ((nflag_emit(is) .eq. 0) .and. &
       &   (npbnd(1, is) .eq. 2 .or. npbnd(2, is) .eq. 2 .or. npbnd(3, is) .eq. 2)) then
!
!-------------------- assign velocity
! - - - - ----------- vf(v): simplified model
!         ----------- vxf(is)
            if (vdtx(is) .ge. 0.0d0) then
                call vflux(vxf(ns:ns + npq(is) - 1), npq(is), &
               &           vdtx(is), peth(is), 1, iconf)
            else
                call vflux(vxf(ns:ns + npq(is) - 1), npq(is), &
               &           -vdtx(is), peth(is), -1, iconf)
                vxf(ns:ne) = -vxf(ns:ne)
            end if
!           --------- check condition code
            if (iconf .ne. 0 .and. myid .eq. 0) then
                write (6, *) 'x:forward injection of species', is, 'is canceled', &
             &             'because iconf =', iconf
            end if
!         ----------- vxb(is)
            if (vdtx(is) .ge. 0.0d0) then
                call vflux(vxb(ns:ns + npq(is) - 1), npq(is), &
               &           vdtx(is), peth(is), -1, iconb)
            else
                call vflux(vxb(ns:ns + npq(is) - 1), npq(is), &
                           -vdtx(is), peth(is), 1, iconb)
                vxb(ns:ne) = -vxb(ns:ne)
            end if
!           --------- check condition code
            if (iconb .ne. 0 .and. myid .eq. 0) then
                write (6, *) 'x:backward injection of species', is, 'is canceled', &
             &             'because iconb =', iconb
            end if
!
!         ----------- vyf(is)
            if (vdty(is) .ge. 0.0d0) then
                call vflux(vyf(ns:ns + npq(is) - 1), npq(is), &
               &           vdty(is), peth(is), 1, iconf)
            else
                call vflux(vyf(ns:ns + npq(is) - 1), npq(is), &
               &           -vdty(is), peth(is), -1, iconf)
                vyf(ns:ne) = -vyf(ns:ne)
            end if
!           --------- check condition code
            if (iconf .ne. 0 .and. myid .eq. 0) then
                write (6, *) 'y:forward injection of species', is, 'is canceled', &
             &             'because iconf =', iconf
            end if
!         ----------- vyb(is)
            if (vdty(is) .ge. 0.0d0) then
                call vflux(vyb(ns:ns + npq(is) - 1), npq(is), &
               &           vdty(is), peth(is), -1, iconb)
            else
                call vflux(vyb(ns:ns + npq(is) - 1), npq(is), &
               &           -vdty(is), peth(is), 1, iconb)
                vyb(ns:ne) = -vyb(ns:ne)
            end if
!           --------- check condition code
            if (iconb .ne. 0 .and. myid .eq. 0) then
                write (6, *) 'y:backward injection of species', is, 'is canceled', &
             &             'because iconb =', iconb
            end if
!
!         ----------- vzf(is)
            if (vdtz(is) .ge. 0.0d0) then
                call vflux(vzf(ns:ns + npq(is) - 1), npq(is), &
               &           vdtz(is), path(is), 1, iconf)
            else
                call vflux(vzf(ns:ns + npq(is) - 1), npq(is), &
               &           -vdtz(is), path(is), -1, iconf)
                vzf(ns:ne) = -vzf(ns:ne)
            end if
!           --------- check condition code
            if (iconf .ne. 0 .and. myid .eq. 0) then
                write (6, *) 'z:forward injection of species', is, 'is canceled', &
             &             'because iconf =', iconf
            end if
!         ----------- vzb(is)
            if (vdtz(is) .ge. 0.0d0) then
                call vflux(vzb(ns:ns + npq(is) - 1), npq(is), &
               &           vdtz(is), path(is), -1, iconb)
            else
                call vflux(vzb(ns:ns + npq(is) - 1), npq(is), &
               &           -vdtz(is), path(is), 1, iconb)
                vzb(ns:ne) = -vzb(ns:ne)
            end if
!           --------- check condition code
            if (iconb .ne. 0 .and. myid .eq. 0) then
                write (6, *) 'z:backward injection of species', is, 'is canceled', &
             &             'because iconb =', iconb
            end if
!
! - - - - ----------- arrange array
!          if(f(is).ne.0.d0) then
!            call arange(vxr(ns),vyr(ns),npq(is),1)
!            call arange(vzr(ns),vyr(ns),npq(is),1)
!          end if
!            call arymod(vxr(ns),vyr(ns),npr(is),npq(is),1,nphi(is))
!            call arymod(vzr(ns),vyr(ns),npr(is),npq(is),1,nphi(is))
!
            call arange(vxf(ns), vyr(ns), npq(is), 1)
            call arange(vyf(ns), vyr(ns), npq(is), 1)
            call arange(vzf(ns), vyr(ns), npq(is), 1)
            call arymod(vxf(ns), vyr(ns), npr(is), npq(is), 1, nphi(is))
            call arymod(vyf(ns), vyr(ns), npr(is), npq(is), 1, nphi(is))
            call arymod(vzf(ns), vyr(ns), npr(is), npq(is), 1, nphi(is))
!
            call arange(vxb(ns), vyr(ns), npq(is), 1)
            call arange(vyb(ns), vyr(ns), npq(is), 1)
            call arange(vzb(ns), vyr(ns), npq(is), 1)
            call arymod(vxb(ns), vyr(ns), npr(is), npq(is), 1, nphi(is))
            call arymod(vyb(ns), vyr(ns), npr(is), npq(is), 1, nphi(is))
            call arymod(vzb(ns), vyr(ns), npr(is), npq(is), 1, nphi(is))

! - - - - ----------- vxr(is)
            call RANN0(vdtx(is), peth(is), vxr(ns:ns + npq(is) - 1), npq(is), icon)
!         ----------- vyr(is)
            call RANN0(vdty(is), peth(is), vyr(ns:ns + npq(is) - 1), npq(is), icon)
!         ----------- vzr(is)
            call RANN0(vdtz(is), path(is), vzr(ns:ns + npq(is) - 1), npq(is), icon)
!

! - - - - ----------- distribute velocity vzr(is) ==> vyr(is),vzr(is)
!          call DVRAU4(iseed1,dranu,npq(is),drawork1,nrawork1,icon)
!          if(icon.ne.0) print*,"[inipin] icon=",icon
!          rad = pi2/nphi(is)
!          iran = 0
!          do m=ns,ns+npq(is)*nphi(is),nphi(is)
!            iran = iran + 1
!            rad0 = pi2*dranu(iran)
!            vs = vxr(m)
!            do l=0,nphi(is)-1
!              vxr(m+l) = vs*dcos(rad0+rad*l)
!              vyr(m+l) = vs*dsin(rad0+rad*l)
!              vzr(m+l) = vzr(m)
!              if(inpf(is).gt.0) then
!                vxf(m+l) = vxf(m)
!                vyf(m+l) = vyf(m)
!                vzf(m+l) = vzf(m)
!              end if
!              if(inpb(is).gt.0) then
!                vxb(m+l) = vxb(m)
!                vyb(m+l) = vyb(m)
!                vzb(m+l) = vzb(m)
!              end if
!            end do
!          end do
!
            call rexch(npr(is), &
           &           d1=vxf(ns:ns + npr(is) - 1), d2=vxb(ns:ns + npr(is) - 1), &
           &           d3=vyf(ns:ns + npr(is) - 1), d4=vyb(ns:ns + npr(is) - 1), &
           &           d5=vzf(ns:ns + npr(is) - 1), d6=vzb(ns:ns + npr(is) - 1))
            call rexch(npr(is), &
           &           d1=vxr(ns:ns + npr(is) - 1), d2=vyr(ns:ns + npr(is) - 1), &
           &           d3=vzr(ns:ns + npr(is) - 1))
!
!
!-------------------- calculate/assign flux
            injct(is) = injct(is)*mltstp
            pdf = abs(npin(is))/(slx*sly*slz)
            pdb = pdf
!
            arrxf(is) = 0.d0
            arrxb(is) = 0.d0
            arryf(is) = 0.d0
            arryb(is) = 0.d0
            arrzf(is) = 0.d0
            arrzb(is) = 0.d0
!
            aa = sqrt(2.d0/pi)*peth(is) &
                &    *exp(-(vdtx(is)/peth(is))**2*0.5d0)
            bb = erf(vdtx(is)/(sqrt(2.d0)*peth(is)))
!         ----------- forward: arrxf(is)
            arrxf(is) = (aa + vdtx(is)*(1.d0 + bb))*0.5d0
            if (myid .eq. 0) print *, "abs(arrxf(is))*pdf", abs(arrxf(is))*pdf
            arrxf(is) = abs(arrxf(is))*pdf*(yu - yl)*(zu - zl)*dt*0.5d0
            if (inpf(is) .le. 0) arrxf(is) = 0.0d0
!         ----------- backward: arrxb(is)
            arrxb(is) = (-aa + vdtx(is)*(1.d0 - bb))*0.5d0
            if (myid .eq. 0) print *, "abs(arrxb(is))*pdb", abs(arrxb(is))*pdb
            arrxb(is) = abs(arrxb(is))*pdb*(yu - yl)*(zu - zl)*dt*0.5d0
            if (inpb(is) .le. 0) arrxb(is) = 0.0d0
!
            aa = sqrt(2.d0/pi)*peth(is) &
                &    *exp(-(vdty(is)/peth(is))**2*0.5d0)
            bb = erf(vdty(is)/(sqrt(2.d0)*peth(is)))
!         ----------- forward: arryf(is)
            arryf(is) = (aa + vdty(is)*(1.d0 + bb))*0.5d0
            if (myid .eq. 0) print *, "abs(arryf(is))*pdf", abs(arryf(is))*pdf
            arryf(is) = abs(arryf(is))*pdf*(zu - zl)*(xu - xl)*dt*0.5d0
            if (inpf(is) .le. 0) arryf(is) = 0.0d0
!         ----------- backward: arryb(is)
            arryb(is) = (-aa + vdty(is)*(1.d0 - bb))*0.5d0
            if (myid .eq. 0) print *, "abs(arryb(is))*pdb", abs(arryb(is))*pdb
            arryb(is) = abs(arryb(is))*pdb*(zu - zl)*(xu - xl)*dt*0.5d0
            if (inpb(is) .le. 0) arryb(is) = 0.0d0
!
            aa = sqrt(2.d0/pi)*path(is) &
                &    *exp(-(vdtz(is)/path(is))**2*0.5d0)
            bb = erf(vdtz(is)/(sqrt(2.d0)*path(is)))
!         ----------- forward: arrzf(is)
            arrzf(is) = (aa + vdtz(is)*(1.d0 + bb))*0.5d0
            if (myid .eq. 0) print *, "abs(arrzf(is))*pdf", abs(arrzf(is))*pdf
            arrzf(is) = abs(arrzf(is))*pdf*(xu - xl)*(yu - yl)*dt*0.5d0
            if (inpf(is) .le. 0) arrzf(is) = 0.0d0
!         ----------- backward: arrzb(is)
            arrzb(is) = (-aa + vdtz(is)*(1.d0 - bb))*0.5d0
            if (myid .eq. 0) print *, "abs(arrzb(is))*pdb", abs(arrzb(is))*pdb
            arrzb(is) = abs(arrzb(is))*pdb*(xu - xl)*(yu - yl)*dt*0.5d0
            if (inpb(is) .le. 0) arrzb(is) = 0.0d0
!
            if (myid .eq. 0) &
           &  write (6, *) 'arrxf(is) =', arrxf(is), 'arrxb(is) =', arrxb(is)
            if (myid .eq. 0) &
           &  write (6, *) 'arryf(is) =', arryf(is), 'arryb(is) =', arryb(is)
            if (myid .eq. 0) &
           &  write (6, *) 'arrzf(is) =', arrzf(is), 'arrzb(is) =', arrzb(is)
!
            if (myid .eq. 0) &
           &  print *, "inipin set for injection from outer boundary: is", is
            cycle
!
!
!-------------------- injection from inner boundary
        else if (nflag_emit(is) .eq. 1) then

!-------------------- assign velocity
            call assign_velocity

!-------------------- calculate/assign flux
            injct(is) = injct(is)*mltstp
            if (curf(is) .gt. 0.0d0 .and. q(is) .ne. 0.0d0) then
                flpf(is) = curf(is)/abs(q(is))
                if (myid .eq. 0) &
               &  print *, "is, flpf from curf: ", is, flpf(is), curf(is)
            end if
!          if(curb(is).gt.0.0d0.and.q(is).ne.0.0d0) then
!            flpb(is) = curb(is)/abs(q(is))
!            if(myid.eq.0) &
!           &  print*,"is, flpb from curb: ",is,flpb(is),curb(is)
!          end if
            if (dnsf(is) .gt. 0.0d0 .or. dnsb(is) .gt. 0.0d0) then
                pdf = dnsf(is)/(dr*dr*dr)
                pdb = dnsb(is)/(dr*dr*dr)
            else
                pdf = npin(is)/(slx*sly*slz)
                pdb = pdf
            end if
!
            do iepl = nns, nne
                fluxf(iepl) = 0.0d0
                fluxb(iepl) = 0.0d0
                ipcplane = ipcpl(iepl)
!
!           --------- forward
                if (curfs(iepl) .gt. 0.0d0 .and. q(is) .ne. 0.0d0) then
                    flpfs(iepl) = curfs(iepl)/abs(q(is))
                    fluxf(iepl) = flpfs(iepl)*dt*0.5d0*peject(iepl)%area
                    if (myid .eq. 0) &
                   &  print *, "iepl, flpfs from curfs: ", iepl, flpfs(iepl), curfs(iepl)
                else if (flpfs(iepl) .gt. 0.0d0) then
                    fluxf(iepl) = flpfs(iepl)*dt*0.5d0*peject(iepl)%area
                else if (flpf(is) .gt. 0.0d0) then
                    fluxf(iepl) = flpf(is)*dt*0.5d0*peject(iepl)%area
                else
                    aa = sqrt(2.d0/pi)*path(is) &
                        &    *exp(-((spa(is) + abvdem(is))/path(is))**2*0.5d0)
                    bb = erf((spa(is) + abvdem(is))/(sqrt(2.d0)*path(is)))
                    fluxf(iepl) = (aa + (spa(is) + abvdem(is))*(1.d0 + bb))*0.5d0
!              if(myid.eq.0) &
!             &  print*,"iepl, particle flux: ",iepl,abs(fluxf(iepl))*pdf
                    fluxf(iepl) = abs(fluxf(iepl))*dt*0.5d0*peject(iepl)%area*pdf
                end if
                if (inpf(is) .le. 0) fluxf(iepl) = 0.0d0
!
!           --------- backward
!            if(flpb(is).gt.0.0d0) then
!              fluxb(iepl) = flpb(is)*dt*0.5d0*peject(iepl)%area
!            else
!              aa = sqrt(2.d0/pi)*path(is) &
!             &    *exp(-((spa(is) - abvdem(is))/path(is))**2*0.5d0)
!              bb = erf((spa(is) - abvdem(is))/(sqrt(2.d0)*path(is)))
!              fluxb(iepl) = (-aa + (spa(is) - abvdem(is))*(1.d0 - bb))*0.5d0
!              if(myid.eq.0) &
!             &  print*,"iepl, particle flux: ",abs(fluxb(iepl))*dnsb(is)
!              fluxb(iepl) = abs(fluxb(iepl))*dt*0.5d0*peject(iepl)%area*pdb
!            endif
!            if(inpb(is).le.0) fluxb(iepl) = 0.0d0
!
!           ---------
                if (nemd(iepl) .eq. 0) then
                    fluxf(iepl) = 0.0d0
                end if
!
                if (fluxf(iepl) .gt. 0.0d0 .and. abs(fluxf(iepl)) .gt. 1.0d-12) then
                    write (6, *) ' plane', iepl, ', myid', myid, &
                   &           ': fluxf =', fluxf(iepl)
                else
                    fluxf(iepl) = 0.0d0
                end if
!
            end do
!
!
!-------------------- photoelectron emission
        else if (nflag_emit(is) .eq. 2) then

!-------------------- assign velocity
            call assign_velocity

!-------------------- calculate/assign flux
            injct(is) = injct(is)*mltstp
            if (curf(is) .gt. 0.0d0 .and. q(is) .ne. 0.0d0) then
                flpf(is) = curf(is)/abs(q(is))
                if (myid .eq. 0) &
               &  print *, "is, flpf from curf: ", is, flpf(is), curf(is)
            end if
!          if(curb(is).gt.0.0d0.and.q(is).ne.0.0d0) then
!            flpb(is) = curb(is)/abs(q(is))
!            if(myid.eq.0) &
!           &  print*,"is, flpb from curb: ",is,flpb(is),curb(is)
!          end if
            if (dnsf(is) .gt. 0.0d0 .or. dnsb(is) .gt. 0.0d0) then
                pdf = dnsf(is)/(dr*dr*dr)
                pdb = dnsb(is)/(dr*dr*dr)
            else
                pdf = npin(is)/(slx*sly*slz)
                pdb = pdf
            end if
!
            do iepl = nns, nne
                fluxf(iepl) = 0.0d0
                fluxb(iepl) = 0.0d0
                ipcplane = ipcpl(iepl)
!
!           --------- forward
                if (curfs(iepl) .gt. 0.0d0 .and. q(is) .ne. 0.0d0) then
                    flpfs(iepl) = curfs(iepl)/abs(q(is))
                    fluxf(iepl) = flpfs(iepl)*dt*0.5d0*peject(iepl)%area
                    if (myid .eq. 0) &
                   &  print *, "iepl, flpfs from curfs: ", iepl, flpfs(iepl), curfs(iepl)
                else if (flpfs(iepl) .gt. 0.0d0) then
                    fluxf(iepl) = flpfs(iepl)*dt*0.5d0*peject(iepl)%area
                else if (flpf(is) .gt. 0.0d0) then
                    fluxf(iepl) = flpf(is)*dt*0.5d0*peject(iepl)%area
                else
                    aa = sqrt(2.d0/pi)*path(is) &
                        &    *exp(-((spa(is) + abvdem(is))/path(is))**2*0.5d0)
                    bb = erf((spa(is) + abvdem(is))/(sqrt(2.d0)*path(is)))
                    fluxf(iepl) = (aa + (spa(is) + abvdem(is))*(1.d0 + bb))*0.5d0
!              if(myid.eq.0) &
!             &  print*,"iepl, particle flux: ",iepl,abs(fluxf(iepl))*pdf
                    fluxf(iepl) = abs(fluxf(iepl))*dt*0.5d0*peject(iepl)%area*pdf
                end if
                if (inpf(is) .le. 0) fluxf(iepl) = 0.0d0
!
!           --------- backward
!            if(flpb(is).gt.0.0d0) then
!              fluxb(iepl) = flpb(is)*dt*0.5d0*peject(iepl)%area
!            else
!              aa = sqrt(2.d0/pi)*path(is) &
!             &    *exp(-((spa(is) - abvdem(is))/path(is))**2*0.5d0)
!              bb = erf((spa(is) - abvdem(is))/(sqrt(2.d0)*path(is)))
!              fluxb(iepl) = (-aa + (spa(is) - abvdem(is))*(1.d0 - bb))*0.5d0
!              if(myid.eq.0) &
!             &  print*,"iepl, particle flux: ",abs(fluxb(iepl))*dnsb(is)
!              fluxb(iepl) = abs(fluxb(iepl))*dt*0.5d0*peject(iepl)%area*pdb
!            endif
!            if(inpb(is).le.0) fluxb(iepl) = 0.0d0
!
!           ---------
                if (nemd(iepl) .eq. 0) then
                    fluxf(iepl) = 0.0d0
                end if
!
!           ---- decrease according to incident angle
                if (geotype(ipcplane) .eq. 0 .or. geotype(ipcplane) .eq. 1) then
                    if (abs(nemd(iepl)) .eq. 1) then
                        fluxf(iepl) = &
                       &  sign(1.0d0, dble(nemd(iepl)))*fluxf(iepl)*dcos(thetax)
                    else if (abs(nemd(iepl)) .eq. 2) then
                        fluxf(iepl) = &
                       &  sign(1.0d0, dble(nemd(iepl)))*fluxf(iepl)*dcos(thetay)
                    else if (abs(nemd(iepl)) .eq. 3) then
                        fluxf(iepl) = &
                       &  sign(1.0d0, dble(nemd(iepl)))*fluxf(iepl)*dcos(thetaz)
                    end if
                else if (geotype(ipcplane) .eq. 2) then
                    align = cylinder(ipcplane)%align
                    if (abs(nemd(iepl)) .eq. align) then
                        if (abs(nemd(iepl)) .eq. 1) then
                            fluxf(iepl) = &
                           &  sign(1.0d0, dble(nemd(iepl)))*fluxf(iepl)*dcos(thetax)
                        else if (abs(nemd(iepl)) .eq. 2) then
                            fluxf(iepl) = &
                           &  sign(1.0d0, dble(nemd(iepl)))*fluxf(iepl)*dcos(thetay)
                        else if (abs(nemd(iepl)) .eq. 3) then
                            fluxf(iepl) = &
                           &  sign(1.0d0, dble(nemd(iepl)))*fluxf(iepl)*dcos(thetaz)
                        end if
                    else
                        if (align .eq. 1) then
                            fluxf(iepl) = fluxf(iepl)*dsin(thetax)
                        else if (align .eq. 2) then
                            fluxf(iepl) = fluxf(iepl)*dsin(thetay)
                        else if (align .eq. 3) then
                            fluxf(iepl) = fluxf(iepl)*dsin(thetaz)
                        end if
                    end if
                end if
!
                if (fluxf(iepl) .gt. 0.0d0 .and. abs(fluxf(iepl)) .gt. 1.0d-12) then
                    write (6, *) ' plane', iepl, ', myid', myid, &
                   &           ': fluxf =', fluxf(iepl)
                else
                    fluxf(iepl) = 0.0d0
                end if
!
            end do
!
!         ------ from solid surface
            if (zl .le. zssurf .and. zu .gt. zssurf) then
                if (flpf(is) .gt. 0.0d0) then
                    fluxub = flpf(is)*(xu - xl)*(yu - yl)*dt*0.5d0
                    fluxub = fluxub*cos(thetaz)
                else
                    aa = sqrt(2.d0/pi)*path(is) &
                        &    *exp(-((spa(is) + abvdem(is))/path(is))**2*0.5d0)
                    bb = erf((spa(is) + abvdem(is))/(sqrt(2.d0)*path(is)))
                    fluxub = (aa + (spa(is) + abvdem(is))*(1.d0 + bb))*0.5d0
!              if(myid.eq.0) &
!             &  print*,"iepl, particle flux: ",iepl,abs(fluxub)*pdf
                    fluxub = abs(fluxub)*(xu - xl)*(yu - yl)*dt*0.5d0*pdf
                    fluxub = fluxub*cos(thetaz)
                end if
                if (inpf(is) .le. 0) then
                    fluxub = 0.0d0
                end if
            else
                fluxub = 0.0d0
            end if

            if (fluxub .gt. 0.0d0 .and. abs(fluxub) .gt. 1.0d-12) then
                write (6, *) ' zssurf myid', myid, &
               &           ': fluxub =', fluxub
            else
                fluxub = 0.0d0
            end if
!
!-------------------- secondary electron emission
        else if (nflag_emit(is) .eq. 3) then

!-------------------- assign velocity
            call assign_velocity

!-------------------- calculate/assign flux
            injct(is) = injct(is)*mltstp
!
            do iepl = nns, nne
                fluxf(iepl) = 0.0d0
                fluxb(iepl) = 0.0d0
                ipcplane = ipcpl(iepl)
!
                if (peject(iepl)%area .gt. 0.0d0) then
                    write (6, *) ' plane', iepl, ', myid', myid, &
                        ': area =', peject(iepl)%area, &
                        ': tarea =', peject(iepl)%tarea
                end if
!
            end do

!--------------------
        end if

    end do ISL2
!**************************************** end of species loop 2

!-------------------- calculate flux for virtual impinging particles
    call MPI_Barrier(MCW, ierr)
    ISL3: do is = 1, nspec
        do ipc = 1, npc
            fluxv(is, ipc) = 0.0d0
            aa = sqrt(2.d0/pi)*path(is) &
                &      *exp(-((spa(is) + abvdem(is))/path(is))**2*0.5d0)
            bb = erf((spa(is) + abvdem(is))/(sqrt(2.d0)*path(is)))
            if (nflag_emit(is) .lt. 0) then
                fluxv(is, ipc) = (aa + (spa(is) + abvdem(is))*(1.d0 + bb))*0.5d0
                fluxv(is, ipc) = abs(fluxv(is, ipc))*(wp(is)**2)/qm(is) &
                                &                   *imarea(ipc)*dt*0.5d0
                if (myid .eq. 0) write (6, *) ' is,ipc', is, ipc, ': fluxv =', fluxv(is, ipc)
            end if
        end do
    end do ISL3

    return

contains

    subroutine assign_velocity
!
!-------------------- vnormf
        if (path(nspec + 1) .gt. 0.0d0) then
            nprfrac1 = nint(npr(is)*mainprop)
            nprfrac2 = npr(is) - nprfrac1
            if (inpf(is) .gt. 0) call vflux(vnormf(ns:ns + nprfrac1 - 1), &
           &                             nprfrac1, spa(is) + abvdem(is), &
           &                             path(is), 1, iconf)
            if (iconf .ne. 0) then
                inpf(is) = 0
                write (6, *) 'forward injection of species', is, 'is canceled', &
             &             'because iconf =', iconf
            end if
            if (nprfrac2 .gt. 0) then
                if (inpf(is) .gt. 0) call vflux(vnormf(ns + nprfrac1: &
               &                                    ns + nprfrac1 + nprfrac2 - 1), &
               &                             nprfrac2, spa(is) + abvdem(is), &
               &                             path(nspec + 1), 1, iconf)
                if (iconf .ne. 0) then
                    inpf(is) = 0
                    write (6, *) 'forward injection of species', is, 'is canceled', &
                 &             'because iconf =', iconf
                end if
            end if
        else
            if (inpf(is) .gt. 0) call vflux(vnormf(ns:ns + npq(is) - 1), &
           &                             npq(is), spa(is) + abvdem(is), &
           &                             path(is), 1, iconf)
            if (iconf .ne. 0) then
                inpf(is) = 0
                write (6, *) 'forward injection of species', is, 'is canceled', &
             &             'because iconf =', iconf
            end if
        end if

! - - - - ----------- random exchange
        call rexch(npr(is), d1=vnormf(ns:ns + npr(is) - 1))
        call rexch(npr(is), d1=vnormf(ns:ns + npr(is) - 1))

!-------------------- vtangf
        if (peth(nspec + 1) .gt. 0.0d0) then
            nprfrac1 = nint(npr(is)*mainprop)
            nprfrac2 = npr(is) - nprfrac1
            call RANN0(spe(is), peth(is), vtangf(ns:ns + nprfrac1 - 1), &
           &           nprfrac1, icon)
            if (nprfrac2 .gt. 0) then
                call RANN0(spe(is), peth(nspec + 1), &
               &           vtangf(ns + nprfrac1:ns + nprfrac1 + nprfrac2 - 1), &
               &           nprfrac2, icon)
            end if
        else
            call RANN0(spe(is), peth(is), vtangf(ns:ns + npq(is) - 1), &
           &           npq(is), icon)
        end if

!-------------------- runiform
        m = ns
        do while (m .lt. ns + npr(is))
            call RANU0(dranu(1:2), 2, icon)
            if (dranu(2) .le. dranu(1)) then
                runiform(m) = dranu(1)
                m = m + 1
            end if
        end do

!-------------------- vnormc,vtangc,psicos
        m = ns
        do while (m .lt. ns + npr(is))
            call RANU0(dranu(1:2), 2, icon)
            psi1 = pi*(dranu(1) - 0.5d0)
            if (dranu(2) .le. dcos(psi1)) then
                vtangc(m) = vnormf(m)*dsin(psi1)
                vnormc(m) = vnormf(m)*dcos(psi1)
                m = m + 1
            end if
        end do
        m = ns
        do while (m .lt. ns + npr(is))
            call RANU0(dranu(1:2), 2, icon)
            psi1 = pi*(dranu(1) - 0.5d0)
            if (dranu(2) .le. dcos(psi1)) then
                psicos(m) = psi1
                m = m + 1
            end if
        end do

!         ----------- vnorms, vtangs1, vtangs2
        m = ns
        do while (m .lt. ns + npr(is))
            call RANU0(dranu(1:3), 3, icon)
            psi1 = pi*(dranu(1) - 0.5d0)
            if (dranu(2) .le. dcos(psi1)*dsin(psi1)*2.0d0) then
                vtangs1(m) = vxf(m)*dsin(psi1)
                vnorms(m) = vxf(m)*dcos(psi1)
                psi2 = 2.0d0*pi*dranu(3)
                vtangs2(m) = vtangs1(m)*dsin(psi2)
                vtangs1(m) = vtangs1(m)*dcos(psi2)
                m = m + 1
            end if
        end do

        return
    end subroutine

end subroutine
