#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine positn(ps, ustep)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   P O S I T N
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .     this subroutine updates particle position.           .
!   .     the leap-frog scheme is used.                        .
!   ............................................................
!
!
!-------------------- parameter and common block
      use oh_type
      use paramt
      use allcom
      ! use detect_collision
      use m_boundary_base
      use particle_collision
      implicit none
!
      integer(kind=8) :: m, mm, ns, ne, nee
      integer(kind=8) :: nprs, npre
      integer(kind=8) :: isfluxt(inpc, ispec)
      integer(kind=8) :: influxsum
      integer(kind=4) :: i, j, k
      integer(kind=4) :: i1, j1, k1
      integer(kind=4) :: is, ustep, ibdy, ipc
      integer(kind=4) :: xl, xu, yl, yu, zl, zu
      integer(kind=4) :: ngx, ngy, ngz
      integer(kind=4) :: nb(3)
      integer(kind=4) :: ps, rid, nud
      integer(kind=4) :: nyield
      integer(kind=4) :: icon, iran
      integer(kind=4) :: addr0, addr1, addr2
!  integer(kind=4) :: oh3_map_region_to_node
      real(kind=8) :: dxl, dxu, dyl, dyu, dzl, dzu
      real(kind=8) :: xlocal, ylocal, zlocal
      real(kind=8) :: x, y, z
      real(kind=8) :: x1, y1, z1, z2, xy1, xz1, yz1, xz2, yz2
      real(kind=8) :: tslx, tsly, tslz, rustep
      real(kind=8) :: v1, v2, v3, v4, v5, v6, v7, v8
      real(kind=8) :: eex, eey, eez, bbx, bby, bbz, boris, vxt, vyt, vzt
      real(kind=8) :: disp1, disp2, disp3
      real(kind=8) :: radsq(inpc)
      real(kind=8) :: vxx, vyy, vzz, vxy, vyx, vyz, vzy, vzx, vxz
      real(kind=8) :: vxymyx, vyzmzy, vzxmxz
      real(kind=8) :: discrim, v2norm, v2normi, paramt1, paramt2, llmt, ulmt
      real(kind=8) :: perg, pergfrac, pemaxinvh(ispec), costhi, tfrac
      real(kind=8) :: xpast0, ypast0, zpast0, xpast1, ypast1, zpast1
      real(kind=8) :: yield, weightr, delta
!
      real(kind=8) :: xmove, ymove, zmove, xsepa, ysepa, termA, termB
      real(kind=8) :: tfrac1, tfrac2
      real(kind=8) :: xpast1a, ypast1a, xpast1b, ypast1b
      real(kind=8) :: xpast2a, ypast2a, xpast2b, ypast2b
      real(kind=8) :: cosbetav, sinbetav
      integer(kind=4) :: chkflg

      real(kind=8) :: gustep

      type(t_CollisionRecord) :: record
      type(t_BoundaryList) :: boundaries

!--------------------
      xl = sdoms(1, 1, sdid(ps) + 1); xu = sdoms(2, 1, sdid(ps) + 1)
      yl = sdoms(1, 2, sdid(ps) + 1); yu = sdoms(2, 2, sdid(ps) + 1)
      zl = sdoms(1, 3, sdid(ps) + 1); zu = sdoms(2, 3, sdid(ps) + 1)
      ngx = xu - xl; ngy = yu - yl; ngz = zu - zl
      dxl = xl; dxu = xu; dyl = yl; dyu = yu; dzl = zl; dzu = zu

      ! call collision_set_local_range(dxl-1, dxu+1, dyl-1, dyu+1, dzl-1, dzu+1)
      ! call collision_clear_all_faces
      ! call collision_add_ssurf3
      boundaries = create_simple_collision_boundaries(sdoms(1:2, 1:3, sdid(ps) + 1))

!--------------------
      pemaxinvh(1:nspec) = 0.5d0/pemax(1:nspec)

!-------------------- for three dimensional system
!     --------------- zero clear rho
!      rho(:,:,:,:,ps) = 0.0d0

!--------------------
      do ipc = 1, npc
          if (geotype(ipc) .eq. 2) then
              radsq(ipc) = cylinder(ipc)%radius*cylinder(ipc)%radius
          else if (geotype(ipc) .eq. 3) then
              radsq(ipc) = sphere(ipc)%radius*sphere(ipc)%radius
          end if
      end do

!============================== species loop
      nee = pbase(ps)
      npre = 0
      do is = 1, isse
          nprs = npre + 1
          npre = npre + npr(is)
      end do
      rustep = 1.0d0/ustep
      tslx = 2.0d0*slx
      tsly = 2.0d0*sly
      tslz = 2.0d0*slz
      ISL1: do is = 1, nspec
          ns = nee + 1
          ne = nee + totalp(is, ps)
          nee = ne
!       -------------
          weightr = q(is)/q(isse)
          delta = weightr*deltaemax(is)
!       ------------- inner loop
          do m = ns, ne
            if (pbuf(m)%preside < 0) then
                nphgram(pbuf(m)%nid + 1, is, ps) = &
                    nphgram(pbuf(m)%nid + 1, is, ps) - 1
                pbuf(m)%nid = -1
                pbuf(m)%preside = 0
            else if (pbuf(m)%preside == OH_PCL_INJECTED) then
                call RANU0(dranu, 1, icon)
                gustep = ustep*dranu(1)
                pbuf(m)%preside = 0
            else ! if pbuf(m)%preside == OH_PCL_ALIVE
                gustep = ustep
            end if
!
              if (pbuf(m)%nid .eq. -1) cycle
!
!         ----------- update particle positions
              xmove = pbuf(m)%vx*ustep
              ymove = pbuf(m)%vy*ustep
              zmove = pbuf(m)%vz*ustep
              pbuf(m)%x = pbuf(m)%x + xmove
              pbuf(m)%y = pbuf(m)%y + ymove
              pbuf(m)%z = pbuf(m)%z + zmove

!         ----------- internal boundary treatment 0
              gustep = ustep
              block
                double precision :: p1(3), p2(3), v(3)
                integer :: iboundary
                type(t_CollisionRecord) :: tmp_record

                v(:) = [pbuf(m)%vx, pbuf(m)%vy, pbuf(m)%vz]
                p2(:) = [pbuf(m)%x, pbuf(m)%y, pbuf(m)%z]
                p1(:) = p2(:) - v(:)*gustep
                
                record%is_collided = .false.
                record%t = 100.0d0  ! "0 < t < 1" should be established when collided.
                do iboundary = 1, boundaries%nboundaries
                    tmp_record = boundaries%boundaries(iboundary)%check_collision(p1, p2)
                    if (tmp_record%is_collided .and. tmp_record%t < record%t) then
                        record = tmp_record
                    end if
                end do

                ! record = boundaries%check_collision(p1(:), p2(:))
              end block
              if (record%is_collided) then
                  pbuf(m)%x = record%position(1)
                  pbuf(m)%y = record%position(2)
                  pbuf(m)%z = record%position(3)
                  pbuf(m)%preside = OH_PCL_TO_BE_ACCUMULATED
              end if
! #include "defsurf.fnc"

!         ----------- external boundary treatment
              nud = 14
              if (pbuf(m)%x .lt. dxl) then
                  nud = nud - 1
              else if (pbuf(m)%x .ge. dxu) then
                  nud = nud + 1
              end if
              if (pbuf(m)%y .lt. dyl) then
                  nud = nud - 3
              else if (pbuf(m)%y .ge. dyu) then
                  nud = nud + 3
              end if
              if (pbuf(m)%z .lt. dzl) then
                  nud = nud - 9
              else if (pbuf(m)%z .ge. dzu) then
                  nud = nud + 9
              end if
              if (nud .ne. 14) then
                  rid = nborps(nud, is, ps)
                  pbuf(m)%nid = rid
                  nphgram(sdid(ps) + 1, is, ps) = nphgram(sdid(ps) + 1, is, ps) - 1
                  if (rid .ne. -1) then
                      nphgram(rid + 1, is, ps) = nphgram(rid + 1, is, ps) + 1
                  end if
                  if (pbuf(m)%x .lt. 0.0d0) then
                      if (npbnd(1, is) .eq. 0) then
                          pbuf(m)%x = pbuf(m)%x + slx
                      else if (npbnd(1, is) .eq. 1) then
                          pbuf(m)%x = -pbuf(m)%x
                          pbuf(m)%vx = -pbuf(m)%vx
                      else
                          gcount(1)%nesc(is) = gcount(1)%nesc(is) + 1
                      end if
                  else if (pbuf(m)%x .ge. slx) then
                      if (npbnd(1, is) .eq. 0) then
                          pbuf(m)%x = pbuf(m)%x - slx
                      else if (npbnd(1, is) .eq. 1) then
                          pbuf(m)%x = tslx - pbuf(m)%x
                          pbuf(m)%vx = -pbuf(m)%vx
                      else
                          gcount(1)%nesc(is) = gcount(1)%nesc(is) + 1
                      end if
                  end if
!
                  if (pbuf(m)%y .lt. 0.0d0) then
                      if (npbnd(2, is) .eq. 0) then
                          pbuf(m)%y = pbuf(m)%y + sly
                      else if (npbnd(2, is) .eq. 1) then
                          pbuf(m)%y = -pbuf(m)%y
                          pbuf(m)%vy = -pbuf(m)%vy
                      else
                          gcount(1)%nesc(is) = gcount(1)%nesc(is) + 1
                      end if
                  else if (pbuf(m)%y .ge. sly) then
                      if (npbnd(2, is) .eq. 0) then
                          pbuf(m)%y = pbuf(m)%y - sly
                      else if (npbnd(2, is) .eq. 1) then
                          pbuf(m)%y = tsly - pbuf(m)%y
                          pbuf(m)%vy = -pbuf(m)%vy
                      else
                          gcount(1)%nesc(is) = gcount(1)%nesc(is) + 1
                      end if
                  end if
!
                  if (pbuf(m)%z .lt. 0.0d0) then
                      if (npbnd(3, is) .eq. 0) then
                          pbuf(m)%z = pbuf(m)%z + slz
                      else if (npbnd(3, is) .eq. 1) then
                          pbuf(m)%z = -pbuf(m)%z
                          pbuf(m)%vz = -pbuf(m)%vz
                      else
                          gcount(1)%nesc(is) = gcount(1)%nesc(is) + 1
                      end if
                  else if (pbuf(m)%z .ge. slz) then
                      if (npbnd(3, is) .eq. 0) then
                          pbuf(m)%z = pbuf(m)%z - slz
                      else if (npbnd(3, is) .eq. 1) then
                          pbuf(m)%z = tslz - pbuf(m)%z
                          pbuf(m)%vz = -pbuf(m)%vz
                      else
                          gcount(1)%nesc(is) = gcount(1)%nesc(is) + 1
                      end if
                  end if
              end if

!         ----------- internal boundary treatment 1
!          if(line_mode.ne.1.and.pbuf(m)%nid.ne.-1) then
              if (line_mode .ne. 1) then
                  do ipc = 1, npc
                      if (pbuf(m)%nid .eq. -1) cycle
                      if ((geotype(ipc) .eq. 0 .or. geotype(ipc) .eq. 1) .and. &
                     &   (pbuf(m)%x .ge. xlpc(ipc) .and. pbuf(m)%x .le. xupc(ipc) .and. &
                     &    pbuf(m)%y .ge. ylpc(ipc) .and. pbuf(m)%y .le. yupc(ipc) .and. &
                     &    pbuf(m)%z .ge. zlpc(ipc) .and. pbuf(m)%z .le. zupc(ipc))) then
                          nphgram(pbuf(m)%nid + 1, is, ps) = nphgram(pbuf(m)%nid + 1, is, ps) - 1
                          pbuf(m)%nid = -1
                          gcount(1)%influx(ipc, is) = gcount(1)%influx(ipc, is) + 1
                          gcount(1)%infhist(pbuf(m)%pid, ipc, is) = &
                         &  gcount(1)%infhist(pbuf(m)%pid, ipc, is) + 1
!
                      else if (geotype(ipc) .eq. 2) then
                          if (cylinder(ipc)%align .eq. 1) then
                              disp1 = pbuf(m)%y - cylinder(ipc)%axis(1)
                              disp2 = pbuf(m)%z - cylinder(ipc)%axis(2)
                              if (pbuf(m)%x .ge. xlpc(ipc) .and. pbuf(m)%x .le. xupc(ipc) .and. &
                             &   disp1*disp1 + disp2*disp2 .lt. radsq(ipc)) then
                                  nphgram(pbuf(m)%nid + 1, is, ps) = &
                                 &  nphgram(pbuf(m)%nid + 1, is, ps) - 1
                                  pbuf(m)%nid = -1
                                  gcount(1)%influx(ipc, is) = gcount(1)%influx(ipc, is) + 1
                                  gcount(1)%infhist(pbuf(m)%pid, ipc, is) = &
                                 &  gcount(1)%infhist(pbuf(m)%pid, ipc, is) + 1
                              end if
                          else if (cylinder(ipc)%align .eq. 2) then
                              disp1 = pbuf(m)%z - cylinder(ipc)%axis(1)
                              disp2 = pbuf(m)%x - cylinder(ipc)%axis(2)
                              if (pbuf(m)%y .ge. ylpc(ipc) .and. pbuf(m)%y .le. yupc(ipc) .and. &
                             &   disp1*disp1 + disp2*disp2 .lt. radsq(ipc)) then
                                  nphgram(pbuf(m)%nid + 1, is, ps) = &
                                &  nphgram(pbuf(m)%nid + 1, is, ps) - 1
                                  pbuf(m)%nid = -1
                                  gcount(1)%influx(ipc, is) = gcount(1)%influx(ipc, is) + 1
                                  gcount(1)%infhist(pbuf(m)%pid, ipc, is) = &
                                 &  gcount(1)%infhist(pbuf(m)%pid, ipc, is) + 1
                              end if
                          else if (cylinder(ipc)%align .eq. 3) then
                              disp1 = pbuf(m)%x - cylinder(ipc)%axis(1)
                              disp2 = pbuf(m)%y - cylinder(ipc)%axis(2)
                              if (pbuf(m)%z .ge. zlpc(ipc) .and. pbuf(m)%z .le. zupc(ipc) .and. &
                             &   disp1*disp1 + disp2*disp2 .lt. radsq(ipc)) then
                                  nphgram(pbuf(m)%nid + 1, is, ps) = &
                                 &  nphgram(pbuf(m)%nid + 1, is, ps) - 1
                                  pbuf(m)%nid = -1
                                  gcount(1)%influx(ipc, is) = gcount(1)%influx(ipc, is) + 1
                                  gcount(1)%infhist(pbuf(m)%pid, ipc, is) = &
                                 &  gcount(1)%infhist(pbuf(m)%pid, ipc, is) + 1
                              end if
                          end if
                      else if (geotype(ipc) .eq. 3) then
                          disp1 = pbuf(m)%x - sphere(ipc)%center(1)
                          disp2 = pbuf(m)%y - sphere(ipc)%center(2)
                          disp3 = pbuf(m)%z - sphere(ipc)%center(3)
                          if (disp1*disp1 + disp2*disp2 + disp3*disp3 .lt. radsq(ipc)) then
                              nphgram(pbuf(m)%nid + 1, is, ps) = &
                             &  nphgram(pbuf(m)%nid + 1, is, ps) - 1
                              pbuf(m)%nid = -1
                              gcount(1)%influx(ipc, is) = gcount(1)%influx(ipc, is) + 1
                              gcount(1)%infhist(pbuf(m)%pid, ipc, is) = &
                             &  gcount(1)%infhist(pbuf(m)%pid, ipc, is) + 1
                          end if
                      end if
                  end do
              end if
          end do

      end do ISL1
!============================== species loop end

!-------------------- diagnostics
!    if(myid.eq.0) then
!      if(intfoc.ne.0) then
!        influxsum = 0
!        do ipc=1,npc
!          influx(ipc) = 0
!          do is=1,nspec
!            influx(ipc) = influx(ipc) + isfluxt(ipc,is)
!          end do
!          influxsum = influxsum + influx(ipc)
!        end do
!        if(istep.eq.0.or.mod(istep-1,intfoc).eq.0) then
!          open(90,file='influx',position='append')
!          open(91,file='isflux',position='append')
!          open(92,file='nesc',position='append')
!        end if
!        write(90,*) t, (influx(ipc),ipc=1,npc), influxsum
!        write(91,*) t, ((isflux(ipc,is),ipc=1,npc),is=1,nspec)
!        write(92,*) t, (nesc(is),is=1,nspec)
!        if(mod(istep,intfoc).eq.0.or.istep.eq.nstep) then
!          close(90)
!          close(91)
!          close(92)
!        end if
!      end if
!    end if

      call boundaries%destroy
      return
end subroutine
