#ifndef ssurf
    if (pbuf(m)%z .le. zssurf) then
        tfrac = (pbuf(m)%z - zssurf)/pbuf(m)%vz
        pbuf(m)%x = pbuf(m)%x - pbuf(m)%vx*tfrac
        pbuf(m)%y = pbuf(m)%y - pbuf(m)%vy*tfrac
        pbuf(m)%z = zssurf
        pbuf(m)%preside = OH_PCL_TO_BE_ACCUMULATED
    end if
#endif
#if ssurf==1
    if (pbuf(m)%z .le. zssurf) then
        tfrac = (pbuf(m)%z - zssurf)/pbuf(m)%vz
        pbuf(m)%x = pbuf(m)%x - pbuf(m)%vx*tfrac
        pbuf(m)%y = pbuf(m)%y - pbuf(m)%vy*tfrac
        pbuf(m)%z = zssurf
        pbuf(m)%preside = OH_PCL_TO_BE_ACCUMULATED
    end if
#elif ssurf==2
    if (pbuf(m)%z .le. zssurf) then
        xsepa = pbuf(m)%x - xbowlc
        ysepa = pbuf(m)%y - ybowlc
        zsepa = pbuf(m)%z - zbowlc
        if (xsepa*xsepa + ysepa*ysepa + zsepa*zsepa .ge. rbwlsq) then
            pbuf(m)%preside = OH_PCL_TO_BE_ACCUMULATED
        end if
    end if
    xsepa = pbuf(m)%x - xdomec
    ysepa = pbuf(m)%y - ydomec
    zsepa = pbuf(m)%z - zdomec
    if (xsepa*xsepa + ysepa*ysepa + zsepa*zsepa .le. rdomsq) then
        pbuf(m)%preside = OH_PCL_TO_BE_ACCUMULATED
    end if
#elif ssurf==3
    if (boundary_type == "rectangle-hole" &
        .or. boundary_type == "cylinder-hole" &
        .or. boundary_type == "hyperboloid-hole") then
        if (pbuf(m)%z <= zssurf+zsbuf) then
            pbuf(m)%preside = OH_PCL_TO_BE_ACCUMULATED
        end if
    end if
    if (pbuf(m)%z <= zssurf+zsbuf) then
        if (.not. &
            ((xlrechole(1) <= pbuf(m)%x .and. pbuf(m)%x <= xurechole(1) .and. &
              ylrechole(1) <= pbuf(m)%y .and. pbuf(m)%y <= yurechole(1) .and. &
              zlrechole(1) <= pbuf(m)%z .and. pbuf(m)%z <= zurechole(1)) .or. &
             (xlrechole(2) <= pbuf(m)%x .and. pbuf(m)%x <= xurechole(2) .and. &
              ylrechole(2) <= pbuf(m)%y .and. pbuf(m)%y <= yurechole(2) .and. &
              zlrechole(2) <= pbuf(m)%z .and. pbuf(m)%z <= zurechole(2)))) then
            pbuf(m)%preside = OH_PCL_TO_BE_ACCUMULATED
        end if
    end if
#endif
