  module lemses
!
    use hdf5io
    use wrapper
    implicit none
!
  contains

  subroutine arange(aa,bb,nam,id)
!
    implicit none
!args
    integer(kind=4) :: id
    integer(kind=8) :: nam
    real(kind=8) :: aa(nam),bb(nam)
!vars
    integer(kind=8) :: i,j,k,l, i0, js, ip, namc
    integer(kind=4) :: iia(1024)
!-------------------- 
      if(id.le.1) return
      if(id.ge.1024) then
        print *,'NUMBER OF MAXWELLIAN EXCEEDS THE LIMIT 1024'
        id = 1024
      end if
      call qsort(aa,nam)
      ip = nam/id
      namc = ip*id
      if(namc.lt.nam) then
        print *,'ERROR IN ARANGE, NUMBER OF ELEMENTS NOT CONSISTANT'
        return
      end if
      do i = 1,id
        iia(i) = ip*(i - 1) + 1
      end do
      call rexch(int(id,8),i1=iia)
      do i = 1,nam
        bb(i) = aa(i)
      end do
      js = 1
      do k = 1,id
        j = js
        i = iia(k)
        i0 = i
        do l = 1,ip
          aa(i) = bb(j)
          i = i + 1
          j = j + id
        end do
        call rexch(ip,d1=aa(i0))
        js = js + 1
      end do
!
    return
  end subroutine arange


  subroutine arymod(aa,bb,in,nam,is,id)
!
    implicit none
!args
    integer(kind=8) :: in, nam
    integer(kind=4) :: is, id
    real(kind=8) :: aa(in),bb(in)
!vars
    integer(kind=8) :: inc, i,j
!-------------------- 
      inc = nam*id
      if(in.lt.inc) then
        print*,'ERROR IN ARYMOD, NUMBER OF ELEMENTS INCONSISTANT'
        print*,'IN,INC,NAM,ID', in, inc, nam, id
        return
      end if
      do i=1,nam
        bb(i) = aa(i)
      end do
      i = is
      do j=1,nam
        aa(i) = bb(j)
        i = i + id
      end do
!
    return
  end subroutine arymod


  subroutine mtsatn(m,dt,nstep,w,attn,damp)
!
    implicit none
!args
    integer(kind=4) :: m, nstep, i
    real(kind=8) :: dt, w, attn, damp
!vars
    real(kind=8) :: theta
!-------------------- 
      attn = 0.0d0
      if(mod(m,2).eq.0) then
        do i=1,m/2
          theta = w*dt*(float(i) - 0.5d0)
          attn = attn + cos(theta)
        end do
        attn = attn*2.0d0/float(m)
      else
        attn = 1.0d0
        do i=1,(m-1)/2
          theta = w*dt*float(i)
          attn = attn + 2.0d0*cos(theta)
        end do
        attn = attn/float(m)
      end if
      damp = attn**(nstep/m)
      attn = log(attn)/(dt*m)
!
    return
  end subroutine mtsatn


  subroutine qsort(x,n)
!
    implicit none
!args
    integer(kind=8) :: n
    real(kind=8) :: x(n)
!vars
    integer(kind=4),parameter :: ipmax=30
    integer(kind=8) :: il(ipmax),ir(ipmax),ip,i,j,k
    real(kind=8) :: p,t
!-------------------- 
      ip     = 1
      il(ip) = 1
      ir(ip) = n
!
   10 continue
      if(ip.gt.0) then
        if(il(ip).ge.ir(ip)) then
          ip = ip - 1
        else
          i  = il(ip)
          j  = ir(ip)
          k  = (i + j)/2
          p = max(min(x(i),x(j)),min(x(j),x(k)),min(x(k),x(i)))
   20     continue
          if(i.le.j) then
   30       continue
            if(x(i).lt.p) then
              i    = i + 1
              goto 30
            end if
   40       continue
            if(x(j).gt.p) then
              j    = j - 1
              goto 40
            end if
            if(i.le.j) then
              t    = x(i)
              x(i) = x(j)
              x(j) = t
              i    = i + 1
              j    = j - 1
            end if
            goto 20
          end if
          if((j-il(ip)).lt.(ir(ip)-i)) then
            il(ip+1) = il(ip)
            ir(ip+1) = j
            il(ip  ) = i
          else
            il(ip+1) = i
            ir(ip+1) = ir(ip)
            ir(ip  ) = j
          end if
          ip = ip + 1
          if(ip.gt.ipmax) then
            print*,'STUCK LEVEL OVER FLOW!!! (QSORT)'
            stop
          end if
        end if
        goto 10
      end if
!
    return
  end subroutine qsort


  subroutine rexch(n,i1,i2,i3,i4,r1,r2,r3,r4, &
 &                   d1,d2,d3,d4,d5,d6,d7,d8)
!
    implicit none
!args
    integer(kind=8),intent(in) :: n
    integer(kind=4),optional,intent(inout) :: i1(n),i2(n),i3(n),i4(n)
    real(kind=4),optional,intent(inout) :: r1(n),r2(n),r3(n),r4(n)
    real(kind=8),optional,intent(inout) :: d1(n),d2(n),d3(n),d4(n)
    real(kind=8),optional,intent(inout) :: d5(n),d6(n),d7(n),d8(n)
!vars
    integer(kind=8) :: nexch, i, n1,n2
    integer(kind=4) :: j, icon
    integer(kind=4) :: ia
    real(kind=4) :: ra
    real(kind=8) :: da
    real(kind=8) :: dranu(n)
!-------------------- 
      nexch = int(n/2)
      do j=1,4
        call ranu0(dranu,n,icon)
        if(icon.ne.0) print*,'[REXCH] RANU0: icon=',icon
        do i=1,nexch
          n1 = int(n*dranu(i*2-1)) + 1
          n2 = int(n*dranu(i*2)) + 1
          if(n1.le.0) cycle
          if(n2.le.0) cycle
          if(n1.gt.n) cycle
          if(n2.gt.n) cycle
!
          if(present(i1)) then
            ia = i1(n1)
            i1(n1) = i1(n2)
            i1(n2) = ia
          end if
!
          if(present(i2)) then
            ia = i2(n1)
            i2(n1) = i2(n2)
            i2(n2) = ia
          end if
!
          if(present(i3)) then
            ia = i3(n1)
            i3(n1) = i3(n2)
            i3(n2) = ia
          end if
!
          if(present(i4)) then
            ia = i4(n1)
            i4(n1) = i4(n2)
            i4(n2) = ia
          end if
!
          if(present(r1)) then
            ra = r1(n1)
            r1(n1) = r1(n2)
            r1(n2) = ra
          end if
!
          if(present(r2)) then
            ra = r2(n1)
            r2(n1) = r2(n2)
            r2(n2) = ra
          end if
!
          if(present(r3)) then
            ra = r3(n1)
            r3(n1) = r3(n2)
            r3(n2) = ra
          end if
!
          if(present(r4)) then
            ra = r4(n1)
            r4(n1) = r4(n2)
            r4(n2) = ra
          end if
!
          if(present(d1)) then
            da = d1(n1)
            d1(n1) = d1(n2)
            d1(n2) = da
          end if
!
          if(present(d2)) then
            da = d2(n1)
            d2(n1) = d2(n2)
            d2(n2) = da
          end if
!
          if(present(d3)) then
            da = d3(n1)
            d3(n1) = d3(n2)
            d3(n2) = da
          end if
!
          if(present(d4)) then
            da = d4(n1)
            d4(n1) = d4(n2)
            d4(n2) = da
          end if
!
          if(present(d5)) then
            da = d5(n1)
            d5(n1) = d5(n2)
            d5(n2) = da
          end if
!
          if(present(d6)) then
            da = d6(n1)
            d6(n1) = d6(n2)
            d6(n2) = da
          end if
!
          if(present(d7)) then
            da = d7(n1)
            d7(n1) = d7(n2)
            d7(n2) = da
          end if
!
          if(present(d8)) then
            da = d8(n1)
            d8(n1) = d8(n2)
            d8(n2) = da
          end if
!
        end do
      end do

    return
  end subroutine rexch


  subroutine ringds(ar,na,vt,bbb)
!
    implicit none
!args
    integer(kind=8) :: na
    real(kind=8) :: ar(na), vt, bbb
!vars
    integer(kind=8) :: j
    real(kind=8) :: ppii,sq2,dv,ds,a,cc,v1,v2,avm,fv1,fv2,dvd
!-------------------- 
      ppii = 3.14159265358979d0
      sq2 = 1.41421356237309d0
      dv = vt
      ds = vt**2/2.0d0/na
      a = 1.0d0/vt**2
      cc = 1.0d0/(1.0d0 - bbb)
      v1 = vt
      v2 = v1
      j = 1
!
   60 continue
      avm = -a*v1*v1
      fv1 = cc*v1*(dexp(avm) - dexp(avm/bbb))
      dvd = ds/fv1
      v1 = v1 - dvd
      if(v1.le.0.0) go to 70
      ar(j) = v1 + 0.5*dvd
      j = j + 1
      go to 60
!
   70 continue
      avm = -a*v2*v2
      fv2 = cc*v2*(dexp(avm) - dexp(avm/bbb))
      dvd = ds/fv2
      v2 = v2 + dvd
      ar(j) = v2 - 0.5*dvd
      j = j + 1
      if(j.le.na) go to 70
      call rexch(na,d1=ar)
!
    return
  end subroutine ringds


  subroutine vflux(ar,na,v0,vt,isl,icon)
!
    implicit none
!args
    integer(kind=4) :: isl, icon
    integer(kind=8) :: na
    real(kind=8) :: ar(na), v0, vt
!vars
    integer(kind=8) :: i,j
    real(kind=8) :: pi = 3.14159265358979d0
    real(kind=8) :: aa,bb,ts,dd,vm,ds,v1,v2,cc,ee,dvd
!-------------------- 
      aa = dsqrt(2.0d0/pi)*vt*exp(-(v0/vt)**2*0.5d0)
      bb = erf(v0/(vt*dsqrt(2.0d0)))
      if(isl.eq. 1) ts = (aa + v0*(1.d0 + bb))*0.5d0
      if(isl.eq.-1) ts = (-aa + v0*(1.d0 - bb))*0.5d0
      dd = v0**2 + 4.0d0*vt**2
      if(v0.ge.0.d0) then
        if(isl.eq. 1) vm = (v0 + dsqrt(dd))*0.5d0
        if(isl.eq.-1) vm = -2.0d0*vt**2/(v0 + dsqrt(dd))
      else
        if(isl.eq.-1) vm = (v0 - dsqrt(dd))*0.5d0
        if(isl.eq. 1) vm =  2.0d0*vt**2/(v0 - dsqrt(dd))
      end if
!
      ds = dabs(ts)/na
      v1 = vm
      v2 = vm
      i = 1
      cc = 1.0d0/(dsqrt(2.0d0*pi)*vt)
!
   10 continue
      ee = cc*v1*dexp(-((v1 - v0)/vt)**2*0.5d0)
      if(ee.eq.0.d0) go to 20
      dvd = ds/ee
      v1 = v1 - dvd
      if(isl.eq. 1.and.v1.le.0.0) go to 20
      if(isl.eq.-1.and.v1.ge.0.0) go to 20
      ar(i) = v1 + 0.5d0*dvd
      i = i + 1
      go to 10
!
   20 continue
      ee = cc*v2*dexp(-((v2 - v0)/vt)**2*0.5d0)
      if(ee.eq.0.d0) go to 30
      dvd = ds/ee
      v2 = v2 + dvd
      ar(i) = v2 - 0.5d0*dvd
      i = i + 1
      if(i.le.na) go to 20
!
   30 continue
      do j=i,na
      ar(j) = 0.0d0
      end do
!
   40 continue
      icon = 0
      if(i.lt.na) icon = i
!
      call rexch(na,d1=ar)
!
    return
  end subroutine vflux


  subroutine vinit(ar,na,v0,vt,xlos)
!
    implicit none
!args
    integer(kind=8) :: na
    real(kind=8) :: ar(na), v0, vt, xlos
!vars
    integer(kind=8) :: n,i,j
    real(kind=8) :: ppii,sq2,vd,nph,ds,x,dx,dvp,axlos,dv,a
    real(kind=8) :: xx,xxm,ts,sqvc,v1,v2,avm,fv1,dvd,fv2
!-------------------- 
      ppii = 3.14159265358979d0
      sq2 = 1.41421356237309d0
      if(xlos.ne.0.) go to  20
      vd = sq2*vt
      nph = na/2
      ds = dsqrt(ppii)/na
      n = 0
      x = 0
      do i=1,nph
        dx = ds/exp(-x**2)
        x = x + dx
        dvp = vd*(x - 0.5*dx)
        n = n + 1
        ar(n) = v0 + dvp
        n = n + 1
        ar(n) = v0 - dvp
      end do
!
   10 continue
      call rexch(na,d1=ar)
    return
!
   20 continue
      axlos = abs(xlos)
      if(xlos.gt.0.) then
        vt = v0*sq2/sqrt(axlos)
      else
        v0 = vt*sqrt(axlos)/sq2
      end if
      dv = vt
      a = 1.0/(dv*dv)
      xx = 0.5*(axlos + 1.0)
      xxm = -xx
      ts = 0.5*(a**xxm)*gamma(xx)
      ds = ts/na
      sqvc = axlos*0.5/a
      v1 = sqrt(sqvc)
      v2 = v1
      j = 1
!
   60 continue
      avm = -a*v1*v1
      fv1 = v1**axlos*exp(avm)
      dvd = ds/fv1
      v1 = v1 - dvd
      if(v1.le.0.0) go to 70
      ar(j) = v1 + 0.5*dvd
      j = j + 1
      go to 60
!
   70 continue
      avm = -a*v2*v2
      fv2 = v2**axlos*exp(avm)
      dvd = ds/fv2
      v2 = v2 + dvd
      ar(j) = v2 - 0.5*dvd
      j = j + 1
      if(j.le.na) go to 70
      call rexch(na,d1=ar)
!
    return
  end subroutine vinit

  end module lemses
