  module ssl2ex
!
    implicit real*8 (a-h,o-z)
!
  contains

  subroutine dalu(a,k,n,epsz,ip,is,vw,icon)
!
!             DALU                LEVEL=1        DATE=86.05.28
!*    *** ALU    ***   *   *********************************************
!     *                                                                *
!     *   SSL2,SSL2/VP !OPYRIGHT FUJITSU LIMITED 1979,1982,1983        *
!     *                                                                *
!     *   A22-11-0202    ALU,DALU,QALU               VERSION-1         *
!     *                                                                *
!     *   AUTHOR....J.MIKAMI,H.MURAOKA,                                *
!     *             K.ITO               1976                           *
!     *                                                                *
!     *   'ALU' DE!OMPOSES THE NON SINGULAR REAL MATRIX BY !ROUT'S     *
!     *   METHOD(LU-DE!OMPOSITION).                                    *
!     *                                                                *
!     *                                                                *
!     *   USAGE                                                        *
!     *        !ALL ALU(A,K,N,EPSZ,IP,IS,VW,I!ON)                      *
!     *          A   ....GIVEN NON SINGULAR REAL MATRIX.               *
!     *              ....RESULTANT LU-DE!OMPOSED MATRIX L AND U.       *
!     *                    ARRAY A IN!LUDES LOWER TRIANGULAR MATRIX L  *
!     *                    ON LOWER TRIANGULAR AREA AND  UNIT UPPER    *
!     *                    TRIANGULAR MATRIX U ON UPPER TRIANGULAR     *
!     *                    AREA WITH NO DIAGONAL ELEMENTS OF U.        *
!     *                    2 DIMENSIONAL ARRAY AS A(K,N).              *
!     *          K   ....GIVEN ADJUSTABLE DIMENSION FOR ARRAY A.       *
!     *          N   ....GIVEN ORDER OF MATRIX A.                      *
!     *          EPSZ....GIVEN MAXIMUM RELATIVE ERROR.                 *
!     *                    IF ZERO ASSIGNED,BE ASSUMED AS IF DEFAULT   *
!     *                    VALUE WAS ASSIGNED(UNIT ROUND OFF*16).      *
!     *          IP  ....RESULTANT TRANSPOSITION VECTOR WHICH          *
!     *                    REPRESENTS ROW-EXCHANGING BY PARTIAL        *
!     *                    PIVOTING.                                   *
!     *                    1 DIMENSIONAL ARRAY,SIZE IS N.              *
!     *          IS  ....RESULTANT SIGN OF DETERMINANT A.              *
!     *          VW  ....AUXILIARY 1 DIMENSIONAL ARRAY,SIZE IS N.      *
!     *          ICON....RESULTANT CONDITION CODE.                     *
!     *                                                                *
!     *   METHOD                                                       *
!     *        CROUT'S METHOD                                          *
!     *                                                                *
!     *   SLAVE SUBROUTINE                                             *
!     *        AMACH                                                   *
!     *                                                                *
!     *   REMARK                                                       *
!     *        IF ASSIGNED 10**(-S) TO EPSZ, IT MEANS THAT IF          *
!     *        CANCELLATION MORE THAN S OCCURED,SO THAT S DIGITS WERE  *
!     *        LOSSED,THEN THE VALUE OF PIVOT IS ASSUMED TO BE ZERO.   *
!     *                                                                *
!     *                                                                *
!     *                                                                *
!     ******************************************************************
!
    implicit REAL*8 (A-H,O-Z)
    dimension a(k,n),vw(n),ip(n)
!    dimension mcode(6)
    real*8    sum
!      data mcode/2ha2,2h2-,2h11,2h-0,2h20,2h2 /
!
!    ------------------------------------------------------------------
!     entry point.
!     set condition code to default value.
!    ------------------------------------------------------------------
      icon=0
      ip(1)=1
!    ------------------------------------------------------------------
!     check arguments
!    ------------------------------------------------------------------
      if(k.lt.n.or.epsz.lt.0.0d0) go to 10
      if(n.gt.1) go to 1000
      if(n.eq.1) go to 15
!    ------------------------------------------------------------------
!     error  icon=30000
!    ------------------------------------------------------------------
   10 icon=30000
      go to 8000
   15 if(a(1,1)) 8000,20,8000
!    ------------------------------------------------------------------
!     error  icon=2000
!    ------------------------------------------------------------------
   20 icon=20000
      go to 8000
!    ------------------------------------------------------------------
!     determine scaling factor for each row.
!    ------------------------------------------------------------------
 1000 do 25 i=1,n
      vw(i)=dabs(a(i,1))
   25 continue
      do 30 j=2,n
      do 27 i=1,n
      if(vw(i).lt.dabs(a(i,j))) vw(i)=dabs(a(i,j))
   27 continue
   30 continue
!    ------------------------------------------------------------------
!     check if scaling factor is zero.
!    ------------------------------------------------------------------
      ogmax=0.0d0
      do 40 i=1,n
      temp1=vw(i)
      if(temp1.eq.0.0d0) go to 20
      if(temp1.gt.ogmax) ogmax=temp1
      vw(i)=1.0d0/temp1
   40 continue
!    ------------------------------------------------------------------
!     set eps to epsz. if epsz is equal to zero,set eps to default
!     value.
!    ------------------------------------------------------------------
      eps=epsz
      if(epsz.eq.0.0d0) eps=dmach(epsz)*16.0d0
      is=1
!.   ------------------------------------------------------------------
!.    determine the pivot of 1st column.
!.   ------------------------------------------------------------------
      og=a(1,1)
      s=dabs(og*vw(1))
      do 80 i=2,n
      if(s.ge.dabs(a(i,1)*vw(i))) go to 80
      og=a(i,1)
      s=dabs(og*vw(i))
      ip(1)=i
   80 continue
!.   ------------------------------------------------------------------
!.    begin triangular decomposition  with partial pivoting.
!.   ------------------------------------------------------------------
 1100 do 180 j=2,n
!.   ------------------------------------------------------------------
!.    exchange rows for partial pivoting.
!.   ------------------------------------------------------------------
      loctm1=j-1
      if(ip(j-1).eq.(j-1)) go to 91
      is=-is
      m=ip(j-1)
      vw(m)=vw(j-1)
      do 90 i=1,loctm1
      l=j-i
      t=a(loctm1,l)
      a(loctm1,l)=a(m,l)
      a(m,l)=t
   90 continue
   91 do 100 i=1,loctm1
      if(ip(i).eq.i) go to 100
      m=ip(i)
      t=a(i,j)
      a(i,j)=a(m,j)
      a(m,j)=t
  100 continue
!.   ------------------------------------------------------------------
!.    check if the pivot is zero.
!.   ------------------------------------------------------------------
      t=a(j-1,j-1)
      if(t.eq.0.0d0.or.dabs(t).lt.dabs(og*eps).or. &
     &   (og.eq.0.0d0.and.dabs(t).le.(ogmax*eps))) go to 20
      vw(j-1)=1.0d0/t
!.   ------------------------------------------------------------------
!.    lu-decomposition by columnwise
!.    at first stage,calculate matrix u entrys(:u(i,j))
!.   ------------------------------------------------------------------
 1200 a(1,j)=a(1,j)*vw(1)
      ik=1
      if(j.eq.2) go to 1300
      do 130 i=2,loctm1
      sum=0.0d0
      m=i-1
      if(ik.eq.1) m=1
      i900=i-1
      do 120 ii=1,i900
      sum=a(i,m )*a(m ,j)+sum
  120 m=m+ik
      ik=-ik
      a(i,j)=(a(i,j)-sum)*vw(i)
  130 continue
!.   ------------------------------------------------------------------
!.    at stage after calculating j-th column u elements,
!.    calculate j-th column l elements(:l(i,j)).
!.   ------------------------------------------------------------------
 1300 s=-1.0d0
      do 170 i=j,n
      sum=0.0d0
      u=a(i,j)
      m=j-1
      if(ik.eq.1) m=1
      i901=j-1
      do 140 ii=1,i901
      sum=a(i,m )*a(m ,j)+sum
  140 m=m+ik
      a(i,j)=u-sum
      t=dabs(a(i,j)*vw(i))
      if(s.ge.t) go to 160
      s=t
      og=u
      ip(j)=i
  160 ik=-ik
  170 continue
  180 continue
!.   ------------------------------------------------------------------
!.    check if diagonal element a(n,n) is relative zero or
!.    absolute zero.
!.   ------------------------------------------------------------------
 1400 t=a(n,n)
      if(t.eq.0.0d0.or.dabs(t).lt.dabs(og*eps).or. &
     &   (og.eq.0.0d0.and.dabs(t).le.(ogmax*eps))) go to 20
 8000 if(icon.ne.0) print*, "ERROR@DALU: icon =", icon

    return
  end subroutine dalu


  subroutine dlax(a,k,n,b,epsz,isw,is,vw,ip,icon)
!
!             DLAX                LEVEL=1        DATE=86.05.28
!*    *** LAX    ***   *   *********************************************
!     *                                                                *
!     *   SSL2,SSL2/VP COPYRIGHT FUJITSU LIMITED 1979,1983             *
!     *                                                                *
!     *   A22-11-0101    LAX,DLAX,QLAX            VERSION-1            *
!     *                                                                *
!     *   AUTHOR....J.MIKAMI,H.MURAOKA,                                *
!     *             K.ITO               1976                           *
!     *                                                                *
!     *   'LAX' SOLVES THE LINEAR EQUATIONS BY CROUT'S METHOD.         *
!     *   AX=B,WHERE A IS REGULAR MATRIX WITH REAL COEFFICIENT.        *
!     *                                                                *
!     *   USAGE                                                        *
!     *        CALL LAX(A,K,N,B,EPSZ,ISW,IS,VW,IP,ICON)                *
!     *          A   ....GIVEN REGULAR COEFFICIENT MATRIX.             *
!     *              ....RESULTANT LU-DECOMPOSED ENTRYS.               *
!     *                    2 DIMENSIONAL ARRAY AS A(K,N),SIZE IS N.    *
!     *          K   ....GIVEN ADJUSTABLE DIMENSION FOR ARRAY A.       *
!     *          N   ....GIVEN ORDER OF MATRIX A.                      *
!     *              ....RESULTANT SOLUTION VECTOR.                    *
!     *          B   ....GIVEN CONSTANT VECTOR.                        *
!     *              ....RESULTANT SOLUTION VECTOR.                    *
!     *                    1 DIMENSIONAL ARRAY,SIZE IS N.              *
!     *          EPSZ....GIVEN VALUE WHICH IS REFERED WHEN  RELATIVE   *
!     *                    ZERO IS DETECTED,POSITIVE.                  *
!     *                    DEFAULT VALUE IS SETTED IF ZERO ASSIGNED    *
!     *          ISW ....GIVEN CONTROL INFORMATION.                    *
!     *                    IF 1,SOLVE EQUATIONS ENTIRLY.               *
!     *                    IF 2,SOLVE EQUATIONS WITH LAST              *
!     *                    LU-DECOMPOSED ENTRYS.                       *
!     *          IS  ....RESULTANT SIGN OF DETERMINANT A.              *
!     *          IP  ....AUXILIARY 1 DIMENSIONED ARRAY,SIZE IS N.      *
!     *                    TRANSPOSITION VECTOR WHICH REPRESENTS       *
!     *                    ROW-EXCHANGING BY PARTIAL PIVOTING.         *
!     *          ICON....RESULTANT CONDITION CODE.                     *
!     *                                                                *
!     *   METHOD                                                       *
!     *        CROUT'S METHOD                                          *
!     *                                                                *
!     *   SLAVE SUBROUTINE                                             *
!     *        ALU,LUX,AMACH                                           *
!     *                                                                *
!     *   REFERENCE CODE                                               *
!     *                                                                *
!     *   REMARK                                                       *
!     *        WHEN PARTIAL PIVOTING,DO ROW EQUILIBRATED PARTIAL       *
!     *        PIVOTING                                                *
!     *                                                                *
!     *                                                                *
!     *                                                                *
!     ******************************************************************
!
    implicit REAL*8 (A-H,O-Z)
    dimension a(k,n),b(n),vw(n),ip(n)
!    dimension mcode(6)
!    data mcode/2ha2,2h2-,2h11,2h-0,2h10,2h1 /
!
      if(isw.eq.1) go to 1000
      if(isw.eq.2) go to 1100
      icon=30000
      go to 8000
 1000 call dalu(a,k,n,epsz,ip,is,vw,icon)
      if(icon.ne.0) go to 8000
 1100 call dlux(b,a,k,n,1,ip,icon)
 8000 if(icon.ne.0) print*, "ERROR@DLAX: icon =", icon
!
    return
  end subroutine dlax


  subroutine dluiv(fa,k,n,ip,icon)
!
!             DLUIV               LEVEL=1        DATE=86.05.28
!*    *** LUIV   ***   *   *********************************************
!     *                                                                *
!     *   SSL2,SSL2/VP COPYRIGHT FUJITSU LIMITED 1979,1982,1983        *
!     *                                                                *
!     *   A22-11-0602    LUIV,DLUIV,QLUIV              VERSION-1       *
!     *                                                                *
!     *   AUTHOR ... J.MIKAMI,K.ITO     1976                           *
!     *                                                                *
!     *                                                                *
!     *   'LUIV' INVERT NON SINGULAR REAL MATRIX A, WHEN LU-DECOMPOSED *
!     *   MATRIX L AND U ARE GIVEN.                                    *
!     *                                                                *
!     *                                                                *
!     *   USAGE                                                        *
!     *        CALL LUIV(FA,K,N,IP,ICON)                               *
!     *          FA  ....GIVEN LU-DECOMPOSED MATRIX L AND U.           *
!     *                    ARRAY FA INCLUDES LOWER TRIANGULAR MATRIX L *
!     *                    ON LOWER TRIANGULAR AREA  AND UNIT UPPER    *
!     *                    TRIANGULAR MATRIX U ON UPPER TRIANGULAR     *
!     *                    AREA WITH NO DIAGONAL ELEMENTS OF U.        *
!     *              ....RESULTANT INVERSE MATRIX OF MATRIX A.         *
!     *                    2 DIMENSIONAL ARRAY AS FA(K,N).             *
!     *          K   ....GIVEN ADJUSTABLE DIMENSION FOR ARRAY FA.      *
!     *          N   ....GIVEN ORDER OF MATRIX A.                      *
!     *          IP  ....GIVEN TRANSPOSITION VECTOR WHICH REPRESENTS   *
!     *                    ROW-EXCHANGING BY PARTIAL PIVOTING.         *
!     *                    1 DIMENSIONAL ARRAY,SIZE IS N.              *
!     *          ICON....RESULTANT CONDITION CODE.                     *
!     *                                                                *
!     *   METHOD                                                       *
!     *        INVERT NON SINGULAR REAL MATRIX A, COMPUTING AS BELOW   *
!     *        INVERSION OF U * INVERSION OF L * PERMUTATION MATRIX P. *
!     *                                                                *
!     *   SLAVE SUBROUTINE                                             *
!     *        NONE.                                                   *
!     *                                                                *
!     *                                                                *
!     *                                                                *
!     ******************************************************************
!
    implicit REAL*8 (A-H,O-Z)
    dimension  fa(k,n),ip(n)
!    dimension  mcode(6)
    real*8     sum
!    data mcode/2ha2,2h2-,2h11,2h-0,2h60,2h2 /
!
!    ------------------------------------------------------------------
!     entry point.
!     set condition code to default value.
!    ------------------------------------------------------------------
      icon=0
!    ------------------------------------------------------------------
!     check arguments
!    ------------------------------------------------------------------
      if(n.le.0.or.k.lt.n) goto 9100
      do 10 i=1,n,1
      if(fa(i,i).eq.0.0d0) goto 9000
   10 continue
!    ------------------------------------------------------------------
!     compute the inverse matrix of l.
!    ------------------------------------------------------------------
 1000 nn=n-1
      if(nn) 70,70,20
   20 do  60  j=1,nn,1
      ik=1
      fa(j,j)=1.0d0/fa(j,j)
      jp1=j+1
      do  50   i=jp1,n
      sum=0.0d0
      m=i-1
      if(ik.eq.1) m=j
      i900=i-1
      do 40 ij=j,i900
      sum=fa(i,m )*fa(m ,j)+sum
   40 m=m+ik
      ik=-ik
      fa(i,j)=-sum/fa(i,i)
   50 continue
   60 continue
   70 fa(n,n)=1.0d0/fa(n,n)
!    ------------------------------------------------------------------
!     compute the inverse matrix of u.
!    ------------------------------------------------------------------
 1100 if(n.eq.1) go to 8000
      if(n.eq.2) go to 120
      do  110 jj=3,n,1
      j=n-jj+3
      j1=j-1
      fa(j1,j)=-fa(j1,j)
      ik=-1
      iz=j-2
      do  100 ii=1,iz
      i=iz-ii+1
      sum=0.0d0
      m=j-1
      if(ik.eq.1) m=i+1
      i902=j-1
      i901=i+1
      do 80 ij=i901,i902
      sum=fa(i,m )*fa(m ,j)+sum
   80 m=m+ik
      ik=-ik
      fa(i,j)=-fa(i,j)-sum
  100 continue
  110 continue
  120 fa(1,2)=-fa(1,2)
!    ------------------------------------------------------------------
!     compute matrix b , where b is multiplication of the inverse
!     matrix of l by the inverse matrix of u.
!    ------------------------------------------------------------------
 1200 do  150  j=1,nn,1
      ik=1
      do  130  i=j,nn
      sum=0.0d0
      m=n
      if(ik.eq.1) m=i+1
      i903=i+1
      do 125 ij=i903,n
      sum=fa(i,m )*fa(m ,j)+sum
  125 m=m+ik
      ik=-ik
      fa(i,j)=fa(i,j)+sum
  130 continue
      jp1=j+1
      ik=1
      do  140  i=1,j
      sum=0.0d0
      m=n
      if(ik.eq.1) m=jp1
      do 135 ij=jp1,n
      sum=fa(i,m )*fa(m ,jp1)+sum
  135 m=m+ik
      ik=-ik
      fa(i,jp1)=sum
  140 continue
  150 continue
!    ------------------------------------------------------------------
!     now,matrix b has been computed.
!     next,compute matrix which is multiplication of matrix b by
!     permutation matrix p.
!     to be exact,change each column drived by transposition vector.
!     at first,check if transposition vector is right.
!    ------------------------------------------------------------------
 1300 do  230  ii=1,n,1
      i=n-ii+1
      if(ip(i).lt.i.or.ip(i).gt.n) goto 9100
      if(ip(i).eq.i) goto 230
      m=ip(i)
      do  220  l=1,n,1
      t=fa(l,i)
      fa(l,i)=fa(l,m)
      fa(l,m)=t
  220 continue
  230 continue
!    ------------------------------------------------------------------
!
!    ------------------------------------------------------------------
 8000 if(icon.ne.0) print*, "ERROR@DLUIV: icon =", icon
    return
!
 9000 icon=20000
      goto 8000
 9100 icon=30000
      goto 8000
!
  end subroutine dluiv


  subroutine dlux(b,a,k,n,isw,ip,icon)
!
!             DLUX                LEVEL=1        DATE=86.05.28
!*    *** LUX    ***   *   *********************************************
!     *                                                                *
!     *   SSL2,SSL2/VP COPYRIGHT FUJITSU LIMITED 1979,1982,1983        *
!     *                                                                *
!     *   A22-11-0302    LUX,DLUX,QLUX                VERSION-1        *
!     *                                                                *
!     *   AUTHOR....J.MIKAMI,H.MURAOKA,                                *
!     *             K.ITO               1976                           *
!     *                                                                *
!     *   'LUX' SOLVES THE LINEAR EQUATIONS LUX=PB.  WHERE L AND U ARE *
!     *   LOWER AND UNIT UPPER TRIANGULAR REAL MATRICES,X IS           *
!     *   SOLUTION VECTOR,P IS PERMUTATION MATRIX AND B IS REAL        *
!     *   CONSTANT VECTOR.                                             *
!     *                                                                *
!     *   USAGE                                                        *
!     *        CALL LUX(B,FA,K,N,ISW,IP,ICON)                          *
!     *          B   ....GIVEN CONSTANT VECTOR.                        *
!     *              ....RESULTANT SOLUTION VECTOR.                    *
!     *                    1 DIMENSIONAL ARRAY,SIZE IS N.              *
!     *          FA  ....GIVEN LU-DECOMPOSED(FACTORIZED) MATRIX WHICH  *
!     *                    CONSIST OF MATRIX L AND U.                  *
!     *                    2 DIMENSIONAL ARRAY AS FA(K,N).             *
!     *          K   ....GIVEN ADJUSTABLE DIMENSION FOR ARRAY FA.      *
!     *          N   ....GIVEN ORDER OF THE MARIX A.                   *
!     *          ISW ....GIVEN CONTROL INFORMATION.                    *
!     *                    IF 1,SOLVE LUX=PB.                          *
!     *                    IF 2,SOLVE LY=PB.                           *
!     *                    IF 3,SOLVE UZ=B.                            *
!     *          IP  ....GIVEN TRANSPOSITION VECTOR WHICH              *
!     *                    REPRESENTS ROW-EXCHANGING BY PARTIAL        *
!     *                    PIVOTING.                                   *
!     *                    1 DIMENSIONAL ARRAY,SIZE IS N.              *
!     *          ICON....RESULTANT CONDITION CODE.                     *
!     *                                                                *
!     *   METHOD                                                       *
!     *        BACK SUBSTITUTION,FORWARD SUBSTITUTION.                 *
!     *                                                                *
!     *   SLAVE SUBROUTINE                                             *
!     *        NONE.                                                   *
!     *                                                                *
!     *                                                                *
!     *                                                                *
!     ******************************************************************
!    ------------------------------------------------------------------
!     for further simplicity,we substitute array fa by array a
!     in below algorithm
!    ------------------------------------------------------------------
!
    implicit REAL*8 (A-H,O-Z)
    dimension  a(k,n),b(n),ip(n)
!    dimension  mcode(6)
    real*8     sum
!    data mcode/2ha2,2h2-,2h11,2h-0,2h30,2h2 /
!
!    ------------------------------------------------------------------
!     entry point
!     set condition code to default value.
!    ------------------------------------------------------------------
      icon=0
!    ------------------------------------------------------------------
!     check arguments
!    ------------------------------------------------------------------
      if(n.le.0.or.k.lt.n.or.isw.le.0.or.isw.ge.4) goto 9000
      do 10 i=1,n
      if(ip(i).lt.i.or.ip(i).gt.n) goto 9000
   10 continue
!.   ------------------------------------------------------------------
!.    if isw ~= 3 then compute ly=pb by forward substitution.
!.   ------------------------------------------------------------------
 1000 if(isw.eq.3) go to 1100
      if(a(1,1).ne.0.0d0) go to 15
   14 icon=20000
      go to 8000
   15 if(ip(1).eq.1) go to 16
      m=ip(1)
      t=b(1)
      b(1)=b(m)
      b(m)=t
   16 b(1)=b(1)/a(1,1)
      if(n.eq.1) go to 1100
      m=1
      ik=1
      do 40 i=2,n
      if(a(i,i).eq.0.0d0) icon=20000
      if(icon.ne.0) go to 40
      w=a(i,i)
      if(ip(i).eq.i) go to 20
      m=ip(i)
      t=b(i)
      b(i)=b(m)
      b(m)=t
   20 sum=0.0d0
      m=i-1
      if(ik.eq.1) m=1
      i900=i-1
      do 30 ij=1,i900
      sum=a(i,m )*b(m )+sum
   30 m=m+ik
      ik=-ik
      b(i)=(b(i)-sum)/w
   40 continue
!.   ------------------------------------------------------------------
!.    check if any diagonal elements are zero.
!.   ------------------------------------------------------------------
      if(icon.ne.0) go to 8000
!    ------------------------------------------------------------------
!     if isw~=2 then compute uz=b by back substitution.
!    ------------------------------------------------------------------
 1100 if(isw.eq.2.or.n.eq.1) go to 8000
      ik=-1
      do 60 j=2,n
      i=n-j+1
      sum=0.0d0
      m=n
      if(ik.eq.1) m=i+1
      i901=i+1
      do 50 ij=i901,n
      sum=a(i,m )*b(m )+sum
   50 m=m+ik
      ik=-ik
      b(i)=b(i)-sum
   60 continue
!    ------------------------------------------------------------------
!
!    ------------------------------------------------------------------
 8000 if(icon.ne.0) print*, "ERROR@DLUX: icon =", icon
    return
!
 9000 icon=30000
      goto 8000
!
  end subroutine dlux


  function dmach(eps)
!
!             DMACH               LEVEL=1        DATE=86.05.28
!*    *** DMACH  ***   *   *********************************************
!     *                                                                *
!     *   SSL2,SSL2/VP COPYRIGHT FUJITSU LIMITED 1979,1983             *
!     *                                                                *
!     *   SLAVE          DMACH                         VERSION-1       *
!     *                                                                *
!     *   ---------- BASE = 16 , DIGIT = 14 ----------                 *
!     *                                                                *
!     *   AUTHOR....J.MIKAMI            1976.7                         *
!     *                                                                *
!     *   'DMACH' SETS THE UNIT ROUND OFF FOR SINGLE PRECISION         *
!     *           ARITHMETIC.                                          *
!     *                                                                *
!     *   USAGE(FUNCTION SUBPROGRAM)                                   *
!     *        DMACH(EPS)                                              *
!     *             EPS...GIVEN, DUMMY ARGUMENT.                       *
!     *                                                                *
!     *   REMARK                                                       *
!     *        DMACH IS MACHIN DEPENDENT ROUTINE.                      *
!     *                                                                *
!     ******************************************************************
!
    implicit REAL*8 (A-H,O-Z)
    real*8        dmach,eps,o
!   integer       i(2)
!   equivalence (i(1) , o)  
!   data  i/$3cc00000,$00000000/
!
      dmach   = 2.220446049250313e-016
!      write(6,'(z16)') dmach
!
    return
  end function dmach


  subroutine dmav(a, k, m, n, x, y, icon)
!
    implicit REAL*8 (A-H,O-Z)
    integer  k, m, n, icon
    real*8   a(k, n), x(n), y(m)
    integer  i, j
    real*8   sum
!    dimension  mcode(6)
!
      icon = 0
      if (m.lt.1 .or. n.eq.0 .or. k.lt.m) then
        icon = 30000
        goto 999
      endif
!
      if (n.gt.0) then
        do 10 i = 1, m
          sum = 0.0d0
          do 20 j = 1, n
            sum = sum + a(i, j) * x(j)
   20     continue
          y(i) = sum
   10   continue
      else
        do 30 i = 1, m
          sum = 0.0d0
          do 40 j = 1, -n
            sum = sum + a(i, j) * x(j)
   40     continue
          y(i) = y(i) - sum
   30   continue
      endif
!
  999 if(icon.ne.0) print*, "ERROR@DMAV: icon =", icon
!
    return
  end subroutine dmav


  end module ssl2ex
