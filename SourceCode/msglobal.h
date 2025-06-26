c	  global parameters:
c	  nzmax: max. interface index;
c	  lmax: max. no of total homogeneous layers (lmax <= nzmax-2);
c	  nrmax: max. no of traces;
c	  nfmax: max. no of frequency samples.
c
      integer nzmax,lmax,nrmax,nfmax
      parameter(nzmax=502,lmax=500,nrmax=200,nfmax=2048)
c
c	  atmospheric parameters
c
      double precision roair,vpair,km2m
      parameter(roair=0.1300d+01,vpair=0.3318d+03,km2m=1.0d+03)
c
      integer lp,nno(nzmax)
      double precision hp(nzmax)
      common /sublayer/ hp,lp,nno
c
c     zrec: receiver depth
c     lzrec: sublayer no of receiver
c
      integer lzrec
      double precision zrec
      common /receiver/ zrec,lzrec
c
c	  filter for free surface reflection
c
      integer surfilter
      common /surface/ surfilter
c       
c     model parameter:
c     n0: number of homogeneous layers
c
      integer n0
      double precision h(lmax),ro(lmax),vp(lmax),vs(lmax)
      double precision qp(lmax),qs(lmax)
      double complex cmp(lmax),cms(lmax)
      common /model/ h,ro,vp,vs,qp,qs,cmp,cms,n0
c
c	  suppression of time-domain aliasing effect
c	  with the complex frequency technique
c
      integer ialias
      double precision fimag
      common /alias/ fimag,ialias
c
c	  original model parameters
c
      integer l0
      double precision z1(lmax),z2(lmax),ro1(lmax),ro2(lmax)
      double precision vp1(lmax),vp2(lmax),vs1(lmax),vs2(lmax)
      double precision qp1(lmax),qp2(lmax),qs1(lmax),qs2(lmax)
      common /model0/z1,z2,ro1,ro2,vp1,vp2,vs1,vs2,qp1,qp2,qs1,qs2,l0
c
      double complex accair,kpair
      double complex acc(lmax),kp(lmax),ks(lmax),cla(lmax),cmu(lmax)
      common /cpara/ accair,kpair,acc,kp,ks,cla,cmu
c
c     source parameters
c
      integer ls,ms,ics
      integer kpower(6)
      double precision zs
      double precision sfct(6)
      common /source/ zs,sfct,ls,ms,ics,kpower

