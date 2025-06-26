      subroutine msmatrix(a,m,k,z,n)
      implicit none
c
      integer m,n
      double precision k,z
      double complex a(m,m)
c
	include 'msglobal.h'
c
      double complex c0,c1,c2
	parameter(c0=(0.d0,0.d0),c1=(1.d0,0.d0),c2=(2.d0,0.d0))
c
      double complex ck,ck2,cx,cz
      double complex gamma,alfa
c
      if(m.ne.4)then
        print *,'error in matrix: n != 4'
        return
      endif
      ck=dcmplx(k,0.d0)
      ck2=dcmplx(k*k,0.d0)
c
      a(1,1)=kp(n)
      a(1,2)=-kp(n)
      a(1,3)=ck
      a(1,4)=ck
      a(2,1)=c2*cmu(n)*ck2-acc(n)
      a(2,2)=a(2,1)
      a(2,3)=c2*ck*cmu(n)*ks(n)
      a(2,4)=-a(2,3)
      a(3,1)=ck
      a(3,2)=ck
      a(3,3)=ks(n)
      a(3,4)=-ks(n)
      a(4,1)=c2*ck*cmu(n)*kp(n)
      a(4,2)=-a(4,1)
      a(4,3)=c2*cmu(n)*ck2-acc(n)
      a(4,4)=a(4,3)
c
      return
      end
