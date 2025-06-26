      subroutine msmatinv(a,m,k,z,n)
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
      integer i,j
      double complex ck,ck2,cx,cz
      double complex gamma,alfa,psi
c
      if(m.ne.4)then
        print *,'error in matinv: n != 4'
        return
      endif
      ck=dcmplx(k,0.d0)
      ck2=dcmplx(k*k,0.d0)
c
      a(1,1)=(acc(n)-c2*ck2*cmu(n))/(c2*acc(n)*kp(n))
      a(1,2)=-c1/(c2*acc(n))
      a(1,3)=ck*cmu(n)/acc(n)
      a(1,4)=ck/(c2*acc(n)*kp(n))
      a(2,1)=-a(1,1)
      a(2,2)=a(1,2)
      a(2,3)=a(1,3)
      a(2,4)=-a(1,4)
      a(3,1)=ck*cmu(n)/acc(n)
      a(3,2)=ck/(c2*acc(n)*ks(n))
      a(3,3)=(acc(n)-c2*ck2*cmu(n))/(c2*acc(n)*ks(n))
      a(3,4)=a(1,2)
      a(4,1)=a(3,1)
      a(4,2)=-a(3,2)
      a(4,3)=-a(3,3)
      a(4,4)=a(1,2)
c
      return
      end
