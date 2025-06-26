      subroutine mshkpsv(hk,f,k,z,n)
      implicit none
c
	include 'msglobal.h'
c
      integer n
      double precision f,k,z
      double complex hk(4,4)
c
      double complex c1,c2,c4
	parameter(c1=(1.d0,0.d0),c2=(2.d0,0.d0),c4=(4.d0,0.d0))
c
	double complex ck,ck2,cx,cx2,cz,cep,cem,cch,csh
      double complex eta,sigma,gamma,kp2,ks2,cz2
c
      ck=dcmplx(k,0.d0)
      ck2=dcmplx(k*k,0.d0)
      cz=dcmplx(z,0.d0)
	cz2=dcmplx(0.5d0*z*z,0.d0)
	cx=dcmplx(k*z,0.d0)
	cx2=dcmplx(k*k*z*z,0.d0)
	kp2=kp(n)*kp(n)
	ks2=ks(n)*ks(n)
	eta=cla(n)/(cla(n)+c2*cmu(n))
	sigma=cmu(n)/(cla(n)+cmu(n))
	gamma=(cla(n)+cmu(n))/(cla(n)+c2*cmu(n))
c
c	propagatior matrix for P-SV waves with |kp*z| and |ks*z| << 1.
c
	hk(1,1)=c1+(kp2-c2*ck2*gamma)*cz2
	hk(1,2)=cz/(cla(n)+c2*cmu(n))
	hk(1,3)=eta*cx
	hk(1,4)=ck*gamma/cmu(n)*cz2
	hk(2,1)=-acc(n)*cz
	hk(2,2)=c1+(kp2-c2*ck2*gamma)*cz2
	hk(2,3)=c2*cmu(n)*ck*(ck2+ks2)*gamma*cz2
	hk(2,4)=cx
	hk(3,1)=-cx
	hk(3,2)=-ck*gamma/cmu(n)*cz2
	hk(3,3)=c1+(ks2+c2*ck2*gamma)*cz2
	hk(3,4)=cz/cmu(n)
	hk(4,1)=-c2*cmu(n)*ck*(ck2+ks2)*gamma*cz2
	hk(4,2)=-eta*cx
	hk(4,3)=(c4*ck2*cmu(n)*gamma-acc(n))*cz
	hk(4,4)=c1+(ks2+c2*ck2*gamma)*cz2
c
	return
	end
