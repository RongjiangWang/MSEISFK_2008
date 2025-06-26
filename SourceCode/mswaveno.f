      subroutine mswaveno(f,k)
      implicit none
c
      double precision f,k
c
	include 'msglobal.h'
c
	double precision pi
	parameter(pi=3.14159265358979d0)
c
	integer n
      double complex cvpair,comega,ck,cro,cvp,cvs
c
	if(f.gt.0.d0)then
        comega=dcmplx(2.d0*pi*f,2.d0*pi*fimag)
	else
	  comega=(0.d0,0.d0)
	endif
      ck=dcmplx(k,0.d0)
	accair=dcmplx(roair,0.d0)*comega*comega
	cvpair=dcmplx(vpair,0.d0)
	kpair=cdsqrt((ck+comega/cvpair)*(ck-comega/cvpair))
	do n=1,n0
        cro=dcmplx(ro(n),0.d0)
        cvp=dcmplx(vp(n),0.d0)*cmp(n)
        kp(n)=cdsqrt((ck+comega/cvp)*(ck-comega/cvp))
        if(vs(n).gt.0.d0)then
          cvs=dcmplx(vs(n),0.d0)*cms(n)
          cmu(n)=cro*cvs*cvs
          ks(n)=cdsqrt((ck+comega/cvs)*(ck-comega/cvs))
	  else
          cmu(n)=(0.d0,0.d0)
          ks(n)=(0.d0,0.d0)
        endif
        cla(n)=cro*cvp*cvp-2.d0*cmu(n)
        acc(n)=dcmplx(ro(n),0.d0)*comega*comega
	enddo
c
      return
      end
