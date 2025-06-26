      subroutine msqmodel(f0,fr,fi)
      implicit none
c
c	calculate q based on the constant q model
c
c	f0 = reference frequency of given q
c	fr = frequency
c	fi = imaginary frequency (f_alias)
c
	double precision f0,fr,fi
c
	include 'msglobal.h'
c
	integer n
c
	do n=1,n0
	  if(qp(n).gt.0.d0)then
	    cmp(n)=dcmplx(1.d0,0.5d0/qp(n))
	  else
	    cmp(n)=(1.d0,0.d0)
	  endif
	  if(qs(n).gt.0.d0)then
	    cms(n)=dcmplx(1.d0,0.5d0/qs(n))
	  else
	    cms(n)=(1.d0,0.d0)
	  endif
	enddo
c
	return
	end
