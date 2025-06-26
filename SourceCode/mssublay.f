	subroutine mssublay(resolut,f,ierr)
	implicit none
c
	integer ierr
	double precision f,resolut(3)
c
	include 'msglobal.h'
c
c	work space
c
	double precision eps
	parameter(eps=1.0d-01)
	integer i,i0,i1,i2,l
	double precision dh,dro,dvp,dvs,dqp,dqs,z,dz
c
	n0=0
c
	do l=1,l0-1
	  dz=z2(l)-z1(l)
	  dvp=2.d0*dabs(vp2(l)-vp1(l))/(vp2(l)+vp1(l))
	  if(vs1(1).gt.0.d0)then
	    dvs=2.d0*dabs(vs2(l)-vs1(l))/(vs2(l)+vs1(l))
	    i2=idnint(2.d0*f*dz*dmax1(dvp/(vp2(l)+vp1(l)),
     &         dvs/(vs2(l)+vs1(l)))/eps)
	  else
	    i2=idnint(2.d0*f*dz*dvp/(vp2(l)+vp1(l))/eps)
	    dvs=0.d0
	  endif
	  dro=2.d0*dabs(ro2(l)-ro1(l))/(ro2(l)+ro1(l))
	  i1=idnint(dmax1(dro/resolut(1),
     &       dvp/resolut(2),dvs/resolut(3)))
	  i0=max0(1,min0(i1,i2))
	  dro=(ro2(l)-ro1(l))/dz
	  dvp=(vp2(l)-vp1(l))/dz
	  dvs=(vs2(l)-vs1(l))/dz
	  dqp=(qp2(l)-qp1(l))/dz
	  dqs=(qs2(l)-qs1(l))/dz
	  dh=dz/dble(i0)
	  do i=1,i0
	    n0=n0+1
	    if(n0.ge.lmax)then
	      ierr=1
	      return
	    endif
	    h(n0)=dh
	    z=(dble(i)-0.5d0)*dh
	    ro(n0)=ro1(l)+dro*z
	    vp(n0)=vp1(l)+dvp*z
	    vs(n0)=vs1(l)+dvs*z
	    qp(n0)=qp1(l)+dqp*z
	    qs(n0)=qs1(l)+dqs*z
	    if(qp(n0).gt.0.d0)then
	      cmp(n0)=dcmplx(1.d0,1.d0/qp(n0))
	    else
	      cmp(n0)=(1.d0,0.d0)
	    endif
	    if(qs(n0).gt.0.d0)then
	      cms(n0)=dcmplx(1.d0,1.d0/qs(n0))
	    else
	      cms(n0)=(1.d0,0.d0)
	    endif
	  enddo
	enddo
c
c	last layer is halfspace
c
	n0=n0+1
	h(n0)=0.d0
	ro(n0)=ro1(l0)
	vp(n0)=vp1(l0)
	vs(n0)=vs1(l0)
	qp(n0)=qp1(l0)
	qs(n0)=qs1(l0)
	if(qp(n0).gt.0.d0)then
	  cmp(n0)=dcmplx(1.d0,1.d0/qp(n0))
	else
	  cmp(n0)=(1.d0,0.d0)
	endif
	if(qs(n0).gt.0.d0)then
	  cms(n0)=dcmplx(1.d0,1.d0/qs(n0))
	else
	  cms(n0)=(1.d0,0.d0)
	endif
c
	write(*,'(7a)')' no ',' thick(km) ','  vp(km/s) ',
     &    '  vs(km/s) ',' ro(g/cm^3)','    qp   ','    qs'
	do i=1,n0
	  dh=h(i)*1.d-05
	  dvp=vp(i)*1.d-05
	  dvs=vs(i)*1.d-05
	  write(*,1001)i,dh,dvp,dvs,ro(i),qp(i),qs(i)
	enddo
1001	format(i4,4f11.4,2f8.1)
	ierr=0
	return
	end
