      subroutine mskern(y,f,k,hydro,eps)
      implicit none
c
c	calculation of response to p-sv source
c       y(6): solution vector (complex)
c       k: wave number
c       f: frequency. 2*pi*f = cyclar frequency
c       eps: relative accuracy
c
	include 'msglobal.h'
c
      double precision f,k,eps
      double complex y(6)
	logical hydro
c
c     work space
c
      double precision eps0,epsp
      parameter(eps0=1.0d-03,epsp=1.0d-06)
c
      integer i,j,l,n,lup,llw,key,lground
      double precision pi,h0,fac,omi,exponent
      double precision wa(8,8),wb(8)
      double complex ca,cb,ck2,ch0,pwave,swave,cfac
      double complex y0(4,2),c0(4,2),c1(4,2),b(2)
      double complex y1(4,2),yup(4,2),ylw(4,2)
      double complex ma(4,4),mai(4,4),hk(4,4),hkw(2,2)
      double complex orth(2,2),coef(2,2),cnorm(2)
c
      do i=1,6
        y(i)=(0.d0,0.d0)
      enddo
c
      pi=4.d0*datan(1.d0)
	ck2=dcmplx(k*k,0.d0)
	if(f.eq.0.d0)return
c
	call mswaveno(f,k)
c
c===============================================================================
c
c     matrix propagation from surface to source
c
      do j=1,2
        do i=1,4
          c0(i,j)=(0.d0,0.d0)
          yup(i,j)=(0.d0,0.d0)
	    y0(i,j)=(0.d0,0.d0)
        enddo
      enddo
      omi=2.d0*pi*f
	lground=lp+1
      do l=ls+1,lp
        if(vs(nno(l-1)).le.0.d0.and.vs(nno(l)).gt.0.d0)then
          lground=l
	    goto 10
        endif
      enddo
10	continue
c
c     yup(2,1): the starting solution vectors
c
      lup=1
c
	if(surfilter.eq.1)then
	  yup(1,1)=kp(nno(lup))
	  yup(2,1)=-acc(nno(lup))
	else
	  yup(1,1)=kpair
	  yup(2,1)=-accair
c	  yup(1,1)=(1.d0,0.d0)
	endif
      if(lup.eq.lzrec)call cmemcpy(yup,y0,8)
c
      do l=lup+1,ls
	  n=nno(l-1)
        call mshkwa(hkw,f,k,hp(l-1),n)
        call caxcb(hkw,yup(1,1),2,2,1,c0(1,1))
        call cmemcpy(c0(1,1),yup(1,1),2)
c
	  if(l.gt.lzrec)then
	    cnorm(1)=cdexp(-kp(n)*dcmplx(hp(l-1),0.d0))
	    y0(1,1)=y0(1,1)*cnorm(1)
	    y0(2,1)=y0(2,1)*cnorm(1)
	  else if(l.eq.lzrec)then
	    call cmemcpy(yup(1,1),y0(1,1),2)
	  endif
      enddo
c
c===============================================================================
c
c     matrix propagation from half-space to sea ground
c
      do i=1,4
        do j=1,2
          c0(i,j)=(0.d0,0.d0)
          ylw(i,j)=(0.d0,0.d0)
        enddo
      enddo
c
c     determination of starting sublayer for half-space
c
      llw=lp
	if(llw.lt.lground)then
	  ylw(1,1)=kp(nno(llw))
	  ylw(2,1)=acc(nno(llw))
        if(llw.gt.ls.and.llw.eq.lzrec)call cmemcpy(ylw,y0,8)
	  goto 300
	endif
c
c
c     c0(4,2): 2 coefficient vectors in the half-space
c
	n=nno(llw)
      if(vs(n).le.0.d0)then
c
c       the lowest layer is fluid
c
	  ylw(1,1)=kp(n)
	  ylw(2,1)=acc(n)
        ylw(3,2)=(1.d0,0.d0)
      else
c
c       the lowest layer is solid
c
        c0(2,1)=(1.d0,0.d0)
        c0(4,2)=(1.d0,0.d0)
        call msmatrix(ma,4,k,0.d0,n)
        call caxcb(ma,c0,4,4,2,ylw)
      endif
      if(llw.gt.ls.and.llw.eq.lzrec)call cmemcpy(ylw,y0,8)
c
      do l=llw-1,lground,-1
        h0=hp(l)
	  ch0=dcmplx(h0,0.d0)
        n=nno(l)
	  if((cdabs(kp(n))+cdabs(ks(n)))*h0.le.eps0)then
	    call mshkpsv(hk,f,k,-h0,n)
	    call caxcb(hk,ylw,4,4,2,y1)
	    call cmemcpy(y1,ylw,8)
	  else
c
c         determination of propagation matrix
c
          call msmatinv(mai,4,k,0.d0,n)
          call msmatrix(ma,4,k,-h0,n)
          call caxcb(mai,ylw,4,4,2,c0)
          pwave=cdexp(-kp(n)*ch0)
          swave=cdexp(-ks(n)*ch0)
c
c         normalization of all modes
c
	    do j=1,2
            cnorm(j)=(0.d0,0.d0)
            do i=1,4
              cnorm(j)=cnorm(j)+c0(i,j)*dconjg(c0(i,j))
            enddo
	    enddo
	    cfac=(1.d0,0.d0)/cdsqrt(cnorm(1)*cnorm(2))
c
c         orthogonalization of the p-sv modes
c
          orth(1,1)=c0(4,2)*cfac
          orth(1,2)=-c0(2,2)*cfac
          orth(2,1)=-c0(4,1)*cfac
          orth(2,2)=c0(2,1)*cfac
          call caxcb(c0,orth,4,2,2,c1)
          if(l.lt.lzrec)then
c
c	      additional normalization to avoid overflow
c
	      do i=1,2
	        orth(i,1)=orth(i,1)*pwave
	        orth(i,2)=orth(i,2)*swave
	      enddo
            call caxcb(y0,orth,4,2,2,y1)
            call cmemcpy(y1,y0,8)
          endif
c
	    c1(1,1)=c1(1,1)*pwave*pwave
c	    c1(2,1)=c1(2,1)
	    c1(3,1)=c1(3,1)*pwave*swave
          c1(4,1)=(0.d0,0.d0)
c
          c1(1,2)=c1(1,2)*swave*pwave
          c1(2,2)=(0.d0,0.d0)
          c1(3,2)=c1(3,2)*swave*swave
c         c1(4,2)=c1(4,2)
c
          call caxcb(ma,c1,4,4,2,ylw)
	  endif
        if(l.gt.ls.and.l.eq.lzrec)call cmemcpy(ylw,y0,8)
      enddo
c
c===============================================================================
c
c     matrix propagation from sea ground to source
c
      ca=ylw(4,1)
      cb=ylw(4,2)
      do i=1,4
        ylw(i,1)=cb*ylw(i,1)-ca*ylw(i,2)
        ylw(i,2)=(0.d0,0.d0)
      enddo
      if(lground.gt.ls.and.lground.le.lzrec)then
        do i=1,4
          y0(i,1)=cb*y0(i,1)-ca*y0(i,2)
          y0(i,2)=(0.d0,0.d0)
        enddo
      endif
300	continue
c
      do l=min0(lground,llw)-1,ls,-1
	  n=nno(l)
        call mshkwa(hkw,f,k,-hp(l),n)
        call caxcb(hkw,ylw(1,1),2,2,1,c0(1,1))
	  call cmemcpy(c0(1,1),ylw(1,1),2)
c
	  if(l.lt.lzrec)then
	    cnorm(1)=cdexp(-kp(n)*dcmplx(hp(l),0.d0))
	    do i=1,4
	      y0(i,1)=y0(i,1)*cnorm(1)
	    enddo
	  else if(l.gt.ls.and.l.eq.lzrec)then
	    call cmemcpy(ylw(1,1),y0(1,1),2)
	  endif
      enddo
c
c===============================================================================
c
c     conditions on the source surface
c
c
c     source function
c
      do i=1,2
        if(kpower(i).eq.1)then
          fac=k
        else
          fac=1.d0
        endif
        b(i)=dcmplx(sfct(i)*fac,0.d0)
        coef(i,1)=yup(i,1)
        coef(i,2)=-ylw(i,1)
      enddo
      key=0
      call cgemp(wa,wb,coef,b,2,1,1.d-99,key)
      if(key.eq.0)then
        print *,'warning in mskern: anormal exit from cgemp!'
        return
      endif
      if(lzrec.le.ls)then
        do i=1,4
          y(i)=b(1)*y0(i,1)
        enddo
      else
        do i=1,4
          y(i)=b(2)*y0(i,1)
        enddo
      endif
      if(hydro)then
c
c       in case of hydrophone, y(1)=pressure=-stress(z)
c
	  y(1)=-y(2)
        y(2)=(0.d0,0.d0)
      endif
      return
      end
