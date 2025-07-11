      subroutine mslayer(ierr)
      implicit none
c
      integer ierr
c
	include 'msglobal.h'
c
      integer l,n,li,lp0
      double precision zswap
      double precision z(nzmax),z0(nzmax)
c
      lp0=1
      z0(lp0)=0.d0
      do n=1,l0-1
        lp0=lp0+1
        z0(lp0)=z0(lp0-1)+h(n)
      enddo
      lp0=lp0+1
      z0(lp0)=zrec
      lp0=lp0+1
      z0(lp0)=zs
c
c     sort the z0-profile
c
      do l=1,lp0-1
        do li=l+1,lp0
          if(z0(li).lt.z0(l))then
            zswap=z0(l)
            z0(l)=z0(li)
            z0(li)=zswap
          endif
        enddo
      enddo
c
c     delete duplicates
c
      lp=1
      z(lp)=0.d0
      do l=2,lp0
        if(z0(l).gt.z(lp))then
          hp(lp)=z0(l)-z(lp)
          lp=lp+1
          z(lp)=z0(l)
        endif
      enddo
      hp(lp)=0.d0
c
c     determine ls,lzrec
c
      do l=1,lp
        if(z(l).eq.zs)ls=l
        if(z(l).eq.zrec)lzrec=l
      enddo
c
c     determine layer no of each depth
c
      li=1
      zswap=h(1)
      nno(1)=1
      do l=2,lp
        if(z(l).ge.zswap.and.li.lt.n0)then
          li=li+1
          zswap=zswap+h(li)
        endif
        nno(l)=li
      enddo
c
      return
      end
