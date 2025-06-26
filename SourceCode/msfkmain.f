      program msfkmain
      implicit none
c
	include 'msglobal.h'
c
c     work space
c
      integer i,j,l,lf,nf,ip,np,ierr
	integer runtime,time
	integer ifselect(2),unit(2)
      double precision pi2,f,f1,f2,df,p,p1,p2,dp
	double precision k,rr,z,energy,qartif
	double precision a(2),resolut(3)
      double complex cy(6)
	logical hydro
      character*50 inputfile,fseries(2)
	character comments*180
c
c     read input file
c
      print *,'#####################################################'
      print *,'#                                                   #'
      print *,'#               Welcome to the program              #'
      print *,'#                                                   #'
      print *,'#                                                   #'
      print *,'#   M   M   SSSS  EEEEE  III   SSSS  FFFFF  K  K    #'
      print *,'#   MM MM  S      E       I   S      F      K K     #'
      print *,'#   M M M   SSS   EEEE    I    SSS   FFFF   KK      #'
      print *,'#   M   M      S  E       I       S  F      K K     #'
      print *,'#   M   M  SSSS   EEEEE  III  SSSS   F      K  K    #'
      print *,'#                                                   #'
      print *,'#                  (Version 2008)                   #'
      print *,'#                                                   #'
      print *,'#                                                   #'
      print *,'#                      by                           #'
      print *,'#                 Rongjiang Wang                    #'
      print *,'#              (wang@gfz-potsdam.de)                #'
      print *,'#                                                   #'
      print *,'#           GeoForschungsZentrum Potsdam            #'
      print *,'#             Last modified: Agu 2008               #'
      print *,'#####################################################'
      print *,'                          '
	write(*,'(a,$)')'the input data file is '
      read(*,'(a)')inputfile
c
      runtime=time()
      pi2=8.d0*datan(1.d0)
c
      open(10,file=inputfile,status='old')
c
c	source parameters
c	=================
c
	call getdata(10,comments)
      read(comments,*)zs
      zs=km2m*zs
c
c	receiver parameters
c	===================
c
	call getdata(10,comments)
      read(comments,*)zrec
      zrec=km2m*zrec
	call getdata(10,comments)
      read(comments,*)f1,f2,df
      nf=idnint((f2-f1)/df)+1
	call getdata(10,comments)
      read(comments,*)p1,p2,dp
      np=idnint((p2-p1)/dp)+1
	call getdata(10,comments)
      read(comments,*)qartif
	if(qartif.gt.0.d0)then
	  qartif=0.5d0/qartif
	else
	  qartif=0.d0
	endif
c
c	output files
c	============
c
	call getdata(10,comments)
      read(comments,*)(ifselect(i),i=1,2)
	call getdata(10,comments)
      read(comments,*)(fseries(i),i=1,2)
c
c	global model parameters
c	=======================
c
	call getdata(10,comments)
	read(comments,*)surfilter
	call getdata(10,comments)
	read(comments,*)(resolut(i),i=1,3)
	call getdata(10,comments)
	read(comments,*)l
	do i=1,3
	  resolut(i)=1.d-02*resolut(i)
	enddo
c
c	multilayered model parameters
c	=============================
c
	do i=1,l
	  call getdata(10,comments)
	  read(comments,*)j,h(i),vp(i),vs(i),ro(i),qp(i),qs(i)
	  h(i)=km2m*h(i)
	  vp(i)=km2m*vp(i)
	  vs(i)=km2m*vs(i)
	enddo
c
c	end of inputs
c	=============
c
	close(10)
c
c	determine upper und lower parameter values of each layer
c
      l0=1
      z1(l0)=0.d0
      do i=2,l
        if(h(i).gt.h(i-1))then
          z1(l0)=h(i-1)
          vp1(l0)=vp(i-1)
          vs1(l0)=vs(i-1)
          ro1(l0)=ro(i-1)
          qp1(l0)=qp(i-1)
          qs1(l0)=qs(i-1)
c
          z2(l0)=h(i)
          vp2(l0)=vp(i)
          vs2(l0)=vs(i)
          ro2(l0)=ro(i)
          qp2(l0)=qp(i)
          qs2(l0)=qs(i)
          l0=l0+1
        else
          z1(l0)=h(i)
          vp1(l0)=vp(i)
          vs1(l0)=vs(i)
          ro1(l0)=ro(i)
          qp1(l0)=qp(i)
          qs1(l0)=qs(i)
        endif
      enddo
      z1(l0)=h(l)
      vp1(l0)=vp(l)
      vs1(l0)=vs(l)
      ro1(l0)=ro(l)
      qp1(l0)=qp(l)
      qs1(l0)=qs(l)
c
c     construction of sublayers at the cutoff frequency
c
	write(*,*)'the multi-layered model at the cuttoff frequency:'
	call mssublay(resolut,f2,ierr)
	if(ierr.eq.1)then
	  print *,'the layer index (lmax) too small defined!'
	  goto 500
	endif
      call mslayer(ierr)
c
      ms=0
      ics=1
      sfct(1)=-1.0d0/(pi2*ro(nno(ls))*vp(nno(ls))**2)
      kpower(1)=0
      sfct(2)=0.d0
      kpower(2)=0
	do i=3,6
	  sfct(i)=0.d+00
	  kpower(i)=0
	enddo
	if(nno(lzrec).eq.1)then
	  hydro=.true.
	else
	  hydro=.false.
	endif
c
	write(*,*)'  no  layer-no       z'
      z=0.d0
      do j=1,lp
        if(j.eq.lzrec.and.j.eq.ls)then
          write(*,1001)j,nno(j),z,' Here is source and receiver!'
        else if(j.eq.lzrec)then
          write(*,1001)j,nno(j),z,' Here is receiver!'
        else if(j.eq.ls)then
          write(*,1001)j,nno(j),z,' Here is source!'
        else
          write(*,1001)j,nno(j),z
        endif
        z=z+hp(j)
	enddo
c
c     open output file
c
	do i=1,2
        if(ifselect(i).eq.1)then
	    unit(i)=20+i
	    open(unit(i),file=fseries(i),status='unknown')
	    write(unit(i),'(a,$)')'  [s/km]-[Hz]'
	    do lf=0,nf-2
	      f=f1+dble(lf)*df
	      write(unit(i),'(f12.4,$)')f
	    enddo
	    write(unit(i),'(f12.4)')f2
	  endif
      enddo
c
	do ip=0,np-1
	  p=p1+dble(ip)*dp
        do i=1,2
          if(ifselect(i).eq.1)then
            write(unit(i),'(E13.4,$)')p
          endif
        enddo
        call msqmodel(f2,f,fimag)
        do lf=0,nf-1
          f=f1+dble(lf)*df
	    fimag=-qartif*f
c
c	    p in s/km, k in 1/m
c
	    k=pi2*f*p/km2m
	    if(f.eq.0.d0.or.p.eq.0.d0)then
	      do i=1,2
	        a(i)=0.d0
	      enddo
	    else
	      call mskern(cy,f,k,hydro,1.d-16)
	      do i=1,2
	        a(i)=pi2*f*cdabs(cy(2*i-1))
	      enddo
	    endif
	    if(lf.lt.nf-1)then
	      do i=1,2
	        if(ifselect(i).eq.1)then
	          write(unit(i),'(E12.4,$)')a(i)
	        endif
            enddo
	    else
            do i=1,2
	        if(ifselect(i).eq.1)then
	          write(unit(i),'(E12.4)')a(i)
	        endif
	      enddo
	    endif
	  enddo
	enddo
c
	do i=1,2
	  if(ifselect(i).eq.1)close(unit(i))
	enddo
c
      runtime=time()-runtime
      call timeprt(runtime)
c
1001  format(2i7,E12.4,a)
1002	format(a,f12.4,a,$)
 500  stop
      end
