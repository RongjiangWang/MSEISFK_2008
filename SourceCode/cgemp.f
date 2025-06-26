      subroutine cgemp(a,b,ca,cb,n,m,eps,key)
      implicit none
c-------------------------------------------------------------------------------
c     subroutine cgemp to solve linear complex equation system                 i
c       a: work space real matrix (2n,2n)                                      i
c       b: work space real matrix (2n,m)                                       i
c       ca: complex coefficient matrix(n,n);                                   i
c       cb: complex right-hand matrix(n,n1) by input,                          i
c           complex solution matrix(n,n1) by return;                           i
c       eps: control constant;                                                 i
c       key: the main term of a column                                         i
c            smaller than eps, key=0: anormal return,                          i
c            else key=1: normal return.                                        i
c-------------------------------------------------------------------------------
      integer n,m,key
      double precision eps
      double precision a(2*n,2*n),b(2*n,m)
      double complex ca(n,n),cb(n,m)
c
      integer i,j
c
      do i=1,n
        do j=1,n
          a(i,j)=dreal(ca(i,j))
          a(n+i,j)=dimag(ca(i,j))
          a(i,n+j)=-a(n+i,j)
          a(n+i,n+j)=a(i,j)
        enddo
        do j=1,m
          b(i,j)=dreal(cb(i,j))
          b(n+i,j)=dimag(cb(i,j))
        enddo
      enddo
      call gemp(a,b,2*n,m,eps,key)
      do i=1,n
        do j=1,m
          cb(i,j)=dcmplx(b(i,j),b(n+i,j))
        enddo
      enddo
c
      return
      end
