      program metropolis

      integer i,j,l,k,q,p
      parameter (N=20000,k=1000,p=100)
      real x(N),y(24*N)


      open(1,file='data_L30.dat',status='old')
      open(2,file='prova_term.dat',status='unknown')

      do i=1,24*N
        read (unit=1, fmt=*) y(i)
      enddo
      do i=1,N
        read (unit=1, fmt=*) x(i)
      enddo
      
      do i=1,N
        write (unit=2, fmt=*) i,x(i)
      end do
      close(1)
      close(2)

      end program
