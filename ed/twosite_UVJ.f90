implicit none
integer N,a,b
real U,V,J,mu

print*,'N=? U=? V=? J=? mu=?'
read*,N,U,V,J,mu
mu=mu+V/2*(2*N-2)+U/2

open(10,file='ed.in')

write(10,*) 4*N, 8*N, 10*N**2-8*N

do a=1,N*2
  write(10,'(2i4,2f8.3)') a,a+N*2,-1d0,0d0
  write(10,'(2i4,2f8.3)') a+N*2,a,-1d0,0d0
end do

do a=1,N*4
  write(10,'(2i4,2f8.3)') a,a,-mu,0d0
end do

do a=1,N
  write(10,'(4i4,2f8.3)') a,a+N,a+N,a,U,0d0
  write(10,'(4i4,2f8.3)') a+2*N,a+N+2*N,a+N+2*N,a+2*N,U,0d0
  do b=a+1,N
    write(10,'(4i4,2f8.3)') a,b,b,a,V,0d0
    write(10,'(4i4,2f8.3)') a,b+N,b+N,a,V,0d0
    write(10,'(4i4,2f8.3)') a+N,b,b,a+N,V,0d0
    write(10,'(4i4,2f8.3)') a+N,b+N,b+N,a+N,V,0d0
    write(10,'(4i4,2f8.3)') a+2*N,b+2*N,b+2*N,a+2*N,V,0d0
    write(10,'(4i4,2f8.3)') a+2*N,b+N+2*N,b+N+2*N,a+2*N,V,0d0
    write(10,'(4i4,2f8.3)') a+N+2*N,b+2*N,b+2*N,a+N+2*N,V,0d0
    write(10,'(4i4,2f8.3)') a+N+2*N,b+N+2*N,b+N+2*N,a+N+2*N,V,0d0

    write(10,'(4i4,2f8.3)') a,b,b,a,-J,0d0
    write(10,'(4i4,2f8.3)') a,b+N,b,a+N,-J,0d0
    write(10,'(4i4,2f8.3)') a+N,b,b+N,a,-J,0d0
    write(10,'(4i4,2f8.3)') a+N,b+N,b+N,a+N,-J,0d0
    write(10,'(4i4,2f8.3)') a+2*N,b+2*N,b+2*N,a+2*N,-J,0d0
    write(10,'(4i4,2f8.3)') a+2*N,b+N+2*N,b+2*N,a+N+2*N,-J,0d0
    write(10,'(4i4,2f8.3)') a+N+2*N,b+2*N,b+N+2*N,a+2*N,-J,0d0
    write(10,'(4i4,2f8.3)') a+N+2*N,b+N+2*N,b+N+2*N,a+N+2*N,-J,0d0

    write(10,'(4i4,2f8.3)') a,a+N,b+N,b,J,0d0
    write(10,'(4i4,2f8.3)') b,b+N,a+N,a,J,0d0
    write(10,'(4i4,2f8.3)') a+2*N,a+N+2*N,b+N+2*N,b+2*N,J,0d0
    write(10,'(4i4,2f8.3)') b+2*N,b+N+2*N,a+N+2*N,a+2*N,J,0d0
  end do
end do

close(10)

end
