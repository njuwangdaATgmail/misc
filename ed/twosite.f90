implicit none
integer N,i,j
real U,mu

print*,'su(N), N=? U=? mu=?'
read*,N,U,mu

open(10,file='ed.in')
write(10,*) 2*N, 4*N, N*(N-1)
do i=1,N
  write(10,'(2i4,2f8.3)') i,i+N,-1d0,0d0
  write(10,'(2i4,2f8.3)') i+N,i,-1d0,0d0
end do
do i=1,N
  write(10,'(2i4,2f8.3)') i,i,U/2*(1-N)-mu,0d0
  write(10,'(2i4,2f8.3)') i+N,i+N,U/2*(1-N)-mu,0d0
end do
do i=1,N
  do j=i+1,N
    write(10,'(4i4,2f8.3)') i,j,j,i,U,0d0
    write(10,'(4i4,2f8.3)') i+N,j+N,j+N,i+N,U,0d0
  end do
end do
close(10)

end
