implicit none
integer L
real(8) U,v,mu
integer Nv,k1,k2,k3,k4,i
integer Ly_model
real(8) U_model,vfpi_model

print*,'Ly,U,vfpi=?'
read*,Ly_model,U_model,vfpi_model

v=vfpi_model/Ly_model
U=U_model/Ly_model

i=Ly_model/4
select case (mod(Ly_model,4))
  case (0)
    L=2*i+1
  case (1)
    L=2*i
  case (2)
    L=2*i+1
  case (3)
    L=2*i+2
  case default
    stop 'error occurs'
end select

mu=U/2*L

Nv=0
do k1=1,L; do k2=1,L; do k3=1,L; do k4=1,L
  if(k1+k2==k3+k4)then
    Nv=Nv+1
  end if
end do; end do; end do; end do

open(10,file='ed.in')
write(10,'(3i8)') L*2,L*2,Nv
do i=1,L
  write(10,'(2i8,2f8.3)') i,i,v*(2*i-L-1)-mu,0d0
end do
do i=1,L
  write(10,'(2i8,2f8.3)') i+L,i+L,-v*(2*i-L-1)-mu,0d0
end do

do k1=1,L; do k2=1,L; do k3=1,L; do k4=1,L
  if(k1+k2==k3+k4)then
    write(10,'(4i8,2f8.3)') k1,k2+L,k3+L,k4,U,0d0
  end if
end do; end do; end do; end do

close(10)

end

