program ssh_j
! set up the input file 'ed.in' for the SSH-J model
! the model parameters are: L,t1,t2,J1,J2,boundary,twist
! accross the boundary, hop=hop*boundary*exp(1i*2*pi*twist)
!                        J = J *boundary
  implicit none
  
  real,parameter :: pi=3.1415926

  ! model parameters
  integer L
  real t1,t2,J1,J2,U
  real boundary,twist    
  
  ! internal variables
  integer i,ii
  real re,im,J

  print*,'L,t1,t2,J1,J2,U,boundary,twist=?'
  read*,L,t1,t2,J1,J2,U,boundary,twist

  open(10,file='ed.in')
  write(10,*) 2*L,4*L,7*L

! hopping terms
  do i=1,L
    ii=mod(i,L)+1
    if(mod(i,2)==0)then
      re=t1;im=0.0
    else
      re=t2;im=0.0
    end if
    if(i==L)then
      re=re*boundary*cos(2*pi*twist)
      im=re*boundary*sin(2*pi*twist)
    end if
    write(10,'(2i4,2f6.2)') i,     ii,   re, im
    write(10,'(2i4,2f6.2)') ii,    i,    re, -im
    write(10,'(2i4,2f6.2)') i+L,   ii+L, re, im
    write(10,'(2i4,2f6.2)') ii+L,  i+L,  re, -im
 end do

! interaction terms
  do i=1,L
    ii=mod(i,L)+1
    if(mod(i,2)==0)then
      J=J1
    else
      J=J2
    end if
    if(i==L)then
      J=J*boundary
    end if
    write(10,'(4i4,2f6.2)') i,    ii+L, ii,   i+L,  J/2,  0.0
    write(10,'(4i4,2f6.2)') i+L,  ii,   ii+L, i,    J/2,  0.0
    write(10,'(4i4,2f6.2)') i,    ii,   ii,   i,    J/4,  0.0
    write(10,'(4i4,2f6.2)') i+L,  ii+L, ii+L, i+L,  J/4,  0.0
    write(10,'(4i4,2f6.2)') i,    ii+L, ii+L, i,    -J/4, 0.0
    write(10,'(4i4,2f6.2)') i+L,  ii,   ii,   i+L,  -J/4, 0.0
    write(10,'(4i4,2f6.2)') i,    i+L,  i+L,  i,    U,    0.0
  end do

  close(10)

end program
