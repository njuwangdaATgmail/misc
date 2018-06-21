IMPLICIT NONE
INTEGER nsec,n,sec,i,dmax,flag,j,nflv,flv,sec1
REAL(8) re,im
INTEGER, ALLOCATABLE :: d(:),table(:,:)
COMPLEX(8), ALLOCATABLE :: h(:,:)
INTEGER, ALLOCATABLE :: c(:,:,:)

OPEN(10,file='ham.dat')
READ(10,*) nsec
ALLOCATE(d(nsec))
DO sec=1,nsec
  READ(10,*) i,d(sec)
END DO
dmax=maxval(d)
ALLOCATE(table(nsec,dmax))
n=0
DO sec=1,nsec
  DO i=1,d(sec)
    n=n+1
    table(sec,i)=n
  END DO
  !PRINT*,sec,':',table(sec,1:d(sec))
END DO
print*,'nhilbert=',n
ALLOCATE(h(n,n))
h=0d0

READ(10,*)
flag=0
DO WHILE(flag==0)
  READ(10,*,iostat=flag) sec,i,j,re,im
  h(table(sec,i),table(sec,j))=cmplx(re,im)
END DO
close(10)
!DO i=1,n
!  do j=1,n
!    if(abs(h(i,j))>1e-6)PRINT*,i,j,h(i,j)
!  end do
!end do
open(10,file='hmat.dat')
write(10,'(2f18.8)') h
close(10)

do nflv=1,16
  if(2**nflv==n)exit
end do
print*,'nflv=',nflv
allocate(c(n,n,nflv))
c=0
open(10,file='anniop.dat')
flag=0
do while(flag==0)
  read(10,*,iostat=flag) flv,sec,i,sec1,j,c(table(sec,i),table(sec1,j),flv)
end do
close(10)
open(10,file='cmat.dat')
write(10,'(1i8)') c
close(10)
END

