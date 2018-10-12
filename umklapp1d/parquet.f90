MODULE workspace
  COMPLEX(8), PARAMETER :: uniti=(0d0,1d0)
  REAL(8), PARAMETER :: pi=acos(-1d0)
  INTEGER :: Nw = 100                       ! even number
  INTEGER :: Niter_BS = 1000
  INTEGER :: Niter_parquet = 100
  REAL(8) :: eps = 1d-4
  REAL(8) :: mixing = 0.2d0
  REAL(8) :: T = 1d0
  REAL(8) :: Vbare(8)          ! bare interactions
  INTEGER :: outputlevel
  COMPLEX(8), ALLOCATABLE :: Self(:)        ! (-Nw:Nw)
  COMPLEX(8), ALLOCATABLE :: Selfnew(:)     ! (-Nw:Nw)
  COMPLEX(8), ALLOCATABLE :: G(:)           ! (-Nw:Nw)
  COMPLEX(8), ALLOCATABLE :: Xpp(:,:)       ! (-Nw:Nw,-Nw:Nw)
  COMPLEX(8), ALLOCATABLE :: Xph(:,:)       ! (-Nw:Nw,-Nw:Nw)
  COMPLEX(8), ALLOCATABLE :: Rp(:,:,:,:)    ! (-Nw:Nw,-Nw:Nw,-Nw:Nw,8), 1-4(up-up), 5-8(up-dn)
  COMPLEX(8), ALLOCATABLE :: Rc(:,:,:,:)
  COMPLEX(8), ALLOCATABLE :: Rd(:,:,:,:)
  COMPLEX(8), ALLOCATABLE :: Ip(:,:,:,:)
  COMPLEX(8), ALLOCATABLE :: Ic(:,:,:,:)
  COMPLEX(8), ALLOCATABLE :: Id(:,:,:,:)
END MODULE

SUBROUTINE BSp()
  USE workspace
  IMPLICIT NONE
  COMPLEX(8) Rnew(-Nw:Nw,-Nw:Nw,-Nw:Nw,8)
  REAL(8) errorbar
  INTEGER iter,i,j,k,i2
  errorbar=1d0
  iter=0
  DO WHILE(errorbar>eps.and.iter<Niter_BS)
    Rnew=0d0
    DO i=-Nw+1,Nw-1,2
      DO j=-Nw+1,Nw-1,2
        DO k=-Nw,Nw,2
          DO i2=-Nw+1,Nw-1,2
            Rnew(i,j,k,1)=Rnew(i,j,k,1)+Ip(i,i2,k,2)*Xpp(i2,k)*(Ip(i2,j,k,1)+Rp(i2,j,k,1)) &
              &                        +Ip(i,i2,k,1)*Xpp(i2,k)*(Ip(i2,j,k,2)+Rp(i2,j,k,2))
            Rnew(i,j,k,2)=Rnew(i,j,k,2)+Ip(i,i2,k,2)*Xpp(i2,k)*(Ip(i2,j,k,2)+Rp(i2,j,k,2)) &
              &                        +Ip(i,i2,k,1)*Xpp(i2,k)*(Ip(i2,j,k,1)+Rp(i2,j,k,1))
            Rnew(i,j,k,3)=Rnew(i,j,k,3)+Ip(i,i2,k,4)*Xpp(i2,k)*(Ip(i2,j,k,3)+Rp(i2,j,k,3)) &
              &                        +Ip(i,i2,k,3)*Xpp(i2,k)*(Ip(i2,j,k,4)+Rp(i2,j,k,4))
            Rnew(i,j,k,4)=Rnew(i,j,k,4)+Ip(i,i2,k,4)*Xpp(i2,k)*(Ip(i2,j,k,4)+Rp(i2,j,k,4)) &
              &                        +Ip(i,i2,k,3)*Xpp(i2,k)*(Ip(i2,j,k,3)+Rp(i2,j,k,3))
            Rnew(i,j,k,5)=Rnew(i,j,k,5)+Ip(i,i2,k,6)*Xpp(i2,k)*(Ip(i2,j,k,5)+Rp(i2,j,k,5)) &
              &                        +Ip(i,i2,k,5)*Xpp(i2,k)*(Ip(i2,j,k,6)+Rp(i2,j,k,6))
            Rnew(i,j,k,6)=Rnew(i,j,k,6)+Ip(i,i2,k,6)*Xpp(i2,k)*(Ip(i2,j,k,6)+Rp(i2,j,k,6)) &
              &                        +Ip(i,i2,k,5)*Xpp(i2,k)*(Ip(i2,j,k,5)+Rp(i2,j,k,5))
            Rnew(i,j,k,7)=Rnew(i,j,k,7)+Ip(i,i2,k,8)*Xpp(i2,k)*(Ip(i2,j,k,7)+Rp(i2,j,k,7)) &
              &                        +Ip(i,i2,k,7)*Xpp(i2,k)*(Ip(i2,j,k,8)+Rp(i2,j,k,8))
            Rnew(i,j,k,8)=Rnew(i,j,k,8)+Ip(i,i2,k,8)*Xpp(i2,k)*(Ip(i2,j,k,8)+Rp(i2,j,k,8)) &
              &                        +Ip(i,i2,k,7)*Xpp(i2,k)*(Ip(i2,j,k,7)+Rp(i2,j,k,7))
          END DO
        END DO
      END DO
    END DO
    Rnew=Rnew*T
    errorbar=maxval(abs(Rnew-Rp))
    Rp=Rnew !Rp+mixing*(Rnew-Rp)
    iter=iter+1
    IF(outputlevel==0) PRINT*,'@BSp:',iter,errorbar
  END DO
END SUBROUTINE

SUBROUTINE BSc()
  USE workspace
  IMPLICIT NONE
  COMPLEX(8) Rnew(-Nw:Nw,-Nw:Nw,-Nw:Nw,8)
  REAL(8) errorbar
  INTEGER iter,i,j,k,i2
  errorbar=1d0
  iter=0
  DO WHILE(errorbar>eps.and.iter<Niter_BS)
    Rnew=0d0
    DO i=-Nw+1,Nw-1,2
      DO j=-Nw+1,Nw-1,2
        DO k=-Nw,Nw,2
          DO i2=-Nw+1,Nw-1,2
            Rnew(i,j,k,1)=Rnew(i,j,k,1)+Ic(i,i2,k,4)*Xph(i2,k)*(Ic(i2,j,k,1)+Rc(i2,j,k,1)) &
              &                        +Ic(i,i2,k,1)*Xph(i2,k)*(Ic(i2,j,k,4)+Rc(i2,j,k,4))
            Rnew(i,j,k,2)=Rnew(i,j,k,2)+Ic(i,i2,k,2)*Xph(i2,k)*(Ic(i2,j,k,2)+Rc(i2,j,k,2)) &
              &                        +Ic(i,i2,k,3)*Xph(i2,k)*(Ic(i2,j,k,3)+Rc(i2,j,k,3))
            Rnew(i,j,k,3)=Rnew(i,j,k,3)+Ic(i,i2,k,2)*Xph(i2,k)*(Ic(i2,j,k,3)+Rc(i2,j,k,3)) &
              &                        +Ic(i,i2,k,3)*Xph(i2,k)*(Ic(i2,j,k,2)+Rc(i2,j,k,2))
            Rnew(i,j,k,4)=Rnew(i,j,k,4)+Ic(i,i2,k,4)*Xph(i2,k)*(Ic(i2,j,k,4)+Rc(i2,j,k,4)) &
              &                        +Ic(i,i2,k,1)*Xph(i2,k)*(Ic(i2,j,k,1)+Rc(i2,j,k,1))
            Rnew(i,j,k,5)=Rnew(i,j,k,5)+Ic(i,i2,k,8)*Xph(i2,k)*(Ic(i2,j,k,5)+Rc(i2,j,k,5)) &
              &                        +Ic(i,i2,k,5)*Xph(i2,k)*(Ic(i2,j,k,8)+Rc(i2,j,k,8))
            Rnew(i,j,k,6)=Rnew(i,j,k,6)+Ic(i,i2,k,6)*Xph(i2,k)*(Ic(i2,j,k,6)+Rc(i2,j,k,6)) &
              &                        +Ic(i,i2,k,7)*Xph(i2,k)*(Ic(i2,j,k,7)+Rc(i2,j,k,7))
            Rnew(i,j,k,7)=Rnew(i,j,k,7)+Ic(i,i2,k,6)*Xph(i2,k)*(Ic(i2,j,k,7)+Rc(i2,j,k,7)) &
              &                        +Ic(i,i2,k,7)*Xph(i2,k)*(Ic(i2,j,k,6)+Rc(i2,j,k,6))
            Rnew(i,j,k,8)=Rnew(i,j,k,8)+Ic(i,i2,k,8)*Xph(i2,k)*(Ic(i2,j,k,8)+Rc(i2,j,k,8)) &
              &                        +Ic(i,i2,k,5)*Xph(i2,k)*(Ic(i2,j,k,5)+Rc(i2,j,k,5))
          END DO
        END DO
      END DO
    END DO
    Rnew=Rnew*T
    errorbar=maxval(abs(Rnew-Rc))
    Rc=Rnew !Rc+mixing*(Rnew-Rc)
    iter=iter+1
    IF(outputlevel==0) PRINT*,'@BSc:',iter,errorbar
  END DO
END SUBROUTINE

SUBROUTINE BSd()
  USE workspace
  IMPLICIT NONE
  COMPLEX(8) Rnew(-Nw:Nw,-Nw:Nw,-Nw:Nw,8)
  REAL(8) errorbar
  INTEGER iter,i,j,k,i2
  errorbar=1d0
  iter=0
  DO WHILE(errorbar>eps.and.iter<Niter_BS)
    Rnew=0d0
    DO i=-Nw+1,Nw-1,2
      DO j=-Nw+1,Nw-1,2
        DO k=-Nw,Nw,2
          DO i2=-Nw+1,Nw-1,2
            Rnew(i,j,k,1)=Rnew(i,j,k,1)+Id(i,i2,k,1)*Xph(i2,k)*(Id(i2,j,k,1)+Rd(i2,j,k,1)) &
              &                        +Id(i,i2,k,3)*Xph(i2,k)*(Id(i2,j,k,3)+Rd(i2,j,k,3)) &
              &                        +Id(i,i2,k,5)*Xph(i2,k)*(Id(i2,j,k,5)+Rd(i2,j,k,5)) &
              &                        +Id(i,i2,k,7)*Xph(i2,k)*(Id(i2,j,k,7)+Rd(i2,j,k,7))
            Rnew(i,j,k,2)=Rnew(i,j,k,2)+Id(i,i2,k,2)*Xph(i2,k)*(Id(i2,j,k,4)+Rd(i2,j,k,4)) &
              &                        +Id(i,i2,k,4)*Xph(i2,k)*(Id(i2,j,k,2)+Rd(i2,j,k,2)) &
              &                        +Id(i,i2,k,6)*Xph(i2,k)*(Id(i2,j,k,8)+Rd(i2,j,k,8)) &
              &                        +Id(i,i2,k,8)*Xph(i2,k)*(Id(i2,j,k,6)+Rd(i2,j,k,6))
            Rnew(i,j,k,3)=Rnew(i,j,k,3)+Id(i,i2,k,1)*Xph(i2,k)*(Id(i2,j,k,3)+Rd(i2,j,k,3)) &
              &                        +Id(i,i2,k,3)*Xph(i2,k)*(Id(i2,j,k,1)+Rd(i2,j,k,1)) &
              &                        +Id(i,i2,k,5)*Xph(i2,k)*(Id(i2,j,k,7)+Rd(i2,j,k,7)) &
              &                        +Id(i,i2,k,7)*Xph(i2,k)*(Id(i2,j,k,5)+Rd(i2,j,k,5))
            Rnew(i,j,k,4)=Rnew(i,j,k,4)+Id(i,i2,k,2)*Xph(i2,k)*(Id(i2,j,k,2)+Rd(i2,j,k,2)) &
              &                        +Id(i,i2,k,4)*Xph(i2,k)*(Id(i2,j,k,4)+Rd(i2,j,k,4)) &
              &                        +Id(i,i2,k,6)*Xph(i2,k)*(Id(i2,j,k,6)+Rd(i2,j,k,6)) &
              &                        +Id(i,i2,k,8)*Xph(i2,k)*(Id(i2,j,k,8)+Rd(i2,j,k,8))
            Rnew(i,j,k,5)=Rnew(i,j,k,5)+Id(i,i2,k,5)*Xph(i2,k)*(Id(i2,j,k,1)+Rd(i2,j,k,1)) &
              &                        +Id(i,i2,k,7)*Xph(i2,k)*(Id(i2,j,k,3)+Rd(i2,j,k,3)) &
              &                        +Id(i,i2,k,1)*Xph(i2,k)*(Id(i2,j,k,5)+Rd(i2,j,k,5)) &
              &                        +Id(i,i2,k,3)*Xph(i2,k)*(Id(i2,j,k,7)+Rd(i2,j,k,7))
            Rnew(i,j,k,6)=Rnew(i,j,k,6)+Id(i,i2,k,6)*Xph(i2,k)*(Id(i2,j,k,4)+Rd(i2,j,k,4)) &
              &                        +Id(i,i2,k,8)*Xph(i2,k)*(Id(i2,j,k,2)+Rd(i2,j,k,2)) &
              &                        +Id(i,i2,k,2)*Xph(i2,k)*(Id(i2,j,k,8)+Rd(i2,j,k,8)) &
              &                        +Id(i,i2,k,4)*Xph(i2,k)*(Id(i2,j,k,6)+Rd(i2,j,k,6))
            Rnew(i,j,k,7)=Rnew(i,j,k,7)+Id(i,i2,k,5)*Xph(i2,k)*(Id(i2,j,k,3)+Rd(i2,j,k,3)) &
              &                        +Id(i,i2,k,7)*Xph(i2,k)*(Id(i2,j,k,1)+Rd(i2,j,k,1)) &
              &                        +Id(i,i2,k,1)*Xph(i2,k)*(Id(i2,j,k,7)+Rd(i2,j,k,7)) &
              &                        +Id(i,i2,k,3)*Xph(i2,k)*(Id(i2,j,k,5)+Rd(i2,j,k,5))
            Rnew(i,j,k,8)=Rnew(i,j,k,8)+Id(i,i2,k,6)*Xph(i2,k)*(Id(i2,j,k,2)+Rd(i2,j,k,2)) &
              &                        +Id(i,i2,k,8)*Xph(i2,k)*(Id(i2,j,k,4)+Rd(i2,j,k,4)) &
              &                        +Id(i,i2,k,2)*Xph(i2,k)*(Id(i2,j,k,6)+Rd(i2,j,k,6)) &
              &                        +Id(i,i2,k,4)*Xph(i2,k)*(Id(i2,j,k,8)+Rd(i2,j,k,8))
          END DO
        END DO
      END DO
    END DO
    Rnew=-Rnew*T  ! "-" comes from the fermion loop
    errorbar=maxval(abs(Rnew-Rd))
    Rd=Rnew !Rd+mixing*(Rnew-Rd)
    iter=iter+1
    IF(outputlevel==0) PRINT*,'@BSd:',iter,errorbar
  END DO
END SUBROUTINE

SUBROUTINE Dyson()
  USE workspace
  IMPLICIT NONE
  INTEGER i,j,k,i2
  COMPLEX(8) Ga(-Nw:Nw,-Nw:Nw,-Nw:Nw,8)  ! 1PI
  Ga=0d0
  DO i=-Nw+1,Nw-1,2
    DO j=-Nw+1,Nw-1,2
      DO k=-Nw,Nw,2
        Ga(i,j,k,:)=-Vbare(:)
        DO i2=-Nw+1,Nw-1,2
          Ga(i,j,k,1)=Ga(i,j,k,1)-Vbare(2)*Xpp(i2,k)*(Ip(i2,j,k,1)+Rp(i2,j,k,1))*T &
            &                    -Vbare(1)*Xpp(i2,k)*(Ip(i2,j,k,2)+Rp(i2,j,k,2))*T
          Ga(i,j,k,2)=Ga(i,j,k,2)-Vbare(2)*Xpp(i2,k)*(Ip(i2,j,k,2)+Rp(i2,j,k,2))*T &
              &                  -Vbare(1)*Xpp(i2,k)*(Ip(i2,j,k,1)+Rp(i2,j,k,1))*T
          Ga(i,j,k,3)=Ga(i,j,k,3)-Vbare(4)*Xpp(i2,k)*(Ip(i2,j,k,3)+Rp(i2,j,k,3))*T &
              &                  -Vbare(3)*Xpp(i2,k)*(Ip(i2,j,k,4)+Rp(i2,j,k,4))*T
          Ga(i,j,k,4)=Ga(i,j,k,4)-Vbare(4)*Xpp(i2,k)*(Ip(i2,j,k,4)+Rp(i2,j,k,4))*T &
              &                  -Vbare(3)*Xpp(i2,k)*(Ip(i2,j,k,3)+Rp(i2,j,k,3))*T
          Ga(i,j,k,5)=Ga(i,j,k,5)-Vbare(6)*Xpp(i2,k)*(Ip(i2,j,k,5)+Rp(i2,j,k,5))*T &
              &                  -Vbare(5)*Xpp(i2,k)*(Ip(i2,j,k,6)+Rp(i2,j,k,6))*T
          Ga(i,j,k,6)=Ga(i,j,k,6)-Vbare(6)*Xpp(i2,k)*(Ip(i2,j,k,6)+Rp(i2,j,k,6))*T &
              &                  -Vbare(5)*Xpp(i2,k)*(Ip(i2,j,k,5)+Rp(i2,j,k,5))*T
          Ga(i,j,k,7)=Ga(i,j,k,7)-Vbare(8)*Xpp(i2,k)*(Ip(i2,j,k,7)+Rp(i2,j,k,7))*T &
              &                  -Vbare(7)*Xpp(i2,k)*(Ip(i2,j,k,8)+Rp(i2,j,k,8))*T
          Ga(i,j,k,8)=Ga(i,j,k,8)-Vbare(8)*Xpp(i2,k)*(Ip(i2,j,k,8)+Rp(i2,j,k,8))*T &
              &                  -Vbare(7)*Xpp(i2,k)*(Ip(i2,j,k,7)+Rp(i2,j,k,7))*T
        END DO
      END DO
    END DO
  END DO
  Selfnew=0d0
  DO i=-Nw+1,Nw-1,2
    DO j=-Nw+1,Nw-1,2
      k=i+j
      IF(abs(k)<=Nw) Selfnew(i)=Selfnew(i)+(Ga(i,j,k,1)+Ga(i,j,k,4) &
        &                                  -Ga(i,i,k,8)-Ga(i,i,k,6) &
        &                                  -Ga(i,i,k,4)-Ga(i,i,k,2))*G(j)
    END DO
  END DO
  Selfnew=Selfnew*T
END SUBROUTINE

SUBROUTINE getX()
  USE workspace
  IMPLICIT NONE
  INTEGER i,j,k
  DO i=-Nw+1,Nw-1,2
    G(i)=1d0/(uniti*i*pi*T-Self(i))
  END DO
  Xpp=0d0
  Xph=0d0
  DO i=-Nw+1,Nw-1,2
    DO k=-Nw,Nw,2
      j=k-i
      IF(abs(j)<=Nw) Xpp(i,k)=G(i)*G(j)
      j=k+i
      IF(abs(j)<=Nw) Xph(i,k)=G(i)*G(j)
    END DO
  END DO
END SUBROUTINE

SUBROUTINE projector()
  USE workspace
  IMPLICIT NONE
  INTEGER i,j,k,i2,j2,k2
  
  Ip=0d0
  Ic=0d0
  Id=0d0
  DO i=-Nw+1,Nw-1,2
    DO j=-Nw+1,Nw-1,2
      DO k=-Nw,Nw,2
        Ip(i,j,k,:)=-Vbare(:)
        Ic(i,j,k,:)=-Vbare(:)
        Id(i,j,k,:)=-Vbare(:)
      END DO
    END DO
  END DO

  DO i=-Nw+1,Nw-1,2
    DO j=-Nw+1,Nw-1,2
      DO k=-Nw,Nw,2

        i2=i
        j2=j
        k2=k-i-j
        IF(abs(k2)<=Nw) Ic(i2,j2,k2,:)=Ic(i2,j2,k2,:)+Rp(i,j,k,:)

        i2=i
        j2=k-j
        k2=j-i
        IF(abs(j2)<=Nw.and.abs(k2)<=Nw) Id(i2,j2,k2,:)=Id(i2,j2,k2,:)+Rp(i,j,k,:)

        i2=i
        j2=j
        k2=k+i+j
        IF(abs(k2)<=Nw) Ip(i2,j2,k2,:)=Ip(i2,j2,k2,:)+Rc(i,j,k,:)

        i2=i
        j2=k+i
        k2=j-i
        IF(abs(j2)<=Nw.and.abs(k2)<=Nw) Id(i2,j2,k2,:)=Id(i2,j2,k2,:)+Rc(i,j,k,:)

        i2=i
        j2=k+i
        k2=k+i+j
        IF(abs(j2)<=Nw.and.abs(k2)<=Nw) Ip(i2,j2,k2,:)=Ip(i2,j2,k2,:)+Rd(i,j,k,:)

        i2=i
        j2=k+i
        k2=j-i
        IF(abs(j2)<=Nw.and.abs(k2)<=Nw) Ic(i2,j2,k2,:)=Ic(i2,j2,k2,:)+Rd(i,j,k,:)

      END DO
    END DO
  END DO

END SUBROUTINE

SUBROUTINE init()
  USE workspace
  IMPLICIT NONE
  OPEN(10,FILE='parquet.in')
  READ(10,*) T
  Vbare=0d0
  READ(10,*) Vbare(7)   !umklapp, all others are zero.
  READ(10,*) Nw
  READ(10,*) Niter_BS
  READ(10,*) Niter_parquet
  READ(10,*) eps
  READ(10,*) mixing
  READ(10,*) outputlevel
  CLOSE(10)
  ALLOCATE(Self(-Nw:Nw),Selfnew(-Nw:Nw),G(-Nw:Nw))
  ALLOCATE(Xpp(-Nw:Nw,-Nw:Nw))
  ALLOCATE(Xph(-Nw:Nw,-Nw:Nw))
  ALLOCATE(Rp(-Nw:Nw,-Nw:Nw,-Nw:Nw,8))
  ALLOCATE(Rc(-Nw:Nw,-Nw:Nw,-Nw:Nw,8))
  ALLOCATE(Rd(-Nw:Nw,-Nw:Nw,-Nw:Nw,8))
  ALLOCATE(Ip(-Nw:Nw,-Nw:Nw,-Nw:Nw,8))
  ALLOCATE(Ic(-Nw:Nw,-Nw:Nw,-Nw:Nw,8))
  ALLOCATE(Id(-Nw:Nw,-Nw:Nw,-Nw:Nw,8))
  Self=0d0
  Rp=0d0
  Rc=0d0
  Rd=0d0
END SUBROUTINE

SUBROUTINE output()
  USE workspace
  IMPLICIT NONE
  INTEGER i,j,k,loc(4)
  OPEN(10,FILE='self.dat')
  DO i=1,Nw,2
    WRITE(10,'(3f12.6)') i*pi*T,Self(i)
  END DO
  CLOSE(10)
  OPEN(10,FILE='Ipfull.dat')
  OPEN(11,FILE='Icfull.dat')
  OPEN(12,FILE='Idfull.dat')
  OPEN(20,FILE='Imax.dat')
  DO i=-Nw+1,Nw-1,2
    DO j=-Nw+1,Nw-1,2
      DO k=-Nw,Nw,2
        WRITE(10,'(19f12.6)') i*pi*T,j*pi*T,k*pi*T,Ip(i,j,k,:)
        WRITE(11,'(19f12.6)') i*pi*T,j*pi*T,k*pi*T,Ic(i,j,k,:)
        WRITE(12,'(19f12.6)') i*pi*T,j*pi*T,k*pi*T,Id(i,j,k,:)
      END DO
    END DO
  END DO
  loc=maxloc(abs(Ip))
  loc(1:3)=loc(1:3)-Nw-1
  WRITE(20,'(1x,a,3f12.6,1i6,2f12.6)') 'P:',loc(1:3)*pi*T,loc(4),Ip(loc(1),loc(2),loc(3),loc(4))
  loc=maxloc(abs(Ic))
  loc(1:3)=loc(1:3)-Nw-1
  WRITE(20,'(1x,a,3f12.6,1i6,2f12.6)') 'C:',loc(1:3)*pi*T,loc(4),Ic(loc(1),loc(2),loc(3),loc(4))
  loc=maxloc(abs(Id))
  loc(1:3)=loc(1:3)-Nw-1
  WRITE(20,'(1x,a,3f12.6,1i6,2f12.6)') 'D:',loc(1:3)*pi*T,loc(4),Id(loc(1),loc(2),loc(3),loc(4))
  CLOSE(10)
  CLOSE(11)
  CLOSE(12)
  CLOSE(20)
END SUBROUTINE

SUBROUTINE exact()
  USE workspace
  IMPLICIT NONE
  INTEGER i
  REAL(8) U,Z
  COMPLEX(8) iw,Gexact
  U=Vbare(7)
  Z=exp(U/T)+exp(-U/T)+14d0
  OPEN(10,FILE='selfexact.dat')
  DO i=1,Nw-1,2
    iw=uniti*i*pi*T
    Gexact=((1d0/(iw-U)+1d0/(iw+U))*(2+exp(-U/T)+exp(U/T))/2+12d0/iw)/Z
    WRITE(10,'(3f12.6)') i*pi*T,iw-1d0/Gexact
  END DO
  CLOSE(10)
END SUBROUTINE

PROGRAM parquet
  USE workspace
  IMPLICIT NONE
  INTEGER iter
  REAL(8) errorbar
  CALL init()
  iter=0
  errorbar=1d0
  DO WHILE(errorbar>eps.and.iter<Niter_parquet)
    CALL getX()
    CALL projector()
    CALL BSp()
    CALL BSc()
    CALL BSd()
    CALL projector()
    CALL Dyson()
    errorbar=maxval(abs(Selfnew-Self))
    Self=Self+mixing*(Selfnew-Self)
    iter=iter+1
    PRINT'(1a,1i6,3f12.6)','@parquet:',iter,errorbar,self(1)
  END DO
  CALL output()
  CALL exact()
END PROGRAM
