MODULE workspace
  COMPLEX(8), PARAMETER :: uniti=(0d0,1d0)
  REAL(8), ALLOCATABLE :: Vbare(:)          ! bare interactions
  INTEGER :: Nw = 100                       ! even number
  REAL(8) :: eps = 1d-4
  INTEGER :: Niter_BS = 1000
  INTEGER :: Niter_parquet = 100
  REAL(8), ALLOCATABLE :: wgrid(:)          ! (-Nw:Nw)*pi*T
  COMPLEX(8), ALLOCATABLE :: Self(:)        ! (-Nw:Nw)
  COMPLEX(8), ALLOCATABLE :: Selfnew(:)     ! (-Nw:Nw)
  COMPLEX(8), ALLOCATABLE :: G(:)           ! (-Nw:Nw)
  COMPLEX(8), ALLOCATABLE :: Xpp(:,:,:)     ! (-Nw:Nw,-Nw:Nw)
  COMPLEX(8), ALLOCATABLE :: Xph(:,:,:)     ! (-Nw:Nw,-Nw:Nw)
  COMPLEX(8), ALLOCATABLE :: Rp(:,:,:,:)    ! (-Nw:Nw,-Nw:Nw,-Nw:Nw,8), 1-4(up-up), 5-8(up-dn)
  COMPLEX(8), ALLOCATABLE :: Rc(:,:,:,:)
  COMPLEX(8), ALLOCATABLE :: Rd(:,:,:,:)
  COMPLEX(8), ALLOCATABLE :: Ip(:,:,:,:)
  COMPLEX(8), ALLOCATABLE :: Ic(:,:,:,:)
  COMPLEX(8), ALLOCATABLE :: Id(:,:,:,:)
END MODULE

SUBROUTINE BSp()
  USE workspace
  COMPLEX(8) Rnew(-Nw:Nw,-Nw:Nw,-Nw:Nw,8)
  REAL(8) :: errorbar = 1d0
  INTEGER :: iter = 0
  INTEGER i,j,k,i2
  DO WHILE(errorbar>eps.and.iter<Niter_BS)
    Rnew=0d0
    DO i=-Nw+1,Nw-1,2
      DO j=-Nw+1,Nw-1,2
        DO k=-Nw,Nw,2
          DO i2=-Nw+1,Nw-1,2
            Rnew(i,j,k,1)=Rnew(i,j,k,1)+Ip(i,i2,k,2)*Xpp(i2,k)*(Ip(i2,j,k,1)+Rp(i2,j,k,1)) &
              &                        +Ip(i,i2,k,1)*Xpp(i2,k)*(Ip(i2,j,k,2)+Rp(i2,j,k,2))
            Rnew(i,j,k,2)=Rnew(i,j,k,2)+Ip(i,i2,k,2)*Xpp(i2,k)*(Ip(i2,j,k,1)+Rp(i2,j,k,1)) &
              &                        +Ip(i,i2,k,1)*Xpp(i2,k)*(Ip(i2,j,k,2)+Rp(i2,j,k,2))
            Rnew(i,j,k,3)=Rnew(i,j,k,3)+Ip(i,i2,k,4)*Xpp(i2,k)*(Ip(i2,j,k,3)+Rp(i2,j,k,3)) &
              &                        +Ip(i,i2,k,3)*Xpp(i2,k)*(Ip(i2,j,k,4)+Rp(i2,j,k,4))
            Rnew(i,j,k,4)=Rnew(i,j,k,4)+Ip(i,i2,k,4)*Xpp(i2,k)*(Ip(i2,j,k,3)+Rp(i2,j,k,3)) &
              &                        +Ip(i,i2,k,3)*Xpp(i2,k)*(Ip(i2,j,k,4)+Rp(i2,j,k,4))
            Rnew(i,j,k,5)=Rnew(i,j,k,5)+Ip(i,i2,k,6)*Xpp(i2,k)*(Ip(i2,j,k,5)+Rp(i2,j,k,5)) &
              &                        +Ip(i,i2,k,5)*Xpp(i2,k)*(Ip(i2,j,k,6)+Rp(i2,j,k,6))
            Rnew(i,j,k,6)=Rnew(i,j,k,6)+Ip(i,i2,k,6)*Xpp(i2,k)*(Ip(i2,j,k,5)+Rp(i2,j,k,5)) &
              &                        +Ip(i,i2,k,5)*Xpp(i2,k)*(Ip(i2,j,k,6)+Rp(i2,j,k,6))
            Rnew(i,j,k,7)=Rnew(i,j,k,7)+Ip(i,i2,k,8)*Xpp(i2,k)*(Ip(i2,j,k,7)+Rp(i2,j,k,7)) &
              &                        +Ip(i,i2,k,7)*Xpp(i2,k)*(Ip(i2,j,k,8)+Rp(i2,j,k,8))
            Rnew(i,j,k,8)=Rnew(i,j,k,8)+Ip(i,i2,k,8)*Xpp(i2,k)*(Ip(i2,j,k,7)+Rp(i2,j,k,7)) &
              &                        +Ip(i,i2,k,7)*Xpp(i2,k)*(Ip(i2,j,k,8)+Rp(i2,j,k,8))
          END DO
        END DO
      END DO
    END DO
    Rnew=Rnew*T
    errorbar=maxval(abs(Rnew-Rp))
    Rp=Rp+mixing*(Rnew-Rp)
    iter=iter+1
    PRINT*,'@BSp:',iter,errorbar
  END DO
END SUBROUTINE

SUBROUTINE BSc()
  USE workspace
  COMPLEX(8) Rnew(-Nw:Nw,-Nw:Nw,-Nw:Nw,8)
  REAL(8) :: errorbar = 1d0
  INTEGER :: iter = 0
  INTEGER i,j,k,i2
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
    Rc=Rc+mixing*(Rnew-Rc)
    iter=iter+1
    PRINT*,'@BSc:',iter,errorbar
  END DO
END SUBROUTINE

SUBROUTINE BSd()
  USE workspace
  COMPLEX(8) Rnew(-Nw:Nw,-Nw:Nw,-Nw:Nw,8)
  REAL(8) :: errorbar = 1d0
  INTEGER :: iter = 0
  INTEGER i,j,k,i2
  DO WHILE(errorbar>eps.and.iter<Niter_BS)
    Rnew=0d0
    DO i=-Nw+1,Nw-1,2
      DO j=-Nw+1,Nw-1,2
        DO k=-Nw,Nw,2
          DO i2=-Nw+1,Nw-1,2
            Rnew(i,j,k,1)=Rnew(i,j,k,1)+Id(i,i2,k,1)*Xph(i2,k)*(Id(i2,j,k,1)+Rd(i2,j,k,1)) &
              &                        +Id(i,i2,k,3)*Xph(i2,k)*(Id(i2,j,k,3)+Rd(i2,j,k,3))
              &                        +Id(i,i2,k,5)*Xph(i2,k)*(Id(i2,j,k,5)+Rd(i2,j,k,5)) &
              &                        +Id(i,i2,k,7)*Xph(i2,k)*(Id(i2,j,k,7)+Rd(i2,j,k,7))
            Rnew(i,j,k,2)=Rnew(i,j,k,2)+Id(i,i2,k,2)*Xph(i2,k)*(Id(i2,j,k,4)+Rd(i2,j,k,4)) &
              &                        +Id(i,i2,k,4)*Xph(i2,k)*(Id(i2,j,k,2)+Rd(i2,j,k,2))
              &                        +Id(i,i2,k,6)*Xph(i2,k)*(Id(i2,j,k,8)+Rd(i2,j,k,8)) &
              &                        +Id(i,i2,k,8)*Xph(i2,k)*(Id(i2,j,k,6)+Rd(i2,j,k,6))
            Rnew(i,j,k,3)=Rnew(i,j,k,3)+Id(i,i2,k,1)*Xph(i2,k)*(Id(i2,j,k,3)+Rd(i2,j,k,3)) &
              &                        +Id(i,i2,k,3)*Xph(i2,k)*(Id(i2,j,k,1)+Rd(i2,j,k,1))
              &                        +Id(i,i2,k,5)*Xph(i2,k)*(Id(i2,j,k,7)+Rd(i2,j,k,7)) &
              &                        +Id(i,i2,k,7)*Xph(i2,k)*(Id(i2,j,k,5)+Rd(i2,j,k,5))
            Rnew(i,j,k,4)=Rnew(i,j,k,4)+Id(i,i2,k,2)*Xph(i2,k)*(Id(i2,j,k,2)+Rd(i2,j,k,2)) &
              &                        +Id(i,i2,k,4)*Xph(i2,k)*(Id(i2,j,k,4)+Rd(i2,j,k,4))
              &                        +Id(i,i2,k,6)*Xph(i2,k)*(Id(i2,j,k,6)+Rd(i2,j,k,6)) &
              &                        +Id(i,i2,k,8)*Xph(i2,k)*(Id(i2,j,k,8)+Rd(i2,j,k,8))
            Rnew(i,j,k,5)=Rnew(i,j,k,5)+Id(i,i2,k,5)*Xph(i2,k)*(Id(i2,j,k,1)+Rd(i2,j,k,1)) &
              &                        +Id(i,i2,k,7)*Xph(i2,k)*(Id(i2,j,k,3)+Rd(i2,j,k,3))
              &                        +Id(i,i2,k,1)*Xph(i2,k)*(Id(i2,j,k,5)+Rd(i2,j,k,5)) &
              &                        +Id(i,i2,k,3)*Xph(i2,k)*(Id(i2,j,k,7)+Rd(i2,j,k,7))
            Rnew(i,j,k,6)=Rnew(i,j,k,6)+Id(i,i2,k,6)*Xph(i2,k)*(Id(i2,j,k,4)+Rd(i2,j,k,4)) &
              &                        +Id(i,i2,k,8)*Xph(i2,k)*(Id(i2,j,k,2)+Rd(i2,j,k,2))
              &                        +Id(i,i2,k,2)*Xph(i2,k)*(Id(i2,j,k,8)+Rd(i2,j,k,8)) &
              &                        +Id(i,i2,k,4)*Xph(i2,k)*(Id(i2,j,k,6)+Rd(i2,j,k,6))
            Rnew(i,j,k,7)=Rnew(i,j,k,7)+Id(i,i2,k,5)*Xph(i2,k)*(Id(i2,j,k,3)+Rd(i2,j,k,3)) &
              &                        +Id(i,i2,k,7)*Xph(i2,k)*(Id(i2,j,k,1)+Rd(i2,j,k,1))
              &                        +Id(i,i2,k,1)*Xph(i2,k)*(Id(i2,j,k,7)+Rd(i2,j,k,7)) &
              &                        +Id(i,i2,k,3)*Xph(i2,k)*(Id(i2,j,k,5)+Rd(i2,j,k,5))
            Rnew(i,j,k,8)=Rnew(i,j,k,8)+Id(i,i2,k,6)*Xph(i2,k)*(Id(i2,j,k,2)+Rd(i2,j,k,2)) &
              &                        +Id(i,i2,k,8)*Xph(i2,k)*(Id(i2,j,k,4)+Rd(i2,j,k,4))
              &                        +Id(i,i2,k,2)*Xph(i2,k)*(Id(i2,j,k,6)+Rd(i2,j,k,6)) &
              &                        +Id(i,i2,k,4)*Xph(i2,k)*(Id(i2,j,k,8)+Rd(i2,j,k,8))
          END DO
        END DO
      END DO
    END DO
    Rnew=Rnew*T
    errorbar=maxval(abs(Rnew-Rd))
    Rd=Rd+mixing*(Rnew-Rd)
    iter=iter+1
    PRINT*,'@BSd:',iter,errorbar
  END DO
END SUBROUTINE


SUBROUTINE Dyson()
  USE workspace
  INTEGER i,j,k
  Selfnew=0d0
  DO i=-Nw+1,Nw-1,2
    DO k=-Nw,Nw,2
      j=i+k
      IF(abs(j)<=Nw) Selfnew(i)=Selfnew(i)+(Ic(i,i,k,8)+Rc(i,i,k,8) &
        &                                  +Ic(i,i,k,6)+Rc(i,i,k,6) &
        &                                  +Ic(i,i,k,4)+Rc(i,i,k,4) &
        &                                  +Ic(i,i,k,2)+Rc(i,i,k,2))*G(j)
    END DO
    DO j=-Nw+1,Nw-1,2
      Selfnew(i)=Selfnew(i)+(Ic(i,j,0,4)+Rc(i,j,0,4) &
        &                   +Ic(i,j,0,1)+Rc(i,j,0,1))*G(j)
    END DO
  END DO
  Selfnew=Selfnew*T
END SUBROUTINE

SUBROUTINE getG()
  USE workspace
  INTEGER i
  DO i=-Nw+1:Nw-1,2
    G(i)=1d0/(uniti*wgrid(i)-Self(i))
  END DO
END SUBROUTINE

SUBROUTINE getX()
  USE workspace
  INTEGER i,j,k
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

PROGRAM parquet
END PROGRAM
