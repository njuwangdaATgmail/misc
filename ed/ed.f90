MODULE mod_ed
! modified 5/2/2013
 
  INTEGER,PARAMETER :: dp=kind(1d0)
  
  TYPE complex_mat
    COMPLEX(dp), ALLOCATABLE :: mat(:,:)
  END TYPE
 
  TYPE subh
    INTEGER nsize
    INTEGER, ALLOCATABLE :: base(:)
    REAL(dp), ALLOCATABLE :: en(:)
    COMPLEX(dp), ALLOCATABLE :: mat(:,:)
  END TYPE


  LOGICAL :: calc_spec=.false.
  INTEGER nflavor    ! include site, spin, orbit... 
  REAL(dp) eps
  INTEGER nhilbert   ! =2**nflavor
  INTEGER nsector    ! number of sector
  INTEGER, ALLOCATABLE :: state(:,:)   ! size of nflavor-by-nhilbert, stores 0 and 1 to mark the state
  INTEGER, ALLOCATABLE :: sector(:)    ! size of nhilbert, point to the sector # of each state
  COMPLEX(dp), ALLOCATABLE :: hop(:,:)
  COMPLEX(dp), ALLOCATABLE :: v(:,:,:,:)
  TYPE(complex_mat), ALLOCATABLE :: psi(:,:) ! size of nsector-by-nflavor
  TYPE(subh), ALLOCATABLE :: h(:)  ! size of nsector
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! declare some private subroutines/functions.
  PRIVATE state2idx,evolve,zeigen  

CONTAINS

  FUNCTION state2idx(sta)
    IMPLICIT NONE
    INTEGER sta(nflavor),i,state2idx
    state2idx=1
    DO i=1,nflavor
      state2idx=state2idx+sta(i)*2**(i-1)
    END DO
  END FUNCTION

  SUBROUTINE checkcomment(fp)
    IMPLICIT NONE
    INTEGER fp
    CHARACTER str
    DO WHILE(.true.)
      READ(fp,*) str
      IF(str/='#'.AND.str/='!')THEN
        BACKSPACE(fp)
        RETURN
      END IF
    END DO
  END SUBROUTINE

  SUBROUTINE ed()
    IMPLICIT NONE
    INTEGER i,j,k,l,sec,info,nhop,nv
    REAL(8) re,im,emin,emin2
    
    eps=1d-8

    OPEN(10,FILE='ed.in')
    READ(10,*) nflavor,nhop,nv
    nhilbert=2**nflavor 
    ALLOCATE(state(nflavor,nhilbert),hop(nflavor,nflavor),v(nflavor,nflavor,nflavor,nflavor),sector(nhilbert))
    DO sec=1,nhop
      CALL checkcomment(10)
      READ(10,*) i,j,re,im
      hop(i,j)=hop(i,j)+cmplx(re,im)
    END DO
    DO sec=1,nv
      CALL checkcomment(10)
      READ(10,*) i,j,k,l,re,im
      v(i,j,k,l)=v(i,j,k,l)+cmplx(re,im)
    END DO
    CLOSE(10)
    
    ! set value to state(:,:)
    state(:,:)=0
    DO i=1,nhilbert
      k=i-1
      j=1
      DO WHILE(.true.)
        state(j,i)=mod(k,2)
        k=(k-state(j,i))/2
        j=j+1
        IF(k==0)EXIT
      END DO
    END DO
    
    sector(:)=0  ! initial sector, sector(i)=0 means state i dose not belong to any sector.
    DO i=1,nhilbert
      sec=sector(i)
      IF(sec==0)THEN
        sec=maxval(sector(1:i))+1
        CALL evolve(i,sec,.false.)   ! find all state k belong this sector, their sector(k)=sec
      END IF
    END DO
    nsector=maxval(sector)
    
    ! set hamiltonian of each sector
    ALLOCATE(h(nsector))
    sector(:)=0
    DO i=1,nhilbert
      sec=sector(i)
      IF(sec==0)THEN
        sec=maxval(sector(1:i))+1
        CALL evolve(i,sec,.true.)
      END IF
    END DO

    OPEN(10,FILE='ham.dat')
    WRITE(10,*) nsector,'    :nsector'
    DO sec=1,nsector
      WRITE(10,'(2i6)') sec,h(sec)%nsize
    END DO
    WRITE(10,*) '--------------------------'
    DO sec=1,nsector
      DO i=1,h(sec)%nsize
        DO j=1,h(sec)%nsize
          IF(abs(h(sec)%mat(i,j))<1d-6)CYCLE
          WRITE(10,'(3i6,2f15.8)')sec,i,j,real(h(sec)%mat(i,j)),aimag(h(sec)%mat(i,j))
        END DO
      END DO
    END DO
    CLOSE(10)

    OPEN(10,FILE='basis.dat')
    DO sec=1,nsector
      WRITE(10,*) 'sec=',sec,':'
      WRITE(10,'(1I12)') h(sec)%base
      WRITE(10,*) '---------------------------'
    END DO
    CLOSE(10)

    !PRINT*,nsector
    OPEN(10,FILE='en.dat')
    DO i=1,nsector
      !print*,h(i)%nsize
      !print'(8f8.3)',h(i)%mat
      !WRITE(10,*) 'sec=',i
      CALL zeigen(h(i)%nsize,h(i)%mat,h(i)%en)
      WRITE(10,'(1f12.6)')h(i)%en
      !WRITE(10,*) '----------------'
    end do

    !========================================
    !   output ground state wave function 
    !========================================
    emin=100000
    i=0
    DO sec=1,nsector
      IF(emin>h(sec)%en(1))THEN
        emin=h(sec)%en(1)
        i=sec
      END IF
    END DO
    print*,'emin=',emin
    
    emin2=100000
    DO sec=1,nsector
      IF(emin2>h(sec)%en(1).AND.abs(h(sec)%en(1)-emin)>1e-6)THEN
        emin2=h(sec)%en(1)
      END IF
    END DO
    print*,'gap=',emin2-emin
    
    OPEN(10,FILE='gswf.dat')
    DO j=1,h(i)%nsize
      WRITE(10,*) h(i)%base(j),real(h(i)%mat(j,1)),aimag(h(i)%mat(j,1))
    END DO
    CLOSE(10)
 


    CALL loc_operator() 

  END SUBROUTINE

  SUBROUTINE evolve(idx0_,sec,seth)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: idx0_,sec
    LOGICAL, INTENT(IN) :: seth
    INTEGER i,j,k,l,n,idx0,idx,ipiv,ii,jj
    LOGICAL gameover
    INTEGER x0(nflavor),x(nflavor)
    
    idx0=idx0_
    x0=state(:,idx0)
    sector(idx0)=sec
    n=1
    gameover=.false.
   
    DO WHILE(.not.gameover)

      DO i=1,nflavor
        DO j=1,nflavor
          IF(abs(hop(i,j))<eps)CYCLE
          x=x0
          IF(x(j)/=1)CYCLE
          x(j)=0
          IF(x(i)/=0)CYCLE
          x(i)=1
          idx=state2idx(x)
          IF(sector(idx)==0)THEN
            n=n+1
            sector(idx)=sec
          END IF
          !print*,sec,i,j,n
        END DO
      END DO

      DO i=1,nflavor
        DO j=1,nflavor
          DO k=1,nflavor
            DO l=1,nflavor
            IF(abs(v(i,j,k,l))<eps)CYCLE
              x=x0
              IF(x(l)/=1)CYCLE
              x(l)=0
              IF(x(k)/=1)CYCLE
              x(k)=0
              IF(x(j)/=0)CYCLE
              x(j)=1
              IF(x(i)/=0)CYCLE
              x(i)=1
              idx=state2idx(x)
              IF(sector(idx)==0)THEN
                n=n+1
                sector(idx)=sec
              END IF
            END DO
          END DO
        END DO
      END DO

      sector(idx0)=-1    ! -1 means I have completed the search starting from idx0
      gameover=.true.
      DO i=1,nhilbert
        IF(sector(i)==sec)THEN
          idx0=i
          x0=state(:,idx0)
          gameover=.false.
          EXIT
        END IF
      END DO
    
    END DO  ! DO WHILE(.not.gameover)


    IF(.not.seth)THEN
      WHERE(sector==-1)
        sector=sec
      END WHERE
      RETURN
    ELSE
      WHERE(sector==-1)
        sector=0
      END WHERE
    END IF
    
    
    h(sec)%nsize=n
    ALLOCATE(h(sec)%base(n),h(sec)%mat(n,n),h(sec)%en(n))
    h(sec)%mat(:,:)=0d0
    
    idx0=idx0_
    x0=state(:,idx0)
    sector(idx0)=sec
    n=1
    jj=1
    h(sec)%base(1)=idx0
    gameover=.false.

    DO WHILE(.not.gameover)
      
      DO i=1,nflavor
        DO j=1,nflavor
          IF(abs(hop(i,j))<eps)CYCLE
          x=x0
          IF(x(j)/=1)CYCLE
          ipiv=sum(x(1:j-1))
          x(j)=0
          IF(x(i)/=0)CYCLE
          ipiv=ipiv+sum(x(1:i-1))
          x(i)=1
          idx=state2idx(x)
          IF(sector(idx)==0)THEN
            n=n+1
            h(sec)%base(n)=idx
            sector(idx)=sec
          END IF
         
          DO ii=1,n   ! if idx has been found before, we need find its ii
            IF(h(sec)%base(ii)==idx)EXIT
          END DO
          IF(mod(ipiv,2)==0)THEN
            h(sec)%mat(ii,jj)=h(sec)%mat(ii,jj)+hop(i,j)
          ELSE
            h(sec)%mat(ii,jj)=h(sec)%mat(ii,jj)-hop(i,j)
          END IF
        END DO
      END DO

      DO i=1,nflavor
        DO j=1,nflavor
          DO k=1,nflavor
            DO l=1,nflavor
              IF(abs(v(i,j,k,l))<eps)CYCLE
              x=x0
              IF(x(l)/=1)CYCLE
              ipiv=sum(x(1:l-1))
              x(l)=0
              IF(x(k)/=1)CYCLE
              ipiv=ipiv+sum(x(1:k-1))
              x(k)=0
              IF(x(j)/=0)CYCLE
              ipiv=ipiv+sum(x(1:j-1))
              x(j)=1
              IF(x(i)/=0)CYCLE
              ipiv=ipiv+sum(x(1:i-1))
              x(i)=1
              idx=state2idx(x)
              IF(sector(idx)==0)THEN
                n=n+1
                h(sec)%base(n)=idx
                sector(idx)=sec
              END IF
              DO ii=1,n
                IF(h(sec)%base(ii)==idx)EXIT
              END DO
              IF(mod(ipiv,2)==0)THEN
                h(sec)%mat(ii,jj)=h(sec)%mat(ii,jj)+v(i,j,k,l)
              ELSE
                h(sec)%mat(ii,jj)=h(sec)%mat(ii,jj)-v(i,j,k,l)
              END IF
            END DO
          END DO
        END DO
      END DO

      sector(idx0)=-1
      gameover=.true.
      DO i=1,nhilbert
        IF(sector(i)==sec)THEN
          idx0=i
          x0=state(:,idx0)
          gameover=.false.
          EXIT
        END IF
      END DO
      DO jj=1,n
        IF(h(sec)%base(jj)==idx0)EXIT
      END DO

    END DO

    WHERE(sector==-1)
      sector=sec
    END WHERE


  END SUBROUTINE

  SUBROUTINE loc_operator()
    IMPLICIT NONE
    INTEGER flv,sec,sec1,d0,d1,n,m,x0(nflavor),sgn,idx1
    
    OPEN(10,FILE='anniop.dat')
    DO flv=1,nflavor

      DO sec=1,nsector
        d0=h(sec)%nsize
        
        DO sec1=1,nsector
          d1=h(sec1)%nsize

          DO n=1,d0

            x0=state(:,h(sec)%base(n))
            IF(x0(flv)==1)CYCLE
            sgn=sum(x0(1:flv-1))
            x0(flv)=1
            idx1=state2idx(x0)
            DO m=1,d1
              IF(h(sec1)%base(m)==idx1)EXIT
            END DO
            IF(m==d1+1)CYCLE
            
            WRITE(10,'(6i6)') flv,sec,n,sec1,m,(-1)**sgn

          END DO

        END DO

      END DO
    END DO

    CLOSE(10)
  END SUBROUTINE

            
#ifdef USELESS
! used for hyb-ctqmc

  SUBROUTINE loc_operator()
    IMPLICIT NONE
    INTEGER flv,sec,sec1,i,j,n,m,d0,d1,idx1,sgn,x0(nflavor),nw
    REAL(dp) psi,emin,Z,eta,wmin,wmax,beta
    REAL(dp),ALLOCATABLE :: w(:),spec(:,:)
    
    IF(calc_spec)THEN
      OPEN(30,FILE='spec.in')
      READ(30,*) nw
      READ(30,*) wmin,wmax
      READ(30,*) beta
      READ(30,*) eta
      CLOSE(30)
      
      ALLOCATE(w(nw),spec(nw,nflavor))
      spec=0d0
      DO i=1,nw
        w(i)=wmin+(i-0.5d0)*(wmax-wmin)/nw
      END DO
    END IF 
    

    OPEN(20,FILE='local.dat')
    WRITE(20,*) 'nsector=',nsector
    WRITE(20,*) 'begin enloc'
    WRITE(20,*) '#sec n <sec,n|H|sec,n>'
    emin=100000
    DO sec=1,nsector
      emin=min(emin,minval(h(sec)%en))
    END DO
    DO sec=1,nsector
      DO n=1,h(sec)%nsize
        WRITE(20,'(2i4,1f20.10)') sec,n,h(sec)%en(n)-emin
      END DO
    END DO
    WRITE(20,*) 'end enloc'
    WRITE(20,*)
    WRITE(20,*) 'begin psiloc'
    WRITE(20,*) '#sec0 n0 flv sec1 n1 <sec0,n0|flv|sec1,n1>'
    DO flv=1,nflavor
      DO sec=1,nsector
        d0=h(sec)%nsize
        DO sec1=1,nsector
          d1=h(sec1)%nsize
          DO i=1,d0
            DO j=1,d1
              ! <sec,i|c|sec1,j> = <sec,i|sec,n><sec,n|c|sec1,m><sec1,m|sec1,j>
              psi=0d0
              DO n=1,d0
                x0=state(:,h(sec)%base(n))
                IF(x0(flv)==1)CYCLE
                sgn=sum(x0(1:flv-1))
                x0(flv)=1
                idx1=state2idx(x0)
                DO m=1,d1
                  IF(h(sec1)%base(m)==idx1)EXIT
                END DO
                IF(m==d1+1)CYCLE
                IF(mod(sgn,2)==0)THEN
                  psi=psi+h(sec)%mat(n,i)*h(sec1)%mat(m,j)
                ELSE
                  psi=psi-h(sec)%mat(n,i)*h(sec1)%mat(m,j)
                END IF
              END DO
              
              IF(abs(psi)>eps)WRITE(20,'(5i4,1f20.10)') sec,i,flv,sec1,j,psi
              
              IF(abs(psi)>eps.AND.calc_spec)THEN
                z=psi**2*(exp(-beta*(h(sec)%en(i)-emin))+exp(-beta*(h(sec1)%en(j)-emin)))  
                spec(:,flv)=spec(:,flv)+z/((w+h(sec)%en(i)-h(sec1)%en(j))**2+eta**2)      

              END IF                

            END DO
          END DO
        END DO
      END DO
    END DO
    WRITE(20,*) 'end psiloc'
    CLOSE(20)
    
    
    IF(calc_spec)THEN
      spec=spec/sum(spec*(wmax-wmin)/nw)
      OPEN(30,FILE='spec.dat')
      DO i=1,nw
        WRITE(30,'(1f12.6)',advance='no')w(i)
        DO j=1,nflavor
          WRITE(30,'(1f12.6)',advance = 'no')spec(i,j)
        END DO
        WRITE(30,*)
      END DO
      CLOSE(30)
    END IF
  
  END SUBROUTINE

#endif
        
  SUBROUTINE zeigen(n,a,v)
    IMPLICIT NONE
    INTEGER n,info
    REAL(8) rwork(3*n),v(n)
    COMPLEX(8) work(10*n),a(n,n)
    IF(n<=0)STOP 'N<=0 is invalid, inside DEIGEN'
    IF(n==1)THEN
      v(1)=a(1,1);a(1,1)=1d0;RETURN
    END IF
    CALL zheev('v','u',n,a,n,v,work,10*n,rwork,info)
    IF(info/=0)STOP 'ERROR @ ZHEEV, inside ZEIGEN'
  END SUBROUTINE

      
  

  
END MODULE

PROGRAM main
  USE mod_ed
  CALL ed()
END 
