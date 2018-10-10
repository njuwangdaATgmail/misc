MODULE workspace
  COMPLEX(8), ALLOCATABLE :: Self(:)        ! (Nw)
  COMPLEX(8), ALLOCATABLE :: Xpp(:,:,:)     ! (Nw,Nv)
  COMPLEX(8), ALLOCATABLE :: Rp(:,:,:,:)    ! (Nw,Nw,Nv,8), 1-4(up-up), 5-8(up-dn)
  COMPLEX(8), ALLOCATABLE :: Rc(:,:,:,:)    ! (Nw,Nw,Nv,8), 1-4(up-up), 5-8(up-dn)
  COMPLEX(8), ALLOCATABLE :: Rd(:,:,:,:)    ! (Nw,Nw,Nv,8), 1-4(up-up), 5-8(up-dn)
  COMPLEX(8), ALLOCATABLE :: Ip(:,:,:,:)    ! (Nw,Nw,Nv,8), 1-4(up-up), 5-8(up-dn)
  COMPLEX(8), ALLOCATABLE :: Ic(:,:,:,:)    ! (Nw,Nw,Nv,8), 1-4(up-up), 5-8(up-dn)
  COMPLEX(8), ALLOCATABLE :: Id(:,:,:,:)    ! (Nw,Nw,Nv,8), 1-4(up-up), 5-8(up-dn)
END MODULE

SUBROUTINE BS()
  USE workspace
END SUBROUTINE

SUBROUTINE Dyson()
  USE workspace
END SUBROUTINE

SUBROUTINE p2c()
  USE workspace
END SUBROUTINE

SUBROUTINE p2d()
  USE workspace
END SUBROUTINE

SUBROUTINE c2p()
  USE workspace
END SUBROUTINE

SUBROUTINE c2d()
  USE workspace
END SUBROUTINE

SUBROUTINE d2p()
  USE workspace
END SUBROUTINE

SUBROUTINE d2c()
  USE workspace
END SUBROUTINE

PROGRAM parquet
END PROGRAM
