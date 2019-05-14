program lanczos_test
  USE COMMON_VARS
  USE MATVEC_PRODUCT
  USE SPARSE_MATRIX
  USE LANCZOS
  USE SETUP
  USE GF_NORMAL
  !
  !
  implicit none
  !Matrix:
  real(8),allocatable    :: Hmat(:,:),Hup(:,:),Hdw(:,:)
  integer                :: Nprint,N,Nloc
  real(8),allocatable    :: Eval(:),Eval_old(:)
  real(8),allocatable    :: Evec(:,:),Eps(:)
  !Lanczos
  integer                :: neigen,nitermax,niter,nblock,Nlanc
  integer                :: i,j,k
  integer                :: iup,idw,jup,jdw,ie
  real(8)                :: Vps,Deps
  complex(8),allocatable :: DeltaMats(:),DeltaReal(:)
  !
  !
  !< setup global dimension variables (in COMMON_VARS)
  Ns  = 7
  Nup = 4
  Ndw = 3
  !
  DimUp  = binomial(Ns,Nup)
  DimDw  = binomial(Ns,Ndw)
  Dim    = DimUp*DimDw
  !
  print*,"Dimensions =",DimUp,DimDw,Dim

  !< number of required eigenstates (1 == groundstate only)
  Neigen = 1
  !
  !< Build Hamiltonian, i.e. a symmetric sparse matrix.
  !The result is saved in spH0,spH0up,spH0dw data structure.
  call Build_spH(Nup,Ndw)
  !
  !
  !< Solve the Problem by exact diagonalization (Lapack) if Dim is small enough
  if(Dim < 2048)then
     allocate(Hmat(Dim,Dim))
     allocate(Hup(DimUp,DimUp))
     allocate(Hdw(DimDw,DimDw))
     call sp_dump_matrix(spH0,Hmat)
     call sp_dump_matrix(spH0up,Hup)
     call sp_dump_matrix(spH0dw,Hdw)
     !
     !< H = H_diag + H_dw * 1_up + 1_dw * H_up
     Hmat = Hmat + kronecker_product(Hdw,eye(DimUp))
     Hmat = Hmat + kronecker_product(eye(DimDw),Hup)
     !
     !< allocate Eigen-vectors (which hold Hmat on input)
     allocate(Evec(Dim,Dim))
     Evec = Hmat
     !< allocate Eigen-values (all)
     allocate(Eval(Dim))
     !< keep the old required Neigen eigen-values for error check later on
     allocate(Eval_old(5*Neigen))
     !
     !< perform the diagonalization
     write(*,"(A)")"LAPACK Diagonalization:"
     call eigh(Evec,Eval)
     !
     !< Write and save result:
     do i=1,5*Neigen
        write(*,*)Eval(i)
     end do
     !
     open(10,file="lapack_Eval_D.dat")
     do i=1,5*Neigen
        write(10,*)Eval(i)
     enddo
     close(10)
     !
     open(20,file="lapack_Evec1_D.dat")
     do i=1,Dim
        write(20,*)Evec(i,1)
     enddo
     close(20)
     !
     Eval_old=Eval(1:5*Neigen)
     !
     deallocate(Eval,Evec)
     print*,""
  endif
  !
  !
  !
  !< Solve the Problem by Lanczos diagonalization
  !
  !< Set the maximum iterations.
  Nitermax=min(Dim,512)
  !
  !< allocate Eigen-vector, we keep only 1
  allocate(Evec(Dim,1)); Evec=0d0
  !
  !< allocate eigen-values, we search only 1
  allocate(Eval(1)); Eval=0d0
  !
  !< perform the diagonalization
  write(*,"(A)")"LANCZOS DIAG :"
  call sp_lanczos(HxV_spH0,Dim,Nitermax,Eval(1),Evec(:,1))
  print*,""
  !
  !< Write and save result:
  do i=1,Neigen
     write(*,*)Eval(i)
  end do
  !
  print*,""
  if(Dim<2048)then
     do i=1,Neigen
        write(*,"(A,ES28.15)")"Error |E_lapack - E_lanczos|:",abs(Eval(i)-Eval_old(i))     
     end do
  endif
  !
  open(100,file="lanczos_Eval_D.dat")
  do i=1,Neigen
     write(100,*)Eval(i)
  end do
  close(100)
  !
  open(100,file="lanczos_Evec1_D.dat")
  do i=1,Dim
     write(100,*)Evec(i,1)
  end do
  close(100)
  !
  call Delete_spH()
  !
  !
  !
  call build_gf_lanczos(nup,ndw,Eval(1),Evec(:,1))
  !
  open(100,file="impGmats_iw.lanc")
  do i=1,size(wmats)
     write(100,*)wmats(i),dimag(impGmats(i)),dreal(impGmats(i))
  enddo
  close(100)
  !
  open(100,file="impGreal_wr.lanc")
  do i=1,size(wreal)
     write(100,*)wreal(i),dimag(impGreal(i)),dreal(impGreal(i))
  enddo
  close(100)
  !
  !
  deallocate(Eval,Evec)  
  print*,""
  !
  !
  !
  print*,"TESTING AGAINST NON-INTERACTING GF: G0(z) = (z - Delta(z))^-1; Delta(z) = sum Vi^2/(z-ei)"
  allocate(Eps(Ns-1))
  Vps  = 1d0/sqrt(dble(Ns-1))
  deps = 2d0/(Ns-1)
  do ie=1,Ns-1
     eps(ie) = -1d0 + (ie-1)*deps
  enddo
  !
  allocate(DeltaMats(size(wmats)))
  DeltaMats=0d0
  do i=1,size(wmats)
     do ie=1,Ns-1
        DeltaMats(i)=DeltaMats(i) + Vps**2/(dcmplx(0d0,wmats(i))-eps(ie))
     enddo
     DeltaMats(i) = dcmplx(0d0,wmats(i)) - DeltaMats(i)
  enddo
  DeltaMats = 1d0/DeltaMats
  !
  allocate(DeltaReal(size(wreal)))
  DeltaReal=0d0
  do i=1,size(wreal)
     do ie=1,Ns-1
        DeltaReal(i)=DeltaReal(i) + Vps**2/(dcmplx(wreal(i),eta)-eps(ie))
     enddo
     DeltaReal(i) = dcmplx(wreal(i),0.01d0) - DeltaReal(i)
  enddo
  DeltaReal = 1d0/DeltaReal

  open(100,file="impGmats_iw.test")
  do i=1,size(wmats)
     write(100,*)wmats(i),dimag(impGmats(i)),dreal(DeltaMats(i))
  enddo
  close(100)
  open(100,file="impGreal_wr.test")
  do i=1,size(wreal)
     write(100,*)wreal(i),dimag(impGreal(i)),dreal(DeltaReal(i))
  enddo
  close(100)
  !
  print*,"Test Matsubara:",sum(abs(impGmats-DeltaMats))/size(wmats)
  print*,"Test RealAxis:",sum(abs(impGreal-DeltaReal))/size(wreal)
end program lanczos_test






