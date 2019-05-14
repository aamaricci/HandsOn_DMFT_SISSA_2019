program lanczos_test
  USE COMMON_VARS
  USE MATVEC_PRODUCT
  USE SPARSE_MATRIX
  USE LANCZOS
  USE SETUP
  !
  !
  implicit none
  !Matrix:
  real(8),allocatable            :: Hmat(:,:),Hup(:,:),Hdw(:,:)
  integer                        :: Nprint,N,Nloc
  real(8),allocatable            :: Eval(:),Eval_old(:)
  real(8),allocatable            :: Evec(:,:)
  !Lanczos
  integer                        :: neigen,nitermax,niter,nblock,Nlanc
  integer                        :: i,j,k
  integer                        :: iup,idw,jup,jdw
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
     !< allocate Eigen-values (all)
     allocate(Eval(Dim))
     !< keep the old required Neigen eigen-values for error check later on
     allocate(Eval_old(5*Neigen))
     !
     !< perform the diagonalization
     write(*,"(A)")"LAPACK Diagonalization:"
     call eigh(Hmat,Eval)
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
        write(20,*)Hmat(i,1)
     enddo
     close(20)
     !
     Eval_old=Eval(1:5*Neigen)
     !
     deallocate(Eval,Hmat)
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
  do i=1,Nprint
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
  deallocate(Eval,Evec)  
  print*,""
  !
  !
  !
end program lanczos_test






