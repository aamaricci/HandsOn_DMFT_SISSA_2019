program diag_test
  USE COMMON_VARS
  USE SPARSE_MATRIX
  USE BUILD_HAM
  USE SETUP
  !
  !
  implicit none
  !Matrix:
  real(8),allocatable            :: Hmat(:,:),Hup(:,:),Hdw(:,:)
  integer                        :: Neigen,N,Nloc
  real(8),allocatable            :: Eval(:)
  real(8),allocatable            :: Evec(:,:)
  integer                        :: i
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
  Neigen = 5
  !
  !< Build Hamiltonian, i.e. a symmetric sparse matrix.
  !The result is saved in spH0,spH0up,spH0dw data structure.
  call Build_spH(Nup,Ndw)
  !
  !
  !< Solve the Problem by exact diagonalization (Lapack) if Dim is small enough
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
  !
  !< allocate Eigen-values (all)
  allocate(Eval(Dim))
  !
  !< perform the diagonalization
  write(*,"(A)")"LAPACK Diagonalization:"
  call eigh(Evec,Eval)
  !
  !< Write and save result:
  do i=1,Neigen
     write(*,*)Eval(i)
  end do
  !
  open(10,file="lapack_Eval_D.dat")
  do i=1,Neigen
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
  deallocate(Eval,Evec)
  print*,""
end program diag_test






