program diag_test
  USE COMMON_VARS
  USE MATVEC_PRODUCT
  USE SPARSE_MATRIX
  USE LANCZOS
  USE SETUP
  !
  !
  implicit none
  !Matrix:
  real(8),allocatable :: Hmat(:,:),Hup(:,:),Hdw(:,:)
  integer             :: Nprint,N,Nloc,isector,unit
  real(8),allocatable :: Eval(:),e0(:,:)
  real(8),allocatable :: Evec(:,:)
  !Lanczos
  real(8)             :: Egs
  integer             :: i,j,k
  integer             :: iup,idw,jup,jdw
  !
  !
  !< setup global dimension variables (in COMMON_VARS)
  !< Number of spins: keep this low for better visualization
  Ns  = 7
  Nup = 4
  Ndw = 3
  !< set default output unit. 6==std.output in fortran-esque
  unit=6
  write(unit,*)"Ns     =",Ns
  write(unit,*)"Ntot   =",2*Ns
  write(unit,*)"2**Ns  =",2**Ns
  write(unit,*)"2**Ntot=",2**(2*Ns)
  write(unit,*)""
  write(unit,*)""
  !
  DimUp  = binomial(Ns,Nup)
  DimDw  = binomial(Ns,Ndw)
  Dim    = DimUp*DimDw
  !
  print*,"Dimensions =",DimUp,DimDw,Dim
  !
  !< Build Hamiltonian, i.e. a symmetric sparse matrix.
  !The result is saved in spH0,spH0up,spH0dw data structure.
  call Build_spH(Nup,Ndw)
  !
  !< Solve the Problem by exact diagonalization (Lapack)
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
  !< allocate Eigen-values (all)
  allocate(Eval(Dim));Eval=0d0
  !
  !< perform the diagonalization
  !< Hmat will contains Eigen-vectors after eigh
  write(*,"(A)")"LAPACK Diagonalization:"
  call eigh(Hmat,Eval)
  !
  !< Write and save result:
  do i=1,5
     write(*,*)Eval(i)
  end do
  !
  deallocate(Hmat,Eval,Hup,Hdw)
  print*,""
  !
  call Delete_spH()
  !

  !
end program diag_test






