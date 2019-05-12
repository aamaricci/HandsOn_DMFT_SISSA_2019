program testLANCZOS_D
  USE COMMON_VARS
  USE SPARSE_MATRIX
  USE SETUP
  USE BUILD_HAM    
  implicit none
  !Matrix:
  real(8),allocatable            :: Hmat(:,:),Hup(:,:),Hdw(:,:)
  !
  !
  !< setup global dimension variables (in COMMON_VARS)
  Ns  = 6
  Nup = 3
  Ndw = 3
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
  !< Dump the sparse matrix (that are small enough here) to dense matrices
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
  !< Analyze and Print the sparse matrices.
  call sp_spy_matrix(spH0,"spH0")
  call sp_spy_matrix(spH0up,"spH0up")
  call sp_spy_matrix(spH0dw,"spH0dw")
  call sp_spy_matrix(Hmat,"Hmat")
  !
end program testLANCZOS_D






