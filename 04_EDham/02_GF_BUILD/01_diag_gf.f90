program diag_test
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
  real(8),allocatable    :: Hup(:,:),Hdw(:,:)
  integer                :: Nprint,N,Nloc,isector,unit,Nsectors
  real(8),allocatable    :: e0(:,:),Eps(:)
  integer                :: i,j,k,ie
  integer                :: iup,idw,jup,jdw
  real(8)                :: Egs
  real(8)                :: Vps,Deps
  complex(8),allocatable :: DeltaMats(:),DeltaReal(:)
  !
  !
  !< setup global dimension variables (in COMMON_VARS)
  !< Number of spins: keep this low for better visualization
  Ns    = 7
  Nsectors = (Ns+1)*(Ns+1)
  !< set default output unit. 6==std.output in fortran-esque
  unit=6
  write(unit,*)"Ns     =",Ns
  write(unit,*)"Ntot   =",2*Ns
  write(unit,*)"2**Ns  =",2**Ns
  write(unit,*)"2**Ntot=",2**(2*Ns)
  write(unit,*)""
  write(unit,*)""
  !
  !< Allocate espace data structure for all sectors:
  call init_espace(Ns)
  allocate(e0(0:Ns,0:Ns))
  !
  isector=0
  do Nup=0,Ns
     do Ndw=0,Ns
        isector=isector+1
        write(unit,"(A,I6,A3,I2,A1,I2,A)")"--- sector",isector,"  (",Nup,",",Ndw,")---- "
        !
        DimUp  = binomial(Ns,Nup)
        DimDw  = binomial(Ns,Ndw)
        Dim    = DimUp*DimDw
        !
        print*,"Dimensions =",DimUp,DimDw,Dim
        call setup_espace(Nup,Ndw,Dim)
        !
        !< Build Hamiltonian, i.e. a symmetric sparse matrix.
        !The result is saved in spH0,spH0up,spH0dw data structure.
        call Build_spH(Nup,Ndw)
        !
        !< Solve the Problem by exact diagonalization (Lapack)
        allocate(Hup(DimUp,DimUp))
        allocate(Hdw(DimDw,DimDw))
        call sp_dump_matrix(spH0,espace(Nup,Ndw)%M)
        call sp_dump_matrix(spH0up,Hup)
        call sp_dump_matrix(spH0dw,Hdw)
        !
        !< H = H_diag + H_dw * 1_up + 1_dw * H_up
        espace(Nup,Ndw)%M = espace(Nup,Ndw)%M + kronecker_product(Hdw,eye(DimUp))
        espace(Nup,Ndw)%M = espace(Nup,Ndw)%M + kronecker_product(eye(DimDw),Hup)
        !
        !< perform the diagonalization
        !< Hmat will contains Eigen-vectors after eigh
        write(*,"(A)",advance='no')"LAPACK Diagonalization:"
        call eigh(espace(Nup,Ndw)%M,espace(Nup,Ndw)%e)
        print*,espace(Nup,Ndw)%e(1)

        e0(Nup,Ndw)=minval(espace(Nup,Ndw)%e)
        !
        deallocate(Hup,Hdw)
        print*,""
        !
        call Delete_spH()
        !
     enddo
  enddo

  !< Get the GS energy and shift all the energies of the spectrum
  egs=minval(e0)
  forall(Nup=0:Ns,Ndw=0:Ns)espace(Nup,Ndw)%e = espace(Nup,Ndw)%e - egs
  !
  !< Get the partition function Z
  zeta_function=0d0
  do Nup=0,Ns
     do Ndw=0,Ns
        Dim = size(espace(Nup,Ndw)%e)
        do i=1,Dim
           zeta_function=zeta_function+exp(-beta*espace(Nup,Ndw)%e(i))
        enddo
     enddo
  enddo
  !
  write(unit,*)""
  write(unit,*)""
  call sleep(2)
  !
  call build_gf_exact()
  !
  open(100,file="impGmats_iw.diag")
  do i=1,size(wmats)
     write(100,*)wmats(i),dimag(impGmats(i)),dreal(impGmats(i))
  enddo
  close(100)
  !
  open(100,file="impGreal_wr.diag")
  do i=1,size(wreal)
     write(100,*)wreal(i),dimag(impGreal(i)),dreal(impGreal(i))
  enddo
  close(100)
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
     DeltaReal(i) = dcmplx(wreal(i),eta) - DeltaReal(i)
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
end program diag_test






