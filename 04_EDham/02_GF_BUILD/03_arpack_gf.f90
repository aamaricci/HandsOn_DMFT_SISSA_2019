program arpack_test
  USE COMMON_VARS
  USE MATVEC_PRODUCT
  USE SPARSE_MATRIX
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
  Neigen = 5
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
     allocate(Eval_old(Neigen))
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
     Eval_old=Eval(1:Neigen)
     !
     deallocate(Eval,Evec)
     print*,""
  endif
  !
  !
  !< Solve the Problem by Arpack diagonalization
  !< Set the maximum iterations.
  Nitermax=min(Dim,512)
  !< Set the block size for Block Lanczos.
  Nblock=min(10*Neigen,Dim)
  !
  !< allocate Eigen-vector, we keep only 1
  allocate(Evec(Dim,Neigen)); Evec=0d0
  !
  !< allocate eigen-values, we search only 1
  allocate(Eval(Neigen));Eval=0d0
  !
  !< perform the diagonalization
  write(*,"(A)")"ARPACK DIAG :"
  call sp_arpack(HxV_spH0,Dim,Neigen,Nblock,Nitermax,Eval,Evec)
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
        write(*,"(I3,A,ES28.15)")i," Error |E_lapack - E_arpack|:",abs(Eval(i)-Eval_old(i))     
     end do
  endif
  !
  open(100,file="arpack_Eval_D.dat")
  do i=1,Neigen
     write(100,*)Eval(i)
  end do
  close(100)
  !
  open(100,file="arpack_Evec1_D.dat")
  do i=1,Dim
     write(100,*)Evec(i,1)
  end do
  close(100)
  !
  call Delete_spH()
  !
  !
  call build_gf_lanczos(nup,ndw,Eval(1),Evec(:,1))
  !
  open(100,file="impGmats_iw.arpack")
  do i=1,size(wmats)
     write(100,*)wmats(i),impGmats(i)
  enddo
  close(100)
  !
  open(100,file="impGreal_wr.arpack")
  do i=1,size(wreal)
     write(100,*)wreal(i),impGreal(i)
  enddo
  close(100)


  deallocate(Eval,Evec)
  print*,""


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
        DeltaReal(i)=DeltaReal(i) + Vps**2/(dcmplx(wreal(i),0.01d0)-eps(ie))
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

contains



  subroutine sp_arpack(MatVec,Ns,Neigen,Nblock,Nitermax,eval,evec,v0)
    !Interface to Matrix-Vector routine:
    interface
       subroutine MatVec(Nloc,vin,vout)
         integer                 :: Nloc
         real(8),dimension(Nloc) :: vin
         real(8),dimension(Nloc) :: vout
       end subroutine MatVec
    end interface
    !Arguments
    integer                   :: Ns
    integer                   :: Neigen
    integer                   :: Nblock
    integer                   :: Nitermax
    real(8)                   :: eval(Neigen)
    real(8)                   :: evec(Ns,Neigen)
    real(8),optional          :: v0(ns)
    !Dimensions:
    integer                   :: maxn,maxnev,maxncv,ldv
    integer                   :: n,nconv,ncv,nev
    !Arrays:
    real(8),allocatable       :: ax(:),d(:,:)
    real(8),allocatable       :: resid(:)
    real(8),allocatable       :: workl(:),workd(:)
    real(8),allocatable       :: v(:,:)
    logical,allocatable       :: select(:)
    integer                   :: iparam(11)
    integer                   :: ipntr(11)
    !Control Vars:
    integer                   :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
    logical                   :: rvec,verb
    integer                   :: i
    real(8)                   :: sigma
    real(8)                   :: tol_
    character                 :: bmat  
    character(len=2)          :: which_
    real(8),external          :: dnrm2
    !
    which_='SA'
    tol_=0d0
    verb=.true.
    !
    maxn   = Ns
    maxnev = Neigen
    maxncv = Nblock
    ldv    = Ns
    if(maxncv>Ns)maxncv=Ns
    n      = maxn
    nev    = maxnev
    ncv    = maxncv
    bmat   = 'I'
    maxitr = Nitermax
    ! 
    allocate(ax(n))
    allocate(d(ncv,2))
    allocate(resid(n))
    allocate(workl(ncv*(ncv+8)))
    allocate(workd(3*n))
    allocate(v(ldv,ncv))
    allocate(select(ncv))
    ax     =0d0
    d      =0d0
    resid  =0d0
    workl  =0d0
    workd  =0d0
    v      =0d0
    lworkl = ncv*(ncv+8)
    info   = 1
    ido    = 0
    ishfts    = 1
    mode1     = 1
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode1
    if(present(v0))then
       resid=v0
    else
       resid=1d0
    endif
    resid=resid/sqrt(dot_product(resid,resid))
    !
    !MAIN LOOP, REVERSE COMMUNICATION
    do
       call dsaupd(ido,bmat,n,which_,nev,tol_,resid,ncv,v,ldv,&
            iparam,ipntr,workd,workl,lworkl,info)
       if(ido/=-1.AND.ido/=1)then
          exit
       end if
       !  Perform matrix vector multiplication
       !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
       call MatVec(N,workd(ipntr(1)),workd(ipntr(2)) )
    end do
    !
    !POST PROCESSING:
    if(info<0)then
       write(*,'(a,i6)')'Error in DSAUPD, info = ', info
       stop
    else
       rvec = .true.
       call dseupd(rvec,'All',select,d,v,ldv,sigma,bmat,n,which_,&
            nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)
       do j=1,neigen
          eval(j)=d(j,1)
          do i=1,ns
             evec(i,j)=v(i,j)
          enddo
       enddo
       !
       !  Compute the residual norm
       !    ||  A*x - lambda*x ||
       !  for the NCONV accurately computed eigenvalues and 
       !  eigenvectors.  (iparam(5) indicates how many are 
       !  accurate to the requested tolerance)
       if(ierr/=0)then
          write(*,'(a,i6)')'Error with DSEUPD (get Evec), ierr = ',ierr
       end if
       !
       if(info==1) then
          write(*,'(a)' ) ' '
          write(*,'(a)' ) '  Maximum number of iterations reached.'
       elseif(info==3) then
          write(*,'(a)' ) ' '
          write(*,'(a)' ) '  No shifts could be applied during implicit '&
               //'Arnoldi update, try increasing NCV.'
       end if
    endif
  end subroutine sp_arpack

end program arpack_test






