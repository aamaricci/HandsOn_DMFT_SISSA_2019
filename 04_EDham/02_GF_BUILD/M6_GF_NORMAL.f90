MODULE GF_NORMAL
  USE COMMON_VARS
  USE SPARSE_MATRIX
  USE LANCZOS
  USE SETUP
  USE MATVEC_PRODUCT
  implicit none

  private

  real(8),parameter           :: pi=acos(-1d0),beta=500d0,eps=0.01d0
  integer,parameter           :: Lmats=2048
  integer,parameter           :: Lreal=2000
  real(8),dimension(Lmats)    :: wmats
  real(8),dimension(Lreal)    :: wreal
  complex(8),dimension(Lmats) :: impGmats
  complex(8),dimension(Lreal) :: impGreal

  public :: build_gf_lanczos
  public :: impGmats,wmats
  public :: impGreal,wreal


contains


  subroutine build_gf_lanczos(Nup,Ndw,state_e,state_cvec)
    integer,intent(in)               :: Nup,Ndw
    real(8)                          :: state_e
    real(8),dimension(:)             :: state_cvec
    integer                          :: Mup,Mdw
    !
    integer                          :: idim,idimUP,idimDW
    integer                          :: jdim,jdimUP,jdimDW
    !
    integer                          :: i,iup,idw,ICup,ICdw
    integer                          :: j,jup,jdw,JCup,JCdw
    integer,dimension(Ns)            :: ibup,ibdw
    real(8),allocatable,dimension(:) :: vvinit
    real(8),allocatable,dimension(:) :: alfa_(:),beta_(:)
    !
    real(8)                          :: sgn,norm2,norm0
    integer                          :: Nitermax,Nlanc,vecDim
    !
    type(sector_map)                 :: HI(2),HJ(2)
    !
    wmats = pi/beta*(2*arange(1,Lmats)-1)
    wreal = linspace(-5d0,5d0,Lreal)    
    !
    !< Build the IN sector (Nup,Ndw). All variables are labelled with i-prefix
    iDimUp = get_sector_dimension(Ns,Nup)
    iDimDw = get_sector_dimension(Ns,Ndw)
    iDim   = iDimUp*iDimDw
    !
    call build_sector(Nup,Ndw,HI)
    !
    !
    !===============================================
    !ADD ONE PARTICLE spin UP:
    !< Build the OUT sector (Mup,Mdw): all variables are labelled with j-prefix
    Mup = Nup + 1
    Mdw = Ndw
    !
    jDimUp = get_sector_dimension(Ns,Mup)
    jDimDw = get_sector_dimension(Ns,Mdw)
    jDim   = jDimUp*jDimDw
    !
    write(*,"(A,2I6)")' add particle:',Mup,Mdw
    !
    !< Allocate the C^+|gs> = |vvinit> vector
    allocate(vvinit(jdim)) ; vvinit=0d0
    !
    call build_sector(Mup,Mdw,HJ)
    do i=1,iDim
       !< get the decomposition of i into iup + idw*2*DimUp
       iup = iup_index(i,iDimUp)
       idw = idw_index(i,iDimUp)
       !
       !< get the corresponding component in the Fock space I_up + I_dw*2**Ns
       Icup = HI(1)%map(iup)
       Icdw = HI(2)%map(idw)
       !
       !< get the partial binary decomposition
       ibup  = bdecomp(Icup,Ns)
       !
       !< check occupation at impurity
       if(ibup(1)/=0)cycle
       !
       !< apply C^+
       call cdg(1,Icup,Jcup,sgn)
       !
       !< retrieve the indices of the out state in the out sector J
       jup = binary_search(HJ(1)%map,Jcup)
       jdw = idw
       j   = jup + (jdw-1)*jDimUp
       !
       !< save the vector component
       vvinit(j) = sgn*state_cvec(i)
    enddo
    !< get rid of the J-sector (we do not need it anymore)
    call delete_sector(HJ)
    !
    !< normalize the |vvinit> vector
    norm2=dot_product(vvinit,vvinit)
    vvinit=vvinit/sqrt(norm2)
    !
    nlanc=min(jDim,300)
    allocate(alfa_(nlanc),beta_(nlanc))
    alfa_=0d0
    beta_=0d0
    !
    !< Build Hamiltonian Matrix of the J-sector
    call Build_spH(Mup,Mdw)
    !< perform tri-diagonalization of the H-matrix, stored in alfas,betas
    call sp_lanczos_tridiag(HxV_spH0,vvinit,alfa_,beta_)
    !< delete Hamiltonian sparse Matrix:   
    call Delete_spH()
    !< Store contribution to GF:
    call add_to_lanczos_gf_normal(norm2,state_e,alfa_,beta_,1)
    !
    deallocate(alfa_,beta_)
    if(allocated(vvinit))deallocate(vvinit)          
    !
    !
    !
    !REMOVE ONE PARTICLE:
    !< Build the OUT sector (Mup,Mdw): all variables are labelled with j-prefix
    Mup = Nup - 1
    Mdw = Ndw
    !
    jDimUp = get_sector_dimension(Ns,Mup)
    jDimDw = get_sector_dimension(Ns,Mdw)
    jDim   = jDimUp*jDimDw
    !
    write(*,"(A,2I6)")' del particle:',Mup,Mdw
    !
    !< Allocate the C|gs> = |vvinit> vector
    allocate(vvinit(jdim)) ; vvinit=0d0
    !
    call build_sector(Mup,Mdw,HJ)
    do i=1,iDim
       !< get the decomposition of i into iup + idw*2*DimUp
       iup = iup_index(i,iDimUp)
       idw = idw_index(i,iDimUp)
       !
       !< get the corresponding component in the Fock space I_up + I_dw*2**Ns
       Icup = HI(1)%map(iup)
       Icdw = HI(2)%map(idw)
       !
       !< get the partial binary decomposition
       ibup  = bdecomp(Icup,Ns)
       !
       !< check occupation at impurity
       if(ibup(1)/=1)cycle
       !
       !< apply C^+
       call c(1,Icup,Jcup,sgn)
       !
       !< retrieve the indices of the out state in the out sector J
       jup = binary_search(HJ(1)%map,Jcup)
       jdw = idw
       j   = jup + (jdw-1)*jDimUp
       !
       !< save the vector component
       vvinit(j) = sgn*state_cvec(i)
    enddo
    !< get rid of the J-sector (we do not need it anymore)
    call delete_sector(HJ)
    !
    !< normalize the |vvinit> vector
    norm2=dot_product(vvinit,vvinit)
    vvinit=vvinit/sqrt(norm2)
    !
    nlanc=min(jDim,300)
    allocate(alfa_(nlanc),beta_(nlanc))
    alfa_=0d0
    beta_=0d0
    !
    !< Build Hamiltonian Matrix of the J-sector
    call Build_spH(Mup,Mdw)
    !< perform tri-diagonalization of the H-matrix, stored in alfas,betas
    call sp_lanczos_tridiag(HxV_spH0,vvinit,alfa_,beta_)
    !< delete Hamiltonian sparse Matrix:   
    call Delete_spH()
    !< Store contribution to GF:
    call add_to_lanczos_gf_normal(norm2,state_e,alfa_,beta_,-1)
    !
    deallocate(alfa_,beta_)
    if(allocated(vvinit))deallocate(vvinit)          
    !===============================================
    !
    !
    call delete_sector(HI)
    !
    return
  end subroutine build_gf_lanczos



  subroutine add_to_lanczos_gf_normal(vnorm2,Ei,alanc,blanc,isign)
    real(8)                                    :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
    Nlanc = size(alanc)
    !
    pesoBZ = vnorm2
    !
    !< diagonalize the tridiagonal H-matrix: gets excitations energies and amplitudes Z(1,j)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       !
       !< K-L sum Matsubara
       do i=1,Lmats
          iw=dcmplx(0d0,wmats(i))
          impGmats(i)=impGmats(i) + peso/(iw-isign*de)
       enddo
       !< K-L sum RealAxis
       do i=1,Lreal
          iw=dcmplx(wreal(i),eps)
          impGreal(i)=impGreal(i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf_normal


END MODULE GF_NORMAL
