MODULE GF_NORMAL
  USE COMMON_VARS
  USE SPARSE_MATRIX
  USE LANCZOS
  USE SETUP
  USE MATVEC_PRODUCT
  implicit none

  private


  real(8),dimension(Lmats)    :: wmats
  real(8),dimension(Lreal)    :: wreal
  complex(8),dimension(Lmats) :: impGmats
  complex(8),dimension(Lreal) :: impGreal

  public :: build_gf_lanczos
  public :: build_gf_exact
  public :: impGmats,wmats
  public :: impGreal,wreal



  integer                          :: iDim,iDimUp,iDimDw,inup,indw
  integer                          :: jDim,jDimUp,jDimDw,jnup,jndw
  !
  integer                          :: i,iup,idw,ICup,ICdw
  integer                          :: j,jup,jdw,JCup,JCdw
  !

contains


  subroutine build_gf_exact()
    complex(8)            :: op_mat
    complex(8)            :: spectral_weight
    real(8)               :: sgn_cdg,sgn_c
    integer,dimension(Ns) :: ibup,ibdw
    integer               :: li,rj
    integer               :: m,i,j,r,k,p
    real(8)               :: Ei,Ej
    real(8)               :: expterm,peso,de,w0
    complex(8)            :: iw
    type(sector_map)      :: HI(2),HJ(2)
    !
    wmats = pi/beta*(2*arange(1,Lmats)-1)
    wreal = linspace(-5d0,5d0,Lreal)
    !
    UP: do iNup=0,Ns
       DW: do iNdw=0,Ns
          !
          jNup = iNup + 1       !arrival sector has one more electron with spin UP (non-magnetic solution)
          jNdw = iNdw
          if(jNup>Ns)cycle      !can not have more electrons UP than Ns
          !
          !< Get dimension of the starting and arriving sectors
          iDimUp  = get_sector_dimension(Ns,iNup)
          iDimDw  = get_sector_dimension(Ns,iNdw)
          iDim    = iDimUp*iDimDw
          !
          jDimUp  = get_sector_dimension(Ns,jNup)
          jDimDw  = get_sector_dimension(Ns,jNdw)
          jDim    = jDimUp*jDimDw
          !
          call build_sector(iNup,iNdw,HI)
          call build_sector(jNup,jNdw,HJ)
          !
          do i=1,iDim          !loop over the states in the i-th sect.
             do j=1,jDim       !loop over the states in the j-th sect.
                op_mat=0.d0
                expterm=exp(-beta*espace(iNup,iNdw)%e(i))+exp(-beta*espace(jNup,jNdw)%e(j))
                if(expterm < 1d-9)cycle
                !
                do li=1,idim              !loop over the component of |I> (IN state!)
                   !< get the decomposition of i into iup + idw*2*DimUp
                   iup = iup_index(li,iDimUp)
                   idw = idw_index(li,iDimUp)
                   !
                   !< get the corresponding component in the Fock space I_up + I_dw*2**Ns
                   Icup = HI(1)%map(iup)
                   Icdw = HI(2)%map(idw)
                   !
                   !< get the partial binary decomposition
                   ibup  = bdecomp(Icup,Ns)
                   !
                   !< check occupation at impurity
                   if(ibup(1)==1)cycle
                   !< apply C^+
                   call cdg(1,Icup,Jcup,sgn_cdg)
                   !
                   !< retrieve the indices of the out state in the out sector J
                   jup = binary_search(HJ(1)%map,Jcup)
                   jdw = idw
                   rj   = jup + (jdw-1)*jDimUp
                   !
                   op_mat=op_mat + espace(jNup,jNdw)%M(rj,j)*sgn_cdg*espace(iNup,iNdw)%M(li,i)
                enddo
                !
                Ei=espace(iNup,iNdw)%e(i)
                Ej=espace(jNup,jNdw)%e(j)
                de=Ej-Ei
                peso=expterm/zeta_function
                spectral_weight=peso*abs(op_mat)**2
                !
                do m=1,Lmats
                   iw=dcmplx(0d0,wmats(m))
                   impGmats(m)=impGmats(m)+spectral_weight/(iw-de)
                enddo
                !
                do m=1,Lreal
                   iw=dcmplx(wreal(m),eta)
                   impGreal(m)=impGreal(m)+spectral_weight/(iw-de)
                enddo
                !
             enddo
          enddo
          !
          call Delete_sector(HI)
          call Delete_sector(HJ)
          !          
       enddo DW
    enddo UP
  end subroutine build_gf_exact






  subroutine build_gf_lanczos(Nup,Ndw,state_e,state_cvec)
    integer,intent(in)               :: Nup,Ndw
    real(8)                          :: state_e
    real(8),dimension(:)             :: state_cvec
    integer                          :: Mup,Mdw

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
          iw=dcmplx(wreal(i),eta)
          impGreal(i)=impGreal(i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf_normal


END MODULE GF_NORMAL
