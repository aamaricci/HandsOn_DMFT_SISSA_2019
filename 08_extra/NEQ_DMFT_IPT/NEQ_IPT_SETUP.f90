module NEQ_IPT_SETUP
  USE NEQ_INPUT_VARS  
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE SF_CONSTANTS
  USE SF_IOTOOLS
  USE SF_SPECIAL
  USE SF_LINALG, only: zeros
  USE SF_MISC, only: assert_shape
  USE DMFT_GK
  USE DMFT_CTRL_VARS
  USE DMFT_FFTGF
  ! USE DMFT_FFTAUX
  ! USE DMFT_MISC
  private


  interface neq_ipt_setup_kb_contour_gf
     module procedure :: neq_continue_equilibirum_lattice_main
     module procedure :: neq_continue_equilibirum_lattice_Nso
     !
     module procedure :: neq_continue_equilibirum_bethe_main
     module procedure :: neq_continue_equilibirum_bethe_Nso
     !
     module procedure :: neq_continue_equilibirum_lattice_nonint_main
     module procedure :: neq_continue_equilibirum_lattice_nonint_Nso
  end interface neq_ipt_setup_kb_contour_gf


  interface neq_setup_contour_gf
     module procedure :: neq_setup_contour_gf_main
     module procedure :: neq_setup_contour_gf_Nso
  end interface neq_setup_contour_gf


  interface neq_setup_contour_sigma
     module procedure :: neq_setup_contour_sigma_main
     module procedure :: neq_setup_contour_sigma_Nso
  end interface neq_setup_contour_sigma


  interface neq_setup_contour_dgf
     module procedure :: neq_setup_contour_dgf_main
     module procedure :: neq_setup_contour_dgf_Nso
  end interface neq_setup_contour_dgf


  public :: neq_ipt_setup_kb_contour_gf     ! read G0 and continue the equilibrium GF,Sigma,G0 to the t=0 contour.


  public :: neq_setup_contour_gf
  public :: neq_setup_contour_sigma
  public :: neq_setup_contour_dgf





contains




  subroutine neq_continue_equilibirum_bethe_main(params,g0,dg0,g,self,wband)
    type(kb_contour_params)             :: params
    type(kb_contour_gf)                 :: g0
    type(kb_contour_dgf)                :: dg0
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    real(8),optional                    :: wband
    real(8)                             :: wband_
    real(8)                             :: wm,res,ims
    logical                             :: bool
    integer                             :: i,j,k,ik,unit,len,N,L,Niw
    complex(8)                          :: zeta
    complex(8)                          :: Self_gtr
    complex(8),allocatable,dimension(:) :: GxG0
    real(8)                             :: Scoeff(2),Gcoeff(4)
    !
    wband_=1d0;if(present(wband))wband_=wband
    !
    call check_kb_contour_gf(params,G0," neq_continue_equilibirum_bethe_main")
    call check_kb_contour_gf(params,G," neq_continue_equilibirum_bethe_main")
    call check_kb_contour_gf(params,Self," neq_continue_equilibirum_bethe_main")
    call check_kb_contour_gf(params,dG0," neq_continue_equilibirum_bethe_main")
    !
    N = params%Nt
    L = params%Ntau
    Niw= params%Niw
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    call neq_setup_contour_gf(params,g0,reg(g0file)) !read G0 from file:
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call neq_setup_contour_sigma(params,g0,self)
    !
    !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
    do i=1,Niw
       zeta  = xi*params%wm(i)+xmu - Self%iw(i)
       G%iw(i) = gfbethe(wm,zeta,wband_)
    enddo
    ! Gcoeff      = tail_coeff_glat(U,0.5d0,0d0,0d0)
    call fft_gf_iw2tau(G%iw,G%mats,beta)!,Gcoeff)     !get G(tau)
    G%ret(1,1)  = -xi                !get G^R(0,0)=-xi
    G%less(1,1) = -xi*G%mats(L)      !get G^<(0,0)= xi*G^M(0-)
    G%lmix(1,0:L)=-xi*G%mats(L:0:-1) !get G^\lmix(0,tau)=-xi*G(beta-tau>0)
    !
    !Derivatives
    allocate(GxG0(0:L))
    !
    !get d/dt G0^R = 0.d0
    dG0%ret(1)  = zero
    !
    !get d/dt G0^< = -xi(-xi)int_0^beta G^\lmix * G0^\rmix
    do k=0,L
       GxG0(k)=G%lmix(1,k)*conjg(G0%lmix(1,L-k))
    end do
    dG0%less(1) = -xi*(-xi)*params%dtau*kb_trapz(GxG0(0:),0,L)
    !
    !get d/dt G0^\lmix = -xi*int_0^beta G0^\lmix * G0^M
    dG0%lmix(0:)= zero
    do j=0,L
       do k=0,j
          GxG0(k)=G%lmix(1,k)*G0%mats(k+L-j)
       end do
       dG0%lmix(j)=dG0%lmix(j)+xi*params%dtau*kb_trapz(GxG0(0:),0,j)
       do k=j,L
          GxG0(k)=G%lmix(1,k)*G0%mats(k-j)
       end do
       dG0%lmix(j)=dG0%lmix(j)-xi*params%dtau*kb_trapz(GxG0(0:),j,L) 
    enddo
    deallocate(GxG0)
    return
  end subroutine neq_continue_equilibirum_bethe_main

  subroutine neq_continue_equilibirum_bethe_Nso(params,g0,dg0,g,self,wband)
    type(kb_contour_params)             :: params
    type(kb_contour_gf)                 :: g0(:,:,:,:) ![Nspin,Nspin,Norb,Norb]
    type(kb_contour_dgf)                :: dg0(:,:,:,:) !as g0
    type(kb_contour_gf)                 :: g(:,:,:,:)   !as g0
    type(kb_contour_gf)                 :: self(:,:,:,:) !as g0
    real(8),optional                    :: wband
    real(8)                             :: wband_
    real(8)                             :: wm,res,ims
    logical                             :: bool
    integer                             :: Nspin,Norb,ispin,jspin,iorb,jorb
    integer                             :: i,j,k,ik,unit,len,N,L,Niw
    complex(8)                          :: zeta
    complex(8)                          :: Self_gtr
    complex(8),allocatable,dimension(:) :: GxG0
    real(8)                             :: Scoeff(2),Gcoeff(4)
    !
    wband_=1d0;if(present(wband))wband_=wband
    !
    Nspin = size(g0,1)
    Norb  = size(g0,3)
    !
    call check_kb_contour_gf(params,G0," neq_continue_equilibirum_bethe_main")
    call check_kb_contour_gf(params,G," neq_continue_equilibirum_bethe_main")
    call check_kb_contour_gf(params,Self," neq_continue_equilibirum_bethe_main")
    call check_kb_contour_gf(params,dG0," neq_continue_equilibirum_bethe_main")
    !
    N = params%Nt
    L = params%Ntau
    Niw= params%Niw
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    call neq_setup_contour_gf(params,g0,reg(g0file)) !read G0 from file:
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call neq_setup_contour_sigma(params,g0,self)
    !
    !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
    allocate(GxG0(0:L))
    !
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Niw
             zeta  = xi*params%wm(i)+xmu - Self(ispin,ispin,iorb,iorb)%iw(i)
             G(ispin,ispin,iorb,iorb)%iw(i) = gfbethe(wm,zeta,wband_)
          enddo
          !
          ! Gcoeff      = tail_coeff_glat(U,0.5d0,0d0,0d0)
          call fft_gf_iw2tau(G(ispin,ispin,iorb,iorb)%iw,G(ispin,ispin,iorb,iorb)%mats(0:),beta)!,Gcoeff)     !get G(tau)
          !
          G(ispin,ispin,iorb,iorb)%ret(1,1)   = -xi                                    !get G^R(0,0)=-xi
          G(ispin,ispin,iorb,iorb)%less(1,1)  = -xi*G(ispin,ispin,iorb,iorb)%mats(L)   !get G^<(0,0)= xi*G^M(0-)
          G(ispin,ispin,iorb,iorb)%lmix(1,0:L)= -xi*G(ispin,ispin,iorb,iorb)%mats(L:0:-1)!L-i) !get G^\lmix(0,tau)=-xi*G(beta-tau>0)
          ! forall(i=0:L)G%lmix(1,i)=-xi*G%mats(L-i) !get G^\lmix(0,tau)=-xi*G(beta-tau>0)
          !
          !Derivatives
          !get d/dt G0^R = 0.d0
          dG0(ispin,ispin,iorb,iorb)%ret(1)  = zero
          !
          !get d/dt G0^< = -xi(-xi)int_0^beta G^\lmix * G0^\rmix
          do k=0,L
             GxG0(k)=G(ispin,ispin,iorb,iorb)%lmix(1,k)*conjg(G0(ispin,ispin,iorb,iorb)%lmix(1,L-k))
          end do
          dG0(ispin,ispin,iorb,iorb)%less(1) = -xi*(-xi)*params%dtau*kb_trapz(GxG0(0:),0,L)
          !
          !get d/dt G0^\lmix = -xi*int_0^beta G0^\lmix * G0^M
          dG0(ispin,ispin,iorb,iorb)%lmix(0:)= zero
          do j=0,L
             do k=0,j
                GxG0(k)=G(ispin,ispin,iorb,iorb)%lmix(1,k)*G0(ispin,ispin,iorb,iorb)%mats(k+L-j)
             end do
             dG0(ispin,ispin,iorb,iorb)%lmix(j)=dG0(ispin,ispin,iorb,iorb)%lmix(j)+xi*params%dtau*kb_trapz(GxG0(0:),0,j)
             do k=j,L
                GxG0(k)=G(ispin,ispin,iorb,iorb)%lmix(1,k)*G0(ispin,ispin,iorb,iorb)%mats(k-j)
             end do
             dG0(ispin,ispin,iorb,iorb)%lmix(j)=dG0(ispin,ispin,iorb,iorb)%lmix(j)-xi*params%dtau*kb_trapz(GxG0(0:),j,L) 
          enddo
       enddo
    enddo
    !
    deallocate(GxG0)
    return
  end subroutine neq_continue_equilibirum_bethe_Nso








  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the  equilibrium to Keldysh contour
  !+-------------------------------------------------------------------+
  subroutine neq_continue_equilibirum_lattice_main(params,g0,gk,dgk,g,self,Hk,Wtk)
    type(kb_contour_params)             :: params
    type(kb_contour_gf)                 :: g0
    type(kb_contour_gf)                 :: gk(:)
    type(kb_contour_dgf)                :: dgk(size(gk))
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    complex(8)                          :: Hk(size(gk))
    real(8)                             :: Wtk(size(gk))
    complex(8),allocatable              :: Gkmats(:)
    integer                             :: i,j,ik,N,L,Niw,Lk
    !
    Lk=size(gk)
    !
    N = params%Nt
    L = params%Ntau
    Niw= params%Niw
    !
    call check_kb_contour_gf(params,G0," neq_continue_equilibirum_lattice")
    call check_kb_contour_gf(params,G," neq_continue_equilibirum_lattice")
    call check_kb_contour_gf(params,Self," neq_continue_equilibirum_lattice")
    do ik=1,Lk
       call check_kb_contour_gf(params,Gk(ik)," neq_continue_equilibirum_lattice")
       call check_kb_contour_gf(params,dGk(ik)," neq_continue_equilibirum_lattice")
    enddo
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    call neq_setup_contour_gf(params,g0,trim(g0file)) !read G0 from file:
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call neq_setup_contour_sigma(params,g0,self)
    !
    !INITIALIZE THE GREENS FUNCTIONS G_loc, G_k & dG_k ^{x=M,<,R,\lmix}
    allocate(Gkmats(Niw))
    do ik=1,Lk
       Gkmats = one/(xi*params%wm(:) + xmu - Hk(ik) - Self%iw(:))
       call neq_setup_contour_gf(params,Gk(ik),gmats=Gkmats)
       !
       call neq_setup_contour_dgf(params,Hk(ik),Gk(ik),Self,dGk(ik))
       !
    enddo
    call sum_kb_contour_gf(params,Wtk(:),Gk(:),G)
    !
    deallocate(Gkmats)
  end subroutine neq_continue_equilibirum_lattice_main

  subroutine neq_continue_equilibirum_lattice_Nso(params,g0,gk,dgk,g,self,Hk,Wtk)
    type(kb_contour_params) :: params
    type(kb_contour_gf)     :: g0(:,:,:,:)   ![Nspin,Nspin,Norb,Norb]
    type(kb_contour_gf)     :: gk(:,:,:,:,:) ![Lk,Nspin,Nspin,Norb,Norb]
    type(kb_contour_dgf)    :: dgk(:,:,:,:,:) !as gk
    type(kb_contour_gf)     :: g(:,:,:,:)     !as g0
    type(kb_contour_gf)     :: self(:,:,:,:)  !as g0
    complex(8)              :: Hk(:,:,:)      ![Nso,Nso,Nk]
    real(8)                 :: Wtk(:)         ![Nk]
    complex(8),allocatable  :: Gkmats(:,:,:,:,:),Smats(:,:,:,:,:) ![Nspin,Nspin,Norb,Norb,Niw]
    integer                 :: Nspin,Norb,Lk
    integer                 :: ispin,jspin,iorb,jorb
    integer                 :: i,j,ik,N,L,Niw
    !
    Nspin = size(g0,1)
    Norb  = size(g0,3)
    !
    Lk    = size(Hk,3)
    !
    N = params%Nt
    L = params%Ntau
    Niw= params%Niw
    !
    call assert_shape_kb_contour_gf(g0,[Nspin,Nspin,Norb,Norb],"neq_continue_equilibirum_lattice_Nso","G0")
    call assert_shape_kb_contour_gf(g,[Nspin,Nspin,Norb,Norb],"neq_continue_equilibirum_lattice_Nso","G")
    call assert_shape_kb_contour_gf(self,[Nspin,Nspin,Norb,Norb],"neq_continue_equilibirum_lattice_Nso","Self")
    call assert_shape_kb_contour_gf(gk,[Lk,Nspin,Nspin,Norb,Norb],"neq_continue_equilibirum_lattice_Nso","Gk")
    call assert_shape_kb_contour_gf(dgk,[Lk,Nspin,Nspin,Norb,Norb],"neq_continue_equilibirum_lattice_Nso","dGk")
    call assert_shape(Hk,[Nspin*Norb,Nspin*Norb,Lk],"neq_continue_equilibirum_lattice_Nso","Hk")
    call assert_shape(Wtk,[Lk],"neq_continue_equilibirum_lattice_Nso","Wtk")
    !
    call check_kb_contour_gf(params,G0,"neq_continue_equilibirum_lattice_Nso","G0")
    call check_kb_contour_gf(params,G,"neq_continue_equilibirum_lattice_Nso","G")
    call check_kb_contour_gf(params,Self,"neq_continue_equilibirum_lattice_Nso","Self")
    call check_kb_contour_gf(params,Gk,"neq_continue_equilibirum_lattice_Nso","Gk")
    call check_kb_contour_gf(params,dGk,"neq_continue_equilibirum_lattice_Nso","dGk")
    ! enddo
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    call neq_setup_contour_gf(params,g0,reg(g0file)) !read G0 from file:
    !
    !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
    call neq_setup_contour_sigma(params,g0,self)
    !
    allocate(Smats(Nspin,Nspin,Norb,Norb,Niw))
    forall(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)Smats(ispin,jspin,iorb,jorb,:) = Self(ispin,jspin,iorb,jorb)%iw(:)
    !
    !INITIALIZE THE GREENS FUNCTIONS G_loc, G_k & dG_k ^{x=M,<,R,\lmix}
    allocate(Gkmats(Nspin,Nspin,Norb,Norb,Niw))
    do ik=1,Lk
       call add_ctrl_var(beta,"beta")
       call add_ctrl_var(xmu,"xmu")
       call dmft_gk_matsubara(Hk(:,:,ik),Wtk(ik),Gkmats,Smats)
       !
       call neq_setup_contour_gf(params,Gk(ik,:,:,:,:),gmats=Gkmats)
       !
       call neq_setup_contour_dgf(params,Hk(:,:,ik),Gk(ik,:,:,:,:),Self,dGk(ik,:,:,:,:))
       !
    enddo
    call sum_kb_contour_gf(params,Wtk,Gk,G)
    !
    deallocate(Gkmats)
  end subroutine neq_continue_equilibirum_lattice_Nso








  !+-------------------------------------------------------------------+
  !PURPOSE: obtain and continue the  equilibrium to Keldysh contour
  ! NON INTERACTING CASE (FOR PRELIMINARY STUDIES)
  !+-------------------------------------------------------------------+
  subroutine neq_continue_equilibirum_lattice_nonint_main(params,gk,dgk,g,Hk,Wtk)
    type(kb_contour_params)             :: params
    type(kb_contour_gf)                 :: gk(:)
    type(kb_contour_dgf)                :: dgk(size(gk))
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: Szero
    complex(8)                          :: Hk(size(gk))
    real(8)                             :: Wtk(size(gk))
    complex(8),allocatable              :: Gkmats(:)
    integer                             :: i,j,ik,N,L,Niw,Lk
    !
    Lk=size(gk)
    !
    N = params%Nt
    L = params%Ntau
    Niw= params%Niw
    !    
    call check_kb_contour_gf(params,G," neq_continue_equilibirum_nonint")
    do ik=1,Lk
       call check_kb_contour_gf(params,Gk(ik)," neq_continue_equilibirum_nonint")
       call check_kb_contour_gf(params,dGk(ik)," neq_continue_equilibirum_nonint")
    enddo
    !
    !INITIALIZE THE GREENS FUNCTIONS G_loc, G_k & dG_k ^{x=M,<,R,\lmix}
    allocate(Gkmats(Niw))
    call allocate_kb_contour_gf(params,Szero) !zero Self-Energy function
    do ik=1,Lk
       Gkmats = one/(xi*params%wm(:) + xmu - Hk(ik))
       call neq_setup_contour_gf(params,Gk(ik),gmats=Gkmats)
       !
       call neq_setup_contour_dgf(params,Hk(ik),Gk(ik),Szero,dGk(ik))
       !
    enddo
    call sum_kb_contour_gf(params,Wtk(:),Gk(:),G)
    !
    call deallocate_kb_contour_gf(Szero)
    deallocate(Gkmats)
  end subroutine neq_continue_equilibirum_lattice_nonint_main

  subroutine neq_continue_equilibirum_lattice_nonint_Nso(params,gk,dgk,g,Hk,Wtk)
    type(kb_contour_params)         :: params
    type(kb_contour_gf)             :: gk(:,:,:,:,:) ![Lk,Nspin,Nspin,Norb,Norb]
    type(kb_contour_dgf)            :: dgk(:,:,:,:,:) !as gk
    type(kb_contour_gf)             :: g(:,:,:,:)     !as g0
    type(kb_contour_gf),allocatable :: Szero(:,:,:,:)  !as g0
    complex(8)                      :: Hk(:,:,:)      ![Nso,Nso,Nk]
    real(8)                         :: Wtk(:)         ![Nk]
    complex(8),allocatable          :: Gkmats(:,:,:,:,:) ![Nspin,Nspin,Norb,Norb,Niw]
    integer                         :: Nspin,Norb,Lk
    integer                         :: ispin,jspin,iorb,jorb
    integer                         :: i,j,ik,N,L,Niw
    !
    Nspin = size(g,1)
    Norb  = size(g,3)
    !
    Lk    = size(Hk,3)
    !
    N = params%Nt
    L = params%Ntau
    Niw= params%Niw
    !
    call assert_shape_kb_contour_gf(g,[Nspin,Nspin,Norb,Norb],"neq_continue_equilibirum_lattice_nonint_Nso","G")
    call assert_shape_kb_contour_gf(gk,[Lk,Nspin,Nspin,Norb,Norb],"neq_continue_equilibirum_lattice_nonint_Nso","Gk")
    call assert_shape_kb_contour_gf(dgk,[Lk,Nspin,Nspin,Norb,Norb],"neq_continue_equilibirum_lattice_nonint_Nso","dGk")
    call assert_shape(Hk,[Nspin*Norb,Nspin*Norb,Lk],"neq_continue_equilibirum_lattice_nonint_Nso","Hk")
    call assert_shape(Wtk,[Lk],"neq_continue_equilibirum_lattice_nonint_Nso","Wtk")
    !
    call check_kb_contour_gf(params,G,"neq_continue_equilibirum_lattice_nonint_Nso","G")
    call check_kb_contour_gf(params,Gk,"neq_continue_equilibirum_lattice_nonint_Nso","Gk")
    call check_kb_contour_gf(params,dGk,"neq_continue_equilibirum_lattice_nonint_Nso","dGk")
    !
    !
    !INITIALIZE THE GREENS FUNCTIONS G_loc, G_k & dG_k ^{x=M,<,R,\lmix}
    allocate(Gkmats(Nspin,Nspin,Norb,Norb,Niw))
    allocate(Szero(Nspin,Nspin,Norb,Norb))
    call allocate_kb_contour_gf(params,Szero)
    do ik=1,Lk
       call add_ctrl_var(beta,"beta")
       call add_ctrl_var(xmu,"xmu")
       call dmft_gk_matsubara(Hk(:,:,ik),Wtk(ik),Gkmats,zeros(Nspin,Nspin,Norb,Norb,Niw))
       !
       call neq_setup_contour_gf(params,Gk(ik,:,:,:,:),gmats=Gkmats)
       !
       call neq_setup_contour_dgf(params,Hk(:,:,ik),Gk(ik,:,:,:,:),Szero,dGk(ik,:,:,:,:))
       !
    enddo
    call sum_kb_contour_gf(params,Wtk,Gk,G)
    !
    call deallocate_kb_contour_gf(Szero)
    deallocate(Gkmats)
  end subroutine neq_continue_equilibirum_lattice_nonint_Nso

























  !+-------------------------------------------------------------------+
  !PURPOSE: setup initial conditions for k-resolved GF
  !this procedure evaluate the derivative of a GF 
  !and set it to the S-K contour.
  !It should be clear that H MUST contain the Hartree term if any in
  !the kernel (i.e. is K--> Self-energy)
  subroutine neq_setup_contour_dgf_main(params,Hk,Gk,Self,dGk)
    type(kb_contour_params)             :: params
    complex(8)                          :: Hk
    type(kb_contour_gf)                 :: Gk
    type(kb_contour_gf)                 :: Self
    type(kb_contour_dgf)                :: dGk
    !
    integer                             :: i,j,k,L,Niw
    real(8)                             :: nk
    complex(8),allocatable,dimension(:) :: SxG
    real(8),dimension(:),allocatable    :: ftau
    !
    L = params%Ntau
    Niw= params%Niw
    !
    call check_kb_contour_gf(params,Gk,"neq_setup_contour_dgf_main")
    call check_kb_contour_gf(params,Self,"neq_setup_contour_dgf_main")
    call check_kb_contour_gf(params,dGk,"neq_setup_contour_dgf_main")
    !
    !Derivatives
    allocate(SxG(0:L))
    !
    !get d/dt G_k^R = -i H(k,0)G_k^R
    dGk%ret(1)  = -xi*Hk*Gk%ret(1,1)
    !
    !get d/dt G_k^< = -i H(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
    dGk%less(1) = -xi*Hk*Gk%less(1,1)
    do k=0,L
       SxG(k)=self%lmix(1,k)*conjg(Gk%lmix(1,L-k))
    end do
    dGk%less(1) = dGk%less(1) - xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,L)
    !
    !get d/dt G_k^\lmix = -xi*H(k,0)*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
    dGk%lmix(0:)= -xi*Hk*Gk%lmix(1,0:)
    do j=0,L
       do k=0,j
          SxG(k)=self%lmix(1,k)*Gk%mats(k+L-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
       do k=j,L
          SxG(k)=self%lmix(1,k)*Gk%mats(k-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,L) 
    enddo
  end subroutine neq_setup_contour_dgf_main


  subroutine neq_setup_contour_dgf_Nso(params,Hk,Gk,Self,dGk)
    type(kb_contour_params)                 :: params
    complex(8),dimension(:,:)               :: Hk    ![Nso,Nso] is fixed
    type(kb_contour_gf),dimension(:,:,:,:)  :: Gk    ![Nspin,Nspin,Norb,Norb]
    type(kb_contour_gf),dimension(:,:,:,:)  :: Self  !as G
    type(kb_contour_dgf),dimension(:,:,:,:) :: dGk   !as G
    !
    integer                                 :: N,L,Niw
    integer                                 :: Nspin,Norb,Nso
    integer                                 :: ispin,jspin,kspin
    integer                                 :: iorb,jorb,korb    
    integer                                 :: ii,jj,kk,j
    integer                                 :: s
    complex(8)                              :: HxG
    complex(8),allocatable,dimension(:)     :: SxG
    !
    N = params%Nt ;     if( N > 1) stop "ERROR neq_setup_contour_dG_Nso: N > 1"
    L = params%Ntau
    Niw= params%Niw
    !
    Nspin = size(Gk,1)
    Norb  = size(Gk,3)
    Nso   = Nspin*Norb
    !
    call check_kb_contour_gf(params,Gk,"neq_setup_contour_dG_Nso")
    call check_kb_contour_gf(params,Self,"neq_setup_contour_dG_Nso")
    call check_kb_contour_gf(params,dGk,"neq_setup_contour_dG_Nso")
    if(any( shape(Hk) /= [Nso,Nso] ) ) stop "ERROR neq_setup_contour_dG_Nso: shape[Hk] != [Nso][Nso]"
    !
    !Derivatives
    allocate(SxG(0:L))
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                ii = iso_indx(ispin,iorb)
                jj = iso_indx(jspin,jorb)
                !
                !Derivatives
                !
                !get d/dt G_k^R = -i H(k,0)G_k^R
                HxG=zero
                do kspin=1,Nspin
                   do korb=1,Norb
                      kk = iso_indx(kspin,korb)
                      HxG = HxG - xi*Hk(ii,kk)*Gk(kspin,jspin,korb,jorb)%ret(1,1)
                   enddo
                enddo
                dGk(ispin,jspin,iorb,jorb)%ret(1)  = HxG
                ! dGk%ret(1)  = -xi*Hk*Gk%ret(1,1)
                !
                !
                !get d/dt G_k^< = -i H(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
                HxG=zero
                do kspin=1,Nspin
                   do korb=1,Norb
                      kk = iso_indx(kspin,korb)
                      HxG = HxG - xi*Hk(ii,kk)*Gk(kspin,jspin,korb,jorb)%less(1,1)
                   enddo
                enddo
                ! dGk%less(1) = -xi*Hk*Gk%less(1,1)
                !
                SxG(0:)=zero
                do s=0,L
                   do kspin=1,Nspin
                      do korb=1,Norb
                         SxG(s)=SxG(s)+Self(ispin,kspin,iorb,korb)%lmix(1,s)*conjg( Gk(kspin,jspin,korb,jorb)%lmix(1,L-s) )
                      enddo
                   enddo
                enddo
                dGk(ispin,jspin,iorb,jorb)%less(1) = HxG - xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,L) 
                ! do k=0,L
                !    SxG(k)=self%lmix(1,k)*conjg(Gk%lmix(1,L-k))
                ! end do
                ! dGk%less(1) = dGk%less(1) - xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,L)
                !
                !
                !get d/dt G_k^\lmix = -xi*H(k,0)*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
                SxG(0:)=zero
                do kspin=1,Nspin
                   do korb=1,Norb
                      kk = iso_indx(kspin,korb)
                      SxG(0:) = SxG(0:) - xi*Hk(ii,kk)*Gk(kspin,jspin,korb,jorb)%lmix(1,0:)
                   enddo
                enddo
                dGk(ispin,jspin,iorb,jorb)%lmix(0:)= SxG(0:)
                !dGk%lmix(0:)= -xi*Hk*Gk%lmix(1,0:)
                !
                do j=0,L
                   !
                   SxG(0:)=zero
                   do s=0,j
                      do kspin=1,Nspin
                         do korb=1,Norb
                            SxG(s)=SxG(s)+Self(ispin,kspin,iorb,korb)%lmix(1,s)*Gk(kspin,jspin,korb,jorb)%mats(s+L-j)
                         enddo
                      enddo
                   end do
                   dGk(ispin,jspin,iorb,jorb)%lmix(j)=dGk(ispin,jspin,iorb,jorb)%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
                   ! do k=0,j
                   !    SxG(k)=self%lmix(1,k)*Gk%mats(k+L-j)
                   ! end do
                   ! dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
                   !
                   !
                   SxG(0:)=zero
                   do s=j,L
                      do kspin=1,Nspin
                         do korb=1,Norb
                            SxG(s)=SxG(s)+Self(ispin,kspin,iorb,korb)%lmix(1,s)*Gk(kspin,jspin,korb,jorb)%mats(s-j)
                         enddo
                      enddo
                   enddo
                   dGk(ispin,jspin,iorb,jorb)%lmix(j)=dGk(ispin,jspin,iorb,jorb)%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,L)
                   ! do k=j,L
                   !    SxG(k)=self%lmix(1,k)*Gk%mats(k-j)
                   ! end do
                   ! dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,L) 
                enddo
                !
             enddo
          enddo
       enddo
    enddo
    !
  contains
    !
    function iso_indx(ispin,iorb) result(iso)
      integer :: ispin
      integer :: iorb
      integer :: iso
      iso = iorb + (ispin-1)*Norb
    end function iso_indx
    !
  end subroutine neq_setup_contour_dgf_Nso






  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################









  subroutine neq_setup_contour_gf_main(params,g,file,gmats)
    type(kb_contour_params)             :: params
    type(kb_contour_gf)                 :: g
    character(len=*),optional           :: file
    complex(8),dimension(:),optional    :: gmats
    logical                             :: bool
    integer                             :: i,unit,N,L,Niw
    real(8),dimension(:),allocatable    :: wm,ftau
    !
    N = params%Nt
    L = params%Ntau
    Niw= params%Niw
    !
    call check_kb_contour_gf(params,G,"neq_setup_equilibrium_gf_main") 
    !
    if(.not.present(file).AND..not.present(gmats))stop "ERROR neq_setup_equilibrium_gf_main: can not setup G, no file or gmats given"

    if(present(file))then
       inquire(file=trim(file),exist=bool)
       if(.not.bool)stop "ERROR neq_setup_equilibrium_gf_main: file to read not present!"
       call read_array(reg(file),g%iw(:))
       !
    elseif(present(gmats))then
       call assert_shape(gmats,[Niw],"neq_setup_equilibrium_gf_main","gmats")
       g%iw = gmats
       !
    end if
    !
    !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
    allocate(ftau(0:Niw))
    call fft_gf_iw2tau(g%iw,ftau(0:),params%beta)
    call fft_extract_gtau(ftau(0:),g%mats(0:))
    g%less(1,1) = -xi*g%mats(L)
    g%ret(1,1)  = -xi
    g%lmix(1,0:L)=-xi*g%mats(L:0:-1)
    !
    deallocate(ftau)
    !
  end subroutine neq_setup_contour_gf_main

  subroutine neq_setup_contour_gf_Nso(params,g,file,gmats)
    type(kb_contour_params)                     :: params
    type(kb_contour_gf)                         :: g(:,:,:,:) ![Nspin][:][Norb][:]
    character(len=*),optional                   :: file
    complex(8),dimension(:,:,:,:,:),optional    :: gmats
    complex(8),dimension(:,:,:,:,:),allocatable :: Gread
    logical                                     :: bool
    integer                                     :: i,unit,N,L,Niw,Nlen,Nspin,Norb
    integer                                     :: ispin,jspin,iorb,jorb
    real(8),dimension(:),allocatable            :: ftau
    !
    N = params%Nt
    L = params%Ntau
    Niw= params%Niw
    !
    Nspin = size(G,1)
    Norb  = size(G,3)
    !
    call assert_shape_kb_contour_gf(G,[Nspin,Nspin,Norb,Norb],"neq_setup_equilibrium_gf_Nso","G")
    call check_kb_contour_gf(params,G,"neq_setup_equilibrium_gf_Nso","G") 
    !
    if(.not.present(file).AND..not.present(gmats))stop "ERROR neq_setup_equilibrium_gf_Nso: can not setup G, no file or gmats given"
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                g(ispin,jspin,iorb,jorb)=zero
             enddo
          enddo
       enddo
    enddo

    if(present(file))then
       inquire(file=trim(file),exist=bool)
       if(.not.bool)stop "ERROR neq_setup_equilibrium_gf_Nso: file to read not present!"
       allocate(Gread(Nspin,Nspin,Norb,Norb,Niw))
       call read_array(reg(file),Gread)
       do ispin=1,Nspin
          ! do jspin=1,Nspin
          jspin=ispin
          do iorb=1,Norb
             ! do jorb=1,Norb
             jorb=iorb
             g(ispin,jspin,iorb,jorb)%iw(:)=Gread(ispin,jspin,iorb,jorb,:)
             ! enddo
          enddo
          ! enddo
       enddo
       deallocate(Gread)
       !
    elseif(present(gmats))then
       call assert_shape(gmats,[Nspin,Nspin,Norb,Norb,Niw],"neq_setup_equilibrium_gf_Nso","gmats")
       do ispin=1,Nspin !do jspin=1,Nspin
          jspin=ispin
          do iorb=1,Norb !do jorb=1,Norb
             jorb=iorb
             g(ispin,jspin,iorb,jorb)%iw(:)=Gmats(ispin,jspin,iorb,jorb,:)
          enddo !enddo
       enddo !enddo
       !
    end if
    !
    !
    !INITIALIZE THE GF G^{x=M,<,R,\lmix}
    allocate(ftau(0:Niw))
    do ispin=1,Nspin !do jspin=1,Nspin
       jspin=ispin
       do iorb=1,Norb !do jorb=1,Norb
          jorb=iorb
          call fft_gf_iw2tau(g(ispin,jspin,iorb,jorb)%iw,ftau(0:),params%beta)
          call fft_extract_gtau(ftau,g(ispin,jspin,iorb,jorb)%mats)
          ! call fft_gf_iw2tau(g(ispin,jspin,iorb,jorb)%iw,g(ispin,jspin,iorb,jorb)%mats(0:),params%beta)
          !
          g(ispin,jspin,iorb,jorb)%less(1,1)  = -xi*g(ispin,jspin,iorb,jorb)%mats(L)
          g(ispin,jspin,iorb,jorb)%ret(1,1)   = -xi
          g(ispin,jspin,iorb,jorb)%lmix(1,0:L)= -xi*g(ispin,jspin,iorb,jorb)%mats(L:0:-1)!L-i)
          !
       enddo !enddo
    enddo !enddo
    !
    deallocate(ftau)
    !
  end subroutine neq_setup_contour_gf_Nso





  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################






  !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix} (this step depends on the imp. solv.)
  ! self^M(0,0) = -*U0*U0*G0(tau)*G0(-tau)*G0(tau)
  ! self^<(0,0) = i^3*U*U*G0(0-)*G0(0+)*G0(0-)
  ! self^>(0,0) = i^3*U*U*G0(0+)*G0(0-)*G0(0+)
  ! self^\lmix(0,t) = i^3*U*U0*G0(-t)*G0(t)*G0(-t)
  ! self^R(0,0) = self^> - self^<
  subroutine neq_setup_contour_sigma_main(params,g0,self)
    type(kb_contour_params)             :: params
    type(kb_contour_gf)                 :: g0
    type(kb_contour_gf)                 :: self
    logical                             :: bool
    integer                             :: i,j,k,ik,unit,len,N,L,Niw,Lk
    complex(8)                          :: Self_gtr
    real(8),dimension(:),allocatable    :: ftau,stau
    !
    N = params%Nt
    L = params%Ntau
    Niw= params%Niw
    !
    call check_kb_contour_gf(params,G0,"neq_setup_equilibrium_sigma_main") 
    call check_kb_contour_gf(params,Self,"neq_setup_equilibrium_sigma_main") 
    !
    !
    allocate(ftau(0:Niw),stau(0:Niw))
    !
    !Get G0(tau)
    call fft_gf_iw2tau(g0%iw,ftau(0:),params%beta)
    !
    !Build Sigma(tau) in IPT
    do j=0,Niw
       stau(j) = Ui*Ui*ftau(j)*ftau(Niw-j)*ftau(j)
    end do
    !
    !Get Sigma(iw)
    call fft_extract_gtau(stau,Self%mats)
    call fft_gf_tau2iw(Self%iw,stau,params%beta)
    if(Ui==0d0)Self%iw=zero
    !
    !Get Sigma on the contour
    Self%less(1,1)  = (xi**3)*U*U*g0%mats(L)*g0%mats(0)*g0%mats(L)
    Self_gtr        = (xi**3)*U*U*g0%mats(0)*g0%mats(L)*g0%mats(0)
    Self%lmix(1,0:) = (xi**3)*U*Ui*g0%mats(L:0:-1)*g0%mats(0:L)*g0%mats(L:0:-1)
    Self%ret(1,1)   = Self_gtr - Self%less(1,1)
    ! do j=0,L
    !    Self%lmix(1,j)=(xi**3)*U*Ui*g0%mats(L-j)*g0%mats(j)*g0%mats(L-j)
    ! end do
    ! Self%ret(1,1) = Self_gtr - Self%less(1,1)
    deallocate(ftau,stau)
    !
    return
  end subroutine neq_setup_contour_sigma_main

  subroutine neq_setup_contour_sigma_Nso(params,g0,self)
    type(kb_contour_params)             :: params
    type(kb_contour_gf)                 :: g0(:,:,:,:)
    type(kb_contour_gf)                 :: self(size(g0,1),size(g0,2),size(g0,3),size(g0,4))
    integer                             :: i,j,k,ik,unit,len,N,L,Niw,Nspin,Norb,ispin,iorb
    complex(8)                          :: Self_gtr
    real(8),dimension(:),allocatable    :: ftau,stau
    !
    N = params%Nt
    L = params%Ntau
    Niw= params%Niw
    !
    Nspin = size(g0,1)
    Norb  = size(g0,3)
    !
    call check_kb_contour_gf(params,G0,"neq_setup_equilibrium_sigma_main") 
    call check_kb_contour_gf(params,Self,"neq_setup_equilibrium_sigma_main") 
    !
    write(*,"(A)") "Get Contour Sigma at t=t`=0 using IPT. Diagonal case only."
    !
    allocate(ftau(0:Niw),stau(0:Niw))
    !
    do ispin=1,Nspin
       do iorb=1,Norb
          !
          !Get G0(tau)
          call fft_gf_iw2tau(g0(ispin,ispin,iorb,iorb)%iw,ftau(0:),params%beta)
          !
          !Build Sigma(tau) in IPT
          do j=0,Niw
             stau(j) = Ui*Ui*ftau(j)*ftau(Niw-j)*ftau(j)
          end do
          !
          !Get Sigma(iw)
          call fft_extract_gtau(stau,Self(ispin,ispin,iorb,iorb)%mats)
          call fft_gf_tau2iw(Self(ispin,ispin,iorb,iorb)%iw,stau,params%beta)
          if(Ui==0d0)Self(ispin,ispin,iorb,iorb)%iw=zero
          !
          !Get Sigma on the contour
          Self(ispin,ispin,iorb,iorb)%less(1,1)  = (xi**3)*U*U*g0(ispin,ispin,iorb,iorb)%mats(L)*g0(ispin,ispin,iorb,iorb)%mats(0)*g0(ispin,ispin,iorb,iorb)%mats(L)
          Self_gtr                               = (xi**3)*U*U*g0(ispin,ispin,iorb,iorb)%mats(0)*g0(ispin,ispin,iorb,iorb)%mats(L)*g0(ispin,ispin,iorb,iorb)%mats(0)
          Self(ispin,ispin,iorb,iorb)%lmix(1,0:) = (xi**3)*U*Ui*g0(ispin,ispin,iorb,iorb)%mats(L:0:-1)*g0(ispin,ispin,iorb,iorb)%mats(0:L)*g0(ispin,ispin,iorb,iorb)%mats(L:0:-1)
          Self(ispin,ispin,iorb,iorb)%ret(1,1)   = Self_gtr - Self(ispin,ispin,iorb,iorb)%less(1,1)
       enddo
    enddo
    !
    deallocate(ftau,stau)
    !
  end subroutine neq_setup_contour_sigma_Nso






  subroutine fft_extract_gtau(g,gred)
    real(8),dimension(0:) :: g
    real(8),dimension(0:) :: gred
    integer               :: N,Nred
    integer               :: i,ip
    real(8)               :: p,mismatch
    N   =size(g)-1
    Nred=size(gred)-1
    gred(0)=g(0)
    gred(Nred)=g(N)
    mismatch=dble(N)/dble(Nred)
    do i=1,Nred-1
       p=dble(i)*mismatch
       ip=int(p)
       gred(i)=g(ip)
    enddo
  end subroutine fft_extract_gtau








end module NEQ_IPT_SETUP





! subroutine neq_continue_equilibirum_(g0,gk,dgk,g,self,Hk,Wtk,params)
!   type(kb_contour_gf)                 :: g0
!   type(kb_contour_gf)                 :: gk(:)
!   type(kb_contour_dgf)                :: dgk(size(gk))
!   type(kb_contour_gf)                 :: g
!   type(kb_contour_gf)                 :: self
!   type(kb_contour_params)             :: params
!   real(8)                             :: Hk(size(gk)),Wtk(size(gk))
!   real(8)                             :: wm,res,ims
!   logical                             :: bool
!   integer                             :: i,j,k,ik,unit,len,N,L,Lf,Lk
!   complex(8)                          :: zeta
!   complex(8)                          :: Self_gtr
!   complex(8),allocatable,dimension(:) :: SxG
!   real(8)                             :: Scoeff(2),Gcoeff(4)
!   real(8),dimension(:),allocatable    :: ftau,stau
!   Lk=size(gk)
!   if(.not.g0%status)stop "init_functions: g0 is not allocated"
!   if(.not.g%status)stop "init_functions: g is not allocated"
!   do ik=1,Lk
!      if(.not.gk(ik)%status)stop "init_functions: gk(ik) is not allocated"
!      if(.not.dgk(ik)%status)stop "init_functions: dgk(ik) is not allocated"
!   enddo
!   if(.not.self%status)stop "init_functions: self is not allocated"
!   if(.not.params%status)stop "init_functions: params is not allocated"
!   !
!   N = params%Nt
!   L = params%Ntau
!   Lf= params%Niw
!   !
!   !CHECK IF G0(IW) IS AVAILABLE OR START FROM THE NON-INTERACTING SOLUTION
!   inquire(file=trim(g0file),exist=bool)
!   if(bool)then
!      write(*,"(A)")"Reading initial G0(iw) from file "//reg(g0file)
!      i = file_length(reg(g0file))
!      if(i/=Lf)then
!         print*,"init_equilibrium_weiss_field: Lfreq in +g0file does not correspond",i
!         stop
!      endif
!      unit = free_unit()
!      open(unit,file=reg(g0file),status='old')
!      do i=1,Lf
!         read(unit,*)wm,ims,res
!         g0%iw(i) = dcmplx(res,ims)
!      enddo
!      close(unit)
!   else
!      write(*,"(A)")"Start from Non-interacting G0(iw)"
!      do i=1,Lf
!         wm    = pi/beta*dble(2*i-1)
!         zeta  = dcmplx(0.d0,wm)
!         g0%iw(i) = sum_overk_zeta(zeta,Hk,Wtk)
!      enddo
!   endif
!   !
!   !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
!   ! Gcoeff = tail_coeff_glat(U,0.5d0,0d0,0d0)
!   ! call fft_gf_iw2tau(g0%iw,g0%mats(0:),params%beta,Gcoeff)
!   ! ! call fftgf_iw2tau(g0%iw,g0%mats(0:),params%beta)
!   allocate(ftau(0:Lf),stau(0:Lf))
!   call fft_gf_iw2tau(g0%iw,ftau(0:),params%beta)
!   call fft_extract_gtau(ftau,g0%mats)
!   g0%less(1,1) = -xi*g0%mats(L)
!   g0%ret(1,1)  = -xi
!   forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
!   !
!   !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
!   ! !(this step depends on the imp. solv.)
!   ! ! self^M(0,0) = -*U0*U0*G0(tau)*G0(-tau)*G0(tau)
!   ! ! self^<(0,0) = i^3*U*U*G0(0-)*G0(0+)*G0(0-)
!   ! ! self^>(0,0) = i^3*U*U*G0(0+)*G0(0-)*G0(0+)
!   ! ! self^\lmix(0,t) = i^3*U*U0*G0(-t)*G0(t)*G0(-t)
!   ! ! self^R(0,0) = self^> - self^<    
!   do j=0,Lf
!      stau(j) = Ui*Ui*ftau(j)*ftau(Lf-j)*ftau(j)
!   end do
!   call fft_extract_gtau(stau,Self%mats)
!   call fft_gf_tau2iw(Self%iw,stau,beta)
!   if(Ui==0d0)Self%iw=zero
!   Self%less(1,1)=(xi**3)*U*U*g0%mats(L)*g0%mats(0)*g0%mats(L)
!   Self_gtr      =(xi**3)*U*U*g0%mats(0)*g0%mats(L)*g0%mats(0)
!   do j=0,L
!      Self%lmix(1,j)=(xi**3)*U*Ui*g0%mats(L-j)*g0%mats(j)*g0%mats(L-j)
!   end do
!   Self%ret(1,1) = Self_gtr - Self%less(1,1)
!   deallocate(ftau,stau)
!   !
!   G%mats=0.d0
!   G%iw = zero
!   do ik=1,Lk 
!      call neq_setup_initial_conditions(gk(ik),dgk(ik),self,hk(ik),params)
!      G%mats(0:)  = G%mats(0:)  + wtk(ik)*gk(ik)%mats(0:)
!      G%iw(:)     = G%iw(:)     + wtk(ik)*gk(ik)%iw(:)
!      G%ret(1,1)  = G%ret(1,1)  + wtk(ik)*gk(ik)%ret(1,1)
!      G%less(1,1) = G%less(1,1) + wtk(ik)*gk(ik)%less(1,1)
!      G%lmix(1,0:)= G%lmix(1,0:)+ wtk(ik)*gk(ik)%lmix(1,0:)
!   enddo
!   return
! end subroutine neq_continue_equilibirum_


! subroutine neq_continue_equilibirum_bethe(g0,dg0,g,self,params,wband)
!   type(kb_contour_gf)                 :: g0
!   type(kb_contour_dgf)                :: dg0
!   type(kb_contour_gf)                 :: g
!   type(kb_contour_gf)                 :: self
!   type(kb_contour_params)             :: params
!   real(8),optional                    :: wband
!   real(8)                             :: wband_
!   real(8)                             :: wm,res,ims
!   logical                             :: bool
!   integer                             :: i,j,k,ik,unit,len,N,L,Lf
!   complex(8)                          :: zeta
!   complex(8)                          :: Self_gtr
!   complex(8),allocatable,dimension(:) :: GxG0
!   real(8)                             :: Scoeff(2),Gcoeff(4)
!   wband_=1d0;if(present(wband))wband_=wband
!   if(.not.g0%status)stop "init_functions: g0 is not allocated"
!   if(.not.dg0%status)stop "init_functions: dg0 is not allocated"
!   if(.not.g%status)stop "init_functions: g is not allocated"
!   if(.not.self%status)stop "init_functions: self is not allocated"
!   if(.not.params%status)stop "init_functions: params is not allocated"
!   N = params%Nt
!   L = params%Ntau
!   Lf= params%Niw
!   !
!   !CHECK IF G0(IW) IS AVAILABLE OR START FROM THE NON-INTERACTING SOLUTION
!   inquire(file=trim(g0file),exist=bool)
!   if(bool)then
!      write(*,"(A)")"Reading initial G0(iw) from file "//reg(g0file)
!      unit = free_unit()
!      open(unit,file=reg(g0file),status='old')
!      i = file_length(reg(g0file))
!      if(i/=Lf)then
!         print*,"init_equilibrium_weiss_field: Lfreq in +g0file does not correspond",i
!         stop
!      endif
!      do i=1,Lf
!         read(unit,*)wm,ims,res
!         g0%iw(i) = dcmplx(res,ims)
!      enddo
!      close(unit)
!   else
!      write(*,"(A)")"Start from Non-interacting G0(iw)"
!      do i=1,Lf
!         wm    = pi/beta*dble(2*i-1)
!         zeta  = xi*wm
!         g0%iw(i) = gfbethe(wm,zeta,wband)
!      enddo
!   endif
!   !
!   !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
!   Gcoeff = tail_coeff_glat(U,0.5d0,0d0,0d0)
!   call fft_gf_iw2tau(g0%iw,g0%mats(0:),params%beta,Gcoeff)
!   g0%less(1,1) = -xi*g0%mats(L)
!   g0%ret(1,1)  = -xi
!   forall(i=0:L)g0%lmix(1,i)=-xi*g0%mats(L-i)
!   !
!   !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
!   !(this step depends on the imp. solv.)
!   ! self^M(0,0) = -*U0*U0*G0(tau)*G0(-tau)*G0(tau)
!   ! self^<(0,0) = i^3*U*U*G0(0-)*G0(0+)*G0(0-)
!   ! self^>(0,0) = i^3*U*U*G0(0+)*G0(0-)*G0(0+)
!   ! self^\lmix(0,t) = i^3*U*U0*G0(-t)*G0(t)*G0(-t)
!   ! self^R(0,0) = self^> - self^<
!   do j=0,L
!      Self%mats(j) = Ui*Ui*g0%mats(j)*g0%mats(L-j)*g0%mats(j)
!   end do
!   Scoeff  = tail_coeff_sigma(Ui,0.5d0)
!   call fft_sigma_tau2iw(Self%iw,Self%mats(0:),beta,Scoeff)
!   Self%iw = xi*dimag(self%iw) !!ACTHUNG: imposing half-filling symmetry
!   Self%less(1,1)=(xi**3)*U*U*g0%mats(L)*g0%mats(0)*g0%mats(L)
!   Self_gtr      =(xi**3)*U*U*g0%mats(0)*g0%mats(L)*g0%mats(0)
!   do j=0,L
!      Self%lmix(1,j)=(xi**3)*U*Ui*g0%mats(L-j)*g0%mats(j)*g0%mats(L-j)
!   end do
!   Self%ret(1,1) = Self_gtr - Self%less(1,1)
!   !
!   !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
!   do i=1,Lf
!      wm    = pi/beta*dble(2*i-1)
!      zeta  = xi*wm - Self%iw(i)
!      G%iw(i) = gfbethe(wm,zeta,wband)
!   enddo
!   Gcoeff      = tail_coeff_glat(U,0.5d0,0d0,0d0)
!   call fft_gf_iw2tau(G%iw,G%mats,beta,Gcoeff)     !get G(tau)
!   G%ret(1,1)  = -xi                               !get G^R(0,0)=-xi
!   G%less(1,1) = -xi*G%mats(L)                  !get G^<(0,0)= xi*G^M(0-)
!   forall(i=0:L)G%lmix(1,i)=-xi*G%mats(L-i) !get G^\lmix(0,tau)=-xi*G(beta-tau>0)
!   !Derivatives
!   allocate(GxG0(0:L))
!   !
!   !get d/dt G0^R = 0.d0
!   dG0%ret(1)  = zero
!   !get d/dt G0^< = -xi(-xi)int_0^beta G^\lmix * G0^\rmix
!   do k=0,L
!      GxG0(k)=G%lmix(1,k)*conjg(G0%lmix(1,L-k))
!   end do
!   dG0%less(1) = -xi*(-xi)*params%dtau*kb_trapz(GxG0(0:),0,L) 
!   !get d/dt G0^\lmix = -xi*int_0^beta G0^\lmix * G0^M
!   dG0%lmix(0:)= zero
!   do j=0,L
!      do k=0,j
!         GxG0(k)=G%lmix(1,k)*G0%mats(k+L-j)
!      end do
!      dG0%lmix(j)=dG0%lmix(j)+xi*params%dtau*kb_trapz(GxG0(0:),0,j)
!      do k=j,L
!         GxG0(k)=G%lmix(1,k)*G0%mats(k-j)
!      end do
!      dG0%lmix(j)=dG0%lmix(j)-xi*params%dtau*kb_trapz(GxG0(0:),j,L) 
!   enddo
!   deallocate(GxG0)
!   return
! end subroutine neq_continue_equilibirum_bethe



! subroutine neq_setup_initial_conditions(Gk,dGk,Self,Hk,params)
!   type(kb_contour_gf)                 :: Gk,Self
!   type(kb_contour_dgf)                :: dGk
!   real(8)                             :: Hk
!   type(kb_contour_params)             :: params
!   integer                             :: i,j,k,L,Niw
!   real(8)                             :: nk
!   complex(8)                          :: epsk
!   complex(8),allocatable,dimension(:) :: SxG
!   real(8),dimension(:),allocatable    :: ftau
!   !
!   L = params%Ntau
!   Niw= params%Niw
!   !
!   do i=1,Niw
!      Gk%iw(i) = one/(xi*params%wm(i) - Hk - Self%iw(i))
!   enddo
!   ! call fft_gf_iw2tau(Gk%iw,Gk%mats(0:),beta)
!   allocate(ftau(0:Niw))
!   call fft_gf_iw2tau(Gk%iw,ftau(0:),beta)          !get G_k(tau)
!   call fft_extract_gtau(ftau,Gk%mats)
!   Gk%less(1,1) = -xi*Gk%mats(L)                 !get G^<_k(0,0)= xi*G^M_k(0-)
!   Gk%ret(1,1)  = -xi                               !get G^R_k(0,0)=-xi
!   forall(i=0:L)Gk%lmix(1,i)=-xi*Gk%mats(L-i) !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
!   !
!   !Derivatives
!   allocate(SxG(0:L))
!   !get d/dt G_k^R = -i e(k,0)G_k^R
!   dGk%ret(1)  = -xi*Hk*Gk%ret(1,1)
!   !get d/dt G_k^< = -i e(k,0)G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
!   do k=0,L
!      SxG(k)=self%lmix(1,k)*conjg(Gk%lmix(1,L-k))
!   end do
!   dGk%less(1) = -xi*Hk*Gk%less(1,1)-&
!        xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,L) 
!   !get d/dt G_k^\lmix = -xi*e(k,0)*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
!   dGk%lmix(0:)= -xi*Hk*Gk%lmix(1,0:)
!   do j=0,L
!      do k=0,j
!         SxG(k)=self%lmix(1,k)*Gk%mats(k+L-j)
!      end do
!      dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
!      do k=j,L
!         SxG(k)=self%lmix(1,k)*Gk%mats(k-j)
!      end do
!      dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,L) 
!   enddo
! end subroutine neq_setup_initial_conditions



! !+-------------------------------------------------------------------+
! !PURPOSE: setup the Weiss Field G0 for the next time-step
! !+-------------------------------------------------------------------+
! subroutine neq_setup_weiss_field(g0,params)
!   type(kb_contour_gf)                   :: g0
!   type(kb_contour_params)               :: params
!   integer                               :: i,j,k,N,L
!   complex(8),allocatable,dimension(:)   :: SxG
!   if(.not.g0%status)stop "init_g0: g0 is not allocated"
!   if(.not.params%status)stop "init_g0: params is not allocated"
!   !
!   N = params%Nt
!   L = params%Ntau
!   !
!   select case(N)
!   case(1)
!      return
!   case(2)
!      !GUESS G0 AT THE NEXT STEP, GIVEN THE INITIAL CONDITIONS
!      do j=1,N
!         g0%ret(N,j) =g0%ret(1,1)
!         g0%less(N,j)=g0%less(1,1)
!      end do
!      do i=1,N-1
!         g0%less(i,N)=g0%less(1,1)
!      end do
!      do j=0,L
!         g0%lmix(N,j)=g0%lmix(1,j)
!      end do
!   case default
!      !EXTEND G0 FROM THE [N-1,N-1] TO THE [N,N] SQUARE TO START DMFT
!      !USING QUADRATIC EXTRAPOLATION
!      do k=1,N-1
!         g0%less(N,k)=2.d0*g0%less(N-1,k)-g0%less(N-2,k)
!         g0%less(k,N)=2.d0*g0%less(k,N-1)-g0%less(k,N-2)
!      end do
!      g0%less(N,N)=2.d0*g0%less(N-1,N-1)-g0%less(N-2,N-2)
!      !
!      do k=0,L
!         g0%lmix(N,k)=2.d0*g0%lmix(N-1,k)-g0%lmix(N-2,k)
!      end do
!      !
!      g0%ret(N,N)=-xi
!      do k=1,N-2
!         g0%ret(N,k)=2.d0*g0%ret(N-1,k)-g0%ret(N-2,k)
!      end do
!      g0%ret(N,N-1)=0.5d0*(g0%ret(N,N)+g0%ret(N,N-2))
!   end select
! end subroutine neq_setup_weiss_field



