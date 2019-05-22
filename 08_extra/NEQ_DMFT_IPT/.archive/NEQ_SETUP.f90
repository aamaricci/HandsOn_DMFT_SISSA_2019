module NEQ_SETUP
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE NEQ_INPUT_VARS  
  USE SF_CONSTANTS
  USE SF_IOTOOLS
  USE SF_SPECIAL
  USE DMFT_FFTGF
  USE DMFT_FFTAUX
  USE DMFT_MISC
  private


  interface neq_setup_contour
     module procedure :: neq_setup_contour_gf
     module procedure :: neq_setup_contour_sigma
     module procedure :: neq_setup_contour_dGk
     module procedure :: neq_setup_contour_dG0
  end interface neq_setup_contour
  public :: neq_setup_contour



contains


  !+-------------------------------------------------------------------+
  !PURPOSE: setup initial conditions for:
  ! - Green's functions
  ! - Self-energy
  ! - derivative of the k-resolved GF: dG_k
  ! - derivative of the Weiss Field: dG0  
  !+-------------------------------------------------------------------+
  subroutine neq_setup_contour_gf(giw,green,params)
    complex(8),dimension(:)          :: giw
    type(kb_contour_gf)              :: green
    type(kb_contour_params)          :: params
    logical                          :: bool
    integer                          :: i,j,k,ik,unit,len,N,L,Lf,Lk
    real(8),allocatable,dimension(:) :: gtau
    !
    if(.not.params%status)stop "neq_push_to_contour_gf ERROR: params is not allocated"
    if(.not.green%status)stop "neq_push_to_contour_gf ERROR: Green is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    if(size(giw)/=Lf)stop "neq_push_to_contour_gf ERROR: size(Giw) != Niw "
    !
    !Push the G(iw) component
    Green%iw = Giw
    !Push the G(tau) --> G.mats component
    allocate(gtau(0:Lf))
    call fft_gf_iw2tau(Green%iw,gtau(0:),params%beta)
    call fft_extract_gtau(gtau,green%mats(0:))
    !Push the G^<(0,0) component
    Green%less(1,1)  = -xi*Green%mats(L)
    !Push the G^Ret(0,0) component
    Green%ret(1,1)   = -xi
    !Push the G^\lmix(0,tau) component
    Green%lmix(1,0:L)= -xi*Green%mats(L:0:-1)
  end subroutine neq_setup_contour_gf

  subroutine neq_setup_contour_sigma(shf,siw,sigma,params)
    complex(8)                       :: shf
    complex(8),dimension(:)          :: siw
    type(kb_contour_sigma)           :: sigma
    type(kb_contour_params)          :: params
    logical                          :: bool
    integer                          :: i,j,k,ik,unit,len,N,L,Lf,Lk
    real(8),allocatable,dimension(:) :: stau
    !
    if(.not.params%status)stop "neq_push_to_contour_gf ERROR: params is not allocated"
    if(.not.sigma%status)stop "neq_push_to_contour_gf ERROR: Sigma is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    !
    if(size(siw)/=Lf)stop "neq_push_to_contour_gf ERROR: size(Siw) != Niw "
    !
    !Push the Sigma(iw) component
    Sigma%hf(1)  = Shf
    Sigma%reg%iw = Siw
    !Push the Sigma(tau) --> Sigma.mats component
    allocate(Stau(0:Lf))
    call fft_sigma_iw2tau(Sigma%reg%iw,stau(0:),params%beta)
    call fft_extract_gtau(Stau,Sigma%reg%mats(0:))
    !Push the G^<(0,0) component
    Sigma%reg%less(1,1)  = -xi*Sigma%reg%mats(L)
    !Push the G^Ret(0,0) component
    Sigma%reg%ret(1,1)   = xi*(Sigma%reg%mats(0)+Sigma%reg%mats(L))  !OK
    !Push the G^\lmix(0,tau) component
    Sigma%reg%lmix(1,0:L)= -xi*Sigma%reg%mats(L:0:-1)
  end subroutine neq_setup_contour_sigma

  subroutine neq_setup_contour_dGk(Gk,dGk,Self,Hk,params)
    type(kb_contour_gf)                 :: Gk
    type(kb_contour_dgf)                :: dGk
    type(kb_contour_sigma)              :: Self
    complex(8)                          :: Hk,Hk_hf
    type(kb_contour_params)             :: params
    integer                             :: i,j,k,L,Lf
    complex(8),allocatable,dimension(:) :: SxG
    !
    L = params%Ntau
    Lf= params%Niw
    !
    allocate(SxG(0:L))
    !
    Hk_hf = Hk - Self%hf(1)
    !
    !Derivatives
    !get d/dt G_k^R = -i [e(k,0)-Sigma_HF(0)]G_k^R
    dGk%ret(1)  = -xi*Hk_hf*Gk%ret(1,1)
    !
    !get d/dt G_k^< = -i [e(k,0)-Sigma_HF(0)]G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
    SxG(0:)=zero
    do k=0,L
       SxG(k)=self%reg%lmix(1,k)*conjg(Gk%lmix(1,L-k))
    end do
    dGk%less(1) = -xi*Hk_hf*Gk%less(1,1)-xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,L) 
    !
    !get d/dt G_k^\lmix = -xi*[e(k,0)-Sigma_HF(0)]*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
    dGk%lmix(0:)= -xi*Hk_hf*Gk%lmix(1,0:)
    do j=0,L
       SxG(0:)=zero
       do k=0,j
          SxG(k)=Self%reg%lmix(1,k)*Gk%mats(k+L-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
       SxG(0:)=zero
       do k=j,L
          SxG(k)=Self%reg%lmix(1,k)*Gk%mats(k-j)
       end do
       dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,L) 
    enddo
    !
    deallocate(SxG)
  end subroutine neq_setup_contour_dGk

  subroutine neq_setup_contour_dG0(g0,dg0,g,self,params,wband)
    type(kb_contour_gf)                 :: g0
    type(kb_contour_dgf)                :: dg0
    type(kb_contour_gf)                 :: g
    type(kb_contour_sigma)              :: self
    type(kb_contour_params)             :: params
    real(8),optional                    :: wband
    real(8)                             :: wband_
    !
    integer                             :: i,j,k,ik,len,N,L,Lf
    complex(8),allocatable,dimension(:) :: GxG0
    !
    wband_=1d0;if(present(wband))wband_=wband
    if(.not.g0%status)stop "init_functions: g0 is not allocated"
    if(.not.dg0%status)stop "init_functions: dg0 is not allocated"
    if(.not.g%status)stop "init_functions: g is not allocated"
    if(.not.self%status)stop "init_functions: self is not allocated"
    if(.not.params%status)stop "init_functions: params is not allocated"
    !
    N = params%Nt
    L = params%Ntau
    Lf= params%Niw
    allocate( GxG0(0:L) )
    !
    !
    !Derivatives d/dt G0:
    !get d/dt G0^R = 0.d0
    dG0%ret(1)  = zero
    !get d/dt G0^< = -xi(-xi)int_0^beta G^\lmix * G0^\rmix
    GxG0(0:)=zero
    do k=0,L
       GxG0(k)=G%lmix(1,k)*conjg(G0%lmix(1,L-k))
    end do
    dG0%less(1) = -xi*(-xi)*params%dtau*kb_trapz(GxG0(0:),0,L) 
    !get d/dt G0^\lmix = -xi*int_0^beta G0^\lmix * G0^M
    dG0%lmix(0:)= zero    
    do j=0,L
       GxG0(0:)=zero
       do k=0,j
          GxG0(k)=G%lmix(1,k)*G0%mats(k+L-j)
       end do
       dG0%lmix(j)=dG0%lmix(j)+xi*params%dtau*kb_trapz(GxG0(0:),0,j)
       GxG0(0:)=zero
       do k=j,L
          GxG0(k)=G%lmix(1,k)*G0%mats(k-j)
       end do
       dG0%lmix(j)=dG0%lmix(j)-xi*params%dtau*kb_trapz(GxG0(0:),j,L) 
    enddo
    !
    deallocate(GxG0)
    return
  end subroutine neq_setup_contour_dG0



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



end module NEQ_SETUP









! !+-------------------------------------------------------------------+
! !PURPOSE: Read equilibrium solution and initialize the corresponding
! ! function.
! ! The code expects to read Sigma_reg(tau), that is the regular part 
! ! of the self-energy in imaginary time. Regular here means deprived 
! ! of the Hartree-Fock or Hartree-Fock-Bogoliubov term (1st order).
! ! The latter term is reconstructed analytically from the knowledge 
! ! of the static observables in the header of the self-energy files:
! ! NORMAL SIGMA:
! ! # U_i n [with U_i = interactino strenght at equilibrium, n=density
! ! ....
! !+-------------------------------------------------------------------+
! subroutine read_equilibrium_sigma_normal(self,params)
!   type(kb_contour_sigma)           :: self
!   type(kb_contour_params)          :: params
!   real(8)                          :: tau
!   logical                          :: bool
!   complex(8)                       :: zeta
!   integer                          :: i,j,k,ik,unit,Len,N,L,Lf,Ltau
!   real(8)                          :: u_,dens
!   real(8),dimension(0:1)           :: Scoeff
!   real(8),allocatable,dimension(:) :: Stau !dummy function for FFT to tau
!   !
!   if(.not.self%status)  stop "read_equilibrium_sigma: sigma is not allocated"
!   if(.not.params%status)stop "read_equilibrium_sigma: params is not allocated"
!   !
!   N = params%Nt
!   L = params%Ntau
!   Lf= params%Niw
!   !
!   inquire(file=reg(sigma_file),exist=bool)
!   if(bool)then
!      write(*,"(A)")"read_equilibrium_sigma_normal: reading Sigma(tau) from file "//reg(sigma_file)
!      dens = read_header(reg(sigma_file))
!      Len  = file_length(reg(sigma_file))
!      Ltau = Len-1-1           !-1 for the header, -1 for the 0
!      if(Ltau < L)stop "read_equilibrium_sigma_normal error: Ltau < params%Ntau"
!      allocate(Stau(0:Ltau))
!      unit = free_unit()
!      open(unit,file=reg(sigma_file),status='old')
!      do i=0,Ltau
!         read(unit,*)tau,Stau(i)
!      enddo
!      close(unit)
!      call fft_extract_gtau(Stau(0:),Self%reg%mats(0:))       !<=== Get Sigma(tau)= xtract(tmp_selt(tau_))
!      self%hf(1) = Ui*(dens-0.5d0)                            !<=== Get Sigma_HF  = Re[Sigma(iw-->infty)]
!      Scoeff  = tail_coeff_sigma(Ui,dens)
!      call fft_sigma_tau2iw(Self%reg%iw,Stau(0:),beta,Scoeff) !<=== Get Sigma(iw) = fft(tmp_self(tau_))
!      deallocate(Stau)
!   else
!      write(*,"(A)")"read_equilibrium_sigma_normal: start from Hartree-Fock Sigma(iw)=Ui*(n-1/2)"
!      dens=0.5d0
!      self%hf(1)   = Ui*(dens-0.5d0)
!      self%reg%iw  = zero
!      self%reg%mats= zero
!   endif
!   self%reg%less(1,1)  = -xi*self%reg%mats(L)                     !OK
!   self%reg%ret(1,1)   =  xi*(self%reg%mats(0)+self%reg%mats(L))  !OK
!   self%reg%lmix(1,0:L)= -xi*self%reg%mats(L:0:-1)                !small errors near 0,beta
! end subroutine read_equilibrium_sigma_normal




! !+-------------------------------------------------------------------+
! !PURPOSE: setup initial conditions for k-resolved GF
! !+-------------------------------------------------------------------+
! subroutine neq_setup_initial_conditions_normal(Gk,dGk,Self,Hk,params)
!   type(kb_contour_gf)                 :: Gk
!   type(kb_contour_dgf)                :: dGk
!   type(kb_contour_sigma)              :: Self
!   complex(8)                          :: Hk
!   ! complex(8)                          :: Hk_hf
!   type(kb_contour_params)             :: params
!   integer                             :: i,j,k,L,Lf
!   complex(8),allocatable,dimension(:) :: SxG
!   ! real(8),allocatable,dimension(:)    :: Stau !dummy function for FFT to tau
!   ! complex(8),allocatable,dimension(:) :: Siw  !dummy function in iw_n
!   !
!   L = params%Ntau
!   Lf= params%Niw
!   ! allocate(Siw(Lf),Stau(0:Lf))
!   !
!   !This is already done using neq_push_to_contour new procedures:
!   !so we pass a Gk which is already continued to the contour.
!   ! Siw   = self%reg%iw + self%hf(1)
!   ! Hk_hf = Hk - self%hf(1)
!   ! do i=1,Lf
!   !    Gk%iw(i) = one/(xi*params%wm(i) - Hk - Siw(i))
!   ! enddo
!   ! call fft_gf_iw2tau(Gk%iw,Stau(0:),beta) !get G_k(tau)
!   ! call fft_extract_gtau(Stau,Gk%mats)     !
!   ! Gk%less(1,1)  = -xi*Gk%mats(L)          !get G^<_k(0,0)= xi*G^M_k(0-)
!   ! Gk%ret(1,1)   = -xi                     !get G^R_k(0,0)=-xi
!   ! Gk%lmix(1,0:L)= -xi*Gk%mats(L:0:-1)     !get G^\lmix_k(0,tau)=xi*G_k(tau<0)=-xi*G_k(beta-tau>0)
!   !
!   !
!   allocate(SxG(0:L))
!   !
!   !Derivatives
!   !get d/dt G_k^R = -i [e(k,0)-Sigma_HF(0)]G_k^R
!   dGk%ret(1)  = -xi*Hk_hf*Gk%ret(1,1)
!   !
!   !get d/dt G_k^< = -i [e(k,0)-Sigma_HF(0)]G_k^< -xi(-xi)int_0^beta S^\lmix*G_k^\rmix
!   SxG(0:)=zero
!   do k=0,L
!      SxG(k)=self%reg%lmix(1,k)*conjg(Gk%lmix(1,L-k))
!   end do
!   dGk%less(1) = -xi*Hk_hf*Gk%less(1,1)-xi*(-xi)*params%dtau*kb_trapz(SxG(0:),0,L) 
!   !
!   !get d/dt G_k^\lmix = -xi*[e(k,0)-Sigma_HF(0)]*G_k^\lmix - xi*int_0^beta G_k^\lmix*G_k^M
!   dGk%lmix(0:)= -xi*Hk_hf*Gk%lmix(1,0:)
!   do j=0,L
!      SxG(0:)=zero
!      do k=0,j
!         SxG(k)=Self%reg%lmix(1,k)*Gk%mats(k+L-j)
!      end do
!      dGk%lmix(j)=dGk%lmix(j)+xi*params%dtau*kb_trapz(SxG(0:),0,j)
!      SxG(0:)=zero
!      do k=j,L
!         SxG(k)=Self%reg%lmix(1,k)*Gk%mats(k-j)
!      end do
!      dGk%lmix(j)=dGk%lmix(j)-xi*params%dtau*kb_trapz(SxG(0:),j,L) 
!   enddo
!   !
!   ! deallocate(Siw,Stau,SxG)
!   deallocate(SxG)
! end subroutine neq_setup_initial_conditions_normal





! !+-------------------------------------------------------------------+
! !PURPOSE: continue the NORMAL equilibrium to Keldysh contour
! !+-------------------------------------------------------------------+
! subroutine neq_continue_equilibirum_normal_lattice(g0,gk,dgk,g,self,Hk,Wtk,params)
!   type(kb_contour_gf)                 :: g0
!   type(kb_contour_gf)                 :: gk(:)
!   type(kb_contour_dgf)                :: dgk(size(gk))
!   type(kb_contour_gf)                 :: g
!   type(kb_contour_sigma)              :: self
!   complex(8)                          :: Hk(size(gk))
!   real(8)                             :: Wtk(size(gk))
!   type(kb_contour_params)             :: params
!   logical                             :: bool
!   integer                             :: i,j,k,ik,unit,len,N,L,Lf,Lk
!   real(8),allocatable,dimension(:)    :: G0tau !dummy function for FFT to tau
!   complex(8),allocatable,dimension(:) :: Siw  !dummy function in iw_n
!   !
!   Lk=size(gk)
!   N = params%Nt
!   L = params%Ntau
!   Lf= params%Niw
!   allocate( Siw(Lf),G0tau(0:Lf) )
!   !
!   if(.not.g0%status)stop "init_functions: g0 is not allocated"
!   if(.not.g%status)stop "init_functions: g is not allocated"
!   if(.not.self%status)stop "init_functions: self is not allocated"
!   do ik=1,Lk
!      if(.not.gk(ik)%status)stop "init_functions: gk(ik) is not allocated"
!      if(.not.dgk(ik)%status)stop "init_functions: dgk(ik) is not allocated"
!   enddo
!   if(.not.params%status)stop "init_functions: params is not allocated"
!   !
!   !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
!   call read_equilibrium_sigma_normal(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
!   !
!   !INITIALIZE THE LOCAL GREEN'S FUNCTION Gloc^{x=M,<,R,\lmix}
!   g=zero
!   do ik=1,Lk
!      call neq_setup_initial_conditions_normal(gk(ik),dgk(ik),self,hk(ik),params)
!   enddo
!   call sum_kb_contour_gf(gk(:),wtk(:),g,params)
!   !
!   !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
!   Siw   = self%reg%iw + self%hf(1)       !Sigma = Sigma_reg + SigmaHF
!   ! g0%iw = one/( one/g%iw + Siw )       !Dyson: G0^-1 = Gloc^-1 + Sigma
!   g0%iw = one/( one/g%iw + self%reg%iw ) !Dyson: tildeG0^-1 = Gloc^-1 + Sigma_reg = Gloc^-1 + Sigma - Sigma_HF
!   call fft_gf_iw2tau(g0%iw,G0tau(0:),params%beta)
!   call fft_extract_gtau(G0tau,g0%mats)
!   g0%less(1,1)  = -xi*g0%mats(L)
!   g0%ret(1,1)   = -xi
!   g0%lmix(1,0:L)= -xi*g0%mats(L:0:-1)
!   !
!   deallocate(G0tau,Siw)
! end subroutine neq_continue_equilibirum_normal_lattice


! !+-----------------------------------------------------------------------------+!
! !PURPOSE: AUXILIARY ROUTINES:
! !  - read_header     : read the header of the seed file for the Self-energy S(tau)
! !  - fft_extract_gtau: extract a G(tau) from a dense to a sparse set.
! !+-----------------------------------------------------------------------------+!
! function read_header(file)  result(obs)
!   character(len=*) :: file
!   integer          :: unit
!   real(8)          :: u_,obs
!   unit = free_unit()
!   open(unit,file=file,status='old')
!   read(unit,*)u_,obs
!   close(unit)
!   if(u_/=Ui)then
!      print*,"read_header error: U_eq in "//file//" different from the input:",Ui
!      stop
!   endif
!   write(*,"(4A)")"Header of the file:",reg(txtfy(u_))," ",reg(txtfy(obs))
! end function read_header



! subroutine neq_continue_equilibirum_normal_bethe(g0,dg0,g,self,params,wband)
!   type(kb_contour_gf)                 :: g0
!   type(kb_contour_dgf)                :: dg0
!   type(kb_contour_gf)                 :: g
!   type(kb_contour_sigma)              :: self
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
!   real(8),allocatable,dimension(:)    :: Gtau
!   complex(8),allocatable,dimension(:) :: Siw
!   !
!   wband_=1d0;if(present(wband))wband_=wband
!   if(.not.g0%status)stop "init_functions: g0 is not allocated"
!   if(.not.dg0%status)stop "init_functions: dg0 is not allocated"
!   if(.not.g%status)stop "init_functions: g is not allocated"
!   if(.not.self%status)stop "init_functions: self is not allocated"
!   if(.not.params%status)stop "init_functions: params is not allocated"
!   !
!   N = params%Nt
!   L = params%Ntau
!   Lf= params%Niw
!   allocate( Siw(Lf), Gtau(0:L), GxG0(0:L) )
!   !
!   !INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix}
!   call read_equilibrium_sigma_normal(self,params)            !<== get Sigma^{x=iw,tau,M,<,R,\lmix}
!   !
!   !INITIALIZE THE GREEN'S FUNCTION G^{x=M,<,R,\lmix}
!   Siw = self%reg%iw + self%hf(1)
!   do i=1,Lf
!      wm    = pi/beta*dble(2*i-1)
!      zeta  = xi*wm - Siw(i)
!      G%iw(i) = gfbethe(wm,zeta,wband)
!   enddo
!   !Gcoeff      = tail_coeff_glat(Ui,dens,xmu,hloc=0d0)
!   call fft_gf_iw2tau(G%iw,Gtau(0:),beta,Gcoeff)     !get G(tau)
!   call fft_extract_gtau(Gtau,G%mats)
!   G%ret(1,1)   = -xi                !get G^R(0,0)=-xi
!   G%less(1,1)  = -xi*G%mats(L)      !get G^<(0,0)= xi*G^M(0-)
!   G%lmix(1,0:L)= -xi*G%mats(L:0:-1) !get G^\lmix(0,tau)=-xi*G(beta-tau>0)
!   !
!   !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
!   g0%iw = one/(one/g%iw + self%reg%iw) !Dyson: tideG0^-1 = Gloc^-1 + Sigma_reg
!   call fft_gf_iw2tau(g0%iw,Gtau(0:),params%beta)
!   call fft_extract_gtau(Gtau,g0%mats)
!   g0%less(1,1)  = -xi*g0%mats(L)
!   g0%ret(1,1)   = -xi
!   g0%lmix(1,0:L)= -xi*g0%mats(L:0:-1)
!   !
!   !Derivatives d/dt G0:
!   !get d/dt G0^R = 0.d0
!   dG0%ret(1)  = zero
!   !get d/dt G0^< = -xi(-xi)int_0^beta G^\lmix * G0^\rmix
!   GxG0(0:)=zero
!   do k=0,L
!      GxG0(k)=G%lmix(1,k)*conjg(G0%lmix(1,L-k))
!   end do
!   dG0%less(1) = -xi*(-xi)*params%dtau*kb_trapz(GxG0(0:),0,L) 
!   !get d/dt G0^\lmix = -xi*int_0^beta G0^\lmix * G0^M
!   dG0%lmix(0:)= zero

!   do j=0,L
!      GxG0(0:)=zero
!      do k=0,j
!         GxG0(k)=G%lmix(1,k)*G0%mats(k+L-j)
!      end do
!      dG0%lmix(j)=dG0%lmix(j)+xi*params%dtau*kb_trapz(GxG0(0:),0,j)
!      GxG0(0:)=zero
!      do k=j,L
!         GxG0(k)=G%lmix(1,k)*G0%mats(k-j)
!      end do
!      dG0%lmix(j)=dG0%lmix(j)-xi*params%dtau*kb_trapz(GxG0(0:),j,L) 
!   enddo
!   !
!   deallocate(Siw,Gtau,GxG0)
!   return
! end subroutine neq_continue_equilibirum_normal_bethe
