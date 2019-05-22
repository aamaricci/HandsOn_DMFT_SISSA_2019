subroutine neq_push_equilibrium_weiss_main(params,g0,gweiss)
  type(kb_contour_params)          :: params
  type(kb_contour_gf)              :: g0
  complex(8),dimension(:)          :: gweiss
  real(8)                          :: wm,res,ims
  logical                          :: bool
  integer                          :: i,unit,N,L,Niw,Nlen
  real(8),dimension(:),allocatable :: wm,ftau
  !
  N = params%Nt
  L = params%Ntau
  Niw= params%Niw
  !
  call check_dimension_kb_contour(params,G0,"neq_setup_equilibrium_weiss_main") 
  !
  call assert_shape(Gweiss,[Niw],"neq_setup_equilibrium_weiss_main","Gweiss")
  write(*,"(A)")"Get G0(iw) from input array Gweiss"
  g0%iw = gweiss
  !
  !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
  allocate(ftau(0:Niw))
  call fft_gf_iw2tau(g0%iw,ftau(0:),params%beta)
  call fft_extract_gtau(ftau,g0%mats)
  !
  g0%less(1,1)  = -xi*g0%mats(L)
  g0%ret(1,1)   = -xi
  g0%lmix(1,0:L)= -xi*g0%mats(L:0:-1)!L-i)
  ! forall(i=0:L)  g0%lmix(1,i)=-xi*g0%mats(L-i)
  !
  deallocate(ftau)
  !
  weiss_status = .true.
  !
end subroutine neq_push_equilibrium_weiss_main

subroutine neq_push_equilibrium_weiss_Nso(params,g0,gweiss)
  type(kb_contour_params)          :: params
  type(kb_contour_gf)              :: g0(:,:,:,:) ![Nspin][:][Norb][:]
  complex(8),dimension(:,:,:,:,:)  :: gweiss
  real(8)                          :: wm,res,ims
  logical                          :: bool
  integer                          :: i,unit,N,L,Niw,Nlen,Nspin,Norb
  integer                          :: ispin,jspin,iorb,jorb
  real(8),dimension(:),allocatable :: ftau
  !
  N = params%Nt
  L = params%Ntau
  Niw= params%Niw
  !
  Nspin = size(G0,1)
  Norb  = size(G0,3)
  !
  call check_dimension_kb_contour(params,G0,"neq_setup_equilibrium_weiss_main") 
  !
  forall(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)g0(ispin,jspin,iorb,jorb)=zero
  !
  call assert_shape(Gweiss,[Nspin,Nspin,Norb,Norb,Niw],"neq_setup_equilibrium_weiss_main","Gweiss")
  write(*,"(A)")"Get diagonal G0(iw) from input array Gweiss"
  forall(ispin=1:Nspin,iorb=1:Norb)g0(ispin,ispin,iorb,iorb)%iw=Gweiss(ispin,ispin,iorb,iorb,:)
  !
  !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
  allocate(ftau(0:Niw))
  do ispin=1,Nspin
     do iorb=1,Norb
        call fft_gf_iw2tau(g0(ispin,ispin,iorb,iorb)%iw,ftau(0:),params%beta)
        call fft_extract_gtau(ftau,g0(ispin,ispin,iorb,iorb)%mats)
        !
        g0(ispin,ispin,iorb,iorb)%less(1,1)  = -xi*g0(ispin,ispin,iorb,iorb)%mats(L)
        g0(ispin,ispin,iorb,iorb)%ret(1,1)   = -xi
        g0(ispin,ispin,iorb,iorb)%lmix(1,0:L)= -xi*g0(ispin,ispin,iorb,iorb)%mats(L:0:-1)!L-i)
        !
     enddo
  enddo
  !
  deallocate(ftau)
  !
  weiss_status = .true.
  !
end subroutine neq_push_equilibrium_weiss_Nso





subroutine neq_setup_equilibrium_weiss_main(params,g0,gweiss)
  type(kb_contour_params)             :: params
  type(kb_contour_gf)                 :: g0
  complex(8),dimension(:),optional    :: gweiss
  real(8)                             :: wm,res,ims
  logical                             :: bool
  integer                             :: i,unit,N,L,Niw,Nlen
  real(8),dimension(:),allocatable    :: wm,ftau
  !
  N = params%Nt
  L = params%Ntau
  Niw= params%Niw
  !
  call check_dimension_kb_contour(params,G0,"neq_setup_equilibrium_weiss_main") 
  !
  if(present(gweiss))then       !If user pass a specific weiss field function:
     call assert_shape(Gweiss,[Niw],"neq_setup_equilibrium_weiss_main","Gweiss")
     write(*,"(A)")"Get G0(iw) from input array Gweiss"
     g0%iw = gweiss
     !
  else                          !Read it from file. If file does exist exit:
     !
     ! inquire(file=trim(g0file),exist=bool)
     ! if(.not.bool)stop "ERROR neq_setup_equilibrium_weiss_main: g0file not present!"
     !
     ! write(*,"(A)")"Get G0(iw) from file w,ImG0,ReG0: "//reg(g0file)
     ! Nlen = file_length(reg(g0file))
     ! if(Nlen/=Niw)stop "ERROR neq_setup_equilibrium_weiss_main: len(g0file)!=Niw"
     ! open(free_unit(unit),file=reg(g0file),status='old')
     ! do i=1,Niw
     !    read(unit,*)wm,ims,res
     !    g0%iw(i) = dcmplx(res,ims)
     ! enddo
     ! close(unit)
     call read_array(reg(g0file),Gweiss)
     ! else        
     !    write(*,"(A)")"Start from Non-interacting G0(iw)"
     !    do i=1,Niw
     !       wm    = pi/beta*dble(2*i-1)
     !       zeta  = dcmplx(0.d0,wm)
     !       g0%iw(i) = sum_overk_zeta(zeta,Hk,Wtk)
     !    enddo
  endif
  !
  !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
  allocate(ftau(0:Niw))
  call fft_gf_iw2tau(g0%iw,ftau(0:),params%beta)
  call fft_extract_gtau(ftau,g0%mats)
  !
  g0%less(1,1)  = -xi*g0%mats(L)
  g0%ret(1,1)   = -xi
  g0%lmix(1,0:L)= -xi*g0%mats(L:0:-1)!L-i)
  ! forall(i=0:L)  g0%lmix(1,i)=-xi*g0%mats(L-i)
  !
  deallocate(ftau)
  !
end subroutine neq_setup_equilibrium_weiss_main

subroutine neq_setup_equilibrium_weiss_Nso(params,g0,gweiss)
  type(kb_contour_params)                  :: params
  type(kb_contour_gf)                      :: g0(:,:,:,:) ![Nspin][:][Norb][:]
  complex(8),dimension(:,:,:,:,:),optional :: gweiss
  real(8)                                  :: wm,res,ims
  logical                                  :: bool
  integer                                  :: i,unit,N,L,Niw,Nlen,Nspin,Norb
  integer                                  :: ispin,jspin,iorb,jorb
  real(8),dimension(:),allocatable         :: ftau
  !
  N = params%Nt
  L = params%Ntau
  Niw= params%Niw
  !
  Nspin = size(G0,1)
  Norb  = size(G0,3)
  !
  call check_dimension_kb_contour(params,G0,"neq_setup_equilibrium_weiss_main") 
  !
  forall(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)g0(ispin,jspin,iorb,jorb)=zero
  !
  if(present(gweiss))then       !If user pass a specific weiss field function:     
     call assert_shape(Gweiss,[Nspin,Nspin,Norb,Norb,Niw],"neq_setup_equilibrium_weiss_main","Gweiss")
     write(*,"(A)")"Get diagonal G0(iw) from input array Gweiss"
     forall(ispin=1:Nspin,iorb=1:Norb)g0(ispin,ispin,iorb,iorb)%iw=Gweiss(ispin,ispin,iorb,iorb,:)
     !
  else                          !Read it from file. If file does exist exit:
     !
     ! inquire(file=trim(g0file),exist=bool)
     ! if(.not.bool)stop "ERROR neq_setup_equilibrium_weiss_main: g0file not present!"
     ! write(*,"(A)")"Get G0(iw) from file w,ImG0,ReG0: "//reg(g0file)
     ! Nlen = file_length(reg(g0file))
     ! if(Nlen/=Niw)stop "ERROR neq_setup_equilibrium_weiss_main: len(g0file)!=Niw"
     ! open(free_unit(unit),file=reg(g0file),status='old')
     ! do i=1,Niw
     !    read(unit,*)wm,ims,res
     !    g0%iw(i) = dcmplx(res,ims)
     ! enddo
     ! close(unit)
     !
     call read_array(reg(g0file),Gweiss)
     !     
  endif
  !
  !INITIALIZE THE WEISS FIELD G0^{x=M,<,R,\lmix}
  allocate(ftau(0:Niw))
  do ispin=1,Nspin
     do iorb=1,Norb
        call fft_gf_iw2tau(g0(ispin,ispin,iorb,iorb)%iw,ftau(0:),params%beta)
        call fft_extract_gtau(ftau,g0(ispin,ispin,iorb,iorb)%mats)
        !
        g0(ispin,ispin,iorb,iorb)%less(1,1)  = -xi*g0(ispin,ispin,iorb,iorb)%mats(L)
        g0(ispin,ispin,iorb,iorb)%ret(1,1)   = -xi
        g0(ispin,ispin,iorb,iorb)%lmix(1,0:L)= -xi*g0(ispin,ispin,iorb,iorb)%mats(L:0:-1)!L-i)
        !
     enddo
  enddo
  deallocate(ftau)
  !
end subroutine neq_setup_equilibrium_weiss_Nso





!INITIALIZE THE SELF-ENERGY SELF^{x=M,<,R,\lmix} (this step depends on the imp. solv.)
! self^M(0,0) = -*U0*U0*G0(tau)*G0(-tau)*G0(tau)
! self^<(0,0) = i^3*U*U*G0(0-)*G0(0+)*G0(0-)
! self^>(0,0) = i^3*U*U*G0(0+)*G0(0-)*G0(0+)
! self^\lmix(0,t) = i^3*U*U0*G0(-t)*G0(t)*G0(-t)
! self^R(0,0) = self^> - self^<
subroutine neq_setup_equilibrium_sigma_main(params,g0,self)
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
  call check_dimension_kb_contour(params,G0,"neq_setup_equilibrium_sigma_main") 
  call check_dimension_kb_contour(params,Self,"neq_setup_equilibrium_sigma_main") 
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
end subroutine neq_setup_equilibrium_sigma_main

subroutine neq_setup_equilibrium_sigma_Nso(params,g0,self)
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
  call check_dimension_kb_contour(params,G0,"neq_setup_equilibrium_sigma_main") 
  call check_dimension_kb_contour(params,Self,"neq_setup_equilibrium_sigma_main") 
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
end subroutine neq_setup_equilibrium_sigma_Nso
