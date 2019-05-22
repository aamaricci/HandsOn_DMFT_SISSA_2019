!----------------------------------------------------------------------------
!  This subroutine solves a Volterra integral equation of the second kind,
!              G(t,t') = Q(t,t') + (K*G)(t,t')
!  for t=n*dt or t'=n*dt, using 2^nd *implicit* Runge-Kutta method.
!----------------------------------------------------------------------------
!Q = Q(t,t'); K = contour_gf
subroutine vie_kb_contour_gf_q(Q,K,G,params)
  type(kb_contour_gf),intent(in)      :: Q
  type(kb_contour_gf),intent(in)      :: K
  type(kb_contour_gf),intent(inout)   :: G
  type(kb_contour_params),intent(in)  :: params
  integer                             :: N,L
  real(8)                             :: dt,dtau
  integer                             :: i,j,s,itau,jtau
  complex(8),dimension(:),allocatable :: KxG
  !
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  allocate(KxG(0:max(N,L)))
  !
  ! Ret component
  ! G^R(t,t') - \int_{t'}^t K^R(t,s)*G^R(s,t')ds = Q^R(t,t')
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  G%ret(N,N)=Q%ret(N,N)
  do j=1,N-1
     G%ret(N,j)=Q%ret(N,j)
     KxG(0:)=zero
     do s=j,N-1
        KxG(s)=K%ret(N,s)*G%ret(s,j)
     eNd do
     G%ret(N,j)=G%ret(N,j) + dt*kb_half_trapz(KxG(0:),j,N-1)
     G%ret(N,j)=G%ret(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
  end do
  !
  ! Lmix component
  ! G^\lmix(t,tau') - \int_0^t K^R(t,s)*G^\lmix(s,tau')ds
  !    = Q^\lmix(t,tau') + \int_0^\beta K^\lmix(t,s)*G^M(s,tau')ds
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  do jtau=0,L
     G%lmix(N,jtau)=Q%lmix(N,jtau)
     !
     KxG(0:)=zero
     do s=0,jtau
        KxG(s)=K%lmix(N,s)*G%mats(s+L-jtau)
     end do
     G%lmix(N,jtau)=G%lmix(N,jtau)-dtau*kb_trapz(KxG(0:),0,jtau)
     KxG(0:)=zero
     do s=jtau,L
        KxG(s)=K%lmix(N,s)*G%mats(s-jtau)
     end do
     G%lmix(N,jtau)=G%lmix(N,jtau)+dtau*kb_trapz(KxG(0:),jtau,L)
     !
     KxG(0:)=zero
     do s=1,N-1
        KxG(s)=K%ret(N,s)*G%lmix(s,jtau)
     end do
     G%lmix(N,jtau)=G%lmix(N,jtau) + dt*kb_half_trapz(KxG(0:),1,N-1)
     !
     G%lmix(N,jtau)=G%lmix(N,jtau)/(1.d0-0.5d0*dt*K%ret(N,N))
  end do
  !
  ! Less component
  ! G^<(t,t') - \int_0^t K^{R}(t,s)*G^{<}(s,t')ds
  !    = Q^<(t,t') - i\int_0^\beta K^\lmix(t,s)*G^\rmix(s,t')ds
  !      + \int_0^{t'} K^<(t,s)*G^A(s,t')ds
  ! G^<(t_{N},t_{j})
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  do j=1, N-1
     G%less(N,j)=Q%less(N,j)
     !
     KxG(0:)=zero
     do s=0,L
        KxG(s)=K%lmix(N,s)*get_rmix(G,s,j,L)
     enddo
     G%less(N,j)=G%less(N,j)-xi*dtau*kb_trapz(KxG(0:),0,L)
     !
     KxG(0:)=zero
     do s=1,j
        KxG(s)=K%less(N,s)*get_adv(G,s,j)
     enddo
     G%less(N,j)=G%less(N,j)+dt*kb_trapz(KxG(0:),1,j)
     !
     KxG(0:)=zero
     do s=1,N-1
        KxG(s)=K%ret(N,s)*G%less(s,j)
     enddo
     G%less(N,j)=G%less(N,j) + dt*kb_half_trapz(KxG(0:),1,N-1)
     !
     G%less(N,j)=G%less(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
  end do
  !
  ! G^<(t_{i},t_{N}) <= Hermite conjugate
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! if (.not.G%anomalous) then
  do i=1,N-1
     G%less(i,N) = -conjg(G%less(N,i))
  end do
  ! else
  !    do i=1,N-1
  !       G%less(i,N) = G%less(N,i)+G%ret(N,i)
  !    end do
  ! endif
  !
  ! G^{<}(t_{n},t_{n}) <= Diagonal
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  G%less(N,N)=Q%less(N,N)
  !
  KxG(0:)=zero
  do s=0,L
     KxG(s)=K%lmix(N,s)*get_rmix(G,s,N,L)
  end do
  G%less(N,N)=G%less(N,N)-xi*dtau*kb_trapz(KxG(0:),0,L)
  !
  KxG(0:)=zero
  do s=1,N
     KxG(s)=K%less(N,s)*get_adv(G,s,N)
  end do
  G%less(N,N)=G%less(N,N)+dt*kb_trapz(KxG(0:),1,N)
  !
  KxG(0:)=zero
  do s=1,N-1
     KxG(s)=K%ret(N,s)*G%less(s,N)
  end do
  G%less(N,N)=G%less(N,N) + dt*kb_half_trapz(KxG(0:),1,N-1)
  !
  G%less(N,N)=G%less(N,N)/(1.d0-0.5d0*dt*K%ret(N,N))
  !
  deallocate(KxG)
end subroutine vie_kb_contour_gf_q





!Q=Q(t,t'); K=contour_sigma
subroutine vie_kb_contour_gf_sigma(Q,Self,G,params)
  type(kb_contour_gf),intent(in)      :: Q
  type(kb_contour_sigma),intent(in)   :: Self
  type(kb_contour_gf),intent(inout)   :: G
  type(kb_contour_params),intent(in)  :: params
  type(kb_contour_gf)                 :: K
  integer                             :: N,L
  real(8)                             :: dt,dtau
  integer                             :: i,j,s,itau,jtau
  complex(8),dimension(:),allocatable :: KxG
  !
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  allocate(KxG(0:max(N,L)))
  call allocate_kb_contour_gf(K,params)
  K = Self%reg
  !
  ! Ret component
  ! G^R(t,t') - \int_{t'}^t K^R(t,s)*G^R(s,t')ds = Q^R(t,t')
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  G%ret(N,N)=Q%ret(N,N)
  do j=1,N-1
     G%ret(N,j)=Q%ret(N,j)
     KxG(0:)=zero
     do s=j,N-1
        KxG(s)=K%ret(N,s)*G%ret(s,j)
     eNd do
     G%ret(N,j)=G%ret(N,j) + dt*kb_half_trapz(KxG(0:),j,N-1)
     G%ret(N,j)=G%ret(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
  end do
  !
  ! Lmix component
  ! G^\lmix(t,tau') - \int_0^t K^R(t,s)*G^\lmix(s,tau')ds
  !    = Q^\lmix(t,tau') + \int_0^\beta K^\lmix(t,s)*G^M(s,tau')ds
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  do jtau=0,L
     G%lmix(N,jtau)=Q%lmix(N,jtau)
     !
     KxG(0:)=zero
     do s=0,jtau
        KxG(s)=K%lmix(N,s)*G%mats(s+L-jtau)
     end do
     G%lmix(N,jtau)=G%lmix(N,jtau)-dtau*kb_trapz(KxG(0:),0,jtau)
     KxG(0:)=zero
     do s=jtau,L
        KxG(s)=K%lmix(N,s)*G%mats(s-jtau)
     end do
     G%lmix(N,jtau)=G%lmix(N,jtau)+dtau*kb_trapz(KxG(0:),jtau,L)
     !
     KxG(0:)=zero
     do s=1,N-1
        KxG(s)=K%ret(N,s)*G%lmix(s,jtau)
     end do
     G%lmix(N,jtau)=G%lmix(N,jtau) + dt*kb_half_trapz(KxG(0:),1,N-1)
     !
     G%lmix(N,jtau)=G%lmix(N,jtau)/(1.d0-0.5d0*dt*K%ret(N,N))
  end do
  !
  ! Less component
  ! G^<(t,t') - \int_0^t K^{R}(t,s)*G^{<}(s,t')ds
  !    = Q^<(t,t') - i\int_0^\beta K^\lmix(t,s)*G^\rmix(s,t')ds
  !      + \int_0^{t'} K^<(t,s)*G^A(s,t')ds
  ! G^<(t_{N},t_{j})
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  do j=1, N-1
     G%less(N,j)=Q%less(N,j)
     !
     KxG(0:)=zero
     do s=0,L
        KxG(s)=K%lmix(N,s)*get_rmix(G,s,j,L)
     enddo
     G%less(N,j)=G%less(N,j)-xi*dtau*kb_trapz(KxG(0:),0,L)
     !
     KxG(0:)=zero
     do s=1,j
        KxG(s)=K%less(N,s)*get_adv(G,s,j)
     enddo
     G%less(N,j)=G%less(N,j)+dt*kb_trapz(KxG(0:),1,j)
     !
     KxG(0:)=zero
     do s=1,N-1
        KxG(s)=K%ret(N,s)*G%less(s,j)
     enddo
     G%less(N,j)=G%less(N,j) + dt*kb_half_trapz(KxG(0:),1,N-1)
     !
     G%less(N,j)=G%less(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
  end do
  !
  ! G^<(t_{i},t_{N}) <= Hermite conjugate
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! if (.not.G%anomalous) then
  do i=1,N-1
     G%less(i,N) = -conjg(G%less(N,i))
  end do
  ! else
  !    do i=1,N-1
  !       G%less(i,N) = G%less(N,i)+G%ret(N,i)
  !    end do
  ! endif
  !
  ! G^{<}(t_{n},t_{n}) <= Diagonal
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  G%less(N,N)=Q%less(N,N)
  !
  KxG(0:)=zero
  do s=0,L
     KxG(s)=K%lmix(N,s)*get_rmix(G,s,N,L)
  end do
  G%less(N,N)=G%less(N,N)-xi*dtau*kb_trapz(KxG(0:),0,L)
  !
  KxG(0:)=zero
  do s=1,N
     KxG(s)=K%less(N,s)*get_adv(G,s,N)
  end do
  G%less(N,N)=G%less(N,N)+dt*kb_trapz(KxG(0:),1,N)
  !
  KxG(0:)=zero
  do s=1,N-1
     KxG(s)=K%ret(N,s)*G%less(s,N)
  end do
  G%less(N,N)=G%less(N,N) + dt*kb_half_trapz(KxG(0:),1,N-1)
  !
  G%less(N,N)=G%less(N,N)/(1.d0-0.5d0*dt*K%ret(N,N))
  !
  deallocate(KxG)
  call deallocate_kb_contour_gf(K)
end subroutine vie_kb_contour_gf_sigma




!Q=delta;K=contour_gf
subroutine vie_kb_contour_gf_delta(K,G,params)
  type(kb_contour_gf),intent(in)      :: K
  type(kb_contour_gf),intent(inout)   :: G
  type(kb_contour_params),intent(in)  :: params
  integer                             :: N,L
  real(8)                             :: dt,dtau
  integer                             :: i,j,s,itau,jtau
  complex(8),dimension(:),allocatable :: KxG
  !
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  allocate(KxG(0:max(N,L)))
  !
  ! Ret component
  ! G^R(t,t') - \int_{t'}^t K^R(t,s)*G^R(s,t')ds = Q^R(t,t')
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  G%ret(N,N)=one                !Q%ret(N,N)
  do j=1,N-1
     G%ret(N,j)=zero            !Q%ret(N,j)
     !
     KxG(0:)=zero
     do s=j,N-1
        KxG(s)=K%ret(N,s)*G%ret(s,j)
     eNd do
     G%ret(N,j)=G%ret(N,j) + dt*kb_half_trapz(KxG(0:),j,N-1)
     !
     G%ret(N,j)=G%ret(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
  end do
  !
  ! Lmix component
  ! G^\lmix(t,tau') - \int_0^t K^R(t,s)*G^\lmix(s,tau')ds
  !    = Q^\lmix(t,tau') + \int_0^\beta K^\lmix(t,s)*G^M(s,tau')ds
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  do jtau=0,L
     G%lmix(N,jtau)=zero        !Q%lmix(N,jtau)
     !
     KxG(0:)=zero
     do s=0,jtau
        KxG(s)=K%lmix(N,s)*G%mats(s+L-jtau)
     end do
     G%lmix(N,jtau)=G%lmix(N,jtau)-dtau*kb_trapz(KxG(0:),0,jtau)
     !
     KxG(0:)=zero
     do s=jtau,L
        KxG(s)=K%lmix(N,s)*G%mats(s-jtau)
     end do
     G%lmix(N,jtau)=G%lmix(N,jtau)+dtau*kb_trapz(KxG(0:),jtau,L)
     !
     KxG(0:)=zero
     do s=1,N-1
        KxG(s)=K%ret(N,s)*G%lmix(s,jtau)
     end do
     G%lmix(N,jtau)=G%lmix(N,jtau) + dt*kb_half_trapz(KxG(0:),1,N-1)
     !
     G%lmix(N,jtau)=G%lmix(N,jtau)/(1.d0-0.5d0*dt*K%ret(N,N))
  end do
  !
  ! Less component
  ! G^<(t,t') - \int_0^t K^{R}(t,s)*G^{<}(s,t')ds
  !    = Q^<(t,t') - i\int_0^\beta K^\lmix(t,s)*G^\rmix(s,t')ds
  !      + \int_0^{t'} K^<(t,s)*G^A(s,t')ds
  ! G^<(t_{N},t_{j})
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  do j=1, N-1
     G%less(N,j)=zero           !Q%less(N,j)
     !
     KxG(0:)=zero
     do s=0,L
        KxG(s)=K%lmix(N,s)*get_rmix(G,s,j,L)
     enddo
     G%less(N,j)=G%less(N,j)-xi*dtau*kb_trapz(KxG(0:),0,L)
     !
     KxG(0:)=zero
     do s=1,j
        KxG(s)=K%less(N,s)*get_adv(G,s,j)
     enddo
     G%less(N,j)=G%less(N,j)+dt*kb_trapz(KxG(0:),1,j)
     !
     KxG(0:)=zero
     do s=1,N-1
        KxG(s)=K%ret(N,s)*G%less(s,j)
     enddo
     G%less(N,j)=G%less(N,j) + dt*kb_half_trapz(KxG(0:),1,N-1)
     !
     G%less(N,j)=G%less(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
  end do
  !
  ! G^<(t_{i},t_{N}) <= Hermite conjugate
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! if (.not.G%anomalous) then
  do i=1,N-1
     G%less(i,N) = -conjg(G%less(N,i))
  end do
  ! else
  !    do i=1,N-1
  !       G%less(i,N) = G%less(N,i)+G%ret(N,i)
  !    end do
  ! endif
  !
  ! G^{<}(t_{n},t_{n}) <= Diagonal
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  G%less(N,N)=zero              !Q%less(N,N)
  !
  KxG(0:)=zero
  do s=0,L
     KxG(s)=K%lmix(N,s)*get_rmix(G,s,N,L)
  end do
  G%less(N,N)=G%less(N,N)-xi*dtau*kb_trapz(KxG(0:),0,L)
  !
  KxG(0:)=zero
  do s=1,N
     KxG(s)=K%less(N,s)*get_adv(G,s,N)
  end do
  G%less(N,N)=G%less(N,N)+dt*kb_trapz(KxG(0:),1,N)
  !
  KxG(0:)=zero
  do s=1,N-1
     KxG(s)=K%ret(N,s)*G%less(s,N)
  end do
  G%less(N,N)=G%less(N,N) + dt*kb_half_trapz(KxG(0:),1,N-1)
  !
  G%less(N,N)=G%less(N,N)/(1.d0-0.5d0*dt*K%ret(N,N))
  !
  deallocate(KxG)
end subroutine vie_kb_contour_gf_delta
