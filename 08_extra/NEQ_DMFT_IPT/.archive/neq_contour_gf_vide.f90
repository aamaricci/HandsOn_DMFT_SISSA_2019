!----------------------------------------------------------------------------
!  This subroutine solves a Volterra integro-differential equation of 
!  the second kind,
!              [i*d/dt-h(t)]G(t,t') = delta(t,t') + (K*G)(t,t'),
!  for t=n*dt or t'=n*dt, using 2^nd implicit Runge-Kutta method.
!----------------------------------------------------------------------------
!K = contour_GF
!H = H_0 + Sigma_HF
subroutine vide_kb_contour_gf_K(H,K,G,dG,dG_new,params)
  complex(8),dimension(:)             :: H
  type(kb_contour_gf),intent(in)      :: K
  type(kb_contour_gf),intent(inout)   :: G
  type(kb_contour_dgf),intent(inout)  :: dG
  type(kb_contour_dgf)                :: dG_new
  type(kb_contour_params),intent(in)  :: params
  integer                             :: N,L
  real(8)                             :: dt,dtau
  integer                             :: i,j,itau,jtau,s
  complex(8),dimension(:),allocatable :: KxG
  complex(8)                          :: dG_less
  !
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  allocate(KxG(0:max(N,L)))
  !
  !TYPE: X=RET,LESS,LMIX (R,<,\LMIX)
  !i d/dt G^x(t,:)  = h(t)G^x(t,:) + Q^x(t,:) + int_{0,t'}^{t}K^R(t,s)G^x(s,:)ds
  !
  !
  !ret component
  ! d/dt G^R(t,:) = -i*h(t)*G^R(t,:) -i*delta(t,:) - i\int_{:}^t K^R(t,s)*G^R(s,:)ds
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
  G%ret(N,N)=-xi
  dG_new%ret(N)=-xi*H(N)*G%ret(N,N)
  do j=1,N-1
     G%ret(N,j)=G%ret(N-1,j) + 0.5d0*dt*dG%ret(j)
     !
     KxG(0:)=zero
     do s=j,N-1
        KxG(s)=K%ret(N,s)*G%ret(s,j)
     enddo
     dG_new%ret(j)=-xi*dt*kb_half_trapz(KxG(0:),j,N-1)
     !
     G%ret(N,j)=G%ret(N,j) + 0.5d0*dt*dG_new%ret(j)
     G%ret(N,j)=G%ret(N,j)/(1.d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
     !
     dG_new%ret(j)=dG_new%ret(j) - xi*H(N)*G%ret(N,j) - 0.5d0*xi*dt*K%ret(N,N)*G%ret(N,j)
  enddo
  !
  !Lmix component
  !d/dt G^\lmix(t,:) = -i*H(t)*G^\lmix(t,:) -i\int_0^\beta K^\lmix(t,s)*G^M(s,:)ds -i\int_0^t K^R(t,s)G^\lmix(s,:)ds
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
  do jtau=0,L
     G%lmix(N,jtau)=G%lmix(N-1,jtau)+0.5d0*dt*dG%lmix(jtau)
     !
     KxG(0:)=zero
     do s=0,jtau
        KxG(s)=K%lmix(N,s)*G%mats(s+L-jtau)
     end do
     dG_new%lmix(jtau)=xi*dtau*kb_trapz(KxG(0:),0,jtau)
     KxG(0:)=zero
     do s=jtau,L
        KxG(s)=K%lmix(N,s)*G%mats(s-jtau)
     end do
     dG_new%lmix(jtau)=dG_new%lmix(jtau)-xi*dtau*kb_trapz(KxG(0:),jtau,L)!<= add -iQ(t)
     !
     KxG(0:)=zero
     do s=1,N-1
        KxG(s)=K%ret(N,s)*G%lmix(s,jtau)
     end do
     dG_new%lmix(jtau)=dG_new%lmix(jtau)-xi*dt*kb_half_trapz(KxG(0:),1,N-1)
     !
     G%lmix(N,jtau)=G%lmix(N,jtau) + 0.5d0*dt*dG_new%lmix(jtau)
     G%lmix(N,jtau)=G%lmix(N,jtau)/(1.d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
     !
     dG_new%lmix(jtau)=dG_new%lmix(jtau)-xi*H(N)*G%lmix(N,jtau)-0.5d0*xi*dt*K%ret(N,N)*G%lmix(N,jtau)
  end do
  !
  !Less component
  !d/dt G^<(t,:) = -i*H(t)*G^<(t,:) -i*[ (-i)*\int_0^\beta K^\lmix(t,s)*G^\rmix(s,:)ds + \int_0^{:}K^<(t,s)*G^A(s,:)ds ]
  !                                 -i*\int_0^t K^R(t,s)*G^<(s,:)ds
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
  ! G^<(t_{N},t_{j}), d/dt G^<(t_{N},t_{j}) <== lower-right triangle
  do j=1,N-1
     G%less(N,j)=G%less(N-1,j) + 0.5d0*dt*dG%less(j)
     !
     KxG(0:)=zero
     do s=0,L
        KxG(s)=K%lmix(N,s)*get_rmix(G,s,j,L)
     end do
     dG_new%less(j)=-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
     !
     KxG(0:)=zero
     do s=1,j
        KxG(s)=K%less(N,s)*get_adv(G,s,j)
     end do
     dG_new%less(j)=dG_new%less(j)-xi*dt*kb_trapz(KxG(0:),1,j)!<= -iQ(t)
     !
     KxG(0:)=zero
     do s=1,N-1
        KxG(s)=K%ret(N,s)*G%less(s,j)
     end do
     dG_new%less(j)=dG_new%less(j)-xi*dt*kb_half_trapz(KxG(0:),1,N-1)
     !
     G%less(N,j)=G%less(N,j) + 0.5d0*dt*dG_new%less(j)
     G%less(N,j)=G%less(N,j)/(1.d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
     !
     dG_new%less(j)=dG_new%less(j)-xi*H(N)*G%less(N,j)-0.5d0*xi*dt*K%ret(N,N)*G%less(N,j)
  end do
  !
  ! G^<(t_{i},t_{N}), d/dt G^<(t_{i},t_{N}) <== upper left triangle
  ! Hermitian conjugate G
  ! if (.not.G%anomalous) then
  do i=1,N-1
     G%less(i,N)=-conjg(G%less(N,i))
  end do
  ! else
  !    do i=1,N-1
  !       G%less(i,N)=G%less(N,i)+G%ret(N,i)
  !    end do
  ! endif
  !
  ! d/dt G^<(t_{N-1},t_{N})
  dG_less=-xi*H(N-1)*G%less(N-1,N)
  !
  KxG(0:)=zero
  do s=0,L
     KxG(s)=K%lmix(N-1,s)*get_rmix(G,s,N,L)
  end do
  dG_less=dG_less-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
  !
  KxG(0:)=zero
  do s=1,N
     KxG(s)=K%less(N-1,s)*get_adv(G,s,N)
  end do
  dG_less=dG_less-xi*dt*kb_trapz(KxG(0:),1,N)
  !
  KxG(0:)=zero
  do s=1,N-1
     KxG(s)=K%ret(N-1,s)*G%less(s,N)
  end do
  dG_less=dG_less-xi*dt*kb_trapz(KxG(0:),1,N-1)
  !
  !G^<(N,N), d/dt G^<(N,N)
  G%less(N,N)=G%less(N-1,N)+0.5d0*dt*dG_less
  !
  KxG(0:)=zero
  do s=0,L
     KxG(s)=K%lmix(N,s)*get_rmix(G,s,N,L)
  end do
  dG_new%less(N)=-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
  !
  KxG(0:)=zero
  do s=1,N
     KxG(s)=K%less(N,s)*get_adv(G,s,N)
  end do
  dG_new%less(N)=dG_new%less(N)-xi*dt*kb_trapz(KxG(0:),1,N)
  !
  KxG(0:)=zero
  do s=1,N-1
     KxG(s)=K%ret(N,s)*G%less(s,N)
  end do
  dG_new%less(N)=dG_new%less(N)-xi*dt*kb_half_trapz(KxG(0:),1,N-1)
  !
  G%less(N,N)=G%less(N,N)+0.5d0*dt*dG_new%less(N)
  G%less(N,N)=G%less(N,N)/(1.d0+0.5d0*xi*dt*H(N)+0.25d0*xi*dt**2*K%ret(N,N))
  !
  dG_new%less(N)=dG_new%less(N)-xi*H(N)*G%less(N,N)-0.5d0*xi*dt*K%ret(N,N)*G%less(N,N)
  !
  deallocate(KxG)
end subroutine vide_kb_contour_gf_K

!K = contour_Sigma
!H = H_0
subroutine vide_kb_contour_gf_Sigma(H0,Sigma,G,dG,dG_new,params)
  complex(8),dimension(:)             :: H0
  type(kb_contour_sigma),intent(in)   :: Sigma
  type(kb_contour_gf),intent(inout)   :: G
  type(kb_contour_dgf),intent(inout)  :: dG
  type(kb_contour_dgf)                :: dG_new
  type(kb_contour_params),intent(in)  :: params
  complex(8),dimension(size(H0))      :: H
  type(kb_contour_gf)                 :: K
  integer                             :: N,L
  real(8)                             :: dt,dtau
  integer                             :: i,j,itau,jtau,s
  complex(8),dimension(:),allocatable :: KxG
  complex(8)                          :: dG_less
  !
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  call allocate_kb_contour_gf(K,params)
  K = Sigma%reg
  H = H0 - Sigma%hf(:)
  !
  call vide_kb_contour_gf_K(H,K,G,dG,dG_new,params)
  !
  call deallocate_kb_contour_gf(K)
end subroutine vide_kb_contour_gf_Sigma
