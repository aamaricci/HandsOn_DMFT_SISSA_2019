!----------------------------------------------------------------------------
!  This subroutine solves a Volterra integro-differential equation of 
!  the second kind,
!              [i*d/dt-h(t)]G(t,t') = delta(t,t') + (K*G)(t,t'),
!  for t=n*dt or t'=n*dt, using 2^nd implicit Runge-Kutta method.
!
! TYPE: X=RET,LESS,LMIX (R,<,\LMIX)
!i d/dt G^x(t,:)  = h(t)G^x(t,:) + Q^x(t,:) + int_{t'}^{t}K^R(t,s)G^x(s,:)ds
!
!Q^x contains further integrals of the form K^a*G*b as obtained from Langreth
!rule's decomposition of the convolution.
!----------------------------------------------------------------------------
subroutine vide_kb_contour_gf_main(params,H,K,G,dG,dG_new)
  type(kb_contour_params),intent(in)  :: params
  complex(8),dimension(:)             :: H ![Ntime]
  type(kb_contour_gf),intent(in)      :: K
  type(kb_contour_gf),intent(inout)   :: G
  type(kb_contour_dgf),intent(inout)  :: dG
  type(kb_contour_dgf)                :: dG_new
  integer                             :: N,L
  real(8)                             :: dt,dtau
  integer                             :: i,j,itau,jtau,s
  complex(8),dimension(:),allocatable :: KxG
  complex(8)                          :: dG_less,A
  !
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  if(N==1) stop "ERROR vide_kb_contour_gf_main: called with N=1"
  !
  if(size(H)/=params%Ntime) stop "ERROR vide_kb_contour_gf_main: size(H) != Ntime"
  call check_kb_contour_gf(params,K,"vide_kb_contour_gf_main")
  call check_kb_contour_gf(params,G,"vide_kb_contour_gf_main")
  call check_kb_contour_gf(params,dG,"vide_kb_contour_gf_main")
  call check_kb_contour_gf(params,dG_new,"vide_kb_contour_gf_main")
  !
  allocate(KxG(0:max(N,L)))
  !
  !TYPE: X=RET,LESS,LMIX (R,<,\LMIX)
  !I D/DT G^X(T,:)  = H(T)G^X(T,:) + Q^X(T,:) + INT_{0,T'}^{T}K^R(T,S)G^X(S,:)DS
  !
  !Ret component
  ! d/dt G^R(t,:) = -i*h(t)*G^R(t,:) -i*delta(t,:) - i\int_{:}^t K^R(t,s)*G^R(s,:)ds
  G%ret(N,N)=-xi
  dG_new%ret(N)=-xi*H(N)*G%ret(N,N)
  do j=1,N-1
     G%ret(N,j)=G%ret(N-1,j) + 0.5d0*dt*dG%ret(j)
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
  do jtau=0,L
     G%lmix(N,jtau)=G%lmix(N-1,jtau)+0.5d0*dt*dG%lmix(jtau)
     do s=0,jtau
        KxG(s)=K%lmix(N,s)*G%mats(s+L-jtau)
     end do
     dG_new%lmix(jtau)=xi*dtau*kb_trapz(KxG(0:),0,jtau)
     do s=jtau,L
        KxG(s)=K%lmix(N,s)*G%mats(s-jtau)
     end do
     dG_new%lmix(jtau)=dG_new%lmix(jtau)-xi*dtau*kb_trapz(KxG(0:),jtau,L)!<= add -iQ(t)
     !
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
  !
  ! G^<(t_{N},t_{j}), d/dt G^<(t_{N},t_{j}) <== lower-right triangle
  do j=1,N-1
     G%less(N,j)=G%less(N-1,j) + 0.5d0*dt*dG%less(j)
     do s=0,L
        KxG(s)=K%lmix(N,s)*conjg(G%lmix(j,L-s))
     end do
     dG_new%less(j)=-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
     do s=1,j
        KxG(s)=K%less(N,s)*conjg(G%ret(j,s))
     end do
     dG_new%less(j)=dG_new%less(j)-xi*dt*kb_trapz(KxG(0:),1,j)!<= -iQ(t)
     !
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
  do i=1,N-1
     G%less(i,N)=-conjg(G%less(N,i))
  end do
  !
  ! d/dt G^<(t_{N-1},t_{N})
  dG_less=-xi*H(N-1)*G%less(N-1,N)
  do s=0,L
     KxG(s)=K%lmix(N-1,s)*conjg(G%lmix(N,L-s))
  end do
  dG_less=dG_less-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
  do s=1,N
     KxG(s)=K%less(N-1,s)*conjg(G%ret(N,s))
  end do
  dG_less=dG_less-xi*dt*kb_trapz(KxG(0:),1,N)
  do s=1,N-1
     KxG(s)=K%ret(N-1,s)*G%less(s,N)
  end do
  dG_less=dG_less-xi*dt*kb_trapz(KxG(0:),1,N-1)
  !
  !G^<(N,N), d/dt G^<(N,N)
  !d/dt G <== d/dt G_new
  G%less(N,N)=G%less(N-1,N)+0.5d0*dt*dG_less
  do s=0,L
     KxG(s)=K%lmix(N,s)*conjg(G%lmix(N,L-s))
  end do
  dG_new%less(N)=-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
  !
  do s=1,N
     KxG(s)=K%less(N,s)*conjg(G%ret(N,s))
  end do
  dG_new%less(N)=dG_new%less(N)-xi*dt*kb_trapz(KxG(0:),1,N)
  !
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
  !
end subroutine vide_kb_contour_gf_main



subroutine vide_kb_contour_gf_Nso(params,H,K,G,dG,dG_new)
  type(kb_contour_params),intent(in)  :: params
  complex(8),dimension(:,:,:)             :: H           ![Nspin*Norb][Nspin*Norb][Nt]
  type(kb_contour_gf),intent(in)          :: K(:,:,:,:)  ![Nspin][Nspin][Norb][Norb]
  type(kb_contour_gf),intent(inout)       :: G(:,:,:,:)     !as K
  type(kb_contour_dgf),intent(inout)      :: dG(:,:,:,:)    !as K
  type(kb_contour_dgf)                    :: dG_new(:,:,:,:)!as K

  integer                                 :: N,L,Nspin,Norb,Nso
  real(8)                                 :: dt,dtau
  integer                                 :: i,j,itau,jtau,s,i1,i2,ii,ik
  integer                                 :: ispin,jspin,iorb,jorb
  integer                                 :: kspin,korb
  complex(8),dimension(:),allocatable     :: KxG
  complex(8)                              :: dG_less
  complex(8),dimension(:,:),allocatable   :: Amat,Bmat,Gmat,dGmat
  !
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  if(N==1) stop "ERROR vide_kb_contour_gf_Nso: called with N=1"
  !
  if(size(H)/=params%Ntime) stop "ERROR vide_kb_contour_gf_main: size(H) != Ntime"
  call check_kb_contour_gf(params,K,"vide_kb_contour_gf_main")
  call check_kb_contour_gf(params,G,"vide_kb_contour_gf_main")
  call check_kb_contour_gf(params,dG,"vide_kb_contour_gf_main")
  call check_kb_contour_gf(params,dG_new,"vide_kb_contour_gf_main")
  !
  allocate(KxG(0:max(N,L)))
  !
  allocate(Amat(Nso,Nso));Amat=zero
  allocate(Bmat(Nso,Nso));Bmat=zero
  allocate(Gmat(Nso,Nso));Gmat=zero
  allocate(dGmat(Nso,Nso));dGmat=zero




  !TYPE: X=RET,LESS,LMIX (R,<,\LMIX)
  !I D/DT G^X(T,:)  = H(T)G^X(T,:) + Q^X(T,:) + INT_{0,T'}^{T}K^R(T,S)G^X(S,:)DS
  !
  !Ret component
  ! d/dt G^R(t,:) = -i*h(t)*G^R(t,:) -i*delta(t,:) - i\int_{:}^t K^R(t,s)*G^R(s,:)ds
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !TIP t_{N}, t`_{N}
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              i1 = iso_indx(ispin,iorb)
              i2 = iso_indx(jspin,jorb)
              G(ispin,jspin,iorb,jorb)%ret(N,N)=-xi
              do kspin=1,Nspin
                 do korb=1,Norb
                    ik = iso_indx(kspin,korb)
                    dG_new(ispin,jspin,iorb,jorb)%ret(N)=dG_new(ispin,jspin,iorb,jorb)%ret(N)-&
                         xi*H(i1,ik,N)*G(kspin,jspin,korb,jorb)%ret(N,N)
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  ! G%ret(N,N)=-xi
  ! dG_new%ret(N)=-xi*H(N)*G%ret(N,N)
  !
  !VERTICAL INTERVAL t_{N}, t`_{j, j=1,...,N-1}
  Amat=zero
  Bmat=zero
  Gmat=zero
  dGmat=zero
  jtime_loop1: do j=1,N-1
     !
     Amat = zeye(Nso)
     !
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 i1 = iso_indx(ispin,iorb)
                 i2 = iso_indx(jspin,jorb)
                 !
                 !Add G^R(t_{N-1},t_j) + dt/2* d_tG^ret(t_{N-1},j)
                 Bmat(i1,i2) = G(ispin,jspin,iorb,jorb)%ret(N-1,j) + 0.5d0*dt*dG(ispin,jspin,iorb,jorb)%ret(j)
                 ! G%ret(N,j)=G%ret(N-1,j) + 0.5d0*dt*dG%ret(j)


                 !Add -xi*K^R(t_N,s)*G^R(s,t_j)
                 KxG(0:)=zero
                 do s=j,N-1
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%ret(N,s)*G(kspin,jspin,korb,jorb)%ret(s,j)
                       enddo
                    enddo
                 enddo
                 dGmat(i1,i2) = -xi*dt*kb_half_trapz(KxG(0:),j,N-1)
                 Bmat(i1,i2) = Bmat(i1,i2) + 0.5d0*dt*dGmat(i1,i2)
                 ! do s=j,N-1
                 !    KxG(s)=K%ret(N,s)*G%ret(s,j)
                 ! enddo
                 ! dG_new%ret(j)=-xi*dt*kb_half_trapz(KxG(0:),j,N-1)
                 ! G%ret(N,j)=G%ret(N,j) + 0.5d0*dt*dG_new%ret(j)

                 Amat(i1,i2) = Amat(i1,i2) + 0.5d0*xi*dt*H(i1,i2,N) + 0.25d0*xi*dt**2*K(ispin,jspin,iorb,jorb)%ret(N,N)
                 ! G%ret(N,j)=G%ret(N,j)/(1.d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
              enddo
           enddo
        enddo
     enddo
     !
     call inv(Amat)
     Gmat = matmul(Amat,Bmat)
     !Update derivative d_t G^R(t,:) as: (dGmat already holds all the terms not included below)
     !d/dt G^R(t,:) = -i*h(t)*G^R(t,:) -i*delta(t,:) - i\int_{:}^t K^R(t,s)*G^R(s,:)ds
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 i1 = iso_indx(ispin,iorb)
                 i2 = iso_indx(jspin,jorb)
                 !
                 do kspin=1,Nspin
                    do korb=1,Norb
                       dGmat(i1,i2) = dGmat(i1,i2) - 0.5d0*xi*dt*K(ispin,kspin,iorb,korb)%ret(N,N)*G(kspin,jspin,korb,jorb)%ret(N,j)
                    enddo
                 enddo
                 !
                 do kspin=1,Nspin
                    do korb=1,Norb
                       ik = iso_indx(kspin,korb)
                       dGmat(i1,i2) = dGmat(i1,i2) - xi*H(i1,ik,N)*G(kspin,jspin,korb,jorb)%ret(N,j)
                    enddo
                 enddo
                 !
                 dG_new(ispin,jspin,iorb,jorb)%ret(j) = dGmat(i1,i2) 
                 G(ispin,jspin,iorb,jorb)%ret(N,j)=Gmat(i1,i2)
                 !
              enddo
           enddo
        enddo
     enddo
     ! dG_new%ret(j)=dG_new%ret(j) - xi*H(N)*G%ret(N,j) - 0.5d0*xi*dt*K%ret(N,N)*G%ret(N,j)
     !
  enddo jtime_loop1


  !
  !Lmix component
  !d/dt G^\lmix(t,:) = -i*H(t)*G^\lmix(t,:) -i\int_0^\beta K^\lmix(t,s)*G^M(s,:)ds -i\int_0^t K^R(t,s)G^\lmix(s,:)ds
  Amat=zero
  Bmat=zero
  Gmat=zero
  dGmat=zero
  jtau_loop: do jtau=0,L
     !
     Amat = zeye(Nso)
     !
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 i1 = iso_indx(ispin,iorb)
                 i2 = iso_indx(jspin,jorb)
                 !
                 !Add G^lmix(t_{N-1},tau_j) + dt/2* d_tG^lmix(t_{N-1},tau_j)
                 Bmat(i1,i2) = G(ispin,jspin,iorb,jorb)%lmix(N-1,jtau)+0.5d0*dt*dG(ispin,jspin,iorb,jorb)%lmix(jtau)
                 !G%lmix(N,jtau)=G%lmix(N-1,jtau)+0.5d0*dt*dG%lmix(jtau)
                 !
                 !Add Q^lmix(t_{N},tau_j) = K^lmix(N,stau)*G^mats(stau-jtau) stau<jtau
                 KxG(0:)=zero
                 do s=0,jtau
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%lmix(N,s)*G(kspin,jspin,korb,jorb)%mats(s+L-jtau)
                       enddo
                    enddo
                 enddo
                 dGmat(i1,i2) = xi*dtau*kb_trapz(KxG(0:),0,jtau)                 
                 ! do s=0,jtau
                 !    KxG(s)=K%lmix(N,s)*G%mats(s+L-jtau)
                 ! end do
                 ! dG_new%lmix(jtau)=xi*dtau*kb_trapz(KxG(0:),0,jtau)
                 !
                 !
                 !Add Q^lmix(t_{N},tau_j) = K^lmix(N,stau)*G^mats(stau-jtau) stau>jtau
                 KxG(0:)=zero
                 do s=jtau,L
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%lmix(N,s)*G(kspin,jspin,korb,jorb)%mats(s-jtau)
                       enddo
                    enddo
                 enddo
                 dGmat(i1,i2)=dGmat(i1,i2)-xi*dtau*kb_trapz(KxG(0:),jtau,L)!<= add -iQ(t)
                 ! do s=jtau,L
                 !    KxG(s)=K%lmix(N,s)*G%mats(s-jtau)
                 ! end do
                 ! dG_new%lmix(jtau)=dG_new%lmix(jtau)-xi*dtau*kb_trapz(KxG(0:),jtau,L)!<= add -iQ(t)
                 !
                 !
                 !Add -xi*K^R(t_N,s)*G^lmix(s,tau_j)
                 KxG(0:)=zero
                 do s=1,N-1
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%ret(N,s)*G(kspin,jspin,korb,jorb)%lmix(s,jtau)
                       enddo
                    enddo
                 enddo
                 dGmat(i1,i2)=dGmat(i1,i2)-xi*dt*kb_half_trapz(KxG(0:),1,N-1)
                 ! do s=1,N-1
                 !    KxG(s)=K%ret(N,s)*G%lmix(s,jtau)
                 ! end do
                 ! dG_new%lmix(jtau)=dG_new%lmix(jtau)-xi*dt*kb_half_trapz(KxG(0:),1,N-1)
                 !
                 Bmat(i1,i2) = Bmat(i1,i2) + 0.5d0*dt*dGmat(i1,i2)
                 !G%lmix(N,jtau)=G%lmix(N,jtau) + 0.5d0*dt*dG_new%lmix(jtau)
                 !
                 Amat(i1,i2) = Amat(i1,i2) + 0.5d0*xi*dt*H(i1,i2,N) + 0.25d0*xi*dt**2*K(ispin,jspin,iorb,jorb)%ret(N,N)
                 ! G%lmix(N,jtau)=G%lmix(N,jtau)/(1.d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
                 !
              enddo
           enddo
        enddo
     enddo
     !
     call inv(Amat)
     Gmat = matmul(Amat,Bmat)
     !Update derivative d_t G^\lmix(t,:) as: (dGmat already holds all the terms not included below)
     !d/dt G^\lmix(t,:) = -i*H(t)*G^\lmix(t,:) -i\int_0^\beta K^\lmix(t,s)*G^M(s,:)ds -i\int_0^t K^R(t,s)G^\lmix(s,:)ds
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 i1 = iso_indx(ispin,iorb)
                 i2 = iso_indx(jspin,jorb)
                 !
                 do kspin=1,Nspin
                    do korb=1,Norb
                       dGmat(i1,i2) = dGmat(i1,i2) - &
                            0.5d0*xi*dt*K(ispin,kspin,iorb,korb)%ret(N,N)*G(kspin,jspin,korb,jorb)%lmix(N,jtau)
                    enddo
                 enddo
                 !
                 do kspin=1,Nspin
                    do korb=1,Norb
                       ik = iso_indx(kspin,korb)
                       dGmat(i1,i2) = dGmat(i1,i2) - xi*H(i1,ik,N)*G(kspin,jspin,korb,jorb)%lmix(N,jtau)
                    enddo
                 enddo
                 !
                 dG_new(ispin,jspin,iorb,jorb)%lmix(jtau) = dGmat(i1,i2) 
                 G(ispin,jspin,iorb,jorb)%lmix(N,jtau)=Gmat(i1,i2)
                 !
              enddo
           enddo
        enddo
     enddo
     !dG_new%lmix(jtau)=dG_new%lmix(jtau)-xi*H(N)*G%lmix(N,jtau)-0.5d0*xi*dt*K%ret(N,N)*G%lmix(N,jtau)
     !
  enddo jtau_loop
  !
  !
  !
  !Less component
  !d/dt G^<(t,:) = -i*H(t)*G^<(t,:) -i*[ (-i)*\int_0^\beta K^\lmix(t,s)*G^\rmix(s,:)ds + \int_0^{:}K^<(t,s)*G^A(s,:)ds ]
  !                                 -i*\int_0^t K^R(t,s)*G^<(s,:)ds
  !
  ! G^<(t_{N},t_{j}), d/dt G^<(t_{N},t_{j}) <== lower-right triangle
  Amat=zero
  Bmat=zero
  Gmat=zero
  dGmat=zero
  jtime_loop2: do j=1,N-1
     !
     Amat = zeye(Nso)
     !
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 i1 = iso_indx(ispin,iorb)
                 i2 = iso_indx(jspin,jorb)
                 !
                 Bmat(i1,i2) = G(ispin,jspin,iorb,jorb)%less(N-1,j) + 0.5d0*dt*dG(ispin,jspin,iorb,jorb)%less(j)
                 !G%less(N,j)=G%less(N-1,j) + 0.5d0*dt*dG%less(j)
                 !
                 KxG(0:)=zero
                 do s=0,L
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%lmix(N,s)*conjg( G(kspin,jspin,korb,jorb)%lmix(j,L-s) ) !rmix <-- lmix
                       enddo
                    enddo
                 enddo
                 dGmat(i1,i2)=-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
                 ! do s=0,L
                 !    KxG(s)=K%lmix(N,s)*conjg(G%lmix(j,L-s))
                 ! end do
                 ! dG_new%less(j)=-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
                 !
                 KxG(0:)=zero
                 do s=1,j
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%less(N,s)*conjg( G(kspin,jspin,korb,jorb)%ret(j,s) ) !adv <-- ret
                       enddo
                    enddo
                 enddo
                 dGmat(i1,i2)=dGmat(i1,i2)-xi*dt*kb_trapz(KxG(0:),1,j) !<= -iQ(t)
                 ! do s=1,j
                 !    KxG(s)=K%less(N,s)*conjg(G%ret(j,s))
                 ! end do
                 ! dG_new%less(j)=dG_new%less(j)-xi*dt*kb_trapz(KxG(0:),1,j)!<= -iQ(t)
                 !
                 KxG(0:)=zero
                 do s=1,N-1
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%ret(N,s)*G(kspin,jspin,korb,jorb)%less(s,j)
                       enddo
                    enddo
                 enddo
                 dGmat(i1,i2)=dGmat(i1,i2)-xi*dt*kb_half_trapz(KxG(0:),1,N-1)
                 ! do s=1,N-1
                 !    KxG(s)=K%ret(N,s)*G%less(s,j)
                 ! end do
                 ! dG_new%less(j)=dG_new%less(j)-xi*dt*kb_half_trapz(KxG(0:),1,N-1)
                 !
                 Bmat(i1,i2) = Bmat(i1,i2) + 0.5d0*dt*dGmat(i1,i2)
                 ! G%less(N,j)=G%less(N,j) + 0.5d0*dt*dG_new%less(j)
                 !
                 Amat(i1,i2) = Amat(i1,i2) + 0.5d0*xi*dt*H(i1,i2,N) + 0.25d0*xi*dt**2*K(ispin,jspin,iorb,jorb)%ret(N,N)
                 ! G%less(N,j)=G%less(N,j)/(1.d0 + 0.5d0*xi*dt*H(N) + 0.25d0*xi*dt**2*K%ret(N,N))
                 !

              enddo
           enddo
        enddo
     enddo
     !
     call inv(Amat)
     Gmat = matmul(Amat,Bmat)
     !Update derivative d_t G^<(t,:) as: (dGmat already holds all the terms not included below)
     !d/dt G^<(t,:) = -i*H(t)*G^<(t,:) -i*[ (-i)*\int_0^\beta K^\lmix(t,s)*G^\rmix(s,:)ds + \int_0^{:}K^<(t,s)*G^A(s,:)ds ]
     !                                 -i*\int_0^t K^R(t,s)*G^<(s,:)ds
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 i1 = iso_indx(ispin,iorb)
                 i2 = iso_indx(jspin,jorb)
                 !
                 do kspin=1,Nspin
                    do korb=1,Norb
                       dGmat(i1,i2) = dGmat(i1,i2) - 0.5d0*xi*dt*K(ispin,kspin,iorb,korb)%ret(N,N)*G(kspin,jspin,korb,jorb)%less(N,j)
                    enddo
                 enddo
                 !
                 do kspin=1,Nspin
                    do korb=1,Norb
                       ik = iso_indx(kspin,korb)
                       dGmat(i1,i2) = dGmat(i1,i2) - xi*H(i1,ik,N)*G(kspin,jspin,korb,jorb)%less(N,j)
                    enddo
                 enddo
                 !
                 dG_new(ispin,jspin,iorb,jorb)%less(j) = dGmat(i1,i2) 
                 G(ispin,jspin,iorb,jorb)%less(N,j)=Gmat(i1,i2)
                 !
              enddo
           enddo
        enddo
     enddo
     ! dG_new%less(j)=dG_new%less(j)-xi*H(N)*G%less(N,j)-0.5d0*xi*dt*K%ret(N,N)*G%less(N,j)
     !
  enddo jtime_loop2
  !
  !
  ! G^<(t_{i},t_{N}), d/dt G^<(t_{i},t_{N}) <== upper left triangle
  ! Hermitian conjugate G
  do i=1,N-1
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 G(ispin,jspin,iorb,jorb)%less(i,N)=-conjg(G(ispin,jspin,iorb,jorb)%less(N,i))
              enddo
           enddo
        enddo
     enddo
  enddo
  ! do i=1,N-1
  !    G%less(i,N)=-conjg(G%less(N,i))
  ! end do
  !
  !

  ! d/dt G^<(t_{N-1},t_{N})
  Amat=zero
  Bmat=zero
  Gmat=zero
  dGmat=zero
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              i1 = iso_indx(ispin,iorb)
              i2 = iso_indx(jspin,jorb)

              do kspin=1,Nspin
                 do korb=1,Norb
                    ik=iso_indx(kspin,korb)
                    dGmat(i1,i2) = dGmat(i1,i2)-xi*H(i1,ik,N-1)*G(kspin,jspin,korb,jorb)%less(N-1,N)
                 enddo
              enddo
              ! dG_less=-xi*H(N-1)*G%less(N-1,N)
              !
              KxG(0:)=zero
              do s=0,L
                 do kspin=1,Nspin
                    do korb=1,Norb
                       KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%lmix(N-1,s)*conjg( G(kspin,jspin,korb,jorb)%lmix(N,L-s) )!rmix <-- lmix
                    enddo
                 enddo
              enddo
              dGmat(i1,i2) = dGmat(i1,i2)-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
              ! do s=0,L
              !    KxG(s)=K%lmix(N-1,s)*conjg(G%lmix(N,L-s))
              ! end do
              ! dG_less=dG_less-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
              !
              KxG(0:)=zero
              do s=1,N
                 do kspin=1,Nspin
                    do korb=1,Norb
                       KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%less(N-1,s)*conjg( G(kspin,jspin,korb,jorb)%ret(N,s) ) !adv <-- ret
                    enddo
                 enddo
              enddo
              dGmat(i1,i2) = dGmat(i1,i2)-xi*dt*kb_trapz(KxG(0:),1,N)
              ! do s=1,N
              !    KxG(s)=K%less(N-1,s)*conjg(G%ret(N,s))
              ! end do
              ! dG_less=dG_less-xi*dt*kb_trapz(KxG(0:),1,N)
              !
              KxG(0:)=zero
              do s=1,N-1
                 do kspin=1,Nspin
                    do korb=1,Norb
                       KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%ret(N-1,s)*G(kspin,jspin,korb,jorb)%less(s,N)
                    enddo
                 enddo
              enddo
              ! do s=1,N-1
              !    KxG(s)=K%ret(N-1,s)*G%less(s,N)
              ! end do
              ! dG_less=dG_less-xi*dt*kb_trapz(KxG(0:),1,N-1)
              !
              !
              !
              !G^<(N,N), d/dt G^<(N,N)
              !
              Bmat(i1,i2) = G(ispin,jspin,iorb,jorb)%less(N-1,N)+0.5d0*dt*dGmat(i1,i2)
              !G%less(N,N)=G%less(N-1,N)+0.5d0*dt*dG_less
              !
              !
              KxG(0:)=zero
              do s=0,L
                 do kspin=1,Nspin
                    do korb=1,Norb
                       KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%lmix(N,s)*conjg( G(kspin,jspin,korb,jorb)%lmix(N,L-s) ) !get_rmix(G,s,N,L)
                    enddo
                 enddo
              enddo
              dGmat(i1,i2)=-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L) !this one reset dGmat
              ! do s=0,L
              !    KxG(s)=K%lmix(N,s)*conjg(G%lmix(N,L-s))
              ! end do
              ! dG_new%less(N)=-xi*(-xi)*dtau*kb_trapz(KxG(0:),0,L)
              !
              !
              KxG(0:)=zero
              do s=1,N
                 do kspin=1,Nspin
                    do korb=1,Norb
                       KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%less(N,s)*conjg( G(kspin,jspin,korb,jorb)%ret(N,s) ) !get_adv(G,s,N)
                    enddo
                 enddo
              enddo
              dGmat(i1,i2)=dGmat(i1,i2) - xi*dt*kb_trapz(KxG(0:),1,N)
              ! do s=1,N
              !    KxG(s)=K%less(N,s)*conjg(G%ret(N,s))
              ! end do
              ! dG_new%less(N)=dG_new%less(N)-xi*dt*kb_trapz(KxG(0:),1,N)
              !
              !
              KxG(0:)=zero
              do s=1,N-1
                 do kspin=1,Nspin
                    do korb=1,Norb
                       KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%ret(N,s)*G(kspin,jspin,korb,jorb)%less(s,N)
                    enddo
                 enddo
              enddo
              dGmat(i1,i2)=dGmat(i1,i2) - xi*dt*kb_half_trapz(KxG(0:),1,N-1)
              ! do s=1,N-1
              !    KxG(s)=K%ret(N,s)*G%less(s,N)
              ! end do
              ! dG_new%less(N)=dG_new%less(N)-xi*dt*kb_half_trapz(KxG(0:),1,N-1)
              !
              Bmat(i1,i2) = Bmat(i1,i2) + 0.5d0*dt*dGmat(i1,i2)
              ! G%less(N,N)=G%less(N,N)+0.5d0*dt*dG_new%less(N)
              !
              Amat(i1,i2) = Amat(i1,i2) + 0.5d0*xi*dt*H(i1,i2,N) + 0.25d0*xi*dt**2*K(ispin,jspin,iorb,jorb)%ret(N,N)
              ! G%less(N,N)=G%less(N,N)/(1.d0+0.5d0*xi*dt*H(N)+0.25d0*xi*dt**2*K%ret(N,N))
              !
           enddo
        enddo
     enddo
  enddo
  !
  !
  call inv(Amat)
  Gmat = matmul(Amat,Bmat)
  !Update derivative d_t G^<(t,:) as: (dGmat already holds all the terms not included below)
  !d/dt G^<(t,:) = -i*H(t)*G^<(t,:) -i*[ (-i)*\int_0^\beta K^\lmix(t,s)*G^\rmix(s,:)ds + \int_0^{:}K^<(t,s)*G^A(s,:)ds ]
  !                                 -i*\int_0^t K^R(t,s)*G^<(s,:)ds
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              i1 = iso_indx(ispin,iorb)
              i2 = iso_indx(jspin,jorb)
              !
              do kspin=1,Nspin
                 do korb=1,Norb
                    dGmat(i1,i2) = dGmat(i1,i2) - 0.5d0*xi*dt*K(ispin,kspin,iorb,korb)%ret(N,N)*G(kspin,jspin,korb,jorb)%less(N,N)
                 enddo
              enddo
              !
              do kspin=1,Nspin
                 do korb=1,Norb
                    ik = iso_indx(kspin,korb)
                    dGmat(i1,i2) = dGmat(i1,i2) - xi*H(i1,ik,N)*G(kspin,jspin,korb,jorb)%less(N,N)
                 enddo
              enddo
              !
              dG_new(ispin,jspin,iorb,jorb)%less(N) = dGmat(i1,i2) 
              G(ispin,jspin,iorb,jorb)%less(N,N)=Gmat(i1,i2)
              !
           enddo
        enddo
     enddo
  enddo
  ! dG_new%less(N)=dG_new%less(N)-xi*H(N)*G%less(N,N)-0.5d0*xi*dt*K%ret(N,N)*G%less(N,N)
  !
  !
  deallocate(KxG)
  deallocate(Amat,Bmat,Gmat,dGmat)
  !
contains
  !
  function iso_indx(ispin,iorb) result(iso)
    integer :: ispin,iorb
    integer :: iso
    iso = iorb + (ispin-1)*Norb
  end function iso_indx
  !
end subroutine vide_kb_contour_gf_Nso
