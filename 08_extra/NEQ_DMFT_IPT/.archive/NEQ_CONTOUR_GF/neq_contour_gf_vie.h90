!----------------------------------------------------------------------------
!  This subroutine solves a Volterra integral equation of the second kind,
!              G(t,t')+(K*G)(t,t')=Q(t,t')
!  for t=n*dt or t'=n*dt, using 2^nd *implicit* Runge-Kutta method.
!----------------------------------------------------------------------------
subroutine vie_kb_contour_gf_main(params,G,K,Q)
  type(kb_contour_params),intent(in)      :: params
  type(kb_contour_gf),intent(inout)       :: G
  type(kb_contour_gf),intent(in)          :: K
  type(kb_contour_gf),intent(in)          :: Q
  integer                                 :: N,L
  real(8)                                 :: dt,dtau
  integer                                 :: i,j,s,itau,jtau
  complex(8),dimension(:),allocatable     :: KxG
  !
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  call check_kb_contour_gf(params,G," vie_kb_contour_gf_main")
  call check_kb_contour_gf(params,K," vie_kb_contour_gf_main") 
  call check_kb_contour_gf(params,Q," vie_kb_contour_gf_main")
  !
  allocate(KxG(0:max(N,L)))
  !
  ! Ret component
  ! G^R(t,t') - \int_{t'}^t K^R(t,s)*G^R(s,t')ds = Q^R(t,t')
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  G%ret(N,N)=Q%ret(N,N)
  do j=1,N-1
     G%ret(N,j)=Q%ret(N,j)
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
     !Adds up to Q (r.h.s.)
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
  ! G^<(t_{N},t_{j}), j=1,...,N-1
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  do j=1, N-1
     G%less(N,j)=Q%less(N,j)
     !
     KxG(0:)=zero
     do s=0,L
        KxG(s)=K%lmix(N,s)*conjg(G%lmix(j,L-s))
     enddo
     G%less(N,j)=G%less(N,j)-xi*dtau*kb_trapz(KxG(0:),0,L)
     !
     KxG(0:)=zero
     do s=1,j
        KxG(s)=K%less(N,s)*conjg(G%ret(j,s))
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
  do i=1,N-1
     G%less(i,N) = -conjg(G%less(N,i))
  end do
  !
  ! G^{<}(t_{n},t_{n}) <= Diagonal
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  G%less(N,N)=Q%less(N,N)
  !
  KxG(0:)=zero
  do s=0,L
     KxG(s)=K%lmix(N,s)*conjg(G%lmix(N,L-s))
  end do
  G%less(N,N)=G%less(N,N)-xi*dtau*kb_trapz(KxG(0:),0,L)
  !
  KxG(0:)=zero
  do s=1,N
     KxG(s)=K%less(N,s)*conjg(G%ret(N,s))
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
  !
end subroutine vie_kb_contour_gf_main


subroutine vie_kb_contour_gf_Nso(params,G,K,Q)
  type(kb_contour_params),intent(in)     :: params
  type(kb_contour_gf),intent(inout)      :: G(:,:,:,:) ![Nspin][Nspin][Norb][Norb]
  type(kb_contour_gf),intent(in)         :: K(size(G,1),size(G,2),size(G,3),size(G,4))
  type(kb_contour_gf),intent(in)         :: Q(size(G,1),size(G,2),size(G,3),size(G,4))
  integer                                :: N,L
  integer                                :: Nspin,Norb,Nso
  real(8)                                :: dt,dtau
  integer                                :: i,j,s,itau,jtau,i1,i2,ii
  integer                                :: ispin,jspin,iorb,jorb
  integer                                :: kspin,korb
  complex(8),dimension(:),allocatable    :: KxG
  complex(8),dimension(:,:),allocatable  :: Amat,Bmat,Gmat
  !
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  call check_kb_contour_gf(params,G," vie_kb_contour_gf_Nso")
  call check_kb_contour_gf(params,K," vie_kb_contour_gf_Nso") 
  call check_kb_contour_gf(params,Q," vie_kb_contour_gf_Nso")
  !
  allocate(KxG(0:max(N,L)))
  !
  allocate(Amat(Nso,Nso));Amat=zero
  allocate(Bmat(Nso,Nso));Bmat=zero
  allocate(Gmat(Nso,Nso));Gmat=zero
  !
  ! Ret component
  ! G^R(t,t') - \int_{t'}^t K^R(t,s)*G^R(s,t')ds = Q^R(t,t')
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !TIP   t_{N}, t`_{N}
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              G(ispin,jspin,iorb,jorb)%ret(N,N)=Q(ispin,jspin,iorb,jorb)%ret(N,N)
           enddo
        enddo
     enddo
  enddo
  !
  !VERTICAL INTERVAL  t_{N}, t`_{j, j=1,...,N-1}
  Amat=zero
  Bmat=zero
  Gmat=zero

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
                 Bmat(i1,i2)=Q(ispin,jspin,iorb,jorb)%ret(N,j)
                 !G%ret(N,j)=Q%ret(N,j)
                 !
                 KxG(0:)=zero
                 do s=j,N-1
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s) + K(ispin,kspin,iorb,korb)%ret(N,s)*G(kspin,jspin,korb,jorb)%ret(s,j)
                       enddo
                    enddo
                 enddo
                 Bmat(i1,i2) = Bmat(i1,i2) + dt*kb_half_trapz(KxG(0:),j,N-1)
                 ! do s=j,N-1
                 !    KxG(s)=K%ret(N,s)*G%ret(s,j)
                 ! end do
                 !G%ret(N,j)=G%ret(N,j) + dt*kb_half_trapz(KxG(0:),j,N-1)
                 !
                 Amat(i1,i2) = Amat(i1,i2) - 0.5d0*dt*K(ispin,jspin,iorb,jorb)%ret(N,N)
                 !G%ret(N,j)=G%ret(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
              enddo
           enddo
        enddo
     enddo
     !
     call inv(Amat)
     Gmat = matmul(Amat,Bmat)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 i1 = iso_indx(ispin,iorb)
                 i2 = iso_indx(jspin,jorb)
                 G(ispin,jspin,iorb,jorb)%ret(N,j)=Gmat(i1,i2)
              enddo
           enddo
        enddo
     enddo
     !
  end do jtime_loop1


  !
  ! Lmix component
  ! G^\lmix(t,tau') - \int_0^t K^R(t,s)*G^\lmix(s,tau')ds
  !    = Q^\lmix(t,tau') + \int_0^\beta K^\lmix(t,s)*G^M(s,tau')ds
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Amat=zero
  Bmat=zero
  Gmat=zero
  jtau_loop: do jtau=0,L
     !
     Amat=zeye(Nso)
     !
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 i1 = iso_indx(ispin,iorb)
                 i2 = iso_indx(jspin,jorb)
                 !
                 Bmat(i1,i2)=Q(ispin,jspin,iorb,jorb)%lmix(N,jtau)                 
                 !G%lmix(N,jtau)=Q%lmix(N,jtau)
                 !
                 !Adds up to Q (r.h.s.)
                 KxG(0:)=zero           
                 do s=0,jtau
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%lmix(N,s)*G(kspin,jspin,korb,jorb)%mats(s+L-jtau)
                       enddo
                    enddo
                 enddo
                 Bmat(i1,i2)=Bmat(i1,i2)-dtau*kb_trapz(KxG(0:),0,jtau)
                 ! do s=0,jtau
                 !    KxG(s)=K%lmix(N,s)*G%mats(s+L-jtau)
                 ! end do
                 ! G%lmix(N,jtau)=G%lmix(N,jtau)-dtau*kb_trapz(KxG(0:),0,jtau)
                 !
                 KxG(0:)=zero
                 do s=jtau,L
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%lmix(N,s)*G(kspin,jspin,korb,jorb)%mats(s-jtau)
                       enddo
                    enddo
                 enddo
                 Bmat(i1,i2)=Bmat(i1,i2)+dtau*kb_trapz(KxG(0:),jtau,L)
                 ! do s=jtau,L
                 !    KxG(s)=K%lmix(N,s)*G%mats(s-jtau)
                 ! end do
                 ! G%lmix(N,jtau)=G%lmix(N,jtau)+dtau*kb_trapz(KxG(0:),jtau,L)
                 !
                 !
                 !Adds up to K*G^lmix (l.h.s.)
                 KxG(0:)=zero
                 do s=1,N-1
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%ret(N,s)*G(kspin,jspin,korb,jorb)%lmix(s,jtau)
                       enddo
                    enddo
                 enddo
                 Bmat(i1,i2) = Bmat(i1,i2) + dt*kb_half_trapz(KxG(0:),1,N-1)
                 ! do s=1,N-1
                 !    KxG(s)=K%ret(N,s)*G%lmix(s,jtau)
                 ! end do
                 ! G%lmix(N,jtau)=G%lmix(N,jtau) + dt*kb_half_trapz(KxG(0:),1,N-1)                 
                 !
                 Amat(i1,i2) = Amat(i1,i2) - 0.5d0*dt*K(ispin,jspin,iorb,jorb)%ret(N,N)
                 ! G%lmix(N,jtau)=G%lmix(N,jtau)/(1.d0-0.5d0*dt*K%ret(N,N))
              enddo
           enddo
        enddo
     enddo
     !
     call inv(Amat)
     Gmat = matmul(Amat,Bmat)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 i1 = iso_indx(ispin,iorb)
                 i2 = iso_indx(jspin,jorb)
                 G(ispin,jspin,iorb,jorb)%lmix(N,jtau)=Gmat(i1,i2)
              enddo
           enddo
        enddo
     enddo
     !
  end do jtau_loop

  !
  ! Less component
  ! G^<(t,t') - \int_0^t K^{R}(t,s)*G^{<}(s,t')ds
  !    = Q^<(t,t') - i\int_0^\beta K^\lmix(t,s)*G^\rmix(s,t')ds
  !      + \int_0^{t'} K^<(t,s)*G^A(s,t')ds
  ! G^<(t_{N},t_{j}), j=1,...,N-1
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  Amat=zero
  Bmat=zero
  Gmat=zero
  jtime_loop2: do j=1,N-1
     !
     Amat=zeye(Nso)
     !
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 i1 = iso_indx(ispin,iorb)
                 i2 = iso_indx(jspin,jorb)
                 !
                 Bmat(i1,i2)=Q(ispin,jspin,iorb,jorb)%less(N,j)
                 !G%less(N,j)=Q%less(N,j)
                 !
                 KxG(0:)=zero
                 do s=0,L
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%lmix(N,s)*conjg( G(kspin,jspin,korb,jorb)%lmix(j,L-s) ) !rmix <-- lmix
                       enddo
                    enddo
                 enddo
                 Bmat(i1,i2) = Bmat(i1,i2) - xi*dtau*kb_trapz(KxG(0:),0,L)
                 ! do s=0,L
                 !    KxG(s)=K%lmix(N,s)*conjg(G%lmix(j,L-s))
                 ! enddo
                 ! G%less(N,j)=G%less(N,j)-xi*dtau*kb_trapz(KxG(0:),0,L)
                 !
                 KxG(0:)=zero
                 do s=1,j
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%less(N,s)*conjg( G(kspin,jspin,korb,jorb)%ret(j,s) ) !adv <-- ret
                       enddo
                    enddo
                 enddo
                 Bmat(i1,i2) = Bmat(i1,i2) + dt*kb_trapz(KxG(0:),1,j)
                 ! do s=1,j
                 !    KxG(s)=K%less(N,s)*conjg(G%ret(j,s))
                 ! enddo
                 ! G%less(N,j)=G%less(N,j)+dt*kb_trapz(KxG(0:),1,j)
                 !
                 !Add to K*G^less
                 KxG(0:)=zero
                 do s=1,N-1
                    do kspin=1,Nspin
                       do korb=1,Norb
                          KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%ret(N,s)*G(kspin,jspin,korb,jorb)%less(s,j)
                       enddo
                    enddo
                 enddo
                 Bmat(i1,i2) = Bmat(i1,i2) + dt*kb_half_trapz(KxG(0:),1,N-1)
                 ! do s=1,N-1
                 !    KxG(s)=K%ret(N,s)*G%less(s,j)
                 ! enddo
                 ! G%less(N,j)=G%less(N,j) + dt*kb_half_trapz(KxG(0:),1,N-1)
                 !
                 Amat(i1,i2) = Amat(i1,i2) - 0.5d0*dt*K(ispin,jspin,iorb,jorb)%ret(N,N)
                 ! G%less(N,j)=G%less(N,j)/(1.d0-0.5d0*dt*K%ret(N,N))
              enddo
           enddo
        enddo
     enddo
     !
     !
     call inv(Amat)
     Gmat = matmul(Amat,Bmat)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 i1 = iso_indx(ispin,iorb)
                 i2 = iso_indx(jspin,jorb)
                 G(ispin,jspin,iorb,jorb)%less(N,j)=Gmat(i1,i2)
              enddo
           enddo
        enddo
     enddo
     !
  enddo jtime_loop2

  !
  ! G^<(t_{i},t_{N}) <= Hermite conjugate
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              do i=1,N-1
                 G(ispin,jspin,iorb,jorb)%less(i,N) = -conjg(G(ispin,jspin,iorb,jorb)%less(N,i))
              end do
           enddo
        enddo
     enddo
  enddo
  ! do i=1,N-1
  !    G%less(i,N) = -conjg(G%less(N,i))
  ! end do
  !
  !
  !
  ! G^{<}(t_{n},t_{n}) <= Diagonal
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Amat=zero
  Bmat=zero
  Gmat=zero
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              i1 = iso_indx(ispin,iorb)
              i2 = iso_indx(jspin,jorb)
              !
              Bmat(i1,i2)=Q(ispin,jspin,iorb,jorb)%less(N,N)
              !G%less(N,N)=Q%less(N,N)
              !
              !Add up to Q
              KxG(0:)=zero
              do s=0,L
                 do kspin=1,Nspin
                    do korb=1,Norb
                       KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%lmix(N,s)*conjg( G(kspin,jspin,korb,jorb)%lmix(N,L-s) ) !rmix <-- lmix
                    enddo
                 enddo
              enddo
              Bmat(i1,i2) = Bmat(i1,i2) - xi*dtau*kb_trapz(KxG(0:),0,L)
              ! do s=0,L
              !    KxG(s)=K%lmix(N,s)*conjg(G%lmix(N,L-s))
              ! end do
              ! G%less(N,N)=G%less(N,N)-xi*dtau*kb_trapz(KxG(0:),0,L)
              !
              KxG(0:)=zero
              do s=1,N
                 do kspin=1,Nspin
                    do korb=1,Norb
                       KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%less(N,s)*conjg( G(kspin,jspin,korb,jorb)%ret(N,s) ) !adv <-- ret
                    enddo
                 enddo
              enddo
              Bmat(i1,i2) = Bmat(i1,i2) + dt*kb_trapz(KxG(0:),1,N)
              ! do s=1,N
              !    KxG(s)=K%less(N,s)*conjg(G%ret(N,s))
              ! end do
              ! G%less(N,N)=G%less(N,N)+dt*kb_trapz(KxG(0:),1,N)
              !
              !Add up to K*G^<
              KxG(0:)=zero
              do s=1,N-1
                 do kspin=1,Nspin
                    do korb=1,Norb
                       KxG(s)=KxG(s)+K(ispin,kspin,iorb,korb)%ret(N,s)*G(kspin,jspin,korb,jorb)%less(s,N)
                    enddo
                 enddo
              enddo
              Bmat(i1,i2) = Bmat(i1,i2) + dt*kb_half_trapz(KxG(0:),1,N-1)
              ! do s=1,N-1
              !    KxG(s)=K%ret(N,s)*G%less(s,N)
              ! end do
              ! G%less(N,N)=G%less(N,N) + dt*kb_half_trapz(KxG(0:),1,N-1)
              !
              Amat(i1,i2) = Amat(i1,i2) - 0.5d0*dt*K(ispin,jspin,iorb,jorb)%ret(N,N)
              ! G%less(N,N)=G%less(N,N)/(1.d0-0.5d0*dt*K%ret(N,N))
           enddo
        enddo
     enddo
  enddo
  !
  call inv(Amat)
  Gmat = matmul(Amat,Bmat)
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              i1 = iso_indx(ispin,iorb)
              i2 = iso_indx(jspin,jorb)
              G(ispin,jspin,iorb,jorb)%less(N,N)=Gmat(i1,i2)
           enddo
        enddo
     enddo
  enddo
  !
  deallocate(KxG)
  deallocate(Amat,Bmat,Gmat)
  !
contains
  !
  function iso_indx(ispin,iorb) result(iso)
    integer :: ispin,iorb
    integer :: iso
    iso = iorb + (ispin-1)*Norb
  end function iso_indx
  !
end subroutine vie_kb_contour_gf_Nso


