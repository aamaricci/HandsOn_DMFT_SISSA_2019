!C(t,t')=(A*B)(t,t'), with t=t_max && t'=0,t_max
subroutine convolute_kb_contour_gf_main(params,A,B,C,dcoeff,ccoeff)
  type(kb_contour_params)             :: params
  type(kb_contour_gf)                 :: A,B,C
  real(8),optional                    :: dcoeff
  complex(8),optional                 :: ccoeff
  integer                             :: N,L
  real(8)                             :: dt,dtau
  complex(8),dimension(:),allocatable :: AxB    
  integer                             :: i,j,k,itau,jtau
  complex(8),dimension(:,:),allocatable :: Amat,Bmat,Cmat
  integer                               :: Ntot
  !
  N   = params%Nt      !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  if(N==1) stop "ERROR convolute_kb_contour_gf_main: called with N=1"
  !
  call check_kb_contour_gf(params,A,"convolute_kb_contour_gf_main") 
  call check_kb_contour_gf(params,B,"convolute_kb_contour_gf_main")
  call check_kb_contour_gf(params,C,"convolute_kb_contour_gf_main")
  !
  allocate(AxB(0:max(L,N)))
  !
  !Ret. component
  !C^R(t,t')=\int_{t'}^t ds A^R(t,s)*B^R(s,t')
  C%ret(N,1:N)=zero
  do j=1,N           !for all t' in 0:t_max
     AxB(0:)  = zero !AxB should be set to zero before integration
     do k=j,N        !store the convolution between t'{=j} and t{=N}
        AxB(k) = A%ret(N,k)*B%ret(k,j)
     enddo
     C%ret(n,j) = C%ret(n,j) + dt*kb_trapz(AxB(0:),j,N)
  enddo
  !
  !Lmix. component
  !C^\lmix(t,tau')=\int_0^{beta} ds A^\lmix(t,s)*B^M(s,tau')
  !                  +\int_0^{t} ds A^R(t,s)*B^\lmix(s,tau')
  !               = I1 + I2
  C%lmix(N,0:L)=zero
  do jtau=0,L
     !I1:
     AxB(0:) = zero
     !break the integral I1 in two parts to take care of the 
     !sign of (tau-tau').
     do k=0,jtau
        AxB(k)=A%lmix(N,k)*B%mats(L+k-jtau)
     end do
     C%lmix(N,jtau)=C%lmix(N,jtau)-dtau*kb_trapz(AxB(0:),0,jtau)
     do k=jtau,L
        AxB(k)=A%lmix(N,k)*B%mats(k-jtau)
     end do
     C%lmix(n,jtau)=C%lmix(n,jtau)+dtau*kb_trapz(AxB(0:),jtau,L)
     !
     !I2:
     AxB(0:) = zero
     do k=1,N
        AxB(k) = A%ret(N,k)*B%lmix(k,jtau)
     enddo
     C%lmix(N,jtau) = C%lmix(N,jtau) + dt*kb_trapz(AxB(0:),1,N)
  enddo
  !
  !Less component
  !C^<(t,t')=-i\int_0^{beta} ds A^\lmix(t,s)*B^\rmix(s,t')
  !             +\int_0^{t'} ds A^<(t,s)*B^A(s,t')
  !             +\int_0^{t} ds A^R(t,s)*B^<(s,t')
  ! (t,t')=>(N,j) <==> Vertical side, no tip (j=1,N-1)
  do j=1,N-1
     C%less(N,j)=zero
     do k=0,L
        AxB(k)=A%lmix(N,k)*conjg(B%lmix(j,L-k))
     end do
     C%less(N,j)=C%less(N,j)-xi*dtau*kb_trapz(AxB(0:),0,L)
     !
     do k=1,j
        AxB(k)=A%less(N,k)*conjg(B%ret(j,k))
     end do
     C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,j)
     !
     do k=1,N
        AxB(k)=A%ret(N,k)*B%less(k,j)
     end do
     C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,N)
  end do
  !
  ! (t,t')=>(i,N) <==> Horizontal side, w/ tip (i=1,N)
  do i=1,N
     C%less(i,N)=zero
     do k=0,L
        AxB(k)=A%lmix(i,k)*conjg(B%lmix(n,L-k))
     end do
     C%less(i,N)=C%less(i,N)-xi*dtau*kb_trapz(AxB(0:),0,L)
     !
     do k=1,N
        AxB(k)=A%less(i,k)*conjg(B%ret(N,k))
     end do
     C%less(i,N)=C%less(i,N)+dt*kb_trapz(AxB(0:),1,N)
     !
     do k=1,i
        AxB(k)=A%ret(i,k)*B%less(k,N)
     end do
     C%less(i,N)=C%less(i,N)+dt*kb_trapz(AxB(0:),1,i)
  end do
  !    
  if(present(dcoeff))then
     C%lmix(N,0:L) = dcoeff*C%lmix(N,0:L)
     C%less(N,1:N-1)=dcoeff*C%less(N,1:N-1)
     C%less(1:N,N)=dcoeff*C%less(1:N,N)
     C%ret(N,1:N)=dcoeff*C%ret(N,1:N)
  endif
  if(present(ccoeff))then
     C%lmix(N,0:L) = ccoeff*C%lmix(N,0:L)
     C%less(N,1:N-1)=ccoeff*C%less(N,1:N-1)
     C%less(1:N,N)=ccoeff*C%less(1:N,N)
     C%ret(N,1:N)=ccoeff*C%ret(N,1:N)
  endif
  deallocate(AxB)
end subroutine convolute_kb_contour_gf_main


subroutine convolute_kb_contour_gf_Nso(params,A,B,C,dcoeff,ccoeff)
  type(kb_contour_params)             :: params
  type(kb_contour_gf)                 :: A(:,:,:,:) ![Nspin][Nspin][Norb][Norb]
  type(kb_contour_gf)                 :: B(size(A,1),size(A,2),size(A,3),size(A,4))
  type(kb_contour_gf)                 :: C(size(A,1),size(A,2),size(A,3),size(A,4))
  real(8),optional                    :: dcoeff
  complex(8),optional                 :: ccoeff
  integer                             :: N,L,Nspin,Norb,Nso
  real(8)                             :: dt,dtau
  complex(8),dimension(:),allocatable :: AxB
  integer                             :: i,j,s,itau,jtau,i1,i2,ii
  integer                             :: ispin,jspin,iorb,jorb
  integer                             :: kspin,korb
  !
  N   = params%Nt      !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dt  = params%dt
  dtau= params%dtau
  !
  if(N==1) stop "ERROR convolute_kb_contour_gf_gf_Nso: called with N=1"
  !
  Nspin = size(A,1)
  Norb  = size(A,3)
  Nso   = Nspin*Norb
  !
  call check_kb_contour_gf(params,A,"convolute_kb_contour_gf_Nso") 
  call check_kb_contour_gf(params,B,"convolute_kb_contour_gf_Nso")
  call check_kb_contour_gf(params,C,"convolute_kb_contour_gf_Nso")
  !
  !
  allocate(AxB(0:max(L,N)));AxB=zero
  !
  ispin_loop: do ispin=1,Nspin
     jspin_loop: do jspin=1,Nspin
        iorb_loop: do iorb=1,Norb
           jorb_loop: do jorb=1,Norb
              !
              !Ret. component
              !C^R(t,t')=\int_{t'}^t ds A^R(t,s)*B^R(s,t')
              C(ispin,jspin,iorb,jorb)%ret(N,1:N)=zero
              do j=1,N           !for all t' in 0:t_max
                 AxB(0:)  = zero !AxB should be set to zero before integration
                 do s=j,N        !store the convolution between t'{=j} and t{=N}
                    do kspin=1,Nspin
                       do korb=1,Norb
                          AxB(s) = AxB(s) + A(ispin,kspin,iorb,korb)%ret(N,s)*B(kspin,jspin,korb,jorb)%ret(s,j)
                       enddo
                    enddo
                 enddo
                 C(ispin,jspin,iorb,jorb)%ret(n,j) = C(ispin,jspin,iorb,jorb)%ret(n,j) + dt*kb_trapz(AxB(0:),j,N)
              enddo
              !
              !Lmix. component
              !C^\lmix(t,tau')=\int_0^{beta} ds A^\lmix(t,s)*B^M(s,tau')
              !                  +\int_0^{t} ds A^R(t,s)*B^\lmix(s,tau')
              !               = I1 + I2
              C(ispin,jspin,iorb,jorb)%lmix(N,0:L)=zero
              do jtau=0,L
                 !I1:
                 AxB(0:) = zero
                 do s=0,jtau
                    do kspin=1,Nspin
                       do korb=1,Norb
                          AxB(s) = AxB(s) + A(ispin,kspin,iorb,korb)%lmix(N,s)*B(kspin,jspin,korb,jorb)%mats(L+s-jtau)
                       enddo
                    enddo
                 enddo
                 C(ispin,jspin,iorb,jorb)%lmix(N,jtau)=C(ispin,jspin,iorb,jorb)%lmix(N,jtau)-dtau*kb_trapz(AxB(0:),0,jtau)
                 do s=jtau,L
                    do kspin=1,Nspin
                       do korb=1,Norb
                          AxB(s) = AxB(s) + A(ispin,kspin,iorb,korb)%lmix(N,s)*B(kspin,jspin,korb,jorb)%mats(s-jtau)
                       enddo
                    enddo
                 enddo
                 C(ispin,jspin,iorb,jorb)%lmix(n,jtau)=C(ispin,jspin,iorb,jorb)%lmix(n,jtau)+dtau*kb_trapz(AxB(0:),jtau,L)
                 !
                 !I2:
                 AxB(0:) = zero
                 do s=1,N
                    do kspin=1,Nspin
                       do korb=1,Norb
                          AxB(s) = AxB(s) + A(ispin,kspin,iorb,korb)%ret(N,s)*B(kspin,jspin,korb,jorb)%lmix(s,jtau)
                       enddo
                    enddo
                 enddo
                 C(ispin,jspin,iorb,jorb)%lmix(N,jtau) = C(ispin,jspin,iorb,jorb)%lmix(N,jtau) + dt*kb_trapz(AxB(0:),1,N)
              enddo
              !
              !Less component
              !C^<(t,t')=-i\int_0^{beta} ds A^\lmix(t,s)*B^\rmix(s,t')
              !             +\int_0^{t'} ds A^<(t,s)*B^A(s,t')
              !             +\int_0^{t} ds A^R(t,s)*B^<(s,t')
              ! (t,t')=>(N,j) <==> Vertical side, with no tip (j=1,N-1)
              C(ispin,jspin,iorb,jorb)%less(N,1:N)=zero
              do j=1,N-1
                 !1.
                 AxB(0:) = zero
                 do s=0,L
                    do kspin=1,Nspin
                       do korb=1,Norb
                          AxB(s) = AxB(s) + A(ispin,kspin,iorb,korb)%lmix(N,s)*conjg( B(kspin,jspin,korb,jorb)%lmix(j,L-s) ) !rmix <-- lmix
                       enddo
                    enddo
                 enddo
                 C(ispin,jspin,iorb,jorb)%less(N,j)=C(ispin,jspin,iorb,jorb)%less(N,j)-xi*dtau*kb_trapz(AxB(0:),0,L)
                 !
                 !2.
                 AxB(0:) = zero
                 do s=1,j
                    do kspin=1,Nspin
                       do korb=1,Norb
                          AxB(s) = AxB(s) + A(ispin,kspin,iorb,korb)%less(N,s)*conjg( B(kspin,jspin,korb,jorb)%ret(j,s) ) !adv <-- ret
                       enddo
                    enddo
                 enddo
                 C(ispin,jspin,iorb,jorb)%less(N,j)=C(ispin,jspin,iorb,jorb)%less(N,j)+dt*kb_trapz(AxB(0:),1,j)
                 !
                 !3.
                 AxB(0:) = zero
                 do s=1,N
                    do kspin=1,Nspin
                       do korb=1,Norb
                          AxB(s) = AxB(s) + A(ispin,kspin,iorb,korb)%ret(N,s)*B(kspin,jspin,korb,jorb)%less(s,j)
                       enddo
                    enddo
                 enddo
                 C(ispin,jspin,iorb,jorb)%less(N,j)=C(ispin,jspin,iorb,jorb)%less(N,j)+dt*kb_trapz(AxB(0:),1,N)
              end do
              !
              !
              ! (t,t')=>(i,N) <==> Horizontal side, w/ tip (i=1,N)
              do i=1,N
                 C(ispin,jspin,iorb,jorb)%less(i,N)=zero
                 !
                 AxB(0:) = zero
                 do s=0,L
                    do kspin=1,Nspin
                       do korb=1,Norb
                          AxB(s) = AxB(s) + A(ispin,kspin,iorb,korb)%lmix(i,s)*conjg( B(kspin,jspin,korb,jorb)%lmix(n,L-s) )
                       enddo
                    enddo
                 enddo
                 C(ispin,jspin,iorb,jorb)%less(i,N)=C(ispin,jspin,iorb,jorb)%less(i,N)-xi*dtau*kb_trapz(AxB(0:),0,L)
                 !
                 AxB(0:) = zero
                 do s=1,N
                    do kspin=1,Nspin
                       do korb=1,Norb
                          AxB(s) = AxB(s) + A(ispin,kspin,iorb,korb)%less(i,s)*conjg( B(kspin,jspin,korb,jorb)%ret(N,s) )
                       enddo
                    enddo
                 enddo
                 C(ispin,jspin,iorb,jorb)%less(i,N)=C(ispin,jspin,iorb,jorb)%less(i,N)+dt*kb_trapz(AxB(0:),1,N)
                 !
                 AxB(0:) = zero
                 do s=1,i
                    do kspin=1,Nspin
                       do korb=1,Norb
                          AxB(s) = AxB(s) + A(ispin,kspin,iorb,korb)%ret(i,s)*B(kspin,jspin,korb,jorb)%less(s,N)
                       enddo
                    enddo
                 enddo
                 C(ispin,jspin,iorb,jorb)%less(i,N)=C(ispin,jspin,iorb,jorb)%less(i,N)+dt*kb_trapz(AxB(0:),1,i)
              end do
              ! C(ispin,jspin,iorb,jorb)%less(1:N-1,N)=zero
              !
              !
              if(present(dcoeff))then
                 C(ispin,jspin,iorb,jorb)%lmix(N,0:L)  = dcoeff*C(ispin,jspin,iorb,jorb)%lmix(N,0:L)
                 C(ispin,jspin,iorb,jorb)%less(N,1:N-1)= dcoeff*C(ispin,jspin,iorb,jorb)%less(N,1:N-1)
                 C(ispin,jspin,iorb,jorb)%less(1:N,N)  = dcoeff*C(ispin,jspin,iorb,jorb)%less(1:N,N)
                 C(ispin,jspin,iorb,jorb)%ret(N,1:N)   = dcoeff*C(ispin,jspin,iorb,jorb)%ret(N,1:N)
              endif
              if(present(ccoeff))then
                 C(ispin,jspin,iorb,jorb)%lmix(N,0:L)  = ccoeff*C(ispin,jspin,iorb,jorb)%lmix(N,0:L)
                 C(ispin,jspin,iorb,jorb)%less(N,1:N-1)= ccoeff*C(ispin,jspin,iorb,jorb)%less(N,1:N-1)
                 C(ispin,jspin,iorb,jorb)%less(1:N,N)  = ccoeff*C(ispin,jspin,iorb,jorb)%less(1:N,N)
                 C(ispin,jspin,iorb,jorb)%ret(N,1:N)   = ccoeff*C(ispin,jspin,iorb,jorb)%ret(N,1:N)
              endif
              !
              !
           enddo jorb_loop
        enddo iorb_loop
     enddo jspin_loop
  enddo ispin_loop
  !
  deallocate(AxB)
  !
end subroutine convolute_kb_contour_gf_Nso
