!C(t,t')=A(t,t') + B(t,t'), with t=t_max && t'=0,t_max
! when called multiple times it sums up to a given array. 
subroutine sum_kb_contour_gf_main(params,A,ak,B,bk,C,iaddup)
  type(kb_contour_params)           :: params
  type(kb_contour_gf)               :: A,B
  type(kb_contour_gf),intent(inout) :: C
  real(8)                           :: ak,bk
  integer                           :: i,N,L
  logical,optional                  :: iaddup
  logical                           :: iaddup_
  !
  iaddup_=.false.;if(present(iaddup))iaddup_=iaddup
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  call check_kb_contour_gf(params,A,"sum_kb_contour_gf_main") 
  call check_kb_contour_gf(params,B,"sum_kb_contour_gf_main")
  call check_kb_contour_gf(params,C,"sum_kb_contour_gf_main")
  !
  if(.not.iaddup_)call del_kb_contour_gf(params,C)
  !
  if(N==1)then
     C%mats(0:) = ak*A%mats(0:) + bk*B%mats(0:)
     C%iw(:)    = ak*A%iw(:)    + bk*B%iw(:)
  endif
  C%ret(N,1:N)   = ak*A%ret(N,1:N)   + bk*B%ret(N,1:N)
  C%less(N,1:N)  = ak*A%less(N,1:N)  + bk*B%less(N,1:N)
  C%lmix(N,0:)   = ak*A%lmix(N,0:)   + bk*B%lmix(N,0:)
  !
  !THIS SHOULD NOT BE INVOLVED IN THE CALCULATION:
  C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
  !
end subroutine sum_kb_contour_gf_main

subroutine sum_kb_contour_gf_Nso(params,A,ak,B,bk,C,iaddup)
  type(kb_contour_params)                :: params
  type(kb_contour_gf),dimension(:,:,:,:) :: A,B ![Nspin,Nspin,Norb,Norb]
  type(kb_contour_gf),intent(inout)      :: C(:,:,:,:) !as A
  real(8)                                :: ak,bk
  integer                                :: N,L,Nspin,Norb
  integer                                :: ispin,jspin,iorb,jorb
  logical,optional                       :: iaddup
  logical                                :: iaddup_
  !
  iaddup_=.false.;if(present(iaddup))iaddup_=iaddup
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  Nspin = size(A,1)
  Norb  = size(A,3)
  !
  !
  call check_kb_contour_gf(params,A,"sum_kb_contour_gf_Nso") 
  call check_kb_contour_gf(params,B,"sum_kb_contour_gf_Nso")
  call check_kb_contour_gf(params,C,"sum_kb_contour_gf_Nso")
  !
  if(.not.iaddup_)call del_kb_contour_gf(params,C)
  !
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              if(N==1)then
                 C(ispin,jspin,iorb,jorb)%mats(0:) = ak*A(ispin,jspin,iorb,jorb)%mats(0:) + bk*B(ispin,jspin,iorb,jorb)%mats(0:)
                 C(ispin,jspin,iorb,jorb)%iw(:)    = ak*A(ispin,jspin,iorb,jorb)%iw(:)    + bk*B(ispin,jspin,iorb,jorb)%iw(:)
              endif
              C(ispin,jspin,iorb,jorb)%ret(N,1:N)   = ak*A(ispin,jspin,iorb,jorb)%ret(N,1:N)   + bk*B(ispin,jspin,iorb,jorb)%ret(N,1:N)
              C(ispin,jspin,iorb,jorb)%less(N,1:N)  = ak*A(ispin,jspin,iorb,jorb)%less(N,1:N)  + bk*B(ispin,jspin,iorb,jorb)%less(N,1:N)
              C(ispin,jspin,iorb,jorb)%lmix(N,0:)   = ak*A(ispin,jspin,iorb,jorb)%lmix(N,0:)   + bk*B(ispin,jspin,iorb,jorb)%lmix(N,0:)
              !
              !THIS SHOULD NOT BE INVOLVED IN THE CALCULATION:
              C(ispin,jspin,iorb,jorb)%less(1:N-1,N)= -conjg(C(ispin,jspin,iorb,jorb)%less(N,1:N-1))
              !
           enddo
        enddo
     enddo
  enddo
end subroutine sum_kb_contour_gf_Nso










subroutine sum_kb_contour_gf_recursive_main(params,A,ak,C,iaddup)
  type(kb_contour_gf)               :: A(:)
  real(8)                           :: ak(size(A))
  type(kb_contour_gf),intent(inout) :: C
  type(kb_contour_params)           :: params
  integer                           :: N,L,Na,i
  logical,optional                  :: iaddup
  logical                           :: iaddup_
  !
  iaddup_=.false.;if(present(iaddup))iaddup_=iaddup
  !
  Na=size(A)
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  do i=1,Na
     call check_kb_contour_gf(params,A(i),"sum_kb_contour_gf_resursive_main") 
  enddo
  call check_kb_contour_gf(params,C,"sum_kb_contour_gf_resursive_main") 
  !
  !Reset the result to zero, unless ADDUP is required.
  if(.not.iaddup_)call del_kb_contour_gf(params,C)
  !
  if(N==1)then
     do i=1,Na
        C%mats(0:) = C%mats(0:) + ak(i)*A(i)%mats(0:)
        C%iw(:)    = C%iw(:)    + ak(i)*A(i)%iw(:)
     enddo
  endif
  do i=1,Na
     C%ret(N,1:N)   = C%ret(N,1:N)  + ak(i)*A(i)%ret(N,1:N)
     C%less(N,1:N)  = C%less(N,1:N) + ak(i)*A(i)%less(N,1:N)
     C%lmix(N,0:)   = C%lmix(N,0:)  + ak(i)*A(i)%lmix(N,0:)
  enddo
  !
  C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
  !
end subroutine sum_kb_contour_gf_recursive_main

subroutine sum_kb_contour_gf_recursive_Nso(params,A,ak,C,iaddup)
  type(kb_contour_gf)               :: A(:,:,:,:,:) ![Nspin,Nspin,Norb,Norb,:]
  real(8)                           :: ak(size(A))  !
  type(kb_contour_gf),intent(inout) :: C(:,:,:,:)   ![Nspin,Nspin,Norb,Norb]
  type(kb_contour_params)           :: params
  integer                           :: N,L,Na,i
  integer                           :: Nspin,Norb
  integer                           :: ispin,jspin,iorb,jorb
  logical,optional                  :: iaddup
  logical                           :: iaddup_
  !
  iaddup_=.false.;if(present(iaddup))iaddup_=iaddup
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  Nspin = size(A,1)
  Norb  = size(A,3)
  Na    = size(A,5)
  !
  do i=1,Na
     call check_kb_contour_gf(params,A(:,:,:,:,i),"sum_kb_contour_gf_resursive_main") 
  enddo
  call check_kb_contour_gf(params,C,"sum_kb_contour_gf_resursive_main") 
  !
  !Reset the result to zero, unless ADDUP is required.
  if(.not.iaddup_)call del_kb_contour_gf(params,C)
  !
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              if(N==1)then
                 do i=1,Na
                    C(ispin,jspin,iorb,jorb)%mats(0:) = C(ispin,jspin,iorb,jorb)%mats(0:) + ak(i)*A(ispin,jspin,iorb,jorb,i)%mats(0:)
                    C(ispin,jspin,iorb,jorb)%iw(:)    = C(ispin,jspin,iorb,jorb)%iw(:)    + ak(i)*A(ispin,jspin,iorb,jorb,i)%iw(:)
                 enddo
              endif
              do i=1,Na
                 C(ispin,jspin,iorb,jorb)%ret(N,1:N)   = C(ispin,jspin,iorb,jorb)%ret(N,1:N)  + ak(i)*A(ispin,jspin,iorb,jorb,i)%ret(N,1:N)
                 C(ispin,jspin,iorb,jorb)%less(N,1:N)  = C(ispin,jspin,iorb,jorb)%less(N,1:N) + ak(i)*A(ispin,jspin,iorb,jorb,i)%less(N,1:N)
                 C(ispin,jspin,iorb,jorb)%lmix(N,0:)   = C(ispin,jspin,iorb,jorb)%lmix(N,0:)  + ak(i)*A(ispin,jspin,iorb,jorb,i)%lmix(N,0:)
              enddo
              !
              C(ispin,jspin,iorb,jorb)%less(1:N-1,N)= -conjg(C(ispin,jspin,iorb,jorb)%less(N,1:N-1))
              !
           enddo
        enddo
     enddo
  enddo
end subroutine sum_kb_contour_gf_recursive_Nso









subroutine sum_kb_contour_dgf_main(params,A,ak,B,bk,C,iaddup)
  type(kb_contour_params)           :: params
  type(kb_contour_dgf)               :: A,B
  type(kb_contour_dgf),intent(inout) :: C
  real(8)                           :: ak,bk
  integer                           :: i,N,L
  logical,optional                  :: iaddup
  logical                           :: iaddup_
  !
  iaddup_=.false.;if(present(iaddup))iaddup_=iaddup
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  call check_kb_contour_gf(params,A,"sum_kb_contour_dgf_main") 
  call check_kb_contour_gf(params,B,"sum_kb_contour_dgf_main")
  call check_kb_contour_gf(params,C,"sum_kb_contour_dgf_main")
  !
  if(.not.iaddup_)call del_kb_contour_gf(params,C)
  !
  C%ret(1:N)   = ak*A%ret(1:N)   + bk*B%ret(1:N)
  C%less(1:N)  = ak*A%less(1:N)  + bk*B%less(1:N)
  C%lmix(0:)   = ak*A%lmix(0:)   + bk*B%lmix(0:)
  !
end subroutine sum_kb_contour_dgf_main

subroutine sum_kb_contour_dgf_Nso(params,A,ak,B,bk,C,iaddup)
  type(kb_contour_params)                 :: params
  type(kb_contour_dgf),dimension(:,:,:,:) :: A,B ![Nspin,Nspin,Norb,Norb]
  type(kb_contour_dgf),intent(inout)      :: C(:,:,:,:)
  real(8)                                 :: ak,bk
  integer                                 :: N,L,Nspin,Norb
  integer                                 :: ispin,jspin,iorb,jorb
  logical,optional                        :: iaddup
  logical                                 :: iaddup_
  !
  iaddup_=.false.;if(present(iaddup))iaddup_=iaddup
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  Nspin = size(A,1)
  Norb  = size(A,3)
  !
  !
  call check_kb_contour_gf(params,A,"sum_kb_contour_dgf_Nso") 
  call check_kb_contour_gf(params,B,"sum_kb_contour_dgf_Nso")
  call check_kb_contour_gf(params,C,"sum_kb_contour_dgf_Nso")
  !
  if(.not.iaddup_)call del_kb_contour_gf(params,C)
  !
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              C(ispin,jspin,iorb,jorb)%ret(1:N)   = ak*A(ispin,jspin,iorb,jorb)%ret(1:N)   + bk*B(ispin,jspin,iorb,jorb)%ret(1:N)
              C(ispin,jspin,iorb,jorb)%less(1:N)  = ak*A(ispin,jspin,iorb,jorb)%less(1:N)  + bk*B(ispin,jspin,iorb,jorb)%less(1:N)
              C(ispin,jspin,iorb,jorb)%lmix(0:)   = ak*A(ispin,jspin,iorb,jorb)%lmix(0:)   + bk*B(ispin,jspin,iorb,jorb)%lmix(0:)
           enddo
        enddo
     enddo
  enddo
end subroutine sum_kb_contour_dgf_Nso








