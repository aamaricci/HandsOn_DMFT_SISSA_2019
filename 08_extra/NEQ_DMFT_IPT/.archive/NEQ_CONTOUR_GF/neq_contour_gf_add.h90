!C(t,t')=A(t,t') + B(t,t'), with t=t_max && t'=0,t_max
subroutine add_kb_contour_gf_main(A,B,C,params)
  type(kb_contour_gf)     :: A,B,C
  type(kb_contour_params) :: params
  integer                 :: N,L
  !
  N   = params%Ntime
  L   = params%Ntau
  !
  call check_kb_contour_gf(params,A,"add_kb_contour_gf_main") 
  call check_kb_contour_gf(params,B,"add_kb_contour_gf_main")
  call check_kb_contour_gf(params,C,"add_kb_contour_gf_main")
  !
  C%less(:,:) = A%less(:,:) + B%less(:,:)
  C%ret(:,:)  = A%ret(:,:)  + B%ret(:,:)
  C%lmix(:,0:)= A%lmix(:,0:)+ B%lmix(:,0:)
  C%mats(0:)  = A%mats(0:)  + B%mats(0:)
  C%iw(:)     = A%iw(:)     + B%iw(:)
  !
end subroutine  add_kb_contour_gf_main

subroutine add_kb_contour_gf_Nso(A,B,C,params)
  type(kb_contour_gf),dimension(:,:,:,:) :: A,B,C ![Nspin][Nspin][Norb][Norb]
  type(kb_contour_params)                :: params
  integer                                :: N,L
  integer                                :: Nspin,Norb
  integer                                :: ispin,jspin,iorb,jorb
  !
  Nspin = size(A,1)
  Norb  = size(A,3)
  !
  N   = params%Ntime
  L   = params%Ntau
  !
  call check_kb_contour_gf(params,A,"add_kb_contour_gf_Nso") 
  call check_kb_contour_gf(params,B,"add_kb_contour_gf_Nso")
  call check_kb_contour_gf(params,C,"add_kb_contour_gf_Nso")
  !
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              C(ispin,jspin,iorb,jorb)%less(:,:) = A(ispin,jspin,iorb,jorb)%less(:,:) + B(ispin,jspin,iorb,jorb)%less(:,:)
              C(ispin,jspin,iorb,jorb)%ret(:,:)  = A(ispin,jspin,iorb,jorb)%ret(:,:)  + B(ispin,jspin,iorb,jorb)%ret(:,:)
              C(ispin,jspin,iorb,jorb)%lmix(:,0:)= A(ispin,jspin,iorb,jorb)%lmix(:,0:)+ B(ispin,jspin,iorb,jorb)%lmix(:,0:)
              C(ispin,jspin,iorb,jorb)%mats(0:)  = A(ispin,jspin,iorb,jorb)%mats(0:)  + B(ispin,jspin,iorb,jorb)%mats(0:)
              C(ispin,jspin,iorb,jorb)%iw(:)     = A(ispin,jspin,iorb,jorb)%iw(:)     + B(ispin,jspin,iorb,jorb)%iw(:)
           enddo
        enddo
     enddo
  enddo
  !
end subroutine  add_kb_contour_gf_Nso







subroutine add_kb_contour_dgf_main(A,B,C,params)
  type(kb_contour_dgf)    :: A,B,C
  type(kb_contour_params) :: params
  integer                 :: N,L
  !
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  call check_kb_contour_gf(params,A,"add_kb_contour_dgf_main") 
  call check_kb_contour_gf(params,B,"add_kb_contour_dgf_main")
  call check_kb_contour_gf(params,C,"add_kb_contour_dgf_main")
  !
  C%less(:) = A%less(:) + B%less(:)
  C%ret(:)  = A%ret(:)  + B%ret(:)  
  C%lmix(0:)= A%lmix(0:)+ B%lmix(0:)
  !
end subroutine add_kb_contour_dgf_main

subroutine add_kb_contour_dgf_Nso(A,B,C,params)
  type(kb_contour_dgf),dimension(:,:,:,:) :: A,B,C ![Nspin][Nspin][Norb][Norb]
  type(kb_contour_params)                :: params
  integer                                :: N,L
  integer                                :: Nspin,Norb
  integer                                :: ispin,jspin,iorb,jorb
  !
  Nspin = size(A,1)
  Norb  = size(A,3)
  !
  N   = params%Ntime
  L   = params%Ntau
  !
  call check_kb_contour_gf(params,A,"add_kb_contour_dgf_Nso") 
  call check_kb_contour_gf(params,B,"add_kb_contour_dgf_Nso")
  call check_kb_contour_gf(params,C,"add_kb_contour_dgf_Nso")
  !
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              C(ispin,jspin,iorb,jorb)%less(:) = A(ispin,jspin,iorb,jorb)%less(:) + B(ispin,jspin,iorb,jorb)%less(:)
              C(ispin,jspin,iorb,jorb)%ret(:)  = A(ispin,jspin,iorb,jorb)%ret(:)  + B(ispin,jspin,iorb,jorb)%ret(:)
              C(ispin,jspin,iorb,jorb)%lmix(0:)= A(ispin,jspin,iorb,jorb)%lmix(0:)+ B(ispin,jspin,iorb,jorb)%lmix(0:)
           enddo
        enddo
     enddo
  enddo
  !
end subroutine  add_kb_contour_dgf_Nso











! subroutine add_kb_contour_gf_recursive(params,A,C)
!   type(kb_contour_params) :: params
!   type(kb_contour_gf)     :: A(:)
!   type(kb_contour_gf)     :: C
!   integer                 :: i,k,Na,N,L
!   logical                 :: checkA,checkC
!   !
!   Na=size(A)
!   N   = params%Ntime    !<== work with the TOTAL size of the contour
!   L   = params%Ntau
!   !
!   do i=1,Na
!      call check_kb_contour_gf(params,A(i),"add_kb_contour_gf_recursive") 
!   enddo
!   call check_kb_contour_gf(params,C,"add_kb_contour_gf_recursive")
!   !
!   C=zero
!   do i=1,Na
!      C%less(:,:) = C%less(:,:) + A(i)%less(:,:)
!      C%ret(:,:)  = C%ret(:,:)  + A(i)%ret(:,:)
!      C%lmix(:,0:)= C%lmix(:,0:)+ A(i)%lmix(:,0:)
!      C%mats(0:)  = C%mats(0:)  + A(i)%mats(0:)
!      C%iw(:)     = C%iw(:)     + A(i)%iw(:)
!      !
!   enddo
! end subroutine  add_kb_contour_gf_recursive
