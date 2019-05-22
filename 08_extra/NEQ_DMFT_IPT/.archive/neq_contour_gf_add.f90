!C(t,t')=A(t,t') + B(t,t'), with t=t_max && t'=0,t_max
!t_max_index==N
subroutine add_kb_contour_gf_simple(A,B,C,params)
  type(kb_contour_gf)     :: A,B,C
  type(kb_contour_params) :: params
  integer                 :: N,L
  logical                 :: checkA,checkB,checkC
  N   = params%Ntime   !<== work with the TOTAL size of the contour
  L   = params%Ntau
  if(  (.not.A%status).OR.&
       (.not.B%status).OR.&
       (.not.C%status))stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
  checkA=check_dimension_kb_contour(A,N,L) 
  checkB=check_dimension_kb_contour(B,N,L)
  checkC=check_dimension_kb_contour(C,N,L)
  !
  C%less(:,:) = A%less(:,:) + B%less(:,:)
  C%ret(:,:)  = A%ret(:,:)  + B%ret(:,:)
  C%lmix(:,0:)= A%lmix(:,0:)+ B%lmix(:,0:)
  C%mats(0:)  = A%mats(0:)  + B%mats(0:)
  C%iw(:)     = A%iw(:)     + B%iw(:)
end subroutine  add_kb_contour_gf_simple
!
subroutine add_kb_contour_gf_recursive(A,C,params)
  type(kb_contour_gf)     :: A(:)
  type(kb_contour_gf)     :: C
  type(kb_contour_params) :: params
  integer                 :: i,Na,N,L
  logical                 :: checkA,checkC
  !
  Na=size(A)
  N   = params%Ntime    !<== work with the TOTAL size of the contour
  L   = params%Ntau
  !
  do i=1,Na
     if(  (.not.A(i)%status) )stop "contour_gf/add_kb_contour_gf: A(i) not allocated"
  enddo
  if(  (.not.C%status) )stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
  do i=1,Na
     checkA=check_dimension_kb_contour(A(i),N,L) 
  enddo
  checkC=check_dimension_kb_contour(C,N,L)
  !
  C=zero
  do i=1,Na
     C%less(:,:) = C%less(:,:) + A(i)%less(:,:)
     C%ret(:,:)  = C%ret(:,:)  + A(i)%ret(:,:)
     C%lmix(:,0:)= C%lmix(:,0:)+ A(i)%lmix(:,0:)
     C%mats(0:)  = C%mats(0:)  + A(i)%mats(0:)
     C%iw(:)     = C%iw(:)     + A(i)%iw(:)
  enddo
end subroutine  add_kb_contour_gf_recursive

subroutine add_kb_contour_gf_delta_d(A,B,C,params)
  type(kb_contour_gf)     :: A,C
  real(8),dimension(:)    :: B
  type(kb_contour_params) :: params
  integer                 :: N,L,i
  logical                 :: checkA,checkC
  !
  N   = params%Ntime
  L   = params%Ntau
  !
  if(  (.not.A%status).OR.&
       (.not.C%status))stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
  checkA=check_dimension_kb_contour(A,N,L) 
  checkC=check_dimension_kb_contour(C,N,L)
  !
  if(size(B)<N)stop "contour_gf/add_kb_contour_gf: size(B) < N"
  !
  C = A
  do i=1,N
     C%ret(i,i) = C%ret(i,i) + B(i)
  enddo
  C%mats(0) = C%mats(0) + B(1)
end subroutine  add_kb_contour_gf_delta_d

subroutine add_kb_contour_gf_delta_c(A,B,C,params)
  type(kb_contour_gf)     :: A,C
  complex(8),dimension(:) :: B
  type(kb_contour_params) :: params
  integer                 :: N,L,i
  logical                 :: checkA,checkC
  !
  N   = params%Ntime
  L   = params%Ntau
  !
  if(  (.not.A%status).OR.&
       (.not.C%status))stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
  checkA=check_dimension_kb_contour(A,N,L) 
  checkC=check_dimension_kb_contour(C,N,L)
  !
  if(size(B)<N)stop "contour_gf/add_kb_contour_gf: size(B) < N"
  !
  C = A
  do i=1,N
     C%ret(i,i) = C%ret(i,i) + B(i)
  enddo
  C%mats(0) = C%mats(0) + B(1)
end subroutine  add_kb_contour_gf_delta_c




subroutine add_kb_contour_dgf(A,B,C,params)
  type(kb_contour_dgf)    :: A,B,C
  type(kb_contour_params) :: params
  integer                 :: N,L
  logical                 :: checkA,checkB,checkC
  if(  (.not.A%status).OR.&
       (.not.B%status).OR.&
       (.not.C%status))stop "contour_gf/add_kb_contour_gf: A,B,C not allocated"
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  checkA=check_dimension_kb_contour(A,N,L)
  checkB=check_dimension_kb_contour(B,N,L)
  checkC=check_dimension_kb_contour(C,N,L)
  !
  C%less(:) = A%less(:) + B%less(:)
  C%ret(:)  = A%ret(:)  + B%ret(:)  
  C%lmix(0:)= A%lmix(0:)+ B%lmix(0:)
end subroutine add_kb_contour_dgf
