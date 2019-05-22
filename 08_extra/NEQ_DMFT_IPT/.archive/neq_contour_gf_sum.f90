! performs the sum a*A + b*B and stores it in C along the perimeter
! can be called multiple times to add up. 
subroutine sum_kb_contour_gf_simple(A,ak,B,bk,C,params,iaddup)
  type(kb_contour_gf)               :: A,B
  type(kb_contour_gf),intent(inout) :: C
  real(8)                           :: ak,bk
  type(kb_contour_params)           :: params
  integer                           :: N,L
  logical,optional                  :: iaddup
  logical                           :: iaddup_
  iaddup_=.false.;if(present(iaddup))iaddup_=iaddup
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  if(  (.not.A%status).OR.&
       (.not.B%status).OR.&
       (.not.C%status))stop "contour_gf/addup_kb_contour_gf: G or Gk not allocated"
  !
  if(.not.iaddup_)call del_kb_contour_gf(C,params)
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
  ! if(.not.C%anomalous)then
  C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
  ! else
  !    C%less(1:N-1,N)= C%less(N,1:N-1)+C%ret(N,1:N-1)
  ! endif
end subroutine sum_kb_contour_gf_simple

subroutine sum_kb_contour_gf_recursive(A,ak,C,params,iaddup)
  type(kb_contour_gf)               :: A(:)
  real(8)                           :: ak(size(A))
  type(kb_contour_gf),intent(inout) :: C
  type(kb_contour_params)           :: params
  integer                           :: N,L,Na,i
  logical,optional                  :: iaddup
  logical                           :: iaddup_
  iaddup_=.false.;if(present(iaddup))iaddup_=iaddup
  !
  Na=size(A)
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  do i=1,Na
     if((.not.A(i)%status) )stop "contour_gf/sum_kb_contour_gf: G or Gk not allocated"
  enddo
  if( (.not.C%status) )stop "contour_gf/addup_kb_contour_gf: G or Gk not allocated"
  if(.not.iaddup_)call del_kb_contour_gf(C,params)
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
  ! if(.not.C%anomalous)then
  C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
  ! else
  !    C%less(1:N-1,N)= C%less(N,1:N-1)+C%ret(N,1:N-1)
  ! endif
end subroutine sum_kb_contour_gf_recursive

subroutine sum_kb_contour_gf_delta_d(A,ak,B,bk,C,params,iaddup)
  type(kb_contour_gf)               :: A
  real(8),dimension(:)              :: B
  type(kb_contour_gf),intent(inout) :: C
  real(8)                           :: ak,bk
  type(kb_contour_params)           :: params
  integer                           :: N,L
  logical,optional                  :: iaddup
  logical                           :: iaddup_
  iaddup_=.false.;if(present(iaddup))iaddup_=iaddup
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  if(  (.not.A%status).OR.&
       (.not.C%status))stop "contour_gf/addup_kb_contour_gf: G or Gk not allocated"
  !
  if(.not.iaddup_)call del_kb_contour_gf(C,params)
  !
  if(N==1)then 
     C%mats(0)  = ak*A%mats(0)  + B(1)
     C%mats(1:) = ak*A%mats(1:)
  endif
  C%ret(N,1:N)   = ak*A%ret(N,1:N)
  C%less(N,1:N)  = ak*A%less(N,1:N)
  C%lmix(N,0:)   = ak*A%lmix(N,0:)
  !
  ! if(.not.C%anomalous)then
  C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
  ! else
  !    C%less(1:N-1,N)= C%less(N,1:N-1)+C%ret(N,1:N-1)
  ! endif
  C%ret(N,N) = C%ret(N,N) + B(N)
end subroutine sum_kb_contour_gf_delta_d

subroutine sum_kb_contour_gf_delta_c(A,ak,B,bk,C,params,iaddup)
  type(kb_contour_gf)               :: A
  complex(8),dimension(:)           :: B
  type(kb_contour_gf),intent(inout) :: C
  real(8)                           :: ak,bk
  type(kb_contour_params)           :: params
  integer                           :: N,L
  logical,optional                  :: iaddup
  logical                           :: iaddup_
  iaddup_=.false.;if(present(iaddup))iaddup_=iaddup
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  if(  (.not.A%status).OR.&
       (.not.C%status))stop "contour_gf/addup_kb_contour_gf: G or Gk not allocated"
  !
  if(.not.iaddup_)call del_kb_contour_gf(C,params)
  !
  if(N==1)then 
     C%mats(0)  = ak*A%mats(0)  + B(1)
     C%mats(1:) = ak*A%mats(1:)
  endif
  C%ret(N,1:N)   = ak*A%ret(N,1:N) 
  C%less(N,1:N)  = ak*A%less(N,1:N)
  C%lmix(N,0:)   = ak*A%lmix(N,0:)
  !
  ! if(.not.C%anomalous)then
  C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
  ! else
  !    C%less(1:N-1,N)= C%less(N,1:N-1)+C%ret(N,1:N-1)
  ! endif
  C%ret(N,N) = C%ret(N,N) + B(N)
end subroutine sum_kb_contour_gf_delta_c











!======= DEL ======= 
! reset the function along the actual perimeter to zero:
subroutine del_kb_contour_gf(G,params)
  type(kb_contour_gf)     :: G
  type(kb_contour_params) :: params
  integer                 :: N,L
  if(  (.not.G%status)) stop "contour_gf/addup_kb_contour_gf: G or Gk not allocated"
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  if(N==1)then
     G%iw = zero
     G%mats = 0d0
  endif
  G%ret(N,1:N)   = zero
  G%less(N,1:N)  = zero
  G%lmix(N,0:)   = zero
end subroutine del_kb_contour_gf

subroutine del_kb_contour_dgf(dG,params)
  type(kb_contour_dgf)    :: dG
  type(kb_contour_params) :: params
  integer                 :: N,L
  if(  (.not.dG%status))stop "contour_gf/del_kb_contour_dgf: dG not allocated"
  N   = params%Nt                 !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  dG%less(1:N) = zero
  dG%ret(1:N)  = zero
  dG%lmix(0:)  = zero
end subroutine del_kb_contour_dgf
