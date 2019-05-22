!+-----------------------------------------------------------------------------+!
!PURPOSE: define some special operations on the GF:
! + get_adv : obtain the Adv component from G
! + get_rmix: obtain the RightMix component from G
! + get_gtr : obtain the Gtr component from G
! + get_bar : obtain the conjugate component in the Nambu space from G (dubbed bar)
!+-----------------------------------------------------------------------------+!
function get_adv(G,i,j) result (adv)
  implicit none
  type(kb_contour_gf),intent(in) :: G
  integer :: i,j
  complex(8) :: adv
  ! if (.not.G%anomalous) then
  adv=conjg(G%ret(j,i))
  return
  ! else
  !    adv=G%ret(j,i)
  !    return
  ! endif
end function get_adv


function get_rmix(G,i,j,L) result (rmix)
  implicit none
  type(kb_contour_gf),intent(in) :: G
  integer :: i,j,L
  complex(8) :: rmix
  ! if (.not.G%anomalous) then
  rmix=conjg(G%lmix(j,L-i))     !ACTHUNG!!! shouldn't it be -G*(beta-tau)?
  return
  ! else
  !   rmix=G%lmix(j,i)
  !   return
  ! endif
end function get_rmix


function get_gtr(G,i,j) result (gtr)
  implicit none
  type(kb_contour_gf),intent(in) :: G
  integer :: i,j
  complex(8) :: gtr
  ! if (.not.G%anomalous) then
  gtr=G%less(i,j)+G%ret(i,j)
  if (i < j)gtr=G%less(i,j)-conjg(G%ret(j,i))
  return
  ! else
  !   gtr=G%less(j,i)
  !   return
  ! endif
end function get_gtr


! subroutine get_bar_gf(A,B,params)
!   type(kb_contour_gf)                 :: A,B
!   type(kb_contour_params)             :: params
!   integer                             :: N,L,Lf
!   real(8)                             :: dt,dtau
!   integer                             :: i,j,k,itau,jtau
!   logical                             :: checkA,checkB,checkC
!   if(  (.not.A%status).OR.&
!        (.not.B%status))stop "contour_gf/get_bar_kb_contour_gf: A,B,C not allocated"
!   N   = params%Nt      !<== work with the ACTUAL size of the contour
!   L   = params%Ntau
!   Lf  = params%Niw
!   dt  = params%dt
!   dtau= params%dtau
!   checkA=check_dimension_kb_contour(A,params%Ntime,L) 
!   checkB=check_dimension_kb_contour(B,params%Ntime,L)
!   !
!   if (.not.B%anomalous) then
!      if(N==1)then
!         ! Matsubara imaginary time component
!         A%mats(0:L) = B%mats(L:0:-1)
!         ! Matsubara frequencies component
!         A%iw = -conjg(B%iw)
!      endif
!      !Ret. component
!      do j=1,N
!         A%ret(n,j) = -conjg(B%ret(n,j))
!      enddo
!      !Lmix. component
!      do jtau=0,L
!         A%lmix(N,jtau) = -conjg(B%lmix(N,L-jtau))
!      enddo
!      !Less component
!      do j=1,N-1
!         A%less(N,j)=-B%less(j,N)+conjg(B%ret(N,j))
!      end do
!      ! (t,t')=>(i,N) <==> Horizontal side, w/ tip (i=1,N)
!      do j=1,N-1
!         A%less(j,N)=-B%less(N,j)-B%ret(N,j)
!      end do
!      !
!   else
!      !
!      if(N==1)then
!         ! Matsubara imaginary time component
!         A%mats(0:L) = B%mats(0:L)
!         ! Matsubara frequencies component
!         A%iw = B%iw
!      endif
!      !Ret. component
!      do j=1,N
!         A%ret(n,j) =  conjg(B%ret(n,j))
!      enddo
!      !Lmix. component
!      do jtau=0,L
!         A%lmix(N,jtau) =  B%lmix(N,L-jtau)
!      enddo
!      !Less component
!      do j=1,N-1
!         A%less(N,j)=-conjg(B%less(j,N))
!      end do
!      ! (t,t')=>(i,N) <==> Horizontal side, w/ tip (i=1,N)
!      do j=1,N-1
!         B%less(N,j)=-conjg(B%less(j,N))
!      end do
!      !
!   endif
! end subroutine get_bar_gf


! subroutine get_bar_sigma(A,B,params)
!   type(kb_contour_sigma)              :: A,B
!   type(kb_contour_params)             :: params
!   integer                             :: N,L,Lf
!   real(8)                             :: dt,dtau
!   integer                             :: i,j,k,itau,jtau
!   logical                             :: checkA,checkB,checkC
!   if(  (.not.A%status).OR.&
!        (.not.B%status))stop "contour_gf/get_bar_kb_contour_gf: A,B not allocated"
!   N   = params%Nt      !<== work with the ACTUAL size of the contour
!   L   = params%Ntau
!   Lf  = params%Niw
!   dt  = params%dt
!   dtau= params%dtau
!   checkA=check_dimension_kb_contour(A,params%Ntime,L) 
!   checkB=check_dimension_kb_contour(B,params%Ntime,L)
!   if (.not.B%reg%anomalous) then
!      if(N==1)then
!         ! Matsubara imaginary time component
!         A%reg%mats(0:L) = B%reg%mats(L:0:-1)
!         ! Matsubara frequencies component
!         A%reg%iw = -conjg(B%reg%iw)
!      endif
!      !Ret. component
!      do j=1,N
!         A%reg%ret(n,j) = -conjg(B%reg%ret(n,j))
!      enddo
!      !Lmix. component
!      do jtau=0,L
!         A%reg%lmix(N,jtau) = -conjg(B%reg%lmix(N,L-jtau))
!      enddo
!      !Less component
!      do j=1,N-1
!         A%reg%less(N,j)=-B%reg%less(j,N)+conjg(B%reg%ret(N,j))
!      end do
!      ! (t,t')=>(i,N) <==> Horizontal side, w/ tip (i=1,N)
!      do j=1,N-1
!         A%reg%less(j,N)=-B%reg%less(N,j)-B%reg%ret(N,j)
!      end do
!      !HF components:
!      A%hf(:) = -conjg(B%hf(:))
!      !
!   else
!      !
!      if(N==1)then
!         ! Matsubara imaginary time component
!         A%reg%mats(0:L) = B%reg%mats(0:L)
!         ! Matsubara frequencies component
!         A%reg%iw = B%reg%iw
!      endif
!      !Ret. component
!      do j=1,N
!         A%reg%ret(n,j) =  conjg(B%reg%ret(n,j))
!      enddo
!      !Lmix. component
!      do jtau=0,L
!         A%reg%lmix(N,jtau) =  B%reg%lmix(N,L-jtau)
!      enddo
!      !Less component
!      do j=1,N-1
!         A%reg%less(N,j)=-conjg(B%reg%less(j,N))
!      end do
!      ! (t,t')=>(i,N) <==> Horizontal side, w/ tip (i=1,N)
!      do j=1,N-1
!         B%reg%less(N,j)=-conjg(B%reg%less(j,N))
!      end do
!      !HF components:
!      A%hf(:) = B%hf(:)
!      !	
!   endif
! end subroutine get_bar_sigma




!+-----------------------------------------------------------------------------+!
!PURPOSE: define some standard operations useful in the code such as
! + equalities among contour GF.
! + product of contour GF times a scalar.
!+-----------------------------------------------------------------------------+!
!
!EQUALITIES:
subroutine kb_contour_gf_equality_(G1,C)
  type(kb_contour_gf),intent(inout) :: G1
  complex(8),intent(in)             :: C
  G1%less(:,:)  = C
  G1%ret(:,:)   = C
  G1%lmix(:,0:) = C
  G1%mats(0:)   = C
  G1%iw(:)      = C
end subroutine kb_contour_gf_equality_
!
subroutine kb_contour_gf_equality__(G1,G2)
  type(kb_contour_gf),intent(inout) :: G1
  type(kb_contour_gf),intent(in)    :: G2
  G1%less(:,:)  = G2%less(:,:)
  G1%ret(:,:)   = G2%ret(:,:)
  G1%lmix(:,0:) = G2%lmix(:,0:)
  G1%mats(0:)   = G2%mats(0:)
  G1%iw(:)      = G2%iw(:)
end subroutine kb_contour_gf_equality__
!
subroutine kb_contour_dgf_equality_(dG1,C)
  type(kb_contour_dgf),intent(inout) :: dG1
  complex(8),intent(in)             :: C
  dG1%less(:)  = C
  dG1%ret(:)   = C
  dG1%lmix(0:) = C
end subroutine kb_contour_dgf_equality_
!
subroutine kb_contour_dgf_equality__(dG1,dG2)
  type(kb_contour_dgf),intent(inout) :: dG1
  type(kb_contour_dgf),intent(in)    :: dG2
  dG1%less(:)  = dG2%less(:)
  dG1%ret(:)   = dG2%ret(:)
  dG1%lmix(0:) = dG2%lmix(0:)
end subroutine kb_contour_dgf_equality__
!
!
!PRODUCT times SCALAR:
function kb_contour_gf_scalarL_d(C,G) result(F)
  real(8),intent(in)             :: C
  type(kb_contour_gf),intent(in) :: G
  type(kb_contour_gf)            :: F
  F%less(:,:) = C*G%less(:,:)
  F%ret(:,:)  = C*G%ret(:,:)
  F%lmix(:,0:)= C*G%lmix(:,0:)
  F%mats(0:)  = C*G%mats(0:)
  F%iw(:)     = C*G%iw(:)
end function kb_contour_gf_scalarL_d
!
function kb_contour_dgf_scalarL_d(C,dG) result(dF)
  real(8),intent(in)             :: C
  type(kb_contour_dgf),intent(in) :: dG
  type(kb_contour_dgf)            :: dF
  dF%less(:) = C*dG%less(:)
  dF%ret(:)  = C*dG%ret(:)
  dF%lmix(0:)= C*dG%lmix(0:)
end function kb_contour_dgf_scalarL_d
!
function kb_contour_gf_scalarL_c(C,G) result(F)
  complex(8),intent(in)          :: C
  type(kb_contour_gf),intent(in) :: G
  type(kb_contour_gf)            :: F
  F%less(:,:) = C*G%less(:,:)
  F%ret(:,:)  = C*G%ret(:,:)
  F%lmix(:,0:)= C*G%lmix(:,0:)
  F%mats(0:)  = C*G%mats(0:)
  F%iw(:)     = C*G%iw(:)
end function kb_contour_gf_scalarL_c
!
function kb_contour_dgf_scalarL_c(C,dG) result(dF)
  complex(8),intent(in)           :: C
  type(kb_contour_dgf),intent(in) :: dG
  type(kb_contour_dgf)            :: dF
  dF%less(:) = C*dG%less(:)
  dF%ret(:)  = C*dG%ret(:)
  dF%lmix(0:)= C*dG%lmix(0:)
end function kb_contour_dgf_scalarL_c
!
function kb_contour_gf_scalarR_d(G,C) result(F)
  real(8),intent(in)             :: C
  type(kb_contour_gf),intent(in) :: G
  type(kb_contour_gf)            :: F
  F%less(:,:) = G%less(:,:)*C
  F%ret(:,:)  = G%ret(:,:)*C
  F%lmix(:,0:)= G%lmix(:,0:)*C
  F%mats(0:)  = G%mats(0:)*C
  F%iw(:)     = G%iw(:)*C
end function kb_contour_gf_scalarR_d
!
function kb_contour_dgf_scalarR_d(dG,C) result(dF)
  real(8),intent(in)             :: C
  type(kb_contour_dgf),intent(in) :: dG
  type(kb_contour_dgf)            :: dF
  dF%less(:) = dG%less(:)*C
  dF%ret(:)  = dG%ret(:)*C
  dF%lmix(0:)= dG%lmix(0:)*C
end function kb_contour_dgf_scalarR_d
!
function kb_contour_gf_scalarR_c(G,C) result(F)
  complex(8),intent(in)             :: C
  type(kb_contour_gf),intent(in) :: G
  type(kb_contour_gf)            :: F
  F%less(:,:) = G%less(:,:)*C
  F%ret(:,:)  = G%ret(:,:)*C
  F%lmix(:,0:)= G%lmix(:,0:)*C
  F%mats(0:)  = G%mats(0:)*C
  F%iw(:)     = G%iw(:)*C
end function kb_contour_gf_scalarR_c
!
function kb_contour_dgf_scalarR_c(dG,C) result(dF)
  complex(8),intent(in)             :: C
  type(kb_contour_dgf),intent(in) :: dG
  type(kb_contour_dgf)            :: dF
  dF%less(:) = dG%less(:)*C
  dF%ret(:)  = dG%ret(:)*C
  dF%lmix(0:)= dG%lmix(0:)*C
end function kb_contour_dgf_scalarR_c
