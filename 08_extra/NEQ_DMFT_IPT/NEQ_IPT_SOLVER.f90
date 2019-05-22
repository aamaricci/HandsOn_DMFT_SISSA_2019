module NEQ_IPT_SOLVER
  USE NEQ_INPUT_VARS
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  implicit none
  private

  public  :: neq_solve_ipt

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Solve with the 2^nd IPT sigma functions
  !+-------------------------------------------------------------------+
  subroutine neq_solve_ipt(params,G0,Sigma)
    type(kb_contour_gf)                   :: G0
    type(kb_contour_gf)                   :: Sigma
    type(kb_contour_params)               :: params
    integer                               :: N,L
    complex(8),dimension(:,:),allocatable :: G0_gtr,Sigma_gtr,G0_rmix
    integer                               :: i,j,itau
    !
    N   = params%Nt                 !<== work with the ACTUAL size of the contour
    L   = params%Ntau
    !
    allocate(G0_gtr(N,N),Sigma_gtr(N,N),G0_rmix(0:L,N))
    do j=1,N
       G0_gtr(N,j)=G0%less(N,j)+G0%ret(N,j)
    end do
    do i=1,N-1
       G0_gtr(i,N)=G0%less(i,n)-conjg(G0%ret(N,i))
    end do
    do j=0,L
       G0_rmix(j,N)  = conjg(G0%lmix(N,L-j))
    enddo

    !Vertical edge
    do j=1,N
       Sigma%less(N,j)= U*U*G0%less(N,j)*G0_gtr(j,N)*G0%less(N,j)
       Sigma_gtr(N,j) = U*U*G0_gtr(N,j)*G0%less(j,N)*G0_gtr(N,j)
    end do
    !Horizontal edge
    do i=1,N-1
       Sigma%less(i,N)= U*U*G0%less(i,N)*G0_gtr(N,i)*G0%less(i,N)
       Sigma_gtr(i,N) = U*U*G0_gtr(i,N)*G0%less(N,i)*G0_gtr(i,N)
    end do
    !Imaginary time edge:
    forall(i=0:L)Sigma%lmix(N,i)  = U*Ui*G0%lmix(N,i)*G0_rmix(i,N)*G0%lmix(N,i)
    forall(j=1:N)Sigma%ret(N,j) = Sigma_gtr(N,j) - Sigma%less(N,j)
  end subroutine neq_solve_ipt


end module NEQ_IPT_SOLVER
