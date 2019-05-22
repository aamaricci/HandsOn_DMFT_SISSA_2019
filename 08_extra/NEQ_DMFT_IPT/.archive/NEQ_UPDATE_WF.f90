!###############################################################
!     PURPOSE  : Constructs some functions used in other places. 
!     AUTHORS  : Adriano Amaricci
!###############################################################
module NEQ_UPDATE_WF
  USE MATRIX
  USE NEQ_VARS_GLOBAL
  implicit none
  private

  public :: neq_update_weiss_field

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Update the Weiss Fields G0^{>,<,Ret} using Dyson equation
  !+-------------------------------------------------------------------+
  subroutine neq_update_weiss_field
    integer :: M,i,j,k,itau,jtau,NN
    real(8) :: R,deg
    real(8) :: w,A,An
    complex(8),dimension(:,:),allocatable :: locGgtr,Sgtr
    complex(8),dimension(:,:),allocatable :: Uno,GammaRet
    !
    complex(8),dimension(:,:),allocatable :: mat_Delta
    complex(8),dimension(:,:),allocatable :: mat_Gamma
    complex(8),dimension(:,:),allocatable :: mat_G0,mat_Sigma,mat_locG
    !
    type(keldysh_contour_gf),save         :: G0_old
    !
    real(8),dimension(Nstep) :: nt
    if(G0_old%status.EQV..false.)call allocate_keldysh_contour_gf(G0_old,Nstep)
    G0_old=G0    

    !=======Component by component inversion==========================
    select case(fupdate)
    case default
       !Matrix update, from testKELDYSHMATGF4
       call msg("update WF with direct inversion of keldysh matrix.")
       !Build Gloc matrix
       allocate(mat_locG(2*nstep,2*nstep))
       mat_locG  = build_keldysh_matrix_gf(locG,Nstep)
       !Build Sigma matrix
       allocate(mat_Sigma(2*Nstep,2*Nstep))
       mat_Sigma = build_keldysh_matrix_gf(Sigma,Nstep)
       !Get G0 matrix:
       allocate(mat_G0(2*nstep,2*nstep))
       mat_locG  = mat_locG*dt**2
       call matrix_inverse(mat_locG)
       mat_G0 = mat_locG + mat_Sigma
       mat_G0 = mat_G0*dt**2
       call matrix_inverse(mat_G0)
       G0%less =  mat_G0(1:Nstep,Nstep+1:2*Nstep)/2.d0
       G0%ret  =  mat_G0(1:Nstep,1:Nstep)
       deallocate(mat_locG,mat_Sigma,mat_G0)

    case(1)  
       !this one can be updated or modified using component by component 
       !inversion which should give better handling of divergence issues.
       call msg("update WF with component-by-component inversion.")
       allocate(locGgtr(Nstep,Nstep),Sgtr(Nstep,Nstep))
       locGgtr = locG%less  + locG%ret  - conjg(transpose(locG%ret))
       Sgtr    = Sigma%less + Sigma%ret - conjg(transpose(Sigma%ret))
       ! !G0ret = [\11 + Gret * Sret]^-1 * Gret
       allocate(Uno(Nstep,Nstep),GammaRet(Nstep,Nstep))
       Uno=zero  ; forall(i=1:nstep)Uno(i,i)=One/dt
       GammaRet = Uno+matmul(locG%ret,Sigma%ret)*dt
       GammaRet = GammaRet*dt**2
       call matrix_inverse(GammaRet)
       G0%ret = matmul(GammaRet,locG%ret)*dt
       !### COMMENTING THIS LINE THE RESULTS ARE IDENTICAL WITH THE THREE METHODS OF UPDATE ###
       !forall(i=1:nstep)G0%ret(i,i)=-xi !???
       !#####################################################################################
       !G0less = GammaR^-1 * Gless * GammaA^-1  -  gR * Sless * gA
       G0%less = matmul(GammaRet,matmul(locG%less,conjg(transpose(GammaRet)))*dt)*dt -&
            matmul(G0%ret,matmul(Sigma%less,conjg(transpose(G0%ret)))*dt)*dt
       ! !G0gtr  = GammaR^-1 * Ggtr * GammaA^-1   -  gR * Sgtr * gA
       ! G0%gtr(0:nstep,0:nstep)  = matmul(GammaRet(0:nstep,0:nstep),matmul(locG%gtr(0:nstep,0:nstep),&
       !      conjg(transpose(GammaRet(0:nstep,0:nstep))))*dt)*dt  -&
       !      matmul(G0ret(0:nstep,0:nstep),matmul(Sigma%gtr(0:nstep,0:nstep),G0adv(0:nstep,0:nstep))*dt)*dt



       ! case(2)
       !    call msg("update with method 2: inversion and matrix-multiplication")
       !    !Build Gloc matrix
       !    allocate(mat_locG(2*nstep,2*nstep))
       !    mat_locG = build_keldysh_matrix_gf(locG,nstep)
       !    !Build Sigma matrix
       !    allocate(mat_Sigma(2*nstep,2*nstep))
       !    mat_Sigma = build_keldysh_matrix_gf(Sigma,nstep)
       !    !Allocate space for other matrices:
       !    allocate(mat_Delta(2*nstep,2*nstep))
       !    allocate(mat_Gamma(2*nstep,2*nstep))
       !    allocate(mat_G0(2*nstep,2*nstep))
       !    mat_Delta= zero ; forall(i=1:2*nstep)mat_Delta(i,i)=One/dt
       !    mat_Gamma= mat_Delta + matmul(mat_Sigma,mat_locG)*dt
       !    mat_Gamma= mat_Gamma*dt**2
       !    call matrix_inverse(mat_Gamma)
       !    mat_G0  =  matmul(mat_locG,mat_Gamma)*dt
       !    G0%less =  mat_G0(1:Nstep,Nstep+1:2*Nstep)/2.d0
       !    G0%ret  =  mat_G0(1:Nstep,1:Nstep)
       !    deallocate(mat_locG,mat_Sigma,mat_G0,mat_Delta,mat_Gamma)
    end select


    G0%less = weight*G0%less + (1.d0-weight)*G0_old%less
    G0%ret  = weight*G0%ret  + (1.d0-weight)*G0_old%ret

    !Save data:
    if(plot3D)call plot_keldysh_contour_gf(G0,time,reg(plot_dir)//"/G0")

  end subroutine neq_update_weiss_field



  !********************************************************************
  !********************************************************************
  !********************************************************************


  ! function build_keldysh_matrix_gf(G,N) result(matG)
  !   type(keldysh_contour_gf)      :: G
  !   complex(8),dimension(2*N,2*N) :: matG
  !   integer                       :: i,j,N
  !   matG=zero
  !   forall(i=1:N,j=1:N)
  !      matG(i,j)         = step(t(i)-t(j))*G%gtr(i,j) + step(t(j)-t(i))*G%less(i,j)
  !      matG(i,N+1+j)     =-G%less(i,j)
  !      matG(N+1+i,j)     = G%gtr(i,j)
  !      matG(N+1+i,N+1+j) =-(step(t(i)-t(j))*G%less(i,j)+ step(t(j)-t(i))*G%gtr(i,j))
  !   end forall
  ! end function build_keldysh_matrix_gf

  function build_keldysh_matrix_gf(G,N) result(matG)
    type(keldysh_contour_gf)      :: G
    complex(8),dimension(2*N,2*N) :: matG
    integer                       :: i,j,N
    matG=zero
    forall(i=1:N,j=1:N)
       matG(i,j)      = G%ret(i,j)        !G_11=G^R
       matG(N+i,N+j)  = conjg(G%ret(j,i)) !G_22=G^A=[G^R]^+
       matG(i,N+j)    = 2.d0*G%less(i,j)  !G_12=2G^<
    end forall
  end function build_keldysh_matrix_gf



  !********************************************************************
  !********************************************************************
  !********************************************************************









end module NEQ_UPDATE_WF
