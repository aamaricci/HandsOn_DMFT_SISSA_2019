MODULE NEQ_CONTOUR_GF
  USE NEQ_CONTOUR
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS
  USE SF_LINALG,    only: zeye,inv
  USE SF_MISC, only: assert_shape
  implicit none
  private

  !>TODO:
  !for intrinsic operations look here:
  !https://stackoverflow.com/questions/38230199/fortran-function-to-overload-multiplication-between-derived-types-with-allocatab
  !
  !add an assert_shape procedure:



  ! KADANOFF-BAYM CONTOUR GREEN'S FUNCTIONS:
  type :: kb_contour_gf
     complex(8),dimension(:,:),allocatable :: less
     complex(8),dimension(:,:),allocatable :: ret
     complex(8),dimension(:,:),allocatable :: lmix
     real(8),dimension(:),allocatable      :: mats
     complex(8),dimension(:),allocatable   :: iw
     logical                               :: status=.false.
  end type kb_contour_gf
  !
  ! KADANOFF-BAYM CONTOUR GREEN'S FUNCTIONS DERIVATIVE
  type :: kb_contour_dgf
     complex(8),dimension(:),allocatable :: less,gtr
     complex(8),dimension(:),allocatable :: ret
     complex(8),dimension(:),allocatable :: lmix
     logical                             :: status=.false.
  end type kb_contour_dgf
  !
  !
  ! ALLOCATION/DEALLOCATION ROUTINES:
  interface allocate_kb_contour_gf
     module procedure :: allocate_kb_contour_gf_main
     module procedure :: allocate_kb_contour_gf_Nso
     module procedure :: allocate_kb_contour_gf_Nlk_main
     module procedure :: allocate_kb_contour_gf_Nlk_Nso
     !
     module procedure :: allocate_kb_contour_dgf_main
     module procedure :: allocate_kb_contour_dgf_Nso
     module procedure :: allocate_kb_contour_dgf_Nlk_main
     module procedure :: allocate_kb_contour_dgf_Nlk_Nso
  end interface allocate_kb_contour_gf
  !
  interface deallocate_kb_contour_gf
     module procedure :: deallocate_kb_contour_gf_main
     module procedure :: deallocate_kb_contour_gf_Nso
     module procedure :: deallocate_kb_contour_gf_Nlk_main
     module procedure :: deallocate_kb_contour_gf_Nlk_Nso
     !
     module procedure :: deallocate_kb_contour_dgf_main
     module procedure :: deallocate_kb_contour_dgf_Nso
     module procedure :: deallocate_kb_contour_dgf_Nlk_main
     module procedure :: deallocate_kb_contour_dgf_Nlk_Nso
  end interface deallocate_kb_contour_gf
  !
  !CHECK:
  interface check_kb_contour_gf
     module procedure :: check_kb_contour_gf_main
     module procedure :: check_kb_contour_gf_Nso
     module procedure :: check_kb_contour_gf_Nlk_main
     module procedure :: check_kb_contour_gf_Nlk_Nso
     !
     module procedure :: check_kb_contour_dgf_main
     module procedure :: check_kb_contour_dgf_Nso
     module procedure :: check_kb_contour_dgf_Nlk_main
     module procedure :: check_kb_contour_dgf_Nlk_Nso
  end interface check_kb_contour_gf
  !
  ! ADD (TOTAL DOMAIN) ROUTINES:
  interface add_kb_contour_gf
     module procedure :: add_kb_contour_gf_main
     module procedure :: add_kb_contour_gf_Nso
     module procedure :: add_kb_contour_gf_Nlk_main
     module procedure :: add_kb_contour_gf_Nlk_Nso
     !
     module procedure :: add_kb_contour_dgf_main
     module procedure :: add_kb_contour_dgf_Nso
     module procedure :: add_kb_contour_dgf_Nlk_main
     module procedure :: add_kb_contour_dgf_Nlk_Nso
  end interface add_kb_contour_gf
  !
  ! SUM (PERIMETER) ROUTINES:
  interface sum_kb_contour_gf
     module procedure :: sum_kb_contour_gf_main
     module procedure :: sum_kb_contour_gf_Nso
     module procedure :: sum_kb_contour_gf_Nlk_main
     module procedure :: sum_kb_contour_gf_Nlk_Nso
     !
     module procedure :: sum_kb_contour_gf_recursive_main
     module procedure :: sum_kb_contour_gf_recursive_Nso
     !
     module procedure :: sum_kb_contour_dgf_main
     module procedure :: sum_kb_contour_dgf_Nso
     module procedure :: sum_kb_contour_dgf_Nlk_main
     module procedure :: sum_kb_contour_dgf_Nlk_Nso
  end interface sum_kb_contour_gf
  !
  ! DEL (PERIMETER) ROUTINES:
  interface del_kb_contour_gf
     module procedure :: del_kb_contour_gf_main
     module procedure :: del_kb_contour_gf_Nso
     module procedure :: del_kb_contour_gf_Nlk_main
     module procedure :: del_kb_contour_gf_Nlk_Nso
     !
     module procedure :: del_kb_contour_dgf_main
     module procedure :: del_kb_contour_dgf_Nso
     module procedure :: del_kb_contour_dgf_Nlk_main
     module procedure :: del_kb_contour_dgf_Nlk_Nso
  end interface del_kb_contour_gf
  !
  ! CONVOLUTION:
  interface convolute_kb_contour_gf
     module procedure ::  convolute_kb_contour_gf_main
     module procedure ::  convolute_kb_contour_gf_Nso_
     module procedure ::  convolute_kb_contour_gf_Nlk_main
     module procedure ::  convolute_kb_contour_gf_Nlk_Nso
  end interface convolute_kb_contour_gf
  !
  !VIE SOLVER:
  interface vie_kb_contour_gf
     module procedure :: vie_kb_contour_gf_main
     module procedure :: vie_kb_contour_gf_Nso_
     module procedure :: vie_kb_contour_gf_Nlk_main
     module procedure :: vie_kb_contour_gf_Nlk_Nso
  end interface vie_kb_contour_gf
  !
  !VIDE SOLVER:
  interface vide_kb_contour_gf
     module procedure :: vide_kb_contour_gf_main
     module procedure :: vide_kb_contour_gf_Nso_
     module procedure :: vide_kb_contour_gf_Nk_main
     module procedure :: vide_kb_contour_gf_Nk_main_Kscalar
     module procedure :: vide_kb_contour_gf_Nk_Nso
     module procedure :: vide_kb_contour_gf_Nk_Nso_Kscalar
  end interface vide_kb_contour_gf
  !
  !SAVE:
  interface save_kb_contour_gf
     module procedure :: save_kb_contour_gf_main
     module procedure :: save_kb_contour_gf_Nso
     module procedure :: save_kb_contour_dgf_main
     module procedure :: save_kb_contour_dgf_Nso
  end interface save_kb_contour_gf
  !
  !READ:
  interface read_kb_contour_gf
     module procedure :: read_kb_contour_gf_main
     module procedure :: read_kb_contour_gf_Nso
     module procedure :: read_kb_contour_dgf_main
     module procedure :: read_kb_contour_dgf_Nso
  end interface read_kb_contour_gf
  !
  !PLOT:
  interface plot_kb_contour_gf
     module procedure :: plot_kb_contour_gf_main
     module procedure :: plot_kb_contour_gf_Nso
     module procedure :: plot_kb_contour_dgf_main
     module procedure :: plot_kb_contour_dgf_Nso
  end interface plot_kb_contour_gf
  !
  !ASSERT_SHAPE:
  interface assert_shape_kb_contour_gf
     module procedure :: kb_contour_gf_assert_shape_N1
     module procedure :: kb_contour_gf_assert_shape_N2
     module procedure :: kb_contour_gf_assert_shape_N3
     module procedure :: kb_contour_gf_assert_shape_N4
     module procedure :: kb_contour_gf_assert_shape_N5
     module procedure :: kb_contour_gf_assert_shape_N6
     module procedure :: kb_contour_gf_assert_shape_N7
     !
     module procedure :: kb_contour_dgf_assert_shape_N1
     module procedure :: kb_contour_dgf_assert_shape_N2
     module procedure :: kb_contour_dgf_assert_shape_N3
     module procedure :: kb_contour_dgf_assert_shape_N4
     module procedure :: kb_contour_dgf_assert_shape_N5
     module procedure :: kb_contour_dgf_assert_shape_N6
     module procedure :: kb_contour_dgf_assert_shape_N7
  end interface assert_shape_kb_contour_gf
  !
  !EXTRAPOLATION:
  interface extrapolate_kb_contour_gf
     module procedure :: extrapolate_kb_contour_gf_main
     module procedure :: extrapolate_kb_contour_gf_Nso
  end interface extrapolate_kb_contour_gf
  !
  !INTEGRATION ROUTINES:
  interface kb_trapz
     module procedure :: kb_trapz_d
     module procedure :: kb_trapz_c
  end interface kb_trapz

  interface kb_half_trapz
     module procedure :: kb_half_trapz_d
     module procedure :: kb_half_trapz_c
  end interface kb_half_trapz
  !
  ! !OVERLOAD OPERATORS  
  interface assignment(=)
     module procedure :: kb_contour_gf_equality_gf_scalar_main
     module procedure :: kb_contour_gf_equality_gf_scalar_Nso
     module procedure :: kb_contour_gf_equality_gf_gf_main
     module procedure :: kb_contour_gf_equality_gf_gf_Nso
     !
     module procedure :: kb_contour_dgf_equality_dgf_scalar_main
     module procedure :: kb_contour_dgf_equality_dgf_scalar_Nso
     module procedure :: kb_contour_dgf_equality_dgf_dgf_main
     module procedure :: kb_contour_dgf_equality_dgf_dgf_Nso
  end interface assignment(=)
  !
  ! interface operator(+)
  !    module procedure :: kb_contour_gf_sum_main
  !    module procedure :: kb_contour_gf_sum_Nso
  ! end interface operator(+)
  !
  ! interface operator(.add.)
  !    module procedure :: kb_contour_gf_add_main
  !    module procedure :: kb_contour_gf_add_Nso
  ! end interface operator(.add.)
  !
  ! interface operator(*)
  !    module procedure ::  kb_contour_gf_ScalarProdLeft_d_main
  !    module procedure ::  kb_contour_gf_ScalarProdLeft_c_main
  !    module procedure ::  kb_contour_gf_ScalarProdLeft_d_Nso
  !    module procedure ::  kb_contour_gf_ScalarProdLeft_c_Nso
  !    !
  !    module procedure ::  kb_contour_dgf_ScalarProdLeft_d_main
  !    module procedure ::  kb_contour_dgf_ScalarProdLeft_c_main
  !    module procedure ::  kb_contour_dgf_ScalarProdLeft_d_Nso
  !    module procedure ::  kb_contour_dgf_ScalarProdLeft_c_Nso
  ! end interface operator(*)
  !
  ! interface operator(.star.)
  !    module procedure :: f_convolute_kb_contour_gf_main
  !    ! module procedure :: f_convolute_kb_contour_gf_Nso
  ! end interface operator(.star.)


  public :: kb_contour_gf
  public :: kb_contour_dgf
  !
  public :: set_kb_contour_gf_time
  !
  public :: allocate_kb_contour_gf
  public :: deallocate_kb_contour_gf
  !
  public :: check_kb_contour_gf
  !
  public :: add_kb_contour_gf
  !
  public :: sum_kb_contour_gf
  public :: del_kb_contour_gf
  !
  public :: convolute_kb_contour_gf,convolute_kb_contour_gf_Nso
  !
  public :: vie_kb_contour_gf,vie_kb_contour_gf_Nso
  !
  public :: vide_kb_contour_gf,vide_kb_contour_gf_Nso
  !
  public :: assert_shape_kb_contour_gf
  public :: save_kb_contour_gf
  public :: inquire_kb_contour_gf
  public :: read_kb_contour_gf
  public :: plot_kb_contour_gf
  !
  public :: extrapolate_kb_contour_gf
  !

  !
  public :: kb_trapz
  public :: kb_half_trapz
  !
  public :: assignment(=)
  ! public :: operator(*)
  ! public :: operator(+)
  ! public :: operator(.add.)
  ! public :: operator(.star.)

  integer,save :: kb_contour_N=0
  integer,save :: kb_contour_L=0
  integer,save :: kb_contour_Niw=0



  integer :: Nlk,Nlat,Nspin,Norb
  integer :: ilk,ilat,ispin,iorb
  integer :: jlk,jlat,jspin,jorb
  integer :: klk,klat,kspin,korb


contains





  !======= ALLOCATE/DEALLOCATE ======= 
  include "NEQ_CONTOUR_GF/neq_contour_gf_allocate.h90"


  !======= CHECK DIMENSION ======= 
  include "NEQ_CONTOUR_GF/neq_contour_gf_check.h90"


  !======= ADD !!C(t,t')=A(t,t') + B(t,t') with t,t'=0,t_max
  include "NEQ_CONTOUR_GF/neq_contour_gf_add.h90"

  !======= SUM !! C(t,t')=A(t,t') + B(t,t'), with t=t_max && t'=0,t_max
  include "NEQ_CONTOUR_GF/neq_contour_gf_sum.h90"


  !======= DEL !! equals G(t,t') to zero along the actual perimeter
  include "NEQ_CONTOUR_GF/neq_contour_gf_del.h90"


  !======= CONVOLUTION !!C(t,t')=(A*B)(t,t'), with t=t_max && t'=0,t_max
  include "NEQ_CONTOUR_GF/neq_contour_gf_convolute.h90"  


  !======= VOLTERRA INTEGRAL EQUATION  !!solve G(t,t')+(K*G)(t,t')=Q(t,t') along the contour edge
  include "NEQ_CONTOUR_GF/neq_contour_gf_vie.h90" 


  !======= VOLTERRA INTEGRO-DIFFERENTIAL EQUATION !! solve [i*d/dt-h(t)]G(t,t') = delta(t,t') + (K*G)(t,t') along the contour edge
  include "NEQ_CONTOUR_GF/neq_contour_gf_vide.h90"


  !======= MISCELLANEOUS !!save,read,plot kb contour gf
  include "NEQ_CONTOUR_GF/neq_contour_gf_misc.h90"  


  !======= EXTRAPOLATION !!extrapolate a KB contour gf from a given time to the next one:
  include "NEQ_CONTOUR_GF/neq_contour_gf_extrapolate.h90"


  !======= OPERATIONS + & * ======= 
  include "NEQ_CONTOUR_GF/neq_contour_gf_operations.h90"




  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################







  subroutine set_kb_contour_gf_time(params)
    type(kb_contour_params) :: params
    kb_contour_N   = params%Nt
    kb_contour_L   = params%Ntau
    kb_contour_Niw = params%Niw
  end subroutine set_kb_contour_gf_time





  ! +get_Adv
  ! +get_Gtr
  ! +get_Rmix
  ! +get_Bar
  !+-----------------------------------------------------------------------------+!
  !PURPOSE: define some special operations on the GF:
  ! + get_adv : obtain the Adv component from G
  ! + get_rmix: obtain the RightMix component from G
  ! + get_gtr : obtain the Gtr component from G
  ! + get_bar : obtain the conjugate component in the Nambu space from G (dubbed bar)
  !+-----------------------------------------------------------------------------+!
  function get_adv(G,i,j) result (adv)
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



















  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################




  !----------------------------------------------------------------------------
  !  This function calculates the sum
  !    \sum_{k=i}^{j} w_{k}^{i,j} f_{k}.
  !    w_{k}^{i,j} = 1/2  for k=i,j
  !                = 1    for i<k<j
  !                      OR
  !                = 1    for i<k<=j (half-edge)
  !----------------------------------------------------------------------------
  function kb_trapz_d(f,ia,ib) result(sum)
    real(8),dimension(0:),intent(in) :: f
    integer,intent(in)              :: ia, ib
    integer                         :: k
    real(8)                         :: sum
    sum=0.d0
    if(ia==ib)then
       return
    else
       sum=sum+0.5d0*f(ia)
       do k=ia+1,ib-1
          sum=sum+f(k)
       end do
       sum=sum+0.5d0*f(ib)
    end if
  end function kb_trapz_d

  function kb_trapz_c(f,ia,ib,origin) result(sum)
    complex(8),dimension(0:),intent(in) :: f
    integer,intent(in)                 :: ia, ib
    integer                            :: k
    complex(8)                         :: sum
    logical,optional,intent(in)        :: origin
    logical                            :: w0
    real(8)                            :: w

    w0=.true.
    if(present(origin)) w0=origin

    w = 0.5d0
    if(.not.w0) w = 1d0

    sum=zero
    if(ia==ib)then
       return
    else
       sum=sum+w*f(ia)
       do k=ia+1,ib-1
          sum=sum+f(k)
       end do
       sum=sum+w*f(ib)
    end if
  end function kb_trapz_c

  function kb_half_trapz_d(f,ia,ib) result(sum)
    real(8),dimension(0:),intent(in) :: f
    integer,intent(in)              :: ia, ib
    integer                         :: k
    real(8)                         :: sum
    sum=0.5d0*f(ia)
    do k=ia+1,ib
       sum=sum+f(k)
    end do
  end function kb_half_trapz_d
  !
  function kb_half_trapz_c(f,ia,ib) result(sum)
    complex(8),dimension(0:),intent(in) :: f
    integer,intent(in)                 :: ia, ib
    integer                            :: k
    complex(8)                         :: sum
    sum=0.5d0*f(ia)
    do k=ia+1,ib
       sum=sum+f(k)
    end do
  end function kb_half_trapz_c


END MODULE NEQ_CONTOUR_GF






! !======= GET OTHER COMPONENTS ======= 
! function gtr_kb_contour_gf(G,N) result(Ggtr)
!   type(kb_contour_gf)       :: G
!   integer                   :: N
!   complex(8),dimension(N,N) :: Ggtr
!   integer                   :: i,j
!   Ggtr = zero
!   forall(i=1:N,j=1:N,i>=j)Ggtr(i,j) = G%less(i,j) + G%ret(i,j)
!   forall(i=1:N,j=1:N,i<j)Ggtr(i,j)=-conjg(Ggtr(j,i))
! end function gtr_kb_contour_gf

! function rmix_kb_contour_gf(G,N,L) result(Grmix)
!   type(kb_contour_gf)         :: G
!   integer                     :: N,L
!   complex(8),dimension(0:L,N) :: Grmix
!   integer                     :: itau,j
!   Grmix = zero
!   forall(itau=0:L,j=1:N)Grmix(itau,j) = conjg(G%lmix(j,L-itau))
! end function rmix_kb_contour_gf
