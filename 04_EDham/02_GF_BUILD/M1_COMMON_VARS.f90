!< This module contains a number of global variables and procedures used elsewhere:
! - eigh : interface to Lapack diagonalization routines
! - str  : a smart string handler
! - eye  : diagonal unit matrix
! - free_unit: return the first non-open unit
MODULE COMMON_VARS
  implicit none
  !
  !< the module is declared PRIVATE: only public declared variable exit from this module.
  !if no public variables or procedures are defined the module becomes almost useless.
  private
  !
  !< relevant public variables shared across the codes.
  integer,public :: Dim                !Sector total dimension
  integer,public :: DimUp              !sub-sector up dimension
  integer,public :: DimDw              !sub-sector dw dimnsion
  integer,public :: Nup                !Number of up spin
  integer,public :: Ndw                !Number of dw spin
  integer,public :: Ns                 !Number of levels per spin
  !
  !
  !fortran-esque: define a unique interface for different, similar procedure.
  ! eigh interface either diagonalization (dsyev) or tri-diagonalization from Lapack.
  ! The choice is automatically done by the compiler at run time depending on the inputs
  ! (different among the procedures and checked at compilation time)
  interface eigh
     module procedure deigh_simple
     module procedure deigh_tridiag
  end interface eigh
  !
  !< str interface to different string manipulation variable-to-string
  interface str
     module procedure str_i_to_ch
     module procedure str_i_to_ch_pad
     module procedure str_r_to_ch
     module procedure str_c_to_ch
     module procedure str_l_to_ch
     module procedure str_ch_to_ch
  end interface str
  !
  !ED RELATED TYPES & PROCEDURES
  !< sector-to-fock space maps
  type sector_map
     integer,dimension(:),allocatable :: map
     logical                          :: status=.false.
  end type sector_map
  !
  !< interface to allocate maps
  interface map_allocate
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate
  !
  !< interface to deallocate maps
  interface map_deallocate
     module procedure :: map_deallocate_scalar
     module procedure :: map_deallocate_vector
  end interface map_deallocate
  !
  !< PUBLIC all procedures
  public :: eigh
  public :: eye
  public :: kronecker_product
  public :: str
  public :: free_unit
  public :: linspace
  public :: arange
  !

  public :: sector_map
  public :: map_allocate
  public :: map_deallocate



contains




  subroutine map_allocate_scalar(H,N)
    type(sector_map) :: H
    integer          :: N
    if(H%status) call map_deallocate_scalar(H)
    allocate(H%map(N))
    H%status=.true.
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:)       :: H
    integer,dimension(size(H))          :: N
    integer                             :: i
    do i=1,size(H)
       call map_allocate_scalar(H(i),N(i))
    enddo
  end subroutine map_allocate_vector

  subroutine map_deallocate_scalar(H)
    type(sector_map) :: H
    if(.not.H%status)then
       write(*,*) "WARNING map_deallocate_scalar: H is not allocated"
       return
    endif
    if(allocated(H%map))deallocate(H%map)
    H%status=.false.
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:) :: H
    integer                       :: i
    do i=1,size(H)
       call map_deallocate_scalar(H(i))
    enddo
  end subroutine map_deallocate_vector



  function eye(n) result(A)
    integer, intent(in) :: n
    real(8)             :: A(n, n)
    integer             :: i
    A = 0d0
    do i = 1, n
       A(i,i) = 1d0
    end do
  end function eye




  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function linspace(start,stop,num,istart,iend,mesh) result(array)
    real(8)          :: start,stop,step,array(num)
    integer          :: num,i
    logical,optional :: istart,iend
    logical          :: startpoint_,endpoint_
    real(8),optional :: mesh
    if(num<0)stop "linspace: N<0, abort."
    startpoint_=.true.;if(present(istart))startpoint_=istart
    endpoint_=.true.;if(present(iend))endpoint_=iend
    if(startpoint_.AND.endpoint_)then
       if(num<2)stop "linspace: N<2 with both start and end points"
       step = (stop-start)/real(num-1,8)
       forall(i=1:num)array(i)=start + real(i-1,8)*step
    elseif(startpoint_.AND.(.not.endpoint_))then
       step = (stop-start)/real(num,8)
       forall(i=1:num)array(i)=start + real(i-1,8)*step
    elseif(.not.startpoint_.AND.endpoint_)then
       step = (stop-start)/real(num,8)
       forall(i=1:num)array(i)=start + real(i,8)*step
    else
       step = (stop-start)/real(num+1,8)
       forall(i=1:num)array(i)=start + real(i,8)*step
    endif
    if(present(mesh))mesh=step
  end function linspace




  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function arange(start,num,iend) result(array)
    integer          :: start,array(num)
    integer          :: num,i
    logical,optional :: iend
    logical          :: endpoint_
    if(num<0)stop "arange: N<0, abort."
    endpoint_=.true.;if(present(iend))endpoint_=iend
    if(endpoint_)then
       forall(i=1:num)array(i)=start+i-1
    else
       forall(i=1:num-1)array(i)=start+i-1
    end if
  end function arange


  !-------------------------------------------------------------------------------------------
  !PURPOSE:  eigenvalue/-vector problem for real symmetric/complex hermitian matrices:
  !-------------------------------------------------------------------------------------------
  subroutine deigh_simple(M,E,jobz,uplo)
    real(8),dimension(:,:),intent(inout) :: M ! M v = E v/v(i,j) = ith component of jth vec.
    real(8),dimension(:),intent(inout)   :: E ! eigenvalues
    character(len=1),optional            :: jobz,uplo
    character(len=1)                     :: jobz_,uplo_
    integer                              :: i,j,n,lda,info,lwork
    real(8),dimension(:),allocatable     :: work
    real(8),dimension(1)                 :: lwork_guess
    jobz_='V';if(present(jobz))jobz_=jobz
    uplo_='U';if(present(uplo))uplo_=uplo
    lda = max(1,size(M,1))
    n   = size(M,2)
    call assert_shape(M,[n,n],"eigh","M")
    Call dsyev(jobz_,uplo_,n,M,lda,E,lwork_guess,-1,info)
    if (info /= 0) then
       print*, "dsyevd returned info =", info
       if (info < 0) then
          print*, "the", -info, "-th argument had an illegal value"
       else
          print*, "the algorithm failed to compute an eigenvalue while working"
          print*, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
          print*, "through", mod(info, n+1)
       end if
       stop 'error deigh: 1st call dsyev'
    end if
    lwork=lwork_guess(1)
    allocate(work(lwork))
    call dsyev(jobz_,uplo_,n,M,lda,E,work,lwork,info)
    if (info /= 0) then
       print*, "dsyevd returned info =", info
       if (info < 0) then
          print*, "the", -info, "-th argument had an illegal value"
       else
          print*, "the algorithm failed to compute an eigenvalue while working"
          print*, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
          print*, "through", mod(info, n+1)
       end if
       stop 'error deigh: 2nd call dsyev'
    end if
    deallocate(work)
  end subroutine deigh_simple




  subroutine deigh_tridiag(D,U,Ev,Irange,Vrange)
    real(8),dimension(:)                :: d
    real(8),dimension(max(1,size(d)-1)) :: u
    real(8),dimension(:,:),optional     :: Ev
    integer,dimension(2),optional       :: Irange
    integer,dimension(2),optional       :: Vrange
    !
    integer                             :: n
    ! = 'N':  Compute eigenvalues only;
    ! = 'V':  Compute eigenvalues and eigenvectors.
    character                           :: jobz_
    ! = 'A': all eigenvalues will be found
    ! = 'V': all eigenvalues in the half-open interval (VL,VU] will be found.
    ! = 'I': the IL-th through IU-th eigenvalues will be found.
    ! For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and
    ! DSTEIN are called
    character                           :: range_
    real(8)                             :: vl,vu 
    integer                             :: il,iu
    real(8),dimension(size(d))          :: w
    real(8),dimension(:,:),allocatable  :: z
    integer,dimension(:),allocatable    :: isuppz
    integer                             :: lwork,liwork
    real(8),dimension(:),allocatable    :: work
    integer,dimension(:),allocatable    :: iwork
    real(8)                             :: abstol = 0d0 !tolerance on approximation error of eigenvalues
    integer                             :: m
    integer                             :: ldz
    integer                             :: info
    !
    n = size(d)
    !
    jobz_ = 'N'; if(present(Ev))jobz_='V'
    !
    range_= 'A'
    m     = n
    if(present(Irange).AND.present(Vrange))stop "EIGH_T: Irange & Vrange both present"
    if(present(Irange))then
       range_ = 'I'
       il = Irange(1)
       iu = Irange(2)
       m  = iu-il+1
    endif
    if(present(Vrange))then
       range_ = 'V'
       vl     = Vrange(1)
       vu     = Vrange(2)
       m      = n
    endif
    !
    if(present(Ev))then
       if(any(shape(Ev)/=[n,m]))stop "EIGH_T ERROR: wrong dimension in Ev"
    endif
    !
    ldz = n
    !
    allocate(z(n,m), isuppz(2*max(1,m)))
    allocate(work(1),iwork(1))
    call dstevr(jobz_,range_,n,d,u,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,-1,iwork,-1,info)
    if(info/=0) then
       print*, "dstevr returned info =", info
       if (info < 0) then
          print*, "the", -info, "-th argument had an illegal value"
       else
          print*, "Internal error"
       end if
       stop 'error dstevr: 1st call dstevr'
    end if
    !
    lwork  = work(1)
    liwork = iwork(1)
    deallocate(work,iwork);allocate(work(lwork),iwork(liwork))
    !
    call dstevr(jobz_,range_,n,d,u,vl,vu,il,iu,abstol,m,w,z,ldz,isuppz,work,lwork,iwork,liwork,info)
    if(info/=0) then
       print*, "dstevr returned info =", info
       if (info < 0) then
          print*, "the", -info, "-th argument had an illegal value"
       else
          print*, "Internal error"
       end if
       stop 'error dstevr: 2nd call dstevr'
    end if
    !
    d  = w
    if(present(Ev)) Ev = z
    !
    deallocate(work,iwork,isuppz,z)
    !
  end subroutine deigh_tridiag



  subroutine assert_shape(A, shap, routine, matname)
    real(8),intent(in)  :: A(:,:)
    integer,intent(in)  :: shap(:)
    character(len=*)    :: routine, matname
    if(any(shape(A) /= shap)) then
       print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
       print*, "Shape should be ", shap
       stop "Aborting due to illegal matrix operation"
    end if
  end subroutine assert_shape


  function free_unit(n) result(unit_)
    integer,optional :: n
    integer          :: unit_,ios
    logical          :: opened
    unit_=100
    do 
       unit_=unit_+1
       INQUIRE(unit=unit_,OPENED=opened,iostat=ios)
       if(.not.opened.AND.ios==0)exit 
       if(unit_>900) stop "ERROR free_unit: no unit free smaller than 900. Possible BUG"
    enddo
    if(present(n))n=unit_
  end function free_unit








  function kronecker_product(A,B) result(AxB)
    real(8),intent(in) :: A(:,:), B(:,:)
    integer            :: i,j
    integer            :: rowA,colA
    integer            :: rowB,colB
    real(8)            :: AxB(size(A,1)*size(B,1),size(A,2)*size(B,2))
    AxB = 0
    rowA=size(A,1) ; colA=size(A,2)
    rowB=size(B,1) ; colB=size(B,2)
    forall(i=1:rowA,j=1:colA)
       AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
    end forall
  end function kronecker_product







  function str_i_to_ch(i4) result(string)
    integer                      :: i4
    character(len=:),allocatable :: string
    character(len=16)            :: string_
    call i4_to_s_left(i4,string_)
    string=trim(adjustl(trim(string_)))
  end function str_i_to_ch

  function str_i_to_ch_pad(i4,Npad) result(string)
    integer                      :: i4
    integer                      :: Npad
    character(len=:),allocatable :: string
    character(len=Npad)          :: string_pad
    call i4_to_s_zero(i4,string_pad)
    string=trim(adjustl(trim(string_pad)))
  end function str_i_to_ch_pad

  function str_r_to_ch(r8,d) result(string)
    real(8)                      :: r8
    integer,optional             :: d
    integer                      :: w_,d_
    character(len=:),allocatable :: string
    character(len=:),allocatable :: string_
    d_=6 ;if(present(d))d_=d
    w_ = get_w_(r8,d_)
    allocate(character(len=w_) :: string_)
    call r8_to_s_left(r8,string_,d_)
    string=trim(adjustl(trim(string_)))
  end function str_r_to_ch


  function str_c_to_ch(c,d) result(string)
    complex(8)                   :: c
    integer,optional             :: d
    integer                      :: w_,d_
    character(len=:),allocatable :: string
    character(len=:),allocatable :: sre,sim
    real(8)                      :: re,im
    d_=6 ;if(present(d))d_=d
    re=dreal(c)
    w_ = get_w_(re,d_)
    allocate(character(len=w_) :: sre)
    call r8_to_s_left(re,sre,d_)
    !
    im=dimag(c)
    w_ = get_w_(im,d_)
    allocate(character(len=w_) :: sim)
    call r8_to_s_left(im,sim,d_)
    string="("//trim(adjustl(trim(sre)))//","//trim(adjustl(trim(sim)))//")"
  end function str_c_to_ch


  function str_l_to_ch(bool) result(string)
    logical          :: bool
    character(len=1) :: string
    string='F'
    if(bool)string='T'
  end function str_l_to_ch

  function str_ch_to_ch(txt) result(string)
    character(len=*)                             :: txt
    character(len=:),allocatable :: string
    string=trim(adjustl(trim(txt)))
  end function str_ch_to_ch


  function get_w_(r8,d) result(w)
    real(8) :: r8
    integer :: d
    integer :: w
    if(r8==0d0)then
       w=d+4
    else
       w=floor(log10(abs(r8)))
       if(w < -1)then
          w = d + 4 + 4
       else
          w = w + d + 4
       endif
    endif
  end function get_w_







  subroutine i4_to_s_left ( i4, s )
    !! I4_TO_S_LEFT converts an I4 to a left-justified string.
    !  Example:
    !    Assume that S is 6 characters long:
    !        I4  S
    !         1  1
    !        -1  -1
    !         0  0
    !      1952  1952
    !    123456  123456
    !   1234567  ******  <-- Not enough room!
    !  Parameters:
    !    Input, integer ( kind = 4 ) I4, an integer to be converted.
    !    Output, character ( len = * ) S, the representation of the integer.
    !    The integer will be left-justified.  If there is not enough space,
    !    the string will be filled with stars.
    character :: c
    integer   :: i
    integer   :: i4
    integer   :: idig
    integer   :: ihi
    integer   :: ilo
    integer   :: ipos
    integer   :: ival
    character(len=*) ::  s
    s = ' '
    ilo = 1
    ihi = len ( s )
    if ( ihi <= 0 ) then
       return
    end if
    !  Make a copy of the integer.
    ival = i4
    !  Handle the negative sign.
    if ( ival < 0 ) then
       if ( ihi <= 1 ) then
          s(1:1) = '*'
          return
       end if
       ival = -ival
       s(1:1) = '-'
       ilo = 2
    end if
    !  The absolute value of the integer goes into S(ILO:IHI).
    ipos = ihi
    !  Find the last digit of IVAL, strip it off, and stick it into the string.
    do
       idig = mod ( ival, 10 )
       ival = ival / 10
       if ( ipos < ilo ) then
          do i = 1, ihi
             s(i:i) = '*'
          end do
          return
       end if
       call digit_to_ch ( idig, c )
       s(ipos:ipos) = c
       ipos = ipos - 1
       if ( ival == 0 ) then
          exit
       end if
    end do
    !  Shift the string to the left.
    s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
    s(ilo+ihi-ipos:ihi) = ' '
  end subroutine i4_to_s_left


  subroutine r8_to_s_left( r8, s, digits)
    !! R8_TO_S_LEFT writes an R8 into a left justified string.
    !    An R8 is a real ( kind = 8 ) value.
    !    A 'F<len(s)>.DIGITS' format is used with a WRITE statement.
    character(len=12)   :: fmt
    integer             :: i
    real(8)             :: r8
    character(len=*)    :: s
    integer             :: s_length,w_
    integer             :: digits
    s_length = len ( s )
    write(fmt,"(A2,I0,A1,I0,A1)")"(F",s_length,".",digits,")"
    if(r8/=0d0)then
       w_=floor(log10(abs(r8)))
       if(w_<-1)write(fmt,"(A3,I0,A1,I0,A1)")"(ES",s_length,".",digits,")"
    endif
    write ( s, fmt ) r8
    s = trim(adjustl(trim( s )))
  end subroutine r8_to_s_left


  subroutine digit_to_ch(digit,ch)
    !! DIGIT_TO_CH returns the character representation of a decimal digit.
    !    Instead of CHAR, we now use the ACHAR function, which
    !    guarantees the ASCII collating sequence.
    !  Example:
    !    DIGIT   CH
    !    -----  ---
    !      0    '0'
    !      1    '1'
    !    ...    ...
    !      9    '9'
    !     17    '*'
    !  Parameters:
    !    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
    !    Output, character CH, the corresponding character.
    character :: ch
    integer   :: digit
    if ( 0 <= digit .and. digit <= 9 ) then
       ch = achar ( digit + 48 )
    else
       ch = '*'
    end if
  end subroutine digit_to_ch


  subroutine i4_to_s_zero ( intval, s )
    !! I4_TO_S_ZERO converts an I4 to a string, with zero padding.
    !    An I4 is an integer ( kind = 4 ).
    !  Example:
    !    Assume that S is 6 characters long:
    !    INTVAL  S
    !         1  000001
    !        -1  -00001
    !         0  000000
    !      1952  001952
    !    123456  123456
    !   1234567  ******  <-- Not enough room!
    !  Parameters:
    !    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
    !    Output, character ( len = * ) S, the representation of the integer.
    !    The integer will be right justified, and zero padded.
    !    If there is not enough space, the string will be filled with stars.
    implicit none
    character c
    integer ( kind = 4 ) i
    integer ( kind = 4 ) idig
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) intval
    integer ( kind = 4 ) ipos
    integer ( kind = 4 ) ival
    character ( len = * ) s
    s = ' '
    ilo = 1
    ihi = len ( s )
    if ( ihi <= 0 ) then
       return
    end if
    !
    !  Make a copy of the integer.
    !
    ival = intval
    !
    !  Handle the negative sign.
    !
    if ( ival < 0 ) then
       if ( ihi <= 1 ) then
          s(1:1) = '*'
          return
       end if
       ival = -ival
       s(1:1) = '-'
       ilo = 2
    end if
    !
    !  Working from right to left, strip off the digits of the integer
    !  and place them into S(ILO:IHI).
    !
    ipos = ihi
    do while ( ival /= 0 .or. ipos == ihi )
       idig = mod ( ival, 10 )
       ival = ival / 10
       if ( ipos < ilo ) then
          do i = 1, ihi
             s(i:i) = '*'
          end do
          return
       end if
       call digit_to_ch ( idig, c )
       s(ipos:ipos) = c
       ipos = ipos - 1
    end do
    !
    !  Fill the empties with zeroes.
    !
    do i = ilo, ipos
       s(i:i) = '0'
    end do
    return
  end subroutine i4_to_s_zero







end module COMMON_VARS

