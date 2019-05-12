MODULE SETUP
  USE COMMON_VARS
  implicit none
  private
  !
  !
  !< build and delete sector (i.e. construct or erase the S-->F maps)
  public :: build_sector
  public :: delete_sector
  !
  !< sector dimensions and state decomposition I=Iup + Idw*2^Ns
  public :: get_sector_dimension
  public :: iup_index
  public :: idw_index
  !
  !< Binomial and binary_decomposition of an integer
  public :: binomial    
  public :: bdecomp
  !
  !< destruction and creation operators in 2nd quantization
  public :: c,cdg
  !
  !< binary search   
  public :: binary_search




contains



  !##################################################################
  !##################################################################
  !AUXILIARY PROCEDURES - Sectors,Nup,Ndw,DimUp,DimDw,...
  !##################################################################
  !##################################################################
  !
  !< Get sub-sector up or dw dimension (binomial)
  function get_sector_dimension(n,np) result(dim)
    integer :: n,np
    integer :: dim
    dim = binomial(n,np)
  end function get_sector_dimension
  !
  !< return the UP index of an state/integer I = I_up + Idw*2^Ns 
  function iup_index(i,DimUp) result(iup)
    integer :: i
    integer :: DimUp
    integer :: iup
    iup = mod(i,DimUp);if(iup==0)iup=DimUp
  end function iup_index
  !
  !< return the DW index of an state/integer I = I_up + Idw*2^Ns 
  function idw_index(i,DimUp) result(idw)
    integer :: i
    integer :: DimUp
    integer :: idw
    idw = (i-1)/DimUp+1
  end function idw_index






  !##################################################################
  !##################################################################
  !BUILD SECTORS
  !##################################################################
  !##################################################################
  subroutine build_sector(Nups,Ndws,H)
    integer                       :: isector
    type(sector_map),dimension(2) :: H
    integer                       :: Nups,Ndws
    integer                       :: DimUps,DimDws
    integer                       :: iup,idw
    integer                       :: nup_,ndw_
    integer                       :: dim,iud
    !
    !
    DimUps = get_sector_dimension(Ns,Nups)
    DimDws = get_sector_dimension(Ns,Ndws)
    !
    call map_allocate(H,[DimUps,DimDws])
    !UP    
    dim=0
    do iup=0,2**Ns-1
       nup_ = popcnt(iup)
       if(nup_ /= Nups)cycle
       dim  = dim+1
       H(1)%map(dim) = iup
    enddo
    !DW
    dim=0
    do idw=0,2**Ns-1
       ndw_= popcnt(idw)
       if(ndw_ /= Ndws)cycle
       dim = dim+1
       H(2)%map(dim) = idw
    enddo
    !
  end subroutine build_sector

  subroutine delete_sector(isector,H)
    integer                   :: isector
    type(sector_map)          :: H(:)
    call map_deallocate(H)
  end subroutine delete_sector




  !##################################################################
  !##################################################################
  !CREATION / DESTRUCTION OPERATORS
  !##################################################################
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates 
  !   |out>=C_pos|in>  OR  |out>=C^+_pos|in> ; 
  !   the sign of |out> has the phase convention, pos labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c

  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
  end subroutine cdg






  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp


  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  function binomial(n1,n2) result(nchoos)
    integer :: n1,n2
    real(8) :: xh
    integer :: i
    integer :: nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial



  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a value in an array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       !stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search











end module SETUP

