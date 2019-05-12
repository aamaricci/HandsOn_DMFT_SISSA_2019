MODULE AUX_FUNX
  USE COMMON_VARS
  implicit none

  interface print_conf
     module procedure :: print_conf_int
     module procedure :: print_conf_ivec
  end interface print_conf

contains



  !+------------------------------------------------------------------+
  !PURPOSE  : constructs the sectors by storing the map to the 
  !states i\in Hilbert_space from the states count in H_sector.
  !|ImpUP,BathUP>|ImpDW,BathDW >
  !+------------------------------------------------------------------+
  subroutine build_sector_normal(nup,ndw,Hup,Hdw)
    integer          :: nup,ndw
    type(sector_map) :: Hup
    type(sector_map) :: Hdw
    integer          :: i,iup,idw,dim
    integer          :: nup_,ndw_
    integer          :: ivec(Ns),jvec(Ns)
    !
    !< Workout UP sub-sector:
    !< allocate UP map to binomial(Ns nup) = dim_sub_sector(nup)
    call map_allocate(Hup,dim_sector_normal(nup))
    dim=0
    do i=0,2**Ns-1
       nup_ = popcnt(i)         !equivalent to sum(binary_decomposition(i))
       if(nup_ /= nup)cycle     !if this state does not have the required number of UP skip
       dim       = dim+1
       Hup%map(dim) = i
    enddo
    !
    !< Workout DW sub-sector:
    !< allocate DW map to binomial(Ns ndw) = dim_sub_sector(ndw)
    call map_allocate(Hdw,dim_sector_normal(ndw))
    dim=0
    do i=0,2**Ns-1
       ndw_=popcnt(i)
       if(ndw_ /= ndw)cycle
       dim       = dim+1
       Hdw%map(dim) = i
    enddo
    !
    !
    !< USELESS: PRETTY PRINT THE STATES OF THE SECTORS:
    do idw=0,2**Ns-1
       jvec  = bdecomp(idw,Ns)
       ndw_  = sum(jvec)
       if(ndw_ /= ndw)cycle
       do iup=0,2**Ns-1
          ivec  = bdecomp(iup,Ns)
          nup_  = sum(ivec)
          if(nup_ /= nup)cycle
          call print_conf_ud(ivec,jvec)
       enddo
    enddo
    write(*,*)""
    !
  end subroutine build_sector_normal
  !



  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Ntot)
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









  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine a(pos,in,out,fsgn)
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
  end subroutine a



  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm+|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine adg(pos,in,out,fsgn)
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
  end subroutine adg





  !+------------------------------------------------------------------+
  !PURPOSE  : return the dimension of a sector
  !+------------------------------------------------------------------+
  function dim_sector_normal(n) result(dim)
    integer          :: n
    integer          :: dim
    dim = binomial(Ns,n)    !this ensures better evaluation of the dimension
  end function dim_sector_normal



  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a value in an array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
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

  









  !+------------------------------------------------------------------+
  !PURPOSE  : print configuration
  !+------------------------------------------------------------------+
  subroutine print_conf_int(i,Ntot)
    integer           :: dim,i,j,unit_,Ntot
    integer           :: ivec(Ntot)
    character(len=2)  :: fbt
    character(len=16) :: fmt
    unit_=6
    write(fbt,'(I2.2)')Ntot
    fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    ivec = bdecomp(i,Ntot)
    i= bjoin(ivec,Ntot)    
    write(unit_,"(I9,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="no")"> - "
    write(unit_,fmt,advance="yes")i
  end subroutine print_conf_int

  subroutine print_conf_ivec(ivec)
    integer           :: dim,i,j,unit_,Ntot
    integer           :: ivec(:)
    character(len=2)  :: fbt
    character(len=16) :: fmt
    unit_=6
    Ntot = size(ivec)
    write(fbt,'(I2.2)')Ntot
    fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    i= bjoin(ivec,Ntot)    
    write(unit_,"(I9,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="no")"> - "
    write(unit_,fmt,advance="yes")i
  end subroutine print_conf_ivec

  subroutine print_conf_ud(ivec,jvec)
    integer,intent(in) :: ivec(:),jvec(size(ivec))
    integer            :: unit_
    integer            :: i,j,iup,idw,Ntot
    character(len=2)   :: fbt
    character(len=20)  :: fmt
    unit_=6
    Ntot = size(ivec)
    write(fbt,'(I2.2)')Ntot
    ! fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//",1x,B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    iup = bjoin(ivec,Ntot)
    idw = bjoin(jvec,Ntot)
    i = bjoin([ivec,jvec],2*Ntot)
    write(unit_,"(I9,1x,I4,1x,A1)",advance="no")i,iup,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A1,I4,A2)",advance="no")">",idw," |"
    write(unit_,"(10I1)",advance="no")(jvec(j),j=1,Ntot)
    write(unit_,"(A1)",advance="yes")">"
    ! write(unit_,fmt,advance="yes")ibits(i,0,Ntot),ibits(i,Ntot,2*Ntot)
  end subroutine print_conf_ud



  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Ntot) with the binary sequence 
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bjoin(ib,Ntot) result(i)
    integer                 :: Ntot
    integer,dimension(Ntot) :: ib
    integer                 :: i,j
    i=0
    do j=0,Ntot-1
       i=i+ib(j+1)*2**j
    enddo
  end function bjoin




  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor
  !+------------------------------------------------------------------+
  function binomial(n1,n2) result(nchoos)
    real(8) :: xh
    integer :: n1,n2,i
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


END MODULE AUX_FUNX




! function Imap_ud(iup,idw,Ns) result(i)
!   integer :: iup,idw,Ns
!   integer :: i
!   i   = iup + idw*2**Ns
! end function Imap_ud
