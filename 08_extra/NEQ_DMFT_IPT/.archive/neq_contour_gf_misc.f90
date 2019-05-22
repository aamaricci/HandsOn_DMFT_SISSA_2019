! CHECK DIMENSION
function check_dimension_kb_contour_gf(G,params) result(bool)
  type(kb_contour_gf)     :: G
  type(kb_contour_params) :: params
  logical                 :: bool
  integer                 :: N,L
  bool=.false.
  N=params%Ntime              !<== check size at maximum time
  L=params%Ntau
  if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) )&
       stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions less/ret"
  if( size(G%lmix)/=N*(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions lmix"
  if( size(G%mats)/=(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions mats"
  bool=.true.
end function check_dimension_kb_contour_gf
!
function check_dimension_kb_contour_sigma(S,params) result(bool)
  type(kb_contour_sigma)  :: S
  type(kb_contour_params) :: params
  logical                 :: bool
  integer                 :: N,L
  bool=.false.
  N=params%Ntime              !<== check size at maximum time
  L=params%Ntau
  if( (size(S%reg%less)/=N**2) .OR. (size(S%reg%ret)/=N**2) )&
       stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions less/ret"
  if( size(S%reg%lmix)/=N*(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions lmix"
  if( size(S%reg%mats)/=(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions mats"
  if( size(S%hf)/=N )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions HF"
  bool=.true.
end function check_dimension_kb_contour_sigma
!
function check_dimension_kb_contour_gf_(G,N,L) result(bool)
  type(kb_contour_gf)     :: G
  logical                 :: bool
  integer                 :: N,L
  bool=.false.
  if( (size(G%less)/=N**2) .OR. (size(G%ret)/=N**2) )&
       stop "ERROR contour_gf/check_dimension_kb_contour_gf: wrong dimensions less/ret"
  if( size(G%lmix)/=N*(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions lmix"
  if( size(G%mats)/=(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions mats"
  bool=.true.
end function check_dimension_kb_contour_gf_
!
function check_dimension_kb_contour_sigma_(S,N,L) result(bool)
  type(kb_contour_sigma)  :: S
  logical                 :: bool
  integer                 :: N,L
  bool=.false.
  if( (size(S%reg%less)/=N**2) .OR. (size(S%reg%ret)/=N**2) )&
       stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions less/ret"
  if( size(S%reg%lmix)/=N*(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions lmix"
  if( size(S%reg%mats)/=(L+1) )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions mats"
  if( size(S%hf)/=N )stop "contour_gf/check_dimension_kb_contour_gf: wrong dimensions HF"
  bool=.true.
end function check_dimension_kb_contour_sigma_
!
function check_dimension_kb_contour_dgf(dG,params) result(bool)
  type(kb_contour_dgf)    :: dG
  type(kb_contour_params) :: params
  logical                 :: bool
  integer                 :: N,L
  bool=.false.
  N=params%Ntime              !<== check size at maximum time
  L=params%Ntau
  if( (size(dG%less)/=N) .OR. (size(dG%ret)/=N) )&
       stop "ERROR contour_gf/check_dimension_kb_contour_dgf: wrong dimensions less/ret"
  if( size(dG%lmix)/=(L+1) )&
       stop "ERROR contour_gf/check_dimension_kb_contour_dgf: wrong dimensions lmix"
  bool=.true.
end function check_dimension_kb_contour_dgf
!
function check_dimension_kb_contour_dgf_(dG,N,L) result(bool)
  type(kb_contour_dgf)    :: dG
  logical                 :: bool
  integer                 :: N,L
  bool=.false.
  if( (size(dG%less)/=N) .OR. (size(dG%ret)/=N) )&
       stop "ERROR contour_gf/check_dimension_kb_contour_dgf: wrong dimensions less/ret"
  if( size(dG%lmix)/=(L+1) )&
       stop "ERROR contour_gf/check_dimension_kb_contour_dgf: wrong dimensions lmix"
  bool=.true.
end function check_dimension_kb_contour_dgf_



! SAVE
subroutine save_kb_contour_gf(G,file)
  type(kb_contour_gf) :: G
  character(len=*)    :: file
  integer             :: unit,N,L,Lf
  if(.not.G%status)stop "contour_gf/save_kb_contour_gf: G not allocated"
  unit = free_unit()
  open(unit,file=reg(file)//"_dimension.data.neqipt")
  N=size(G%less,1)
  L=size(G%lmix,2)-1
  Lf=size(G%iw)
  write(unit,"(A1,3A4,3X)")"#","N","L","Lf"
  write(unit,"(1X,3I4,3X)")N,L,Lf
  close(unit)
  call save_array(reg(file)//"_less.data.neqipt",G%less(:,:))
  call save_array(reg(file)//"_ret.data.neqipt", G%ret(:,:))
  call save_array(reg(file)//"_lmix.data.neqipt",G%lmix(:,0:))
  call save_array(reg(file)//"_mats.data.neqipt",G%mats(0:))
  call save_array(reg(file)//"_iw.data.neqipt",G%iw(:))
end subroutine save_kb_contour_gf
!
subroutine save_kb_contour_sigma(S,file)
  type(kb_contour_sigma) :: S
  character(len=*)       :: file
  integer                :: unit,i
  if(.not.S%status)stop "contour_gf/save_kb_contour_gf: S not allocated"
  unit = free_unit()
  open(unit,file=reg(file)//"_hfb.data.neqipt")
  do i=1,size(S%hf(:))
     write(unit,*)S%hf(i)
  enddo
  close(unit)
  call save_kb_contour_gf(S%reg,file)
end subroutine save_kb_contour_sigma
!
subroutine save_kb_contour_dgf(dG,file)
  type(kb_contour_dgf) :: dG
  character(len=*)     :: file
  integer              :: unit,N,L
  if(.not.dG%status)stop "contour_gf/save_kb_contour_dgf: dG not allocated"
  unit = free_unit()
  N=size(dG%less)
  L=size(dG%lmix)-1
  open(unit,file=reg(file)//"_dimension.data.neqipt")
  write(unit,"(A1,2A4,2X)")"#","N","L"
  write(unit,"(1X,2I4,2X)")N,L
  close(unit)
  call save_array(reg(file)//"_less.data.neqipt",dG%less(:))
  call save_array(reg(file)//"_ret.data.neqipt", dG%ret(:))
  call save_array(reg(file)//"_lmix.data.neqipt",dG%lmix(0:))
end subroutine save_kb_contour_dgf





! READ
subroutine read_kb_contour_gf(G,file)
  type(kb_contour_gf)  :: G
  character(len=*)     :: file
  logical              :: check
  check = inquire_kb_contour_gf(file)
  if(.not.G%status.OR..not.check)stop "contour_gf/read_kb_contour_gf: G not allocated"
  call read_array(trim(file)//"_less.data.neqipt",G%less(:,:))
  call read_array(trim(file)//"_ret.data.neqipt",G%ret(:,:))
  call read_array(trim(file)//"_lmix.data.neqipt",G%lmix(:,0:))
  call read_array(trim(file)//"_mats.data.neqipt",G%mats(0:))
  call read_array(trim(file)//"_iw.data.neqipt",G%iw(:))
end subroutine read_kb_contour_gf
!
subroutine read_kb_contour_sigma(S,file)
  type(kb_contour_sigma) :: S
  character(len=*)       :: file
  logical                :: check
  integer                :: unit,i,Len
  check = inquire_kb_contour_sigma(file)
  if(.not.S%status.OR..not.check)stop "contour_gf/read_kb_contour_sigma: S not allocated"
  Len  = file_length(reg(file)//"_hfb.data.neqipt")
  if( Len > size(S%hf) ) stop "contour_gf/read_kb_contour_sigma: length(file) > size(S%hf)"
  unit = free_unit()
  open(unit,file=reg(file)//"_hfb.data.neqipt")
  do i=1,Len
     read(unit,*)S%hf(i)
  enddo
  close(unit)
  call read_kb_contour_gf(S%reg,file)
end subroutine read_kb_contour_sigma
!
subroutine read_kb_contour_dgf(dG,file)
  type(kb_contour_dgf) :: dG
  character(len=*)     :: file
  logical              :: check
  check = inquire_kb_contour_dgf(file)
  if(.not.dG%status.OR..not.check)stop "contour_gf/read_kb_contour_dgf: dG not allocated"
  call read_array(trim(file)//"_less.data.neqipt",dG%less(:))
  call read_array(trim(file)//"_ret.data.neqipt",dG%ret(:))
  call read_array(trim(file)//"_lmix.data.neqipt",dG%lmix(0:))
end subroutine read_kb_contour_dgf
!
!
!aux: inquire:
function inquire_kb_contour_gf(file) result(check)
  integer          :: i
  logical          :: check,bool(5)
  character(len=*) :: file
  character(len=16),dimension(5)  :: ctype=([ character(len=5) :: 'less','ret','lmix','mats','iw'])
  check=.true.
  do i=1,5
     inquire(file=reg(file)//"_"//reg(ctype(i))//".data.neqipt",exist=bool(i))
     if(.not.bool(i))inquire(file=reg(file)//"_"//reg(ctype(i))//".data.neqipt.gz",exist=bool(i))
     check=check.AND.bool(i)
  enddo
end function inquire_kb_contour_gf
!
function inquire_kb_contour_sigma(file) result(check)
  integer          :: i
  logical          :: check,bool(6)
  character(len=*) :: file
  character(len=16),dimension(5)  :: ctype=([ character(len=5) :: 'less','ret','lmix','mats','iw'])
  check=.true.
  do i=1,5
     inquire(file=reg(file)//"_"//reg(ctype(i))//".data.neqipt",exist=bool(i))
     if(.not.bool(i))inquire(file=reg(file)//"_"//reg(ctype(i))//".data.neqipt.gz",exist=bool(i))
     check=check.AND.bool(i)
  enddo
  inquire(file=reg(file)//"_hfb.data.neqipt",exist=bool(6))
  check=check.AND.bool(6)
end function inquire_kb_contour_sigma
!
function inquire_kb_contour_dgf(file) result(check)
  integer          :: i
  logical          :: check,bool(3)
  character(len=*) :: file
  character(len=16),dimension(3)  :: ctype=([character(len=5) :: 'less','ret','lmix'])
  check=.true.
  do i=1,3
     inquire(file=reg(file)//"_"//reg(ctype(i))//".data.neqipt",exist=bool(i))
     if(.not.bool(i))inquire(file=reg(file)//"_"//reg(ctype(i))//".data.neqipt.gz",exist=bool(i))
     check=check.AND.bool(i)
  enddo
end function inquire_kb_contour_dgf




! PLOT
subroutine plot_kb_contour_gf(file,G,params)
  character(len=*)        :: file
  type(kb_contour_gf)     :: G
  type(kb_contour_params) :: params
  integer                 :: Nt,i
  if(.not.G%status)stop "contour_gf/plot_kb_contour_gf: G is not allocated" 
  Nt=params%Ntime
  call splot3d(reg(file)//"_less_t_t.data.neqipt",params%t(:Nt),params%t(:Nt),G%less(:Nt,:Nt))
  call splot3d(reg(file)//"_ret_t_t.data.neqipt",params%t(:Nt),params%t(:Nt),G%ret(:Nt,:Nt))
  call splot3d(reg(file)//"_lmix_t_tau.data.neqipt",params%t(:Nt),params%tau(0:),G%lmix(:Nt,0:))
  call splot(reg(file)//"_mats_tau.data.neqipt",params%tau(0:),G%mats(0:))
  call splot(reg(file)//"_mats_iw.data.neqipt",params%wm(:),G%iw(:))
end subroutine plot_kb_contour_gf
!
subroutine plot_kb_contour_sigma(file,S,params)
  character(len=*)        :: file
  type(kb_contour_sigma)  :: S
  type(kb_contour_params) :: params
  integer                 :: Nt,i
  if(.not.S%status)stop "contour_gf/plot_kb_contour_sigma: S is not allocated" 
  Nt=params%Ntime
  call splot3d(reg(file)//"_less_t_t.data.neqipt",params%t(:Nt),params%t(:Nt),S%reg%less(:Nt,:Nt))
  call splot3d(reg(file)//"_ret_t_t.data.neqipt",params%t(:Nt),params%t(:Nt),S%reg%ret(:Nt,:Nt))
  call splot3d(reg(file)//"_lmix_t_tau.data.neqipt",params%t(:Nt),params%tau(0:),S%reg%lmix(:Nt,0:))
  call splot(reg(file)//"_mats_tau.data.neqipt",params%tau(0:),S%reg%mats(0:))
  call splot(reg(file)//"_mats_iw.data.neqipt",params%wm(:),S%reg%iw(:)+S%hf(1))
  call splot(reg(file)//"_hfb_t.data.neqipt",params%t(:Nt),S%hf(:Nt))
end subroutine plot_kb_contour_sigma
!
subroutine plot_kb_contour_dgf(file,dG,params)
  character(len=*)        :: file
  type(kb_contour_dgf)    :: dG
  type(kb_contour_params) :: params
  integer                 :: Nt
  if(.not.dG%status)stop "contour_gf/plot_kb_contour_gf: G is not allocated" 
  Nt=params%Ntime
  call splot(reg(file)//"_less_t.data.neqipt",params%t(:Nt),dG%less(:Nt))
  call splot(reg(file)//"_ret_t.data.neqipt",params%t(:Nt),dG%ret(:Nt))
  call splot(reg(file)//"_lmix_tau.data.neqipt",params%tau(0:),dG%lmix(0:))
end subroutine plot_kb_contour_dgf






! EXTRAPOLATION
!extrapolate a function from a given time to the next one:
subroutine extrapolate_kb_contour_gf(g,params)
  type(kb_contour_gf)     :: g
  type(kb_contour_params) :: params
  integer                 :: i,j,k,N,L
  if(.not.g%status)     stop "extrapolate_kb_contour_gf: g is not allocated"
  if(.not.params%status)stop "extrapolate_kb_contour_gf: params is not allocated"
  N = params%Nt
  L = params%Ntau
  select case(N)
  case(1)
     return
  case(2)
     !GUESS G AT THE NEXT STEP, GIVEN THE INITIAL CONDITIONS
     do j=1,N
        g%ret(N,j) =g%ret(1,1)
        g%less(N,j)=g%less(1,1)
     end do
     do i=1,N-1
        g%less(i,N)=g%less(1,1)
     end do
     do j=0,L
        g%lmix(N,j)=g%lmix(1,j)
     end do
  case default
     !EXTEND G FROM THE [N-1,N-1] TO THE [N,N] SQUARE TO START DMFT
     !USING QUADRATIC EXTRAPOLATION
     do k=1,N-1
        g%less(N,k)=2.d0*g%less(N-1,k)-g%less(N-2,k)
        g%less(k,N)=2.d0*g%less(k,N-1)-g%less(k,N-2)
     end do
     g%less(N,N)=2.d0*g%less(N-1,N-1)-g%less(N-2,N-2)
     !
     do k=0,L
        g%lmix(N,k)=2.d0*g%lmix(N-1,k)-g%lmix(N-2,k)
     end do
     !
     g%ret(N,N)=-xi
     do k=1,N-2
        g%ret(N,k)=2.d0*g%ret(N-1,k)-g%ret(N-2,k)
     end do
     g%ret(N,N-1)=0.5d0*(g%ret(N,N)+g%ret(N,N-2))
  end select
end subroutine extrapolate_kb_contour_gf

subroutine extrapolate_kb_contour_sigma(s,params)
  type(kb_contour_sigma)  :: s
  type(kb_contour_params) :: params
  integer                 :: i,j,k,N,L
  if(.not.s%status)     stop "extrapolate_kb_contour_gf: g is not allocated"
  if(.not.params%status)stop "extrapolate_kb_contour_gf: params is not allocated"
  N = params%Nt
  select case(N)
  case(1)
     return
  case(2)
     !GUESS Shf AT THE NEXT STEP, GIVEN THE INITIAL CONDITIONS
     s%hf(N)=s%hf(1)
  case default
     !EXTEND Shf FROM THE [N-1,N-1] TO THE [N,N] SQUARE TO START DMFT
     !USING QUADRATIC EXTRAPOLATION
     S%hf(N)=2.d0*S%hf(N-1)-S%hf(N-2)
  end select
  call extrapolate_kb_contour_gf(s%reg,params)
end subroutine extrapolate_kb_contour_sigma






