!======= SAVE
subroutine save_kb_contour_gf_main(G,file)
  type(kb_contour_gf) :: G
  character(len=*)    :: file
  integer             :: unit
  call save_array(reg(file)//"_less.data",G%less(:,:))
  call save_array(reg(file)//"_ret.data", G%ret(:,:))
  call save_array(reg(file)//"_lmix.data",G%lmix(:,0:))
  call save_array(reg(file)//"_mats.data",G%mats(0:))
  call save_array(reg(file)//"_iw.data",G%iw(:))
end subroutine save_kb_contour_gf_main
!
subroutine save_kb_contour_dgf_main(dG,file)
  type(kb_contour_dgf) :: dG
  character(len=*)     :: file
  integer              :: unit
  call save_array(reg(file)//"_less.data",dG%less(:))
  call save_array(reg(file)//"_ret.data", dG%ret(:))
  call save_array(reg(file)//"_lmix.data",dG%lmix(0:))
end subroutine save_kb_contour_dgf_main











!======= READ ======= 
subroutine read_kb_contour_gf_main(G,file)
  type(kb_contour_gf)  :: G
  character(len=*)     :: file
  logical              :: check
  check = inquire_kb_contour_gf(file)
  if(.not.G%status.OR..not.check)stop "contour_gf/read_kb_contour_gf: G not allocated"
  call read_array(trim(file)//"_less.data",G%less(:,:))
  call read_array(trim(file)//"_ret.data",G%ret(:,:))
  call read_array(trim(file)//"_lmix.data",G%lmix(:,0:))
  call read_array(trim(file)//"_mats.data",G%mats(0:))
  call read_array(trim(file)//"_iw.data",G%iw(:))
end subroutine read_kb_contour_gf_main
!
subroutine read_kb_contour_dgf_main(dG,file)
  type(kb_contour_dgf) :: dG
  character(len=*)     :: file
  logical              :: check
  check = inquire_kb_contour_dgf(file)
  if(.not.dG%status.OR..not.check)stop "contour_gf/read_kb_contour_dgf: dG not allocated"
  call read_array(trim(file)//"_less.data",dG%less(:))
  call read_array(trim(file)//"_ret.data",dG%ret(:))
  call read_array(trim(file)//"_lmix.data",dG%lmix(0:))
end subroutine read_kb_contour_dgf_main
!
!
function inquire_kb_contour_gf(file) result(check)
  integer          :: i
  logical          :: check,bool(5)
  character(len=*) :: file
  character(len=16),dimension(5)  :: ctype=([ character(len=5) :: 'less','ret','lmix','mats','iw'])
  check=.true.
  do i=1,5
     inquire(file=reg(file)//"_"//reg(ctype(i))//".data",exist=bool(i))
     if(.not.bool(i))inquire(file=reg(file)//"_"//reg(ctype(i))//".data.gz",exist=bool(i))
     check=check.AND.bool(i)
  enddo
end function inquire_kb_contour_gf
!
function inquire_kb_contour_dgf(file) result(check)
  integer          :: i
  logical          :: check,bool(3)
  character(len=*) :: file
  character(len=16),dimension(3)  :: ctype=([character(len=5) :: 'less','ret','lmix'])
  check=.true.
  do i=1,3
     inquire(file=reg(file)//"_"//reg(ctype(i))//".data",exist=bool(i))
     if(.not.bool(i))inquire(file=reg(file)//"_"//reg(ctype(i))//".data.gz",exist=bool(i))
     check=check.AND.bool(i)
  enddo
end function inquire_kb_contour_dgf






!======= PLOT ======= 
subroutine plot_kb_contour_gf_main(params,file,G)
  character(len=*)        :: file
  type(kb_contour_gf)     :: G
  type(kb_contour_params) :: params
  integer                 :: Nt
  if(.not.G%status)stop "contour_gf/plot_kb_contour_gf: G is not allocated" 
  Nt=params%Ntime
  call splot3d(reg(file)//"_less_t_t.neqipt",params%t(:Nt),params%t(:Nt),G%less(:Nt,:Nt))
  call splot3d(reg(file)//"_ret_t_t.neqipt",params%t(:Nt),params%t(:Nt),G%ret(:Nt,:Nt))
  call splot3d(reg(file)//"_lmix_t_tau.neqipt",params%t(:Nt),params%tau(0:),G%lmix(:Nt,0:))
  call splot(reg(file)//"_mats_tau.neqipt",params%tau(0:),G%mats(0:))
  call splot(reg(file)//"_mats_iw.neqipt",params%wm(:),G%iw(:))
end subroutine plot_kb_contour_gf_main
!
subroutine plot_kb_contour_dgf_main(params,file,dG)
  character(len=*)        :: file
  type(kb_contour_dgf)    :: dG
  type(kb_contour_params) :: params
  integer                 :: Nt
  if(.not.dG%status)stop "contour_gf/plot_kb_contour_gf: G is not allocated" 
  Nt=params%Ntime
  call splot(reg(file)//"_less_t.neqipt",params%t(:Nt),dG%less(:Nt))
  call splot(reg(file)//"_ret_t.neqipt",params%t(:Nt),dG%ret(:Nt))
  call splot(reg(file)//"_lmix_tau.neqipt",params%tau(0:),dG%lmix(0:))
end subroutine plot_kb_contour_dgf_main
