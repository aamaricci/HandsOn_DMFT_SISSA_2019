subroutine check_kb_contour_gf_main(params,G,pname)
  type(kb_contour_gf)       :: G
  type(kb_contour_params)   :: params
  integer                   :: N,L,Lf
  character(len=*),optional :: pname
  character(len=100)        :: pname_
  !
  pname_="check_kb_contour_gf_main";if(present(pname))pname_=reg(pname)
  !
  N=params%Ntime              !<== check size at maximum time
  L=params%Ntau
  !
  if(.not.G%status) stop "check_kb_contour_gf_main: G not allocated"
  call assert_shape(G%less,[N,N],reg(pname_),"G%less")
  call assert_shape(G%ret,[N,N],reg(pname_),"G%ret")
  call assert_shape(G%lmix,[N,L+1],reg(pname_),"G%lmix")
  call assert_shape(G%mats,[L+1],reg(pname_),"G%mats")
end subroutine check_kb_contour_gf_main
!
subroutine check_kb_contour_gf_Nso(params,G,pname)
  type(kb_contour_gf)       :: G(:,:,:,:) !![Nspin][Nspin][Norb][Norb]
  type(kb_contour_params)   :: params
  integer                   :: N,L,Lf,Nspin,Norb
  integer                   :: ispin,jspin,iorb,jorb
  character(len=*),optional :: pname
  character(len=100)        :: pname_
  !
  pname_="check_kb_contour_gf_Nso";if(present(pname))pname_=reg(pname)
  !
  Nspin = size(G,1)
  Norb  = size(G,3)
  if(any( shape(G) /= [Nspin,Nspin,Norb,Norb] ) ) stop "ERROR check_kb_contour_gf_Nso: shape[G] != [Nspin][Nspin][Norb][Norb]"
  !
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              call check_kb_contour_gf_main(params,G(ispin,jspin,iorb,jorb),trim(pname_))
           enddo
        enddo
     enddo
  enddo
end subroutine check_kb_contour_gf_Nso





subroutine check_kb_contour_dgf_main(params,dG,pname)
  type(kb_contour_dgf)      :: dG
  type(kb_contour_params)   :: params
  integer                   :: N,L
  character(len=*),optional :: pname
  character(len=100)        :: pname_
  pname_="check_kb_contour_dgf";if(present(pname))pname_=reg(pname)
  !
  N=params%Ntime              !<== check size at maximum time
  L=params%Ntau
  !
  if(.not.dG%status) stop "check_kb_contour_dgf: dG not allocated"
  call assert_shape(dG%less,[N],reg(pname_),"dG%less")
  call assert_shape(dG%ret,[N],reg(pname_),"dG%ret")
  call assert_shape(dG%lmix,[L+1],reg(pname_),"dG%lmix")
end subroutine check_kb_contour_dgf_main

subroutine check_kb_contour_dgf_Nso(params,dG,pname)
  type(kb_contour_dgf)      :: dG(:,:,:,:) !![Nspin][Nspin][Norb][Norb]
  type(kb_contour_params)   :: params
  integer                   :: N,L,Lf,Nspin,Norb
  integer                   :: ispin,jspin,iorb,jorb
  character(len=*),optional :: pname
  character(len=100)        :: pname_
  !
  pname_="check_kb_contour_dgf_Nso";if(present(pname))pname_=reg(pname)
  !
  Nspin = size(dG,1)
  Norb  = size(dG,3)
  if(any( shape(dG) /= [Nspin,Nspin,Norb,Norb] ) ) stop "ERROR check_kb_contour_dgf_Nso: shape[dG] != [Nspin][Nspin][Norb][Norb]"
  !
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              call check_kb_contour_dgf_main(params,dG(ispin,jspin,iorb,jorb),trim(pname_))
           enddo
        enddo
     enddo
  enddo
end subroutine check_kb_contour_dgf_Nso
