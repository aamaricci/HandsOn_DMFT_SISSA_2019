subroutine allocate_kb_contour_gf_main(params,G)
  type(kb_contour_gf)     :: G
  type(kb_contour_params) :: params
  integer                 :: i,j,N,L,Lf
  !
  if(allocated(G%less))deallocate(G%less)
  if(allocated(G%ret)) deallocate(G%ret)
  if(allocated(G%lmix))deallocate(G%lmix)
  if(allocated(G%mats))deallocate(G%mats)
  if(allocated(G%iw))deallocate(G%iw)
  N = params%Ntime            !<== allocate at maximum time
  L = params%Ntau
  Lf= params%Niw
  allocate(G%less(N,N))  ; G%less=zero
  allocate(G%ret(N,N))   ; G%ret=zero
  allocate(G%lmix(N,0:L)); G%lmix=zero
  allocate(G%mats(0:L))  ; G%mats=0d0
  allocate(G%iw(Lf))     ; G%iw=zero
  G%status=.true.
  !
end subroutine allocate_kb_contour_gf_main

subroutine allocate_kb_contour_gf_Nso(params,G)
  type(kb_contour_gf)     :: G(:,:,:,:) ![Nspin][Nspin][Norb][Norb]
  type(kb_contour_params) :: params
  integer                 :: i,j,N,L,Lf,Nspin,Norb,ispin,jspin,iorb,jorb
  Nspin = size(G,1)
  Norb  = size(G,3)
  if(any( shape(G) /= [Nspin,Nspin,Norb,Norb] ) ) stop "ERROR allocate_kb_contour_gf_Nso: shape[G] != [Nspin][Nspin][Norb][Norb]"
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              call allocate_kb_contour_gf_main(params,G(ispin,jspin,iorb,jorb))
           enddo
        enddo
     enddo
  enddo
  !
end subroutine allocate_kb_contour_gf_Nso










!##################################################################
!##################################################################
!##################################################################
!##################################################################









subroutine allocate_kb_contour_dgf_main(params,dG,wgtr)
  type(kb_contour_dgf)    :: dG
  type(kb_contour_params) :: params
  integer                 :: i,j,N,L
  logical,optional        :: wgtr
  if(allocated(dG%less))deallocate(dG%less)
  if(allocated(dG%ret)) deallocate(dG%ret)
  if(allocated(dG%lmix))deallocate(dG%lmix)
  N=params%Ntime           !<== allocate at maximum time
  L=params%Ntau
  allocate(dG%less(N))  ; dG%less=zero
  allocate(dG%ret(N))   ; dG%ret=zero
  allocate(dG%lmix(0:L)); dG%lmix=zero
  if(present(wgtr).AND.wgtr)then
     allocate(dG%gtr(N))
     dG%gtr=zero
  endif
  dG%status=.true.
end subroutine allocate_kb_contour_dgf_main

subroutine allocate_kb_contour_dgf_Nso(params,dG)
  type(kb_contour_dgf)     :: dG(:,:,:,:) ![Nspin][Nspin][Norb][Norb]
  type(kb_contour_params) :: params
  integer                 :: i,j,N,L,Lf,Nspin,Norb,ispin,jspin,iorb,jorb
  Nspin = size(dG,1)
  Norb  = size(dG,3)
  if(any( shape(dG) /= [Nspin,Nspin,Norb,Norb] ) ) stop "ERROR allocate_kb_contour_gf_Nso: shape[dG] != [Nspin][Nspin][Norb][Norb]"
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              call allocate_kb_contour_dgf_main(params,dG(ispin,jspin,iorb,jorb))
           enddo
        enddo
     enddo
  enddo
  !
end subroutine allocate_kb_contour_dgf_Nso







!##################################################################
!##################################################################
!##################################################################
!##################################################################






subroutine deallocate_kb_contour_gf_main(G)
  type(kb_contour_gf) :: G
  if(.not.G%status)stop "contour_gf/deallocate_kb_contour_gf: G not allocated"
  deallocate(G%less,G%ret,G%lmix,G%mats,G%iw)
  G%status=.false.
end subroutine deallocate_kb_contour_gf_main

subroutine deallocate_kb_contour_gf_Nso(G)
  type(kb_contour_gf)     :: G(:,:,:,:) ![Nspin][Nspin][Norb][Norb]
  integer                 :: i,j,N,L,Lf,Nspin,Norb,ispin,jspin,iorb,jorb
  Nspin = size(G,1)
  Norb  = size(G,3)
  if(any( shape(G) /= [Nspin,Nspin,Norb,Norb] ) ) stop "ERROR deallocate_kb_contour_gf_Nso: shape[G] != [Nspin][Nspin][Norb][Norb]"
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              call deallocate_kb_contour_gf_main(G(ispin,jspin,iorb,jorb))
           enddo
        enddo
     enddo
  enddo
  !
end subroutine deallocate_kb_contour_gf_Nso








!##################################################################
!##################################################################
!##################################################################
!##################################################################






subroutine deallocate_kb_contour_dgf_main(dG)
  type(kb_contour_dgf) :: dG
  if(.not.dG%status)stop "contour_gf/deallocate_kb_contour_dgf: dG not allocated"
  deallocate(dG%less,dG%ret,dG%lmix)
  if(allocated(dG%gtr))deallocate(dG%gtr)
  dG%status=.false.
end subroutine deallocate_kb_contour_dgf_main

subroutine deallocate_kb_contour_dgf_Nso(dG)
  type(kb_contour_dgf)     :: dG(:,:,:,:) ![Nspin][Nspin][Norb][Norb]
  integer                 :: i,j,N,L,Lf,Nspin,Norb,ispin,jspin,iorb,jorb
  Nspin = size(dG,1)
  Norb  = size(dG,3)
  if(any( shape(dG) /= [Nspin,Nspin,Norb,Norb] ) ) stop "ERROR deallocate_kb_contour_gf_Nso: shape[dG] != [Nspin][Nspin][Norb][Norb]"
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              call deallocate_kb_contour_dgf_main(dG(ispin,jspin,iorb,jorb))
           enddo
        enddo
     enddo
  enddo
  !
end subroutine deallocate_kb_contour_dgf_Nso
