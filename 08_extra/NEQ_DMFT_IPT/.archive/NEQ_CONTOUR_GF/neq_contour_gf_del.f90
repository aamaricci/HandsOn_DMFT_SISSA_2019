subroutine del_kb_contour_gf_main(params,G)
  type(kb_contour_params) :: params
  type(kb_contour_gf)     :: G
  integer                 :: N,L
  !
  call check_kb_contour_gf(params,G,"del_kb_contour_gf_main")
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  if(N==1)then
     G%iw = zero
     G%mats = 0d0
  endif
  G%ret(N,1:N)   = zero
  G%less(N,1:N)  = zero
  G%lmix(N,0:)   = zero
end subroutine del_kb_contour_gf_main

subroutine del_kb_contour_gf_Nso(params,G)
  type(kb_contour_params) :: params
  type(kb_contour_gf)     :: G(:,:,:,:) ![Nspin][:][Norb][:]
  integer                 :: N,L,Nspin,Norb
  integer                 :: ispin,jspin,iorb,jorb
  !
  call check_kb_contour_gf(params,G,"del_kb_contour_gf_Nso")
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  Nspin = size(G,1)
  Norb  = size(G,3)
  !
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              if(N==1)then
                 G(ispin,jspin,iorb,jorb)%iw = zero
                 G(ispin,jspin,iorb,jorb)%mats = 0d0
              endif
              G(ispin,jspin,iorb,jorb)%ret(N,1:N)   = zero
              G(ispin,jspin,iorb,jorb)%less(N,1:N)  = zero
              G(ispin,jspin,iorb,jorb)%lmix(N,0:)   = zero
           enddo
        enddo
     enddo
  enddo
end subroutine del_kb_contour_gf_Nso





subroutine del_kb_contour_dgf_main(params,dG)
  type(kb_contour_params) :: params
  type(kb_contour_dgf)    :: dG
  integer                 :: N,L
  !
  call check_kb_contour_gf(params,dG,"del_kb_contour_dgf_main")
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  dG%ret(1:N)   = zero
  dG%less(1:N)  = zero
  dG%lmix(0:)   = zero
end subroutine del_kb_contour_dgf_main

subroutine del_kb_contour_dgf_Nso(params,dG)
  type(kb_contour_params) :: params
  type(kb_contour_dgf)    :: dG(:,:,:,:) ![Nspin][:][Norb][:]
  integer                 :: N,L,Nspin,Norb
  integer                 :: ispin,jspin,iorb,jorb
  !
  call check_kb_contour_gf(params,dG,"del_kb_contour_dgf_Nso")
  !
  N   = params%Nt   !<== work with the ACTUAL size of the contour
  L   = params%Ntau
  !
  Nspin = size(dG,1)
  Norb  = size(dG,3)
  !
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              dG(ispin,jspin,iorb,jorb)%ret(1:N)   = zero
              dG(ispin,jspin,iorb,jorb)%less(1:N)  = zero
              dG(ispin,jspin,iorb,jorb)%lmix(0:)   = zero
           enddo
        enddo
     enddo
  enddo
end subroutine del_kb_contour_dgf_Nso
