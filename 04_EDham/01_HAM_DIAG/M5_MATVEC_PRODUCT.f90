MODULE MATVEC_PRODUCT
  USE COMMON_VARS
  USE SPARSE_MATRIX
  USE LANCZOS
  USE SETUP
  implicit none

  private

  public :: Build_SpH
  public :: Delete_SpH
  !
  public :: HxV_spH0


contains




  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine Build_spH(Nups,Ndws)
    integer                       :: Nups,Ndws   
    integer,dimension(Ns)         :: nup,ndw
    integer                       :: Dim
    integer                       :: i,iup,idw
    integer                       :: j,jup,jdw
    integer                       :: m,mup,mdw
    integer                       :: ms
    integer                       :: impi,impi_up,impi_dw
    integer                       :: iorb,jorb,ispin,jspin,ibath
    integer                       :: kp,k1,k2,k3,k4
    integer                       :: alfa,beta,ie
    real(8)                       :: sg1,sg2,sg3,sg4
    real(8)                       :: htmp,htmpup,htmpdw,deps
    !Hamiltonian hard coded paramters:
    real(8)                       :: xmu
    real(8)                       :: Vps,Uloc,Eps(Ns-1)
    !
    type(sector_map),dimension(2) :: Hs
    !
    !< generate Hamiltonian parameters:
    Vps  = 1d0/sqrt(dble(Ns-1))
    deps = 2d0/(Ns-1)
    do ie=1,Ns-1
       eps(ie) = -1d0 + (ie-1)*deps
    enddo
    Uloc = 3d0
    !
    call build_sector(Nups,Ndws,Hs)
    !
    DimUp = get_sector_dimension(Ns,Nups)
    DimDw = get_sector_dimension(Ns,Ndws)
    Dim=DimUp*DimDw
    !
    !
    call sp_init_matrix(spH0,Dim)
    call sp_init_matrix(spH0dw,DimDw)
    call sp_init_matrix(spH0up,DimUp)
    !
    !LOCAL PART
    do i=1,Dim
       iup = iup_index(i,DimUp)
       idw = idw_index(i,DimUp)
       !
       nup  = bdecomp(Hs(1)%map(iup),Ns)
       ndw  = bdecomp(Hs(2)%map(idw),Ns)
       !
       !Diagonal Elements, i.e. local part
       htmp = 0d0
       htmp = htmp + Uloc*(nup(1)*ndw(1))
       htmp = htmp - 0.5d0*Uloc*(nup(1)+ndw(1)) + 0.25d0*uloc
       !
       !> H_Bath: local bath energy contribution.
       !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
       do kp=1,Ns-1
          alfa = 1 + kp
          htmp =htmp + eps(kp)*nup(alfa) !UP
          htmp =htmp + eps(kp)*ndw(alfa) !DW
       enddo
       !
       call sp_insert_element(spH0,htmp,i,i)
    enddo
    !
    ! IMP UP <--> BATH UP
    do iup=1,DimUp
       mup  = Hs(1)%map(iup)
       nup  = bdecomp(mup,Ns)
       do kp=1,Ns-1
          alfa = 1 + kp
          if( (nup(1)==1) .AND. (nup(alfa)==0) )then              
             call c(1,mup,k1,sg1)
             call cdg(alfa,k1,k2,sg2)
             jup = binary_search(Hs(1)%map,k2)
             htmp = vps*sg1*sg2
             call sp_insert_element(spH0up,htmp,iup,jup)
          endif
          if( (nup(1)==0) .AND. (nup(alfa)==1) )then
             call c(alfa,mup,k1,sg1)
             call cdg(1,k1,k2,sg2)
             jup=binary_search(Hs(1)%map,k2)
             htmp = vps*sg1*sg2
             call sp_insert_element(spH0up,htmp,iup,jup)
          endif
       enddo
    enddo
    !
    !IMP DW <--> BATH DW
    do idw=1,DimDw
       mdw  = Hs(2)%map(idw)
       ndw  = bdecomp(mdw,Ns)
       do kp=1,Ns-1
          alfa = 1 + kp
          if( (ndw(1)==1) .AND. (ndw(alfa)==0) )then
             call c(1,mdw,k1,sg1)
             call cdg(alfa,k1,k2,sg2)
             jdw=binary_search(Hs(2)%map,k2)
             htmp=vps*sg1*sg2
             call sp_insert_element(spH0dw,htmp,idw,jdw)
          endif
          !
          !
          if( (ndw(1)==0) .AND. (ndw(alfa)==1) )then
             call c(alfa,mdw,k1,sg1)
             call cdg(1,k1,k2,sg2)
             jdw=binary_search(Hs(2)%map,k2)
             htmp=vps*sg1*sg2
             call sp_insert_element(spH0dw,htmp,idw,jdw)
          endif
       enddo
    enddo
    !
    call Delete_sector(Hs)
    !
  end subroutine Build_spH



  subroutine Delete_spH
    if(spH0%status)call sp_delete_matrix(spH0)
    if(spH0dw%status)call sp_delete_matrix(spH0dw)
    if(spH0up%status)call sp_delete_matrix(spH0up)
  end subroutine Delete_spH




  subroutine HxV_spH0(N,v,Hv)
    integer                         :: N
    real(8),dimension(N)            :: v
    real(8),dimension(N)            :: Hv
    integer                         :: i,iup,idw,j,jup,jdw,jj
    !
    Hv=0d0
    !
    do i=1,N
       do j=1,spH0%row(i)%Size
          Hv(i) = Hv(i) + spH0%row(i)%vals(j)*v(spH0%row(i)%cols(j))
       end do
    end do
    !
    do idw=1,DimDw
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          do jj=1,spH0dw%row(idw)%Size
             j     = iup +  (spH0dw%row(idw)%cols(jj)-1)*DimUp
             Hv(i) = Hv(i) + spH0dw%row(idw)%vals(jj)*V(j)
          enddo
       enddo
    enddo
    !
    do idw=1,DimDw
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp          
          do jj=1,spH0up%row(iup)%Size
             j = spH0up%row(iup)%cols(jj) + (idw-1)*DimUp
             Hv(i) = Hv(i) + spH0up%row(iup)%vals(jj)*V(j)
          enddo
       enddo
    enddo
    !
    return
  end subroutine HxV_spH0







end module MATVEC_PRODUCT







