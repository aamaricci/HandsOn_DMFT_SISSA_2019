MODULE NEQ_MEASURE
  USE NEQ_INPUT_VARS
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE ELECTRIC_FIELD
  USE SF_CONSTANTS
  USE SF_IOTOOLS, only:free_unit
  implicit none

  private

  public :: measure_observables
  public :: measure_dens
  public :: measure_docc
  public :: measure_ekin
  public :: measure_epot
  public :: measure_etot
  public :: measure_current

contains



  !+-------------------------------------------------------------------+
  !PURPOSE: measure some observables and print them
  !+-------------------------------------------------------------------+
  subroutine measure_observables(params,g,self)
    type(kb_contour_gf)     :: g
    type(kb_contour_gf)     :: self
    type(kb_contour_params) :: params
    integer                 :: unit,itime
    real(8)                 :: dens,docc,ekin,epot,etot
    itime = params%Nt
    unit = free_unit()
    open(unit,file="observables.info")
    write(unit,"(8A20)")"time","N","Docc","Ekin","Epot","Etot"
    close(unit)
    open(unit,file="observables.neqipt",position="append")
    dens = measure_dens(params,g,self)
    docc = measure_docc(params,g,self)
    ekin = measure_ekin(params,g,self)
    epot = measure_epot(params,g,self)
    etot = ekin + epot
    write(unit,"(6F20.12)")params%t(itime),dens,docc,ekin,epot,etot
    close(unit)
  end subroutine measure_observables



  !+-------------------------------------------------------------------+
  !PURPOSE: measure current
  !+-------------------------------------------------------------------+
  subroutine measure_current(params,Gk,Vkt,Wtk)
    type(kb_contour_gf)            :: gk(:)
    real(8),dimension(:,:,:)       :: Vkt
    real(8),dimension(:)           :: wtk
    type(kb_contour_params)        :: params
    integer                        :: unit,itime,Lk,ik,i
    real(8),dimension(size(Vkt,3)) :: Jloc
    real(8)                        :: nkt
    !
    Lk=size(gk)
    itime = params%Nt
    !
    if(size(Vkt,1)<params%Ntime)stop "neq_measure_current: dim(Vkt,2) < Ntime"
    if(size(Vkt,2)/=Lk)stop "neq_measure_current: dim(Vkt,3) != Lk"
    !
    unit = free_unit()
    open(unit,file="current.info")
    write(unit,"(8A20)")"time","Jx","Jy","Jz"
    close(unit)
    !
    open(unit,file="current.neqipt",position="append")
    Jloc=0d0
    do ik=1,Lk
       nkt  = dimag(Gk(ik)%less(itime,itime))
       Jloc = Jloc + Wtk(ik)*Nkt*Vkt(itime,ik,:)
    enddo
    write(unit,"(4F20.12)")params%t(itime),(Jloc(i),i=1,size(Jloc))
    close(unit)
  end subroutine measure_current





  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the density at a given istant of time
  ! n(t)=-xi*G^<(t,t)
  !+-------------------------------------------------------------------+
  function measure_dens(params,g,self) result(dens)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: dens
    integer                             :: N
    N = params%Nt
    dens = dimag(G%less(N,N))
  end function measure_dens

  function measure_dens_Nso(params,g,self) result(dens)
    type(kb_contour_gf)                 :: g(Nspin,Nspin,Norb,Norb)
    type(kb_contour_gf)                 :: self(Nspin,Nspin,Norb,Norb)
    type(kb_contour_params)             :: params
    real(8)                             :: dens(Nspin,Norb)
    integer                             :: N,ispin,iorb
    N = params%Nt
    do ispin=1,Nspin
       do iorb=1,Norb
          dens(ispin,iorb) = dimag(G(ispin,ispin,iorb,iorb)%less(N,N))
       enddo
    enddo
  end function measure_dens_Nso



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the double occupancy at a given istant of time
  ! d(t)=n_up(t)*n_do(t)-1/U0*[Self^M*G^M]
  !      n_up(t)*n_do(t)-i/U*[Self^R*G^< + Self^<*G^A + Self^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_docc(params,g,self) result(docc)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: docc
    integer                             :: i,k,j,N,L
    complex(8),dimension(:),allocatable :: SxG
    real(8)                             :: nt
    N = params%Nt
    L = params%Ntau
    !
    nt   = dimag(G%less(N,N))
    allocate(SxG(0:max(N,L)))
    docc = nt**2
    if(N==1)then
       if(ui/=0.d0)then
          do k=0,L
             SxG(k)=Self%mats(L-k)*G%mats(k)
          end do
          docc=docc-1.d0/Ui*params%dtau*kb_trapz(SxG(0:),0,L)
       endif
    else
       if(u/=0.d0)then
          do k=0,L
             SxG(k)=Self%lmix(N,k)*conjg(G%lmix(N,L-k))
          end do
          docc=docc + 1.d0/U*params%dtau*dimag( (-xi)*kb_trapz(SxG(0:),0,L) )
          do k=1,N
             SxG(k)=Self%ret(N,k)*G%less(k,N)
          end do
          docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
          do k=1,N
             SxG(k)=Self%less(N,k)*conjg(G%ret(N,k))
          end do
          docc=docc + 1.d0/U*params%dt*dimag(kb_trapz(SxG(0:),1,N))
       endif
    endif
    deallocate(SxG)
  end function measure_docc


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_ekin(params,g,self) result(ekin)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: ekin
    integer                             :: i,k,j,N,L
    complex(8),dimension(:),allocatable :: Ker
    real(8)                             :: nt
    N = params%Nt
    L = params%Ntau
    !
    allocate(Ker(0:max(N,L)))
    if(N==1)then
       do k=0,L
          Ker(k)=G%mats(L-k)*G%mats(k)
       end do
       ekin = -2.d0*params%dtau*kb_trapz(Ker(0:),0,L)
    else
       do k=0,L
          Ker(k)=G%lmix(N,k)*conjg(G%lmix(N,L-k))
       end do
       ekin=2.d0*params%dtau*dimag( (-xi)*kb_trapz(Ker(0:),0,L) )
       do k=1,N
          Ker(k)=G%ret(N,k)*G%less(k,N)
       end do
       ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
       do k=1,N
          Ker(k)=G%less(N,k)*conjg(G%ret(N,k))
       end do
       ekin=ekin + 2.d0*params%dt*dimag(kb_trapz(Ker(0:),1,N))
    endif
    deallocate(Ker)
  end function measure_ekin



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! U(t)= U*docc(t) - n(t) + 1/4
  !+-------------------------------------------------------------------+
  function measure_epot(params,g,self) result(epot)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: epot,docc,nt
    integer                             :: i,k,j,N,L
    N = params%Nt
    L = params%Ntau
    !
    if(N==1)then
       nt   = measure_dens(params,g,self)
       docc = measure_docc(params,g,self)
       epot = Ui*(docc - nt + 0.25d0)
    else
       nt   = measure_dens(params,g,self)
       docc = measure_docc(params,g,self)
       epot = U*(docc - nt + 0.25d0)
    endif
  end function measure_epot



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_etot(params,g,self) result(etot)
    type(kb_contour_gf)                 :: g
    type(kb_contour_gf)                 :: self
    type(kb_contour_params)             :: params
    real(8)                             :: etot,ekin,epot
    ekin = measure_ekin(params,g,self)
    epot = measure_epot(params,g,self)
    etot = ekin + epot
  end function measure_etot




END MODULE NEQ_MEASURE
