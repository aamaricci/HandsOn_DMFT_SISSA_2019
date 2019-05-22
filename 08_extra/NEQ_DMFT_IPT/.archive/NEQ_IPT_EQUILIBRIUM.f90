module NEQ_IPT_EQUILIBRIUM
  USE NEQ_INPUT_VARS
  USE NEQ_FFT_GF
  USE SF_LINALG, only: inv
  USE SF_CONSTANTS, only: xi,pi,one
  USE SF_ARRAYS, only:arange
  USE SF_SPECIAL, only: bethe_lattice
  USE SF_IOTOOLS, only: free_unit
  implicit none
  private


  !half-filling solver:
  public :: ipt_solve_matsubara

  public :: ipt_measure_potential_energy_matsubara
  public :: ipt_measure_kinetic_energy_matsubara
  public :: ipt_measure_hartree_energy_matsubara
  public :: ipt_measure_dens_matsubara
  public :: ipt_measure_docc_matsubara
  public :: ipt_measure_zeta_matsubara
  !


contains


  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order perturbation theory in Matsubara normal: 
  !+-------------------------------------------------------------------+
  subroutine ipt_solve_matsubara(weiss,sigma)
    complex(8),dimension(:)           :: weiss
    complex(8),dimension(size(weiss)) :: sigma
    real(8),dimension(0:size(weiss))  :: fg0_tau,sigma_tau
    integer                           :: Lf,i,unit
    !
    Lf = size(Weiss)
    !
    call fft_iw2tau(Weiss,fg0_tau(0:),beta)
    forall(i=0:Niw)sigma_tau(i)=Ui*Ui*fg0_tau(i)*fg0_tau(Niw-i)*fg0_tau(i)
    call fft_tau2iw(sigma_tau(0:),Sigma,beta)
    !
    open(free_unit(unit),file="Sigma_tau.ipt")
    do i=0,Niw
       write(unit,*)i*beta/Niw,sigma_tau(i)
    enddo
    close(unit)
    !
  end subroutine ipt_solve_matsubara



  !+-------------------------------------------------------------------+
  !PURPOSE: measure observables in Matsubara formalism
  !+-------------------------------------------------------------------+
  function ipt_measure_dens_matsubara(Green) result(ipt_dens)
    complex(8),dimension(:) :: green
    real(8)                 :: ipt_dens
    ipt_dens = fft_gbeta_minus(Green,beta)
  end function ipt_measure_dens_matsubara



  !PURPOSE: measure renormalization constant zeta
  function ipt_measure_zeta_matsubara(Sigma) result(ipt_zeta)
    complex(8),dimension(:) :: Sigma
    real(8)                 :: ipt_zeta,wm1
    wm1=pi/beta
    ipt_zeta = 1d0 - dimag(Sigma(1))/wm1
    ipt_zeta = 1d0/ipt_zeta
  end function ipt_measure_zeta_matsubara




  !PURPOSE: measure double occupancy
  function ipt_measure_docc_matsubara(Green,Sigma) result(ipt_docc)
    complex(8),dimension(:)           :: Green
    complex(8),dimension(size(Green)) :: Sigma
    real(8)                           :: ipt_docc,epot,ehar
    !
    epot = ipt_measure_potential_energy_matsubara(Green,Sigma)
    ehar = ipt_measure_hartree_energy_matsubara(Green)
    ipt_docc = 0.25d0
    if(Ui > 0d0)ipt_docc = epot/Ui - ehar/Ui
  end function ipt_measure_docc_matsubara





  !PURPOSE: measure potential energy
  function ipt_measure_potential_energy_matsubara(Green,Sigma) result(ipt_Epot)
    complex(8),dimension(:)           :: Green
    complex(8),dimension(size(Green)) :: Sigma
    real(8)                           :: ipt_Epot
    ipt_Epot=sum(dreal(Sigma)*dreal(Green))
    ipt_Epot=ipt_Epot-sum(dimag(Sigma)*dimag(Green))
    ipt_Epot=ipt_Epot/beta*2d0
  end function ipt_measure_potential_energy_matsubara



  !PURPOSE: measure hartree energy term
  function ipt_measure_hartree_energy_matsubara(Green) result(ipt_Ehartree)
    complex(8),dimension(:)           :: Green
    real(8)                           :: ipt_Ehartree,n
    n = fft_gbeta_minus(Green,beta)
    ipt_Ehartree = -Ui*n + Ui*0.25d0 
  end function ipt_measure_hartree_energy_matsubara





  !PURPOSE: measure kinetic energy
  function ipt_measure_kinetic_energy_matsubara(Hk,Wtk,Sigma) result(ipt_Ekin)
    complex(8),dimension(:,:,:)   :: Hk ![Nso][Nso][Lk]
    real(8),dimension(size(Hk),3) :: Wtk
    complex(8),dimension(:)       :: Sigma
    real(8)                       :: ipt_Ekin
    ipt_Ekin = f_ipt_kinetic_normal(Hk(1,1,:),Wtk,Sigma)
  end function ipt_measure_kinetic_energy_matsubara

  function f_ipt_kinetic_normal(Hk,Wtk,Sigma) result(ipt_Ekin)
    complex(8),dimension(:)          :: Hk
    real(8),dimension(size(Hk))      :: Wtk
    complex(8),dimension(:)          :: Sigma
    integer                          :: Lk,No,Liw
    integer                          :: i,ik,iorb
    real(8),dimension(:),allocatable :: wm
    real(8)                          :: Sigma_HF
    complex(8)                       :: Ak,Bk
    complex(8)                       :: Ck,Zk
    complex(8)                       :: Zeta,Gk,Tk
    real(8)                          :: Tail0,Tail1,spin_degeneracy
    !
    real(8)                          :: H0,ipt_Ekin
    !
    No = 1
    Lk = size(Hk)
    Liw= size(Sigma)
    !
    allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    !
    Sigma_HF = dreal(Sigma(Liw))
    !
    H0=0d0
    Zk=1d0
    do ik=1,Lk
       Ak= Hk(ik)
       Bk=-Hk(ik)-Sigma_HF
       do i=1,Liw
          Gk = (xi*wm(i)+xmu)*Zk - Hk(ik) - Sigma(i)
          Gk = 1d0/Gk
          Tk = Zk/(xi*wm(i)) - Bk/(xi*wm(i))**2
          Ck = Ak*(Gk - Tk)
          H0 = H0 + Wtk(ik)*Ck
       enddo
    enddo
    spin_degeneracy=2!3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2.d0*spin_degeneracy
    !
    Tail0=0d0
    Tail1=0d0
    do ik=1,Lk
       Ak= Hk(ik)
       Bk=-Hk(ik)-Sigma_HF
       Ck= Ak*Bk
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*Ak
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*Ck
    enddo
    Tail0=spin_degeneracy*Tail0
    Tail1=spin_degeneracy*Tail1*beta
    ipt_Ekin=H0+Tail0+Tail1
    deallocate(wm)
  end function f_ipt_kinetic_normal














end module NEQ_IPT_EQUILIBRIUM
