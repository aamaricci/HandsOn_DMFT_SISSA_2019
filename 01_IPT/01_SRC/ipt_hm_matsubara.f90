program hmipt_matsubara
  USE IPT_FFTGF
  implicit none
  !
  real(8),parameter    :: pi = acos(-1d0), D=1d0
  complex(8),parameter :: xi = (0d0,1d0),one=(1d0,0d0),zero=(0d0,0d0)
  !
  integer                :: nloop          !dmft loop variables
  real(8)                :: uloc           !local  interaction
  real(8)                :: beta           !inverse temperature
  real(8)                :: dmft_error     !DMFT error threshold
  logical                :: converged,check
  real(8)                :: wmix
  integer                :: i,iloop,L
  complex(8)             :: zeta
  !
  complex(8),allocatable :: sigma(:)
  complex(8),allocatable :: fg0(:),fg0_prev(:)
  complex(8),allocatable :: fg(:)
  !
  real(8),allocatable    :: wm(:),fg0_tau(:),sigma_tau(:)
  real(8)                :: dens,docc,z
  !
  !A fortran way to read input variables:
  namelist/ipt_variable/Nloop,Uloc,Beta,L,dmft_error,wmix
  Nloop=100
  Beta=100d0
  Uloc=2d0
  DMFT_error=1d-5
  L=2*4096
  wmix=0.5d0
  !
  open(999,file="inputIPT.conf")
  read(999,nml=ipt_variable)
  write(*,nml=ipt_variable)     !write on the screen the used input
  close(999)


  !allocate functions:
  allocate(fg(L));fg=zero
  allocate(sigma(L));sigma=zero
  allocate(fg0(L));fg0=zero
  allocate(fg0_prev(L));fg0_prev=zero



  !build freq. array
  allocate(wm(L))
  do i=1,L
     wm(i)  = pi/beta*(2*i-1)
  enddo


  !get or read first sigma 
  call  get_initial_function(Sigma,"Sigma.restart")



  !dmft loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,"(A,i5)",advance="no")"DMFT-loop",iloop
     !
     !     
     !SELF-CONSISTENCY: 
     do i=1,L
        zeta  = xi*wm(i) - sigma(i)
        fg(i) = gfbethe(wm(i),zeta,D)
     enddo
     !
     fg0_prev= fg0                       !
     fg0     = one/(one/fg + sigma)
     if(iloop>1)fg0 = wmix*fg0 + (1.d0-wmix)*fg0_prev !mix to avoid loops
     !
     !IMPURITY SOLVER: fg0-->sigma
     call solve_ipt_mats(sigma)
     !
     dens = ipt_measure_dens_matsubara(fg)
     z    = ipt_measure_zeta_matsubara(sigma,fg0)
     docc = ipt_measure_docc_matsubara(sigma,fg0)
     write(*,"(3F15.9,3x)",advance="no")dens,docc,z

     !Check CONVERGENCE on the Weiss Field:
     converged=check_convergence(fg0)
  enddo



  !WRITE Q.TIES TO FILES:
  open(101,file="G_iw.ipt")
  open(102,file="G0_iw.ipt")
  open(103,file="Sigma_iw.ipt")
  do i=1,L
     write(101,*)wm(i),dimag(fg(i)),dreal(fg(i))
     write(102,*)wm(i),dimag(fg0(i)),dreal(fg0(i))
     write(103,*)wm(i),dimag(sigma(i)),dreal(sigma(i))
  enddo
  close(101)
  close(102)
  close(103)


  open(999,file="observables_last.ipt")
  write(999,"(3F15.9,1x)")dens,docc,z
  close(999)


  open(103,file="Sigma.restart")
  do i=1,L
     write(103,*)sigma(i)
  enddo
  close(103)




contains


  !PURPOSE  : Solve 2^nd order perturbation theory
  subroutine solve_ipt_mats(sigma)
    complex(8),dimension(L) :: sigma
    real(8),dimension(0:L)  :: fg0_tau
    real(8),dimension(0:L)  :: sigma_tau
    call fftgf_iw2tau(fg0,fg0_tau(0:),beta)
    forall(i=0:L)sigma_tau(i)=Uloc*Uloc*fg0_tau(i)*fg0_tau(L-i)*fg0_tau(i)
    call fftgf_tau2iw(sigma_tau(0:),Sigma,beta)
    Sigma=xi*dimag(Sigma)
  end subroutine solve_ipt_mats



  !Get the initial Sigma: either by reading from file or assuming a HF form
  subroutine get_initial_function(self,file)
    complex(8),dimension(:) :: self
    character(len=*)        :: file
    logical                 :: check
    inquire(file=file,exist=check)
    if(check)then
       print*,'Reading sigma'
       open(100,file=file)
       do i=1,size(self)
          read(100,*)self(i)
       enddo
       close(100)
    endif
  end subroutine get_initial_function


  !Get the Bethe Lattice Local GF (Matsubara)
  elemental function gfbethe(w,zeta,d)
    real(8),intent(in)    :: w,d
    complex(8),intent(in) :: zeta
    complex(8)            :: gfbethe,sqroot
    real(8)               :: sq,sig
    sqroot=sqrt(zeta**2-d**2)
    sq=dimag(sqroot)
    sig=w*sq/abs(w*sq)
    gfbethe=2.d0/(zeta+sig*sqroot)
    return
  end function gfbethe



  function ipt_measure_dens_matsubara(Green) result(ipt_dens)
    complex(8),dimension(:)           :: Green
    real(8)                           :: ipt_dens,Gtau(0:size(green))
    call fftgf_iw2tau(Green,Gtau(0:),beta)
    ipt_dens = -Gtau(size(Sigma))
  end function ipt_measure_dens_matsubara

  function ipt_measure_zeta_matsubara(Sigma,Weiss) result(ipt_zeta)
    complex(8),dimension(:)           :: Sigma
    complex(8),dimension(size(Sigma)) :: Weiss
    real(8)                           :: ipt_zeta,wm1
    wm1=pi/beta
    ipt_zeta = 1d0 - dimag(Sigma(1))/wm1
    ipt_zeta = 1d0/ipt_zeta
  end function ipt_measure_zeta_matsubara

  function ipt_measure_docc_matsubara(Sigma,Weiss) result(ipt_docc)
    complex(8),dimension(:)           :: Sigma
    complex(8),dimension(size(Sigma)) :: Weiss
    complex(8),dimension(size(Sigma)) :: Green
    real(8)                           :: ipt_docc,epot,ehar,dens
    !
    Green = one/Weiss - Sigma
    Green = one/Green
    !
    Epot = sum(dreal(Sigma)*dreal(Green))
    Epot = Epot - sum(dimag(Sigma)*dimag(Green))
    Epot = Epot/beta*2d0
    !
    dens = ipt_measure_dens_matsubara(Green)
    ehar = -Uloc*dens + Uloc*0.25d0
    !
    ipt_docc = 0.25d0
    if(uloc > 0d0)ipt_docc = epot/uloc - ehar/uloc
  end function ipt_measure_docc_matsubara





  !Evaluate error and check convergence:
  !Xnew: base array to check convergence
  function check_convergence(Xnew) result(convergence)
    complex(8),intent(in)       :: Xnew(:)
    integer                     :: i,j,Msum
    logical                     :: convergence  
    real(8)                     :: err
    real(8)                     :: M,S
    complex(8),save,allocatable :: Xold(:)
    integer,save                :: success=0,check=1
    Msum=size(Xnew)
    if(.not.allocated(Xold))allocate(Xold(Msum))
    S=0d0 ; M=0d0
    do i=1,Msum
       M=M + abs(Xnew(i)-Xold(i))
       S=S + abs(Xnew(i))
    enddo
    err= M/S
    Xold=Xnew
    !
    if(err < dmft_error)then
       success=success+1
    else
       success=0
    endif
    !
    convergence=.false.
    if(success > 1)convergence=.true.
    !
    if(check>=Nloop)convergence=.true.
    check=check+1
    !
    write(*,"(A,ES15.7)")"error=",err
  end function check_convergence

end program
