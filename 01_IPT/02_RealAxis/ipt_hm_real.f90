program hmipt_realaxis
  implicit none
  !
  real(8),parameter                  :: pi = acos(-1d0), D=1d0
  complex(8),parameter               :: xi = (0d0,1d0),one=(1d0,0d0),zero=(0d0,0d0)
  !
  integer                            :: nloop          !dmft loop variables
  real(8)                            :: uloc           !local  interaction
  real(8)                            :: beta           !inverse temperature
  real(8)                            :: dmft_error     !DMFT error threshold
  logical                            :: converged,check
  real(8)                            :: wmix,eps
  integer                            :: i,iloop,L
  complex(8)                         :: zeta
  !
  complex(8),allocatable             :: sigma(:)
  complex(8),allocatable             :: fg0(:),fg0_prev(:)
  complex(8),allocatable             :: fg(:)
  !
  real(8),allocatable                :: wr(:)
  !
  real(8),dimension(:),allocatable   :: A0m,A0p,P1,P2
  integer,allocatable,dimension(:,:) :: iy_m_ix
  real(8)                            :: mesh


  !
  !A fortran way to read input variables:
  namelist/ipt_variable/Nloop,Uloc,Beta,L,dmft_error,wmix,eps
  Nloop=100
  Beta=100d0
  Uloc=2d0
  L=2*4096
  DMFT_error=1d-5
  wmix=0.5d0
  eps=0.01d0
  !
  open(999,file="inputIPT.conf")
  read(999,nml=ipt_variable)
  write(*,nml=ipt_variable)     !write on the screen the used input
  close(999)


  !allocate functions:
  allocate(fg(L))      ;fg=zero
  allocate(sigma(L))   ;sigma=zero
  allocate(fg0(L))     ;fg0=zero
  allocate(fg0_prev(L));fg0_prev=zero


  !build freq. array
  allocate(wr(L))
  mesh = 10d0/L
  do i=1,L
     wr(i) = -5d0 + (i-1)*mesh
  enddo

  !get or read first sigma 
  call  get_initial_function(Sigma,"Sigma.restart")

  call get_frequency_index

  !dmft loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,"(A,i5,A)",advance="no")"DMFT-loop",iloop," "
     !
     !     
     !SELF-CONSISTENCY:
     do i=1,L
        zeta = dcmplx(wr(i),eps) - sigma(i)
        fg(i) = gfbether(wr(i),zeta,D)
     enddo
     !
     fg0_prev= fg0                                    !
     fg0     = one/(one/fg + sigma)                   !
     if(iloop>1)fg0 = wmix*fg0 + (1.d0-wmix)*fg0_prev !mix to avoid loops
     !
     !IMPURITY SOLVER: fg0-->sigma
     call solve_ipt_real(sigma)
     !
     !Check CONVERGENCE on the Weiss Field:
     converged=check_convergence(fg0)
  enddo



  !WRITE Q.TIES TO FILES:
  open(101,file="G_wreal.ipt")
  open(102,file="G0_wreal.ipt")
  open(103,file="Sigma_wreal.ipt")
  do i=1,L
     write(101,*)wr(i),dimag(fg(i)),dreal(fg(i))
     write(102,*)wr(i),dimag(fg0(i)),dreal(fg0(i))
     write(103,*)wr(i),dimag(sigma(i)),dreal(sigma(i))
  enddo
  close(101)
  close(102)
  close(103)



  open(103,file="Sigma.restart")
  do i=1,L
     write(103,*)sigma(i)
  enddo
  close(103)




contains


  !PURPOSE  : Create an array of the indices y-x for a faster evaluation
  subroutine get_frequency_index()
    integer :: ix,iy,iz
    if(.not.allocated(iy_m_ix))allocate(iy_m_ix(L,L))
    iy_m_ix=0
    do ix=1,L
       do iy=1,L
          iz = iy - ix + L/2 
          if(iz<1 .OR. iz>L) iz=-1 !out of range-> if(iz>-L)
          iy_m_ix(iy,ix)=iz
       enddo
    enddo
    if(.not.allocated(A0m))allocate(A0m(L))
    if(.not.allocated(A0p))allocate(A0p(L))
    if(.not.allocated(P1)) allocate(P1(L))
    if(.not.allocated(P2)) allocate(P2(L))
  end subroutine get_frequency_index


  !PURPOSE  : Solve 2^nd order perturbation theory
  subroutine solve_ipt_real(sigma)
    complex(8),dimension(L) :: sigma
    integer                 :: ix,iy,iz
    real(8)                 :: sum1,sum2
    real(8),dimension(L)    :: reS,imS
    !
    call getAs
    call getPolarization
    !
    do ix=1,L
       sum1=zero
       sum2=zero
       do iy=1,L
          iz= iy_m_ix(iy,ix)
          if(iz>0)then
             sum1=sum1+A0p(iy)*P1(iz)*mesh
             sum2=sum2+A0m(iy)*P2(iz)*mesh
          end if
       enddo
       imS(ix)=-Uloc*Uloc*(sum1+sum2)*pi
    enddo
    reS   = kronig(imS,wr,L)
    sigma = reS + xi*imS
  end subroutine solve_ipt_real


  !PURPOSE  : Get auxiliary array Aplus, Aminus for faster polarization evaluation
  subroutine getAs()
    real(8) :: dos(L)
    dos(:) =-dimag(fg0(:))/pi
    A0p(:) = dos(:)*fermi(wr(:),beta)
    A0m(:) = dos(:)*(1.d0-fermi(wr(:),beta))
  end subroutine getAs


  !PURPOSE  : Get polarization bubbles
  subroutine getPolarization
    integer :: ix,iy,iz    
    P1=zero
    P2=zero
    do ix=1,L
       do iy=1,L
          iz= iy_m_ix(iy,ix)
          if(iz>0)then
             P1(ix)=P1(ix) + A0m(iy)*A0p(iz)*mesh
             P2(ix)=P2(ix) + A0p(iy)*A0m(iz)*mesh
          endif
       enddo
    enddo
  end subroutine getPolarization








  



  !+-------------------------------------------------------------------+
  ! ANCILLARY PROCEDURES:
  !+-------------------------------------------------------------------+
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

  !PURPOSE  : get the hilber transfom of a given "zeta" with bethe dos
  function gfbether(w,zeta,d)
    real(8)               :: w,d
    complex(8)            :: zeta
    complex(8)            :: gfbether,sqroot
    real(8)               :: sig
    if(dreal(zeta)==0.d0)zeta=dcmplx(1.d-8,dimag(zeta))
    sqroot=sqrt(zeta**2-d**2)
    sig=dreal(zeta)/abs(dreal(zeta))
    gfbether=2.d0/(zeta+sig*sqroot)
  end function gfbether

  !PURPOSE  : Evaluate a safe Fermi function
  elemental function fermi(x,beta)
    real(8),intent(in) :: x, beta 
    real(8)            :: fermi
    if(x*beta > 100d0)then
       fermi=0d0
       return
    endif
    fermi = 1d0/(1d0+exp(beta*x))
  end function fermi

  !PURPOSE  : Perform a fast Kramers-K\"onig integration: 
  function kronig(fi,wr,M) result(fr)
    integer              :: i,j,M
    real(8),dimension(M) :: fi,wr,fr
    real(8),dimension(M) :: logo,deriv
    real(8)              :: dh,sum
    dh=wr(2)-wr(1)
    logo=0.d0
    do i=2,M-1
       logo(i) = log( (wr(M)-wr(i))/(wr(i)-wr(1)) )
    enddo
    deriv(1)= (fi(2)-fi(1))/dh
    deriv(M)= (fi(M)-fi(M-1))/dh
    do i=2,M-1
       deriv(i) = (fi(i+1)-fi(i-1))/(2*dh)
    enddo
    fr=0.d0
    do i=1,M
       sum=0.d0
       do j=1,M
          if(i/=j)then
             sum=sum+(fi(j)-fi(i))*dh/(wr(j)-wr(i))
          else
             sum=sum+deriv(i)*dh
          endif
       enddo
       fr(i) = (sum + fi(i)*logo(i))/pi
    enddo
    return
  end function kronig

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





end program hmipt_realaxis
