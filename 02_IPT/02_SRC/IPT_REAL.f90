module IPT_REAL_AXIS
  implicit none
  private

  integer                             :: MM
  real(8),dimension(:),allocatable    :: A0m,A0p,P1,P2
  integer,allocatable,dimension(:,:)  :: iy_m_ix
  complex(8),dimension(:),allocatable :: fg0,sigma
  real(8),dimension(:),allocatable    :: wr
  real(8)                             :: mesh

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : interface function for the impurity solver (SOPT)
  ! half-filling case
  !+-------------------------------------------------------------------+
  function ipt_solve_real(fg0_,wr_) result(sigma_)
    complex(8),dimension(:)          :: fg0_
    complex(8),dimension(size(fg0_)) :: sigma_
    real(8),dimension(size(fg0_))    :: wr_
    MM=size(fg0_)
    if(.not.allocated(wr))allocate(wr(MM))
    if(.not.allocated(fg0))allocate(fg0(MM))
    if(.not.allocated(sigma))allocate(sigma(MM))
    if(.not.allocated(iy_m_ix))call get_frequency_index       
    !
    fg0 = fg0_
    wr  = wr_
    mesh= abs(wr(2)-wr(1))
    !
    call getAs
    call getPolarization
    call Sopt
    !
    sigma_=sigma
  end function ipt_solve_real


  !PURPOSE  : Create an array of the indices y-x for a faster evaluation
  subroutine get_frequency_index()
    integer :: ix,iy,iz
    if(.not.allocated(iy_m_ix))allocate(iy_m_ix(MM,MM))
    iy_m_ix=0
    do ix=1,MM
       do iy=1,MM
          iz = iy - ix + MM/2 
          if(iz<1 .OR. iz>MM) iz=-1 !out of range-> if(iz>-L)
          iy_m_ix(iy,ix)=iz
       enddo
    enddo
    if(.not.allocated(A0m))allocate(A0m(MM))
    if(.not.allocated(A0p))allocate(A0p(MM))
    if(.not.allocated(P1)) allocate(P1(MM))
    if(.not.allocated(P2)) allocate(P2(MM))
  end subroutine get_frequency_index


  !PURPOSE  : Get auxiliary array Aplus, Aminus for faster polarization evaluation
  subroutine getAs
    real(8) :: dos(MM)
    dos(:) =-dimag(fg0(:))/pi
    A0p(:) = dos(:)*fermi(wr(:),beta)
    A0m(:) = dos(:)*(1.d0-fermi(wr(:),beta))
  end subroutine getAs

  !PURPOSE  : Get polarization bubbles
  subroutine getPolarization
    integer :: ix,iy,iz    
    P1=zero
    P2=zero
    do ix=1,MM
       do iy=1,MM
          iz= iy_m_ix(iy,ix)
          if(iz>0)then
             P1(ix)=P1(ix) + A0m(iy)*A0p(iz)*mesh
             P2(ix)=P2(ix) + A0p(iy)*A0m(iz)*mesh
          endif
       enddo
    enddo
  end subroutine getPolarization

  !PURPOSE  : Solve 2^nd order perturbation theory
  subroutine Sopt
    integer :: ix,iy,iz
    real(8) :: sum1,sum2
    real(8),dimension(MM) :: reS,imS
    do ix=1,MM
       sum1=zero
       sum2=zero
       do iy=1,MM
          iz= iy_m_ix(iy,ix)
          if(iz>0)then
             sum1=sum1+A0p(iy)*P1(iz)*mesh
             sum2=sum2+A0m(iy)*P2(iz)*mesh
          end if
       enddo
       imS(ix)=-Uloc(1)*Uloc(1)*(sum1+sum2)*pi
    enddo
    reS = kronig(imS,wr,size(ImS))
    sigma = reS + xi*imS
  end subroutine Sopt





  !+-------------------------------------------------------------------+
  ! ANCILLARY PROCEDURES:
  !+-------------------------------------------------------------------+
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


END MODULE IPT_REAL_AXIS




