! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
MODULE ED_HAMILTONIAN_DIRECT_HxV
  USE ED_HAMILTONIAN_COMMON
  implicit none
  private

  integer                              :: iiup,iidw
  integer                              :: iud,jj
  integer                              :: i,iup,idw
  integer                              :: j,jup,jdw
  integer                              :: m,mup,mdw
  integer                              :: ishift
  integer                              :: isector,jsector
  integer                              :: ms
  integer                              :: impi
  integer                              :: iorb,jorb,ispin,jspin,ibath
  integer                              :: kp,k1,k2,k3,k4
  integer                              :: ialfa,ibeta
  real(8)                              :: sg1,sg2,sg3,sg4
  real(8)                              :: htmp,htmpup,htmpdw
  logical                              :: Jcondition
  integer                              :: Nfoo
  real(8),dimension(:,:,:),allocatable :: diag_hybr ![Nspin,Norb,Nbath]
  real(8),dimension(:,:,:),allocatable :: bath_diag ![Nspin,Norb/1,Nbath]




  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_main
  public  :: directMatVec_orbs
#ifdef _MPI
  public  :: directMatVec_MPI_main
  public  :: directMatVec_MPI_orbs
#endif



contains


  subroutine directMatVec_main(Nloc,vin,Hv)
    integer                             :: Nloc
    real(8),dimension(Nloc)             :: vin
    real(8),dimension(Nloc)             :: Hv
    real(8),dimension(:),allocatable    :: vt,Hvt
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)               :: Nup,Ndw
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization, bath energy
    include "ED_HAMILTONIAN/diag_hybr_bath.f90"
    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "ED_HAMILTONIAN/direct/HxV_local.f90"
    !
    !UP HAMILTONIAN TERMS
    include "ED_HAMILTONIAN/direct/HxV_up.f90"
    !    
    !DW HAMILTONIAN TERMS
    include "ED_HAMILTONIAN/direct/HxV_dw.f90"
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       include "ED_HAMILTONIAN/direct/HxV_non_local.f90"
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_main




  subroutine directMatVec_orbs(Nloc,vin,Hv)
    integer                             :: Nloc
    real(8),dimension(Nloc)             :: vin
    real(8),dimension(Nloc)             :: Hv
    real(8),dimension(:),allocatable    :: vt,Hvt
    integer                             :: isector
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)               :: Nup,Ndw
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization, bath energy
    include "ED_HAMILTONIAN/diag_hybr_bath.f90"
    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "ED_HAMILTONIAN/direct/Orbs/HxV_local.f90" 
    !
    !UP HAMILTONIAN TERMS
    include "ED_HAMILTONIAN/direct/Orbs/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS
    include "ED_HAMILTONIAN/direct/Orbs/HxV_dw.f90"
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_orbs







#ifdef _MPI
  subroutine directMatVec_MPI_main(Nloc,vin,Hv)
    integer                             :: Nloc,N
    real(8),dimension(Nloc)             :: Vin
    real(8),dimension(Nloc)             :: Hv
    real(8),dimension(:),allocatable    :: vt,Hvt
    !
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)               :: Nup,Ndw
    !
    integer                             :: MpiIerr
    integer,allocatable,dimension(:)    :: Counts
    integer,allocatable,dimension(:)    :: Offset
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    !Get diagonal hybridization, bath energy
    include "ED_HAMILTONIAN/diag_hybr_bath.f90"
    !
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    ! if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "ED_HAMILTONIAN/direct_mpi/HxV_local.f90"
    !
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "ED_HAMILTONIAN/direct_mpi/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    allocate(vt(mpiQup*DimDw)) ;vt=0d0
    allocate(Hvt(mpiQup*DimDw));Hvt=0d0
    call vector_transpose_MPI(DimUp,MpiQdw,Vin,DimDw,MpiQup,vt) !Vin^T --> Vt
    include "ED_HAMILTONIAN/direct_mpi/HxV_dw.f90"
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=0d0         !reallocate Vt
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt) !Hvt^T --> Vt
    Hv = Hv + Vt
    deallocate(vt)
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       !
       allocate(vt(N)) ; vt = 0d0
       call allgather_vector_MPI(MpiComm,vin,vt)
       !
       include "ED_HAMILTONIAN/direct_mpi/HxV_non_local.f90"
       !
       deallocate(Vt)
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_MPI_main



  subroutine directMatVec_MPI_Orbs(Nloc,vin,Hv)
    integer                             :: Nloc
    real(8),dimension(Nloc)             :: Vin
    real(8),dimension(Nloc)             :: Hv
    real(8),dimension(:),allocatable    :: vt,Hvt
    !
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)               :: Nup,Ndw
    integer,allocatable,dimension(:)    :: Counts,Displs
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    !Get diagonal hybridization, bath energy
    include "ED_HAMILTONIAN/diag_hybr_bath.f90"
    !
    !    
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    !
    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "ED_HAMILTONIAN/direct_mpi/Orbs/HxV_local.f90"
    !
    !NON-LOCAL TERMS:
    ! 
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "ED_HAMILTONIAN/direct_mpi/Orbs/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    allocate(vt(mpiQup*DimDw)) ;vt=0d0
    allocate(Hvt(mpiQup*DimDw));Hvt=0d0
    call vector_transpose_MPI(DimUp,MpiQdw,Vin,DimDw,MpiQup,vt) !Vin^T --> Vt
    include "ED_HAMILTONIAN/direct_mpi/Orbs/HxV_dw.f90"
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=0d0         !reallocate Vt
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt) !Hvt^T --> Vt
    Hv = Hv + Vt
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_MPI_Orbs
#endif


END MODULE ED_HAMILTONIAN_DIRECT_HXV
