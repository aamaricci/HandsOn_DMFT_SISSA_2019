MODULE ED_HAMILTONIAN_COMMON
  USE SF_MISC,    only: assert_shape
  USE SF_LINALG,  only: kronecker_product,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_SETUP
  implicit none


  !> MPI local variables (shared)
  ! #ifdef _MPI
  !   integer                                   :: MpiComm=MPI_UNDEFINED
  ! #else
  !   integer                                   :: MpiComm=0
  ! #endif
  !   logical                                   :: MpiStatus=.false.
  !   logical                                   :: MpiMaster=.true.
  !   integer                                   :: MpiIerr
  !   integer                                   :: MpiRank=0
  !   integer                                   :: MpiSize=1
  !   integer                                   :: MpiQ=1
  !   integer                                   :: MpiQup=1
  !   integer                                   :: MpiQdw=1
  !   integer                                   :: MpiR=0
  !   integer                                   :: MpiRup=0
  !   integer                                   :: MpiRdw=0
  !   integer                                   :: MpiIstart
  !   integer                                   :: MpiIend
  !   integer                                   :: MpiIshift
  !
  integer                                   :: Dim
  integer                                   :: DimUp
  integer                                   :: DimDw
  integer,allocatable,dimension(:)          :: DimUps
  integer,allocatable,dimension(:)          :: DimDws
  !
  integer                                   :: Hsector=0
  logical                                   :: Hstatus=.false.
  type(sector_map),dimension(:),allocatable :: Hs



contains



  !####################################################################
  !               ALL-2-ALL-V VECTOR MPI TRANSPOSITION 
  !####################################################################
#ifdef _MPI
  subroutine vector_transpose_MPI(nrow,qcol,a,ncol,qrow,b)    
    integer                            :: nrow,ncol,qrow,qcol
    real(8)                            :: a(nrow,qcol)
    real(8)                            :: b(ncol,qrow)
    integer,allocatable,dimension(:,:) :: send_counts,send_offset
    integer,allocatable,dimension(:,:) :: recv_counts,recv_offset
    integer                            :: counts,Ntot
    integer :: i,j,irank,ierr
    !
    counts = Nrow/MpiSize
    Ntot   = Ncol/MpiSize
    if(mod(Ncol,MpiSize)/=0)Ntot=Ntot+1
    !
    allocate(send_counts(0:MpiSize-1,Ntot));send_counts=0
    allocate(send_offset(0:MpiSize-1,Ntot));send_offset=0
    allocate(recv_counts(0:MpiSize-1,Ntot));recv_counts=0
    allocate(recv_offset(0:MpiSize-1,Ntot));recv_offset=0
    !
    do i=1,qcol
       do irank=0,MpiSize-1
          if(irank < mod(Nrow,MpiSize))then
             send_counts(irank,i) = counts+1
          else
             send_counts(irank,i) = counts
          endif
       enddo
    enddo
    !
    do i=1,Ntot
       call MPI_AllToAll(&
            send_counts(:,i),1,MPI_INTEGER,&
            recv_counts(:,i),1,MPI_INTEGER,&
            MpiComm,ierr)
    enddo
    !
    do i=1,Ntot
       do irank=1,MpiSize-1
          send_offset(irank,i) = send_counts(irank-1,i) + send_offset(irank-1,i)
       enddo
    enddo
    !
    !Get the irank=0 elements, i.e. first entries:
    recv_offset(0,1) = 0
    do i=2,Ntot
       recv_offset(0,i) = sum(recv_counts(0,:i-1))
    enddo
    !the rest of the entries:
    do i=1,Ntot
       do irank=1,MpiSize-1
          recv_offset(irank,i) = recv_offset(irank-1,i) + sum(recv_counts(irank-1,:))
       enddo
    enddo
    !
    !
    do j=1,Ntot
       call MPI_AllToAllV(&
            A(:,j),send_counts(:,j),send_offset(:,j),MPI_DOUBLE_PRECISION,&
            B(:,:),recv_counts(:,j),recv_offset(:,j),MPI_DOUBLE_PRECISION,&
            MpiComm,ierr)
    enddo
    !
    call local_transpose(b,ncol,qrow)
    !
    return
  end subroutine vector_transpose_MPI


  subroutine local_transpose(mat,nrow,ncol)
    integer                      :: nrow,ncol
    real(8),dimension(Nrow,Ncol) :: mat
    mat = transpose(reshape(mat,[Ncol,Nrow]))
  end subroutine local_transpose
#endif



end MODULE ED_HAMILTONIAN_COMMON




