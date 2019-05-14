! > BUILD SPARSE HAMILTONIAN of the SECTOR
MODULE ED_HAMILTONIAN_SPARSE_HxV
  USE ED_HAMILTONIAN_COMMON
  implicit none
  private


  !>Sparse Matric constructors
  public :: ed_buildh_main
  public :: ed_buildh_orbs

  !>Sparse Mat-Vec product using stored sparse matrix
  public  :: spMatVec_main
  public  :: spMatVec_orbs
#ifdef _MPI
  public  :: spMatVec_MPI_main
  public  :: spMatVec_MPI_orbs
#endif


  integer                              :: i,iup,idw
  integer                              :: j,jup,jdw
  integer                              :: m,mup,mdw
  integer                              :: ms,iud
  integer                              :: impi
  integer                              :: iorb,jorb,ispin,jspin,ibath
  integer                              :: kp,k1,k2,k3,k4
  integer                              :: ialfa,ibeta,indx
  real(8)                              :: sg1,sg2,sg3,sg4
  real(8)                              :: htmp,htmpup,htmpdw
  logical                              :: Jcondition
  integer                              :: Nfoo
  real(8),dimension(:,:,:),allocatable :: diag_hybr ![Nspin,Norb,Nbath]
  real(8),dimension(:,:,:),allocatable :: bath_diag ![Nspin,Norb/1,Nbath]


contains


  subroutine ed_buildh_main(isector,Hmat)
    integer                              :: isector   
    real(8),dimension(:,:),optional      :: Hmat
    real(8),dimension(:,:),allocatable   :: Htmp_up,Htmp_dw,Hrdx
    integer,dimension(2*Ns_Ud)           :: Indices    ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)      :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                :: Nup,Ndw    ![Ns]
    !
#ifdef _MPI
    if(Mpistatus .AND. MpiComm == MPI_COMM_NULL)return
#endif
    !
    if(.not.Hstatus)stop "ed_buildh_main ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(present(Hmat))&
         call assert_shape(Hmat,[getdim(isector), getdim(isector)],"ed_buildh_main","Hmat")
    !
    !Get diagonal hybridization, bath energy
    include "ED_HAMILTONIAN/diag_hybr_bath.f90"
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0d,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0d,Dim)
       if(Jhflag)then
          call sp_set_mpi_matrix(MpiComm,spH0nd,mpiIstart,mpiIend,mpiIshift)
          call sp_init_matrix(MpiComm,spH0nd,Dim)
       endif
    else
       call sp_init_matrix(spH0d,Dim)
       if(Jhflag)call sp_init_matrix(spH0nd,Dim)
    endif
#else
    call sp_init_matrix(spH0d,Dim)
    if(Jhflag)call sp_init_matrix(spH0nd,Dim)
#endif
    call sp_init_matrix(spH0dws(1),DimDw)
    call sp_init_matrix(spH0ups(1),DimUp)
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN TERMS
    include "ED_HAMILTONIAN/stored/H_local.f90"
    !
    !NON-LOCAL HAMILTONIAN TERMS
    if(jhflag)then
       include "ED_HAMILTONIAN/stored/H_non_local.f90"
    endif
    !
    !UP TERMS
    include "ED_HAMILTONIAN/stored/H_up.f90"
    !
    !DW TERMS
    include "ED_HAMILTONIAN/stored/H_dw.f90"
    !-----------------------------------------------!
    !
    if(present(Hmat))then
       Hmat = 0d0
       allocate(Htmp_up(DimUp,DimUp));Htmp_up=0d0
       allocate(Htmp_dw(DimDw,DimDw));Htmp_dw=0d0
       !
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0d,Hmat)
       else
          call sp_dump_matrix(spH0d,Hmat)
       endif
#else
       call sp_dump_matrix(spH0d,Hmat)
#endif
       !
       if(Jhflag)then
          allocate(Hrdx(Dim,Dim));Hrdx=0d0
#ifdef _MPI
          if(MpiStatus)then
             call sp_dump_matrix(MpiComm,spH0nd,Hrdx)
          else
             call sp_dump_matrix(spH0nd,Hrdx)
          endif
#else
          call sp_dump_matrix(spH0nd,Hrdx)
#endif
          Hmat = Hmat + Hrdx
          deallocate(Hrdx)
       endif
       !
       call sp_dump_matrix(spH0ups(1),Htmp_up)
       call sp_dump_matrix(spH0dws(1),Htmp_dw)
       Hmat = Hmat + kronecker_product(Htmp_dw,eye(DimUp))
       Hmat = Hmat + kronecker_product(eye(DimDw),Htmp_up)
       !
       deallocate(Htmp_up,Htmp_dw)
    endif
    !
    deallocate(diag_hybr,bath_diag)
    return
    !
  end subroutine ed_buildh_main





  subroutine ed_buildh_orbs(isector,Hmat)
    integer                             :: isector   
    real(8),dimension(:,:),optional     :: Hmat
    real(8),dimension(:,:),allocatable  :: Htmp_up,Htmp_dw,Hrdx
    integer,dimension(2*Ns_Ud)          :: Indices    ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)               :: Nup,Ndw    ![Ns]
    !
#ifdef _MPI
    if(Mpistatus .AND. MpiComm == MPI_COMM_NULL)return
#endif
    !
    if(.not.Hstatus)stop "ed_buildh_main ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(present(Hmat))&
         call assert_shape(Hmat,[getdim(isector), getdim(isector)],"ed_buildh_main","Hmat")
    !
    !Get diagonal hybridization, bath energy
    include "ED_HAMILTONIAN/diag_hybr_bath.f90"
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0d,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0d,Dim)
    else
       call sp_init_matrix(spH0d,Dim)
    endif
#else
    call sp_init_matrix(spH0d,Dim)
#endif
    do iud=1,Ns_Ud
       call sp_init_matrix(spH0dws(iud),DimDws(iud))
       call sp_init_matrix(spH0ups(iud),DimUps(iud))
    enddo
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN TERMS
    include "ED_HAMILTONIAN/stored/Orbs/H_local.f90"
    !
    !UP TERMS
    include "ED_HAMILTONIAN/stored/Orbs/H_up.f90"
    !
    !DW TERMS
    include "ED_HAMILTONIAN/stored/Orbs/H_dw.f90"
    !-----------------------------------------------!
    !
    if(present(Hmat))then
       Hmat = 0d0
       allocate(Htmp_up(DimUp,DimUp));Htmp_up=0d0
       allocate(Htmp_dw(DimDw,DimDw));Htmp_dw=0d0
       !
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0d,Hmat)
       else
          call sp_dump_matrix(spH0d,Hmat)
       endif
#else
       call sp_dump_matrix(spH0d,Hmat)
#endif
       do iud=1,Ns_Ud
          allocate(Hrdx(DimDws(iud),DimDws(iud)));Hrdx=0d0
          call sp_dump_matrix(spH0dws(iud),Hrdx)
          call build_Htmp_dw(iud,Hrdx,DimDws(iud),Htmp_dw)
          Hmat = Hmat + kronecker_product(Htmp_dw,eye(DimUp))
          deallocate(Hrdx)
          !
          allocate(Hrdx(DimUps(iud),DimUps(iud)));Hrdx=0d0
          call sp_dump_matrix(spH0ups(iud),Hrdx)
          call build_Htmp_up(iud,Hrdx,DimUps(iud),Htmp_up)
          Hmat = Hmat + kronecker_product(eye(DimDw),Htmp_up)
          deallocate(Hrdx)
          !
       enddo
       !
       deallocate(Htmp_up,Htmp_dw)
    endif
    !
    deallocate(diag_hybr,bath_diag)
    return
    !
  contains
    subroutine build_Htmp_up(iud,H,Dim,Hup)
      integer                        :: iud,Dim,i
      real(8),dimension(Dim,Dim)     :: H
      real(8),dimension(DimUp,DimUp) :: Hup
      if(dim/=DimUps(iud))stop "error in build_Htmp_up"
      if(iud==1)then
         Hup= kronecker_product(H,eye(product(DimUps(2:))))
      else if(iud==Ns_Ud)then
         Hup= kronecker_product(eye(product(DimUps(1:Ns_Ud-1))),H)
      else
         Hup= kronecker_product( &
              kronecker_product( &
              eye(product(DimUps(1:iud-1))), H) , eye(product(DimUps(iud+1:Ns_Ud))) )
      end if
    end subroutine build_Htmp_up
    subroutine build_Htmp_dw(iud,H,Dim,Hdw)
      integer                        :: iud,Dim,i
      real(8),dimension(Dim,Dim)     :: H
      real(8),dimension(DimDw,DimDw) :: Hdw
      if(dim/=DimDws(iud))stop "error in build_Htmp_dw"
      if(iud==1)then
         Hdw= kronecker_product(H,eye(product(DimDws(2:))))
      else if(iud==Ns_Ud)then
         Hdw= kronecker_product(eye(product(DimDws(1:Ns_Ud-1))),H)
      else
         Hdw= kronecker_product( &
              kronecker_product( &
              eye(product(DimDws(1:iud-1))), H) , eye(product(DimDws(iud+1:Ns_Ud))) )
      end if
    end subroutine build_Htmp_dw
  end subroutine ed_buildh_orbs












  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial
  ! - MPI
  !+------------------------------------------------------------------+
  subroutine spMatVec_main(Nloc,v,Hv)
    integer                         :: Nloc
    real(8),dimension(Nloc)         :: v
    real(8),dimension(Nloc)         :: Hv
    real(8)                         :: val
    integer                         :: i,iup,idw,j,jup,jdw,jj
    !
    !
    Hv=0d0
    !
    !Local:
    do i = 1,Nloc
       do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(spH0d%row(i)%cols(j))
       enddo
    enddo
    !
    !DW:
    do iup=1,DimUp
       !
       do idw=1,DimDw
          i = iup + (idw-1)*DimUp
          do jj=1,spH0dws(1)%row(idw)%Size
             jup = iup
             jdw = spH0dws(1)%row(idw)%cols(jj)
             val = spH0dws(1)%row(idw)%vals(jj)
             j     = jup +  (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !
    !UP:
    do idw=1,DimDw
       !
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          do jj=1,spH0ups(1)%row(iup)%Size
             jup = spH0ups(1)%row(iup)%cols(jj)
             jdw = idw
             val = spH0ups(1)%row(iup)%vals(jj)
             j =  jup + (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !
    !Non-Local:
    if(jhflag)then
       do i = 1,Nloc
          do j=1,spH0nd%row(i)%Size
             val   = spH0nd%row(i)%vals(j)
             jj    = spH0nd%row(i)%cols(j)
             Hv(i) = Hv(i) + val*V(jj)
          enddo
       enddo
    endif
    !
  end subroutine spMatVec_main

  subroutine spMatVec_orbs(Nloc,v,Hv)
    integer                    :: Nloc
    real(8),dimension(Nloc)    :: v
    real(8),dimension(Nloc)    :: Hv
    integer                    :: i,iup,idw,j,jup,jdw,jj
    integer                    :: iud
    integer,dimension(2*Ns_Ud) :: Indices,Jndices
    !
    !
    Hv=0d0
    !
    do i = 1,Nloc
       do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(spH0d%row(i)%cols(j))
       enddo
    enddo
    !
    !
    do i=1,Dim
       call state2indices(i,[DimUps,DimDws],Indices)
       do iud=1,Ns_Ud
          !
          !UP:
          iup = Indices(iud)
          do jj=1,spH0ups(iud)%row(iup)%Size
             Jndices = Indices ; Jndices(iud) = spH0ups(iud)%row(iup)%cols(jj)
             call indices2state(Jndices,[DimUps,DimDws],j)
             Hv(i) = Hv(i) + spH0ups(iud)%row(iup)%vals(jj)*V(j)
          enddo
          !
          !DW:
          idw = Indices(iud+Ns_Ud)
          do jj=1,spH0dws(iud)%row(idw)%Size
             Jndices = Indices ; Jndices(iud+Ns_Ud) = spH0dws(iud)%row(idw)%cols(jj)
             call indices2state(Jndices,[DimUps,DimDws],j)
             Hv(i) = Hv(i) + spH0dws(iud)%row(idw)%vals(jj)*V(j)
          enddo
          !
       enddo
    enddo
    !
  end subroutine spMatVec_orbs


#ifdef _MPI
  subroutine spMatVec_mpi_main(Nloc,v,Hv)
    integer                          :: Nloc
    real(8),dimension(Nloc)          :: v
    real(8),dimension(Nloc)          :: Hv
    !
    integer                          :: N
    real(8),dimension(:),allocatable :: vt,Hvt
    real(8),dimension(:),allocatable :: vin
    real(8)                          :: val
    integer                          :: i,iup,idw,j,jup,jdw,jj
    !local MPI
    integer                          :: irank,MpiIerr
    integer,allocatable,dimension(:) :: Counts
    integer,allocatable,dimension(:) :: Offset
    !
    if(MpiComm==Mpi_Comm_Null)return
    if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    Hv=0d0
    do i=1,Nloc                 !==spH0%Nrow
       do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(i)
       end do
    end do
    !
    !Non-local terms.
    !UP part: contiguous in memory.
    do idw=1,MpiQdw
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          hxv_up: do jj=1,spH0ups(1)%row(iup)%Size
             jup = spH0ups(1)%row(iup)%cols(jj)
             jdw = idw
             val = spH0ups(1)%row(iup)%vals(jj)
             j   = jup + (idw-1)*DimUp
             Hv(i) = Hv(i) + val*v(j)
          end do hxv_up
       enddo
    end do
    !
    !DW part: non-contiguous in memory -> MPI transposition
    !Transpose the input vector as a whole:
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    !
    allocate(vt(mpiQup*DimDw)) ;vt=0d0
    allocate(Hvt(mpiQup*DimDw));Hvt=0d0
    call vector_transpose_MPI(DimUp,MpiQdw,v,DimDw,MpiQup,vt)
    Hvt=0d0    
    do idw=1,MpiQup             !<= Transposed order:  column-wise DW <--> UP  
       do iup=1,DimDw           !<= Transposed order:  column-wise DW <--> UP
          i = iup + (idw-1)*DimDw
          hxv_dw: do jj=1,spH0dws(1)%row(iup)%Size
             jup = spH0dws(1)%row(iup)%cols(jj)
             jdw = idw             
             j   = jup + (jdw-1)*DimDw
             val = spH0dws(1)%row(iup)%vals(jj)
             Hvt(i) = Hvt(i) + val*vt(j)
          end do hxv_dw
       enddo
    end do
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ; vt=0d0
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt)
    Hv = Hv + Vt
    deallocate(vt)
    !
    !
    !Non-Local:
    if(jhflag)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       ! 
       allocate(vt(N)) ; vt = 0d0
       call allgather_vector_MPI(MpiComm,v,vt)
       !
       do i=1,Nloc
          matmul: do j=1,spH0nd%row(i)%Size
             Hv(i) = Hv(i) + spH0nd%row(i)%vals(j)*Vt(spH0nd%row(i)%cols(j))
          enddo matmul
       enddo
       deallocate(Vt)
    endif
    !
  end subroutine spMatVec_mpi_main


  subroutine spMatVec_mpi_orbs(Nloc,v,Hv)
    integer                          :: Nloc
    real(8),dimension(Nloc)          :: v
    real(8),dimension(Nloc)          :: Hv
    !
    integer                          :: N
    real(8),dimension(:),allocatable :: vt,Hvt
    real(8)                          :: val
    integer                          :: i,iup,idw,j,jup,jdw,jj
    integer                          :: iiup,iidw
    integer                          :: iud
    integer,dimension(2*Ns_Ud)       :: Indices,Jndices
    !local MPI
    integer                          :: irank
    !
    if(MpiComm==Mpi_Comm_Null)return
    if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    Hv=0d0
    do i=1,Nloc                 !==spH0%Nrow
       do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(i)
       end do
    end do
    !
    !
    !Non-local terms.
    !UP part: contiguous in memory.
    do iidw=1,MpiQdw
       do iiup=1,DimUp
          i = iiup + (iidw-1)*DimUp
          call state2indices(i,[DimUps,DimDws],Indices)
          do iud=1,Ns_Ud
             !
             iup = Indices(iud)
             hxv_up: do jj=1,spH0ups(iud)%row(iup)%Size
                Jndices      = Indices
                Jndices(iud) = spH0ups(iud)%row(iup)%cols(jj)
                call indices2state(Jndices,[DimUps,DimDws],j)
                Hv(i) = Hv(i) + spH0ups(iud)%row(iup)%vals(jj)*v(j)
             end do hxv_up
             !
          enddo
       enddo
    end do
    !
    !DW part: non-contiguous in memory -> MPI transposition
    !Transpose the input vector as a whole:
    !
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    !
    allocate(vt(mpiQup*DimDw)) ;vt=0d0
    allocate(Hvt(mpiQup*DimDw));Hvt=0d0
    call vector_transpose_MPI(DimUp,MpiQdw,v,DimDw,MpiQup,vt)
    Hvt=0d0    
    do iidw=1,MpiQup            !<= Transposed order:  column-wise DW <--> UP  
       do iiup=1,DimDw          !<= Transposed order:  column-wise DW <--> UP  
          i = iiup + (iidw-1)*DimDw
          call state2indices(i,[DimDws,DimUps],Indices)
          do iud=1,Ns_Ud
             !
             iup = Indices(iud)
             hxv_dw: do jj=1,spH0dws(iud)%row(iup)%Size
                Jndices      = Indices
                Jndices(iud) = spH0dws(iud)%row(iup)%cols(jj)
                call indices2state(Jndices,[DimDws,DimUps],j)
                Hvt(i) = Hvt(i) + spH0dws(iud)%row(iup)%vals(jj)*vt(j)
             end do hxv_dw
             !
          enddo
       enddo
    end do
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ; vt=0d0
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt)
    Hv = Hv + vt
  end subroutine spMatVec_mpi_orbs
#endif



end MODULE ED_HAMILTONIAN_SPARSE_HXV
