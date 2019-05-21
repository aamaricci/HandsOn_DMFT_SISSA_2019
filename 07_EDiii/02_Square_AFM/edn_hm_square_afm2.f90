!   Solve the Hubbard model with AFM 2 atoms in the basis 
program ed_hm_square_afm2
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                       :: ip,iloop,ilat,ineq,Lk,Nso,Nlso,ispin,iorb
  logical                                       :: converged
  integer                                       :: Nineq,Nlat
  !Bath:
  integer                                       :: Nb
  real(8),allocatable                           :: Bath_ineq(:,:),Bath_prev(:,:)

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal,Greal_ineq
  !Hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk ![Nlat*Nspin*Norb,Nlat*Nspin*Norb,Nk]
  complex(8),allocatable,dimension(:,:)         :: modelHloc ![Nlat*Nspin*Norb,Nlat*Nspin*Norb]
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc_ineq
  real(8),allocatable,dimension(:)              :: Wtk
  !variables for the model:
  character(len=16)                             :: finput
  real(8)                                       :: ts,wmixing
  integer                                       :: Nktot,Nkx,Nkpath
  logical                                       :: spinsym,neelsym

  integer                                       :: comm,rank
  logical                                       :: master


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput   , "FINPUT" , default='inputED.conf')
  call parse_input_variable(ts     , "TS"     , finput, default=1.d0)
  call parse_input_variable(nkx    , "NKX"    , finput, default=25)
  call parse_input_variable(nkpath , "NKPATH" , finput, default=500)
  call parse_input_variable(wmixing, "WMIXING", finput, default=0.75d0)
  call parse_input_variable(spinsym, "SPINSYM", finput, default=.false.)
  call parse_input_variable(neelsym, "NEELSYM", finput, default=.true.)

  call ed_read_input(trim(finput),comm)
  call set_store_size(1024)

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  Nlat=2
  Nineq=2
  Nso=Nspin*Norb
  Nlso=Nlat*Nso
  if(Norb/=1)stop  "Norb != 1"
  if(Nspin/=2)stop "Nspin != 2"

  if(neelsym)then
     Nineq=1
     write(*,*)"Using Neel symmetry to refold BZ"
     write(*,*)"Using Nineq sites=",Nineq
     write(*,*)"Symmetries used are:"
     write(*,*)"(site=2,l,s)=(site=1,l,-s)"
  endif

  if(spinsym)sb_field=0.d0

  !Allocate Weiss Field:
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb))

  !< Build H(k) and H_loc=sum_k H(k)
  call build_hk("Hk_2d_afm2.dat")
  Hloc = lso2nnn_reshape(modelHloc,Nlat,Nspin,Norb)
  do ip=1,Nineq
     Hloc_ineq(ip,:,:,:,:) = Hloc(ip,:,:,:,:)
  enddo

  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath_ineq(Nineq,Nb))
  allocate(Bath_prev(Nineq,Nb))
  call ed_init_solver(comm,Bath_ineq,Hloc_ineq)
  do ip=1,Nineq
     call break_symmetry_bath(Bath_ineq(ip,:),sb_field,(-1d0)**(ip+1))
  enddo


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !solve the impurity problem:
     call ed_solve(Comm,Bath_ineq,Hloc_ineq)
     !
     !retrieve inequivalent self-energies:
     call ed_get_sigma_matsubara(Smats_ineq,Nineq)
     call ed_get_sigma_real(Sreal_ineq,Nineq)
     !
     !extend them to the lattice using symmetry if this applies
     do ip=1,Nineq
        Smats(ip,:,:,:,:,:) = Smats_ineq(ip,:,:,:,:,:)
        Sreal(ip,:,:,:,:,:) = Sreal_ineq(ip,:,:,:,:,:)
     enddo
     if(neelsym)then
        do ispin=1,2
           Smats(2,ispin,ispin,:,:,:)=Smats(1,3-ispin,3-ispin,:,:,:)
           Sreal(2,ispin,ispin,:,:,:)=Sreal(1,3-ispin,3-ispin,:,:,:)
        enddo
     endif
     !
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Comm,Hk,Wtk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
     !
     !fold to the inequivalent sites
     do ip=1,Nineq
        Gmats_ineq(ip,:,:,:,:,:) = Gmats(ip,:,:,:,:,:)
     enddo
     !
     !
     !
     ! compute the Weiss field (only the Nineq ones)
     call dmft_self_consistency(Comm,Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,SCtype=cg_scheme)
     !
     !
     ! fit baths and mix result with old baths
     call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=1)
     if(spinsym)then
        call spin_symmetrize_bath(Bath_ineq,save=.true.)
     else
        call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=2)
     endif
     !
     !
     ! Mixing:
     if(iloop>1)Bath_ineq = wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq

     ! Convergence
     if(master)converged = check_convergence(Weiss_ineq(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     call Bcast_MPI(Comm,converged)
     !
     if(master)call end_loop
  enddo


  call dmft_gloc_realaxis(Comm,Hk,Wtk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)




contains



  !--------------------------------------------------------------------!
  !PURPOSE: BUILD THE H(k) FOR THE BHZ-AFM MODEL.
  !--------------------------------------------------------------------!
  subroutine build_hk(file)
    character(len=*)                        :: file
    integer                                 :: Npts
    integer                                 :: i,j,k,ik,iorb,jorb
    integer                                 :: ix,iy,iz
    real(8)                                 :: kx,ky,kz
    real(8),dimension(:),allocatable        :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable      :: kpath
    real(8),dimension(2)                    :: bk1,bk2,kvec
    real(8)                                 :: n(Nlso)
    complex(8)                              :: w
    complex(8)                              :: Gmats(Nlso,Nlso,Lmats),Greal(Nlso,Nlso,Lreal)
    complex(8)                              :: Smats(Nlso,Nlso,Lmats),Sreal(Nlso,Nlso,Lreal)
    !
    Nktot=Nkx*Nkx
    write(LOGfile,*)"Using Nk_total="//txtfy(Nktot)
    !
    !
    !>Reciprocal lattice basis vector  
    bk1=  pi*[ 1d0, -1d0 ]
    bk2=2*pi*[ 0d0,  1d0 ]
    call TB_set_bk(bk1,bk2)
    !
    !
    !>Get TB Hamiltonian matrix
    allocate(Hk(Nlso,Nlso,Nktot))
    allocate(Wtk(Nktot))
    allocate(modelHloc(Nlso,Nlso))
    call TB_build_model(Hk,hk_model,Nlso,[Nkx,Nkx])
    Wtk=1.d0/dble(Nktot)
    modelHloc = sum(Hk(:,:,:),dim=3)/Nktot
    where(abs(dreal(modelHloc))<1.d-9)modelHloc=0.d0
    call TB_write_Hloc(modelHloc,"Hloc.dat")
    !
    !
    !solve along the standard path in the 2D BZ.
    Npts=4
    allocate(kpath(Npts,3))
    kpath(1,:)=kpoint_Gamma
    kpath(2,:)=kpoint_X1
    kpath(3,:)=kpoint_M1
    kpath(4,:)=kpoint_Gamma
    call TB_solve_model(Hk_model,Nlso,kpath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=20) :: 'G', 'X', 'M', 'G'],&
         file="Eigenbands_afm2.nint")
  end subroutine build_hk




  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(N,N)     :: hk
    complex(8),dimension(N,N)     :: h0,tk
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = kpoint(1)
    ky = kpoint(2)
    !
    ! Hk =  -t * | 0                  1 + e^ikx(e^ikx + e^iky) |
    !            | 1 + e^-ikx(e^-ikx + e^-iky)   0             |
    !
    hk=zero
    hk(1,2) = -ts*(one+exp(xi*2*kx)+exp(xi*(kx+ky))+exp(xi*(kx-ky)))
    hk(2,1) = -ts*(one+exp(-xi*2*kx)+exp(-xi*(kx+ky))+exp(-xi*(kx-ky)))
  end function hk_model



end program ed_hm_square_afm2



