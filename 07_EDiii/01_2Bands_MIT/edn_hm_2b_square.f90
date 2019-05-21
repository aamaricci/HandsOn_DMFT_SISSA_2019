program ed_hm_2bands
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                       :: iloop,Lk,Nso,Nlso,Nlat
  logical                                       :: converged
  integer                                       :: ispin,ilat!,i,j

  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:)              :: Bath,Bath_prev

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats,Greal

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk
  complex(8),allocatable,dimension(:,:)         :: modelHloc
  complex(8),allocatable,dimension(:,:,:,:)     :: Hloc
  real(8),allocatable,dimension(:)              :: Wtk

  integer,allocatable,dimension(:)              :: ik2ix,ik2iy
  real(8),dimension(2)                          :: bk1,bk2

  !variables for the model:
  integer                                       :: Nk,Nkpath
  real(8)                                       :: ts,wmixing,Mh
  character(len=32)                             :: finput
  character(len=32)                             :: hkfile
  logical                                       :: spinsym
  !


  !Parse additional variables && read Input && read H(k)^2x2
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(mh,"MH",finput,default=0d0)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS","inputED.conf",default=0.5d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call ed_read_input(trim(finput))
  !
  !< Control variables for DMFT_TOOLS
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")


  if(Norb/=2)stop "Wrong setup from input file: Norb=2"
  Nso=Nspin*Norb
  Nlso=Nso



  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  allocate(Hloc(Nspin,Nspin,Norb,Norb));Hloc=zero

  !Build the Hamiltonian on a grid or on a path
  call build_hk(trim(hkfile))
  Hloc = so2nn_reshape(modelHloc,Nspin,Norb)

  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_prev(Nb))
  call ed_init_solver(Bath,Hloc)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(Bath,Hloc)

     call ed_get_sigma_matsubara(Smats)

     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)

     ! compute the Weiss field (only the Nineq ones)
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,SCtype=cg_scheme)

     !Fit the new bath, starting from the old bath + the supplied Weiss
     call ed_chi2_fitgf(Weiss,Bath,ispin=1)

     !MIXING:
     if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop)
     call end_loop
  enddo

  call ed_get_sigma_real(Sreal)
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)



contains




  !--------------------------------------------------------------------!
  !Lattice Hamitonian:
  !--------------------------------------------------------------------!
  function hk_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: kx,ky
    real(8)                         :: dxy
    !
    if(Nlso/=2)stop "hk_model: Nso != 2"
    !
    kx=kpoint(1)
    ky=kpoint(2)
    !
    hk(:,:) = zero
    !
    hk(1,1) = -2d0*ts*(cos(kx)+cos(ky)) + Mh
    hk(2,2) = -2d0*ts*(cos(kx)+cos(ky)) - Mh
    !
  end function hk_model






  !---------------------------------------------------------------------
  !PURPOSE: get model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional             :: file
    integer                               :: i,j,ik
    integer                               :: ix,iy
    real(8)                               :: kx,ky  
    integer                               :: iorb,jorb
    integer                               :: isporb,jsporb
    integer                               :: ispin,jspin
    integer                               :: unit
    complex(8),dimension(Nlso,Nlso,Lmats) :: Gmats,fooSmats
    complex(8),dimension(Nlso,Nlso,Lreal) :: Greal,fooSreal
    real(8),dimension(2)                  :: kvec
    real(8)                               :: blen,area_hex,area_rect,points_in,points_tot
    real(8),allocatable,dimension(:)      :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable    :: kpath

    Lk= Nk*Nk

    write(LOGfile,*)"Build H(k)    :",Lk
    write(LOGfile,*)"# of SO-bands :",Nlso

    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    allocate(wtk(Lk));Wtk=0d0

    call TB_set_bk([pi2,0d0],[0d0,pi2])

    call TB_build_model(Hk, hk_model, Nlso, [Nk,Nk])

    Wtk = 1d0/Lk

    if(present(file))&
         call TB_write_hk(Hk, trim(file), &
         No = Nlso,Nd = Norb,Np = 0,Nineq = 1,&
         Nkvec=[Nk,Nk])
    !
    allocate(modelHloc(Nlso,Nlso))
    modelHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(modelHloc))<1.d-4)modelHloc=0d0

    !path: G X M G
    allocate(kpath(4,2))
    kpath(1,:)=[0d0,0d0]
    kpath(2,:)=[ pi,0d0]
    kpath(3,:)=[ pi, pi]
    kpath(4,:)=[0d0,0d0]
    call TB_solve_model(hk_model,Nlso,kpath,Nkpath,&
         colors_name=[red1,green1,blue1],&
         points_name=[character(len=10) :: "G","X","M", "G"],&
         file="Eigenbands.nint")
    !
  end subroutine build_hk








end program ed_hm_2bands



