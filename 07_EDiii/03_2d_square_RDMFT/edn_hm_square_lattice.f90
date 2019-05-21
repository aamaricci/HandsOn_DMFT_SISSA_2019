program ed_hm_square_lattice
  USE DMFT_ED
  !
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                                       :: iloop,Nb,Lk,Nx,Nso,ip,ik,iq,Nlat,Nineq,ilat,jlat,iw
  logical                                       :: converged
  real(8)                                       :: wband,ts,wmixing
  !Bath:
  real(8),allocatable                           :: Bath_ineq(:,:),Bath_prev(:,:)
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc(:,:,:,:,:),Hloc_ineq(:,:,:,:,:)
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Gmats,Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Greal,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gijmats
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Weiss
  !
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Gmats_ineq,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Greal_ineq,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Weiss_ineq
  character(len=16)                             :: finput
  complex(8),allocatable                        :: Hij(:,:,:) 
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=2.0d0,comment="hopping parameter")
  call parse_input_variable(Nx,"Nx",finput,default=10,comment="Number of points for each side of the lattice")
  !
  call ed_read_input(trim(finput))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  Nso=1

  !The number of sites in the lattice
  Nlat = Nx*Nx
  print*,"Nlat=",Nlat
  Nineq = 1                     !there is only one inequivalent sites --> all Sigmas are identical


  !Allocate Local Fields for all lattice sites:
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  !for the inequivalent sites (1)
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb))



  ! GET THE TIGHT BINDING HAMILTONIAN FOR THE SQUARE LATTICE 
  allocate(Hij(Nlat,Nlat,1))
  Hij(:,:,1) = one*Htb_square_lattice(Nrow=Nx,Ncol=Nx,ts=ts)

  Hloc = zero
  Hloc_ineq = zero


  Smats = zero
  Smats_ineq = zero


  !setup solver
  Nb=get_bath_dimension()
  allocate(bath_ineq(Nineq,Nb))
  allocate(bath_prev(Nineq,Nb))
  !
  call ed_init_solver(Bath_ineq,Hloc_ineq)



  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(Bath_ineq,Hloc_ineq) 
     call ed_get_sigma_matsubara(Smats_ineq,Nineq)
     call ed_get_sigma_real(Sreal_ineq,Nineq)
     Smats=zero
     Sreal=zero
     do ip=1,Nlat
        Smats(ip,:,:,:,:,:) = Smats_ineq(1,:,:,:,:,:)
        Sreal(ip,:,:,:,:,:) = Sreal_ineq(1,:,:,:,:,:)
     enddo


     !Compute the local gfs:
     call dmft_gloc_matsubara(Hij,[1d0],Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
     Gmats_ineq(1,:,:,:,:,:) = Gmats(1,:,:,:,:,:)


     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq)
     else
        call dmft_delta(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq)
     endif

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(bath_ineq,Weiss_ineq,Hloc_ineq)

     !MIXING:
     if(iloop>1)Bath_ineq = wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq

     !Check convergence (if required change chemical potential)
     converged = check_convergence(Weiss_ineq(1,1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)

     call end_loop
  enddo

  !Compute the local gfs:
  call dmft_gloc_realaxis(Hij,[1d0],Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)


  !Get kinetic energy:
  Smats = zero
  call dmft_kinetic_energy(Hij,[1d0],Smats)




contains


  function Htb_square_lattice(Nrow,Ncol,pbc_row,pbc_col,ts) result(H0)
    integer                                :: Nrow
    integer                                :: Ncol
    logical,optional                       :: pbc_row,pbc_col
    logical                                :: pbc_row_,pbc_col_
    real(8),optional                       :: ts
    real(8)                                :: ts_
    real(8),dimension(Nrow*Ncol,Nrow*Ncol) :: H0
    integer                                :: i,jj,row,col,link(4),j,iorb,jorb
    integer                                :: unit
    !
    pbc_row_=.true. ; if(present(pbc_row)) pbc_row_=pbc_row
    pbc_col_=.true. ; if(present(pbc_col)) pbc_col_=pbc_col
    ts_=0.5d0;if(present(ts))ts_=ts
    !
    H0 = 0.d0
    unit=free_unit()
    !+- 2D LATTICE (NROW x NCOL) -+!
    if(Nlat /= Nrow*Ncol) stop "get_lattice_hamiltonian error: Nlat != Nrow*Ncol"
    !THESE ARE STILL GLOBAL VARIABLES...
    do row=0,Nrow-1
       do col=0,Ncol-1
          i=col+ 1 + row*Ncol
          !
          !
          !right hop
          link(1)= i + 1     
          if((col+1)==Ncol) then
             link(1)=0  
             if(pbc_col_)link(1)=1+row*Ncol  
          end if
          !left  hop
          link(3)= i - 1    
          if((col-1)<0)     then
             link(3)=0  
             if(pbc_col_)link(3)=Ncol+row*Ncol
          end if
          !up    hop
          link(2)= i + Ncol 
          if((row+1)==Nrow) then
             link(2)=0  
             if(pbc_row_)link(2)=col+1
          end if
          !down  hop
          link(4)= i - Ncol 
          if((row-1)<0)     then
             link(4)=0  
             if(pbc_row_)link(4)=col+1+(Nrow-1)*Ncol
          end if
          !
          do jj=1,4
             if(link(jj)>0)H0(i,link(jj))=-ts_ !! ts must be negative.
          enddo
          !
       enddo
    enddo
    open(unit,file='Htb_square_lattice.ed')
    do i=1,Nlat
       write(unit,"(5000(F5.2,1x))")(H0(i,j),j=1,Nlat)
    enddo
    close(unit)
  end function Htb_square_lattice







end program ed_hm_square_lattice



