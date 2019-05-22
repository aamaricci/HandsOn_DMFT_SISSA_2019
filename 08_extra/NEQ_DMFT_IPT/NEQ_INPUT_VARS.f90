MODULE NEQ_INPUT_VARS
  USE SF_VERSION
  USE SF_CONSTANTS
  USE SF_PARSE_INPUT
  implicit none

  !GIT VERSION
  include "revision.inc"

  !Gloabl  variables
  integer                 :: Norb          !Norb =# of impurity orbitals
  integer                 :: Nspin         !Nspin=# spin degeneracy (max 2)
  integer                 :: Ntime         !Number of Time steps
  integer                 :: Ntau          !Imaginary time slices
  integer                 :: Niw           !Number of frequencies
  real(8)                 :: dt            !real-time step
  real(8)                 :: dtau          !imag-time step
  real(8)                 :: beta          !inverse temperature
  integer                 :: nloop         !dmft loop variables
  real(8)                 :: Ui            !equilibrium local interaction
  real(8)                 :: U             !non-equilibrium local interaction
  real(8)                 :: xmu           !chemical potential (shifted by Uloc/2)
  real(8)                 :: Vbath         !coupling to the Thermostat
  integer                 :: Lbath         !number of frequency in the bash DOS
  real(8)                 :: Wbath         !Width of the BATH DOS
  real(8)                 :: Walpha        !exponent of the pseudo-gapped bath.
  real(8)                 :: Wgap          !gap of the gapped bath
  character(len=16)       :: bath_type     !choose the shape of the BATH
  real(8)                 :: eps           !broadening
  real(8)                 :: dmft_error     !convergence error threshold
  character(len=16)       :: field_type    !choose the profile of the electric field
  integer                 :: Nsuccess      !number of convergence success
  character(len=32)       :: g0file,sigfile




  !ELECTRIC FIELD VARIABLES:
  real(8)                 :: Efield        !Electric field strength
  real(8),dimension(3)    :: Evect         !Electric field vectors as input
  real(8)                 :: Ton,Toff      !turn on/off time, t0 also center of the pulse
  integer                 :: Ncycles       !Number of cycles in pulsed light packet
  real(8)                 :: omega0        !parameter for the Oscilatting field and Pulsed light
  real(8)                 :: E1            !Electric field strenght for the AC+DC case (tune to resonate)

  ! !type(vect2D)            :: Ak,Ek         !Electric field vector potential and vector
  ! real(8)                 :: Efield        !Electric field strength
  ! real(8)                 :: Ex,Ey         !Electric field vectors as input
  ! real(8)                 :: t0,t1         !turn on/off time, t0 also center of the pulse
  ! integer                 :: Ncycles       !Number of cycles in pulsed light packet
  ! real(8)                 :: omega0        !parameter for the Oscilatting field and Pulsed light
  ! real(8)                 :: E1            !Electric field strenght for the AC+DC case (tune to resonate)



  ! !CONTAINER FOR THE CONTOUR PARAMETERS 
  ! !=========================================================  
  ! type(kb_contour_params) :: cc_params 
  !MOVED INTO NEQ_CONTOUR



contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine neq_read_input(inputFILE)
    character(len=*)               :: inputFILE
    !GLOBAL
    call parse_input_variable(Norb       , "NORB",inputFILe,default=1,comment="Number of impurity orbitals.")
    call parse_input_variable(Nspin      , "NSPIN",inputFILE,default=1,comment="Number of spin degeneracy (max 2)")
    call parse_input_variable(Ntime      , "NTIME" , inputFILE , default      =100 , comment="Number of Real-Time steps")
    call parse_input_variable(Ntau       , "NTAU" , inputFILE , default       =50 , comment="Number of Imag-Time steps")
    call parse_input_variable(Niw        , "NIW" , inputFILE , default      =4096 , comment="Number of Matsubara Frequencies")
    call parse_input_variable(dt         , "DT" , inputFILE , default         =0.1d0 , comment="Real-time step")
    call parse_input_variable(beta       , "BETA" , inputFILE , default       =10d0 , comment="Inverse temperature")
    call parse_input_variable(Ui         , "Ui" , inputFILE , default         =0d0 , comment="equilibrium local interaction")
    call parse_input_variable(U          , "U" , inputFILE , default          =1d0 , comment="non-equilibrium local interaction")
    call parse_input_variable(xmu        , "XMU" , inputFILE , default        =0.d0 , comment="chemical potential")
    call parse_input_variable(eps        , "EPS" , inputFILE , default        =0.01d0 , comment="broadening")
    call parse_input_variable(nloop      , "NLOOP" , inputFILE , default      =30 , comment="Max number of DMFT loop")
    call parse_input_variable(dmft_error , "DMFT_ERROR" , inputFILE , default =1.d-3 , comment="DMFT convergence threshold")
    call parse_input_variable(Nsuccess   , "NSUCCESS" , inputFILE , default   =1 , comment="number of consecutive successes on convergence")
    !BATH
    call parse_input_variable(Vbath      , "VBATH" , inputFILE , default      =0d0 , comment="coupling to the thermostat")
    call parse_input_variable(Lbath      , "LBATH" , inputFILE , default      =1000 , comment="number of frequencies in the thermostat DOS")
    call parse_input_variable(bath_type  , "BATH_TYPE" , inputFILE , default  ='flat' , comment="thermostat DOS type [flat,gauss,bethe,...]")
    call parse_input_variable(wbath      , "WBATH" , inputFILE , default      =20d0 , comment="Width of the thermostat DOS")
    call parse_input_variable(walpha     , "WALPHA" , inputFILE , default     =1d0 , comment="exponent of the pseudo-gapped thermostat")
    call parse_input_variable(wgap       , "WGAP" , inputFILE , default       =5d0 , comment="gap of the gapped thermostat")
    !EFIELD
    call parse_input_variable(field_type , "FIELD_TYPE" , inputFILE , default ='dc' , comment="profile type of the electric field ")
    call parse_input_variable(Efield     , "EFIELD" , inputFILE , default     =0d0 , comment="electric field strength")
    call parse_input_variable(Evect      , "EVECT" , inputFILE , default      =[1d0,0d0,0d0] , comment="electric field direction (normalized)")
    call parse_input_variable(ton        , "TON" , inputFILE , default        =0d0 , comment="turn on time or center of the pulse")
    call parse_input_variable(toff       , "TOFF" , inputFILE , default       =10000d0 , comment="turn off time")
    call parse_input_variable(ncycles    , "NCYCLES" , inputFILE , default    =1 , comment="number of cycles in pulsed light signal ")
    call parse_input_variable(omega0     , "OMEGA0" , inputFILE , default     =acos(-1d0) , comment="parameter for the Oscilatting field and Pulsed light")
    call parse_input_variable(E1         , "E1" , inputFILE , default         =0d0 , comment="Electric field strenght for the AC+DC case (tune to resonate)")
    call parse_input_variable(g0file     , "G0FILE" , inputFILE , default     ="G0.restart" , comment="File with G0(iw) + header of the form: dens, <H>, <H**2>, <HDC>")
    call parse_input_variable(sigfile    , "SIGFILE" , inputFILE , default    ="S.restart" , comment="File with Sigma(iw) + header")
    call save_input_file(inputFILE)
    call scifor_version()
    call code_version(revision)
  end subroutine neq_read_input



  subroutine neq_save_input(inputFILE)
    character(len=*)               :: inputFILE
    call save_input_file(inputFILE)
  end subroutine neq_save_input



end module NEQ_INPUT_VARS

