program testEDNupNdw
  !< use COMMON_VARIABLES module to load global variables 
  USE COMMON_VARS
  !
  !< use AUXILIARY FUNCTIONS module to load some auxiliary procedures: binary_decomposition, sector_buildup, etc...
  USE AUX_FUNX
  !
  !< fortran-esque to force explicit declaration of all variables.
  implicit none
  !
  !< define variables to be used in THIS code. Variables are shared with
  !any procedure contained here.
  integer             :: i,j,k,m,l
  integer             :: Ntot,isector
  integer             :: iup,idw,isz,int
  integer             :: Nup,Ndw,Sz
  integer             :: Dim
  integer,allocatable :: GetSector(:,:),getNup(:),getNdw(:),getSz(:),getNt(:),GetDim(:)
  integer,allocatable :: hmap(:),UPmap(:),DWmap(:)
  integer,allocatable :: ivec(:)
  integer             :: unit
  character(len=6)    :: mode
  character(len=2)    :: fbt,fbt2
  character(len=16)   :: fmt,fmt2
  real(8)             :: sgn1,sgn2
  logical             :: bool
  !
  !< CUSTOM TYPE: a user defined variable/structure which contains
  !the map between the states of a given sector S and those of the Fock space F
  type(sector_map)    :: H,Hup,Hdw
  !
  !  
  !< define some useful parameters:
  !
  !< Number of spins: keep this low for better visualization
  Ns    = 3
  !< Total number of levels for fermions. 1Fermion=2*spin
  Ntot  = 2*Ns
  !< set default output unit. 6==std.output in fortran-esque
  unit=6
  write(unit,*)"Ns     =",Ns
  write(unit,*)"Ntot   =",Ntot
  write(unit,*)"2**Ns  =",2**Ns
  write(unit,*)"2**Ntot=",2**Ntot
  write(unit,*)""
  write(unit,*)""
  !
  !
  !< Build recursively all sectors: This is just a test run.
  !Building a sector de-facto means to allocate and construct the UP-DW maps of the two sub-sectors
  !to the Fock space. Each sub-map, i.e. Hup%map, contains the ordered (!) list of the Fock spaces which correspond
  !to each state of the sector: e.g.
  !S --> F
  !1 --> 11
  !2 --> 13
  !3 --> 14
  !4 --> 19
  !5 --> 21
  !6 --> 22 
  !7 --> 35 
  !8 --> 37 
  !9 --> 38  
  isector=0
  do Nup=0,Ns
     do Ndw=0,Ns
        isector=isector+1
        write(unit,"(A,I6,A3,I2,A1,I2,A)")"--- sector",isector,"  (",Nup,",",Ndw,")---- "
        call build_sector_normal(Nup,Ndw,Hup,Hdw)
        deallocate(Hup%map,Hdw%map)
     enddo
  enddo
  write(unit,*)""
  write(unit,*)""
  call sleep(2)
  !
  !
  !
  !< Test that the given sector (i.e. the Maps S-->F) actually reproduce the states of the Fock space
  Nup=2
  Ndw=1
  write(unit,*)""
  write(unit,"(A,I2,A1,I2,A,I10)")"Test --- sector (",Nup,",",Ndw,")---- ",dim_sector_normal(Nup)*dim_sector_normal(Ndw)
  call build_sector_normal(Nup,Ndw,Hup,Hdw)
  !
  write(unit,*)"I = I_UP + I_DW*2*Ns"
  write(unit,*)"Map: UP",Hup%map
  write(unit,*)"                              X        "
  write(unit,*)"Map: DW",Hdw%map, " * 2^3"
  write(unit,*)""
  do j = 1,dim_sector_normal(Ndw)
     do i = 1,dim_sector_normal(Nup)
        k = Hup%map(i) + 2**Ns*Hdw%map(j)
        write(unit,'(i6)',advance='no')k
     enddo
  enddo
  write(unit,*)""
  deallocate(Hup%map,Hdw%map)
  write(unit,*)""
  write(unit,*)""
  write(unit,*)""
  call sleep(2)
  !
  !
  !
  !< Apply a combination of operators C^dagger_i C_j to hop an electron from i  to J 
  Nup=2
  Ndw=1
  write(unit,"(A,I2,A1,I2,A)")"Cdg_2up * C_3up ---  --- in sector (",Nup,",",Ndw,")---- "
  write(unit,'(A)')"Apply C_up^+_2 Cup_3 on |I_up>|I_dw>"
  !
  !< Allocate a vector containing the binary decomposition of I (the starting state in F)  
  allocate(Ivec(Ns))
  !
  !< build the sector (i.e. construct Hup,Hdw maps)
  call build_sector_normal(Nup,Ndw,Hup,Hdw)  !
  !
  !|m> = |m_up>|m_dw> (m_dw is any here) ==> |I> = |I_up>|I_dw> (I_dw is any here)
  !
  ! Cdg_2up *C_3up|m> = Cdg_2up *C_3up|I> = Cdg_2up *C_3up|I_up>|I_dw>
  ! = |J_up>|J_dw> == |J_up>|I_dw> (J_dw=I_dw no spin dw has been changed)
  ! <== |l_up>|l_dw> = |l>
  !
  m = 2
  I = Hup%map(m)
  write(unit,'(A,I3)')"m_up=",m
  write(unit,'(A,I3,A,I3,A)')"I_up=",I," = UPmap(",m,")"
  !
  !< Binary decomposition of I_up \in Fock_up: Ivec = 101
  ivec = bdecomp(I,Ns) 
  call print_conf(I,Ns)
  !
  !< Before applying the operator check that there is no electron at position 2 and there is one at position 3
  if(ivec(3)==1.AND.ivec(2)==0)then
     !< apply from the right (a first)
     call a(3,i,k,sgn1);call print_conf(k,Ns)
     call adg(2,k,j,sgn2);call print_conf(j,Ns)
  end if
  !< The result of the application of the operators is a fermionic sign:
  write(unit,*)sgn2
  !< retrieve the sector state l_up so that: which corresponds to the end state J
  !use binary search to loop up for the J state in the ordered map Hup (scale logN but saves
  !a lot of memory.
  l = binary_search(Hup%map,J)
  write(unit,'(A,I3)')"J_up=",J
  write(unit,'(A,I3,A,I3,A)')"l_up=",l," = inv_UPmap(",J,")"
  write(unit,*)""
  deallocate(Hup%map,Hdw%map)
  deallocate(Ivec)



end program testEDNupNdw
