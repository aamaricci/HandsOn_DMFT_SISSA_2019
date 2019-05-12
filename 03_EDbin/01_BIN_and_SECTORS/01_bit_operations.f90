program testBIT
  !
  !< use COMMON_VARIABLES module to load global variables 
  USE COMMON_VARS
  !
  !< fortran-esque to force explicit declaration of all variables.
  implicit none
  !
  !< define variables to be used in THIS code. Variables are shared with
  !any procedure contained here.
  integer             :: i,j,k,m,l
  integer             :: Ntot
  integer             :: iup,idw,int
  integer             :: Nup,Ndw,Sz
  integer             :: Dim
  integer             :: unit
  character(len=6)    :: mode
  character(len=2)    :: fbt,fbt2
  character(len=16)   :: fmt,fmt2
  real(8)             :: sgn1,sgn2
  logical             :: bool



  write(*,*) "This code shows some BIT-function from Fortran" 


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
  write(unit,*)"Exercises with BITS:"
  write(unit,*)"" 
  !
  !
  !< internal writing: F-esque trick: write a number to a variable (as a unit). Need to
  !specify the format. I2.2 == Integer, 2spaces. 2 digits. FBT=string of length 2
  !
  !< define a suitable format to pretty print the Integers and their Bit decomposition.
  !I9 = Integer 9spaces. Bn.m is a bit format for integers, n spaces, m digits.
  !remark: Bn.m print write the number right-to-left (as default for bits)
  !
  !Format for 2*Ns digits integers
  write(fbt,'(I2.2)')2*Ns
  fmt="(I9,A,B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
  !
  !Format for Ns digits integers
  write(fbt2,'(I2.2)')Ns
  fmt2="(I9,A,B"//adjustl(trim(fbt2))//"."//adjustl(trim(fbt2))//")"
  !
  !
  !< Let's pick a number: 13 ==> 001101 
  i=13                          !001101
  !
  !< write it in Integer + Bit format
  write(unit,fmt)i,' => ',i
  call sleep(1)
  !
  !
  !
  write(unit,*)""
  write(unit,*)"BTEST: loop over 0:2*Ns-1 bits and check if the corresponding bit is 0 or 1"
  write(unit,fmt)i,' => ',i
  write(unit,"(A)",advance='no')"             " !conform to FMT conform.
  do k=2*Ns-1,0,-1              !we actually loop backward to conform to FMT format
     write(unit,"(L)",advance="no")btest(i,k)
  enddo
  write(unit,*)""
  call sleep(1)
  !
  !
  !< Show IBSET functioning. Set a bit to 1
  write(unit,*)""
  write(unit,*)"IBSET: loop over 0:2*Ns-1 bits and sets each bit to 1, print the resulting integer"
  do k=0,2*Ns-1
     j = ibset(i,k)
     write(unit,fmt)j,' => ',j
  enddo
  write(unit,*)""
  call sleep(1)
  !
  !
  !< Show IBCLR functioning. Set a bit to 0
  write(unit,*)"IBCLR:loop over 0:2*Ns-1 bits and sets each bit to 0, print the resulting integer"
  do k=0,2*Ns-1
     j = ibclr(i,k)
     write(unit,fmt)j,' => ',j
  enddo
  write(unit,*)""
  call sleep(1)
  !
  !
  write(unit,*)"IEOR: returns the bitwise Boolean exclusive-OR of I and K"
  write(unit,fmt)i,' => ',i
  k=2**1+2**3
  write(unit,fmt)k,' => ',k
  write(unit,*)""
  j = ieor(i,k)
  write(unit,fmt)j,' => ',j
  write(unit,*)""
  call sleep(1)
  !
  !
  !
  write(unit,*)"IBITS:extracts a field of length LEN from I, starting from bit position POS and extending left for LEN bits"
  write(unit,fmt)i,' => ',i
  write(unit,fmt2)ibits(i,0,3),' => ',ibits(i,0,3)
  write(unit,fmt2)ibits(i,3,3),' => ',ibits(i,3,3)
  write(unit,*)"5 + 1 * 2^3 = 13"
  write(unit,*)""


  write(unit,*)"DONE..."
end program testBIT
