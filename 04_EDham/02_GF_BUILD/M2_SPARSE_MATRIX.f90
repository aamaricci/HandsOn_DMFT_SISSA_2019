!< THIS MODULE CONTAINS A SPECIAL, FLEXIBLE, DATA STRUCTURE TO CONTAINS AND
! HANDLE A SPARSE MATRIX. ESSENTALLY:
! A SPARSE MATRIX M IS STORED AS A NUMBER OF NROW OF DYNAMIC ROWS, CONTAININNG VALUES & COLUMNS
! OF NON-ZERO ELEMENTS ONLY. THE ROWS ARE DYNAMICS ARRAYS WHICH RESIZE ON REQUEST. ACCESS IS O(1)
MODULE SPARSE_MATRIX
  USE COMMON_VARS, only: free_unit,str
  implicit none
  private
  !
  !< DYNAMICS SPARSE ROW (COMPRESSED STORE ROW)
  type sparse_row_csr
     integer                                   :: size !actual size
     real(8),dimension(:),allocatable          :: vals !values
     integer,dimension(:),allocatable          :: cols !columns
  end type sparse_row_csr
  !
  !< SPARSE MATRIX DATA STRUCTURE
  type sparse_matrix_csr
     type(sparse_row_csr),dimension(:),pointer :: row !an array of rows (number of rows is known a priori
     integer                                   :: Nrow !number of rows
     integer                                   :: Ncol !number of columns
     logical                                   :: status=.false. !bool for allocation 
  end type sparse_matrix_csr
  !
  !
  !< DEFINE SOME OPERATIONS TO BE PERFORMED ON THE SPARSE MATRIX OBJECTS 
  !INIT SPARSE MATRICES 
  interface sp_init_matrix
     module procedure :: sp_init_matrix_csr
  end interface sp_init_matrix
  !
  !DELETE SPARSE MATRIX 
  interface sp_delete_matrix
     module procedure :: sp_delete_matrix_csr
  end interface sp_delete_matrix
  !
  !INSERT ELEMENTS
  interface sp_insert_element
     module procedure :: sp_insert_element_csr
  end interface sp_insert_element
  !
  !LOAD STANDARD MATRIX INTO SPARSE MATRICES
  interface sp_load_matrix
     module procedure :: sp_load_matrix_csr
  end interface sp_load_matrix
  !
  !DUMP SPARSE MATRIX INTO STANDARD MATRIX
  interface sp_dump_matrix
     module procedure :: sp_dump_matrix_csr
  end interface sp_dump_matrix
  !
  !SPY PRINT SPARSE MATRIX
  interface sp_spy_matrix
     module procedure :: sp_spy_matrix_csr
     module procedure :: sp_spy_matrix_dense
  end interface sp_spy_matrix
  !
  !
  !< DECLARE PUBLIC THE NEEDED VARIABLES, OBJECTS AND PROCEDURES
  !CSR Sparse Matrix
  public :: sparse_matrix_csr
  public :: sp_init_matrix      !init the sparse matrix   !checked
  public :: sp_delete_matrix    !delete the sparse matrix !checked
  public :: sp_insert_element   !insert an element        !checked
  public :: sp_load_matrix      !create sparse from array !checked
  public :: sp_dump_matrix      !dump sparse into array   !checked
  public :: sp_spy_matrix       !print sparse matrix nicely !checked
  !
  !
  !
  !< interface to procedures required to dynamically re-allocate and extend each row
  ! i.e. resize the vectors values,columns in each row if another element is to be added.
  interface add_to
     module procedure :: add_to_I
     module procedure :: add_to_D
     module procedure :: add_to_Z
  end interface add_to
  !
  !
  !
  !< SPARSE MATRIX OBJECT TO BE USED IN THE TEST CODE
  type(sparse_matrix_csr),public  :: spH0 !DIAGONAL PART
  type(sparse_matrix_csr),public  :: spH0up !UP-PART
  type(sparse_matrix_csr),public  :: spH0dw !DW-PART




contains       


  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  subroutine sp_init_matrix_csr(sparse,N,N1)
    type(sparse_matrix_csr),intent(inout) :: sparse
    integer                               :: N
    integer,optional                      :: N1
    integer                               :: i
    !
    !put here a delete statement to avoid problems
    if(sparse%status)stop "sp_init_matrix: alreay allocate can not init"
    !
    sparse%Nrow=N
    sparse%Ncol=N 
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%row(N))
    do i=1,N
       sparse%row(i)%size=0
       allocate(sparse%row(i)%vals(0)) !empty array
       allocate(sparse%row(i)%cols(0)) !empty array
    end do
    !
    sparse%status=.true.
    !
  end subroutine sp_init_matrix_csr








  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_matrix_csr(sparse)    
    type(sparse_matrix_csr),intent(inout) :: sparse
    integer                           :: i
    type(sparse_row_csr),pointer          :: row
    !
    if(.not.sparse%status)stop "Warning SPARSE/sp_delete_matrix: sparse not allocated already."
    !
    do i=1,sparse%Nrow
       deallocate(sparse%row(i)%vals)
       deallocate(sparse%row(i)%cols)
       sparse%row(i)%Size  = 0
    enddo
    deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
  end subroutine sp_delete_matrix_csr















  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element_csr(sparse,value,i,j)
    type(sparse_matrix_csr),intent(inout) :: sparse
    real(8),intent(in)                    :: value
    integer,intent(in)                    :: i,j
    !
    type(sparse_row_csr),pointer          :: row
    integer                               :: column,pos
    logical                               :: iadd
    !
    column = j
    !
    row => sparse%row(i)
    !
    iadd = .false.                          !check if column already exist
    if(any(row%cols == column))then         !
       pos = binary_search(row%cols,column) !find the position  column in %cols        
       iadd=.true.                          !set Iadd to true
    endif
    !
    if(iadd)then                            !this column exists so just sum it up       
       row%vals(pos)=row%vals(pos) + value  !add up value to the current one in %vals
    else                                    !this column is new. increase counter and store it 
       ! row%vals = [row%vals,value]
       ! row%cols = [row%cols,column]
       call add_to(row%vals,value)
       call add_to(row%cols,column)
       row%Size = row%Size + 1
    endif
    !
    if(row%Size > sparse%Ncol)stop "sp_insert_element_csr ERROR: row%Size > sparse%Ncol"
    !
  end subroutine sp_insert_element_csr








  !+------------------------------------------------------------------+
  !PURPOSE: load a regular matrix (2dim array) into a sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_load_matrix_csr(matrix,sparse)
    real(8),dimension(:,:),intent(in) :: matrix
    type(sparse_matrix_csr),intent(inout) :: sparse    
    integer                           :: i,j,Ndim1,Ndim2
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)   
    !
    if(sparse%status)call sp_delete_matrix_csr(sparse)
    call sp_init_matrix_csr(sparse,Ndim1,Ndim2)
    !
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=0.d0)call sp_insert_element_csr(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine sp_load_matrix_csr







  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix_csr(sparse,matrix)
    type(sparse_matrix_csr),intent(in)   :: sparse
    real(8),dimension(:,:),intent(inout) :: matrix
    integer                              :: i,j,Ndim1,Ndim2
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_matrix: dimensions error"
    !
    matrix=0.d0
    do i=1,Ndim1
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = matrix(i,sparse%row(i)%cols(j)) + sparse%row(i)%vals(j)
       enddo
    enddo
  end subroutine sp_dump_matrix_csr






  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+
  subroutine sp_print_matrix_csr(sparse,unit,fmt)
    type(sparse_matrix_csr)          :: sparse
    integer,optional             :: unit
    integer                      :: unit_
    integer                      :: i,j,Ns
    character(len=*),optional    :: fmt
    character(len=64)            :: fmt_
    type(sparse_row_csr),pointer     :: row
    integer                      :: count=0
    unit_=6;if(present(unit))unit_=unit
    fmt_='F15.9';if(present(fmt))fmt_=fmt
    write(*,*)"Print sparse matrix (compact mode) ->",unit_
    do i=1,sparse%Nrow
       row => sparse%row(i)
       do j=1,row%Size
          write(unit_,"("//trim(fmt_)//",A1,I0,3X)",advance='no')row%vals(j),',',row%cols(j)
       end do
       write(unit_,*)
    enddo
    write(unit_,*)
  end subroutine sp_print_matrix_csr









  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+  
  subroutine sp_spy_matrix_csr(sparse,header)
    type(sparse_matrix_csr)          :: sparse
    character ( len = * )           :: header
    integer                         :: N1,N2
    character ( len = 255 )         :: command_filename
    integer                         :: command_unit
    character ( len = 255 )         :: data_filename
    integer                         :: data_unit
    integer                         :: i, j
    character ( len = 6 )           :: n1_s,n2_s,n1_i,n2_i
    integer                         :: nz_num
    character ( len = 255 )         :: png_filename
    !
    !  Create data file.
    !
    !
    N1 = sparse%Nrow
    N2 = sparse%Ncol
    data_filename = trim ( header ) // '_data.dat'
    open (unit=free_unit(data_unit), file = data_filename, status = 'replace' )
    nz_num = 0
    do i=1,N1
       do j=1,sparse%row(i)%size
          write(data_unit,'(2x,i6,2x,i6)') sparse%row(i)%cols(j),i
          nz_num = nz_num + 1
       enddo
    enddo
    close(data_unit)
    !
    write(*,"(A,A)")trim(header)," analysis:"
    write(*,"(A,2I12)") " number of non-zero elements, dimensions:",nz_num,N1*N2
    write(*,"(A,F7.2,A1)")" sparse-ness                            :",dble(nz_num)/N1/N2*100,"%"
    write(*,*)""
    !
    !  Create command file.
    !
    command_filename = "plot_"//str(header)//'_commands.gp'
    open(unit = free_unit(command_unit), file = command_filename, status = 'replace' )
    write(command_unit,'(a)') '#unset key'
    write(command_unit,'(a)') 'set terminal postscript eps enhanced color font "Times-Roman,16"'
    write(command_unit,'(a)') 'set output "|ps2pdf -sEPSCrop - '//str(header)//".pdf"//'"'
    write(command_unit,'(a)') 'set size ratio -1'
    write(command_unit,'(a)') 'set xlabel "<--- J --->"'
    write(command_unit,'(a)') 'set ylabel "<--- I --->"'
    write(command_unit,'(a,i6,a)')'set title "',nz_num,' nonzeros for '//str(header)//'"'
    write(command_unit,'(a)') 'set timestamp'
    write(command_unit,'(a)' )'plot [x=1:'//str(N1)//'] [y='//str(N2)//':1] "'//&
         str(data_filename)//'" w p pt 5 ps 0.4 lc rgb "red"'
    close ( unit = command_unit )
    return
  end subroutine sp_spy_matrix_csr


  subroutine sp_spy_matrix_dense(h,header)
    real(8),dimension(:,:)          :: h
    character ( len = * )           :: header
    integer                         :: N1,N2
    character ( len = 255 )         :: command_filename
    integer                         :: command_unit
    character ( len = 255 )         :: data_filename
    integer                         :: data_unit
    integer                         :: i, j
    character ( len = 6 )           :: n1_s,n2_s,n1_i,n2_i
    integer                         :: nz_num
    character ( len = 255 )         :: png_filename
    !
    !  Create data file.
    !
    !
    N1 = size(h,1)
    N2 = size(h,2)
    data_filename = trim ( header ) // '_data.dat'
    open (unit=free_unit(data_unit), file = data_filename, status = 'replace' )
    nz_num = 0
    do i=1,N1
       do j=1,N2
          if(h(i,j)/=0d0)then
             write(data_unit,'(2x,i6,2x,i6)') j,i
             nz_num = nz_num + 1
          endif
       enddo
    enddo
    close(data_unit)
    !
    write(*,"(A,A)")trim(header)," analysis:"
    write(*,"(A,2I12)")" number of non-zero elements, dimensions:",nz_num,N1*N2
    write(*,"(A,F7.2,A1)")" sparse-ness                            :",dble(nz_num)/N1/N2*100,"%"
    write(*,*)""
    !
    !  Create command file.
    !
    command_filename = "plot_"//str(header)//'_commands.gp'
    open(unit = free_unit(command_unit), file = command_filename, status = 'replace' )
    write(command_unit,'(a)') '#unset key'
    write(command_unit,'(a)') 'set terminal postscript eps enhanced color font "Times-Roman,16"'
    write(command_unit,'(a)') 'set output "|ps2pdf -sEPSCrop - '//str(header)//".pdf"//'"'
    write(command_unit,'(a)') 'set size ratio -1'
    write(command_unit,'(a)') 'set xlabel "<--- J --->"'
    write(command_unit,'(a)') 'set ylabel "<--- I --->"'
    write(command_unit,'(a,i6,a)')'set title "',nz_num,' nonzeros for '//str(header)//'"'
    write(command_unit,'(a)') 'set timestamp'
    write(command_unit,'(a)' )'plot [x=1:'//str(N1)//'] [y='//str(N2)//':1] "'//&
         str(data_filename)//'" w p pt 5 ps 0.4 lc rgb "red"'
    close ( unit = command_unit )
    return
  end subroutine sp_spy_matrix_dense









  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       !stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search

  subroutine add_to_I(vec,val)
    integer,dimension(:),allocatable,intent(inout) :: vec
    integer,intent(in)                             :: val  
    integer,dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_I

  subroutine add_to_D(vec,val)
    real(8),dimension(:),allocatable,intent(inout) :: vec
    real(8),intent(in)                             :: val  
    real(8),dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_D

  subroutine add_to_Z(vec,val)
    complex(8),dimension(:),allocatable,intent(inout) :: vec
    complex(8),intent(in)                             :: val  
    complex(8),dimension(:),allocatable               :: tmp
    integer                                           :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_Z

end module SPARSE_MATRIX






