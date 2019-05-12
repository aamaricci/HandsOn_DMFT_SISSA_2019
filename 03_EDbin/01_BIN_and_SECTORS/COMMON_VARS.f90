MODULE COMMON_VARS
  implicit none
  integer :: Ns

  type sector_map
     integer,dimension(:),allocatable :: map
  end type sector_map

  interface map_allocate
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate

contains

  subroutine map_allocate_scalar(H,N)
    type(sector_map) :: H
    integer :: N
    allocate(H%map(N))
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:) :: H
    integer,dimension(size(H))    :: N
    integer :: i
    do i=1,size(H)
       allocate(H(i)%map(N(i)))
    enddo
  end subroutine map_allocate_vector

END MODULE COMMON_VARS

