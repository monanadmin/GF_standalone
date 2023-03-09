module ModVector

  use moddata
  implicit none

  private

  public :: vector_t
  public :: init
  public :: free_memory
  public :: get_size
  public :: get_num_elements
  public :: insert
  public :: vector_put
  public :: get
  public :: remove
  public :: print_all
  

  ! Vector data type
  type :: vector_t
    private
    type(data_t), allocatable, dimension(:) :: vector
    integer :: num_elements
    
  end type vector_t

contains


  ! Initialize vector
  subroutine init(self, vec_max_size)
    implicit none 
    type(vector_t), intent(inout) :: self
    integer, intent(in) :: vec_max_size

    self%num_elements = 0
    allocate(self%vector(vec_max_size))
  end subroutine init


  ! Free the entire list and all data, beginning at SELF
  subroutine free_memory(self)
    type(vector_t), intent(inout) :: self

    if(allocated(self%vector)) deallocate(self%vector)
    self%num_elements = 0
  end subroutine free_memory


  function get_size(self) result(size_vec)
    type(vector_t), intent(in) :: self
    integer :: size_vec
    if(allocated(self%vector)) then
      size_vec = size(self%vector)
    else
      size_vec = 0
    endif
  end function get_size


  function get_num_elements(self) result(num_elements)
    type(vector_t), intent(in) :: self
    integer :: num_elements
    num_elements = self%num_elements
  end function get_num_elements

  ! ToDo ....

  ! Insert a value at end of the vector
  subroutine insert(self, data)
    type(vector_t), intent(inout) :: self
    type(data_t), intent(in) :: data

    self%num_elements = self%num_elements + 1
    self%vector(self%num_elements) = data
  end subroutine insert


  ! Store the encoded data the index_data position
  subroutine vector_put(self, data, data_index)
    type(vector_t), intent(inout) :: self
    type(data_t), intent(in) :: data
    integer, intent(in) :: data_index

    self%vector(data_index) = data
  end subroutine vector_put


  ! Return the DATA stored in data_index
  function get(self, data_index) result(data)
    type(vector_t), intent(inout) :: self
    integer, intent(in) :: data_index
    type(data_t) :: data
    data = self%vector(data_index)
  end function get

  subroutine print_all(self) 
    type(vector_t), intent(inout) :: self
    integer :: data_index
    integer, parameter :: views = 5
    write(*, '(A8)', advance='NO') 'vector = ('
    do data_index = 1, min(self%num_elements, views)
      write(*,'(i8, "," )',advance='NO')  self%vector(data_index)%x
    end do
    if (self%num_elements > views) then
      ! print last elements 
      write(*,'(A5)',advance='NO') ' ... '
      do data_index = max(self%num_elements - views, views +1 ), self%num_elements
        write(*,'(i8, "," )',advance='NO')  self%vector(data_index)%x
      end do
    endif
    write(*, '(A2)', advance='YES') ' )'

  end subroutine print_all


  function remove(self, data_to_remove) result(is_removed)
    type(vector_t), intent(inout) :: self
    type(data_t), intent(in) :: data_to_remove
    integer :: index
    logical :: is_removed

    is_removed = .false.
    if (self%num_elements == 0) then
      print *, '***** trying to remove from an empty vector. Ignoring'
      return
    endif

    do index = 1, self%num_elements
      if (self%vector(index)%x == data_to_remove%x) then
        self%vector(index:self%num_elements-1) = self%vector(index+1:self%num_elements)  ! Bidu
        self%num_elements = self%num_elements -1
        is_removed = .true.
      endif
    enddo

  end function remove




end module ModVector