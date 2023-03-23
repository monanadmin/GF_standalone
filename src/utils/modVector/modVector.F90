module ModVector

  implicit none

  private

  character(len=*), parameter   :: p_source_name = 'modVector.F90'
  
  
  ! Vector inside data type
  type :: data_t
     integer :: index_value
  end type data_t
  ! Vector data type
  
  type :: vector_t
    private
    type(data_t), pointer, dimension(:) :: vector
    integer :: num_elements
    
  end type vector_t

  type(vector_t) :: instance

  ! public :: init_instance
  public :: init
  public :: free_memory
  public :: get_size
  public :: get_num_elements
  public :: insert
  public :: vector_put
  public :: get
  public :: get_index_value
  public :: remove
  public :: print_all
  public :: insert_range


contains

  ! Initialize vector
  subroutine init(vec_max_size)
    implicit none 
    integer, intent(in) :: vec_max_size

    instance%num_elements = 0
    allocate(instance%vector(vec_max_size))
  end subroutine init


  ! Insert a range of values in vector
  subroutine insert_range(start_val, end_val)
    implicit none
     
    integer, intent(in)           :: start_val
    integer, intent(in)           :: end_val
    integer                       :: i, vector_ind
    character(len=*), parameter   :: p_procedure_name = 'insert_range' 
     
    if (end_val < start_val) then
      print*, p_procedure_name, "Erro: end_val", end_val, " menor que start_val", start_val 
      stop
    endif
            
    vector_ind=0      
    do i = start_val, end_val
      vector_ind=vector_ind+1
      instance%vector(vector_ind)%index_value=i
    enddo
    instance%num_elements = vector_ind
  
  end subroutine insert_range


  ! Free the entire list and all data, beginning at SELF
  subroutine free_memory()
    implicit none
    if(associated(instance%vector)) deallocate(instance%vector)
    instance%num_elements = 0
  end subroutine free_memory


  ! get vector max size (by init)
  function get_size() result(size_vec)
    integer :: size_vec
    if(associated(instance%vector)) then
      size_vec = size(instance%vector)
    else
      size_vec = 0
    endif
  end function get_size


  ! get vector actual num of elements
  function get_num_elements() result(num_elements)
    integer :: num_elements
    num_elements = instance%num_elements
  end function get_num_elements


  ! Insert a value at end of the vector
  ! TODO - change data to index_value 
  subroutine insert(data)
    type(data_t), intent(in) :: data

    instance%num_elements = instance%num_elements + 1
    instance%vector(instance%num_elements) = data
  end subroutine insert


  ! Store the encoded data the index_data position
  ! TODO - change data to index_value 
  subroutine vector_put(data, data_index)
    type(data_t), intent(in) :: data
    integer, intent(in) :: data_index

    instance%vector(data_index) = data
  end subroutine vector_put


  ! Return the DATA stored in data_index
  function get(data_index) result(data)
    integer, intent(in) :: data_index
    type(data_t) :: data
    data = instance%vector(data_index)
  end function get
  
  function get_index_value(data_index) result(data)
    integer, intent(in) :: data_index
    integer :: data
    data = instance%vector(data_index)%index_value
  end function get_index_value

  
  ! print all elements of vector
  subroutine print_all() 
    integer :: data_index
    integer, parameter :: views = 5
    write(*, '(A8)', advance='NO') 'vector = ('
    do data_index = 1, min(instance%num_elements, views)
      write(*,'(i8, "," )',advance='NO')  instance%vector(data_index)%index_value
    end do
    if (instance%num_elements > views) then
      ! print last elements 
      write(*,'(A5)',advance='NO') ' ... '
      do data_index = max(instance%num_elements - views, views +1 ), instance%num_elements
        write(*,'(i8, "," )',advance='NO')  instance%vector(data_index)%index_value
      end do
    endif
    write(*, '(A2)', advance='YES') ' )'

  end subroutine print_all


  ! remove element containing data from parameter
  function remove(index_value_param) result(is_removed)
    integer, intent(in) :: index_value_param
    integer :: index
    logical :: is_removed

    is_removed = .false.
    if (instance%num_elements == 0) then
      print *, '***** trying to remove from an empty vector. Ignoring'
      return
    endif

    do index = 1, instance%num_elements
      if (instance%vector(index)%index_value == index_value_param) then
        instance%vector(index:instance%num_elements-1) = instance%vector(index+1:instance%num_elements)  ! Bidu
        instance%num_elements = instance%num_elements -1
        is_removed = .true.
      endif
    enddo

  end function remove

end module ModVector