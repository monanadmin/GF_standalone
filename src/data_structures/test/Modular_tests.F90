! A simple generic linked list test program
program Modular_tests

  use Moddata

  implicit none

  logical all_tests_passed 
  
  integer, parameter :: p_num_of_insertions = 1000000
  integer, parameter :: p_num_of_deletions  = 1000

  real, parameter :: p_insert_max_time_allowed = 0.2
  real, parameter :: p_visit_max_time_allowed  = 0.1
  real, parameter :: p_remove_max_time_allowed = 10.0

  print*, ">>>>> Running Modular Tests"
  
  all_tests_passed = .true.
  all_tests_passed = all_tests_passed .and. test_insert_performance()
  all_tests_passed = all_tests_passed .and. test_list_remove_performance()
  all_tests_passed = all_tests_passed .and. test_vector_insert_performance()
  all_tests_passed = all_tests_passed .and. test_vector_remove_performance()

  if (all_tests_passed) then
    print*, ">>>>> All tests OK !"
    call exit(0)
  else
    print*, ">>>>> Some Tests failed!"
    call exit(-1)
  endif


  contains 


    logical function test_insert_performance() result(test_result)
      use Modlist
      implicit none

      type(list_t), pointer :: ll => null()
      type(list_t), pointer :: list_pointer => null()
      type(data_t), allocatable :: dat_a

      integer :: an_integer
      real time_initial, time_final 


      print*, ">>>>> Running Test List Insert and Visit Performance"
      test_result = .true.

      ! List Insert
      !
      call cpu_time(time_initial)
      allocate(dat_a)
      dat_a%x = 1
      call init(ll, dat_a)
      deallocate(dat_a)
      do an_integer = 2, p_num_of_insertions
        allocate(dat_a)
        dat_a%x = an_integer
        call insert(ll, DATA=dat_a)
        deallocate(dat_a)
      enddo
      call cpu_time(time_final)
      
      ! check performance
      print *, 'Insert Max time = ', p_insert_max_time_allowed, '; Run time = ', time_final-time_initial
      if (time_final-time_initial > p_insert_max_time_allowed) then
        print *, 'Insert Max time exceeded!!!'
        test_result = .false.
        call free_memory(ll)
        return
      endif

      ! check HEAD
      dat_a = get(ll)
      print*, 'Checking head element ...', dat_a%x
      if(dat_a%x /= p_num_of_insertions) then
        print *, 'Head element should be: ', p_num_of_insertions, ' but was', dat_a%x
        test_result = .false.
        call free_memory(ll)
        return
      endif

      ! List Visit test
      call cpu_time(time_initial)
      list_pointer => ll
      do 
        ! do a quick visit
        dat_a = get(list_pointer)
        list_pointer => next(list_pointer)
        if (.not. associated(list_pointer)) exit
      enddo
      call cpu_time(time_final)


      print*, 'Checking last element ...', dat_a%x
      if (dat_a%x /= 1) then
        print *, 'Last element should be 1 but was', dat_a%x
        test_result = .false.
        call free_memory(ll)
        return
      endif


      print *, 'Visit Max time = ', p_visit_max_time_allowed, '; Run time = ', time_final-time_initial
      if (time_final-time_initial > p_visit_max_time_allowed) then
        print *, 'Visit Max time exceeded!!!'
        test_result = .false.
        call free_memory(ll)
        return
      endif

      call free_memory(ll)
      return
    end function test_insert_performance


    logical function test_list_remove_performance() result(test_result)
      use Modlist
      implicit none

      type(list_t), pointer :: ll => null()
      type(list_t), pointer :: list_pointer => null()
      type(data_t) :: dat_a, data_to_remove


      integer :: an_integer
      real :: time_initial, time_final 

      integer :: list_size 
      logical :: is_removed
  
      print*, ">>>>> Running Test List random remove Performance"
      test_result = .true.

      list_size = p_num_of_insertions

      ! Insert values
      dat_a%x = 1
      call init(ll, dat_a)
      do an_integer = 2, p_num_of_insertions
        dat_a%x = an_integer
        call insert(ll, DATA=dat_a)
      enddo

      print *, "removing ", p_num_of_deletions, " elements of list of size ", p_num_of_insertions
      call cpu_time(time_initial)
      
      do an_integer = 1, p_num_of_deletions
        list_pointer => next(ll)
        ! removes using a step
        data_to_remove%x = an_integer * (p_num_of_insertions/p_num_of_deletions) - an_integer
        ! print *, 'removing ', data_to_remove%x
        is_removed = remove(list_pointer, data_to_remove)
        if(is_removed) list_size = list_size -1
      enddo
      call cpu_time(time_final)


      print *, 'Remove Max time = ', p_remove_max_time_allowed, '; Run time = ', time_final-time_initial
      if (time_final-time_initial > p_remove_max_time_allowed) then
        print *, 'Remove Max time exceeded!!!'
        test_result = .false.
        call free_memory(ll)
        return
      endif

      call free_memory(ll)
      return

    end function test_list_remove_performance


    logical function test_vector_insert_performance() result(test_result)
      use ModVector
      implicit none

      type(vector_t) :: vec 
      type(data_t) :: dat_a

      integer :: an_integer, data_index
      real time_initial, time_final 

      print*, ">>>>> Running Test Vector Insert and Visit Performance"
      test_result = .true.

      call cpu_time(time_initial)
      call init(vec, p_num_of_insertions)
      do an_integer = 1, p_num_of_insertions
        dat_a%x = an_integer
        call insert(vec, DATA=dat_a)
      enddo
      call cpu_time(time_final)
      call print_all(vec)
      
      ! check performance
      print *, 'Insert Max time = ', p_insert_max_time_allowed, '; Run time = ', time_final-time_initial
      if (time_final-time_initial > p_insert_max_time_allowed) then
        print *, 'Insert Max time exceeded!!!'
        test_result = .false.
        call free_memory(vec)
        return
      endif

      data_index=1
      dat_a = get(vec,data_index)
      print*, 'Checking Head element ...', dat_a%x
      if (dat_a%x /= 1) then
        print *, 'Head element should be 1 but was', dat_a%x
        test_result = .false.
        call free_memory(vec)
        return
      endif

      data_index=p_num_of_insertions
      dat_a = get(vec,data_index)
      print*, 'Checking last element ...', dat_a%x
      if(dat_a%x /= p_num_of_insertions) then
        print *, 'Last element should be: ', p_num_of_insertions, ' but was', dat_a%x
        test_result = .false.
        call free_memory(vec)
        return
      endif

      ! do a quick visit to every element
      call cpu_time(time_initial)
      do an_integer = 1, get_num_elements(vec)
        dat_a = get(vec,data_index)
      enddo
      call cpu_time(time_final)


      print *, 'Visit Max time = ', p_visit_max_time_allowed, '; Run time = ', time_final-time_initial
      if (time_final-time_initial > p_visit_max_time_allowed) then
        print *, 'Visit Max time exceeded!!!'
        test_result = .false.
        call free_memory(vec)
        return
      endif

      call free_memory(vec)
      return

    end function test_vector_insert_performance

    logical function test_vector_remove_performance() result(test_result)
      use ModVector
      implicit none
      type(vector_t) :: vec
      type(data_t) :: dat_a

      integer :: an_integer
      real :: time_initial, time_final 

      type(data_t) :: data_to_remove
      logical :: dummy

      print*, ">>>>> Running Test vector random remove Performance"
      test_result = .true.
      ! Insert values
      dat_a%x = 1
      call init(vec, p_num_of_insertions)
      do an_integer = 1, p_num_of_insertions
        call insert(vec, DATA=dat_a)
        dat_a%x = an_integer+1
      enddo
      call print_all(vec)
      print *, "removing ", p_num_of_deletions, " elements of vector of size ", p_num_of_insertions
      call cpu_time(time_initial)
      do an_integer = 1, p_num_of_deletions
        ! removes using a step
        data_to_remove%x = an_integer * (p_num_of_insertions/p_num_of_deletions) - an_integer
        ! print*, 'removing ', data_to_remove%x
        dummy = remove(vec, data_to_remove)
      enddo
      call cpu_time(time_final)
      call print_all(vec)

      print *, 'Remove Max time = ', p_remove_max_time_allowed, '; Run time = ', time_final-time_initial
      if (time_final-time_initial > p_remove_max_time_allowed) then
        print *, 'Remove Max time exceeded!!!'
        test_result = .false.
        call free_memory(vec)
        return
      endif

      call free_memory(vec)
      return

    end function test_vector_remove_performance

  
end program Modular_tests