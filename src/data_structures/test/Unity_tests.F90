! A simple generic linked list test program
program Unity_tests

  use Moddata

  implicit none

  integer :: test_errors_sum
  
  print*, ""
  print*, ">>>>> Running Unity Tests"

  test_errors_sum = 0
  test_errors_sum = test_errors_sum + test_list_init()
  test_errors_sum = test_errors_sum + test_list_inserts()
  test_errors_sum = test_errors_sum + test_list_removes()
  test_errors_sum = test_errors_sum + test_vector_init()
  test_errors_sum = test_errors_sum + test_vector_inserts()
  test_errors_sum = test_errors_sum + test_vector_removes()
  
  print *, "" 
  if (test_errors_sum == 0) then
    print*, ">>>>> ALL TESTS OK !"
    call exit(0)
  else
    print*, ">>>>> TESTS FAILED. TOTAL TESTS FAILED = ", test_errors_sum
    call exit(-1)
  endif


  contains 

    ! List Tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !
    integer function test_list_init() result(test_error)
      use Modlist
      implicit none 

      type(list_t), pointer :: ll => null()
      type(data_t), allocatable :: dat_to_insert
      type(data_t) :: dat_test

      print*, ""
      print*, ">>>>> Running Test List Initilize"
      test_error = 0

      allocate(dat_to_insert)
      dat_to_insert%x = 1

      call init(ll, dat_to_insert)
      print *, 'Initializing list with data:', dat_to_insert

      print *, 'Testing head node'      
      dat_test = get(ll)
      if (dat_test%x /= 1) then
        print *, 'Head node data should be: 1 but was', dat_test%x
        test_error = 1
      endif
      
      call free_memory(ll)
      return

    end function


    integer function test_list_inserts() result(test_error)
      use Modlist
      implicit none 

      type(list_t), pointer :: ll => null()
      type(data_t), allocatable :: dat_x
      type(data_t) :: dat_test

      print*, ""
      print*, ">>>>> Running Test List Insert second element"
      test_error = 0

      allocate(dat_x)
      dat_x%x = 1
      call init(ll, dat_x)
      print *, 'Initializing list with data:', dat_x%x
      deallocate(dat_x)

      allocate(dat_x)
      dat_x%x = 2
      call insert(ll, dat_x)
      print *, 'Inserting node with data:', dat_x%x
      deallocate(dat_x)

      print *, 'Testing head node'      
      dat_test = get(ll)
      if (dat_test%x /= 2) then
        print *, '!!!!! TEST FAILED !!!!! Head node data should be: 2 but was', dat_test%x
        call free_memory(ll)
        test_error = 1
        return
      endif
      
      print *, 'Testing second node'
      dat_test = get(next(ll))
      if (dat_test%x /= 1) then
        print *, '!!!!! TEST FAILED !!!!! Second node data should be: 1 but was', dat_test%x
        call free_memory(ll)
        test_error = 1
        return
      endif

      ! Free the list
      call free_memory(ll)

    end function

    integer function test_list_removes() result(test_error)
      use Modlist
      implicit none 

      type(list_t), pointer :: ll, node_curr  => null()
      type(data_t) :: dat_x
      type(data_t) :: dat_test
      integer :: test_value
      logical :: is_removed 

      print*, ""
      print*, ">>>>> Running Test List removes elements"
      test_error = 0

      dat_x%x = 10
      call init(ll, dat_x)
      print *, 'Initializing list with data:', dat_x%x
      dat_x%x = 20
      call insert(ll, dat_x)
      print *, 'Inserting node with data:', dat_x%x
      dat_x%x = 30
      call insert(ll, dat_x)
      print *, 'Inserting node with data:', dat_x%x
      dat_x%x = 40
      call insert(ll, dat_x)
      print *, 'Inserting node with data:', dat_x%x

      print *, 'removes last element'
      dat_x%x = 10
      is_removed = remove(ll, dat_x)

      print *, 'Testing nodes'      
      node_curr => ll
      test_value = 40
      do
        dat_test = get(node_curr)
        print *, 'Checking node ', dat_test%x
        if (dat_test%x /= test_value) then
          print *, '!!!!! TEST FAILED !!!!! Node data should be: ',test_value,' but was', dat_test%x
          call free_memory(ll)
          test_error = 1
          return
        endif
        if (.not. associated(next(node_curr))) exit
        node_curr => next(node_curr)
        test_value = test_value - 10
      enddo

      if(test_value == 10 ) then
        print *, '!!!!! TEST FAILED !!!!! last element not removed'
        call free_memory(ll)
        test_error = 1
        return
      endif
      
      print *, 'removes second element'
      dat_x%x = 30
      is_removed = remove(ll, dat_x)
      node_curr => ll
      dat_test = get(node_curr)
      print *, 'Checking node ', dat_test%x
      test_value = 40
      if (dat_test%x /= test_value) then
        print *, '!!!!! TEST FAILED !!!!! Node data should be: ',test_value,' but was', dat_test%x
        call free_memory(ll)
        test_error = 1
        return
      endif
      node_curr => next(node_curr)
      dat_test = get(node_curr)
      print *, 'Checking node ', dat_test%x
      test_value = 20
      if (dat_test%x /= test_value) then
        print *, '!!!!! TEST FAILED !!!!! Node data should be: ',test_value,' but was', dat_test%x
        call free_memory(ll)
        test_error = 1
        return
      endif

      print *, 'removes first element'
      dat_x%x = 40
      is_removed = remove(ll, dat_x)
      node_curr => ll
      dat_test = get(node_curr)
      print *, 'Checking node ', dat_test%x
      test_value = 20
      if (dat_test%x /= test_value) then
        print *, '!!!!! TEST FAILED !!!!! Node data should be: ',test_value,' but was', dat_test%x
        call free_memory(ll)
        test_error = 1
        return
      endif

      print *, 'removes the only one element '
      dat_x%x = 20
      is_removed = remove(ll, dat_x)
      if (associated(ll)) then
        print *, '!!!!! TEST FAILED !!!!! List should not contains elements'
        call free_memory(ll)
        test_error = 1
        return
      endif

    end function

    ! Vector Tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !
    integer function test_vector_init() result(test_error)
      use ModVector
      implicit none 

      type(vector_t) :: vec
      integer, parameter :: p_vector_size = 1000000

      print*, ""
      print*, ">>>>> Running Test Vector Initilize"
      test_error = 0


      print *, 'Testing get_size = 0'
      if(get_size(vec) /= 0) then
        print *, '!!!!! TEST FAILED !!!!! Size of vector should be: 0 but was', get_size(vec)
        test_error = 1
        call free_memory(vec)
        return
      endif

      print *, 'Initializing Vector with size = ', p_vector_size
      call init(vec, p_vector_size)

      print *, 'Testing get_size = ', p_vector_size
      if(get_size(vec) /= p_vector_size) then
        print *, '!!!!! TEST FAILED !!!!! Size of vector should be: ', p_vector_size,' but was', get_size(vec)
        test_error = 1
        call free_memory(vec)
        return
      endif

      print *, 'Testing num_elements = 0'
      if(get_num_elements(vec) /= 0) then
        print *, '!!!!! TEST FAILED !!!!! Num elements should be: 0 but was', get_num_elements(vec)
        test_error = 1
        call free_memory(vec)
        return
      endif

      call free_memory(vec)

    end function


    integer function test_vector_inserts() result(test_error)
      use ModVector
      implicit none 

      type(vector_t) :: vec
      type(data_t) :: dat_to_insert
      type(data_t) :: dat_test
      integer, parameter :: p_vector_size = 1000000
      
      print*, ""
      print*, ">>>>> Running Test Vector Insert elements"
      test_error = 0

      print *, 'Initializing Vector with size = ', p_vector_size
      call init(vec, p_vector_size)

      dat_to_insert%x = 10
      print *, 'Inserting vector element with: ', dat_to_insert%x
      call insert(vec, dat_to_insert)

      print *, 'Testing num_elements = 1'
      if(get_num_elements(vec) /= 1) then
        print *, '!!!!! TEST FAILED !!!!! Num elements should be: 1 but was', get_num_elements(vec)
        test_error = 1
        call free_memory(vec)
        return
      endif

      print *, 'Testing first element'
      dat_test = get(vec, get_num_elements(vec))
      if(dat_test%x /= 10) then
        print *, '!!!!! TEST FAILED !!!!! First element should be: 10 but was', dat_test%x
        test_error = 1
        call free_memory(vec)
        return
      endif

      dat_to_insert%x = 20
      print *, 'Inserting vector element with: ', dat_to_insert%x
      call insert(vec, dat_to_insert)

      print *, 'Testing num_elements = 2'
      if(get_num_elements(vec) /= 2) then
        print *, '!!!!! TEST FAILED !!!!! Num elements should be: 2 but was', get_num_elements(vec)
        test_error = 1
        call free_memory(vec)
        return
      endif

      print *, 'Testing second element'
      dat_test = get(vec, get_num_elements(vec))
      if(dat_test%x /= 20) then
        print *, '!!!!! TEST FAILED !!!!! First element should be: 20 but was', dat_test%x
        test_error = 1
        call free_memory(vec)
        return
      endif

      call free_memory(vec)

    end function test_vector_inserts


    integer function test_vector_removes() result(test_error)
      use ModVector
      implicit none 

      type(vector_t) :: vec
      type(data_t) :: dat_to_insert, dat_to_remove
      type(data_t) :: dat_test
      integer, parameter :: p_vector_size = 1000000
      integer :: num_elements_test, index_element, test_value
      logical :: dummy
      
      print*, ""
      print*, ">>>>> Running Test Vector Removes elements"
      test_error = 0

      print *, 'Initializing Vector with size = ', p_vector_size
      call init(vec, p_vector_size)

      dat_to_insert%x = 10
      print *, 'Insertting elements in vector:' 
      call insert(vec, dat_to_insert)
      dat_to_insert%x = 20
      call insert(vec, dat_to_insert)
      dat_to_insert%x = 30
      call insert(vec, dat_to_insert)
      dat_to_insert%x = 40
      call insert(vec, dat_to_insert)
      call print_all(vec)

      print *, 'Testing remove last '
      dat_to_remove%x = 40
      dummy = remove(vec, dat_to_remove)
      num_elements_test = 3
      print *, 'Testing num_elements = ', num_elements_test
      if(get_num_elements(vec) /= num_elements_test) then
        print *, '!!!!! TEST FAILED !!!!! Num elements should be: ', num_elements_test, ' but was ', get_num_elements(vec)
        test_error = 1
        call free_memory(vec)
        return
      endif
      print *, 'Testing elements'
      call print_all(vec)
      test_value = 10
      do index_element = 1, 3
        dat_test = get(vec, index_element)
        if(dat_test%x /= test_value) then
          print *, '!!!!! TEST FAILED !!!!! element index', index_element, ' should be: ', test_value, ' but was', dat_test%x
          test_error = 1
          call free_memory(vec)
          return
        endif
        test_value = test_value + 10
      enddo

      print *, 'Testing remove index 2'
      dat_to_remove%x = 20
      dummy = remove(vec, dat_to_remove)
      num_elements_test = 2
      print *, 'Testing num_elements = ', num_elements_test
      if(get_num_elements(vec) /= num_elements_test) then
        print *, '!!!!! TEST FAILED !!!!! Num elements should be: ', num_elements_test, ' but was ', get_num_elements(vec)
        test_error = 1
        call free_memory(vec)
        return
      endif
      print *, 'Testing elements'
      call print_all(vec)
      test_value = 10
      dat_test = get(vec, 1)
      if(dat_test%x /= test_value) then
        print *, '!!!!! TEST FAILED !!!!! First element should be: ', test_value, ' but was', dat_test%x
        test_error = 1
        call free_memory(vec)
        return
      endif
      test_value = 30
      dat_test = get(vec, 2)
      if(dat_test%x /= test_value) then
        print *, '!!!!! TEST FAILED !!!!! First element should be: ', test_value, ' but was', dat_test%x
        test_error = 1
        call free_memory(vec)
        return
      endif

      print *, 'Testing remove index 1'
      dat_to_remove%x = 10
      dummy = remove(vec, dat_to_remove)
      num_elements_test = 1
      print *, 'Testing num_elements = ', num_elements_test
      if(get_num_elements(vec) /= num_elements_test) then
        print *, '!!!!! TEST FAILED !!!!! Num elements should be: ', num_elements_test, ' but was ', get_num_elements(vec)
        test_error = 1
        call free_memory(vec)
        return
      endif
      print *, 'Testing elements'
      call print_all(vec)
      test_value = 30
      dat_test = get(vec, 1)
      if(dat_test%x /= test_value) then
        print *, '!!!!! TEST FAILED !!!!! First element should be: ', test_value, ' but was', dat_test%x
        test_error = 1
        call free_memory(vec)
        return
      endif

      call free_memory(vec)

    end function test_vector_removes

end program Unity_tests