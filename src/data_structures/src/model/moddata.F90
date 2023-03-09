! A derived type for storing data.
module moddata
  implicit none

  private
  public :: data_t
!   public :: data_ptr

  ! Data is stored in data_t
  type :: data_t
     integer :: x
  end type data_t

!   ! A trick to allow us to store pointers in the list
!   type :: data_ptr
!      type(data_t), pointer :: p
!   end type data_ptr
end module moddata