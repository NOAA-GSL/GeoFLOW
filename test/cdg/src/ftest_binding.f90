     subroutine dofprod(ret, a, b, n) bind(C,name="dofprod")
       use, intrinsic :: ISO_C_BINDING
       implicit none
       real (C_FLOAT), intent(inout), dimension(*) :: ret
       real (C_FLOAT), intent(inout), dimension(*) :: a, b
       integer(C_INT), intent(inout)               :: n
       integer                                     :: j

       do j = 1, n 
         ret(j) = a(j) * b(j)
       enddo
     end subroutine dofprod

     subroutine dodprod(ret, a, b, n) bind(C,name="dodprod")
       use, intrinsic :: ISO_C_BINDING
       implicit none
       real (C_DOUBLE), intent(inout), dimension(*) :: ret
       real (C_DOUBLE), intent(inout), dimension(*) :: a, b
       integer(C_INT), intent (inout)               :: n
       integer                                      :: j

       do j = 1, n 
         ret(j) = a(j) * b(j)
       enddo
     end subroutine dodprod

