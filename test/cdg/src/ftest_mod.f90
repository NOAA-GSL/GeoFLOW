module ftest

     private

     public :: doprod, dofprod, dodprod

     interface doprod
       module procedure dofprod, dodprod
     end interface

     contains

     subroutine dofprod(ret, a, b)
       implicit none
       real, intent(out) :: ret
       real, intent (in) :: a, b

       ret = a * b
     end subroutine dofprod

     subroutine dodprod(ret, a, b)
       implicit none
       real*8, intent(out) :: ret
       real*8, intent (in) :: a, b

       ret = a * b
     end subroutine dodprod

end module ftest
