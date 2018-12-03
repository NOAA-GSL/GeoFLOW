c======================================================================================
c Name         : dmxmp
c Date         : 1/1/18 (DLR)
c Copyright    : Copyright 2018. Colorado State University. All rights reserved
c Description  : Hand-unrolled m x m routines to call from Fortran cff_blas 
c                driver routine, for GDOUBLE types
c======================================================================================




c************************************************************************************
c************************************************************************************
c METHOD     : dmxm1
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm1(A, na, B, nb, C, nc)
!$acc routine vector
      real*8     A(na,1), B(1,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i,1) * B(1,j)
        enddo 
      enddo


      return
c
c  end of subroutine dmxm1
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm2
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm2(A, na, B, nb, C, nc)
!$acc routine vector
      real*8     A(na,2), B(2,nc), C(na,nc) 
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i,1) * B(1,j)
     2           + A(i,2) * B(2,j)
        enddo 
      enddo


      return
c
c  end of subroutine dmxm2
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm3
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm3(A, na, B, nb, C, nc)
!$acc routine vector
      real*8     A(na,3), B(3,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i,1) * B(1,j)
     2           + A(i,2) * B(2,j)
     3           + A(i,3) * B(3,j)
        enddo 
      enddo


      return
c
c  end of subroutine dmxm3
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm4
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm4(A, na, B, nb, C, nc)
!$acc routine vector
      real*8     A(na,4), B(4,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i,1) * B(1,j)
     2           + A(i,2) * B(2,j)
     3           + A(i,3) * B(3,j)
     4           + A(i,4) * B(4,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm4
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm5
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm5(A, na, B, nb, C, nc)
!$acc routine vector
      real*8     A(na,5), B(5,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i,1) * B(1,j)
     2           + A(i,2) * B(2,j)
     3           + A(i,3) * B(3,j)
     4           + A(i,4) * B(4,j)
     5           + A(i,5) * B(5,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm5
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm6
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm6(A, na, B, nb, C, nc)
!$acc routine vector
      real*8     A(na,6), B(6,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i,1) * B(1,j)
     2           + A(i,2) * B(2,j)
     3           + A(i,3) * B(3,j)
     4           + A(i,4) * B(4,j)
     5           + A(i,5) * B(5,j)
     6           + A(i,6) * B(6,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm6
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm7
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm7(A, na, B, nb, C, nc)
!$acc routine vector
      real*8     A(na,7), B(7,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i,1) * B(1,j)
     2           + A(i,2) * B(2,j)
     3           + A(i,3) * B(3,j)
     4           + A(i,4) * B(4,j)
     5           + A(i,5) * B(5,j)
     6           + A(i,6) * B(6,j)
     7           + A(i,7) * B(7,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm7
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm8
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm8(A, na, B, nb, C, nc)
      real*8     A(na,8), B(8,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i,1) * B(1,j)
     2           + A(i,2) * B(2,j)
     3           + A(i,3) * B(3,j)
     4           + A(i,4) * B(4,j)
     5           + A(i,5) * B(5,j)
     6           + A(i,6) * B(6,j)
     7           + A(i,7) * B(7,j)
     8           + A(i,8) * B(8,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm8
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm9
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm9(A, na, B, nb, C, nc)
      real*8     A(na,9), B(9,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i,1) * B(1,j)
     2           + A(i,2) * B(2,j)
     3           + A(i,3) * B(3,j)
     4           + A(i,4) * B(4,j)
     5           + A(i,5) * B(5,j)
     6           + A(i,6) * B(6,j)
     7           + A(i,7) * B(7,j)
     8           + A(i,8) * B(8,j)
     9           + A(i,9) * B(9,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm9
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm10
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm10(A, na, B, nb, C, nc)
      real*8     A(na,10), B(10,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm10
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm11
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm11(A, na, B, nb, C, nc)
      real*8     A(na,11), B(11,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm11
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm12
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm12(A, na, B, nb, C, nc)
      real*8     A(na,12), B(12,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm12
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm13
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm13(A, na, B, nb, C, nc)
      real*8     A(na,13), B(13,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm13
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm14
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm14(A, na, B, nb, C, nc)
      real*8     A(na,14), B(14,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
     4           + A(i,14) * B(14,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm14
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm15
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm15(A, na, B, nb, C, nc)
      real*8     A(na,15), B(15,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
     4           + A(i,14) * B(14,j)
     5           + A(i,15) * B(15,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm15
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm16
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm16(A, na, B, nb, C, nc)
      real*8     A(na,16), B(16,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
     4           + A(i,14) * B(14,j)
     5           + A(i,15) * B(15,j)
     6           + A(i,16) * B(16,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm16
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm17
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm17(A, na, B, nb, C, nc)
      real*8     A(na,17), B(17,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
     4           + A(i,14) * B(14,j)
     5           + A(i,15) * B(15,j)
     6           + A(i,16) * B(16,j)
     7           + A(i,17) * B(17,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm17
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm18
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm18(A, na, B, nb, C, nc)
      real*8     A(na,18), B(18,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
     4           + A(i,14) * B(14,j)
     5           + A(i,15) * B(15,j)
     6           + A(i,16) * B(16,j)
     7           + A(i,17) * B(17,j)
     8           + A(i,18) * B(18,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm18
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm19
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm19(A, na, B, nb, C, nc)
      real*8     A(na,19), B(19,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
     4           + A(i,14) * B(14,j)
     5           + A(i,15) * B(15,j)
     6           + A(i,16) * B(16,j)
     7           + A(i,17) * B(17,j)
     8           + A(i,18) * B(18,j)
     9           + A(i,19) * B(19,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm19
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm20
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm20(A, na, B, nb, C, nc)
      real*8     A(na,20), B(20,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
     4           + A(i,14) * B(14,j)
     5           + A(i,15) * B(15,j)
     6           + A(i,16) * B(16,j)
     7           + A(i,17) * B(17,j)
     8           + A(i,18) * B(18,j)
     9           + A(i,19) * B(19,j)
     O           + A(i,20) * B(20,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm20
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm21
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm21(A, na, B, nb, C, nc)
      real*8     A(na,21), B(21,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
     4           + A(i,14) * B(14,j)
     5           + A(i,15) * B(15,j)
     6           + A(i,16) * B(16,j)
     7           + A(i,17) * B(17,j)
     8           + A(i,18) * B(18,j)
     9           + A(i,19) * B(19,j)
     O           + A(i,20) * B(20,j)
     1           + A(i,21) * B(21,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm21
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm22
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm22(A, na, B, nb, C, nc)
      real*8     A(na,22), B(22,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
     4           + A(i,14) * B(14,j)
     5           + A(i,15) * B(15,j)
     6           + A(i,16) * B(16,j)
     7           + A(i,17) * B(17,j)
     8           + A(i,18) * B(18,j)
     9           + A(i,19) * B(19,j)
     O           + A(i,20) * B(20,j)
     1           + A(i,21) * B(21,j)
     2           + A(i,22) * B(22,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm22
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm23
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm23(A, na, B, nb, C, nc)
      real*8     A(na,23), B(23,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
     4           + A(i,14) * B(14,j)
     5           + A(i,15) * B(15,j)
     6           + A(i,16) * B(16,j)
     7           + A(i,17) * B(17,j)
     8           + A(i,18) * B(18,j)
     9           + A(i,19) * B(19,j)
     O           + A(i,20) * B(20,j)
     1           + A(i,21) * B(21,j)
     2           + A(i,22) * B(22,j)
     3           + A(i,23) * B(23,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm23
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxm24
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm24(A, na, B, nb, C, nc)
      real*8     A(na,24), B(24,nc), C(na,nc)
      integer    na, nb, nc
      integer    i, j

c  Do matrix-matrix multiply using cache-blocking:
      do j = 1, nb
        do i = 1, na
          C(i,j) = A(i, 1) * B(1 ,j)
     2           + A(i, 2) * B(2 ,j)
     3           + A(i, 3) * B(3 ,j)
     4           + A(i, 4) * B(4 ,j)
     5           + A(i, 5) * B(5 ,j)
     6           + A(i, 6) * B(6 ,j)
     7           + A(i, 7) * B(7 ,j)
     8           + A(i, 8) * B(8 ,j)
     9           + A(i, 9) * B(9 ,j)
     O           + A(i,10) * B(10,j)
     1           + A(i,11) * B(11,j)
     2           + A(i,12) * B(12,j)
     3           + A(i,13) * B(13,j)
     4           + A(i,14) * B(14,j)
     5           + A(i,15) * B(15,j)
     6           + A(i,16) * B(16,j)
     7           + A(i,17) * B(17,j)
     8           + A(i,18) * B(18,j)
     9           + A(i,19) * B(19,j)
     O           + A(i,20) * B(20,j)
     1           + A(i,21) * B(21,j)
     2           + A(i,22) * B(22,j)
     3           + A(i,23) * B(23,j)
     4           + A(i,24) * B(24,j)
        enddo 
      enddo

      return
c
c  end of subroutine dmxm24
c
      end

