c======================================================================================
c Name         : mxmp
c Date         : 1/1/18 (DLR)
c Copyright    : Copyright 2018. Colorado State University. All rights reserved
c Description  : Hand-unrolled m x m routines to call from Fortran cff_blas driver 
c                routine, for GFLOAT types
c======================================================================================




c************************************************************************************
c************************************************************************************
c METHOD     : fmxm1
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm1(A, na, B, nb, C, nc)
!$acc routine vector
      real*4     A(na,1), B(1,nc), C(na,nc)
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
c  end of subroutine fmxm1
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm2
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm2(A, na, B, nb, C, nc)
!$acc routine vector
      real*4     A(na,2), B(2,nc), C(na,nc) 
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
c  end of subroutine fmxm2
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm3
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm3(A, na, B, nb, C, nc)
!$acc routine vector
      real*4     A(na,3), B(3,nc), C(na,nc)
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
c  end of subroutine fmxm3
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm4
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm4(A, na, B, nb, C, nc)
!$acc routine vector
      real*4     A(na,4), B(4,nc), C(na,nc)
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
c  end of subroutine fmxm4
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm5
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm5(A, na, B, nb, C, nc)
!$acc routine vector
      real*4     A(na,5), B(5,nc), C(na,nc)
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
c  end of subroutine fmxm5
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm6
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm6(A, na, B, nb, C, nc)
!$acc routine vector
      real*4     A(na,6), B(6,nc), C(na,nc)
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
c  end of subroutine fmxm6
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm7
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm7(A, na, B, nb, C, nc)
!$acc routine vector
      real*4     A(na,7), B(7,nc), C(na,nc)
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
c  end of subroutine fmxm7
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm8
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm8(A, na, B, nb, C, nc)
      real*4     A(na,8), B(8,nc), C(na,nc)
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
c  end of subroutine fmxm8
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm9
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm9(A, na, B, nb, C, nc)
      real*4     A(na,9), B(9,nc), C(na,nc)
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
c  end of subroutine fmxm9
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm10
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm10(A, na, B, nb, C, nc)
      real*4     A(na,10), B(10,nc), C(na,nc)
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
c  end of subroutine fmxm10
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm11
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm11(A, na, B, nb, C, nc)
      real*4     A(na,11), B(11,nc), C(na,nc)
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
c  end of subroutine fmxm11
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm12
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm12(A, na, B, nb, C, nc)
      real*4     A(na,12), B(12,nc), C(na,nc)
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
c  end of subroutine fmxm12
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm13
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm13(A, na, B, nb, C, nc)
      real*4     A(na,13), B(13,nc), C(na,nc)
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
c  end of subroutine fmxm13
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm14
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm14(A, na, B, nb, C, nc)
      real*4     A(na,14), B(14,nc), C(na,nc)
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
c  end of subroutine fmxm14
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm15
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm15(A, na, B, nb, C, nc)
      real*4     A(na,15), B(15,nc), C(na,nc)
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
c  end of subroutine fmxm15
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm16
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm16(A, na, B, nb, C, nc)
      real*4     A(na,16), B(16,nc), C(na,nc)
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
c  end of subroutine fmxm16
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm17
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm17(A, na, B, nb, C, nc)
      real*4     A(na,17), B(17,nc), C(na,nc)
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
c  end of subroutine fmxm17
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm18
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm18(A, na, B, nb, C, nc)
      real*4     A(na,18), B(18,nc), C(na,nc)
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
c  end of subroutine fmxm18
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm19
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm19(A, na, B, nb, C, nc)
      real*4     A(na,19), B(19,nc), C(na,nc)
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
c  end of subroutine fmxm19
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm20
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm20(A, na, B, nb, C, nc)
      real*4     A(na,20), B(20,nc), C(na,nc)
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
c  end of subroutine fmxm20
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm21
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm21(A, na, B, nb, C, nc)
      real*4     A(na,21), B(21,nc), C(na,nc)
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
c  end of subroutine fmxm21
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm22
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm22(A, na, B, nb, C, nc)
      real*4     A(na,22), B(22,nc), C(na,nc)
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
c  end of subroutine fmxm22
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm23
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm23(A, na, B, nb, C, nc)
      real*4     A(na,23), B(23,nc), C(na,nc)
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
c  end of subroutine fmxm23
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : fmxm24
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine fmxm24(A, na, B, nb, C, nc)
      real*4     A(na,24), B(24,nc), C(na,nc)
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
c  end of subroutine fmxm24
c
      end

