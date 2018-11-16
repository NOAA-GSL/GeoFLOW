c======================================================================================
c Name         : qmxmp
c Date         : 1/1/18 (DLR)
c Copyright    : Copyright 2018. Colorado State University. All rights reserved
c Description  : Hand-unrolled m x m routines to call from Fortran cff_blas 
c                driver routine, for GQUAD types
c======================================================================================




c************************************************************************************
c************************************************************************************
c METHOD     : qmxm1
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm1(A, na, B, nb, C, nc)
!$acc routine vector
      real*16    A(na,1), B(1,nc), C(na,nc)
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
c  end of subroutine qmxm1
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm2
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm2(A, na, B, nb, C, nc)
!$acc routine vector
      real*16    A(na,2), B(2,nc), C(na,nc) 
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
c  end of subroutine qmxm2
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm3
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm3(A, na, B, nb, C, nc)
!$acc routine vector
      real*16    A(na,3), B(3,nc), C(na,nc)
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
c  end of subroutine qmxm3
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm4
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm4(A, na, B, nb, C, nc)
!$acc routine vector
      real*16    A(na,4), B(4,nc), C(na,nc)
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
c  end of subroutine qmxm4
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm5
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm5(A, na, B, nb, C, nc)
!$acc routine vector
      real*16    A(na,5), B(5,nc), C(na,nc)
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
c  end of subroutine qmxm5
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm6
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm6(A, na, B, nb, C, nc)
!$acc routine vector
      real*16    A(na,6), B(6,nc), C(na,nc)
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
c  end of subroutine qmxm6
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm7
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm7(A, na, B, nb, C, nc)
!$acc routine vector
      real*16    A(na,7), B(7,nc), C(na,nc)
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
c  end of subroutine qmxm7
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm8
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm8(A, na, B, nb, C, nc)
      real*16    A(na,8), B(8,nc), C(na,nc)
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
c  end of subroutine qmxm8
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm9
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm9(A, na, B, nb, C, nc)
      real*16    A(na,9), B(9,nc), C(na,nc)
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
c  end of subroutine qmxm9
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm10
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm10(A, na, B, nb, C, nc)
      real*16    A(na,10), B(10,nc), C(na,nc)
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
c  end of subroutine qmxm10
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm11
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm11(A, na, B, nb, C, nc)
      real*16    A(na,11), B(11,nc), C(na,nc)
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
c  end of subroutine qmxm11
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm12
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm12(A, na, B, nb, C, nc)
      real*16    A(na,12), B(12,nc), C(na,nc)
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
c  end of subroutine qmxm12
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm13
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm13(A, na, B, nb, C, nc)
      real*16    A(na,13), B(13,nc), C(na,nc)
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
c  end of subroutine qmxm13
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm14
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm14(A, na, B, nb, C, nc)
      real*16    A(na,14), B(14,nc), C(na,nc)
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
c  end of subroutine qmxm14
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm15
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm15(A, na, B, nb, C, nc)
      real*16    A(na,15), B(15,nc), C(na,nc)
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
c  end of subroutine qmxm15
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm16
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm16(A, na, B, nb, C, nc)
      real*16    A(na,16), B(16,nc), C(na,nc)
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
c  end of subroutine qmxm16
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm17
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm17(A, na, B, nb, C, nc)
      real*16    A(na,17), B(17,nc), C(na,nc)
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
c  end of subroutine qmxm17
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm18
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm18(A, na, B, nb, C, nc)
      real*16    A(na,18), B(18,nc), C(na,nc)
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
c  end of subroutine qmxm18
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm19
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm19(A, na, B, nb, C, nc)
      real*16    A(na,19), B(19,nc), C(na,nc)
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
c  end of subroutine qmxm19
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm20
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm20(A, na, B, nb, C, nc)
      real*16    A(na,20), B(20,nc), C(na,nc)
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
c  end of subroutine qmxm20
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm21
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm21(A, na, B, nb, C, nc)
      real*16    A(na,21), B(21,nc), C(na,nc)
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
c  end of subroutine qmxm21
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm22
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm22(A, na, B, nb, C, nc)
      real*16    A(na,22), B(22,nc), C(na,nc)
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
c  end of subroutine qmxm22
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm23
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm23(A, na, B, nb, C, nc)
      real*16    A(na,23), B(23,nc), C(na,nc)
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
c  end of subroutine qmxm23
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : qmxm24
c DESCRIPTION: mxm utility  routine that explicitly unrolls by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine qmxm24(A, na, B, nb, C, nc)
      real*16    A(na,24), B(24,nc), C(na,nc)
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
c  end of subroutine qmxm24
c
      end

