!Parallel Program
!Final Project Parallelizing Gaussian Elimination
!Solve AX = C

use omp_lib
implicit none
integer,parameter       :: nt = 16, nflush = 2 * 3 * 1024 * 1024
integer,parameter       :: m = 1024, k = m, n = m + k
double precision        :: flush(nflush) !Flushes caches
integer                 :: i, j, l, maxindex
double precision        :: A(m,m), X(m,k), XE(m,k), B(m,n), C(m,k)
double precision        :: RESIDUAL(m,k), t0, t1, time, num_ops
double precision        :: gflops, maxvalue
double precision        :: D(m,m), E(m,m), temp, reciptemp, if_zero,
largest, alpha

!Set threads
call omp_set_num_threads(nt)

!Initialize Arrays
call random_number(A)
B(1:m,1:m) = A
call random_number(XE) !Exact Solution
C = matmul(A,XE)
B(1:m,(m+1):(m+k)) = C

!Begin program

!Flush caches
!$omp parallel do schedule(static) shared(flush) private(i)
do i = 1, nflush
  flush(i) = dble(i)
enddo
!$omp end parallel do

!Start timer
t0 = omp_get_wtime()

maxvalue = 1/maxval(abs(B))

!Find pivots and zero B

do i = 1,m
  if (i<m)then
    largest = 0.d0
    do j = 1,m
      if (abs(B(j,i)) > abs(largest)) then
        largest = B(j,i)
        maxindex = j
      endif
    enddo

!interchange rows

    do j = 1,n
      temp = B(i,j)
      B(i,j) = B(maxindex,j)
      B(maxindex,j) = temp
    enddo
  endif

  if (i == m) then
    maxindex = m !IF i = m THE PIVOT IS B(i,i)
  endif
!FINDING THE RELATIVE ERROR

 if_zero = abs(B(maxindex,i))*maxvalue

  if (if_zero < 1.d-10) then
    print*,'Matrix is numerically singular'
  endif
!Divide the row by the pivot
  reciptemp = 1/B(i,i)
  do j = 1,n
    B(i,j) = B(i,j)*reciptemp
  enddo

!THE FOLLOWING LOOP ZEROS OUT THE APPROPRIATE PART OF THE AUGEMENTED
!MATRIX

!$omp parallel default(none) shared(B,alpha,i) private(l,j)
!$omp do schedule(runtime)

  do j = 1,m
    if(j .ne. i) then
      alpha = -B(j,i)
      do l = i+1,n
        B(j,l) = B(i,l) * alpha + B(j,l)
      enddo
      B(j,i) = 0.d0
    endif
  enddo

!$omp enddo
!$omp end parallel

Enddo

t1 = omp_get_wtime()
time = t1 - t0 !time in seconds

!ESTIMATED NUMBER OF OPERATIONS AND PERFORMANCE
num_ops = dble(k)*dble(m)*dble(m)*dble(m)/3.d0  !ESTIMATE
gflops = (num_ops/time)*1.d-9
print*,'Estimated performance is ', gflops,' gflops'
print*,'Time = ', time,' seconds'

!Compute the solution
X(1:m,1:k) = B(1:m,(m+1):(m+k))

!FIND AND PRINT ERROR
print*,'m = ',m,' k = ',k,' number of right hand sides'
RESIDUAL = matmul(A,X) - C
print*,'RESIDUAL = ', maxval(abs(RESIDUAL))
print*,'Maximum error ',maxval(abs(X-XE))
print*,'flush(8) + 2 = ',flush(8)+2 !dead code elim
print*,'check check'
end

