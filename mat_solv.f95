!  Matrix solving various routines by Alan Rominger
!   version 0.3
!
!  *** version notes ***
!  0.1  Version was optimized to miminize code length
!  0.2  adding implicit none statements and expanding code in order
!       to be more effective
!  0.3s implimenting assumed size arrays so won't pass the argument
!       n, which is the order of the system.  Was forced to add 
!       interfaces to do this, mat_solv_interfaces.
!       reduce # of flops with some intermediary values.
!  0.3  reverted back to passing system order, because the syntax
!       procedure (name) :: is a 2003 feature and won't work with
!       g95 and probably several others.
!
!   Items I would like to use (some later version)
!    - where, elsewhere, end where
!    - array of pointers
!    - automatic operations on data use

program matrix_solver
  implicit none
  
  procedure (coef_s_x) :: print_system, gauss_alg, gausspp, gausspp_num
  procedure (mat_mat) :: matinv
  
  integer, parameter :: n = 3, m = 2
  double precision :: in_mat(n+m,n+m), s(n+m), x(n+m), invt(n+m,n+m)
  integer :: i,j
 
  1002 Format(1/,7x,'------------',2x,a,2x,'------------')
  
  write(*,*) ' start seed ',rand(5)
  do i = 1,n+m
    do j = 1,n+m
      in_mat(i,j) = rand()
    end do
    if (i<=n) in_mat(i,i+2:n) = 0.0d0
    if (i<=n) in_mat(i,1:i-2) = 0.0d0
    s(i) = rand()*2.0
  end do
  
  call block_alg(in_mat,s,x,n,m)
  write(*,1002) 'Block Thomas Solve'
  call print_system(in_mat,s,x,n+m)
  
  call gausspp(in_mat,s,x,n+m)
  write(*,1002) 'Gauss with Piv'
  call print_system(in_mat,s,x,n+m)

  call gausspp_num(in_mat,s,x,n+m)
  write(*,1002) 'Gauss with numerical correction'
  call print_system(in_mat,s,x,n+m)
  
  call gauss_alg(in_mat,s,x,n+m)
  write(*,1002) 'Gauss alg bare'
  call print_system(in_mat,s,x,n+m)
  
  call block_LT_alg(in_mat,s,x,n,m)
  write(*,1002) 'Block LT solve'
  call print_system(in_mat,s,x,n+m)
  
!   call jacobi_alg(in_mat,s,x,n+m)
!   write(*,1002) 'Jacobi algorithm'
!   call print_system(in_mat,s,x)
  

end program matrix_solver

subroutine print_system(coef,s,x,n)  !  generic interface
  implicit none
  integer, intent(in) :: n
  double precision, dimension(n,n), intent(in) :: coef
  double precision, dimension(n),   intent(in) :: s
  double precision, dimension(n),   intent(out):: x
  integer :: i

  1000 Format(5(F9.5),2('   |  ',F10.5))
  1001 Format(5(/5(F9.5)))  !  For printing matrix
  do i = 1,n
    write(*,1000) coef(i,1:n),s(i),x(i)
  end do
  x=0.
  
end subroutine print_system

!   **********************************************************************************
!  Solves the 2 eqn matrix system (D-CA*B)f=(S2-CA-1S1), e=A*S1-A*Bf where A is lower triangular
subroutine block_LT_alg(coef,s,x,n,m)
  implicit none
  integer, intent(in) :: n,m
  double precision, intent(in) :: coef(n+m,n+m), s(n+m)
  double precision, intent(out) :: x(n+m)
  double precision :: fAS(m,m+1), AnBS(n,m+1), syst(n+m,n+m+1)
  double precision :: pass1(1:n), pass2(1:m-1)
  integer :: i,j
  procedure (coef_s_x) :: LT_solve, gausspp
  syst(1:n+m,n+m+1) = s
  syst(1:n+m,1:n+m) = coef
  do j = 1,m+1
    call LT_solve(syst(1:n,1:n),syst(1:n,n+j),pass1,n)
    AnBS(1:n,j) = pass1
  end do
  do i = 1,m
    do j = 1,m+1
      fAS(i,j) = syst(n+i,n+j) - sum(coef(n+i,1:n)*AnBS(1:n,j))
    end do
  end do
  call gausspp(fAS(1:m,1:m),fAS(1:m,m+1),pass2,m)
  x(n+1:n+m) = pass2
  do i = 1,n
    x(i) = AnBS(i,m+1) - sum(AnBS(i,1:m)*x(n+1:n+m))
  end do
end subroutine block_LT_alg

!   **********************************************************************************
!  Solves the 2 eqn matrix system (D-CA*B)f=(S2-CA-1S1), e=A*S1-A*Bf where A is tridiag
subroutine block_alg(coef,s,x,n,m)
  implicit none
  integer, intent(in) :: n,m
  double precision, intent(in)  :: coef(n+m,n+m), s(n+m)
  double precision, intent(out) :: x(n+m)
  double precision :: fAS(m,m+1), AnBS(n,m+1), syst(n+m,n+m+1)
  integer :: i,j
  procedure (coef_s_x) :: thomas_alg, gausspp
  syst(1:n+m,n+m+1) = s
  syst(1:n+m,1:n+m) = coef
  do j = 1,m+1
    call thomas_alg(syst(1:n,1:n),syst(1:n,n+j),AnBS(1:n,j),n)
  end do
  do i = 1,m
    do j = 1,m+1
      fAS(i,j) = syst(n+i,n+j) - sum(coef(n+i,1:n)*AnBS(1:n,j))
    end do
  end do
  call gausspp(fAS(1:m,1:m),fAS(1:m,m+1),x(n+1:n+m),m)
  do i = 1,n
    x(i) = AnBS(i,m+1) - sum(AnBS(i,1:m)*x(n+1:n+m))
  end do
end subroutine block_alg

!   **********************************************************************************
!  Solves tridiagonal system
subroutine thomas_alg(coef,s,x,n)
  implicit none
  integer, intent(in) :: n
  double precision, dimension(n,n), intent(in) :: coef
  double precision, dimension(n),   intent(in) :: s
  double precision, dimension(n),   intent(out):: x
  double precision, dimension(n) :: alpha, g
  integer :: i,j
  alpha(1) = coef(1,1)
  g(1) = s(1)/alpha(1)
  do i = 2,n
    alpha(i) = coef(i,i)-coef(i,i-1)*(coef(i-1,i)/alpha(i-1))
    g(i) = (s(i)-coef(i,i-1)*g(i-1))/alpha(i)
  end do
  x(n) = g(n)
  do i = n-1, 1, -1
    x(i) = g(i)-(coef(i,i+1)/alpha(i))*x(i+1)
  end do
end subroutine thomas_alg

!   **********************************************************************************
!  Function applies Gauss algorithm as a part of the matrix solving code collection
!    This is adapted specifically to a lower triagonal matrix, the U part minus one diagonal
!    is assumed to be all zero.
subroutine LT_solve(coef,s,x,n)
  implicit none
  integer, intent(in) :: n
  double precision, dimension(n,n), intent(in) :: coef
  double precision, dimension(n),   intent(in) :: s
  double precision, dimension(n),   intent(out):: x
  double precision :: syst(n,n+1)
  integer :: i,j
  syst(1:n,n+1) = s
  syst(1:n,1:n) = coef
  do j = 1,n-1    ! Collumn
    do i = j+1,n    ! Row
      if (syst(i,j) == 0.0d0) cycle
      syst(i,j+1:n+1) = syst(i,j+1:n+1)-syst(j,j+1:n+1)*(syst(i,j)/syst(j,j))  ! row add
    end do
  end do
  x(n) = syst(n,n+1) / syst(n,n)
  do i = n-1,1,-1
    x(i) = (syst(i,n+1) - x(i+1)*syst(i,i+1)) / syst(i,i)
  end do
end subroutine LT_solve

!   **********************************************************************************
subroutine gauss_alg(coef,s,x,n)             !  Gauss Alg
  implicit none
  integer, intent(in) :: n
  double precision, dimension(n,n), intent(in) :: coef
  double precision, dimension(n),   intent(in) :: s
  double precision, dimension(n),   intent(out):: x
  double precision, dimension(n,n+1) :: syst
  integer :: i, j
  double precision :: alpha
  syst(1:n,n+1) = s
  syst(1:n,1:n) = coef
  do j = 1,n      ! Collumn
    if (syst(j,j)==0.0d0) stop ' Gauss failure; singular matrix'
    do i = j+1,n    ! Row
      alpha = syst(i,j)/syst(j,j)
      if (syst(i,j) == 0.0d0) cycle
      syst(i,j+1:n+1) = syst(i,j+1:n+1)-syst(j,j+1:n+1)*alpha  ! row add
    end do
  end do
  do i = n,1,-1
    x(i) = (syst(i,n+1) - sum(x(i+1:n)*syst(i,i+1:n))) / syst(i,i)
  end do
end subroutine gauss_alg

!   **********************************************************************************
!  Gaussian elimination with partial pivoting
subroutine gausspp(coef,s,x,n)
  implicit none
  integer, intent(in) :: n
  double precision, dimension(n,n), intent(in) :: coef
  double precision, dimension(n),   intent(in) :: s
  double precision, dimension(n),   intent(out):: x
  double precision :: rowstore(n+1), syst(n,n+1)
  integer :: i,j, istore
  double precision :: alpha
  syst(1:n,n+1) = s
  syst(1:n,1:n) = coef
  do j=1,n
    istore = maxloc( abs(syst(j:n,j)) ,dim=1) + j-1
    if (istore /= j) then
      rowstore(j:n+1) = syst(j,j:n+1)
      syst(j,j:n+1) = syst(istore,j:n+1)
      syst(istore,j:n+1) = rowstore(j:n+1)
    end if
    if (syst(j,j)==0.0d0) stop ' Gauss failure; singular matrix'
    do i=j+1,n
      alpha = syst(i,j)/syst(j,j)
      if (syst(i,j) == 0.0d0) cycle
      syst(i,j+1:n+1) = syst(i,j+1:n+1) - syst(j,j+1:n+1)*alpha
    end do
  end do
  do i = n,1,-1
    x(i) = (syst(i,n+1) - sum(x(i+1:n)*syst(i,i+1:n))) / syst(i,i)
  end do
end subroutine gausspp


!   **********************************************************************************
!  Gaussian alg with P.P. and numerical correction
recursive subroutine gausspp_num(coef,s,x,n)
  implicit none
  integer, intent(in) :: n
  double precision, dimension(n,n), intent(in) :: coef
  double precision, dimension(n),   intent(in) :: s
  double precision, dimension(n),   intent(out):: x
  double precision :: rowstore(n+1), syst(n,n+1), error(n), x2(n), s_trans(n,n)
  integer :: i,j, istore, itrs
  double precision :: alpha, tol, errval, errcom
  tol = 1.0e-5
  syst(1:n,n+1) = s
  syst(1:n,1:n) = coef
  s_trans = 0.
  do i=1,n
    s_trans(i,i) = 1.
  end do
  do j=1,n-1
    istore = maxloc( abs(syst(j:n,j)) ,dim=1) + j-1
    if (istore /= j) then
      rowstore(j:n+1) = syst(j,j:n+1)
      syst(j,j:n+1) = syst(istore,j:n+1)
      syst(istore,j:n+1) = rowstore(j:n+1)
      rowstore(1:n) = s_trans(j,1:n) ! store s_trans for iteration
      s_trans(j,1:n) = s_trans(istore,1:n)
      s_trans(istore,1:n) = rowstore(1:n)
    end if
    if (syst(j,j)==0.0d0) stop ' Gauss failure; singular matrix'
    do i=j+1,n
      alpha = syst(i,j)/syst(j,j)
      if (syst(i,j) == 0.0d0) cycle
      syst(i,j+1:n+1) = syst(i,j+1:n+1) - syst(j,j+1:n+1)*alpha
      s_trans(i,1:n) = s_trans(i,1:n) - s_trans(j,1:n)*alpha
    end do
  end do
  do i = n,1,-1
    x(i) = (syst(i,n+1) - sum(x(i+1:n)*syst(i,i+1:n))) / syst(i,i)
  end do
  errcom = tol * sum(abs(s))
  errval = 0.
  do i = 1,n
    error(i) = sum(coef(i,1:n)*x(1:n))-s(i)
    errval = errval + abs(error(i))
  end do
  itrs = 0
  if (errval > errcom) then
    itrs = itrs+1
    call gausspp_rec(error,x2,errcom)
    x = x - x2*0.5
!     if (itrs>5) exit
    errval = 0.
    do i = 1,n
      error(i) = sum(coef(i,1:n)*x(1:n))-s(i)
      errval = errval + abs(error(i))
    end do
  end if
  contains
    recursive subroutine gausspp_rec(s1,x,ref)
      double precision, intent(in) :: s1(n)
      double precision, intent(out):: x(n), ref
      double precision, dimension(n) :: x2, error, s2
      double precision :: errval
      do i = n,1,-1
        s2(i) = sum(s_trans(i,1:n)*s1(1:n))
        x(i) = (s2(i) - sum(x(i+1:n)*syst(i,i+1:n))) / syst(i,i)
      end do
      errval = 0.
      do i = 1,n
        error(i) = sum(coef(i,1:n)*x(1:n))-s1(i)
        errval = errval + abs(error(i))
      end do
!      write(*,'(a,2E15.5)') '   call_rec ',errval/n,errcom/n
      if (errval > ref) then
        call gausspp_rec(error,x2,ref)
        x = x - x2
      end if
    end subroutine gausspp_rec
end subroutine gausspp_num

!   **********************************************************************************
!  Matrix inversion by classic method with partial pivoting
subroutine matinv(coef,inver,n)
  implicit none
  integer, intent(in) :: n
  double precision, dimension(n,n), intent(in) :: coef
  double precision, dimension(n,n), intent(out):: inver
  double precision :: syst(n,2*n), rowstore(2*n)
  integer :: istore, i, j
  syst(1:n,1:n) = coef(1:n,1:n)
  syst(1:n,n+1:2*n) = 0.
  forall (i = 1:n) syst(i,i+n) = 1.
  do j=1,n
    istore = maxloc( abs(syst(j:n,j)) ,dim=1) + j-1
    if (istore /= j) then
      rowstore(j:2*n) = syst(j,j:2*n)
      syst(j,j:2*n) = syst(istore,j:2*n)
      syst(istore,j:2*n) = rowstore(j:2*n)
    end if
    if (syst(j,j)==0.0d0) stop ' Inverting failure; singular matrix'
    syst(j,j+1:2*n) = syst(j,j+1:2*n)/syst(j,j)
    do i=j+1,n
      if (syst(i,j) == 0.0d0) cycle
      syst(i,j+1:2*n) = syst(i,j+1:2*n) - syst(j,j+1:2*n)*syst(i,j)
    end do
  end do
  do i = n-1,1,-1
    do j = i+1,n
      if (syst(i,j) == 0.0d0) cycle
      syst(i,n+1:2*n) = syst(i,n+1:2*n) - syst(j,n+1:2*n)*syst(i,j) 
    end do
  end do
  inver = syst(1:n,n+1:2*n)
end subroutine matinv

! !   **********************************************************************************
! subroutine jacobi_alg(coef,s,x,n)             !  Jacobi Alg
!   integer          :: n, i, j, k
!   double precision :: coef(n,n), s(n), x(n), xkm(n), eps=1.e-7, omega
!   omega = 0.1
!   x = 0.1d0; k=0
!   do
!     x = x/sum(x)
!     xkm = x
!     do i = 1,n
!       x(i) = (1.d0-omega)*xkm(i)+omega*(s(i)-sum(x(1:i-1)*coef(i,1:i-1))-sum(xkm(i+1:n)*coef(i,i+1:n)))/coef(i,i)
!     end do
! !    write(*,*) ' x ',x
!     if (maxval(abs((x-xkm*sum(x))/x)) < eps) exit
!     k = k+1;  if (k>500) stop '  Jacobi mat solve failed to converge'
!   end do
! end subroutine jacobi_alg