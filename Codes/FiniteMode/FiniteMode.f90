program main
 implicit none
 integer, parameter :: dp = selected_real_kind(15)
 integer :: i,n
 real(dp) :: dt,t,tf
 real(dp), allocatable :: x1(:),x2(:),x3(:),x4(:),x5(:)

 dt = 0.0001
 n = 2000000
 allocate(x1(n),x2(n),x3(n),x4(n),x5(n))

 x1(1) = 0.540323
 x2(1) = -1.543569
 x3(1) = -0.680421
 x4(1) = -1.185361
 x5(1) = -0.676307
 

 open(1,file = "x1.dat")
 open(2,file = "x2.dat")
 open(3,file = "x3.dat")
 open(4,file = "x4.dat")
 open(5,file = "x5.dat")

 do i = 1,n
  x1(i+1) = x1(i)+dt*(x2(i)*x3(i)+x5(i)*x4(i)-2.0*x2(i)*x5(i))
  x2(i+1) = x2(i)+dt*(x3(i)*x4(i)+x1(i)*x5(i)-2.0*x3(i)*x1(i))
  x3(i+1) = x3(i)+dt*(x4(i)*x5(i)+x2(i)*x1(i)-2.0*x4(i)*x2(i))
  x4(i+1) = x4(i)+dt*(x5(i)*x1(i)+x3(i)*x2(i)-2.0*x5(i)*x3(i))
  x5(i+1) = x5(i)+dt*(x1(i)*x2(i)+x4(i)*x3(i)-2.0*x1(i)*x4(i))
 end do 

 write(1,*) x1
 write(2,*) x2
 write(3,*) x3
 write(4,*) x4
 write(5,*) x5
 
 deallocate(x1,x2,x3,x4,x5)

end program