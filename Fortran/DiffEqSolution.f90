module methods
implicit none

    contains
    ! == A subroutine that reads and outputs data from file to console      ==
    ! == Подпрограмма, выполняющая чтение и вывод исходных данных в консоль ==
    subroutine ReadAndPrintData(F,m,k,v0,h,n,method)
    implicit none
    
        real       :: F, m, k, h, v0
        integer    :: n, method
        
        ! == Reading data from file         ==
        ! == Чтение входных данных из файла ==
        open(1, file='Data.txt')

        read(1,*) F, m, k, v0, h, n, method
        close(1)
        print*, 'File read successfully.'
        
        ! == Outputting read data to the console ==
        ! == Вывод прочитанных данных в консоль  == 
        write(*,'(a5, a8, a8, a10)') 'F', 'm', 'k', 'step'
        write(*,'(f8.2,f8.2,f8.2,f8.4)') F,m,k,h
        write(*,*) ''
    end subroutine ReadAndPrintData
    
    !== The function that performs the Euler's method  ==
    !== Функция, выполняющая вычисление по методу Эйлера ==
    function Euler(F,m,k,h,n) result (v)
        real, allocatable :: v(:)
        real              :: F,m,k,h
        integer           :: i, n
        
        allocate(v(n+1))
        
        ! == Solution via Euler's method ==
        ! == Решение методом Эйлера      == 
        print*, 'Finding velocity via Eulers method...'
        v(1) = 0
        do i=2,n
            v(i) = v(i-1) + h*fxy(F,m,k,v(i-1))
        end do
    end function Euler
    
    !== The function that performs the Runge-Kutta 2nd Order method        ==
    !== Функция, выполняющая вычисление по методу Рунге-Кутты 2-го порядка ==
    function Runge_Kutta_2nd_Order(F,m,k,h,n) result(v)
        real, allocatable :: v(:)
        real              :: F,m,k,h
        real              :: k1, k2
        integer           :: i, n
        
        allocate(v(n+1))
        
        ! == Solution via Runge-Kutta 2nd Order method ==
        ! == Решение методом Рунге-Кутты 2-го порядка  ==
        print*, 'Finding velocity via Runge-Kutta 2nd Order method...'
        v(1) = 0
        do i=2,n
            k1 = h*fxy(F,m,k,v(i-1))
            k2 = h*fxy(F,m,k,(v(i-1)+(k1/2)))
            v(i) = v(i-1) + k2
        end do

    end function Runge_Kutta_2nd_Order
    
    !== The function that performs the Runge-Kutta 4th Order method        ==
    !== Функция, выполняющая вычисление по методу Рунге-Кутты 4-го порядка ==
    function Runge_Kutta_4th_Order(F,m,k,h,n) result(v)
        real, allocatable :: v(:)
        real              :: F,m,k,h
        real              :: k1, k2, k3, k4
        integer           :: i, n
        
        allocate(v(n+1))
        
        ! == Solution via Runge-Kutta 4th Order method ==
        ! == Решение методом Рунге-Кутты 4-го порядка  ==
        print*, 'Finding velocity via Runge-Kutta 4th Order method...'
        v(1) = 0
        do i=2,n
            k1 = h*fxy(F,m,k,v(i-1))
            k2 = h*fxy(F,m,k,(v(i-1)+(k1/2)))
            k3 = h*fxy(F,m,k,(v(i-1)+(k2/2)))
            k4 = h*fxy(F,m,k,(v(i-1)+k3))
            v(i) = v(i-1) + (k1+2*k2+2*k3+k4)/6
        end do
    end function Runge_Kutta_4th_Order
    
    ! == Function declaration   ==
    ! == Описываем нашу функцию == 
    real function fxy(F,m,k,v)
    implicit none
        real, intent(in) :: F,m,k,v
        fxy = (F-k*(v**2))/m
    end function fxy
    
end module methods

module percent
implicit none
    contains
    
    subroutine PercentCounter(F,m,k,h,n,v2)
    implicit none
    real, allocatable :: v1(:), v2(:), t(:)
    real, allocatable :: x(:)
    real, allocatable :: percent(:)
    real :: F, m, k, h, v0
    integer i,n

    allocate(v1(n))
    allocate(x(n))
    allocate(percent(n))
    allocate(t(n))

    t(1) = 0
    do i=1,n
        t(i) = (i-1)*h
        v1(i) = sqrt(F/k) * tanh( abs(((1/2)*log(abs(((1+(sqrt(k/F)*v0)))/(1-(sqrt(k/F)*v0))))) - (sqrt(k*F)/m)*t(i) ))
    end do

    do i=1,n
        x(i) = abs(v1(i) - v2(i))
        percent(i) = 100*x(i)/v2(i)
    end do

    open(3, file='Percents.txt', status='replace')
    do i=1,n
        write(3,'(f12.10)') percent(i)
    end do
    close(3)
    
    deallocate(v1)
    deallocate(x)
    deallocate(percent)
    deallocate(t)

    end subroutine PercentCounter
    
end module percent
    
module OrderOfAccuracy
use methods
use percent
implicit none
      
      contains
    
      subroutine OrderOfAccuracyEstimation()
      real       :: F, m, k, h1, h2, v0
      integer    :: n1, n2, method
      real, allocatable :: v_h1(:), v_h2(:)
      
      open(4, file='DataForOrderEstimation.txt')
      read(4,*) F, m, k, v0, h1, h2, n1, n2, method
      close(4)
      
      allocate(v_h1(n1+1))
      allocate(v_h2(n2+1))
      
      if (method==1) then
        v_h1 = Euler(F,m,k,h1,n1)
        v_h2 = Euler(F,m,k,h2,n2)
        
      else if (method==2) then
        v_h1 = Runge_Kutta_2nd_Order(F,m,k,h1,n1)
        v_h2 = Runge_Kutta_2nd_Order(F,m,k,h2,n2)
    
      else if (method==3) then
        v_h1 = Runge_Kutta_4th_Order(F,m,k,h1,n1)
        v_h2 = Runge_Kutta_4th_Order(F,m,k,h2,n2)
        
      else
        print*, 'Error in selecting method'
        pause
      end if
      
      call PercentCounter(F,m,k,h1,n1,v_h1)
      call PercentCounter(F,m,k,h2,n2,v_h2)
      
      end subroutine 
end module OrderOfAccuracy
    
program main
use methods
use percent
use OrderOfAccuracy
implicit none

real, allocatable :: v(:), t(:)
real              :: F, m, k, h, v0
integer           :: i, n, method

    ! == Calling the subroutine for reading data ==
    ! == Вызываем подпрограмму для чтения данных ==
    call ReadAndPrintData(F,m,k,v0,h,n,method)
    allocate(v(n+1))
    allocate(t(n+1))
    
    t(1) = 0
    do i = 2,n
        t(i) = (i-1)*h
    end do
    
    print*, '====================================='
    print*, ''
    
    ! == Selecting method of solution ==    
    ! == Выбор метода решения         ==
    if (method==1) then
        v = Euler(F,m,k,h,n)
    
    else if (method==2) then
        v = Runge_Kutta_2nd_Order(F,m,k,h,n)
    
    else if (method==3) then
        v = Runge_Kutta_4th_Order(F,m,k,h,n)
    
    else if (method==4) then
        v(1) = 0
        do i=1,n
            v(i) = sqrt(F/k) * tanh( abs(((1/2)*log(abs(((1+(sqrt(k/F)*v0)))/(1-(sqrt(k/F)*v0))))) - (sqrt(k*F)/m)*t(i) ))
        end do
        
    else
        print*, 'Error in selecting method'
        pause
    end if
    
    ! == Saving solution to a new file     ==
    ! == Сохраняем вычисление в новый файл ==
     open(2, file='OutputData.txt', status='replace')
        write(2,'(a12, a12)') 'v(t)', 't'
        
        do i=1,n
            write(2,'(f12.5,f12.5)') t(i)
        end do
        write(2,*) '==========================='
        do i=1,n
            write(2,'(f12.5,f12.5)') v(i)
        end do
        
     close(2)
        
    print*, 'File successfully created.'
    call PercentCounter(F,m,k,h,n,v)
    
    call OrderOfAccuracyEstimation()
    
    deallocate(v)
    deallocate(t)
pause    
end program main
    
