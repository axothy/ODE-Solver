module methods
implicit none

    contains
    ! == A subroutine that reads and outputs data from file to console      ==
    ! == Ïîäïðîãðàììà, âûïîëíÿþùàÿ ÷òåíèå è âûâîä èñõîäíûõ äàííûõ â êîíñîëü ==
    subroutine ReadAndPrintData(F,m,k,v0,h,n,method)
    implicit none
    
        real       :: F, m, k, h, v0
        integer    :: n, method
        
        ! == Reading data from file         ==
        ! == ×òåíèå âõîäíûõ äàííûõ èç ôàéëà ==
        open(1, file='Data.txt')

        read(1,*) F, m, k, v0, h, n, method
        close(1)
        print*, 'File read successfully.'
        
        ! == Outputting read data to the console ==
        ! == Âûâîä ïðî÷èòàííûõ äàííûõ â êîíñîëü  == 
        write(*,'(a5, a8, a8, a10)') 'F', 'm', 'k', 'step'
        write(*,'(f8.2,f8.2,f8.2,f8.4)') F,m,k,h
        write(*,*) ''
    end subroutine ReadAndPrintData
    
    !== The function that performs the Euler's method  ==
    !== Ôóíêöèÿ, âûïîëíÿþùàÿ âû÷èñëåíèå ïî ìåòîäó Ýéëåðà ==
    function Euler(F,m,k,h,n) result (v)
        real, allocatable :: v(:)
        real              :: F,m,k,h
        integer           :: i, n
        
        allocate(v(n+1))
        
        ! == Solution via Euler's method ==
        ! == Ðåøåíèå ìåòîäîì Ýéëåðà      == 
        print*, 'Finding velocity via Eulers method...'
        v(1) = 0
        do i=2,n
            v(i) = v(i-1) + h*fxy(F,m,k,v(i-1))
        end do
    end function Euler
    
    !== The function that performs the Runge-Kutta 2nd Order method        ==
    !== Ôóíêöèÿ, âûïîëíÿþùàÿ âû÷èñëåíèå ïî ìåòîäó Ðóíãå-Êóòòû 2-ãî ïîðÿäêà ==
    function Runge_Kutta_2nd_Order(F,m,k,h,n) result(v)
        real, allocatable :: v(:)
        real              :: F,m,k,h
        real              :: k1, k2
        integer           :: i, n
        
        allocate(v(n+1))
        
        ! == Solution via Runge-Kutta 2nd Order method ==
        ! == Ðåøåíèå ìåòîäîì Ðóíãå-Êóòòû 2-ãî ïîðÿäêà  ==
        print*, 'Finding velocity via Runge-Kutta 2nd Order method...'
        v(1) = 0
        do i=2,n
            k1 = h*fxy(F,m,k,v(i-1))
            k2 = h*fxy(F,m,k,(v(i-1)+(k1/2)))
            v(i) = v(i-1) + k2
        end do

    end function Runge_Kutta_2nd_Order
    
    !== The function that performs the Runge-Kutta 4th Order method        ==
    !== Ôóíêöèÿ, âûïîëíÿþùàÿ âû÷èñëåíèå ïî ìåòîäó Ðóíãå-Êóòòû 4-ãî ïîðÿäêà ==
    function Runge_Kutta_4th_Order(F,m,k,h,n) result(v)
        real, allocatable :: v(:)
        real              :: F,m,k,h
        real              :: k1, k2, k3, k4
        integer           :: i, n
        
        allocate(v(n+1))
        
        ! == Solution via Runge-Kutta 4th Order method ==
        ! == Ðåøåíèå ìåòîäîì Ðóíãå-Êóòòû 4-ãî ïîðÿäêà  ==
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
    ! == Îïèñûâàåì íàøó ôóíêöèþ == 
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
    
program main
use methods
use percent

implicit none

real, allocatable :: v(:), t(:)
real              :: F, m, k, h, v0
integer           :: i, n, method

    ! == Calling the subroutine for reading data ==
    ! == Âûçûâàåì ïîäïðîãðàììó äëÿ ÷òåíèÿ äàííûõ ==
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
    ! == Âûáîð ìåòîäà ðåøåíèÿ         ==
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
    ! == Ñîõðàíÿåì âû÷èñëåíèå â íîâûé ôàéë ==
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
    
    deallocate(v)
    deallocate(t)
pause    
end program main
    
