
program Seepage_Analysis
    implicit none
    integer :: n,i,j,k, X, Y, GRD_PTS, REM, HALF, REM1, Y1
    real :: HEAD
    real, dimension (:,:), allocatable :: H1, H2,H3
    
    interface
        subroutine Matrix_Trans(In1, In2)
        implicit none
        integer :: nrow, ncol, i, j
        real, dimension (:,:), allocatable :: In1, In2
        end subroutine Matrix_Trans      
    
        subroutine Matrix_Multi(In1,In2,Output)
        implicit none
        integer :: Row1, Col, Col2
        real, dimension (:,:), allocatable :: In1, In2, Output
        end subroutine Matrix_Multi
    
        subroutine MatrixInverse(Input1, Output1)
        implicit none    
        integer :: nrow, ncol
            real, dimension(:,:), allocatable :: Input1, Output1
        end subroutine MatrixInverse    
        
        subroutine PrintF(Input1)
        implicit none
        integer :: nrow, ncol
        real, dimension(:,:), allocatable :: Input1
        end subroutine PrintF
    end interface
    
    open(1,file = "E:\Fortran\Exercies\Solution to Navier-Stoke equation\test.txt")
    open(2,file = "E:\Fortran\Exercies\Solution to Navier-Stoke equation\RHS.txt")
    
    print *, "Enter number of X grids (Put odd numbers)"
    read (*,*) X
    print *, "(Enter number of Y grids)-1"
    read (*,*) Y
    print *, "Enter value of head"
    read (*,*) HEAD
    GRD_PTS = Y*(X-1)
    HALF = GRD_PTS/2
    print *, GRD_PTS
    allocate (H1(GRD_PTS, GRD_PTS), H2(GRD_PTS,1), H3(GRD_PTS,1))
    H1 = 0.0
!Writing the left half
    !Writing the first layer
    do i = 1, HALF,Y 
        REM = int((HALF-1)/Y + 1)
        H1(i,i) = -4
        H1(i,i+1) = 2
        if (i ==1) then !leftmost boundary or origin in this case
            H1(i,i+Y) = 2
        elseif (i == 1+Y*(REM-1)) then
            H1(i,i-Y) = 1
        else !in between left and right walls
            H1(i,i-Y) = 1
            H1(i,i+Y) = 1
        end if
    end do
    do k = 2,Y
    if (k == Y) then !Writing the top layer
        do i = k, HALF, Y
        REM = int((HALF-1)/Y + 1) !To find the rightmost cells
            H1(i,i) = -4
            H1(i,i-1) = 1
            if (i==k) then !To find the  leftmost cells
                H1(i,i+Y) = 2 
            elseif (i == k + Y*(REM-1)) then
                H1(i,i-Y) = 1
            else
                H1(i,i+Y) = 1
                H1(i,i-Y) = 1
            end if
         end do
        else !Writing the inner sandwich
        do i = k, HALF, Y
        REM = int((HALF-1)/Y + 1)
            H1(i,i) = -4
            H1(i,i-1) = 1
            H1(i,i+1) = 1
            if (i==k) then
                H1(i,i+Y) = 2
            elseif (i == k + Y*(REM-1)) then
                H1(i,i-Y) = 1
            else
                H1(i,i+Y) = 1
                H1(i,i-Y) = 1
            end if
        end do
        end if 
    end do
!Writing First Half ends
!Writing Second Half
    !Writing the first layer
    do i = HALF+1,GRD_PTS,Y 
        REM1 = int((HALF-1)/Y + 1)
        print *, REM1
        H1(i,i) = -4
        H1(i,i+1) = 2
        if (i == HALF+1) then !leftmost boundary or origin in this case
            H1(i,i+Y) = 1
        elseif (i == HALF+1+Y*(REM1-1)) then !rightmost boundary
            H1(i,i-Y) = 2
        else !in between left and right walls
            H1(i,i-Y) = 1
            H1(i,i+Y) = 1
        end if
    end do
    
    Y1 = HALF+Y
    do k = 2+HALF,Y1
    if (k == Y1) then !Writing the top layer
        do i = k, GRD_PTS, Y
        REM1 = int((HALF-1)/Y + 1)
            H1(i,i) = -4
            H1(i,i-1) = 1
            if (i==k) then
                H1(i,i+Y) = 1
            elseif (i == k + Y*(REM-1)) then
                H1(i,i-Y) = 2
            else
                H1(i,i+Y) = 1
                H1(i,i-Y) = 1
            end if
         end do
        else
        do i = k, GRD_PTS, Y
        REM = int((HALF-1)/Y + 1)
            H1(i,i) = -4
            H1(i,i-1) = 1
            H1(i,i+1) = 1
            if (i==k) then
                H1(i,i+Y) = 1
            elseif (i == k + Y*(REM-1)) then
                H1(i,i-Y) = 2
            else
                H1(i,i+Y) = 1
                H1(i,i-Y) = 1
            end if
        end do
        end if 
    end do
    !END WRITING THE H1 MATRIX
    !START WRITING THE RHS MATRIX I.E. H3
    H3 = 0.0
    !Considering the central zone
    do i = HALF-Y, HALF+Y
        if (i == HALF) then
            H3(i,1) = -HEAD*1.5
        else
            H3(i,1) = -HEAD/2
        end if
    end do
    do i = Y, HALF-Y, Y
        H3(i,1) = -HEAD
    end do
    !End writing RHS matrix
    !Finding the required values
    
    
   
    do i = 1, GRD_PTS
        do j = 1, GRD_PTS
    write(1,'(f6.2)',advance = 'no'), H1(i,j)
        END DO
        write(1,*)
    end do

    do i = 1, GRD_PTS
    write(2,'(f10.4)',advance = 'no'), H3(i,1)
    write(2,*)
    end do   
    
    
    call PrintF(H1)
            

    end program Seepage_Analysis
    
subroutine Matrix_Trans(In1, In2)
    implicit none
    integer :: nrow, ncol, i, j
    real, dimension (:,:), allocatable :: In1, In2
    nrow = size(In1, dim = 1)
    ncol = size(In1, dim = 2)
    allocate (In2(ncol,nrow))

    do j = 1, ncol
        do i = 1,nrow
            In2(j,i) = In1(i,j)
        end do
    end do
    end subroutine Matrix_Trans

subroutine Matrix_Multi(In1,In2,Output)
    implicit none
    integer :: Row1, Col, Col2, i, j, k
    real, dimension (:,:), allocatable :: In1, In2, Output
    Row1 = size(In1, dim = 1)
    Col = size(In1, dim = 2)
    Col2 = size(In2, dim = 2)
    allocate (Output(Row1,Col2))

    do i = 1,Row1
        do j = 1,Col2
            Output(i,j) = 0.0
            do k = 1,Col
                Output(i,j) = Output(i,j) + In1(i,k)*In2(k,j)
            end do
        end do
    end do
    end subroutine Matrix_Multi
    
    
    
subroutine MatrixInverse(Input1, Output1)
    implicit none
    integer :: nrow, ncol
    real, dimension(:,:), allocatable :: Input1, Output1
    integer :: i, j, k
    real :: X11, Xik
    nrow = size(Input1, dim = 1)
    ncol = nrow
    allocate(Output1(nrow, ncol))
    
    !writing an Identity matrix
    do i = 1, nrow
        do j = 1, ncol
            if (i == j) then
                Output1(i,j) = 1
            else
                Output1(i,j) = 0
            end if
        end do
    end do
    
    do k = 1, nrow
        X11 = Input1(k,k)
        do j = 1, ncol
            Input1(k,j) = Input1(k,j)/X11
            Output1(k,j) = Output1(k,j)/X11
        end do
          
        do i = 1, nrow
            if (i /= k) then
                Xik = Input1(i,k)
                do j = 1, ncol
                  Input1(i,j) = Input1(i,j) - Xik * Input1(k,j)
                  Output1(i,j) = Output1(i,j) - Xik * Output1(k,j)
                end do
            end if
        end do
    end do
end subroutine MatrixInverse 

subroutine PrintF(Input1)
    implicit none
    integer :: nrow, ncol
    real, dimension(:,:), allocatable :: Input1
    integer :: i, j
    nrow = size(Input1, dim = 1)
    ncol =size(Input1, dim = 2)
    do i = 1, nrow
        do j = 1, ncol
                write (*, '(F15.3)', advance='no'), Input1(i,j)
            end do
            write(*,*)
        end do
    end subroutine PrintF
