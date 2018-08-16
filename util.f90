module util
    
contains

subroutine show_matrix(mtx)
    real(kind=8) :: mtx(:, :)
    integer :: sh(2), i, j
    sh = shape(mtx)
    do i = 1, sh(1)
        do j = 1, sh(1)
            write (*, '(E)', advance="no") mtx(i, j)
        end do
        write (*, *)
    end do
    write (*, *)
end subroutine
    
end module