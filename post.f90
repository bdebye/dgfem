module post

implicit none

type scalar_field
    real, dimension(3) :: coord
    real :: value
end type

type vector_field
    real, dimension(3) :: coord
    real, dimension(3) :: value
end type

character(len=20), parameter :: POS_HEADER = 'View "field" {'
character(len=2), parameter :: POS_FOOTER = '};'

contains

subroutine save_scalar_field(filename, field)
    character(len=*), intent(in) :: filename
    type(scalar_field), dimension(:), intent(in) :: field
    integer :: i

    open(unit=10, file=filename)
    write (10, '(A)'), POS_HEADER
    do i = 1, size(field)
        write (10, 100), field(i) % coord, field(i) % value
    end do
    write (10, '(A)'), POS_FOOTER
    close(10)

    100 format ('SP(', F, ',', F, ',', F, '){', E, '};')
end subroutine

subroutine save_vector_field(filename, field)
    character(len=*), intent(in) :: filename
    type(vector_field), dimension(:), intent(in) :: field
    integer :: i

    open(unit=10, file=filename)
    write (10, '(A)'), POS_HEADER
    do i = 1, size(field)
        write (10, 100), field(i) % coord, field(i) % value
    end do
    write (10, '(A)'), POS_FOOTER
    close(10)

    100 format ('VP(', F, ',', F, ',', F, '){', E, ',', E, ',', E'};')
end subroutine

end module
