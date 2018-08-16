
module linalg

implicit none

contains

function determin(a)
    use lapack95
    real(kind=8) :: a(:, :), determin
    real(kind=8), allocatable :: b(:, :)
    integer :: sp(2)
    integer :: i

    sp = shape(a)
    allocate(b(sp(1), sp(2)))
    b = a
    call getrf(b)
    determin = 1.0
    do i = 1, sp(1)
        determin = determin * b(i, i)
    end do
    deallocate(b)
end function

subroutine inverse(a)
    use lapack95
    real(kind=8) :: a(:, :)
    integer, allocatable :: ipiv(:)
    integer :: sp(2)

    sp = shape(a)
    allocate(ipiv(sp(1)))

    call getrf(a, ipiv)
    call getri(a, ipiv)
    deallocate(ipiv)
end subroutine

function cross(a, b) result(c)
    real(kind=8) :: a(3), b(3), c(3)
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
end function

function dot(a, b) result(d)
    real(kind=8) :: a(3), b(3), d
    d = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)
end function

function norm(a) result(d)
    real(kind=8) a(3), d
    d = sqrt(a(1) * a(1) + a(2) * a(2) + a(3) * a(3))
end function

function identity(n) result(Id)
    integer :: n, i
    real(kind=8) :: Id(n, n)
    do i = 1, n
        Id(i, i) = 1.0d0
    end do
end function

function zeros(m, n) result(z)
    integer :: m, n
    real(kind=8) :: z(m, n)
    z = 0.0d0
end function

function ones(m, n) result(z)
    integer :: m, n
    real(kind=8) :: z(m, n)
    z = 1.0d0
end function

end module
