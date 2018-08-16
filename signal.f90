module signal
    
contains

function pow(a, n)
    real(kind=8) :: pow, a
    integer :: n, i
    pow = 1.0d0
    do i = 1, n
        pow = pow * a
    end do
end function

function gauss_dd(t, tau) result(A)
    real(kind=8) :: t, A, tau
    A = (((pow(t, 2) - pow(tau, 2)) / pow(tau, 4)) * exp(-pow(t, 2) &
            / pow(tau, 2) / 2.0d0)) / ((- pow(tau, 2)) / pow(tau, 4))
end function

function neumann_pulse(t, tau) result(A)
    real(kind=8) :: t, tau, A
    A = - (t / tau) * exp(-0.5d0 * pow(t / tau, 2))
end function

function gauss_module(t, tau, omega) result(A)
    real(kind=8) :: t, tau, omega
    A = exp(-0.5d0 * pow(t / tau, 2)) * sin(omega * t)
end function

function gauss_differential(t, tau) result(A)
    real(kind=8) :: t, tau, A
    A = - (t / pow(tau, 2)) * exp(-pow(t, 2) / pow(tau, 2) / 2.0d0)
end function

function gauss_module_diff(t, tau, omega) result(A)
    real(kind=8) :: t, tau, omega
    A = omega * cos(omega * t) *exp(-0.5d0 * pow(t / tau, 2)) + sin(omega * t) * gauss_differential(t, tau)
end function
    
end module