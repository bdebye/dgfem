module dgfem
    
use linalg
implicit none

integer :: Ne
integer, parameter :: icomb(2, 6) = [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]
real(kind=8), parameter :: M(4, 4) = [[0.1d0, 0.5d-1, 0.5d-1, 0.5d-1], &
                                      [0.5d-1, 0.1d0, 0.5d-1, 0.5d-1], &
                                      [0.5d-1, 0.5d-1, 0.1d0, 0.5d-1], &
                                      [0.5d-1, 0.5d-1, 0.5d-1, 0.1d0]]

real(kind=8), parameter :: pi = 4.0d0 * atan(1.0d0)
real(kind=8), parameter :: eps0 = 8.854d-012
real(kind=8), parameter :: sig0 = 0.0d0
real(kind=8), parameter :: mu0 = 1.257d-006
real(kind=8), parameter :: c0 = 1.0d0 / sqrt(eps0 * mu0)
real(kind=8), parameter :: Zs = sqrt(mu0 / eps0)
real(kind=8), parameter :: Ys = sqrt(eps0 / mu0)
real(kind=8), parameter :: tau = 1.0d-9
real(kind=8), parameter :: dt = 1d-11
real(kind=8), parameter :: f = 2.0d8
real(kind=8), parameter :: omega = 2.0d0 * pi * f

integer, parameter :: N = 3000
integer :: record_id
real(kind=8), parameter :: record_point(3) = [1.0d0, 2.0d0, 2.0d0]
real(kind=8) :: record_data(N)

integer :: source_elem(20), source_elem_number
integer :: sid1, sid2

type elem
    integer :: eid
    integer :: pid(4)
    integer :: exterior_id(4)
    
    real(kind=8) :: epsilon, mu, sigma
    real(kind=8) :: volumn
    
    real(kind=8) :: He(6)
    real(kind=8) :: Ee(6)
    real(kind=8) :: fee(6)
    real(kind=8) :: fhe(6)
    real(kind=8) :: bee(6)
    real(kind=8) :: bhe(6)
    
    real(kind=8) :: l(6)
    real(kind=8) :: f_area(4)
    
    real(kind=8) :: coord_matrix(4, 4)
    real(kind=8) :: coef_matrix(4, 4)
    
    real(kind=8) :: Tee(6, 6)
    real(kind=8) :: The(6, 6)
    real(kind=8) :: Se(6, 6)
    real(kind=8) :: Ree(6, 6)
    real(kind=8) :: Rhe(6, 6)
    
contains

    procedure :: u => element_u
    procedure :: p => element_p
    procedure :: F => element_F
    procedure :: N => base_value
    procedure :: curl_N => curl_base
    
    procedure :: E_field => element_E_field
    procedure :: H_field => element_H_field
    procedure :: center_coord => element_center_coord
    
    procedure :: center_E_field => element_center_E_field
    procedure :: center_H_field => element_center_H_field
    procedure :: center_Ez => element_center_Ez
    
    procedure :: Se_inte_func => element_Se_inte_func
    
end type

type(elem), allocatable, target :: elements(:)
    
contains

function simplex_id(e, i) result(si)
    integer :: e, i, si
    si = icomb(i, e)
end function

subroutine init_elem(eid)
    use gmsh
    use linalg
    integer :: eid, i, li1, li2, e
    type(elem), pointer :: el
    
    el => elements(eid)
    el % eid = eid
    el % pid = tetrah(eid) % node_number_list
    el % epsilon = eps0
    el % mu = mu0
    el % sigma = 0.0d0
    do i = 1, 4
        el % coord_matrix(1: 3, i) = pts(el % pid(i))
        el % f_area(i) = face_area(eid, i)
        el % exterior_id(i) = adjacent_element(eid, i)
    end do
    el % coord_matrix(4, :) = 1.0d0
    el % volumn = abs(determin(el % coord_matrix)) / 6.0d0
    el % coef_matrix = el % coord_matrix
    call inverse(el % coef_matrix)
    do e = 1, 6
        li1 = simplex_id(e, 1)
        li2 = simplex_id(e, 2)
        el % l(e) = norm(pts(el % pid(li1)) - pts(el % pid(li2)))
    end do
    
    el % bee = 0.0d0
    el % bhe = 0.0d0
    
    call fill_Tee(el)
    call fill_The(el)
    call fill_Se(el)
    call fill_Ree(el)
    call fill_Rhe(el)
end subroutine

subroutine fill_Tee(el)
    type(elem) :: el
    integer :: i, j
    do i = 1, 6
        do j = 1, 6
            el % Tee(i, j) = el % epsilon * el % F(i, j)
        end do
    end do
end subroutine

subroutine fill_The(el)
    type(elem) :: el
    integer :: i, j
    do i = 1, 6
        do j = 1, 6
            el % The(i, j) = el % mu * el % F(i, j)
        end do
    end do
end subroutine

subroutine fill_Ree(el)
    type(elem) :: el
    integer :: i, j
    if (el % sigma == 0.0d0) then
        el % Ree = 0.0d0
    else
        do i = 1, 6
            do j = 1, 6
                el % The(i, j) = el % sigma * el % F(i, j)
            end do
        end do
    end if
end subroutine

subroutine fill_Rhe(el)
    type(elem) :: el
    el % Rhe = 0.0d0
end subroutine

subroutine update_fee(el)
    use gmsh
    type(elem) :: el
    type(quad_point) :: im_data(4)
    real(kind=8) :: n(3), p(3), area, weight
    integer :: i, j, k
    el % fee = 0.0d0
    do j = 1, 4
        if (el % exterior_id(j) == -1) then
            cycle
        else
            im_data = triang_quad_point(el % eid, j)
            n = normal_of_face(el % eid, j)
            area = el % f_area(j)
            do k = 1, 4
                p = im_data(k) % coord
                weight = im_data(k) % W
                do i = 1, 6
                    el % fee(i) = el % fee(i) + area * weight * dot(el % N(i, p), &
                            cross(n, elements(el % exterior_id(j)) % H_field(p)))
                end do
            end do
        end if
    end do
end subroutine

subroutine update_fhe(el)
    use gmsh
    type(elem) :: el
    type(quad_point) :: im_data(4)
    real(kind=8) :: n(3), p(3), area, weight
    integer :: i, j, k
    el % fhe = 0.0d0
    do j = 1, 4
        if (el % exterior_id(j) == -1) then
            cycle
        else
            im_data = triang_quad_point(el % eid, j)
            n = normal_of_face(el % eid, j)
            area = el % f_area(j)
            do k = 1, 4
                p = im_data(k) % coord
                weight = im_data(k) % W
                do i = 1, 6
                    el % fhe(i) = el % fhe(i) + area * weight * dot(cross(n, el % N(i, p)), &
                            elements(el % exterior_id(j)).E_field(p))
                end do
            end do
        end if
    end do
end subroutine

subroutine update_Ee(el)
    use lapack95
    use util
    type(elem) :: el
    real(kind=8) :: A(6, 6)
    integer :: ipiv(6)
    A = el % Tee / dt + 0.5d0 * el % Ree
    call update_fee(el)
    el % Ee = matmul(el % Tee / dt - 0.5d0 * el % Ree, el % Ee) &
                + matmul(el % Se, el % He) + el % fee + el % bee

    call getrf(A, ipiv)
    call getrs(A, ipiv, el % Ee)
end subroutine

subroutine update_He(el)
    use lapack95
    type(elem) :: el
    real(kind=8) :: A(6, 6)
    integer :: ipiv(6)
    A = el % The / dt + 0.5d0 * el % Rhe
    call update_fhe(el)
    el % He = matmul(el % The / dt - 0.5d0 * el % Rhe, el % He) &
                - matmul(el % Se, el % Ee) + el % fhe + el % bhe
    call getrf(A, ipiv)
    call getrs(A, ipiv, el % He)
end subroutine

subroutine impress_source(el, pid1, pid2, dI)
    type(elem) :: el
    integer :: pid1, pid2
    real(kind=8) :: dI
    
    if (all([pid1, pid2] == [el % pid(1), el % pid(2)], 1)) then
        el % bee(1) = - el % l(1) * dI
    else if (all([pid1, pid2] == [el % pid(2), el % pid(1)], 1)) then
        el % bee(1) = el % l(1) * dI
    else if (all([pid1, pid2] == [el % pid(1), el % pid(3)], 1)) then
        el % bee(2) = - el % l(2) * dI
    else if (all([pid1, pid2] == [el % pid(3), el % pid(1)], 1)) then
        el % bee(2) = el % l(2) * dI
    else if (all([pid1, pid2] == [el % pid(1), el % pid(4)], 1)) then
        el % bee(3) = - el % l(3) * dI
    else if (all([pid1, pid2] == [el % pid(4), el % pid(1)], 1)) then
        el % bee(3) = el % l(3) * dI
    else if (all([pid1, pid2] == [el % pid(2), el % pid(3)], 1)) then
        el % bee(4) = - el % l(4) * dI
    else if (all([pid1, pid2] == [el % pid(3), el % pid(2)], 1)) then
        el % bee(4) = el % l(4) * dI
    else if (all([pid1, pid2] == [el % pid(2), el % pid(4)], 1)) then
        el % bee(5) = - el % l(5) * dI
    else if (all([pid1, pid2] == [el % pid(4), el % pid(2)], 1)) then
        el % bee(5) = el % l(5) * dI
    else if (all([pid1, pid2] == [el % pid(3), el % pid(4)], 1)) then
        el % bee(6) = - el % l(6) * dI
    else if (all([pid1, pid2] == [el % pid(4), el % pid(3)], 1)) then
        el % bee(6) = el % l(6) * dI
    end if
end subroutine

#define a_(i) (el.coef_matrix(i, 4))
#define b_(i) (el.coef_matrix(i, 1))
#define c_(i) (el.coef_matrix(i, 2))
#define d_(i) (el.coef_matrix(i, 3))

function element_u(el, i, j)
    real(kind=8) :: element_u(3)
    class(elem) :: el
    integer :: i, j
    element_u = [c_(i) * d_(j) - c_(j) * d_(i), &
         b_(j) * d_(i) - b_(i) * d_(j), &
         b_(i) * c_(j) - b_(j) * c_(i)]
end function

function element_p(el, i, j)
    real(kind=8) :: element_p
    class(elem) :: el
    integer :: i, j
    element_p = b_(i) * b_(j) + c_(i) * c_(j) + d_(i) * d_(j)
end function

function element_F(el, i, j)
    use util
    class(elem) :: el
    integer :: i, j, i1, i2, j1, j2
    real(kind=8) :: element_F
    i1 = simplex_id(i, 1)
    i2 = simplex_id(i, 2)
    j1 = simplex_id(j, 1)
    j2 = simplex_id(j, 2)
    element_F =  el % volumn * el % l(i) * el % l(j) &
                    * (el % p(i2, j2) * M(i1, j1) &
                     - el % p(i2, j1) * M(i1, j2) &
                     - el % p(i1, j2) * M(i2, j1) &
                     + el % p(i1, j1) * M(i2, j2) &
                )
end function

function base_value(el, i, p) result(vector)
    class(elem) :: el
    integer :: i, i1, i2
    real(kind=8) :: p(3), vector(3), Li1, Li2
    i1 = simplex_id(i, 1)
    i2 = simplex_id(i, 2)
    Li1 = a_(i1) + b_(i1) * p(1) + c_(i1) * p(2) + d_(i1) * p(3)
    Li2 = a_(i2) + b_(i2) * p(1) + c_(i2) * p(2) + d_(i2) * p(3)
    vector = el % l(i) * ( &
            Li1 * [b_(i2), c_(i2), d_(i2)] - &
            Li2 * [b_(i1), c_(i1), d_(i1)]   &
        )
end function

function curl_base(el, i, p) result(vector)
    class(elem) :: el
    integer :: i, i1, i2
    real(kind=8) :: p(3), vector(3)
    i1 = simplex_id(i, 1)
    i2 = simplex_id(i, 2)
    vector = 2.0d0 * el % l(i) * el % u(i1, i2)
end function

function element_E_field(el, p) result(Ev)
    class(elem) :: el
    real(kind=8) :: p(3), Ev(3)
    Ev = el % Ee(1) * el % N(1, p) + &
         el % Ee(2) * el % N(2, p) + &
         el % Ee(3) * el % N(3, p) + &
         el % Ee(4) * el % N(4, p) + &
         el % Ee(5) * el % N(5, p) + &
         el % Ee(6) * el % N(6, p)
end function

function element_H_field(el, p) result(Hv)
    class(elem) :: el
    real(kind=8) :: p(3), Hv(3)
    Hv = el % He(1) * el % N(1, p) + &
         el % He(2) * el % N(2, p) + &
         el % He(3) * el % N(3, p) + &
         el % He(4) * el % N(4, p) + &
         el % He(5) * el % N(5, p) + &
         el % He(6) * el % N(6, p)
end function

function element_center_coord(el) result(coord)
    use gmsh
    class(elem) :: el
    real(kind=8) :: coord(3)
    coord = center_of_element(el % eid)
end function

function element_center_E_field(el) result(Ev)
    class(elem) :: el
    real(kind=8) :: Ev(3)
    Ev = el % E_field(el % center_coord())
end function

function element_center_H_field(el) result(Hv)
    class(elem) :: el
    real(kind=8) :: Hv(3)
    Hv = el % H_field(el % center_coord())
end function

function element_center_Ez(el) result(Ez)
    class(elem) :: el
    real(kind=8) :: Ez
    Ez = dot(el % center_E_field(), [0.0d0, 0.0d0, 1.0d0])
end function

function element_Se_inte_func(el, i, j, p) result(d)
    use linalg
    class(elem) :: el
    integer :: i, j
    real(kind=8) :: p(3), d
    d = 0.5d0 * (dot(el % N(i, p), el % curl_N(j, p)) + &
                 dot(el % curl_N(i, p), el % N(j, p)))
end function

subroutine fill_Se(el)
    use gmsh
    type(elem) :: el
    type(quad_point) :: im_data(5)
    integer :: i, j, k
    im_data = tetrah_quad_point(el % eid)
    do i = 1, 6
        do j = 1, 6
            do k = 1, 5
                el % Se(i, j) = el % Se(i, j) + 6.0d0 * el % volumn &
                    * im_data(k) % W &
                    * el % Se_inte_func(i, j, im_data(k) % coord)
            end do
        end do
    end do
end subroutine

subroutine impose_boundary(el, f_local)
    use gmsh
    type(elem) :: el
    type(quad_point) :: im_data(4)
    real(kind=8) :: n(3), ac, area
    real(kind=8) :: p(3), weight
    integer :: i, j, k, f_local
    im_data = triang_quad_point(el % eid, f_local)
    area = el % f_area(f_local)
    n = normal_of_face(el % eid, f_local)
    do k = 1, 4
        weight = im_data(k) % W
        p = im_data(k) % coord
        do i = 1, 6
            do j = 1, 6
                ac = area * weight * dot(cross(n, el % N(i, p)), &
                        cross(n, el % N(j, p)))
                el % Ree(i, j) = el % Ree(i, j) + Ys * ac
                el % Rhe(i, j) = el % Ree(i, j) + Zs * ac
            end do
        end do
    end do
end subroutine

subroutine configure_source()
    use gmsh
    real(kind=8) :: p1(3), p2(3)
    p1 = [2.0d0, 2.0d0, 1.9d0]
    p2 = [2.0d0, 2.0d0, 2.1d0]
    sid1 = pid_from_coord(p1)
    sid2 = pid_from_coord(p2)
    call eid_from_pid([sid1, sid2], source_elem, source_elem_number)
    print *, "Configured current source..."
end subroutine

subroutine impress_current_source(t)
    use signal
    real(kind=8) :: t, dI
    integer :: i
    dI = gauss_module(t - 5.0d0 * tau, tau, omega)
    do i = 1, source_elem_number
        call impress_source(elements(source_elem(i)), sid1, sid2, dI)
    end do
end subroutine

subroutine impose_boundary_condition()
    use gmsh
    integer :: i
    !$omp parallel do
    do i = 1, size(boundary)
        call impose_boundary(elements(boundary(i) % eid), &
                boundary(i) % local_fid)
    end do
    !$omp end parallel do
end subroutine

subroutine init_dgfem()
    use gmsh
    use util
    integer :: eid
    print *, "Initializing elements..."
    Ne = tetrah_n
    allocate(elements(Ne))
    !$omp parallel do
    do eid = 1, tetrah_n
        call init_elem(eid)
    end do
    !$omp end parallel do
    call impose_boundary_condition
end subroutine

subroutine update_E()
    integer :: i
    !$omp parallel do
    do i = 1, Ne
        call update_Ee(elements(i))
    end do
    !$omp end parallel do
end subroutine

subroutine update_H()
    integer :: i
    !$omp parallel do
    do i = 1, Ne
        call update_He(elements(i))
    end do
    !$omp end parallel do
end subroutine

subroutine show_element_info(eid)
    use util
    integer :: eid
    type(elem) :: el
    el = elements(eid)
    print *, "Tee"
    call show_matrix(el % Tee)
    print *, "The"
    call show_matrix(el % The)
    print *, "Se"
    call show_matrix(el % Se)
    print *, "Ree"
    call show_matrix(el % Ree)
    print *, "Rhe"
    call show_matrix(el % Rhe)
    print *, "fee"
    print *, el % fee
    print *, "fhe"
    print *, el % fhe
    print *, "bee"
    print *, el % bee
    print *, "bhe"
    print *, el % bhe
end subroutine


function time_step_requirement(el) result(time_step)
    type(elem) :: el
    real(kind=8) :: time_step
    time_step = 4.0d0 * el % volumn / (c0 * sum(el % f_area)) / (4.0d0 * sqrt(5.0d0) / 3.0d0 + 8.0d0 / 3.0d0)
end function

function maximum_time_step() result(time_step)
    real(kind=8) :: time_step
    integer :: i
    time_step = time_step_requirement(elements(1))
    do i = 1, Ne
        if (time_step_requirement(elements(i)) < time_step) then
            time_step = time_step_requirement(elements(i))
            !print *, i, time_step_requirement(elements(i))
        end if
    end do

end function

subroutine test_element()
    use util
    use gmsh
    integer :: i
    
    do i = 1, Ne
        if (norm(elements(i) % center_E_field()) > 100d0) then
            print *, i
        end if
    end do
end subroutine

subroutine export_electric_field(id)
    character(len=100) :: filename
    character(len=20), parameter :: POS_HEADER = 'View "field" {'
    character(len=2), parameter :: POS_FOOTER = '};'
    real(kind=8) :: coord(3), field(3)
    integer :: i, id

    print *, "Export electric field..."
    write (filename, '(A, I, A)') "view", id, ".pos"
    
    open(unit=10, file=filename)
    write (10, '(A)'), POS_HEADER
    do i = 1, Ne
        coord = elements(i) % center_coord()
        field = elements(i) % center_E_field()
        write (10, 100), coord, field
    end do
    write (10, '(A)'), POS_FOOTER
    close(10)

    100 format ('VP(', F, ',', F, ',', F, '){', E, ',', E, ',', E'};')
end subroutine

subroutine export_Ez_field(id)
    use gmsh
    type(elem) :: el
    character(len=100) :: filename
    character(len=20), parameter :: POS_HEADER = 'View "field" {'
    character(len=2), parameter :: POS_FOOTER = '};'
    real(kind=8) :: p1(3), p2(3), p3(3), p4(3),&
                    v1, v2, v3, v4
    integer :: i, id

    print *, "Export electric field..."
    write (filename, '(A, I, A)') "view", id, ".pos"
    
    open(unit=10, file=filename)
    write (10, '(A)'), POS_HEADER
    do i = 1, Ne
        el = elements(i)
        p1 = pts(el % pid(1))
        p2 = pts(el % pid(2))
        p3 = pts(el % pid(3))
        p4 = pts(el % pid(4))
        v1 = dot(el % E_field(p1), [0.0d0, 0.0d0, 1.0d0])
        v2 = dot(el % E_field(p2), [0.0d0, 0.0d0, 1.0d0])
        v3 = dot(el % E_field(p3), [0.0d0, 0.0d0, 1.0d0])
        v4 = dot(el % E_field(p4), [0.0d0, 0.0d0, 1.0d0])
        write (10, 100), p1, p2, p3, p4, v1, v2, v3, v4
    end do
    write (10, '(A)'), POS_FOOTER
    close(10)

    100 format ('SS(', 11(F, ','), F, '){', 3(E, ','), E, '};')
end subroutine

subroutine record_field(i)
    integer :: i
    real(kind=8) :: ef(3)
    ef = elements(record_id) % E_field(record_point)
    record_data(i) = ef(3)
end subroutine

subroutine output_record()
    integer :: i
    real(kind=8) :: t
    open(unit=10, file="signal.txt")
    do i = 1, N
        t = dt * i
        write (10, *), t, record_data(i)
    end do
    close(10)
end subroutine

function point_in_element(p, el) result(test)
    real(kind=8) :: p(3), c(4)
    type(elem) :: el
    integer :: i
    logical :: test
    do i = 1, 4
        c(i) = a_(i) * p(1) + b_(i) * p(2) + c_(i) * p(3) + d_(i)
    end do
    test = all(c >= 0.0d0 .and. c <= 1.0d0, 1)
end function

function element_contains_point(p) result(eid)
    integer :: i, eid
    real(kind=8) :: p(3)
    eid = -1
    do i = 1, Ne
        if (point_in_element(p, elements(i))) then
            eid = i
            return
        end if
    end do
    print *, "Error, failed to find the corresponding element of record point..."
end function

subroutine configure_signal_record()
    record_id = element_contains_point(record_point)
    print *, "The signal at the point (", record_point, ") will be recorded..."
end subroutine
    
end module

program main
    use gmsh
    use dgfem
    use signal
    real(kind=8) :: t
    call load_mesh("fetd_dipole.msh")
    call init_dgfem
    call configure_source
    call configure_signal_record
    print *, "Maximum time step width: ", maximum_time_step()
    print *, "-----------------------------------------------------------------"
    print *, "Start updating electromagnetic field..."
    print *
    do i = 1, N
        t = i * dt
        print *, i, t
        call impress_current_source(t)
        call update_E
        call update_H
        call record_field(i)
        if (mod(i, 100) == 0) then
            call export_Ez_field(i)
        end if
    end do
    call output_record
    call export_electric_field(-1)
    pause
end program