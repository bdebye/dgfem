module gmsh

implicit none

integer :: node_n
integer :: tetrah_n
integer :: triang_n
integer :: segment_n

type msh_node
    real(kind=8), dimension(3) :: coord
end type

type msh_segment
    integer :: physical_tag
    integer :: geometrical_tag
    integer, dimension(2) :: node_number_list
end type

type msh_triang
    integer :: physical_tag
    integer :: geometrical_tag
    integer, dimension(3) :: node_number_list
end type

type msh_tetrah
    integer :: physical_tag
    integer :: geometrical_tag
    integer, dimension(4) :: node_number_list
end type

type msh_face
    integer :: eid
    integer :: local_fid
end type

type quad_point
    real(kind=8) :: coord(3)
    real(kind=8) :: W
end type

type(msh_node), dimension(:), allocatable :: node
type(msh_segment), dimension(:), allocatable :: segment
type(msh_triang), dimension(:), allocatable :: triang
type(msh_tetrah), dimension(:), allocatable :: tetrah

type(msh_face), dimension(:), allocatable :: boundary

integer, parameter :: SEGMENT_TAG = 1
integer, parameter :: TRIANG_TAG = 2
integer, parameter :: TETRAH_TAG = 4

contains

subroutine count_element_number(filename)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=100) :: buffer
    integer :: count, i, element_n, tmp

    open(unit=10, file=filename)
    do i = 1, 4
        read (10, '(A)'), buffer
    end do
    read (10, *), count
    node_n = count
    do i = 1, count + 1
        read (10, '(A)'), buffer
    end do
    read (10, '(A)'), buffer
    read (10, *), element_n
    segment_n = 0
    triang_n = 0
    tetrah_n = 0
    
    do i = 1, element_n
        read (10, *), count, tmp
        if(tmp == SEGMENT_TAG) then
        segment_n = segment_n + 1
        else if(tmp == TRIANG_TAG) then
            triang_n = triang_n + 1
        else if(tmp == TETRAH_TAG) then
            tetrah_n = tetrah_n + 1
        end if
    end do
    close(10)
end subroutine

subroutine read_node(unit, i)
    integer, intent(in) :: unit, i
    integer :: n
    read (unit, *), n, node(i) % coord
!	print *, node(i) % coord
end subroutine

subroutine bubble_sort(a, n)
  implicit none
  integer :: n, a(n)
  integer i, j, temp
  do i = n - 1, 1, -1 
    do j = 1, i
      if (a(j) > a(j + 1)) then
        temp = a(j)
        a(j) = a(j + 1)
        a(j + 1) = temp
      end if
    end do
  end do
  return
end subroutine

subroutine read_element(unit, segment_i, triang_i, tetrah_i)
    integer, intent(in) :: unit
    integer, intent(inout) :: segment_i, triang_i, tetrah_i
    integer :: tmp, type
    integer, dimension(20) :: buffer
    character(len=50) :: line
    read (unit, '(A)'), line
!	print *, line
    read (line, *), tmp, type
    if(type == SEGMENT_TAG) then 
        read (line, *), buffer(1:7)
        segment_i = segment_i + 1
        segment(segment_i) % physical_tag = buffer(4)
        segment(segment_i) % geometrical_tag = buffer(5)
        segment(segment_i) % node_number_list = buffer(6:7)
        call bubble_sort(segment(segment_i) % node_number_list, 2)
!		print *, segment(segment_i) % node_number_list
    else if (type == TRIANG_TAG) then
        read (line, *), buffer(1:8)
        triang_i = triang_i + 1
        triang(triang_i) % physical_tag = buffer(4)
        triang(triang_i) % geometrical_tag = buffer(5)
        triang(triang_i) % node_number_list = buffer(6:8)
        call bubble_sort(triang(triang_i) % node_number_list, 3)
!		print *, triang(triang_i) % node_number_list
    else if (type == TETRAH_TAG) then
        read (line, *), buffer(1:9)
        tetrah_i = tetrah_i + 1
        tetrah(tetrah_i) % physical_tag = buffer(4)
        tetrah(tetrah_i) % geometrical_tag = buffer(5)
        tetrah(tetrah_i) % node_number_list = buffer(6:9)
        call bubble_sort(tetrah(tetrah_i) % node_number_list, 4)
!		print *, tetrah(tetrah_i) % node_number_list
    end if
end subroutine

subroutine jump_line(unit, n)
    character(len=100) :: buffer
    integer, intent(in) :: unit, n
    integer :: i
    do i = 1, n
        read (unit, '(A)'), buffer
    end do
end subroutine

subroutine search_boundary()
    integer :: i, j, k
    integer :: count = 0
    print *, "Searching mesh boundary..."
    do i = 1, tetrah_n
        do j = 1, triang_n
            do k = 1, 4
                if (all(pid_of_face(i, k) == triang(j) % node_number_list, 1)) then
                    count = count + 1
                end if
            end do
        end do
    end do
    allocate(boundary(count))
    print *, "Totally ", count, " boundary facets..."
    count = 1
    do i = 1, tetrah_n
        do j = 1, triang_n
            do k = 1, 4
                if (all(pid_of_face(i, k) == triang(j) % node_number_list, 1)) then
                    boundary(count) % eid = i
                    boundary(count) % local_fid = k
                    count = count + 1
                end if
            end do
        end do
    end do
end subroutine

subroutine load_mesh(filename)
    character(len=*), intent(in) :: filename
    integer :: i
    integer :: segment_i = 0, triang_i = 0, tetrah_i = 0
    integer :: element_n

    print *, "laoding mesh..."
    call clear_mesh()
    call count_element_number(filename)
    allocate(node(node_n))
    allocate(segment(segment_n))
    allocate(triang(triang_n))
    allocate(tetrah(tetrah_n))

    open(unit=10, file=filename)
    call jump_line(10, 5)
    do i = 1, node_n
        call read_node(10, i)
    end do
    call jump_line(10, 2)
    read (10, *) element_n
    do i = 1, element_n
        call read_element(10, segment_i, triang_i, tetrah_i)
    end do
    close(10)
    call search_boundary()
    print *, "Finished loading mesh..."
end subroutine

subroutine clear_mesh()
    node_n = 0
    segment_n = 0
    triang_n = 0
    tetrah_n = 0
    if(allocated(node)) then
        deallocate(node)
    end if
    if(allocated(segment)) then
        deallocate(segment)
    end if
    if(allocated(triang)) then
        deallocate(triang)
    end if
    if(allocated(tetrah)) then
        deallocate(tetrah)
    end if
    if(allocated(boundary)) then
        deallocate(boundary)
    end if
end subroutine

function pts(pid) result(p)
    integer :: pid
    real(kind=8) :: p(3)
    p = node(pid)%coord 
end function

function find_id(id, list_id) result(pos)
    integer :: i, pos, id
    integer :: list_id(:)
    pos = -1
    do i = 1, size(list_id)
        if (id == list_id(i)) then
            pos = i
            return
        end if
    end do
end function

function count_share(list_a, list_b) result(count)
    integer :: list_a(:), list_b(:)
    integer :: i, count
    count = 0
    do i = 1, size(list_a)
        if (find_id(list_a(i), list_b) > 0) then
            count = count + 1
        end if
    end do
end function

function pid_from_coord(coord) result(pid)
    integer i, pid
    real(kind=8) :: coord(3)
    pid = -1
    do i = 1, node_n
        if (all(abs(node(i) % coord - coord) < 1.0d-9, 1)) then
            pid = i
            return
        end if
    end do
end function

function pid_of_face(eid, local) result(pids)
    integer :: eid, local, pids(3)
    select case(local)
        case(1)
            pids = tetrah(eid) % node_number_list([2, 3, 4])
        case(2)
            pids = tetrah(eid) % node_number_list([1, 3, 4])
        case(3)
            pids = tetrah(eid) % node_number_list([1, 2, 4])
        case(4)
            pids = tetrah(eid) % node_number_list([1, 2, 3])
        case default
            pids = [-1, -1, -1]
    end select
end function

subroutine eid_from_pid(pid, eid, num)
    integer :: pid(:), eid(:)
    integer :: num
    integer :: count = 1, i
    do i = 1, tetrah_n
        if (count_share(pid, tetrah(i) % node_number_list) &
                == size(pid)) then
            if (count > size(eid)) then
                print *, "Error, buffer not enough..."
                exit
            end if
            eid(count) = i
            count = count + 1
        end if
    end do
    num = count - 1
end subroutine

function adjacent_element(eid, local) result(id)
    integer :: eid, local, id, i
    integer :: pids(3)
    id = -1
    pids = pid_of_face(eid, local)
    do i = 1, tetrah_n
        if (count_share(pids, tetrah(i) % node_number_list) == 3 &
            .and. (i .ne. eid)) then
            id = i
            return
        end if
    end do
end function

function normal_of_face(eid, local) result(n)
    use linalg
    integer :: eid, local, pid(3), i
    real(kind=8) :: p(3, 3), v(3)
    real(kind=8) :: n(3)
    pid = pid_of_face(eid, local)
    do i = 1, 3
        p(:, i) = pts(pid(i))
    end do
    v = pts(tetrah(eid) % node_number_list(local))
    n = cross(p(:, 2) - p(:, 1), p(:, 3) - p(:, 2))
    n = n / norm(n)
    if (dot(p(:, 1) - v, n) < 0) then
        n = - n
    end if
end

function center_of_element(eid) result(center)
    integer :: eid, pid(4)
    real(kind=8) :: center(3)
    pid = tetrah(eid) % node_number_list
    center = node(pid(1)) % coord &
           + node(pid(2)) % coord &
           + node(pid(3)) % coord &
           + node(pid(4)) % coord
    center = center / 4.0d0
end function

function center_of_face(eid, local) result(center)
    integer :: eid, local, pid(3)
    real(kind=8) :: center(3)
    pid = pid_of_face(eid, local)
    center = node(pid(1)) % coord &
           + node(pid(2)) % coord &
           + node(pid(3)) % coord
    center = center / 3d0
end function

function face_area(eid, local) result(area)
    use linalg
    integer :: pid(3), eid, local
    real(kind=8) :: area
    pid = pid_of_face(eid, local)
    area = norm(cross(pts(pid(2)) - pts(pid(1)), pts(pid(3)) - pts(pid(1)))) / 2.0d0
end function


function tetrah_quad_point(eid) result(qpoints)
    integer :: eid, pid(4)
    real(kind=8) :: r(3), s(3), t(3), base(3)
    type(quad_point) :: qpoints(5)
    
    pid = tetrah(eid) % node_number_list
    base = pts(pid(1))
    r = pts(pid(2)) - pts(pid(1))
    s = pts(pid(3)) - pts(pid(1))
    t = pts(pid(4)) - pts(pid(1))
    
    qpoints(1) % coord = base + r / 4.0d0 + s / 4.0d0 + t / 4.0d0
    qpoints(2) % coord = base + r / 2.0d0 + s / 6.0d0 + t / 6.0d0
    qpoints(3) % coord = base + r / 6.0d0 + s / 2.0d0 + t / 6.0d0
    qpoints(4) % coord = base + r / 6.0d0 + s / 6.0d0 + t / 2.0d0
    qpoints(5) % coord = base + r / 6.0d0 + s / 6.0d0 + t / 6.0d0
    
    qpoints(1) % W = - 4.0d0 / 30.0d0
    qpoints(2) % W = 9.0d0 / 120.0d0
    qpoints(3) % W = 9.0d0 / 120.0d0
    qpoints(4) % W = 9.0d0 / 120.0d0
    qpoints(5) % W = 9.0d0 / 120.0d0
end function

function triang_quad_point(eid, f_local) result(qpoints)
    integer :: eid, f_local, pid(3)
    real(kind=8) :: r(3), s(3), base(3)
    type(quad_point) :: qpoints(4)
    pid = pid_of_face(eid, f_local)
    
    base = pts(pid(1))
    r = pts(pid(2)) - pts(pid(1))
    s = pts(pid(3)) - pts(pid(1))
    
    qpoints(1) % coord = base + r / 3.0d0 + s / 3.0d0
    qpoints(2) % coord = base + 0.6d0 * r + s / 5.0d0
    qpoints(3) % coord = base + r / 5.0d0 + 0.6d0 * s
    qpoints(4) % coord = base + r / 5.0d0 + s / 5.0d0
    
    qpoints(1) % W = - 9.0d0 / 32.0d0
    qpoints(2) % W = 25.0d0 / 96.0d0
    qpoints(3) % W = 25.0d0 / 96.0d0
    qpoints(4) % W = 25.0d0 / 96.0d0
end function

end module
