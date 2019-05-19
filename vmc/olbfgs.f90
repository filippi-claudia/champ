module olbfgs
implicit none

    private 
    
        integer, parameter :: dp = kind(0.d0)
        real(dp), allocatable :: s(:, :)
        real(dp), allocatable :: y(:, :)
        real(dp), allocatable :: gradient_prev(:)
        real(dp), allocatable :: parms_prev(:)
        real(dp), parameter :: eps = 1.0e-10_dp
        integer :: curvature_index = 0

    public initialize_olbfgs, olbfgs_iteration, update_hessian

contains

    subroutine initialize_olbfgs(num_pars, history_size) 
        integer, intent(in) :: num_pars
        integer, intent(in) :: history_size
        integer :: i

        ! initial setup upon first call
        allocate(s(history_size, num_pars))
        allocate(y(history_size, num_pars))
        allocate(gradient_prev(num_pars))
        allocate(parms_prev(num_pars))

        gradient_prev = (/(0, i=1, num_pars)/)
        parms_prev = (/(0, i=1, num_pars)/)

    end subroutine

    subroutine olbfgs_iteration(parameters, gradient, step_size, iteration)
        real(dp), intent(inout) :: parameters(:)
        real(dp), intent(in) :: gradient(:)
        real(dp), intent(in) :: step_size
        integer, intent(in) :: iteration

        integer :: n
        
        ! local data
        real(dp), allocatable :: p(:)
        allocate(p(n))

        n = size(s(1, :))

        ! compute initial search direction
        p = initial_direction(gradient, iteration)

        ! save previous gradient and parameters value
        gradient_prev = gradient(1:n)
        parms_prev = parameters(1:n)

        ! update parmaters
        parameters(1:n) = parameters(1:n) + step_size * p(1:n)

    end subroutine

    subroutine update_hessian(parameters, gradient)
        real(dp), dimension(:), intent(in) :: parameters
        real(dp), dimension(:), intent(in) :: gradient
        integer :: m, n
        m = size(s(:, 1))
        n = size(s(1, :))

        ! update head of "queue" used to store curvature pairs
        curvature_index = mod(curvature_index , m) + 1

        ! update curvature pairs s, y according to algorithm
        s(curvature_index, :) = parameters - parms_prev
        y(curvature_index, :) = gradient(1:n) - gradient_prev
    end subroutine

    function initial_direction(gradient, iteration)
        real(dp), allocatable :: initial_direction(:) 
        real(dp), intent(in) :: gradient(:)
        integer, intent(in) :: iteration

        ! local variables
        integer :: n, m, i, head
        real(dp) :: alpha, beta, tot
        real(dp), allocatable :: alphas(:)
        n = size(s(1, :)) ! number of parameters
        m = size(s(:, 1)) ! history size

        ! allocate return value
        allocate(initial_direction(n))

        ! allocate local data
        allocate(alphas(min(m, iteration)))

        ! perform hessian approximation using history
        initial_direction = -gradient(1:n)

        ! first of two-loop recursion
        ! TODO: figure out way to use curvature_meow as queue
        do i=min(m, iteration),1,-1
            head = transform_index(i, iteration)
            alpha = dot_product(s(head, :), initial_direction) &
                / dot_product(s(head, :), y(head, :))
            initial_direction = initial_direction - alpha * y(head, :)
            alphas(head) = alpha
        end do

        ! scale search direction
        if (iteration == 1) then
            initial_direction = initial_direction * eps
        else
            tot = 0.0_dp
            do i=1, min(m, iteration)
                head = transform_index(i, iteration)
                tot = tot + dot_product(s(head, :), y(head, :)) &
                    / dot_product(y(head, :), y(head, :))
            end do
            initial_direction = initial_direction * tot/min(m, iteration)
        end if

        ! second of two-loop recursion
        do i=1, min(m, iteration)
            head = transform_index(i, iteration)
            alpha = alphas(head)
            beta = dot_product(y(head, :), initial_direction) &
                / dot_product(y(head, :), s(head, :))
            initial_direction = initial_direction + (alpha - beta) * s(head, :)
        end do

    end function 

    function transform_index(index, iteration)
        integer, intent(in) :: index, iteration
        integer :: n, m, transform_index

        ! get history size, problem dimensionality
        m = size(s(:, 1))
        n = size(s(1, :))

        ! determine true index into s, y
        if (iteration <= m) then
            transform_index = index
        else
            transform_index = mod(index - 1 + curvature_index, m) + 1
        end if

    end function

end module
