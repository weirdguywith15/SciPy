program asteroid_belt_simulation
  implicit none

  ! Constants
  real, parameter :: G = 6.67430E-11        ! Gravitational constant (m^3 kg^-1 s^-2)
  real, parameter :: M_SUN = 1.989E30       ! Mass of the Sun (kg)
  real, parameter :: AU = 1.496E11          ! Astronomical Unit (m)
  real, parameter :: DAY = 86400.0          ! Seconds in a day
  integer, parameter :: N_AST = 100         ! Number of asteroids
  integer, parameter :: N_STEPS = 720       ! Time steps (30 days, ~1 hr per step)
  real, parameter :: DT = DAY / 24.0        ! Time step (1 hour in seconds)

  ! Variables
  real :: x(N_AST), y(N_AST), z(N_AST)      ! Positions (m)
  real :: vx(N_AST), vy(N_AST), vz(N_AST)   ! Velocities (m/s)
  real :: ax(N_AST), ay(N_AST), az(N_AST)   ! Accelerations (m/s^2)
  real :: r, t
  integer :: i, step

  ! Output file
  open(unit=10, file='asteroid_orbits.dat', status='replace')

  ! Initialize asteroid positions and velocities
  call initialize_asteroids(x, y, z, vx, vy, vz)
  
  ! Time loop
  t = 0.0
  do step = 1, N_STEPS
    ! Compute accelerations
    call compute_accelerations(x, y, z, ax, ay, az)
    
    ! Update positions and velocities using RK4
    call rk4_step(x, y, z, vx, vy, vz, ax, ay, az)
    
    ! Output data every 24 steps (once per day)
    if (mod(step, 24) == 0) then
      write(10, '(F10.2, 100(3F15.6))') t/DAY, (x(i)/AU, y(i)/AU, z(i)/AU, i=1,N_AST)
    end if
    
    t = t + DT
  end do

  close(10)
  print *, 'Simulation complete. Data written to asteroid_orbits.dat'

contains

  ! Subroutine to initialize asteroid positions and velocities
  subroutine initialize_asteroids(x, y, z, vx, vy, vz)
    real, intent(out) :: x(:), y(:), z(:), vx(:), vy(:), vz(:)
    real :: r, theta, v_circ
    integer :: i
    
    do i = 1, N_AST
      ! Random radius between 2.0 and 3.5 AU (typical asteroid belt range)
      r = (2.0 + 1.5 * real(i-1)/(N_AST-1)) * AU
      theta = 2.0 * 3.14159 * real(i-1)/(N_AST-1)  ! Spread evenly in angle
      
      ! Position in xy-plane (circular orbits initially)
      x(i) = r * cos(theta)
      y(i) = r * sin(theta)
      z(i) = 0.0  ! Assume planar belt for simplicity
      
      ! Circular velocity: v = sqrt(GM/r)
      v_circ = sqrt(G * M_SUN / r)
      vx(i) = -v_circ * sin(theta)  ! Tangential velocity
      vy(i) = v_circ * cos(theta)
      vz(i) = 0.0
    end do
  end subroutine initialize_asteroids

  ! Subroutine to compute gravitational accelerations
  subroutine compute_accelerations(x, y, z, ax, ay, az)
    real, intent(in) :: x(:), y(:), z(:)
    real, intent(out) :: ax(:), ay(:), az(:)
    real :: r
    integer :: i
    
    do i = 1, N_AST
      r = sqrt(x(i)2 + y(i)2 + z(i)**2)
      ! Acceleration due to Sun (only central force considered)
      ax(i) = -G * M_SUN * x(i) / r**3
      ay(i) = -G * M_SUN * y(i) / r**3
      az(i) = -G * M_SUN * z(i) / r**3
    end do
  end subroutine compute_accelerations

  ! Subroutine to perform one RK4 integration step
  subroutine rk4_step(x, y, z, vx, vy, vz, ax, ay, az)
    real, intent(inout) :: x(:), y(:), z(:), vx(:), vy(:), vz(:)
    real, intent(in) :: ax(:), ay(:), az(:)
    real :: k1x(N_AST), k1y(N_AST), k1z(N_AST)
    real :: k1vx(N_AST), k1vy(N_AST), k1vz(N_AST)
    real :: k2x(N_AST), k2y(N_AST), k2z(N_AST)
    real :: k2vx(N_AST), k2vy(N_AST), k2vz(N_AST)
    real :: k3x(N_AST), k3y(N_AST), k3z(N_AST)
    real :: k3vx(N_AST), k3vy(N_AST), k3vz(N_AST)
    real :: k4x(N_AST), k4y(N_AST), k4z(N_AST)
    real :: k4vx(N_AST), k4vy(N_AST), k4vz(N_AST)
    real :: xt(N_AST), yt(N_AST), zt(N_AST)
    real :: vxt(N_AST), vyt(N_AST), vzt(N_AST)
    real :: axt(N_AST), ayt(N_AST), azt(N_AST)
    integer :: i

    ! k1: initial derivatives
    k1x = vx; k1y = vy; k1z = vz
    k1vx = ax; k1vy = ay; k1vz = az

    ! k2: midpoint with k1
do i = 1, N_AST
      xt(i) = x(i) + 0.5 * DT * k1x(i)
      yt(i) = y(i) + 0.5 * DT * k1y(i)
      zt(i) = z(i) + 0.5 * DT * k1z(i)
    end do
    call compute_accelerations(xt, yt, zt, axt, ayt, azt)
    k2x = vx + 0.5 * DT * k1vx
    k2y = vy + 0.5 * DT * k1vy
    k2z = vz + 0.5 * DT * k1vz
    k2vx = axt; k2vy = ayt; k2vz = azt

    ! k3: midpoint with k2
    do i = 1, N_AST
      xt(i) = x(i) + 0.5 * DT * k2x(i)
      yt(i) = y(i) + 0.5 * DT * k2y(i)
      zt(i) = z(i) + 0.5 * DT * k2z(i)
    end do
    call compute_accelerations(xt, yt, zt, axt, ayt, azt)
    k3x = vx + 0.5 * DT * k2vx
    k3y = vy + 0.5 * DT * k2vy
    k3z = vz + 0.5 * DT * k2vz
    k3vx = axt; k3vy = ayt; k3vz = azt

    ! k4: end point with k3
    do i = 1, N_AST
      xt(i) = x(i) + DT * k3x(i)
      yt(i) = y(i) + DT * k3y(i)
      zt(i) = z(i) + DT * k3z(i)
    end do
    call compute_accelerations(xt, yt, zt, axt, ayt, azt)
    k4x = vx + DT * k3vx
    k4y = vy + DT * k3vy
    k4z = vz + DT * k3vz
    k4vx = axt; k4vy = ayt; k4vz = azt

    ! Update positions and velocities
    do i = 1, N_AST
      x(i) = x(i) + (DT/6.0) * (k1x(i) + 2.0*k2x(i) + 2.0*k3x(i) + k4x(i))
      y(i) = y(i) + (DT/6.0) * (k1y(i) + 2.0*k2y(i) + 2.0*k3y(i) + k4y(i))
      z(i) = z(i) + (DT/6.0) * (k1z(i) + 2.0*k2z(i) + 2.0*k3z(i) + k4z(i))
      vx(i) = vx(i) + (DT/6.0) * (k1vx(i) + 2.0*k2vx(i) + 2.0*k3vx(i) + k4vx(i))
      vy(i) = vy(i) + (DT/6.0) * (k1vy(i) + 2.0*k2vy(i) + 2.0*k3vy(i) + k4vy(i))
      vz(i) = vz(i) + (DT/6.0) * (k1vz(i) + 2.0*k2vz(i) + 2.0*k3vz(i) + k4vz(i))
    end do
  end subroutine rk4_step

end program asteroid_belt_simulation
