! Copyright (c) 2024 Roberto E. Navarro <roberto.navarro@udec.cl>
!                    Departamento de Fisica
!                    Facultad de Ciencias Fisicas y Matematicas
!                    Universidad de Concepcion


module precision_module
  use, intrinsic :: iso_fortran_env, only: dtype => real64
end module precision_module

module particle_pushers
  use precision_module
  implicit none

contains
  function boris_method(v, E, B) result(vnew)
    real(dtype), dimension(3), intent(in) :: v, E, B
    real(dtype), dimension(3) :: vnew

    real(dtype), dimension(3) :: vminus
    real(dtype) :: fp

    vminus = v + E

    fp = 2.0 / (1.0 + dot_product(B,B))

    vnew = vminus + cross_product(vminus, B)
    vnew = vminus + fp * cross_product(vnew, B)

    vnew = vnew + E
  end function boris_method


  ! Calculate the cross product between two 3D vectors `a` and `b`.
  ! Arguments
  ! ---------
  !   a, b : dimension(3) vectors
  !
  ! Returns
  ! -------
  !   dimension(3) vector representing the cross product of `a` and `b`.
  pure function cross_product(a, b) result(cross)
    real(dtype), dimension(3), intent(in) :: a, b
    real(dtype), dimension(3) :: cross

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product
end module particle_pushers


module fields
  use precision_module
  implicit none

contains
  pure function magnetic_dipole(x) result(B)
    real(dtype), dimension(3), intent(in) :: x
    real(dtype), dimension(3) :: B
    real(dtype) :: r

    r = norm2(x)

    ! We assume a magnetic dipole of moment vec(m)=[0,0,1]
    B = [real(dtype):: 0, 0, -1.0/r**3] + 3 * x(3) * x / r**5
  end function magnetic_dipole

end module fields


program main
  use particle_pushers
  use fields
  implicit none

  real(dtype), parameter :: q = 1.0, m = 1.0

  real(dtype), parameter :: dt = 0.01
  integer, parameter :: nsteps = 64000

  real(dtype), parameter :: halfdt = 0.5*dt

  real(dtype), dimension(3) :: x, v, E, B
  integer :: i

  ! Initial conditions
  x = [1.0, 0.0, 0.0]
  v = [0.0, 0.0, 0.09]

  ! Static EM fields (normalized)
  E = [0.0, 0.0, 0.0] * halfdt * q / m

  open(11, file="output.dat", action="write", access="stream")
  do i=1,nsteps
     x = x + halfdt * v         ! x(t)

     B = magnetic_dipole(x) * halfdt * q / m

     v = boris_method(v, E, B)  ! v(t+dt/2)
     x = x + halfdt * v         ! x(t+dt/2)

     write(11) x, v
  end do

  close(11)

end program main
