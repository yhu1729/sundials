! -----------------------------------------------------------------------------
! Programmer(s): Yifan Hu @ UMBC
! -----------------------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2025-2026, Lawrence Livermore National Security,
! University of Maryland Baltimore County, and the SUNDIALS contributors.
! Copyright (c) 2013-2025, Lawrence Livermore National Security
! and Southern Methodist University.
! Copyright (c) 2002-2013, Lawrence Livermore National Security.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! -----------------------------------------------------------------------------
! Driver of radau5.f (E. Hairer and G. Wanner) for the reference files of
! firk_test_radau5_steps (FIRKODE design milestone M2). run_radau5.py
! downloads radau5.f and its linear algebra, compiles them with this driver
! and runs it; none of them is part of SUNDIALS.
!
! Usage: radau5_driver PROBLEM RTOL QUOT OUTFILE STEPS
!
!   PROBLEM  rober (Robertson kinetics to t = 1e11), vdpol (van der Pol in
!            the form of Hairer and Wanner, eps = 1e-6, i.e. mu = 1000, to
!            t = 2) or orego (the Oregonator to t = 360)
!   RTOL     the nominal relative tolerance given to RADAU5, which replaces
!            it by 0.1 RTOL^(2/3)
!   QUOT     ATOL / RTOL, which RADAU5 keeps
!   OUTFILE  the reference file to write
!   STEPS    1 to write the accepted steps (t_n, h_n) to 7 digits, 0 not to
!
! RADAU5 runs with its default parameters (IWORK and WORK zero), the
! analytic Jacobian, scalar tolerances and the initial step 1e-6.
! -----------------------------------------------------------------------------

module radau5_driver_out
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  integer, parameter :: dp = real64
  integer :: out_unit = 10
  logical :: write_steps = .false.
end module radau5_driver_out

program radau5_driver
  use radau5_driver_out
  implicit none
  integer, parameter :: nmax = 3
  integer, parameter :: lwork = nmax*(4*nmax + 12) + 20
  integer, parameter :: liwork = 3*nmax + 20
  character(len=32) :: prob
  character(len=512) :: arg, outfile
  real(dp) :: rtol(1), atol(1), rtol_nom, atol_nom, quot, x, xend, h, h0, expm
  real(dp) :: y(nmax), work(lwork), rpar(1)
  integer :: iwork(liwork), ipar(1), n, idid, i, steps_flag, ios
  ! RADAU5 may overwrite its integer arguments (MLJAC and MUJAC for a full
  ! Jacobian): variables, not constants
  integer :: itol, ijac, mljac, mujac, imas, mlmas, mumas, iout
  external :: fcn_rober, jac_rober, fcn_vdpol, jac_vdpol, fcn_orego, &
    jac_orego, mas_dummy, solout, radau5

  if (command_argument_count() /= 5) then
    write (*, '(a)') 'usage: radau5_driver PROBLEM RTOL QUOT OUTFILE STEPS'
    stop 2
  end if
  call get_command_argument(1, prob)
  call get_command_argument(2, arg)
  read (arg, *, iostat=ios) rtol_nom
  if (ios /= 0) stop 2
  call get_command_argument(3, arg)
  read (arg, *, iostat=ios) quot
  if (ios /= 0) stop 2
  call get_command_argument(4, outfile)
  call get_command_argument(5, arg)
  read (arg, *, iostat=ios) steps_flag
  if (ios /= 0) stop 2

  atol_nom = quot*rtol_nom
  rtol(1) = rtol_nom
  atol(1) = atol_nom
  h0 = 1.0e-6_dp
  h = h0
  x = 0.0_dp
  y = 0.0_dp
  select case (trim(prob))
  case ('rober')
    n = 3
    y(1) = 1.0_dp
    xend = 1.0e11_dp
  case ('vdpol')
    n = 2
    y(1) = 2.0_dp
    xend = 2.0_dp
  case ('orego')
    n = 3
    y(1:3) = [1.0_dp, 2.0_dp, 3.0_dp]
    xend = 360.0_dp
  case default
    write (*, '(2a)') 'unknown problem ', trim(prob)
    stop 2
  end select
  iwork = 0
  work = 0.0_dp
  itol = 0
  ijac = 1
  mljac = n
  mujac = 0
  imas = 0
  mlmas = n
  mumas = 0
  iout = 1

  open (newunit=out_unit, file=trim(outfile), status='replace', &
        action='write', iostat=ios)
  if (ios /= 0) stop 3
  write_steps = (steps_flag == 1)
  write (out_unit, '(2a)') 'problem ', trim(prob)
  write (out_unit, '(a,i0)') 'n ', n
  write (out_unit, '(a,es25.17e3)') 'tend ', xend
  write (out_unit, '(a,es25.17e3)') 'h0 ', h0
  write (out_unit, '(a,es25.17e3)') 'rtol_nominal ', rtol_nom
  write (out_unit, '(a,es25.17e3)') 'atol_nominal ', atol_nom

  select case (trim(prob))
  case ('rober')
    call radau5(n, fcn_rober, x, y, xend, h, rtol, atol, itol, jac_rober, &
                ijac, mljac, mujac, mas_dummy, imas, mlmas, mumas, solout, &
                iout, work, lwork, iwork, liwork, rpar, ipar, idid)
  case ('vdpol')
    call radau5(n, fcn_vdpol, x, y, xend, h, rtol, atol, itol, jac_vdpol, &
                ijac, mljac, mujac, mas_dummy, imas, mlmas, mumas, solout, &
                iout, work, lwork, iwork, liwork, rpar, ipar, idid)
  case ('orego')
    call radau5(n, fcn_orego, x, y, xend, h, rtol, atol, itol, jac_orego, &
                ijac, mljac, mujac, mas_dummy, imas, mlmas, mumas, solout, &
                iout, work, lwork, iwork, liwork, rpar, ipar, idid)
  end select

  ! the tolerances RADAU5 used, 0.1 RTOL^(2/3) and ATOL with the same ratio,
  ! in its own arithmetic (it restores RTOL and ATOL on return)
  expm = 2.0_dp/3.0_dp
  write (out_unit, '(a,es25.17e3)') 'rtol ', 0.1_dp*rtol_nom**expm
  write (out_unit, '(a,es25.17e3)') 'atol ', &
    0.1_dp*rtol_nom**expm*(atol_nom/rtol_nom)
  write (out_unit, '(a,i0)') 'idid ', idid
  write (out_unit, '(a,i0)') 'nfcn ', iwork(14)
  write (out_unit, '(a,i0)') 'njac ', iwork(15)
  write (out_unit, '(a,i0)') 'nstep ', iwork(16)
  write (out_unit, '(a,i0)') 'naccpt ', iwork(17)
  write (out_unit, '(a,i0)') 'nrejct ', iwork(18)
  write (out_unit, '(a,i0)') 'ndec ', iwork(19)
  write (out_unit, '(a,i0)') 'nsol ', iwork(20)
  write (out_unit, '(a,es25.17e3)') 't ', x
  do i = 1, n
    write (out_unit, '(a,i0,1x,es25.17e3)') 'y ', i, y(i)
  end do
  write (out_unit, '(a)') 'end'
  close (out_unit)
  if (idid /= 1) stop 1
end program radau5_driver

! the accepted steps: SOLOUT is called at the initial point (NR = 1) and
! after every accepted step
subroutine solout(nr, xold, x, y, cont, lrc, n, rpar, ipar, irtrn)
  use radau5_driver_out
  implicit none
  integer, intent(in) :: nr, lrc, n, ipar(*)
  integer, intent(inout) :: irtrn
  real(dp), intent(in) :: xold, x, y(n), cont(lrc), rpar(*)
  if (write_steps .and. nr > 1) then
    write (out_unit, '(a,es14.6e3,1x,es14.6e3)') 'step ', x, x - xold
  end if
  irtrn = 0
end subroutine solout

subroutine mas_dummy(n, am, lmas, rpar, ipar)
  use radau5_driver_out, only: dp
  implicit none
  integer, intent(in) :: n, lmas, ipar(*)
  real(dp), intent(inout) :: am(lmas, n)
  real(dp), intent(in) :: rpar(*)
end subroutine mas_dummy

! Robertson kinetics
subroutine fcn_rober(n, x, y, f, rpar, ipar)
  use radau5_driver_out, only: dp
  implicit none
  integer, intent(in) :: n, ipar(*)
  real(dp), intent(in) :: x, y(n), rpar(*)
  real(dp), intent(out) :: f(n)
  f(1) = -0.04_dp*y(1) + 1.0e4_dp*y(2)*y(3)
  f(2) = 0.04_dp*y(1) - 1.0e4_dp*y(2)*y(3) - 3.0e7_dp*y(2)*y(2)
  f(3) = 3.0e7_dp*y(2)*y(2)
end subroutine fcn_rober

subroutine jac_rober(n, x, y, dfy, ldfy, rpar, ipar)
  use radau5_driver_out, only: dp
  implicit none
  integer, intent(in) :: n, ldfy, ipar(*)
  real(dp), intent(in) :: x, y(n), rpar(*)
  real(dp), intent(out) :: dfy(ldfy, n)
  dfy(1, 1) = -0.04_dp
  dfy(1, 2) = 1.0e4_dp*y(3)
  dfy(1, 3) = 1.0e4_dp*y(2)
  dfy(2, 1) = 0.04_dp
  dfy(2, 2) = -1.0e4_dp*y(3) - 6.0e7_dp*y(2)
  dfy(2, 3) = -1.0e4_dp*y(2)
  dfy(3, 1) = 0.0_dp
  dfy(3, 2) = 6.0e7_dp*y(2)
  dfy(3, 3) = 0.0_dp
end subroutine jac_rober

! van der Pol, y1' = y2, y2' = ((1 - y1^2) y2 - y1) / eps
subroutine fcn_vdpol(n, x, y, f, rpar, ipar)
  use radau5_driver_out, only: dp
  implicit none
  integer, intent(in) :: n, ipar(*)
  real(dp), intent(in) :: x, y(n), rpar(*)
  real(dp), intent(out) :: f(n)
  real(dp), parameter :: eps = 1.0e-6_dp
  f(1) = y(2)
  f(2) = ((1.0_dp - y(1)**2)*y(2) - y(1))/eps
end subroutine fcn_vdpol

subroutine jac_vdpol(n, x, y, dfy, ldfy, rpar, ipar)
  use radau5_driver_out, only: dp
  implicit none
  integer, intent(in) :: n, ldfy, ipar(*)
  real(dp), intent(in) :: x, y(n), rpar(*)
  real(dp), intent(out) :: dfy(ldfy, n)
  real(dp), parameter :: eps = 1.0e-6_dp
  dfy(1, 1) = 0.0_dp
  dfy(1, 2) = 1.0_dp
  dfy(2, 1) = (-2.0_dp*y(1)*y(2) - 1.0_dp)/eps
  dfy(2, 2) = (1.0_dp - y(1)**2)/eps
end subroutine jac_vdpol

! the Oregonator
subroutine fcn_orego(n, x, y, f, rpar, ipar)
  use radau5_driver_out, only: dp
  implicit none
  integer, intent(in) :: n, ipar(*)
  real(dp), intent(in) :: x, y(n), rpar(*)
  real(dp), intent(out) :: f(n)
  f(1) = 77.27_dp*(y(2) + y(1)*(1.0_dp - 8.375e-6_dp*y(1) - y(2)))
  f(2) = (y(3) - (1.0_dp + y(1))*y(2))/77.27_dp
  f(3) = 0.161_dp*(y(1) - y(3))
end subroutine fcn_orego

subroutine jac_orego(n, x, y, dfy, ldfy, rpar, ipar)
  use radau5_driver_out, only: dp
  implicit none
  integer, intent(in) :: n, ldfy, ipar(*)
  real(dp), intent(in) :: x, y(n), rpar(*)
  real(dp), intent(out) :: dfy(ldfy, n)
  dfy(1, 1) = 77.27_dp*(1.0_dp - 2.0_dp*8.375e-6_dp*y(1) - y(2))
  dfy(1, 2) = 77.27_dp*(1.0_dp - y(1))
  dfy(1, 3) = 0.0_dp
  dfy(2, 1) = -y(2)/77.27_dp
  dfy(2, 2) = -(1.0_dp + y(1))/77.27_dp
  dfy(2, 3) = 1.0_dp/77.27_dp
  dfy(3, 1) = 0.161_dp
  dfy(3, 2) = 0.0_dp
  dfy(3, 3) = -0.161_dp
end subroutine jac_orego
