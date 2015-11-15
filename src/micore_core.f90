!
! MICO-RE
!   [Minimal Implementation of Cloud Optical Retrieval]
!
! micore_core.f90
!
! Author: Rintaro Okamura
!
! Description:
!   main routine and libraries of MICO-RE cloud optical retrieval code
!
! ChangeLog
!   201511xx : First Version
!


module micore_core
  implicit none
  private

  public :: R_, RD_, R4_, verbose_flag
  public :: get_cmd_args
  public :: micore_retrieval

  ! based on HPARX library
  integer, parameter :: R_  = selected_real_kind(13) ! default precision
  integer, parameter :: RD_ = selected_real_kind(13) ! higher  precision
  integer, parameter :: R4_ = selected_real_kind(6)  ! 4-byte real

  ! flag for verbose mode
  logical, parameter :: verbose_flag = .true.

  ! threshold value for convergence of cost function
  real(R_), parameter :: threshold = 0.05_R_
  ! max # of iteration
  integer, parameter :: max_iter = 9999

contains
  ! get commandline arguments
  !   this code is based on HPARX library
  subroutine get_cmd_args(narg, argv, argmsg)
    integer, intent(inout)    :: narg
    character(*), intent(out) :: argv(:)
    character(*), intent(in)  :: argmsg
    integer :: nn, i

    nn = command_argument_count()

    if (nn < narg) then
      write (*,*) argmsg
      stop
    end if

    do i = 1, nn
      call get_command_argument(i, argv(i))
    end do

    narg = nn
  end subroutine get_cmd_args

  ! get location inside of grid
  !   this code is based on HPARX library
  subroutine grid_idx_loc(xg, x, ix, rat, mext)
    real(R_), intent(in)  :: xg(:) ! vector
    real(R_), intent(in)  :: x     ! a value
    integer,  intent(out) :: ix    ! found grid index (1 <= ix <= size(xg)-1)
    real(R_), intent(out) :: rat   ! location between xg(ix) and xg(ix+1)
    !// xg(ix) <= x < xg(ix+1) if xg is  increasing vector and xg(1) <= x < xg(nx)
    !// xg(ix) >= x > xg(ix+1) if xg is deccreasing vector and xg(1) >= x > xg(nx)
    integer,  intent(in), optional :: mext ! if 1, extraporation (rat will be < 0 or > 1)
    integer  :: i, iz

    ! Too few grid points
    iz = size(xg)
    if (iz <= 1) then
       ix = 1
       rat = 0.0_R_

       ! Increasing vector
    else if (xg(1) < xg(iz)) then
       if (x <= xg(1)) then ! out of lower limit
          ix = 1
          rat = 0.0_R_
          if (present(mext)) then
             if (mext == 1) rat = (x - xg(1)) / (xg(2) - xg(1))
          end if
       else if (x >= xg(iz)) then ! out of upper limit
          ix = iz - 1
          rat = 1.0_R_
          if (present(mext)) then
             if (mext == 1) rat = (x - xg(ix)) / (xg(iz) - xg(ix))
          end if
       else                ! in the grid
          ix = 1
          do
             if (iz <= ix + 1) exit
             i = (ix + iz) / 2
             if (x >= xg(i)) then
                ix = i
             else
                iz = i
             end if
          end do
          rat = (x - xg(ix)) / (xg(ix + 1) - xg(ix))
       end if

       ! Decreasing vector
    else
       if (x >= xg(1)) then ! out of lower limit
          ix = 1
          rat = 0.0_R_
          if (present(mext)) then
             if (mext == 1) rat = (x - xg(1)) / (xg(2) - xg(1))
          end if
       else if (x <= xg(iz)) then ! out of upper limit
          ix = iz - 1
          rat = 1.0_R_
          if (present(mext)) then
             if (mext == 1) rat = (x - xg(ix)) / (xg(iz) - xg(ix))
          end if
       else                ! in the grid
          ix = 1
          do
             if (iz <= ix + 1) exit
             i = (ix + iz) / 2
             if (x <= xg(i)) then
                ix = i
             else
                iz = i
             end if
          end do
          rat = (x - xg(ix)) / (xg(ix + 1) - xg(ix))
       end if
    end if
  end subroutine grid_idx_loc

  ! get coefficients of akima interpolation
  !   y = y0 + c1*x + c2*x^2 + c3*x^3
  !   this code is based on HPARX library
  subroutine akima_coefs(xtab, ytab, x, y, c1, c2, c3)
    real(R_), intent(in)  :: xtab(:)
    real(R_), intent(in)  :: ytab(:)
    real(R_), intent(in)  :: x
    real(R_), intent(out) :: y
    real(R_), intent(out) :: c1, c2, c3
    real(R_) :: dydxtab(0:1), d(-2:2), dx, dy, w0, w1, ax
    integer  :: nx, i, ix

    nx = size(ytab)

    call grid_idx_loc(xtab, x, ix, ax)

    !Trapezoidal difference
    if (ix >= 3 .and. ix <= nx - 3) then
       d(-2:2) = (ytab(ix-1:ix+3) - ytab(ix-2:ix+2)) / (xtab(ix-1:ix+3) - xtab(ix-2:ix+2))
    else if (ix <= 2) then
       if (ix == 2) then
          d(-1:2) = (ytab(ix  :ix+3) - ytab(ix-1:ix+2)) / (xtab(ix  :ix+3) - xtab(ix-1:ix+2))
          d(-2) = 2.0_R_ * d(-1) - d(0)
       else
          d(0:2)  = (ytab(ix+1:ix+3) - ytab(ix  :ix+2)) / (xtab(ix+1:ix+3) - xtab(ix  :ix+2))
          d(-1) = 2.0_R_ * d(0)  - d(1)
          d(-2) = 2.0_R_ * d(-1) - d(0)
       end if
    else
       if (ix == nx - 2) then
          d(-2:1) = (ytab(ix-1:ix+2) - ytab(ix-2:ix+1)) / (xtab(ix-1:ix+2) - xtab(ix-2:ix+1))
          d(2) = 2.0_R_ * d(1) - d(0)
       else !if (ix == nx - 1) then
          d(-2:0) = (ytab(ix-1:ix+1) - ytab(ix-2:ix  )) / (xtab(ix-1:ix+1) - xtab(ix-2:ix  ))
          d(1) = 2.0_R_ * d(0) - d(-1)  
          d(2) = 2.0_R_ * d(1) - d(0)
       end if
    end if

    ! Derivative estimates
    do i = 0, 1
       w1 = abs(d(i+1) - d(i  ))
       w0 = abs(d(i-1) - d(i-2))
       if (w1 + w0 == 0.0_R_) then
          dydxtab(i) = (d(i-1) + d(i)) * 0.5_R_
       else
          dydxtab(i) = (w1 * d(i-1) + w0 * d(i)) / (w1 + w0)
       end if
    end do

    dx = xtab(ix+1) - xtab(ix)
    dy = ytab(ix+1) - ytab(ix)
    c1 = dydxtab(0) * dx
    c2 =  3.0_R_ * dy - (c1 + dydxtab(1)*dx) - c1
    c3 = -2.0_R_ * dy + (c1 + dydxtab(1)*dx)
    y = ytab(ix) + ax * (c1 + ax * (c2 + ax * c3))
  end subroutine akima_coefs

  ! separate look-up table into several vectors
  subroutine separate_lut(lut, lut_refs1, lut_refs2, tau_arr, cder_arr)
    real(R_), intent(in)  :: lut(:,:)
    real(R_), intent(out) :: lut_refs1(:), lut_refs2(:) ! reflectances of lut
    real(R_), intent(out) :: tau_arr(:), cder_arr(:)    ! tau and cder in lut
    integer :: i, lutsize

    lutsize = size(lut, 1)

    do i = 1, lutsize
      tau_arr(i)   = lut(i, 1)
      cder_arr(i)  = lut(i, 2)
      lut_refs1(i) = lut(i, 3)
      lut_refs2(i) = lut(i, 4)
    end do
  end subroutine separate_lut

  ! estimate initial tau and cder by least square method
  function estimate_initial_values(lut_refs1, lut_refs2, tau_arr, cder_arr, obs_ref)
    real(R_), intent(in)  :: lut_refs1(:), lut_refs2(:) ! reflectances of lut
    real(R_), intent(in)  :: tau_arr(:), cder_arr(:)    ! tau and cder in lut
    real(R_), intent(in)  :: obs_ref(:)
    real(R_) :: estimate_initial_values(2)
    integer :: minind

    minind = minloc((lut_refs1 - obs_ref(1))**2 + (lut_refs2 - obs_ref(2))**2, 1)

    estimate_initial_values(1) = tau_arr(minind)
    estimate_initial_values(2) = cder_arr(minind)
  end function estimate_initial_values

  ! estimate reflectances from cloud properties by using look-up table
  function estimate_refs(lut_refs1, lut_refs2, tau_arr, cder_arr, tau, cder)
    real(R_), intent(in)  :: lut_refs1(:), lut_refs2(:) ! reflectances of lut
    real(R_), intent(in)  :: tau_arr(:), cder_arr(:)    ! tau and cder in lut
    real(R_), intent(in)  :: tau, cder
    real(R_) :: estimate_refs(2)

    estimate_refs(1) = 0.0
    estimate_refs(2) = 0.0
  end function estimate_refs

  ! cost function J
  !   J = (R_obs1 - R_est1)^2 + (R_obs2 - R_est2)^2
  function cost_func(obs_ref, est_ref)
    real(R_), intent(in) :: obs_ref(:)
    real(R_), intent(in) :: est_ref(:)
    real(R_) :: cost_func

    cost_func = (obs_ref(1) - est_ref(1))**2 - (obs_ref(2) - est_ref(2))**2
  end function cost_func

  ! update cloud properties
  function update_cloud_properties(tau, cder)
    real(R_), intent(in) :: tau, cder
    real(R_) :: update_cloud_properties(2)

    update_cloud_properties(1) = tau
    update_cloud_properties(2) = cder
  end function update_cloud_properties

  ! main routine of retrieval code
  !   input:
  !     lut: look-up table (tau, cder -> ref1, ref2)
  !       (i, 1): tau  (Cloud Optical Thickness)
  !       (i, 2): cder (Cloud Droplet Effective Radius)
  !       (i, 3): reflectance1
  !       (i, 4): reflectance2
  !     obs_ref: array of observed reflectances
  !   output:
  !     tau: estimated cloud optical thickness
  !     cder: estimated cloud droplet effective radius
  subroutine micore_retrieval(lut, obs_ref, tau, cder, cost_res)
    real(R_), intent(in)  :: lut(:,:)
    real(R_), intent(in)  :: obs_ref(:)
    real(R_), intent(out) :: tau, cder
    real(R_), intent(out) :: cost_res
    real(R_) :: lut_refs1(size(lut, 1)), lut_refs2(size(lut, 1)) ! reflectances of lut
    real(R_) :: tau_arr(size(lut,1)), cder_arr(size(lut, 1)) ! tau and cder in lut
    real(R_) :: est_ref(2) ! estimated reflectances
    real(R_) :: cps(2) ! cloud physical parameters
    real(R_) :: prev_cost
    integer :: i

    ! initialization
    call separate_lut(lut, lut_refs1, lut_refs2, tau_arr, cder_arr)
    cps = estimate_initial_values(lut_refs1, lut_refs2, tau_arr, cder_arr, obs_ref)

    if (verbose_flag) then
      write (*,*) "# 1st-estimated TAU:  ", cps(1)
      write (*,*) "# 1st-estimated CDER: ", cps(2)
    end if

    stop

    ! main loop of optimal estimation
    do i = 1, max_iter
      if (verbose_flag) write (*,*) "# iterate step: ", i

      est_ref = estimate_refs(lut_refs1, lut_refs2, tau_arr, cder_arr, tau, cder)
      if (verbose_flag) then
        write (*,*) "# estimated REF1: ", est_ref(1)
        write (*,*) "# estimated REF2: ", est_ref(2)
      end if

      cost_res = cost_func(obs_ref, est_ref)
      if (cost_res < threshold .or. (prev_cost - cost_res) < threshold) exit
      if (verbose_flag) then
        write (*,*) "# COST: ", cost_res
      end if

      cps = update_cloud_properties(cps(1), cps(2))
      if (verbose_flag) then
        write (*,*) "# estimated TAU:  ", cps(1)
        write (*,*) "# estimated CDER: ", cps(2)
      end if
    end do

    tau  = cps(1)
    cder = cps(2)
  end subroutine micore_retrieval
end module micore_core

