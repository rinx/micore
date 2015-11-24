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
!   20151118 : First Version
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
  real(R_), parameter :: threshold = 0.000000001_R_
  real(R_), parameter :: diff_thre = 0.000000001_R_
  ! max # of iteration
  integer, parameter :: max_iter = 99999

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

  ! simple insert sort
  function insert_sort(dat)
    real(R_) :: dat(:)
    real(R_), allocatable :: insert_sort(:)
    real(R_) :: tmp
    integer  :: i, j

    do i = 2, size(dat)
      tmp = dat(i)
      if (dat(i-1) > tmp) then
        j = i
        do
          dat(j) = dat(j-1)
          j = j - 1
          if (j < 1 .or. dat(j-1) <= tmp) exit
        end do
        dat(j) = tmp
      end if
    end do

    allocate (insert_sort(size(dat)))
    insert_sort(:) = dat(:)
  end function insert_sort

  ! select unique elements from array
  function select_uniq_elems(dat)
    real(R_) :: dat(:)
    real(R_), allocatable :: select_uniq_elems(:)
    integer :: i, cnt

    dat(:) = insert_sort(dat(:))
    cnt = 1
    do i = 2, size(dat)
      if (dat(i-1) /= dat(i)) then
        cnt = cnt + 1
      end if
    end do
    allocate (select_uniq_elems(cnt))
    cnt = 1
    select_uniq_elems(1) = dat(1)
    do i = 2, size(dat)
      if (select_uniq_elems(cnt) /= dat(i)) then
        cnt = cnt + 1
        select_uniq_elems(cnt) = dat(i)
      end if
    end do
  end function select_uniq_elems

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

    call grid_idx_loc(xtab, x, ix, ax, 1)

    ! Trapezoidal difference
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

  ! easy alias for akima_coefs
  function akima_intp(xtab, ytab, x)
    real(R_) :: xtab(:)
    real(R_) :: ytab(:)
    real(R_) :: x
    real(R_) :: akima_intp
    real(R_) :: c1, c2, c3

    call akima_coefs(xtab, ytab, x, akima_intp, c1, c2, c3)
  end function akima_intp

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

  ! 1D-linear algebra solver for limited case (only for 2x2-matrix)
  function solve_1Dlinalg_choles(A, b) result(x)
    real(R_) :: A(2,2)
    real(R_) :: b(2)
    real(R_) :: x(2)
    real(R_) :: L(2,2), z(2)

    L(1,1) = sqrt(A(1,1))
    L(1,2) = 0
    L(2,1) = A(1,2) / L(1,1)
    L(2,2) = sqrt(A(2,2) - L(2,1)**2)

    z(1) = b(1) / L(1,1)
    z(2) = (b(2) - L(2,1) * z(1)) / L(2,2)

    x(2) = z(2) / L(2,2)
    x(1) = (z(1) - L(2,1) * x(2)) / L(1,1)
  end function solve_1Dlinalg_choles

  ! estimate initial tau and cder by least square method
  function estimate_initial_values(lut_refs1, lut_refs2, tau_arr, cder_arr, obs_ref)
    real(R_) :: lut_refs1(:), lut_refs2(:) ! reflectances of lut
    real(R_) :: tau_arr(:), cder_arr(:)    ! tau and cder in lut
    real(R_) :: obs_ref(:)
    real(R_) :: estimate_initial_values(2)
    integer :: minind

    minind = minloc((lut_refs1 - obs_ref(1))**2 + (lut_refs2 - obs_ref(2))**2, 1)

    estimate_initial_values(1) = tau_arr(minind)
    estimate_initial_values(2) = cder_arr(minind)
  end function estimate_initial_values

  ! estimate reflectances from cloud properties by using look-up table
  subroutine estimate_refs(lut_refs1, lut_refs2, tau_arr, cder_arr, tau, cder, est_refs, akic)
    real(R_), intent(in)  :: lut_refs1(:), lut_refs2(:) ! reflectances of lut
    real(R_), intent(in)  :: tau_arr(:), cder_arr(:)    ! tau and cder in lut
    real(R_), intent(in)  :: tau, cder
    real(R_), intent(out) :: est_refs(2)
    real(R_), intent(out) :: akic(12) ! an array of akima coefficients
    real(R_), allocatable :: unq_tau(:), unq_cder(:) ! unique tau and cder
    real(R_), allocatable :: tmp_ref1(:), tmp_ref2(:)
    real(R_), allocatable :: tmp_ref(:,:)
    integer  :: i, j

    ! extract unique values
    allocate (unq_tau(size(select_uniq_elems(tau_arr))))
    allocate (unq_cder(size(select_uniq_elems(cder_arr))))
    unq_tau(:)  = select_uniq_elems(tau_arr(:))
    unq_cder(:) = select_uniq_elems(cder_arr(:))

    allocate(tmp_ref1(size(unq_cder)), tmp_ref2(size(unq_cder)))
    allocate(tmp_ref(size(unq_tau),2))

    ! interpolation with CDER
    do i = 1, size(unq_tau)
      do j = 1, size(unq_cder)
        tmp_ref1(j) = lut_refs1(size(unq_cder) * (i-1) + j)
        tmp_ref2(j) = lut_refs2(size(unq_cder) * (i-1) + j)
      end do
      tmp_ref(i,1) = akima_intp(unq_cder, tmp_ref1, cder)
      tmp_ref(i,2) = akima_intp(unq_cder, tmp_ref2, cder)
    end do

    ! interpolation with TAU
    call akima_coefs(unq_tau, tmp_ref(:,1), tau, est_refs(1), akic(1), akic(2), akic(3))
    call akima_coefs(unq_tau, tmp_ref(:,2), tau, est_refs(2), akic(4), akic(5), akic(6))

    deallocate(tmp_ref1, tmp_ref2, tmp_ref)

    allocate(tmp_ref1(size(unq_tau)), tmp_ref2(size(unq_tau)))
    allocate(tmp_ref(size(unq_cder),2))

    ! interpolation with TAU
    do i = 1, size(unq_cder)
      do j = 1, size(unq_tau)
        tmp_ref1(j) = lut_refs1(size(unq_cder) * (j-1) + i)
        tmp_ref2(j) = lut_refs2(size(unq_cder) * (j-1) + i)
      end do
      tmp_ref(i,1) = akima_intp(unq_tau, tmp_ref1, tau)
      tmp_ref(i,2) = akima_intp(unq_tau, tmp_ref2, tau)
    end do

    ! interpolation with CDER
    ! tmp_ref1 is for dummy
    call akima_coefs(unq_cder, tmp_ref(:,1), cder, tmp_ref1(1), akic(7), akic(8), akic(9))
    call akima_coefs(unq_cder, tmp_ref(:,2), cder, tmp_ref1(2), akic(10), akic(11), akic(12))
    deallocate (unq_tau, unq_cder)
    deallocate(tmp_ref1, tmp_ref2, tmp_ref)
  end subroutine estimate_refs

  ! cost function J
  !   J = (R_obs1 - R_est1)^2 + (R_obs2 - R_est2)^2
  function cost_func(obs_ref, est_ref)
    real(R_) :: obs_ref(:)
    real(R_) :: est_ref(:)
    real(R_) :: cost_func

    cost_func = (obs_ref(1) - est_ref(1))**2 + (obs_ref(2) - est_ref(2))**2
  end function cost_func

  ! update cloud properties
  function update_cloud_properties(obs_ref, est_ref, cps, akic)
    real(R_) :: obs_ref(2), est_ref(2)
    real(R_) :: cps(2)
    real(R_) :: akic(12)
    real(R_) :: update_cloud_properties(2)
    real(R_) :: k(2,2) ! Jacobian matrix
    real(R_) :: refdiff_vec(2) ! difference vector of reflectances

    ! Jacobian matrix
    k(1,1) = akic(1)  + cps(1) * (2 * akic(2)  + cps(1) * 3 * akic(3))  ! (dR_1)/(dtau) at tau = tau_i
    k(1,2) = akic(7)  + cps(2) * (2 * akic(8)  + cps(2) * 3 * akic(9))  ! (dR_1)/(dre)  at re  = re_i
    k(2,1) = akic(4)  + cps(1) * (2 * akic(5)  + cps(1) * 3 * akic(6))  ! (dR_2)/(dtau) at tau = tau_i
    k(2,2) = akic(10) + cps(2) * (2 * akic(11) + cps(2) * 3 * akic(12)) ! (dR_2)/(dre)  at re  = re_i

    refdiff_vec(:) = obs_ref(:) - est_ref(:)

    update_cloud_properties(:) = cps(:) + &
      solve_1Dlinalg_choles(matmul(transpose(k(:,:)), k(:,:)), &
      matmul(transpose(k(:,:)), refdiff_vec(:)))
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
    real(R_) :: akic(12) ! an array of akima coefficients
    real(R_) :: prev_cost = 100.0_R_
    real(R_) :: best_cost = 100.0_R_
    real(R_) :: best_cps(2) = (/0.0_R_, 0.0_R_/)
    integer :: i

    ! initialization
    call separate_lut(lut, lut_refs1, lut_refs2, tau_arr, cder_arr)
    cps = estimate_initial_values(lut_refs1, lut_refs2, tau_arr, cder_arr, obs_ref)

    ! main loop of optimal estimation
    do i = 1, max_iter
      if (verbose_flag) write (*,*) "# iterate step: ", i

      if (verbose_flag) then
        write (*,*) "# estimated TAU:  ", cps(1)
        write (*,*) "# estimated CDER: ", cps(2)
      end if

      call estimate_refs(lut_refs1, lut_refs2, tau_arr, cder_arr, cps(1), cps(2), est_ref, akic)
      if (verbose_flag) then
        write (*,*) "# observed  REF1: ", obs_ref(1)
        write (*,*) "# observed  REF2: ", obs_ref(2)
        write (*,*) "# estimated REF1: ", est_ref(1)
        write (*,*) "# estimated REF2: ", est_ref(2)
      end if

      cost_res = cost_func(obs_ref, est_ref)
      if (verbose_flag) then
        write (*,*) "# COST: ", cost_res
      end if
      if (cost_res < threshold .or. abs(prev_cost - cost_res) < diff_thre) exit
      prev_cost = cost_res

      if (cost_res < best_cost) then
        best_cost = cost_res
        best_cps(:) = cps(:)
      endif

      cps = update_cloud_properties(obs_ref, est_ref, cps, akic)
    end do

    if (cost_res < best_cost) then
      best_cps(:) = cps(:)
    else
      cost_res = best_cost
    endif

    tau  = best_cps(1)
    cder = best_cps(2)
  end subroutine micore_retrieval
end module micore_core

