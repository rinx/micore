!
! MICO-RE
!   [Minimal Implementation of Cloud Optical Retrieval]
!
! micore_main.f90
!
! Author: Rintaro Okamura
!
! Description:
!   An interface of MICO-RE cloud optical retrieval code
!
! Dependency:
!   micore_core.f90
!
! ChangeLog
!   20151127 : First Version
!

program micore
  use micore_core

  implicit none

  integer, parameter :: IUI = 10

  character(256) :: argv(3)
  character(512) :: argmsg
  integer :: narg
  integer :: irecl, ios, i, ivar, cnt

  real(R_), allocatable :: lut(:,:)
  real(R_) :: obs_ref(2)
  real(R_) :: tau, cder
  real(R_) :: cost
  real(R4_) :: wrk(4)

  ! get commandline arguments
  narg = 3
  argmsg = "Usage: micore lutfile obs_ref1 obs_ref2"
  call get_cmd_args(narg, argv, argmsg)
  read (argv(2), *) obs_ref(1)
  read (argv(3), *) obs_ref(2)

  ! debug
  if (verbose_flag) then
    write (*,*) '# OBS_REF1: ', obs_ref(1)
    write (*,*) '# OBS_REF2: ', obs_ref(2)
  end if

  ! read lut-file
  inquire (iolength=irecl) (1.0_R4_, ivar=1, 4)
  open (IUI, file=argv(1), access='direct', form='unformatted', status='old', recl=irecl)

  i = 1
  do
    read (IUI, rec=i, iostat=ios) wrk
    if (ios /= 0) exit
    i = i + 1
  end do

  cnt = i - 1

  allocate (lut(cnt,4))

  do i = 1, cnt
    read (IUI, rec=i, iostat=ios) wrk
    if (ios /= 0) exit
    lut(i,:) = wrk(:)
  end do

  if (verbose_flag) then
    write (*,*) "# Look-up table"
    write (*,*) "# TAU CDER REF1 REF2"
    do i = 1, cnt
      write (*,*) lut(i,:)
    end do
  end if

  ! call main routine of retrieval
  call micore_retrieval(lut, obs_ref, tau, cder, cost)

  write (*,*) 'TAU:  ', tau
  write (*,*) 'CDER: ', cder
  write (*,*) 'COST: ', cost

  ! finalize
  deallocate (lut)
  close(IUI)

  stop

end program micore

