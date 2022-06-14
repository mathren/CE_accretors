! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

module run_star_extras

  use star_lib
  use star_def
  use const_def
  use math_lib

  implicit none

  ! these routines are called by the standard run_star check_model
contains

  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! this is the place to set any procedure pointers you want to change
    ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


    ! the extras functions in this file will not be called
    ! unless you set their function pointers as done below.
    ! otherwise we use a null_ version which does nothing (except warn).

    s% extras_startup => extras_startup
    s% extras_start_step => extras_start_step
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve
    s% how_many_extra_history_columns => how_many_extra_history_columns
    s% data_for_extra_history_columns => data_for_extra_history_columns
    s% how_many_extra_profile_columns => how_many_extra_profile_columns
    s% data_for_extra_profile_columns => data_for_extra_profile_columns

    s% how_many_extra_history_header_items => how_many_extra_history_header_items
    s% data_for_extra_history_header_items => data_for_extra_history_header_items
    s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
    s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

  end subroutine extras_controls


  subroutine extras_startup(id, restart, ierr)
    integer, intent(in) :: id
    logical, intent(in) :: restart
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! use 11-16 to avoid overlap with lxtra(1) and lxtra(2) used in run_binary_extras.f90
    if (.not. restart) then
       s% lxtra(11) = .false.
       s% lxtra(12) = .false.
       s% lxtra(13) = .false.
       s% lxtra(14) = .false.
       s% lxtra(15) = .false.
       s% lxtra(16) = .false.
    end if


  end subroutine extras_startup


  integer function extras_start_step(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_start_step = 0
  end function extras_start_step


  ! returns either keep_going, retry, or terminate.
  integer function extras_check_model(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    real(dp) :: error, atol, rtol
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going
    if (.false. .and. s% star_mass_h1 < 0.35d0) then
       ! stop when star hydrogen mass drops to specified level
       extras_check_model = terminate
       write(*, *) 'have reached desired hydrogen mass'
       return
    end if

    ! ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(11) .eqv. .false.) .and. (s% r(1)/Rsun >= 19.9)) then
       ! absolute and relative tolerance
       atol = 0.5d-2
       rtol = 0.5d-2
       error = abs(s%r(1)/Rsun - 20.0)/ &
            (atol+rtol*max(s%r(1)/Rsun, 20.0))
       if (error > 1.0) then
          extras_check_model = retry
       end if
    end if

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(12) .eqv. .false.) .and. (s% r(1)/Rsun >= 99.0)) then
       ! absolute and relative tolerance
       atol = 1d-2
       rtol = 1d-2
       error = abs(s%r(1)/Rsun - 100.0)/ &
            (atol+rtol*max(s%r(1)/Rsun, 100.0))
       if (error > 1.0) then
          extras_check_model = retry
       end if
    end if

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(13) .eqv. .false.) .and. (s% r(1)/Rsun >= 199.0)) then
       ! absolute and relative tolerance
       atol = 1.5d-2
       rtol = 1.5d-2
       error = abs(s%r(1)/Rsun - 200.0)/ &
            (atol+rtol*max(s%r(1)/Rsun, 200.0))
       if (error > 1.0) then
          extras_check_model = retry
          print*, "error", error
       end if
    end if

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(14) .eqv. .false.) .and. (s% r(1)/Rsun >= 299.0)) then
       ! absolute and relative tolerance
       atol = 1.5d-2
       rtol = 1.5d-2
       error = abs(s%r(1)/Rsun - 300.0)/ &
            (atol+rtol*max(s% r(1)/Rsun, 300.0))
       if (error > 1.0) then
          extras_check_model = retry
       end if
    end if

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(15) .eqv. .false.) .and. (s% r(1)/Rsun >= 499.0)) then
       ! absolute and relative tolerance
       atol = 1d-2
       rtol = 1d-2
       error = abs(s%r(1)/Rsun - 500.0)/ &
            (atol+rtol*max(s%r(1)/Rsun, 500.0))
       if (error > 1.0) then
          extras_check_model = retry
       end if
    end if

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(16) .eqv. .false.) .and. (s% r(1)/Rsun >= 999.0)) then
       ! absolute and relative tolerance
       atol = 1.5d-2
       rtol = 1.5d-2
       error = abs(s%r(1)/Rsun - 1000.0)/ &
            (atol+rtol*max(s%r(1)/Rsun, 1000.0))
       if (error > 1.0) then
          extras_check_model = retry
       end if
    end if



    ! by default, indicate where (in the code) MESA terminated
    if (extras_check_model == terminate) s% termination_code = t_extras_check_model
  end function extras_check_model


  integer function how_many_extra_history_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_columns = 0
  end function how_many_extra_history_columns


  subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! note: do NOT add the extras names to history_columns.list
    ! the history_columns.list is only for the built-in history column options.
    ! it must not include the new column names you are adding here.


  end subroutine data_for_extra_history_columns


  integer function how_many_extra_profile_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_columns = 0
  end function how_many_extra_profile_columns


  subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
    integer, intent(in) :: id, n, nz
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: k
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! note: do NOT add the extra names to profile_columns.list
    ! the profile_columns.list is only for the built-in profile column options.
    ! it must not include the new column names you are adding here.

    ! here is an example for adding a profile column
    !if (n /= 1) stop 'data_for_extra_profile_columns'
    !names(1) = 'beta'
    !do k = 1, nz
    !   vals(k,1) = s% Pgas(k)/s% P(k)
    !end do

  end subroutine data_for_extra_profile_columns


  integer function how_many_extra_history_header_items(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_header_items = 0
  end function how_many_extra_history_header_items


  subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return

    ! here is an example for adding an extra history header item
    ! also set how_many_extra_history_header_items
    ! names(1) = 'mixing_length_alpha'
    ! vals(1) = s% mixing_length_alpha

  end subroutine data_for_extra_history_header_items


  integer function how_many_extra_profile_header_items(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_header_items = 0
  end function how_many_extra_profile_header_items


  subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return

    ! here is an example for adding an extra profile header item
    ! also set how_many_extra_profile_header_items
    ! names(1) = 'mixing_length_alpha'
    ! vals(1) = s% mixing_length_alpha

  end subroutine data_for_extra_profile_header_items


  ! returns either keep_going or terminate.
  ! note: cannot request retry; extras_check_model can do that.
  integer function extras_finish_step(id)
    integer, intent(in) :: id
    real(dp) :: TAMS_h1_treshold
    character (len=strlen) :: fname
    integer :: ierr
    integer :: k
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_finish_step = keep_going

    ! stop if Lnuc == L_phot+L_neu
    ! if (abs(s% L_phot + s% power_nonnuc_neutrinos - s% L_nuc_burn_total) <= 1d-9) then
    !    print *, "reached relax state"
    !    extras_finish_step = terminate
    ! end if

    if (s%lxtra(11) .eqv. .false.) then
       ! save profile for R=10Rsun
       if (s% r(1)/Rsun >= 20) then
          s% lxtra(11) = .true.
          write(fname, fmt="(a11)") '20Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a10)") '20Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a6)") '20Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(12) .eqv. .false.) then
       ! save profile for R=100Rsun
       if (s% r(1)/Rsun >= 100) then
          s% lxtra(12) = .true.
          write(fname, fmt="(a12)") '100Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '100Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '100Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(13) .eqv. .false.) then
       ! save profile for R=200Rsun
       if (s% r(1)/Rsun >= 200) then
          s% lxtra(13) = .true.
          write(fname, fmt="(a12)") '200Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '200Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '200Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(14) .eqv. .false.) then
       ! save profile for R=300Rsun
       if (s% r(1)/Rsun >= 300) then
          s% lxtra(14) = .true.
          write(fname, fmt="(a12)") '300Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '300Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '300Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(15) .eqv. .false.) then
       ! save profile for R=500Rsun
       if (s% r(1)/Rsun >= 500) then
          s% lxtra(15) = .true.
          write(fname, fmt="(a12)") '500Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '500Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '500Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(16) .eqv. .false.) then
       ! save profile for R=1000Rsun
       if (s% r(1)/Rsun >= 1000) then
          s% lxtra(16) = .true.
          write(fname, fmt="(a13)") '1000Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1000Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1000Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    TAMS_h1_treshold = 1d-2

    if ((s% center_h1 < TAMS_h1_treshold) .and. (s% center_he4 < 1.0d-4) .and. (s% center_c12 < 2.0d-2)) then
       write(*,'(g0)') "Single star depleted carbon, terminating from run_star_extras"
       extras_finish_step = terminate
    endif



    if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
  end function extras_finish_step


  subroutine extras_after_evolve(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
  end subroutine extras_after_evolve


end module run_star_extras
