!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
program main
  implicit none
  character(50) :: calc_mode

  read(*,*)calc_mode

  select case(calc_mode)
  case("single-cell")
    call single_cell
!  case("multi-scale")
!    call multi_scale
  case default
    write(*,"(A,2x,A)")"Invalid calc_mode",calc_mode
    stop
  end select
      

end program main
