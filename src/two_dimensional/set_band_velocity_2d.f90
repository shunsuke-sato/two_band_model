!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine set_band_velocity_2d(vk_xy)
  use global_variables_2d
  implicit none
  real(8),intent(out) :: vk_xy(2,-NKx:NKx,-NKy:NKy)
  integer :: ikx,iky

  select case(nband_type)
  case(N_PARABOLIC_BAND)
    do ikx = -NKx,NKx
      vk_xy(1,ix,iy) = kx(ikx)/mass_r
    end do
    do iky = -NKy,NKy
      vk_xy(2,ix,iy) = ky(iky)/mass_r
    end do
  case default
    write(*,"(A,2x,A)")"Invalid nband_type",nband_type
    stop
  end select


end subroutine set_band_velocity_2d



