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
  real(8) :: dipx, dipy, mu_y

  select case(nband_type)
  case(N_PARABOLIC_BAND)
!$omp parallel do private(ikx)
    do ikx = -NKx,NKx
      vk_xy(1,ikx,:) = kx(ikx)/mass_r
    end do
!$omp parallel do private(iky)
    do iky = -NKy,NKy
      vk_xy(2,:,iky) = ky(iky)/mass_r
    end do
  case(N_COS_BAND)
    dipx = 2d0*pi/kx_max
    dipy = 2d0*pi/ky_max
!$omp parallel do private(ikx, iky) collapse(2)
    do ikx = -NKx,NKx
      do iky = -NKy,NKy
        vk_xy(1,ikx,iky) = 0.5d0*band_width*dipx*sin(dipx*kx(ikx))*cos(dipy*ky(iky))
        vk_xy(2,ikx,iky) = 0.5d0*band_width*dipy*cos(dipx*kx(ikx))*sin(dipy*ky(iky))
      end do
    end do
  case(N_COS4_BAND)
    dipx = pi/kx_max
    dipy = pi/ky_max
!$omp parallel do private(ikx, iky) collapse(2)
    do ikx = -NKx,NKx
      do iky = -NKy,NKy
        vk_xy(1,ikx,iky) = 4d0*band_width*dipx&
          *cos(dipx*kx(ikx))**3*sin(dipx*kx(ikx))*cos(dipy*ky(iky))**4
        vk_xy(2,ikx,iky) = 4d0*band_width*dipy&
          *cos(dipy*ky(iky))**3*sin(dipy*ky(iky))*cos(dipx*kx(ikx))**4
      end do
    end do
  case(N_HBN_BAND)
    dipx = 2d0*pi/kx_max
    mu_y = 2d0/(dipx**2*band_width)
!$omp parallel do private(ikx, iky) collapse(2)
    do ikx = -NKx,NKx
      do iky = -NKy,NKy
        vk_xy(1,ikx,iky) = dipx*0.5d0*band_width*sin(dipx*kx(ikx))
        vk_xy(2,ikx,iky) = ky(iky)/mu_y
      end do
    end do
  case default
    write(*,"(A,2x,A)")"Invalid nband_type",nband_type
    stop
  end select

end subroutine set_band_velocity_2d



