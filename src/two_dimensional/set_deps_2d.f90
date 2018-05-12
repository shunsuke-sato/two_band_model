!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine set_deps_2d
  use global_variables_2d
  implicit none
  integer :: ikx,iky

  select case(nband_type)
  case(N_PARABOLIC_BAND)
!$omp parallel do private(ikx, iky) collapse(2)
    do ikx = -NKx,NKx
      do iky = -NKy,NKy
        deps(ikx,iky) = eps_g + 0.5d0/mass_r*(kx(ikx)**2+ky(iky)**2)
      end do
    end do
  case(N_COS_BAND)
!$omp parallel do private(ikx, iky) collapse(2)
    do ikx = -NKx,NKx
      do iky = -NKy,NKy
        deps(ikx,iky) = eps_g - 0.5d0*band_width*(&
          cos(2d0*pi*kx(ikx)/kx_max)*cos(2d0*pi*ky(ikx)/ky_max) -1d0)
      end do
    end do
  case(N_COS4_BAND)
!$omp parallel do private(ikx, iky) collapse(2)
    do ikx = -NKx,NKx
      do iky = -NKy,NKy
        deps(ikx,iky) = eps_g +band_width*&
          (1d0 - cos(pi*kx(ikx)/kx_max)**2*cos(pi*ky(ikx)/ky_max))
      end do
    end do
  case default
    write(*,"(A,2x,A)")"Invalid nband_type",nband_type
    stop
  end select


  return
end subroutine set_deps_2d
