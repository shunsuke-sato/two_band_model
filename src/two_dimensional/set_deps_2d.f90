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
  case default
    write(*,"(A,2x,A)")"Invalid nband_type",nband_type
    stop
  end select


  return
end subroutine set_deps_2d
