!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine set_deps
  use global_variables
  implicit none
  integer :: ikr,ikz

  select case(nband_type)
  case(N_PARABOLIC_BAND)
!$omp parallel do private(ikz, ikr)
    do ikz = -NKz,NKz
      do ikr = 1,NKr
        deps(ikr,ikz) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
      end do
    end do
  case(N_NONPARABOLIC_BAND)
!$omp parallel do private(ikz, ikr)
    do ikz = -NKz,NKz
      do ikr = 1,NKr
        deps(ikr,ikz) = eps_g*sqrt(1d0 + (kr(ikr)**2+kz(ikz)**2)/(mass_r*eps_g)) 
      end do
    end do
  case default
    write(*,"(A,2x,A)")"Invalid nband_type",nband_type
    stop
  end select


  return
end subroutine set_deps
