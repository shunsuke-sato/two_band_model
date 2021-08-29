!---------------------------------------------------!
! Copyright (c) 2020 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine add_impulsive_distortion
  use global_variables
  implicit none
  integer :: ikr,ikz,ik
  complex(8) :: alpha

!$omp parallel do private(ikz)
  do ikz = -NKz,NKz
    kz(ikz) = kz0(ikz) + Act(0)
  end do

  call set_deps

!$omp parallel do private(ikz,ikr,alpha)
  do ikz = -NKz,NKz
    do ikr = 1,NKr
      alpha = piz_vc*kick_impulsive/deps(ikr,ikz)      
      zCt(1,ikr,ikz) = 1d0/sqrt(1d0 + alpha**2)
      zCt(2,ikr,ikz) = -alpha/sqrt(1d0 + alpha**2)
    end do
  end do
 

end subroutine add_impulsive_distortion
