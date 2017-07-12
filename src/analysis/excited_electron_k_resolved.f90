!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine excited_electron_k_resolved(nex_kz)
  use global_variables
  implicit none
  real(8),intent(out) :: nex_kz(-NKz:NKz)
  real(8) :: ss
  integer :: ikz,ikr

!$omp parallel do private(ikz,ss,ikr)
  do ikz = -NKz,NKz

    ss = 0d0
    do ikr = 1,NKr
      ss = ss + abs(zCt(2,ikr,ikz))**2*kr(ikr)
    end do

    nex_kz(ikz) = ss *2d0/((2d0*pi)**3)*(2d0*pi*dkr) 

  end do

end subroutine excited_electron_k_resolved
