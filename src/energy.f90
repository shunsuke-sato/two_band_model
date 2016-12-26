!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine energy(Eex)
  use global_variables
  implicit none
  real(8) :: Eex
  integer :: ikr,ikz

  Eex = 0d0
!$omp parallel
!$omp do private(ikz, ikr) reduction(+:Eex)
  do ikz = -NKz,NKz
  do ikr = 1,NKr
    Eex = Eex+ kz(ikz)*abs(zCt(2,ikr,ikz))**2*kr(ikr)*deps(ikr,ikz)
  end do
  end do
!$omp end parallel
  Eex = Eex*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 

  return
end subroutine energy
