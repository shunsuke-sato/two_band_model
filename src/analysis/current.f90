!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine current(jz_intra,jz_inter)
  use global_variables
  implicit none
  integer :: ikr,ikz
  real(8) :: jz_intra,jz_inter
  complex(8) :: zfact

  call set_ddeps_dkz

  jz_intra = 0d0
  jz_inter = 0d0
!$omp parallel
!$omp do private(ikz, ikr) reduction(+:jz_intra,jz_inter)
  do ikz = -NKz,NKz
  do ikr = 1,NKr
    jz_intra = jz_intra + ddeps_dkz(ikr,ikz)*abs(zCt(2,ikr,ikz))**2*kr(ikr)
    jz_inter = jz_inter + 2d0*real(conjg(zCt(1,ikr,ikz))*zCt(2,ikr,ikz))*kr(ikr)
  end do
  end do
!$omp end parallel
  jz_intra = jz_intra*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
  jz_inter = jz_inter*piz_vc*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 

  return
end subroutine current
