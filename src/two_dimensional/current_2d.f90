!---------------------------------------------------!
! Copyright (c) 2018 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine current_2d(jxyz_intra,jxyz_inter)
  use global_variables_2d
  implicit none
  real(8),intent(out) :: jxyz_intra(3),jxyz_inter(3)
  real(8) :: vk_xy(2,-NKx:NKx,-NKy:NKy)
  real(8) :: tmp
  integer :: ikx,iky
  complex(8) :: zfact

  call set_band_velocity_2d(vk_xy)

  jxyz_intra = 0d0
  tmp = 0d0
  do ikx = -NKx,NKx
  do iky = -Nky,NKz
    jxyz_intra(1:2) = jxyz_intra(1:2) + vk_xy(1:2,ikx,iky)*abs(zCt(2,ikx,iky))**2
    tmp = tmp + 2d0*real(conjg(zCt(1,ikx,iky))*zCt(2,ikx,iky))
  end do
  end do
  tmp = tmp*dkx*dky
  jxyz_intra(:) = jxyz_intra(:)*dkx*dky
  jxyz_inter(1) = pix_vc*tmp*dkx*dky
  jxyz_inter(2) = piy_vc*tmp*dkx*dky
  jxyz_inter(3) = piz_vc*tmp*dkx*dky

  return
end subroutine current_2d
