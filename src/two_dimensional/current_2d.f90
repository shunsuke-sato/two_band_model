!---------------------------------------------------!
! Copyright (c) 2018 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine current_2d(jxy_intra,jxy_inter)
  use global_variables_2d
  implicit none
  real(8),intent(out) :: jxy_intra(2),jxy_inter(2)
  real(8) :: vk_xy(2,-NKx:NKx,-NKy:NKy)
  integer :: ikx,iyk
  real(8) ::
  complex(8) :: zfact

  call set_band_velocity_2d(vk_xy)

  jxy_intra = 0d0
  tmp = 0d0
  do ikz = -NKz,NKz
  do ikr = 1,NKr
    jxy_intra(1:2) = jxy_intra(1:2) + vk_xy(1:2,ikx,iky)*abs(zCt(2,ikx,iky))**2
    tmp = tmp + 2d0*real(conjg(zCt(1,ikx,iky))*zCt(2,ikx,iky))
  end do
  end do
  tmp = tmp*dkx*dky
  jxy_intra(1) = tmp*pix_vc
  jxy_intra(2) = tmp*piy_vc
  jxy_inter = jxy_inter*dkx*dky

  return
end subroutine current_2d
