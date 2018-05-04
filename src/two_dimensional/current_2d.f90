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
  real(8),parameter :: ss = 1d0/(2d0*pi)**2
  real(8) :: vk_xy(2,-NKx:NKx,-NKy:NKy)
  real(8) :: tmp,jx,jy
  integer :: ikx,iky
  complex(8) :: zfact

  call set_band_velocity_2d(vk_xy)

  jx = 0d0; jy = 0d0
  tmp = 0d0
!$omp parallel do private(ikx, iky) reduction(+:tmp,jx,jy) collapse(2)
  do ikx = -NKx,NKx
  do iky = -Nky,NKy
    jx = jx + vk_xy(1,ikx,iky)*abs(zCt(2,ikx,iky))**2
    jy = jy + vk_xy(2,ikx,iky)*abs(zCt(2,ikx,iky))**2
    tmp = tmp + 2d0*real(conjg(zCt(1,ikx,iky))*zCt(2,ikx,iky))
  end do
  end do

  jxyz_intra(1) = 2d0*jx*dkx*dky*ss
  jxyz_intra(2) = 2d0*jy*dkx*dky*ss
  jxyz_intra(3) = 0d0
  jxyz_inter(1) = 2d0*pix_vc*tmp*dkx*dky*ss
  jxyz_inter(2) = 2d0*piy_vc*tmp*dkx*dky*ss
  jxyz_inter(3) = 2d0*piz_vc*tmp*dkx*dky*ss

  return
end subroutine current_2d
