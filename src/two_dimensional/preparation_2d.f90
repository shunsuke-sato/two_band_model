!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine preparation_2d
  use global_variables_2d
  implicit none
  integer :: ikx,iky

  allocate(zCt(2,-NKx:NKx,-Nky:NKy),deps(-NKx:NKx,-NKy:NKy))
  allocate(kx0(-NKx:NKx),kx(-NKx:NKx))
  allocate(ky0(-NKy:NKy),ky(-NKy:NKy))

  zCt = 0d0; zCt(1,:,:) = 1d0
  deps = 0d0

  dkx = kx_max/dble(NKx)
  dky = ky_max/dble(NKy)

  do ikx = -NKx,NKx
    kx0(ikx) = dkx*dble(ikx)
  end do
  kx = kx0

  do iky = -NKy,NKy
    ky0(iky) = dky*dble(iky)
  end do
  ky = ky0

  call set_deps_2d

  return
end subroutine preparation_2d
