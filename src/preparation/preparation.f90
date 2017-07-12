!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine preparation
  use global_variables
  implicit none
  integer :: ikr,ikz,ik


  allocate(zCt(2,NKr,-Nkz:NKz),deps(NKr,-NKz:NKz),deps_int(NKr,-NKz:NKz))
  allocate(ddeps_dkz(NKr,-NKz:NKz))
  allocate(kz0(-NKz:NKz),kz(-NKz:NKz),kr(NKr))
  zCt = 0d0; zCt(1,:,:) = 1d0
  deps = 0d0
  deps_int = 0d0
  

  dkr = kr_max/dble(NKr)
  dkz = kz_max/dble(NKz)

  do ikz = -NKz,NKz
    kz0(ikz) = dkz*dble(ikz)
  end do
  kz = kz0

  do ikr = 1,NKr
    kr(ikr) = dkr*dble(ikr)
  end do



  kz = kz0
  call set_deps

  return
end subroutine preparation
