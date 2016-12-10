!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine prop_Ac_ms
  use global_variables_ms
  implicit none
  integer :: ix
  real(8) :: jz_intra,jz_inter
  real(8) :: jz_t(Nx_L:Nx_R)
  complex(8) :: zfact

  jz_t = 0d0
  jz_t(1:Mx) = jz(1:Mx)

  Az_old = Az
  Az = Az_new
  Az_new = 0d0

  do ix = Nx_L+1,Nx_R-1
    Az_new(ix)=2d0*Az(ix)-Az_old(ix) &
      +(c_light*dt/Hx)**2*(Az(ix+1)-2d0*Az(ix)+Az(ix-1)) &
      -4d0*pi*dt**2*jz_t(ix)
  end do

  return
end subroutine prop_Ac_ms
