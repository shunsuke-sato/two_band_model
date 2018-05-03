!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine single_cell_2d
  use global_variables_2d
  implicit none


  if(Nprocs /= 1)call err_finalize("Parallelization is not supported &
    for single-cell calculation.")

  read(*,*)kx_max,ky_max
  read(*,*)NKx, NKy
  read(*,*)Nt,dt
  read(*,*)envelope_1
  read(*,*)Iwcm2_1,omega_ev_1,tpulse_fs_1,CEP_2pi_1
  read(*,*)dir_pol_1(1),dir_pol_2(1)
  read(*,*)envelope_2
  read(*,*)Iwcm2_2,omega_ev_2,tpulse_fs_2,CEP_2pi_2
  read(*,*)dir_pol_2(1),dir_pol_2(2)
  read(*,*)Tdelay_fs

  call preparation_2d

  call set_deps_2d

end subroutine single_cell_2d
