!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine single_cell_2d
  use global_variables_2d
  implicit none
  integer :: it
  real(8) :: jxyz_intra(3), jxyz_inter(3)

  if(Nprocs /= 1)call err_finalize("Parallelization is not supported &
    for single-cell calculation.")

  read(*,*)kx_max,ky_max
  read(*,*)NKx, NKy
  read(*,*)Nt,dt
  read(*,*)envelope_1
  read(*,*)Iwcm2_1,omega_ev_1,tpulse_fs_1,CEP_2pi_1
  read(*,*)dir_pol_1(1:3)
  read(*,*)envelope_2
  read(*,*)Iwcm2_2,omega_ev_2,tpulse_fs_2,CEP_2pi_2
  read(*,*)dir_pol_2(1:3)
  read(*,*)Tdelay_fs

  call preparation_2d
  call input_Ac_2d

  do it = 0,Nt
    write(*,*)'it=',it,'/',Nt
    call dt_evolve_2d(it)
    call current_2d(jxyz_intra,jxyz_inter)
    jt_xy(it+1,1:2) = jxy_intra(:) + jxy_inter(:)
  end do

  open(21,file='Act_jt.out')
  do it = 0,Nt
    write(21,'(999e26.16e3)')dt*dble(it),Act_xy(it,1:2),jt_xy(it,1:2) &
      ,-0.5d0*(Act_xy(it+1,1:2)-Act_xy(it-1,1:2))/dt
  end do
  close(21)

end subroutine single_cell_2d
