!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine input_Ac_ms
  use global_variables_ms
  implicit none
  integer :: it
  real(8) :: tt

!  allocate(Act(0:Nt+1),Act_dt2(0:Nt+1),jtz(0:Nt+1),jtz_intra(0:Nt+1),jtz_inter(0:Nt+1))

  E0_1=5.338d-9*sqrt(Iwcm2_1)
  omega_1 = omega_ev_1/(2d0*Ry)
  tpulse_1 = tpulse_fs_1/0.02418d0
  E0_2=5.338d-9*sqrt(Iwcm2_2)
  omega_2 = omega_ev_2/(2d0*Ry)
  tpulse_2 = tpulse_fs_2/0.02418d0
  Tdelay = Tdelay_fs/0.02418d0


  

  return
end subroutine input_Ac_ms
