!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine input_Ac_2d
  use global_variables_2d
  implicit none
  integer :: it
  real(8) :: tt
  real(8) :: Es,Up,alpha

  allocate(Act_xyz(-1:Nt+2,3),jt_xyz(0:Nt+1,3))
  allocate(Act_pump(-1:Nt+2),Act_probe(-1:Nt+2) )

  E0_1=5.338d-9*sqrt(Iwcm2_1)
  omega_1 = omega_ev_1/(2d0*Ry)
  tpulse_1 = tpulse_fs_1/0.02418d0
  E0_2=5.338d-9*sqrt(Iwcm2_2)
  omega_2 = omega_ev_2/(2d0*Ry)
  tpulse_2 = tpulse_fs_2/0.02418d0
  Tdelay = Tdelay_fs/0.02418d0

!
  Es = piz_vc**2*E0_1**2/eps_g**3
  Up = E0_1**2/(4d0*mass_r*omega_1**2)
  alpha = (Es+Up)/eps_g

  write(*,"(A,2x,999e26.16e3)")"Stark shift, Es=",Es
  write(*,"(A,2x,999e26.16e3)")"Ponderomotive energy, Up=",Up
  write(*,"(A,2x,999e26.16e3)")"(Es+Up)/Eg=",(Es+Up)/eps_g
  write(*,"(A,2x,999e26.16e3)")"Es/Up=",Es/Up
  write(*,"(A,2x,999e26.16e3)")"E0_1=",E0_1
  write(*,"(A,2x,999e26.16e3)")"E0_2=",E0_2


  Act_xyz = 0d0; jt_xyz=0d0
  Act_pump = 0d0; Act_probe = 0d0

!Pump
  select case(envelope_1)
  case("cos4cos")
    do it = 0,Nt+1
      tt = dt*dble(it)
      if(abs(tt-0.5d0*tpulse_1) < 0.5d0*tpulse_1)then
        Act_pump(it) = -(E0_1/omega_1)*cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1)**4 &
          *sin(omega_1*(tt-0.5d0*tpulse_1)+2d0*pi*CEP_2pi_1)
      end if
    end do

  case("cos2cos")
    do it = 0,Nt+1
      tt = dt*dble(it)
      if(abs(tt-0.5d0*tpulse_1) < 0.5d0*tpulse_1)then
        Act_pump(it) = -(E0_1/omega_1)*cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1)**2 &
          *sin(omega_1*(tt-0.5d0*tpulse_1) +2d0*pi*CEP_2pi_1)
      end if
    end do

  case("cos6cos")
    do it = 0,Nt+1
      tt = dt*dble(it)
      if(abs(tt-0.5d0*tpulse_1) < 0.5d0*tpulse_1)then
        Act_pump(it) = -(E0_1/omega_1)*cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1)**6 &
          *sin(omega_1*(tt-0.5d0*tpulse_1)+2d0*pi*CEP_2pi_1)
      end if
    end do

  case default
    stop "Invalid envelope_1"
  end select

!Probe
  select case(envelope_2)
  case("cos4cos")
    do it = 0,Nt+1
      tt = dt*dble(it)

      if(abs(tt-0.5d0*tpulse_1-Tdelay) < 0.5d0*tpulse_2)then
        Act_probe(it) = -(E0_2/omega_2) &
          *cos(pi*(tt-0.5d0*tpulse_1-Tdelay)/tpulse_2)**4 &
          *sin(omega_2*(tt-0.5d0*tpulse_1-Tdelay)+2d0*pi*CEP_2pi_2)
      end if
    end do


  case("cos2cos")  
    do it = 0,Nt+1
      tt = dt*dble(it)

      if(abs(tt-0.5d0*tpulse_1-Tdelay) < 0.5d0*tpulse_2)then
        Act_probe(it) = -(E0_2/omega_2) &
          *cos(pi*(tt-0.5d0*tpulse_1-Tdelay)/tpulse_2)**2 &
          *sin(omega_2*(tt-0.5d0*tpulse_1-Tdelay)+2d0*pi*CEP_2pi_2)
      end if
    end do

  case("cos6cos")
    do it = 0,Nt+1
      tt = dt*dble(it)

      if(abs(tt-0.5d0*tpulse_1-Tdelay) < 0.5d0*tpulse_2)then
        Act_probe(it) =  -(E0_2/omega_2) &
          *cos(pi*(tt-0.5d0*tpulse_1-Tdelay)/tpulse_2)**6 &
          *sin(omega_2*(tt-0.5d0*tpulse_1-Tdelay)+2d0*pi*CEP_2pi_2)
      end if
    end do

  case default
    stop "Invalid envelope_2"
  end select

  Act_xyz(:,1) = dir_pol_1(1)*Act_pump(:) + dir_pol_2(1)*Act_probe(:)
  Act_xyz(:,2) = dir_pol_1(2)*Act_pump(:) + dir_pol_2(2)*Act_probe(:)
  Act_xyz(:,3) = dir_pol_1(3)*Act_pump(:) + dir_pol_2(3)*Act_probe(:)


  return
end subroutine input_Ac_2d
