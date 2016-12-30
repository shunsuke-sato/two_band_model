!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine input_Ac
  use global_variables
  implicit none
  integer :: it
  real(8) :: tt
  real(8) :: Es,Up,alpha

  allocate(Act(-1:Nt+1),jtz(0:Nt+1),jtz_intra(0:Nt+1),jtz_inter(0:Nt+1))

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


  Act = 0d0

!Pump
  select case(envelope_1)
  case("cos4cos")
    do it = 0,Nt+1
      tt = dt*dble(it)
      if(abs(tt-0.5d0*tpulse_1) < 0.5d0*tpulse_1)then
        Act(it) = -(E0_1/omega_1)*cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1)**4 &
          *sin(omega_1*(tt-0.5d0*tpulse_1))
      end if
    end do

    Act_dt2 = 0d0
    do it = 0,Nt+1
      tt = dt*dble(it) + 0.5d0*dt
      if(abs(tt-0.5d0*tpulse_1) < 0.5d0*tpulse_1)then
        Act_dt2(it) = -(E0_1/omega_1)*cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1)**4 &
          *sin(omega_1*(tt-0.5d0*tpulse_1))
      end if
    end do
  case("cos2cos")
    do it = 0,Nt+1
      tt = dt*dble(it)
      if(abs(tt-0.5d0*tpulse_1) < 0.5d0*tpulse_1)then
        Act(it) = -(E0_1/omega_1)*cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1)**2 &
          *sin(omega_1*(tt-0.5d0*tpulse_1))
      end if
    end do

    Act_dt2 = 0d0
    do it = 0,Nt+1
      tt = dt*dble(it) + 0.5d0*dt
      if(abs(tt-0.5d0*tpulse_1) < 0.5d0*tpulse_1)then
        Act_dt2(it) = -(E0_1/omega_1)*cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1)**2 &
          *sin(omega_1*(tt-0.5d0*tpulse_1))
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
        Act(it) = Act(it) -(E0_2/omega_2) &
          *cos(pi*(tt-0.5d0*tpulse_1-Tdelay)/tpulse_2)**4 &
          *sin(omega_2*(tt-0.5d0*tpulse_1-Tdelay))
      end if
    end do

    Act_dt2 = 0d0
    do it = 0,Nt+1
      tt = dt*dble(it) + 0.5d0*dt

      if(abs(tt-0.5d0*tpulse_1-Tdelay) < 0.5d0*tpulse_2)then
        Act_dt2(it) = Act_dt2(it) -(E0_2/omega_2) &
          *cos(pi*(tt-0.5d0*tpulse_1-Tdelay)/tpulse_2)**4 &
          *sin(omega_2*(tt-0.5d0*tpulse_1-Tdelay))
      end if
    end do
  case("cos2cos")  
    do it = 0,Nt+1
      tt = dt*dble(it)

      if(abs(tt-0.5d0*tpulse_1-Tdelay) < 0.5d0*tpulse_2)then
        Act(it) = Act(it) -(E0_2/omega_2) &
          *cos(pi*(tt-0.5d0*tpulse_1-Tdelay)/tpulse_2)**4 &
          *sin(omega_2*(tt-0.5d0*tpulse_1-Tdelay))
      end if
    end do

    Act_dt2 = 0d0
    do it = 0,Nt+1
      tt = dt*dble(it) + 0.5d0*dt

      if(abs(tt-0.5d0*tpulse_1-Tdelay) < 0.5d0*tpulse_2)then
        Act_dt2(it) = Act_dt2(it) -(E0_2/omega_2) &
          *cos(pi*(tt-0.5d0*tpulse_1-Tdelay)/tpulse_2)**4 &
          *sin(omega_2*(tt-0.5d0*tpulse_1-Tdelay))
      end if
    end do
  case default
    stop "Invalid envelope_1"
  end select


  return
end subroutine input_Ac
