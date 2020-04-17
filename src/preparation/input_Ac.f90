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

  allocate(Act(-1:Nt+2),jtz(0:Nt+1),jtz_intra(0:Nt+1),jtz_inter(0:Nt+1))
  allocate(Act_pump(-1:Nt+2),Act_probe(-1:Nt+2) )

!  E0_1=5.338d-9*sqrt(Iwcm2_1)
  E0_1=E0_1_V_m*a_B*1d-10/(2d0*Ry)
  omega_1 = omega_ev_1/(2d0*Ry)
  tpulse_1 = tpulse_fs_1/0.02418d0
!  E0_2=5.338d-9*sqrt(Iwcm2_2)
  E0_2=E0_2_V_m*a_B*1d-10/(2d0*Ry)
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


  Act = 0d0; jtz=0d0; jtz_intra = 0d0; jtz_inter = 0d0
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

  Act = Act_pump + Act_probe

  if(if_pure_intraband_fields_exist)then
    call input_Ac_intra
  end if

  return

  contains
    subroutine input_Ac_intra
      use global_variables
      implicit none
      integer :: it
      real(8) :: tt, xx

      E0_static = E0_static_V_AA*a_B/ev
      T_duration_static = tpulse_1

      allocate(Act_intra(-1:Nt+2))

     
      do it = -1, Nt+2
        tt = dt*dble(it)
        xx = tt - 0.5d0*T_duration_static
        Act_intra(it) = -E0_static*xx
      end do



    end subroutine input_Ac_intra

end subroutine input_Ac
!--------------------------------------------------------------------------------
