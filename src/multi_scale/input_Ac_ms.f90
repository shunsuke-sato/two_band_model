!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine input_Ac_ms
  use global_variables_ms
  implicit none
  integer :: it,ix
  real(8) :: Xstart
  real(8) :: tt,xx

!  allocate(Act(0:Nt+1),Act_dt2(0:Nt+1),jtz(0:Nt+1),jtz_intra(0:Nt+1),jtz_inter(0:Nt+1))
  allocate(xn(Nx_L:Nx_R), Az(Nx_L:Nx_R), Az_new(Nx_L:Nx_R), Az_old(Nx_L:Nx_R), jz(Mx))
  Az_new = 0d0; Az = 0d0; Az_old = 0d0; jz = 0d0
  do ix = Nx_L,Nx_R
    xn(ix) = Hx*(dble(ix)-0.5d0)
  end do



  E0_1=5.338d-9*sqrt(Iwcm2_1)
  omega_1 = omega_ev_1/(2d0*Ry)
  tpulse_1 = tpulse_fs_1/0.02418d0
  E0_2=5.338d-9*sqrt(Iwcm2_2)
  omega_2 = omega_ev_2/(2d0*Ry)
  tpulse_2 = tpulse_fs_2/0.02418d0
  Tdelay = Tdelay_fs/0.02418d0

  Xstart=5*Hx

! Pump
  select case(envelope_1)
  case("cos4cos")
     do ix=NX_L,0
        xx=xn(ix)
        if(xx > -Xstart-tpulse_1*c_light .and. xx < -Xstart) then
           Az(ix)=-E0_1/omega_1*sin(pi*(xx+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**4 &
                *cos(omega_1*(xx+Xstart+tpulse_1*c_light)/c_light)
        endif
    
        xx=xx-dt*c_light
        if(xx > -Xstart-tpulse_1*c_light+dt*c_light .and. xx < -Xstart+dt*c_light) then
           Az_new(ix)=-E0_1/omega_1*sin(pi*(xx+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**4 &
                *cos(omega_1*(xx+Xstart+tpulse_1*c_light)/c_light)
        end if
     end do
  case("cos2cos")
     do ix=NX_L,0
        xx=xn(ix)
        if(xx > -Xstart-tpulse_1*c_light .and. xx < -Xstart) then
           Az(ix)=-E0_1/omega_1*sin(pi*(xx+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
                *cos(omega_1*(xx+Xstart+tpulse_1*c_light)/c_light)
        endif
    
        xx=xx-dt*c_light
        if(xx > -Xstart-tpulse_1*c_light+dt*c_light .and. xx < -Xstart+dt*c_light) then
           Az_new(ix)=-E0_1/omega_1*sin(pi*(xx+Xstart+tpulse_1*c_light)/(tpulse_1*c_light))**2 &
                *cos(omega_1*(xx+Xstart+tpulse_1*c_light)/c_light)
        end if
     end do
  case default
     stop "Invalid envelope_1"
  end select

! Probe
  select case(envelope_2)
  case("cos4cos")
     do ix=NX_L,0
        xx=xn(ix)
        if(xx > -Xstart-(tpulse_1*0.5d0+Tdelay)*c_light .and. &
             xx < -Xstart-(tpulse_1*0.5d0+Tdelay-tpulse_2)*c_light ) then
           Az(ix)=Az(ix)-E0_2/omega_2*sin(pi*(xx+Xstart+(tpulse_1*0.5d0+Tdelay)*c_light) &
                /(tpulse_2*c_light))**4 &
                *cos(omega_2*(xx+Xstart+(tpulse_1*0.5d0+Tdelay)*c_light)/c_light)
        endif

         
        xx=xx-dt*c_light
        if(xx > -Xstart-(tpulse_1*0.5d0+Tdelay)*c_light+dt*c_light &
             .and. xx < -Xstart-(tpulse_1*0.5d0+Tdelay-tpulse_2)*c_light+dt*c_light ) then
           Az_new(ix)=Az_new(ix)-E0_2/omega_2*sin(pi*(xx+Xstart+(tpulse_1*0.5d0+Tdelay)*c_light)&
                /(tpulse_2*c_light))**4 &
                *cos(omega_2*(xx+Xstart+(tpulse_1*0.5d0+Tdelay)*c_light)/c_light)
        end if
     end do
  case("cos2cos")
     do ix=NX_L,0
        xx=xn(ix)
        if(xx > -Xstart-(tpulse_1*0.5d0+Tdelay)*c_light .and. &
             xx < -Xstart-(tpulse_1*0.5d0+Tdelay-tpulse_2)*c_light ) then
           Az(ix)=Az(ix)-E0_2/omega_2*sin(pi*(xx+Xstart+(tpulse_1*0.5d0+Tdelay)*c_light) &
                /(tpulse_2*c_light))**2 &
                *cos(omega_2*(xx+Xstart+(tpulse_1*0.5d0+Tdelay)*c_light)/c_light)
        endif

         
        xx=xx-dt*c_light
        if(xx > -Xstart-(tpulse_1*0.5d0+Tdelay)*c_light+dt*c_light &
             .and. xx < -Xstart-(tpulse_1*0.5d0+Tdelay-tpulse_2)*c_light+dt*c_light ) then
           Az_new(ix)=Az_new(ix)-E0_2/omega_2*sin(pi*(xx+Xstart+(tpulse_1*0.5d0+Tdelay)*c_light)&
                /(tpulse_2*c_light))**2 &
                *cos(omega_2*(xx+Xstart+(tpulse_1*0.5d0+Tdelay)*c_light)/c_light)
        end if
     end do
  case default
     stop "Invalid envelope_2"
  end select

  return
end subroutine input_Ac_ms
