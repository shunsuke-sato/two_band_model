module global_variables
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0,1d0)
  real(8),parameter :: fs = 0.024189d0
  real(8),parameter :: ev = 27.2114d0
  real(8),parameter :: b_a = 0.529d0
  real(8),parameter :: clight = 137.035999139d0

! material parameters
  real(8),parameter :: tau_z = 1d0, delta_gap = 0d0
! Fermi velocity of graphene PRL 101, 226405 (2008).
  real(8),parameter :: velocity = clight*1.12d6/299792458d0
  integer,parameter :: nkx = 256, nky = 256
  real(8),parameter :: kx_max = 2d0, ky_max = 2d0
  real(8) :: kx(nkx),ky(nky)
  real(8) :: dkx,dky
  complex(8) :: zpsi(2,nkx,nky)

! time propagation
  real(8),parameter :: Tprop = 100d0/fs, dt = 0.08d0
  integer,parameter :: nt = aint(Tprop/dt)+1

  real(8) :: jt(2,0:nt), ac(2,0:nt), ac_dt2(2,0:nt)

end module global_variables
!----------------------------------------------------------------------------------------!
program main
  use global_variables
  implicit none

  call init
  call init_ac

  call time_propagation

end program main
!----------------------------------------------------------------------------------------!
subroutine init
  use global_variables
  implicit none
  integer :: ikx,iky
  real(8) :: kxt, kyt

  do ikx = 1,nkx
    kx(ikx) = -kx_max + 2d0*kx_max*dble(ikx-1)/dble(nkx-1)
  end do
  dkx = 2d0*kx_max/dble(nkx-1)

  do iky = 1,nky
    ky(iky) = -ky_max + 2d0*ky_max*dble(iky-1)/dble(nky-1)
  end do
  dky = 2d0*ky_max/dble(nky-1)

  if(delta_gap == 0d0)then ! Dirac cone
    do ikx = 1, nkx
      kxt = kx(ikx)
      do iky = 1,nky
        kyt = ky(iky)
        zpsi(1,ikx,iky) = -(tau_z*kxt-zI*kyt)/sqrt(2d0*(kxt**2+kyt**2))
        zpsi(2,ikx,iky) = 1d0/sqrt(2d0)
      end do
    end do
  else
    stop 'Kane band is not implemented yet.'
  end if


end subroutine init
!----------------------------------------------------------------------------------------!
subroutine time_propagation
  use global_variables
  implicit none
  real(8) :: jxy(2)
  integer :: it

  open(21,file="current.out")
  call current(jxy(:),0)
  write(21,"(999e26.16e3)")0d0,ac(:,0),jxy
  do it = 0,nt-1
    if(mod(it,max(100,nt/100)) == 0)write(*,*)"it=",it,nt
    call dt_evolve(it)
    call current(jxy(:),it)
    write(21,"(999e26.16e3)")dt*(it+1),ac(:,it+1),jxy
  end do
  close(21)


end subroutine time_propagation
!----------------------------------------------------------------------------------------!
subroutine dt_evolve(it)
  use global_variables
  integer,intent(in) :: it
  integer :: ikx, iky
  real(8) :: kxt, kyt
  complex(8) :: zpsi_t(2,0:4)

  do ikx = 1,nkx
    do iky = 1,nky

      kxt = kx(ikx) + ac(1,it)
      kyt = kx(iky) + ac(2,it)
! RK1
      zpsi_t(:,0) = zpsi(:,ikx,iky)
      zpsi_t(1,1) = 0.5d0*delta_gap*zpsi_t(1,0) &
        + velocity*(tau_z*kxt - zI*kyt)*zpsi_t(2,0) 
      zpsi_t(2,1) = velocity*(tau_z*kxt + zI*kyt)*zpsi_t(1,0)  &
        - 0.5d0*delta_gap*zpsi_t(2,0)

      kxt = kx(ikx) + ac_dt2(1,it)
      kyt = kx(iky) + ac_dt2(2,it)
! RK2
      zpsi_t(:,0) = zpsi(:,ikx,iky) - zI*0.5d0*dt*zpsi_t(:,1)

      zpsi_t(1,2) = 0.5d0*delta_gap*zpsi_t(1,0) &
        + velocity*(tau_z*kxt - zI*kyt)*zpsi_t(2,0) 
      zpsi_t(2,2) = velocity*(tau_z*kxt + zI*kyt)*zpsi_t(1,0)  &
        - 0.5d0*delta_gap*zpsi_t(2,0)

! RK3
      zpsi_t(:,0) = zpsi(:,ikx,iky) - zI*0.5d0*dt*zpsi_t(:,2)

      zpsi_t(1,3) = 0.5d0*delta_gap*zpsi_t(1,0) &
        + velocity*(tau_z*kxt - zI*kyt)*zpsi_t(2,0) 
      zpsi_t(2,3) = velocity*(tau_z*kxt + zI*kyt)*zpsi_t(1,0)  &
        - 0.5d0*delta_gap*zpsi_t(2,0)

      kxt = kx(ikx) + ac(1,it+1)
      kyt = kx(iky) + ac(2,it+1)

! RK4
      zpsi_t(:,0) = zpsi(:,ikx,iky) - zI*dt*zpsi_t(:,3)

      zpsi_t(1,4) = 0.5d0*delta_gap*zpsi_t(1,0) &
        + velocity*(tau_z*kxt - zI*kyt)*zpsi_t(2,0) 
      zpsi_t(2,4) = velocity*(tau_z*kxt + zI*kyt)*zpsi_t(1,0)  &
        - 0.5d0*delta_gap*zpsi_t(2,0)

      zpsi(:,ikx,iky) = zpsi(:,ikx,iky) -zI*dt/6d0*( &
        zpsi_t(:,1) + 2d0*zpsi_t(:,2) + 2d0*zpsi_t(:,3) + zpsi_t(:,4) )
      
    end do
  end do
  

end subroutine dt_evolve
!----------------------------------------------------------------------------------------!
subroutine init_ac
  use global_variables
  implicit none
  integer :: it
  real(8) :: tt,tt2
  real(8) :: E0,omega,tpulse
  real(8) :: E_SD, T_SD

  E0 = 1d7*b_a*1d-10/ev
  omega = 0.19074d0/ev ! 6.5 micron
  tpulse = 2d3/fs

! Source-drain
  E_SD = 1d-4
  T_SD = 10d0/fs

  ac = 0d0
  ac_dt2 = 0d0

  do it = 0,nt

    tt = dt*it
    tt2 = dt*it + 0.5d0*dt

! source-drain
    if(tt<T_SD)then
      ac(1,it) = ac(1,it) - E_SD*(tt/2d0 - T_SD/(2d0*pi)*sin(pi*tt/T_SD))
    else
      ac(1,it) = ac(1,it) - E_SD*(tt - T_SD/2d0)
    end if
    if(tt2<T_SD)then
      ac_dt2(1,it) = ac_dt2(1,it) - E_SD*(tt2/2d0 - T_SD/(2d0*pi)*sin(pi*tt2/T_SD))
    else
      ac_dt2(1,it) = ac_dt2(1,it) - E_SD*(tt2 - T_SD/2d0)
    end if


     if( abs(tt-0.5d0*tpulse)<0.5d0*tpulse )then
        ac(1,it) = ac(1,it) - E0/omega*cos(pi*(tt-0.5d0*tpulse)/tpulse)**2*sin(omega*(tt-0.5d0*tpulse))
     end if

     tt = dt*it + dt/2d0
     if( abs(tt-0.5d0*tpulse)<0.5d0*tpulse )then
        ac_dt2(1,it) = ac_dt2(1,it) - E0/omega*cos(pi*(tt-0.5d0*tpulse)/tpulse)**2*sin(omega*(tt-0.5d0*tpulse))
     end if

     
  end do


end subroutine init_ac
!----------------------------------------------------------------------------------------!
subroutine current(jxy,it)
  use global_variables
  implicit none
  real(8) :: jxy(2)
  integer :: it
  integer :: ikx, iky
  real(8) :: kxt, kyt
  real(8) :: jxt,jyt,jx0,jy0,xx
  complex(8) :: zs

  jxy = 0d0

  do ikx = 1,nkx
    do iky = 1,nky

       kxt = kx(ikx) + ac(1,it)
       kyt = kx(iky) + ac(2,it)
       
       jxt = real(zpsi(1,ikx,iky)*conjg(zpsi(2,ikx,iky)))
       jyt = real(zI*zpsi(1,ikx,iky)*conjg(zpsi(2,ikx,iky)))
       
       zs = -velocity*(tau_z*kxt-zI*kyt)&
            /(delta_gap*0.5d0&
            +sqrt(delta_gap**2*0.25d0+velocity**2*(kxt**2+kyt**2)))
       
       xx = abs(zs)**2
       jx0 =  2d0/(1d0+xx)*real(zs)
       jy0 = -2d0/(1d0+xx)*aimag(zs)
       
       jxy(1) = jxy(1) + jxt - jx0
       jxy(2) = jxy(2) + jyt - jy0
       
     end do
  end do

  jxy = jxy*dkx*dky

end subroutine current
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
