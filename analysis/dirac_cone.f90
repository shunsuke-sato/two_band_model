module global_variables
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0,1d0)
  real(8),parameter :: fs = 0.024189d0

! material parameters
  real(8),parameter :: tau_z = 1d0, delta_gap = 0d0, velocity = 1d0
  integer,parameter :: nkx = 256, nky = 256
  real(8),parameter :: kx_max = 2d0, ky_max = 2d0
  real(8) :: kx(nk),ky(nk)
  complex(8) :: zpsi(2,nkx,nky)

! time propagation
  real(8),parameter :: Tprop = 20d0/fs, dt = 0.08d0
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
    kx(ik) = -kx_max + 2d0*kx_max*dble(ikx-1)/dble(nkx-1)
  end do

  do iky = 1,nky
    ky(ik) = -ky_max + 2d0*ky_max*dble(iky-1)/dble(nky-1)
  end do

  if(delta_gap == 0d0)then ! Dirac cone
    do ikx = 1, nkx
      kxt = kx(ikx)
      do iky = 1,nky
        kyt = ky(iky)
        zpsi(1,ikx,iky) = -(tau_z*kxt-zI*kyt)/sqrt(2d0*(kx**2+ky**2))
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
  integer :: it

  do it = 0,nt-1
    call dt_evolve(it)
  end do


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

  ac = 0d0
  ac_dt2 = 0d0


end subroutine init_ac
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------!
