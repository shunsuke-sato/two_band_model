!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
module global_variables
  implicit none
! mathematical parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zI = (0d0,1d0)
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0

! control parameter
  character(2) :: calc_mode = 'rf' ! rt or ft

! Material parameters
  real(8),parameter :: eps_g = 9d0/(2d0*Ry)  !1.52d0/(2d0*Ry) 
  real(8),parameter :: mass_r = 0.4d0 !1d0/(1d0/0.57d0+1d0/0.067d0)
! Discretization
  integer,parameter :: Nepsk_xy = 512, NKz = 512
  real(8),parameter :: eps_kxy_max = 1d0,eps_kz_max = 1d0
  real(8),parameter :: kz_max = sqrt(2d0*mass_r*eps_kz_max)
  real(8),parameter :: depskxy = eps_kxy_max/Nepsk_xy,dkz = kz_max/Nkz

! Time-propagation
  integer,parameter :: Nt = 10000
  real(8),parameter :: dt = 0.02d0
  complex(8) :: zDj(0:Nt)
  real(8) :: sin2w0_t_t0(0:Nt),cosw0_t_t0(0:Nt)

! Fourier transform
  integer,parameter :: Nw = 2000
  real(8),parameter :: wi = 0d0/(2d0*Ry), wf = 20d0/(2d0*Ry), dw = (wf-wi)/Nw
  complex(8) :: zDj_w(0:Nw)

! Laser parameters
  real(8),parameter :: E0 = 0d-2
  real(8),parameter :: omega0 = 1.55d0/(2d0*Ry)
  real(8),parameter :: T0 = 0d0



end module global_variables
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
program jdos_rt
  use global_variables
  implicit none

  select case(calc_mode)
  case('rt')
    call init_tfunction
    call calc_zDj
    call write_zDj
  case('ft')
    call read_zDj
    call calc_zDj_w
  case('rf')
    call init_tfunction
    call calc_zDj
    call write_zDj
    call calc_zDj_w
  end select

end program jdos_rt
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine init_tfunction
  use global_variables
  implicit none
  real(8) :: tt
  integer :: it

  do it = 0, Nt
    tt = dt*it + T0
    sin2w0_t_t0(it) = sin(2d0*omega0*tt) - sin(2d0*omega0*T0)
    cosw0_t_t0(it)  = cos(omega0*tt) - cos(omega0*T0)
  end do

end subroutine init_tfunction
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine calc_zDj
  use global_variables
  implicit none
  real(8) :: tt, eps_kxy, eps_kz, eps_kxyz, eps_tot, kz, Up
  real(8) :: theta1, theta2
  integer :: it, ikz, iepskxy
  real(8) :: ss1,ss2
  real(8) :: fact_xy(0:Nepsk_xy),fact_z(0:NKz)
  fact_xy = 1d0; fact_z = 1d0
  fact_xy(0) = 0.5d0; fact_z(0) = 0.5d0


  zDj = 0d0
  Up = E0**2/(4d0*mass_r*omega0**2)
  theta2 = E0**2/(8d0*mass_r*omega0**3)

!$omp parallel do private(iepskxy, ikz, eps_kxy, kz, &
!$omp & eps_kz, eps_kxyz, eps_tot,theta1,tt,ss1,ss2) &
!$omp collapse(2) reduction(+:zDj)
  do iepskxy = 0,Nepsk_xy
    do ikz = 0,Nkz
      eps_kxy = depskxy * iepskxy
      kz = dkz * ikz
      eps_kz = 0.5d0*kz**2/mass_r
      eps_kxyz = eps_g + eps_kxy + eps_kz
      eps_tot = eps_kxyz + Up
      theta1 = kz*E0/(mass_r*omega0**2)

      do it = 0,Nt
        tt = dt*it + T0
!== Start: full response ==
        ss1 = eps_tot*(tt-T0) - theta2*sin2w0_t_t0(it)
        ss2 = theta1*cosw0_t_t0(it)
        zDj(it) = zDj(it) + exp(-zI*ss1)*cos(ss2)*fact_xy(iepskxy)*fact_z(ikz)
!== End:   full response ==

!!== Start: full response (weak field limit) ==
!        zDj(it) = zDj(it) + (exp(-zI*eps_tot*(tt-T0)) &
!          +exp(-zI*eps_kxyz*(tt-T0))*( &
!          zI*theta2*sin2w0_t_t0(it) &
!          -0.5d0*(theta1*cosw0_t_t0(it))**2 &
!          ))*fact_xy(iepskxy)*fact_z(ikz)
!!== End:   full response (weak field limit) ==

      end do

    end do
  end do

  zDj = zDj*depskxy*dkz*4d0*pi*mass_r/(2d0*pi)**3


end subroutine calc_zDj
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine write_zDj
  use global_variables
  implicit none
  real(8) :: tt
  integer :: it

  open(20,file='zDj_t.out')
  write(20,"(A,2x,e26.16e3)")'# mass_r =',mass_r
  write(20,"(A,2x,e26.16e3)")'# eps_g  =',eps_g
  write(20,"(A,2x,e26.16e3)")'# omega0 =',omega0
  write(20,"(A,2x,e26.16e3)")'# T0     =',T0
  do it = 0,Nt
    write(20,"(999e26.16e3)")dt*it,zDj(it)
  end do
  close(20)

end subroutine write_zDj
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine read_zDj
  use global_variables
  implicit none
  real(8) :: tt
  integer :: it
  real(8) :: tmp1, tmp2, tmp3

  open(20,file='zDj_t.out')
  read(20,*)
  read(20,*)
  read(20,*)
  read(20,*)
  do it = 0,Nt
    read(20,*)tmp1,tmp2,tmp3
    zDj(it) = tmp2 + zI*tmp3
  end do
  close(20)

end subroutine read_zDj
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine calc_zDj_w
  use global_variables
  implicit none
  real(8) :: tt,ww,xx,ff
  integer :: it,iw
  complex(8) :: zs


  do iw = 0,Nw
    ww = wi + dw*iw
    zs = 0.5d0 * zDj(0)
    do it = 1,Nt
      tt = dt*it
      xx = dble(it)/dble(Nt)
      ff = 1d0 - 3d0*xx**2 + 2d0*xx**3
      zs = zs + zDj(it)*exp(zI*ww*tt)*ff
    end do
    zDj_w(iw) = zs*dt/pi

  end do


  open(20,file='zDj_w.out')
  write(20,"(A,2x,e26.16e3)")'# mass_r =',mass_r
  write(20,"(A,2x,e26.16e3)")'# eps_g  =',eps_g
  write(20,"(A,2x,e26.16e3)")'# omega0 =',omega0
  write(20,"(A,2x,e26.16e3)")'# T0     =',T0
  do iw = 0,Nw
    ww = wi + dw*iw
    ff = 0; if(ww > eps_g) ff = sqrt(2d0*mass_r)**3/(2d0*pi)**2*sqrt(ww-eps_g)
    write(20,"(999e26.16e3)")ww,real(zDj_w(iw)),ff
  end do
  close(20)

end subroutine calc_zDj_w
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
