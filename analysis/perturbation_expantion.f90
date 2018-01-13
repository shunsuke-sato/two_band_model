!---------------------------------------------------!
! Copyright (c) 2018 Shunsuke A. Sato               !
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
!  real(8),parameter :: eps_g = (40.73d0+1.52d0)/(2d0*Ry)  !1.52d0/(2d0*Ry) 
  real(8),parameter :: eps_g = 42d0/(2d0*Ry)  !1.52d0/(2d0*Ry) 
  real(8),parameter :: mass_r = 0.067d0
! Discretization
  integer,parameter :: Nepsk_xy = 128, NKz = 128
  real(8),parameter :: eps_kxy_max = 0.5d0,eps_kz_max = 0.5d0
  real(8),parameter :: kz_max = sqrt(2d0*mass_r*eps_kz_max)
  real(8),parameter :: depskxy = eps_kxy_max/Nepsk_xy,dkz = kz_max/Nkz
  real(8),parameter :: gamma = (0.1d0/(2d0*ry))

! Time-propagation
  integer,parameter :: Nt = 10000
  real(8),parameter :: dt = 0.02d0
  complex(8) :: zDj(0:Nt)
  real(8) :: sin2w0_t_t0(0:Nt),cosw0_t_t0(0:Nt),pow_cosw0_t_t0(0:Nt)

! Fourier transform
  integer,parameter :: Nw = 400
  real(8),parameter :: wi = 36d0/(2d0*Ry), wf = 50d0/(2d0*Ry), dw = (wf-wi)/Nw
  complex(8) :: zDj_w(0:Nw)

! Laser parameters
  real(8),parameter :: E0 = 1d-2
  real(8),parameter :: omega0 = 1.55d0/(2d0*Ry)
  real(8),parameter :: T0 = 0d0



end module global_variables
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
program main
  use global_variables
  implicit none
  integer :: it, ikz, iepskxy, iw, nt_delay
  real(8) :: Up, theta2,pvc
  real(8) :: eps_kxy, eps_kz, eps_kxyz, kz
  real(8) :: ww
  real(8) :: tt
  complex(8)  :: zEw0, zEw_2p, zEw_2m
  complex(8) :: zJ3w(0:Nw),zsigma, zc3

  open(20,file='sigma_tr.out')

  pvc = 1d0
  Up = E0**2/(4d0*mass_r*omega0**2)

  tt = 0d0
  nt_delay = 16

  do it = 0, nt_delay
    tt = it*(2d0*pi/omega0)/dble(nt_delay)
    

  zj3w = 0d0

  do iepskxy = 0, Nepsk_xy
    write(*,*)iepskxy
    do ikz = 0, Nkz
      eps_kxy = depskxy * iepskxy
      kz = dkz * ikz
      eps_kz = 0.5d0*kz**2/mass_r
      eps_kxyz = eps_g + eps_kxy + eps_kz
      theta2  = ( (kz*E0)/(mass_r*omega0) )**2

      do iw = 0,Nw
        ww = wi + dw*dble(iw)
        zEw0 = exp(zI*ww*tt)
        zEw_2p = exp(zI*(ww+2d0*omega0)*tt)
        zEw_2m = exp(zI*(ww-2d0*omega0)*tt)

!! positive frequency 
! fundamental frequency
        zc3 = 0d0
        zc3 = zc3 + Up*zEw0/(ww-eps_kxyz+zi*gamma)
        zc3 = zc3 + (1d0/(2d0*eps_kxyz**2)*theta2-Up/eps_kxyz)*zEw0
        zc3 = zc3 + theta2/4d0*(1d0/(ww-eps_kxyz+omega0+zi*gamma))* &
          (1d0/(ww-eps_kxyz+zi*gamma)-1d0/eps_kxyz)*zEw0
        zc3 = zc3 + theta2/4d0*(1d0/(ww-eps_kxyz-omega0+zi*gamma))* &
          (1d0/(ww-eps_kxyz+zi*gamma)-1d0/eps_kxyz)*zEw0

! + 2*omega0
        zc3 = zc3 -0.25d0*theta2/(ww-eps_kxyz+omega0+zi*gamma)* &
          (1d0/(ww-eps_kxyz+2d0*omega0+zI*gamma)-1d0/eps_kxyz)*zEw_2p
        zc3 = zc3 -0.5d0*Up/(ww-eps_kxyz+2d0*omega0+zi*gamma)*zEw_2p
        zc3 = zc3 -0.5d0*(0.5d0*theta2/eps_kxyz**2 - Up/eps_kxyz)* zEw_2p

! - 2*omega0
        zc3 = zc3 -0.25d0*theta2/(ww-eps_kxyz-omega0+zi*gamma)* &
          (1d0/(ww-eps_kxyz-2d0*omega0+zI*gamma)-1d0/eps_kxyz)*zEw_2m
        zc3 = zc3 -0.5d0*Up/(ww-eps_kxyz-2d0*omega0+zi*gamma)*zEw_2m
        zc3 = zc3 -0.5d0*(0.5d0*theta2/eps_kxyz**2 - Up/eps_kxyz)* zEw_2m
!!
        zc3 = zc3*(-zI*pvc/eps_kxyz)/(ww-eps_kxyz+zI*gamma)
        zJ3w(iw) =  zJ3w(iw) + pvc*zc3

        cycle ! debug
!! negative frequency
        ww = -ww
! fundamental frequency
        zc3 = 0d0
        zc3 = zc3 + Up*zEw0/(ww-eps_kxyz+zi*gamma)
        zc3 = zc3 + (1d0/(2d0*eps_kxyz**2)*theta2-Up/eps_kxyz)*zEw0
        zc3 = zc3 + theta2/4d0*(1d0/(ww-eps_kxyz+omega0+zi*gamma))* &
          (1d0/(ww-eps_kxyz+zi*gamma)-1d0/eps_kxyz)*zEw0
        zc3 = zc3 + theta2/4d0*(1d0/(ww-eps_kxyz-omega0+zi*gamma))* &
          (1d0/(ww-eps_kxyz+zi*gamma)-1d0/eps_kxyz)*zEw0

! + 2*omega0
        zc3 = zc3 -0.25d0*theta2/(ww-eps_kxyz+omega0+zi*gamma)* &
          (1d0/(ww-eps_kxyz+2d0*omega0+zI*gamma)-1d0/eps_kxyz)*zEw_2p
        zc3 = zc3 -0.5d0*Up/(ww-eps_kxyz+2d0*omega0+zi*gamma)*zEw_2p
        zc3 = zc3 -0.5d0*(0.5d0*theta2/eps_kxyz**2 - Up/eps_kxyz)* zEw_2p

! - 2*omega0
        zc3 = zc3 -0.25d0*theta2/(ww-eps_kxyz-omega0+zi*gamma)* &
          (1d0/(ww-eps_kxyz-2d0*omega0+zI*gamma)-1d0/eps_kxyz)*zEw_2m
        zc3 = zc3 -0.5d0*Up/(ww-eps_kxyz-2d0*omega0+zi*gamma)*zEw_2m
        zc3 = zc3 -0.5d0*(0.5d0*theta2/eps_kxyz**2 - Up/eps_kxyz)* zEw_2m

        zc3 = zc3*(-zI*pvc/eps_kxyz)/(ww-eps_kxyz+zI*gamma)
        zJ3w(iw) =  zJ3w(iw) +  conjg(pvc*zc3)


        
      end do


    end do
  end do

  zJ3w = zJ3w*depskxy*dkz*4d0*pi*mass_r/(2d0*pi)**3
  do iw = 0,nw
    ww = wi + dw*dble(iw)
    zEw0 = exp(zI*ww*tt)
    zsigma = zJ3w(iw)/zEw0
    write(20,"(999e26.16e3)")tt,ww,zsigma,zJ3w(iw)

  end do

  write(20,*)

end do

  close(20)

end program main
