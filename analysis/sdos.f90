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
  real(8),parameter :: a_B=0.529177d0,ev=27.2114d0


! Material parameters
  real(8),parameter :: eps_g = 3.2d0/ev  !1.52d0/(2d0*Ry) 
  real(8),parameter :: mass_r = 0.13d0 !1d0/(1d0/0.57d0+1d0/0.067d0)

  integer,parameter :: nband_type_parabolic = 0
  integer,parameter :: nband_type_Kane = 1
  integer,parameter :: nband_type = nband_type_parabolic


! Discretization
  integer,parameter :: NKr = 64, NKz = 64
  real(8),parameter :: kz_max = 0.1d0, kr_max = 0.1d0
  real(8),parameter :: dkz = kz_max/Nkz, dkr = kr_max/Nkr


! Fourier transform
  real(8),parameter :: gamma = 0.1d0/ev
  real(8),parameter :: Tprop = 10d0/gamma
  real(8),parameter :: dt = 10d0
  integer,parameter :: nt = aint(Tprop/dt)+1
  integer,parameter :: Nw = 200
  real(8),parameter :: wi = 2.0d0/ev, wf = 4d0/ev, dw = (wf-wi)/Nw
  real(8) :: delta_Dj_w(0:Nw)

! Laser parameters
  real(8),parameter :: E0 = 1d6*a_B*1d-10/ev

end module global_variables
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
program main
  implicit none

  call calc_delta_Dj_parabolic
  call calc_delta_Dj_kane


end program main
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine calc_delta_Dj_parabolic
  use global_variables
  implicit none
  integer :: iw, ikz, ikr, it
  real(8) :: ww,kz,kr,eps0, c1,c2,ss,tt
  complex(8) :: ztheta, ztheta0

  delta_Dj_w = 0d0

!$omp parallel default(shared), private(iw,ikz,ikr,it,ww,kz,kr,eps0,c1,c2,tt,ztheta,ztheta0,ss) 
!$omp do
  do iw = 0,nw
    write(*,*)'iw=',iw
    ww = wi + dw*iw

    do ikz = -nkz,nkz
      kz = dkz*ikz

      do ikr = 1,nkr
        kr = dkr*ikr

        eps0 = eps_g + 0.5d0*(kz**2 + kr**2)/mass_r
        c1 = kz*E0/(2d0*mass_r)
        c2 = E0**2/(6d0*mass_r)

        ss = 0d0
        do it = 0,nt
          tt = dt*it
          ztheta0 = -zI*(eps0 - ww)*tt -gamma*tt
          ztheta  = ztheta0 -zI*(-c1*tt**2 + c2*tt**3)

          ss = ss + (exp(ztheta) - exp(ztheta0))

        end do
        ss = ss*dt*dkz*dkr*2d0*pi*kr

        delta_Dj_w(iw) = delta_Dj_w(iw) + ss

      end do
    end do
  end do
!$omp end do
!$omp end parallel  

  open(20,file='delta_jdos_parabolic_band.out')
  do iw =0,nw
    ww = wi + dw*iw
    write(20,*)ww,delta_Dj_w(iw)
  end do
  close(20)


end subroutine calc_delta_Dj_parabolic
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine calc_delta_Dj_kane
  use global_variables
  implicit none
  integer :: iw, ikz, ikr, it
  real(8) :: ww,kz,kr,eps0, ss,tt
  real(8) :: a,b,c,d,F0,Ft
  complex(8) :: ztheta, ztheta0

  delta_Dj_w = 0d0

!$omp parallel default(shared), private(iw,ikz,ikr,it,ww,kz,kr,eps0,a,b,c,d,F0,Ft,tt,ztheta,ztheta0,ss) 
!$omp do
  do iw = 0,nw
    write(*,*)'iw=',iw
    ww = wi + dw*iw

    do ikz = -nkz,nkz
      kz = dkz*ikz

      do ikr = 1,nkr
        kr = dkr*ikr

        eps0 = eps_g*sqrt(1d0  + (kz**2 + kr**2)/(mass_r*eps_g))
        a = E0**2/(mass_r*eps_g)
        b = -2d0*kz*E0/(mass_r*eps_g)
        c = 1d0 + (kz**2+kr**2)/(mass_r*eps_g)
        d = (b**2 - 4d0*a*c)/(8d0*a)
        F0 = b*sqrt(c)/(4d0*a)-d*log(abs(b+sqrt(a*c)))/sqrt(a)

        ss = 0d0
        do it = 0,nt
          tt = dt*it
          ztheta0 = -zI*(eps0 - ww)*tt -gamma*tt

          Ft = (2d0*a*tt+b)/(4d0*a)*sqrt(a*tt**2+b*tt+c)-d*log(abs(2*a*tt+b+2*sqrt(a*(a*tt**2+b*tt+c))))/sqrt(a)
          ztheta  = ztheta0 -zI*(Ft - F0)

          ss = ss + (exp(ztheta) - exp(ztheta0))

        end do
        ss = ss*dt*dkz*dkr*2d0*pi*kr

        delta_Dj_w(iw) = delta_Dj_w(iw) + ss

      end do
    end do
  end do
!$omp end do
!$omp end parallel  

  open(20,file='delta_jdos_parabolic_kane.out')
  do iw =0,nw
    ww = wi + dw*iw
    write(20,*)ww,delta_Dj_w(iw)
  end do
  close(20)


end subroutine calc_delta_Dj_kane
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
