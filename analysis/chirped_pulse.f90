!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module global_variables
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0,1d0)

  real(8),parameter :: fs = 0.02419d0
  real(8),parameter :: ev = 27.211d0

  real(8) :: sigma0 = 200d-3/fs*sqrt(2d0)
  real(8) :: omega0 = 40d0/ev

  real(8) :: sigma = 1.333d0/fs*sqrt(2d0)
  real(8) :: beta = (1.55d0/ev)/(1.33/fs)
  real(8) :: sigma0_w, sigma_w



  real(8) :: norm, phi

  integer,parameter :: Nt = 20000
  integer,parameter :: Nw = 200
  real(8) :: Et0(-Nt:Nt), Et(-Nt:Nt),dt

end module global_variables


program main
  use global_variables
  implicit none
  real(8) :: tmp, tt
  integer :: it, iw
  real(8) :: ww, dw
  complex(8) :: exp_c, zE0, zE

  sigma_w = 1d0/sqrt((1d0/(2d0*sigma**2))/((1d0/(2d0*sigma**2))**2+beta**2))
  write(*,"(A)")"Parameters for chirped pulse:"
  write(*,"(A,2x,es26.16e3,A)")"Pulse width = ",sigma*fs/sqrt(2d0)," fs"
  write(*,"(A,2x,es26.16e3,A)")"FWHM        = ",sigma*fs/sqrt(2d0)*2d0*sqrt(2d0*log(2d0))," fs"
  write(*,"(A,2x,es26.16e3,A)")"Photon energy = ",omega0*ev," eV"
  write(*,"(A,2x,es26.16e3,A)")"beta        = ",beta*ev/fs," eV/fs"
  write(*,"(A,2x,es26.16e3,A)")"Band width    = ",sigma_w*ev," eV"

  sigma0 = sqrt(0.5d0/((1d0/2d0/sigma**2)+beta**2/(1d0/2d0/sigma**2)))
  sigma0_w = sqrt(1d0/2d0/sigma0**2)

  write(*,"(A)")"Parameters for un-chirped pulse:"
  write(*,"(A,2x,es26.16e3,A)")"Pulse width   = ",sigma0*fs/sqrt(2d0)," fs"
  write(*,"(A,2x,es26.16e3,A)")"FWHM          = ",sigma0*fs/sqrt(2d0)*2d0*sqrt(2d0*log(2d0))," fs"
  write(*,"(A,2x,es26.16e3,A)")"Photon energy = ",omega0*ev," eV"
  write(*,"(A,2x,es26.16e3,A)")"Band width    = ",sigma0_w*ev," eV"

  norm  = sqrt(sqrt( ((1d0/2d0/sigma**2)**2 + beta**2)/(1d0/2d0/sigma0**2)**2 ))

  dt = sigma*5d0/Nt

  do it = -Nt,Nt
    tt = dt*it
    Et0(it) = zI*exp(zI*omega0*tt)*exp(-0.5d0*tt**2/sigma0**2)
    Et(it) = norm*zI*exp(zI*(omega0+beta*tt)*tt)*exp(-0.5d0*tt**2/sigma**2)
  end do
  open(21,file="chirped_pulse.dat")
  do it = -Nt,Nt
    tt = dt*it
    write(21,"(999e26.16e3)")tt*fs*1d3,Et0(it),Et(it)
  end do
  close(21)

  open(21,file="chirped_pulse_w.dat")
  dw = sigma0_w * 5d0/dble(Nw)
  do iw = -Nw,Nw
    ww = dw*iw + omega0

    zE0 = 0d0
    zE  = 0d0
    do it = -Nt,Nt
      exp_c = exp(zI*ww*it*dt)
      zE0 = zE0 + exp_c*Et0(it)
      zE  = zE  + exp_c*Et(it)
    end do
    zE0 = zE0*dt
    zE = zE*dt

    write(21,"(999e26.16e3)")ww*ev,abs(zE0)**2,abs(zE)**2

  end do
  close(21)

end program main
