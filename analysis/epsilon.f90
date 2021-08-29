program main
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zI = (0d0,1d0)
  integer,parameter :: Nt = 6000, Nw = 500
  real(8),parameter :: wi = 0.1d0/27.211d0, wf = 20d0/27.211d0, dw = (wf-wi)/Nw
  real(8) :: tt(0:Nt),Act(0:Nt),jt(0:Nt)
  real(8) :: dt,ww,xx,ff
  integer :: it,iw
  complex(8) :: zAw,zEw,zJw,zsigma,zeps
  real(8) :: kick_impulsive

  kick_impulsive = 1d-4

  do it = 0,Nt
    read(*,*)tt(it),Act(it),jt(it)
  end do
  dt = tt(1) - tt(0)
  jt(0) = jt(0)*0.5d0

  open(21,file="zeps.out")
  do iw = 0,Nw
    ww = wi + dw*iw
    zAw = 0d0; zJw = 0d0
    do it = 0,Nt
      xx = (tt(it)-tt(0))/(tt(Nt) - tt(0))
      ff = 1d0 -3d0 * xx**2 + 2d0*xx**3
      zAw = zAw + exp(zI*ww*tt(it)) *Act(it)
      zjw = zjw + exp(zI*ww*tt(it)) *jt(it)*ff
    end do
    zEw = -(-zI*ww)*zAw*dt
    zJw = - zJw*dt
!    zsigma = zJw/zEw
    zsigma = zJw/kick_impulsive
    zeps = 1d0 + 4d0*pi*zi*zsigma/ww

    write(21,"(999e26.16e3)")ww,zeps,zsigma

  end do

  close(21)




end program main
