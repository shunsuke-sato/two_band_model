program main
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: ev = 1d0/27.2114d0, fs = 1d0/0.024189d0
  integer,parameter :: nt = 10400, nw = 512
  real(8),parameter :: ww_i = 0.1d0*ev, ww_f = 30d0*ev, dw = (ww_f-ww_i)/nw
  real(8),parameter :: tpulse = 20d0*fs
  real(8) :: tt(0:nt),Act(0:nt),jtz(0:nt),jtz_intra(0:nt),jtz_inter(0:nt),Et(0:nt)
  complex(8) :: zEw, zjw, zjw_intra, zjw_inter, zfact
  real(8) :: window(0:nt)
  real(8) :: dt, ww, xx
  integer :: iw, it


  open(20,file='Act_jtz.out')
  do it = 0,nt
    read(20,*)tt(it),Act(it),jtz(it),jtz_intra(it),jtz_inter(it),Et(it)
  end do
  close(20)
  dt = tt(1) - tt(0)


  window = 0d0
  do it = 0, nt
    xx = tt(it) - 0.5d0*tpulse
    if(abs(xx) < 0.5d0*tpulse)then
      window(it) = cos(pi*xx/tpulse)**4
    end if
  end do

  jtz = jtz * window
  jtz_intra = jtz_intra * window
  jtz_inter = jtz_inter * window

  open(30,file='hhg_spec.out')
  do iw = 0, nw
    ww = ww_i + dw*iw

    zEw = 0d0
    zjw = 0d0
    zjw_intra = 0d0
    zjw_inter = 0d0

    do it = 0, nt
      zfact = exp(zi*ww*tt(it))
      zEw = zEw + Et(it)*zfact
      zjw = zjw + jtz(it)*zfact
      zjw_intra = zjw_intra + jtz_intra(it)*zfact
      zjw_inter = zjw_inter + jtz_inter(it)*zfact
    end do

    zEw = zEw*dt
    zjw = zjw*dt
    zjw_intra = zjw_intra*dt
    zjw_inter = zjw_inter*dt

    
    write(30,"(999e26.16e3)")ww,abs(zEw)**2, abs(zjw)**2, abs(zjw_intra)**2 &
      , abs(zjw_inter)**2

  end do
  close(30)
  


end program main
