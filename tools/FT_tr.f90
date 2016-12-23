program main
  implicit none
  complex(8),parameter :: zI = (0d0,1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)
  integer,parameter :: Nt = 30000, Nw = 200
  real(8),parameter :: wi = 20d0/27.2, wf = 50d0/27.2d0, dw =(wf-wi)/dble(Nw)
  real(8),parameter :: gamma = 0.5d0/27.2,alpha = -1d-1 !- 8.166d-3
  real(8) :: tt(0:Nt+1),Act(0:Nt+1),jt(0:Nt+1),dt,jt_pp(0:Nt+1),jt_pn(0:Nt+1),jt_np(0:Nt+1)
  real(8) :: Act_pp(0:Nt+1),Act_pn(0:Nt+1),Act_np(0:Nt+1)
  integer :: it,iw
  real(8) :: f1,f2,f3,f4,f5,f6,f7,x
  complex(8) :: zs1,zs2, zs3,zs4,zAw(0:Nw),zjw(0:Nw),zjw_intra(0:Nw),zjw_inter(0:Nw)
  real(8) :: ww,tav,tav_pump,tav_probe
  complex(8) :: zsigma,zeps,zsigma0,zeps0,zEw
  integer,parameter :: Nt_delay = 30
  integer :: it_delay
  character(50) :: cit,filename

!  write(*,*)pi
  open(20,file='Pump_only_Act_jtz.out')
  do it = 0,Nt
    read(20,*)tt(it),Act_pn(it),jt_pn(it)
  end do
  close(20)
  open(20,file='Probe_only_Act_jtz.out')
  do it = 0,Nt
    read(20,*)tt(it),Act_np(it),jt_np(it)
  end do
  close(20)
  tav_pump = sum(Act_pn(0:Nt)**2*tt(0:Nt))/sum(Act_pn(0:Nt)**2)
  tav_probe = sum(Act_np(0:Nt)**2*tt(0:Nt))/sum(Act_np(0:Nt)**2)

  open(30,file='zeps_tr.out')

  do it_delay = 1,Nt_delay
     write(cit,"(I3.3)")it_delay
     filename=trim(cit)//"_Act_jtz.out"
  open(20,file=filename)
  do it = 0,Nt
    read(20,*)tt(it),Act_pp(it),jt_pp(it)
  end do
  close(20)
  dt = tt(1) - tt(0)
  jt(:) = jt_pp(:)-jt_pn(:) !+ alpha*Act(:)
  Act(:) = Act_pp(:) - Act_pn(:)

  tav = sum(Act(0:Nt)**2*tt(0:Nt))/sum(Act(0:Nt)**2)

  write(*,*)(tav-tav_pump)*0.02418

  do iw = 0,Nw
!    write(*,*)'iw=',iw,'/',Nw
    ww = wi + dble(iw)*dw

    zs1 = 0d0; zs2 = 0d0
    zs3 = 0d0; zs4 = 0d0
    do it = 0,Nt
!      x = (tt(it)-tt(0))/(tt(Nt)-tt(0))
!      f1 = 1d0 -3d0*x**2 + 2d0*x**3
      if(tt(it) > tav -1d0/0.02418d0)then
        f1 = exp(-gamma*(tt(it)-tav))
        zs1 = zs1 + jt(it)*exp(zI*ww*tt(it))*f1
        zs2 = zs2 + Act(it)*exp(zI*ww*tt(it))*f1
      end if
      if(tt(it) > tav_probe-1d0/0.02418d0)then
        f1 = exp(-gamma*(tt(it)-tav_pump))
        zs3 = zs3 + jt_np(it)*exp(zI*ww*tt(it))*f1
        zs4 = zs4 + Act_np(it)*exp(zI*ww*tt(it))*f1
      end if

    end do
    zsigma = -zs1/(zI*ww*zs2)
    zsigma0= -zs3/(zI*ww*zs4)
    zeps = 1d0 + 4d0*pi*zI*zsigma/ww
    zeps0= 1d0 + 4d0*pi*zI*zsigma0/ww

    write(30,'(999e26.16e3)')tav-tav_pump,ww,zeps,zeps0
  end do
  write(30,*)
  end do
  close(30)


end program main
