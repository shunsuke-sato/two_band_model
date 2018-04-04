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
  integer,parameter :: nkx = 32, nky = 32
  real(8),parameter :: kx_max = 1d0, ky_max = 1d0
  real(8),parameter :: kx_shift = 0d0, ky_shift = 0d0
! Fermi-Dirac distribution
  real(8),parameter :: mu_F = 0d0/ev
  real(8),parameter :: kbT  = 80d0/11604.505d0/ev

  real(8) :: kx(nkx),ky(nky)
  real(8) :: dkx,dky
  complex(8) :: zpsi(2,2,nkx,nky)
  real(8) :: occ(2,nkx,nky), eps(2,nkx,nky)

! time propagation
  real(8),parameter :: Tprop = 2.2d3/fs, dt = 0.08d0
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
  real(8) :: kxt, kyt, theta
  complex(8) :: zs


  do ikx = 1,nkx
    kx(ikx) = -kx_max + 2d0*kx_max*dble(ikx-1)/dble(nkx-1)
  end do
  dkx = 2d0*kx_max/dble(nkx-1)
  kx = kx + kx_shift

  do iky = 1,nky
    ky(iky) = -ky_max + 2d0*ky_max*dble(iky-1)/dble(nky-1)
  end do
  dky = 2d0*ky_max/dble(nky-1)
  ky = ky + ky_shift

  eps = 0d0

  if(delta_gap == 0d0)then ! Dirac cone
    do ikx = 1, nkx
      kxt = kx(ikx)
      do iky = 1,nky
        kyt = ky(iky)

        zs = tau_z*kxt+zI*kyt
        if(zs /= 0d0)then
          theta = -aimag(log(-zs))
        else
          theta = 0d0
        end if

        zpsi(1,1,ikx,iky) = exp(zI*theta)/sqrt(2d0)
        zpsi(2,1,ikx,iky) = 1d0/sqrt(2d0)
        eps(1,ikx,iky) = -velocity*sqrt(tau_z**2*kxt**2 + kyt**2)
        
        zpsi(1,2,ikx,iky) = exp(zI*(theta+pi))/sqrt(2d0)
        zpsi(2,2,ikx,iky) = 1d0/sqrt(2d0)
        eps(2,ikx,iky) = velocity*sqrt(tau_z**2*kxt**2 + kyt**2)
        
      end do
    end do
  else
    stop 'Kane band is not implemented yet.'
  end if
  
! Fermi-Dirac distribution
  do ikx = 1,nkx
    do iky = 1,nky

      occ(1,ikx,iky) = 2d0/(exp((eps(1,ikx,iky)-mu_F)/kbT)+1d0)
      occ(2,ikx,iky) = 2d0/(exp((eps(2,ikx,iky)-mu_F)/kbT)+1d0)

    end do
  end do

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

    if(it*dt > (1d3+10d0+1d0)/fs)then
      call e_h_distribution(it)
      exit
    end if

  end do
  close(21)


end subroutine time_propagation
!----------------------------------------------------------------------------------------!
subroutine dt_evolve(it)
  implicit none
  integer,intent(inout) :: it

  call dt_evolve_Taylor(it)
!  call dt_evolve_Magnus(it)
!  call dt_evolve_Magnus_1st(it)


end subroutine dt_evolve
!!----------------------------------------------------------------------------------------!
!subroutine dt_evolve_Magnus(it)
!  use global_variables
!  implicit none
!  integer,intent(in) :: it
!  real(8) :: acx0, acx1, acx2
!  real(8) :: acy0, acy1, acy2
!  integer :: ikx, iky
!  real(8) :: kxt, kyt,kxt_eff, kyt_eff
!  real(8) :: delta,lambda(2)
!  complex(8) :: zHeff(2,2),zalpha,zx,zvec(2,2),zc(2)
!  real(8) :: const,ss
!
!
!  if(delta_gap /= 0d0)stop 'delta_gap has to be zero in Magnus propagator.'
!
!  acx0 = ac(1,it)
!  acx2 = (ac(1,it)-2d0*ac_dt2(1,it)+ac(1,it+1))/(0.5d0*dt)**2
!  acx1 = (ac(1,it+1)-acx0-0.5d0*dt**2*acx2)/dt
!
!  acy0 = ac(2,it)
!  acy2 = (ac(2,it)-2d0*ac_dt2(2,it)+ac(2,it+1))/(0.5d0*dt)**2
!  acy1 = (ac(2,it+1)-acy0-0.5d0*dt**2*acy2)/dt
!
!  const = -dt**2/6d0*tau_z*velocity**2
!
!
!!$omp parallel default(shared), private(ikx,iky,kxt,kyt,kxt_eff,kyt_eff,ss,delta,zalpha, &
!!$omp & lambda,zx,zvec,zc) 
!!$omp do collapse(2)
!  do ikx = 1,nkx
!    do iky = 1,nky
!      kxt = kx(ikx) + acx0
!      kxt_eff = kxt + 0.5d0*dt*acx1 + dt**2/6d0*acx2
!      kyt = ky(iky) + acy0
!      kyt_eff = kyt + 0.5d0*dt*acy1 + dt**2/6d0*acy2
!      
!      ss = kxt*acy1-acx1*kyt
!!      zHeff(1,1) =  const*ss
!!      zHeff(2,1) =  velocity*(tau_z*kxt_eff + zI*kyt_eff)
!!      zHeff(1,2) =  velocity*(tau_z*kxt_eff - zI*kyt_eff)
!!      zHeff(2,2) = -const*ss
!
!      delta  = const*ss
!      zalpha = velocity*(tau_z*kxt_eff + zI*kyt_eff)
!!
!      if(delta >= 0d0)then
!        
!! vector 1
!        lambda(1) = sqrt(delta**2 + abs(zalpha)**2)
!        zx = zalpha/(lambda(1) + delta)
!        ss = 1d0/sqrt(1d0 + abs(zx)**2)
!        zvec(1,1) = ss
!        zvec(2,1) = ss*zx
!
!! vector 2
!        lambda(2) = -lambda(1)
!        zvec(1,2) = -ss*conjg(zx)
!        zvec(2,2) =  ss
!
!      else
!
!! vector 1
!        lambda(1) = sqrt(delta**2 + abs(zalpha)**2)
!        zx = conjg(zalpha)/(lambda(1) - delta)
!        ss = 1d0/sqrt(1d0 + abs(zx)**2)
!        zvec(1,1) = ss*zx
!        zvec(2,1) = ss
!
!! vector 2
!        lambda(2) = -lambda(1)
!        zvec(1,2) =  ss
!        zvec(2,2) = -ss*conjg(zx)
!
!      end if
!
!! projection
!      zc(1) = conjg(zvec(1,1))*zpsi(1,ikx,iky) + conjg(zvec(2,1))*zpsi(2,ikx,iky)
!      zc(2) = conjg(zvec(1,2))*zpsi(1,ikx,iky) + conjg(zvec(2,2))*zpsi(2,ikx,iky)
!
!! propagation
!      zc(1) = zc(1)*exp(-zI*lambda(1)*dt)
!      zc(2) = zc(2)*exp(-zI*lambda(2)*dt)
!
!! update wavefunction
!      zpsi(:,ikx,iky) = zc(1)*zvec(:,1) + zc(2)*zvec(:,2)
!      
!
!    end do
!  end do
!!$omp end do
!!$omp end parallel  
!
!end subroutine dt_evolve_Magnus
!!----------------------------------------------------------------------------------------!
!subroutine dt_evolve_Magnus_1st(it)
!  use global_variables
!  implicit none
!  integer,intent(in) :: it
!  real(8) :: acx0, acx1, acx2
!  real(8) :: acy0, acy1, acy2
!  integer :: ikx, iky
!  real(8) :: kxt, kyt,kxt_eff, kyt_eff
!  real(8) :: delta,lambda(2)
!  complex(8) :: zHeff(2,2),zalpha,zx,zvec(2,2),zc(2)
!  real(8) :: const,ss
!
!
!  if(delta_gap /= 0d0)stop 'delta_gap has to be zero in Magnus propagator.'
!
!  acx0 = ac(1,it)
!  acx2 = (ac(1,it)-2d0*ac_dt2(1,it)+ac(1,it+1))/(0.5d0*dt)**2
!  acx1 = (ac(1,it+1)-acx0-0.5d0*dt**2*acx2)/dt
!
!  acy0 = ac(2,it)
!  acy2 = (ac(2,it)-2d0*ac_dt2(2,it)+ac(2,it+1))/(0.5d0*dt)**2
!  acy1 = (ac(2,it+1)-acy0-0.5d0*dt**2*acy2)/dt
!
!  delta = 0d0
!
!
!!$omp parallel default(shared), private(ikx,iky,kxt,kyt,kxt_eff,kyt_eff,ss,zalpha, &
!!$omp & lambda,zx,zvec,zc) 
!!$omp do collapse(2)
!  do ikx = 1,nkx
!    do iky = 1,nky
!      kxt = kx(ikx) + acx0
!      kxt_eff = kxt + 0.5d0*dt*acx1 + dt**2/6d0*acx2
!      kyt = ky(iky) + acy0
!      kyt_eff = kyt + 0.5d0*dt*acy1 + dt**2/6d0*acy2
!      
!      zalpha = velocity*(tau_z*kxt_eff + zI*kyt_eff)
!!
!! vector 1
!      lambda(1) = sqrt(delta**2 + abs(zalpha)**2)
!      zx = zalpha/(lambda(1) + delta)
!      ss = 1d0/sqrt(1d0 + abs(zx)**2)
!      zvec(1,1) = ss
!      zvec(2,1) = ss*zx
!
!! vector 2
!      lambda(2) = -lambda(1)
!      zvec(1,2) = -ss*conjg(zx)
!      zvec(2,2) =  ss
!
!
!! projection
!      zc(1) = conjg(zvec(1,1))*zpsi(1,ikx,iky) + conjg(zvec(2,1))*zpsi(2,ikx,iky)
!      zc(2) = conjg(zvec(1,2))*zpsi(1,ikx,iky) + conjg(zvec(2,2))*zpsi(2,ikx,iky)
!
!! propagation
!      zc(1) = zc(1)*exp(-zI*lambda(1)*dt)
!      zc(2) = zc(2)*exp(-zI*lambda(2)*dt)
!
!! update wavefunction
!      zpsi(:,ikx,iky) = zc(1)*zvec(:,1) + zc(2)*zvec(:,2)
!      
!
!    end do
!  end do
!!$omp end do
!!$omp end parallel  
!
!end subroutine dt_evolve_Magnus_1st
!!----------------------------------------------------------------------------------------!
subroutine dt_evolve_Taylor(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ikx, iky
  real(8) :: kxt, kyt
  complex(8) :: zpsi_t(2,2,0:4)

!$omp parallel default(shared), private(ikx,iky,kxt,kyt,zpsi_t) 
!$omp do collapse(2)
  do ikx = 1,nkx
    do iky = 1,nky

      kxt = kx(ikx) + ac(1,it)
      kyt = ky(iky) + ac(2,it)
! RK1
      zpsi_t(:,:,0) = zpsi(:,:,ikx,iky)
      zpsi_t(1,:,1) = 0.5d0*delta_gap*zpsi_t(1,:,0) &
        + velocity*(tau_z*kxt - zI*kyt)*zpsi_t(2,:,0) 
      zpsi_t(2,:,1) = velocity*(tau_z*kxt + zI*kyt)*zpsi_t(1,:,0)  &
        - 0.5d0*delta_gap*zpsi_t(2,:,0)

      kxt = kx(ikx) + ac_dt2(1,it)
      kyt = ky(iky) + ac_dt2(2,it)
! RK2
      zpsi_t(:,:,0) = zpsi(:,:,ikx,iky) - zI*0.5d0*dt*zpsi_t(:,:,1)

      zpsi_t(1,:,2) = 0.5d0*delta_gap*zpsi_t(1,:,0) &
        + velocity*(tau_z*kxt - zI*kyt)*zpsi_t(2,:,0) 
      zpsi_t(2,:,2) = velocity*(tau_z*kxt + zI*kyt)*zpsi_t(1,:,0)  &
        - 0.5d0*delta_gap*zpsi_t(2,:,0)

! RK3
      zpsi_t(:,:,0) = zpsi(:,:,ikx,iky) - zI*0.5d0*dt*zpsi_t(:,:,2)

      zpsi_t(1,:,3) = 0.5d0*delta_gap*zpsi_t(1,:,0) &
        + velocity*(tau_z*kxt - zI*kyt)*zpsi_t(2,:,0) 
      zpsi_t(2,:,3) = velocity*(tau_z*kxt + zI*kyt)*zpsi_t(1,:,0)  &
        - 0.5d0*delta_gap*zpsi_t(2,:,0)

      kxt = kx(ikx) + ac(1,it+1)
      kyt = ky(iky) + ac(2,it+1)

! RK4
      zpsi_t(:,:,0) = zpsi(:,:,ikx,iky) - zI*dt*zpsi_t(:,:,3)

      zpsi_t(1,:,4) = 0.5d0*delta_gap*zpsi_t(1,:,0) &
        + velocity*(tau_z*kxt - zI*kyt)*zpsi_t(2,:,0) 
      zpsi_t(2,:,4) = velocity*(tau_z*kxt + zI*kyt)*zpsi_t(1,:,0)  &
        - 0.5d0*delta_gap*zpsi_t(2,:,0)

      zpsi(:,:,ikx,iky) = zpsi(:,:,ikx,iky) -zI*dt/6d0*( &
        zpsi_t(:,:,1) + 2d0*zpsi_t(:,:,2) + 2d0*zpsi_t(:,:,3) + zpsi_t(:,:,4) )
      
    end do
  end do
!$omp end do
!$omp end parallel  

end subroutine dt_evolve_Taylor
!----------------------------------------------------------------------------------------!
subroutine init_ac
  use global_variables
  implicit none
  integer :: it
  real(8) :: tt,tt2, xx, xx2
  real(8) :: E0,omega,tpulse
  real(8) :: E_SD, T_SD
  real(8) :: gap_Floquet


!  gap_Floquet = 0.1d0/ev
  E0 = 1d7*b_a*1d-10/ev
  omega = 0.19074d0/ev ! 6.5 micron
  tpulse = 1d3/fs

!  E0 = omega*0.5d0/velocity*sqrt( &
!    (gap_Floquet + omega)**2 - omega**2 &
!    )

  write(*,"(A,2x,e26.16e3,A)")"Field strength =",E0*ev/b_a*1d10, "eV/m"

! Source-drain
  E_SD = 1d-8
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

! circular laser
    xx = tt - T_SD - 0.5d0*tpulse
    xx2 = tt2 - T_SD - 0.5d0*tpulse

     if( abs(xx)<0.5d0*tpulse )then
        ac(1,it) = ac(1,it) - E0/omega*cos(pi*xx/tpulse)**4*sin(omega*xx)
        ac(2,it) = ac(2,it) - E0/omega*cos(pi*xx/tpulse)**4*cos(omega*xx)
     end if

     if( abs(xx2)<0.5d0*tpulse )then
        ac_dt2(1,it) = ac_dt2(1,it) - E0/omega*cos(pi*xx2/tpulse)**4*sin(omega*xx2)
        ac_dt2(2,it) = ac_dt2(2,it) - E0/omega*cos(pi*xx2/tpulse)**4*cos(omega*xx2)
     end if
     
  end do


end subroutine init_ac
!----------------------------------------------------------------------------------------!
subroutine current(jxy,it)
  use global_variables
  implicit none
  real(8),intent(out) :: jxy(2)
  integer,intent(in) :: it
  real(8) :: jx,jy
  integer :: ikx, iky
  real(8) :: kxt, kyt
  real(8) :: jxt,jyt,jx0,jy0,xx
  complex(8) :: zs

  jx = 0d0
  jy = 0d0

!$omp parallel default(shared), private(ikx,iky,kxt,kyt,jxt,jyt,zs,jx0,jy0) 
!$omp do reduction(+:jx,jy) collapse(2)
  do ikx = 1,nkx
    do iky = 1,nky

       kxt = kx(ikx) + ac(1,it)
       kyt = ky(iky) + ac(2,it)
       
       jxt = 2d0*occ(1,ikx,iky)*real(zpsi(1,1,ikx,iky)*conjg(zpsi(2,1,ikx,iky))) &
            +2d0*occ(2,ikx,iky)*real(zpsi(1,2,ikx,iky)*conjg(zpsi(2,2,ikx,iky)))
       jyt = 2d0*occ(1,ikx,iky)*real(zI*zpsi(1,1,ikx,iky)*conjg(zpsi(2,1,ikx,iky))) &
            +2d0*occ(2,ikx,iky)*real(zI*zpsi(1,2,ikx,iky)*conjg(zpsi(2,2,ikx,iky)))
       
       zs =  -(tau_z*kxt-zI*kyt)&
            /sqrt(kxt**2+kyt**2)
       

       jx0 =  2d0*0.5d0*real(zs)   *(occ(1,ikx,iky)-occ(2,ikx,iky))
       jy0 =  2d0*0.5d0*real(zI*zs)*(occ(1,ikx,iky)-occ(2,ikx,iky))
       
       jx = jx + jxt - jx0
       jy = jy + jyt - jy0
       
     end do
  end do
!$omp end do
!$omp end parallel  

  jxy(1) = jx*dkx*dky
  jxy(2) = jy*dkx*dky

end subroutine current
!----------------------------------------------------------------------------------------!
subroutine intra_current(jxy,it)
  use global_variables
  implicit none
  real(8),intent(out) :: jxy(2)
  integer,intent(in) :: it
  real(8) :: jx,jy
  real(8) :: jxt,jyt,jx0,jy0
  integer :: ikx, iky
  real(8) :: kxt, kyt
  complex(8) :: zeig_t(2,2), zs
  real(8) :: theta, ovl_t(2,2)

  jx = 0d0
  jy = 0d0

!$omp parallel default(shared), private(ikx,iky,kxt,kyt,jxt,jyt,zs,jx0,jy0,zeig_t,ovl_t,theta) 
!$omp do reduction(+:jx,jy) collapse(2)
  do ikx = 1,nkx
    do iky = 1,nky

       kxt = kx(ikx) + ac(1,it)
       kyt = ky(iky) + ac(2,it)

       zs = tau_z*kxt+zI*kyt
       if(zs /= 0d0)then
         theta = -aimag(log(-zs))
       else
         theta = 0d0
       end if

       zeig_t(1,1) = exp(zI*theta)/sqrt(2d0)
       zeig_t(2,1) = 1d0/sqrt(2d0)

       zeig_t(1,2) = exp(zI*(theta+pi))/sqrt(2d0)
       zeig_t(2,2) = 1d0/sqrt(2d0)

       ovl_t(1,1) = abs(sum(conjg(zeig_t(:,1))*zpsi(:,1,ikx,iky)))**2
       ovl_t(2,1) = abs(sum(conjg(zeig_t(:,2))*zpsi(:,1,ikx,iky)))**2
       ovl_t(1,2) = abs(sum(conjg(zeig_t(:,1))*zpsi(:,2,ikx,iky)))**2
       ovl_t(2,2) = abs(sum(conjg(zeig_t(:,2))*zpsi(:,2,ikx,iky)))**2

       zs =  -(tau_z*kxt-zI*kyt)&
         /sqrt(kxt**2+kyt**2)
       
       jx0 =  2d0*0.5d0*real(zs)   *(occ(1,ikx,iky)-occ(2,ikx,iky))
       jy0 =  2d0*0.5d0*real(zI*zs)*(occ(1,ikx,iky)-occ(2,ikx,iky))
       jxt =  2d0*0.5d0*real(zs)*( &
          (occ(1,ikx,iky)*ovl_t(1,1)+occ(2,ikx,iky)*ovl_t(1,2)) &
         -(occ(1,ikx,iky)*ovl_t(2,1)+occ(2,ikx,iky)*ovl_t(2,2)))
       jyt =  2d0*0.5d0*real(zI*zs)*( &
          (occ(1,ikx,iky)*ovl_t(1,1)+occ(2,ikx,iky)*ovl_t(1,2)) &
         -(occ(1,ikx,iky)*ovl_t(2,1)+occ(2,ikx,iky)*ovl_t(2,2)))

       jx = jx + jxt - jx0
       jy = jy + jyt - jy0


     end do
   end do
!$omp end do
!$omp end parallel  


  jxy(1) = jx*dkx*dky
  jxy(2) = jy*dkx*dky

end subroutine intra_current
!----------------------------------------------------------------------------------------!
subroutine e_h_distribution(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8) :: occ_t(2,nkx,nky)
  character(256) :: cit, cfilename
  integer :: ikx, iky
  real(8) :: kxt, kyt
  complex(8) :: zeig_t(2,2), zs
  real(8) :: theta
  
  do ikx = 1,nkx
    do iky = 1,nky

      kxt = kx(ikx) + ac(1,it)
      kyt = ky(iky) + ac(2,it)

      zs = tau_z*kxt+zI*kyt
      if(zs /= 0d0)then
        theta = -aimag(log(-zs))
      else
        theta = 0d0
      end if

      zeig_t(1,1) = exp(zI*theta)/sqrt(2d0)
      zeig_t(2,1) = 1d0/sqrt(2d0)

      zeig_t(1,2) = exp(zI*(theta+pi))/sqrt(2d0)
      zeig_t(2,2) = 1d0/sqrt(2d0)

      occ_t(1,ikx,iky) = occ(1,ikx,iky)*abs(sum(conjg(zeig_t(:,1))*zpsi(:,1,ikx,iky)))**2 &
                        +occ(2,ikx,iky)*abs(sum(conjg(zeig_t(:,1))*zpsi(:,2,ikx,iky)))**2 

      occ_t(2,ikx,iky) = occ(1,ikx,iky)*abs(sum(conjg(zeig_t(:,2))*zpsi(:,1,ikx,iky)))**2 &
                        +occ(2,ikx,iky)*abs(sum(conjg(zeig_t(:,2))*zpsi(:,2,ikx,iky)))**2 


    end do
  end do

  write(cit,"(I9.9)")it
  cfilename = trim(cit)//"_eh_dist.out"
  open(900,file=cfilename)
  write(900,"(A,2x,I7,I7)")"# nkx,nky=",nkx,nky
  write(900,"(A,2x,2e26.16e3)")"# acx,acy=",ac(:,it)
  write(900,"(A)")"# kx, ky, (-)hole dist., (+)electron dist."
  do ikx = 1,nkx
    do iky = 1,nky

      kxt = kx(ikx) + ac(1,it)
      kyt = ky(iky) + ac(2,it)

      write(900,"(999e26.16e3)")kxt,kyt,occ_t(1,ikx,iky)-2d0,occ_t(2,ikx,iky)
    end do
    write(900,*)
  end do
  

  close(900)
  

end subroutine e_h_distribution
!----------------------------------------------------------------------------------------!
