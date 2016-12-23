!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine single_cell
  use global_variables
  implicit none
  integer :: it,ikz,ikr
  real(8) :: jz_intra,jz_inter,Etz
  integer :: it_delay
  integer,parameter :: Ndelay = 30
  real(8),parameter :: Tdelay_fs_ini=-3d0,Tdelay_fs_fin=3d0
  real(8),parameter :: dTdelay_fs = (Tdelay_fs_fin-Tdelay_fs_ini)/dble(Ndelay)
  real(8) :: Iwcm2_1_t,Iwcm2_2_t
  character(50) :: citer,filename

  if(Nprocs /= 1)call err_finalize("Parallelization is not supported &
    for single-cell calculation.")

  read(*,*)kr_max,kz_max
  read(*,*)NKr, NKz
  read(*,*)Nt,dt
  read(*,*)envelope_1
  read(*,*)Iwcm2_1,omega_ev_1,tpulse_fs_1
  read(*,*)envelope_2
  read(*,*)Iwcm2_2,omega_ev_2,tpulse_fs_2
  read(*,*)Tdelay_fs

  Iwcm2_1_t=Iwcm2_1
  Iwcm2_2_t=Iwcm2_2


  do it_delay = -2,Ndelay
    write(*,*)'it_delay=',it_delay,'/',Ndelay
     if(it_delay == -2)then ! Pump only
        Iwcm2_1=Iwcm2_1_t
        Iwcm2_2=0d0
        Tdelay_fs = 0d0
     else if(it_delay == -1)then ! Probe only
        Iwcm2_1=0d0
        Iwcm2_2=Iwcm2_2_t
        Tdelay_fs = 0d0
     else ! Pump-probe
        Iwcm2_1=Iwcm2_1_t
        Iwcm2_2=Iwcm2_2_t
        Tdelay_fs = Tdelay_fs_ini + dTdelay_fs*dble(it_delay)
     end if



  call preparation
  call input_Ac

  do it = 0,Nt
!    write(*,*)'it=',it,'/',Nt
    call current(jz_intra,jz_inter)
    jtz_intra(it) = jz_intra; jtz_inter(it) = jz_inter
    jtz(it) = jtz_intra(it) + jtz_inter(it)

!=== deps_int, deps ====
    deps_int = deps_int + 0.5d0*(dt*0.5d0)*deps
!$omp parallel
!$omp do private(ikz, ikr)
    do ikz = -NKz,NKz
      kz(ikz) = kz0(ikz) + Act_dt2(it)
      do ikr = 1,NKr
        deps(ikr,ikz) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
      end do
    end do
!$omp end parallel
    deps_int = deps_int + 0.5d0*(dt*0.5d0)*deps
!=== deps_int, deps ====

    Etz = -(Act(it+1)-Act(it))/dt
    call dt_evolve(Etz)

!=== deps_int, deps ====
    deps_int = deps_int + 0.5d0*(dt*0.5d0)*deps
!$omp parallel
!$omp do private(ikz, ikr)
    do ikz = -NKz,NKz
      kz(ikz) = kz0(ikz) + Act(it+1)
      do ikr = 1,NKr
        deps(ikr,ikz) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
      end do
    end do
!$omp end parallel
    deps_int = deps_int + 0.5d0*(dt*0.5d0)*deps
!=== deps_int, deps ====

  end do


  if(it_delay == -2)then
     open(21,file='Pump_only_Act_jtz.out')
     do it = 0,Nt
        write(21,'(999e26.16e3)')dt*dble(it),Act(it),jtz(it),jtz_intra(it),jtz_inter(it),Act_dt2(it)
     end do
     close(21)
  else if(it_delay == -1)then
     open(21,file='Probe_only_Act_jtz.out')
     do it = 0,Nt
        write(21,'(999e26.16e3)')dt*dble(it),Act(it),jtz(it),jtz_intra(it),jtz_inter(it),Act_dt2(it)
     end do
     close(21)
  else
     write(citer,"(I3.3)")it_delay
     filename=trim(citer)//"_Act_jtz.out"
     open(21,file=filename)
     do it = 0,Nt
        write(21,'(999e26.16e3)')dt*dble(it),Act(it),jtz(it),jtz_intra(it),jtz_inter(it),Act_dt2(it)
     end do
     close(21)
  end if

  deallocate(zCt,deps,deps_int)
  deallocate(kz0,kz,kr)
  deallocate(Act,Act_dt2,jtz,jtz_intra,jtz_inter)

  end do
end subroutine single_cell
