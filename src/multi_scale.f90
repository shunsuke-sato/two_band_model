!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine multi_scale
  use global_variables_ms
  implicit none

  if(myrank == 0)then
    read(*,*)kr_max,kz_max
    read(*,*)NKr, NKz
    read(*,*)Nt,dt
    read(*,*)Iwcm2_1,omega_ev_1,tpulse_fs_1
    read(*,*)Iwcm2_2,omega_ev_2,tpulse_fs_2
    read(*,*)Tdelay_fs
    read(*,*)Nx_L,Nx_R,Mx
    read(*,*)Hx
  end if

  call MPI_BCAST(kr_max,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(kz_max,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(NKr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(NKz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(Nt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(dt,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(Iwcm2_1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(omega_ev_1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(tpulse_fs_1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(Iwcm2_2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(omega_ev_2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(tpulse_fs_2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(Tdelay_fs,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(Nx_L,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(Nx_R,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(Mx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(Hx,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  

  call preparation_ms
!  call input_Ac
!
!
!
!  do it = 0,Nt
!    write(*,*)'it=',it,'/',Nt
!    call current(jz_intra,jz_inter)
!    jtz_intra(it) = jz_intra; jtz_inter(it) = jz_inter
!    jtz(it) = jtz_intra(it) + jtz_inter(it)
!
!!=== deps_int, deps ====
!    deps_int = deps_int + 0.5d0*(dt*0.5d0)*deps
!    do ikz = -NKz,NKz
!      kz(ikz) = kz0(ikz) + Act_dt2(it)
!      do ikr = 1,NKr
!        deps(ikr,ikz) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
!      end do
!    end do
!    deps_int = deps_int + 0.5d0*(dt*0.5d0)*deps
!!=== deps_int, deps ====
!
!    Etz = -(Act(it+1)-Act(it))/dt
!    call dt_evolve(Etz)
!
!!=== deps_int, deps ====
!    deps_int = deps_int + 0.5d0*(dt*0.5d0)*deps
!    do ikz = -NKz,NKz
!      kz(ikz) = kz0(ikz) + Act(it+1)
!      do ikr = 1,NKr
!        deps(ikr,ikz) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
!      end do
!    end do
!    deps_int = deps_int + 0.5d0*(dt*0.5d0)*deps
!!=== deps_int, deps ====
!
!  end do
!  
!  open(21,file='Act_jtz.out')
!  do it = 0,Nt
!    write(21,'(999e26.16e3)')dt*dble(it),Act(it),jtz(it),jtz_intra(it),jtz_inter(it),Act_dt2(it)
!  end do
!  close(21)

  if(Myrank == 0) write(*,*)  'This calculation is shutdown successfully!'
  call MPI_FINALIZE(ierr)

end subroutine multi_scale
