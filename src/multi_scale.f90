!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine multi_scale
  use global_variables_ms
  implicit none
  integer :: it

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
  call input_Ac_ms
  jz_store(0,Mx_s:Mx_e)=0d0
  Az_store(0,Mx_s:Mx_e)=Az_new(Mx_s:Mx_e)

  if(myrank == 0)then
    open(20,file="Act_vac.out")
    write(20,"(A)")"# Time (a.u.), Act (front surf.), Act (rear surf.)"
    write(20,"(999e26.16e3)")0d0,Az_new(0),Az_new(Mx+1)
  end if

  do it = 0,Nt
    write(*,*)'it=',it,'/',Nt

    call prop_elec_ms
    call current_ms
    call prop_Ac_ms
    if(myrank == 0)write(20,"(999e26.16e3)")dt*(it+1),Az_new(0),Az_new(Mx+1)
    jz_store(it+1,Mx_s:Mx_e)=jz(Mx_s:Mx_e)
    Az_store(it+1,Mx_s:Mx_e)=Az_new(Mx_s:Mx_e)
    if(mod(it,1000) == 0)call write_field_iter(it)
  end do

  if(Myrank == 0)close(20)
  if(Myrank == 0) write(*,*)  'This calculation is shutdown successfully!'
  call MPI_FINALIZE(ierr)

end subroutine multi_scale
