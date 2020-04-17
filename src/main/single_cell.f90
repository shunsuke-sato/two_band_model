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
  real(8) :: jz_intra,jz_inter,Etz,Eex
  real(8) :: nex1st,nex2nd,nex3rd,nex4th
  real(8),allocatable :: nex0th_kz(:),nex1st_kz(:),nex2nd_kz(:)

  if(Nprocs /= 1)call err_finalize("Parallelization is not supported &
    for single-cell calculation.")

  read(*,*)kr_max,kz_max
  read(*,*)NKr, NKz
  read(*,*)Nt,dt
  read(*,*)envelope_1
  read(*,*)E0_1_V_m,omega_ev_1,tpulse_fs_1,CEP_2pi_1
  read(*,*)envelope_2
  read(*,*)E0_2_V_m,omega_ev_2,tpulse_fs_2,CEP_2pi_2
  read(*,*)Tdelay_fs
  read(*,*)E0_static_V_AA

  call preparation
  call input_Ac

  open(21,file="Eex_nex.out")
  open(22,file="nex_k.out")
  allocate(nex0th_kz(-NKz:NKz),nex1st_kz(-NKz:NKz),nex2nd_kz(-NKz:NKz))

  do it = 0,Nt
    write(*,*)'it=',it,'/',Nt
    call dt_evolve(it)

    
    call current(jz_intra,jz_inter)
    jtz_intra(it+1) = jz_intra; jtz_inter(it+1) = jz_inter
    jtz(it+1) = jtz_intra(it+1) + jtz_inter(it+1)

    if(mod(it,100) == 0 .or. it == Nt)then
       call energy(Eex)
       call excited_electron(nex1st,nex2nd,nex3rd,nex4th,it+1)
       write(21,"(999e26.16e3)")dt*(it+1),Eex,nex1st,nex2nd,nex3rd,nex4th
    end if

    if(mod(it,500) == 0 .or. it == Nt)then
      call excited_electron_k_resolved(nex0th_kz,nex1st_kz,nex2nd_kz,it)
      do ikz = -NKz, NKz,4
        write(22,"(9999e26.16e3)")dt*(it+1),kz0(ikz)+Act(it+1), &
          nex0th_kz(ikz),nex1st_kz(ikz),nex2nd_kz(ikz)
      end do
        write(22,*)
    end if


  end do

  close(21)
  close(22)

  open(21,file='Act_jtz.out')
  do it = 0,Nt
    write(21,'(999e26.16e3)')dt*dble(it),Act(it),jtz(it),jtz_intra(it),jtz_inter(it) &
      ,-0.5d0*(Act(it+1)-Act(it-1))/dt
  end do
  close(21)

end subroutine single_cell
