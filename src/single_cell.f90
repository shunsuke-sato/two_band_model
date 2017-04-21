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
  real(8) :: jz_intra,jz_inter,Etz,Eex,nex


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

  call preparation
  call input_Ac

  open(21,file="Eex_nex.out")

  do it = 0,Nt
    write(*,*)'it=',it,'/',Nt
    call dt_evolve(it)

    
    call current(jz_intra,jz_inter)
    jtz_intra(it+1) = jz_intra; jtz_inter(it+1) = jz_inter
    jtz(it+1) = jtz_intra(it+1) + jtz_inter(it+1)

    if(mod(it,100) == 0 .or. it == Nt)then
       call energy(Eex)
       call excited_electron(nex,it+1)
       write(21,"(999e26.16e3)")dt*(it+1),Eex,nex
    end if

  end do

  close(21)
  open(21,file='Act_jtz.out')
  do it = 0,Nt
    write(21,'(999e26.16e3)')dt*dble(it),Act(it),jtz(it),jtz_intra(it),jtz_inter(it)
  end do
  close(21)

end subroutine single_cell
