!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine preparation_ms
  use global_variables_ms
  implicit none
  integer :: ix,Mx_ave, Mx_rem
  integer :: ikr,ikz,ik

  if(myrank == 0)write(*,"(A)")"== Preparation of the multi-scale calculation starts. =="  
  Mx_ave = Mx/Nprocs; Mx_rem = Mx - Mx_ave*Nprocs
  if(Myrank + 1 <= Mx_rem)then
    Mx_s = Myrank*(Mx_ave + 1) + 1
    Mx_e = (Myrank+1)*(Mx_ave + 1)
  else
    Mx_s = Mx_rem*(Mx_ave + 1) + (Myrank - Mx_rem)*Mx_ave + 1
    Mx_e = Mx_s + Mx_ave -1
  end if

  allocate(zCt(2,NKr,-Nkz:NKz,Mx_s:Mx_e),deps(NKr,-NKz:NKz,Mx_s:Mx_e) &
    ,deps_int(NKr,-NKz:NKz,Mx_s:Mx_e))
  allocate(kz0(-NKz:NKz),kz(-NKz:NKz),kr(NKr))
  zCt = 0d0; zCt(1,:,:,:) = 1d0
  deps = 0d0
  deps_int = 0d0
  

  dkr = kr_max/dble(NKr)
  dkz = kz_max/dble(NKz)

  do ikz = -NKz,NKz
    kz0(ikz) = dkz*dble(ikz)
  end do
  kz = kz0

  do ikr = 1,NKr
    kr(ikr) = dkr*dble(ikr)
  end do

  do ikz = -NKz,NKz
    kz(ikz) = kz0(ikz)
    do ikr = 1,NKr
      deps(ikr,ikz,:) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
    end do
  end do

  allocate(jz_store(0:Nt+1,Mx_s:Mx_e), Az_store(0:Nt+1,Mx_s:Mx_e))

  if(myrank == 0)write(*,"(A)")"== Preparation of the multi-scale calculation ends. =="  
  return
end subroutine preparation_ms
