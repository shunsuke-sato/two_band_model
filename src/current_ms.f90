!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine current_ms
  use global_variables_ms
  implicit none
  integer :: ikr,ikz,im
  real(8) :: jz_intra,jz_inter
  real(8) :: jz_l(Mx)
  complex(8) :: zfact

  jz_l = 0d0
  do im = Mx_s,Mx_e
    jz_intra = 0d0
    jz_inter = 0d0
    do ikz = -NKz,NKz
      kz(ikz) = kz0(ikz) + Az_new(im)
      do ikr = 1,NKr
        jz_intra = jz_intra + kz(ikz)*abs(zCt(2,ikr,ikz,im))**2*kr(ikr)
        zfact = exp(-zI*deps_int(ikr,ikz,im))
        jz_inter = jz_inter &
          + 2d0*real(zfact*conjg(zCt(1,ikr,ikz,im))*zCt(2,ikr,ikz,im))*kr(ikr)
      end do
    end do

    jz_intra = jz_intra/mass_r*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
    jz_inter = jz_inter*piz_vc*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 

    jz_l(im) = jz_intra + jz_inter
  end do

  call MPI_ALLREDUCE(jz_l,jz,Mx,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  return
end subroutine current_ms
