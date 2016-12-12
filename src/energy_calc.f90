!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine energy_calc(Eelec,Eelemag,Etot)
  use global_variables_ms
  implicit none
  real(8),intent(out) :: Eelec,Eelemag,Etot
  real(8) :: Eelec_l,Et,Bt
  integer :: ix,im,ikz,ikr

  Eelec_l = 0d0
  do im = Mx_s,Mx_e
    do ikz = -NKz,NKz
      do ikr = 1,NKr
        Eelec_l = Eelec_l + deps(ikr,ikz,im)*abs(zCt(2,ikr,ikz,im))**2*kr(ikr)
      end do
    end do
  end do
  Eelec_l = Eelec_l * 2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
  call MPI_ALLREDUCE(Eelec_l,Eelec,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  Eelemag = 0d0
  do ix = Nx_L+1,Nx_R-1
    Et = -0.5d0*(Az_new(ix)-Az_old(ix))/dt
    Bt = -0.5d0*(Az(ix+1)-Az(ix-1))*c_light/Hx
    Eelemag = Eelemag + (Et**2+Bt**2)
  end do
  
  Etot = Eelec + Eelemag
!!== Unit should be consistent between Eelec and Eelemag ==!!

end subroutine energy_calc
