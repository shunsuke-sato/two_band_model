!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine excited_electron(nex,it)
  use global_variables
  implicit none
  real(8),intent(out) :: nex
  integer,intent(in) :: it
  real(8) :: lambda_v,lambda_c,theta_p,theta_m,eps_p,eps_m
  real(8) :: Etz0,Etz1,alpha,ss
  complex(8) :: zx,zy
  integer :: ikr,ikz
  complex(8) :: zeig_vec_v(2), zeig_vec_c(2)

  Etz0 = 0.5d0*(Act(it+1)-Act(it-1))/dt
  nex = 0d0

!$omp parallel

!=== deps_int, deps ====
!$omp do private(ikz, ikr)
  do ikz = -NKz,NKz
    kz(ikz) = kz0(ikz) + Act(it)*fact_intra
    do ikr = 1,NKr
      deps(ikr,ikz) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
    end do
  end do
!=== deps_int, deps ====




!$omp do private(ikz,ikr,alpha,lambda_v,lambda_c,zx,zy,ss,zeig_vec_c,zeig_vec_v) reduction(+:nex)
  do ikz = -NKz,NKz
  do ikr = 1,NKr

    alpha = -piz_vc*Etz0/deps(ikr,ikz)
    lambda_v = 0.5d0*(deps(ikr,ikz)-sqrt(deps(ikr,ikz)**2+4d0*alpha**2))
    lambda_c = 0.5d0*(deps(ikr,ikz)+sqrt(deps(ikr,ikz)**2+4d0*alpha**2))
    zx = zi*alpha/(deps(ikr,ikz)-lambda_v)
    zy = zi*alpha/lambda_c

    ss = 1d0/sqrt(1d0+abs(zx)**2)
    zeig_vec_v(1) = ss; zeig_vec_v(2) = zx*ss

    ss = 1d0/sqrt(abs(zy)**2+1d0)
    zeig_vec_c(1) = zy*ss; zeig_vec_c(2) = ss

    zx = sum(conjg(zeig_vec_v(:))*zCt(:,ikr,ikz))
    zy = sum(conjg(zeig_vec_c(:))*zCt(:,ikr,ikz))

    nex = nex+ abs(zy)**2*kr(ikr)
  end do
  end do
!$omp end parallel
  nex = nex*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 

  return
end subroutine excited_electron
