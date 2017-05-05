!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine excited_electron(nex,nex_v,nex_r,it)
  use global_variables
  implicit none
  real(8),intent(out) :: nex,nex_v,nex_r
  integer,intent(in) :: it
  real(8) :: lambda_v,lambda_c,theta_p,theta_m,eps_p,eps_m,dlambda
  real(8) :: Etz0,Etz1,alpha,ss,dEt_dt
  complex(8) :: zx,zy
  integer :: ikr,ikz
  complex(8) :: zeig_vec_v(2), zeig_vec_c(2)
  complex(8) :: zeig_vec_acc_v(2), zeig_vec_acc_c(2)
  real(8) :: lambda_acc_v,lambda_acc_c
  real(8) :: xx, xx_dot, ff, ff_dot, eta, eta_dot
  real(8) :: deps_dot(NKr,-NKz:NKz)

  Etz0 = -0.5d0*(Act(it+1)-Act(it-1))/dt
  dEt_dt = -(Act(it+1)-2d0*Act(it) + Act(it-1))/dt**2
  nex = 0d0
  nex_v = 0d0
  nex_r = 0d0

!$omp parallel

!=== deps_int, deps ====
!$omp do private(ikz, ikr)
  do ikz = -NKz,NKz
    kz(ikz) = kz0(ikz) + Act(it)*fact_intra
    do ikr = 1,NKr
      deps(ikr,ikz) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
      deps_dot(ikr,ikz) = -1d0/mass_r*kz(ikz)*Etz0*fact_intra
    end do
  end do
!=== deps_int, deps ====


!$omp do private(ikz,ikr,alpha,lambda_v,lambda_c,zx,zy,ss,zeig_vec_c,zeig_vec_v, &
!$omp& eta,eta_dot,xx,xx_dot,lambda_acc_v,lambda_acc_c,zeig_vec_acc_v,zeig_vec_acc_c) &
!$omp& reduction(+:nex, nex_v, nex_r)
  do ikz = -NKz,NKz
  do ikr = 1,NKr

    alpha = piz_vc*Etz0/deps(ikr,ikz)
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
    nex_v = nex_v+ abs(zCt(2,ikr,ikz))**2*kr(ikr)

! real carrier
    eta = 2d0*piz_vc*Etz0/deps(ikr,ikz)
    eta_dot = 2d0*piz_vc*(dEt_dt*deps(ikr,ikz) -2d0*Etz0*deps_dot(ikr,ikz) )/deps(ikr,ikz)**3
    xx = eta/(1d0 + sqrt(1d0 + eta**2))
    xx_dot = eta_dot/(1d0 + sqrt(1d0 + eta**2))/sqrt(1d0 + eta**2)

    alpha = xx_dot/(1d0 + xx**2)

    lambda_acc_v = (lambda_v + lambda_c) &
      - sqrt((lambda_c - lambda_v)**2 + 4d0*alpha**2)
    lambda_acc_v = 0.5d0 * lambda_acc_v

    lambda_acc_c = (lambda_v + lambda_c) &
      + sqrt((lambda_c - lambda_v)**2 + 4d0*alpha**2)
    lambda_acc_c = 0.5d0 * lambda_acc_c

    zx = alpha/(lambda_acc_c - lambda_v)
    zy = -alpha/(lambda_c - lambda_acc_v)
    
    zeig_vec_acc_v = zeig_vec_v + zy*zeig_vec_c
    zeig_vec_acc_c = zx*zeig_vec_v + zeig_vec_c

    ss = sum(abs(zeig_vec_acc_v)**2); zeig_vec_acc_v = zeig_vec_acc_v/sqrt(ss)
    ss = sum(abs(zeig_vec_acc_c)**2); zeig_vec_acc_c = zeig_vec_acc_c/sqrt(ss)

    zy = sum(conjg(zeig_vec_acc_c(:))*zCt(:,ikr,ikz))

    nex_r = nex_r+ abs(zy)**2*kr(ikr)


  end do
  end do
!$omp end parallel

  nex = nex*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
  nex_v = nex_v*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
  nex_r = nex_r*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 

  return
end subroutine excited_electron
