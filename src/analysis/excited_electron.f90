!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine excited_electron(nex1st,nex2nd,nex3rd,nex4th,it)
  use global_variables
  implicit none
  real(8),intent(out) :: nex1st,nex2nd,nex3rd,nex4th
  integer,intent(in) :: it
  real(8) :: theta_p,theta_m,eps_p,eps_m,dlambda
  real(8) :: Etz0,Etz1,alpha,ss,dEt_dt,d2Et_dt2
  complex(8) :: zx,zy
  integer :: ikr,ikz
  complex(8) :: zeig_vec_2nd_v(2), zeig_vec_2nd_c(2)
  complex(8) :: zeig_vec_3rd_v(2), zeig_vec_3rd_c(2)
  complex(8) :: zeig_vec_4th_v(2), zeig_vec_4th_c(2)
  real(8) :: lambda_2nd_v, lambda_2nd_c
  real(8) :: lambda_3rd_v, lambda_3rd_c
  real(8) :: lambda_4th_v, lambda_4th_c
  real(8) :: xx, xx_dot, xx_dot2, ff, ff_dot, eta, eta_dot, eta_dot2
  real(8) :: gamma, gamma_dot, gamma_dot2
  real(8) :: gamma_2nd, deps_2nd, eta_2nd, eta_2nd_dot
  real(8) :: yy
  real(8) :: deps_dot(NKr,-NKz:NKz),deps_dot2(NKr,-NKz:NKz)


  nex1st = 0d0
  nex2nd = 0d0
  nex3rd = 0d0
  nex4th = 0d0

  Etz0 = -0.5d0*(Act(it+1)-Act(it-1))/dt
  dEt_dt = -(Act(it+1)-2d0*Act(it) + Act(it-1))/dt**2
  d2Et_dt2 = -(0.5d0*Act(it+2)-Act(it+1)+Act(it-1)-0.5d0*Act(it-2))/dt**3



!$omp parallel 

!=== deps_int, deps ====
!$omp do private(ikz, ikr)
  do ikz = -NKz,NKz
    kz(ikz) = kz0(ikz) + Act(it)*fact_intra
    do ikr = 1,NKr
      deps(ikr,ikz) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
      deps_dot(ikr,ikz) = -1d0/mass_r*kz(ikz)*Etz0*fact_intra
      deps_dot2(ikr,ikz) = 1d0/mass_r*(Etz0**2-kz(ikz)*dEt_dt)*fact_intra
    end do
  end do
!=== deps_int, deps ====


!$omp do private(ikz,ikr,gamma,gamma_dot,gamma_dot2, &
!$omp& eta,eta_dot,eta_dot2,xx,xx_dot,xx_dot2, &
!$omp& lambda_2nd_v,lambda_2nd_c,ss,zeig_vec_2nd_v,zeig_vec_2nd_c, &
!$omp& zx,zy,gamma_2nd,deps_2nd,eta_2nd,yy,& 
!$omp& lambda_3rd_v,lambda_3rd_c,zeig_vec_3rd_v,zeig_vec_3rd_c) &
!$omp& reduction(+:nex1st,nex2nd,nex3rd)
  do ikz = -NKz,NKz
  do ikr = 1,NKr

    gamma      = piz_vc*Etz0/deps(ikr,ikz)
    gamma_dot  = piz_vc*(dEt_dt*deps(ikr,ikz)-Etz0*deps_dot(ikr,ikz))/deps(ikr,ikz)**2
    gamma_dot2 = piz_vc*(d2Et_dt2*deps(ikr,ikz)-Etz0*deps_dot2(ikr,ikz))/deps(ikr,ikz)**2 &
                -gamma_dot*2d0*deps_dot(ikr,ikz)/deps(ikr,ikz)

    eta      = 2d0*gamma/deps(ikr,ikz)
    eta_dot  = 2d0*(gamma_dot*deps(ikr,ikz)-gamma*deps_dot(ikr,ikz))/deps(ikr,ikz)**2
    eta_dot2 = 2d0*(gamma_dot2*deps(ikr,ikz)-gamma*deps_dot2(ikr,ikz))/deps(ikr,ikz)**2 &
              -eta_dot*2d0*deps_dot(ikr,ikz)/deps(ikr,ikz)

    xx      = eta/(1d0+sqrt(1d0+eta**2))
    xx_dot  = eta_dot /((1d0+sqrt(1d0+eta**2))*sqrt(1d0+eta**2))
    xx_dot2 = eta_dot2/((1d0+sqrt(1d0+eta**2))*sqrt(1d0+eta**2))  &
             -eta*xx_dot**2*(2d0+1d0/sqrt(1d0+eta**2)) 

! nex_1st
    nex1st = nex1st+ abs(zCt(2,ikr,ikz))**2*kr(ikr)

! nex_2nd
    lambda_2nd_v = 0.5d0*(deps(ikr,ikz)-sqrt(deps(ikr,ikz)**2+4d0*gamma**2))
    lambda_2nd_c = 0.5d0*(deps(ikr,ikz)+sqrt(deps(ikr,ikz)**2+4d0*gamma**2))

    ss = 1d0/sqrt(1d0+xx**2)
    zeig_vec_2nd_v(1) = ss      ; zeig_vec_2nd_v(2) = zI*xx*ss
    zeig_vec_2nd_c(1) = zI*xx*ss; zeig_vec_2nd_c(2) = ss

    zx = sum(conjg(zeig_vec_2nd_v(:))*zCt(:,ikr,ikz))
    zy = sum(conjg(zeig_vec_2nd_c(:))*zCt(:,ikr,ikz))

    nex2nd = nex2nd+ abs(zy)**2*kr(ikr)

! nex_3rd
    gamma_2nd = xx_dot/(1d0+xx**2)
    deps_2nd  = sqrt(deps(ikr,ikz)**2 + 4d0*gamma**2)
    eta_2nd   = 2d0*gamma_2nd/deps_2nd
    yy        = -eta_2nd/(1d0+sqrt(1d0 + eta_2nd**2))

    lambda_3rd_v = (lambda_2nd_v + lambda_2nd_c) &
      - sqrt((lambda_2nd_c - lambda_2nd_v)**2 + 4d0*gamma_2nd**2)
    lambda_3rd_v = 0.5d0 * lambda_3rd_v

    lambda_3rd_c = (lambda_2nd_v + lambda_2nd_c) &
      + sqrt((lambda_2nd_c - lambda_2nd_v)**2 + 4d0*gamma_2nd**2)
    lambda_3rd_c = 0.5d0 * lambda_3rd_c

    ss = 1d0/sqrt(1d0 + yy**2)
    
    zeig_vec_3rd_v =  zeig_vec_2nd_v*ss       + zeig_vec_2nd_c*ss*yy
    zeig_vec_3rd_c = -zx*zeig_vec_2nd_v*ss*yy + zeig_vec_2nd_c*ss

    zy = sum(conjg(zeig_vec_3rd_c(:))*zCt(:,ikr,ikz))

    nex3rd = nex3rd+ abs(zy)**2*kr(ikr)

  end do
  end do
!$omp end parallel

  nex1st = nex1st*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
  nex2nd = nex2nd*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
  nex3rd = nex3rd*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
  nex4th = nex4th*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 

  return
end subroutine excited_electron
