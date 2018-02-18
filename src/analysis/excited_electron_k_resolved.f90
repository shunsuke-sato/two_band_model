!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine excited_electron_k_resolved(nex0th_kz,nex1st_kz &
                                                ,nex2nd_kz &
                                                ,nex3rd_kz &
                                                ,nex4th_kz &
                                                ,it)
  use global_variables
  implicit none
  real(8),intent(out) :: nex0th_kz(-NKz:NKz)
  real(8),intent(out) :: nex1st_kz(-NKz:NKz)
  real(8),intent(out) :: nex2nd_kz(-NKz:NKz)
  real(8),intent(out) :: nex3rd_kz(-NKz:NKz)
  real(8),intent(out) :: nex4th_kz(-NKz:NKz)
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
  real(8) :: gamma_2nd, deps_2nd, gamma_2nd_dot, deps_2nd_dot
  real(8) :: eta_2nd, eta_2nd_dot
  real(8) :: yy, yy_dot
  real(8) :: gamma_3rd, eta_3rd, deps_3rd, zz
  real(8) :: deps_dot(NKr,-NKz:NKz),deps_dot2(NKr,-NKz:NKz)

  nex0th_kz = 0d0
  nex1st_kz = 0d0
  nex2nd_kz = 0d0  
  nex3rd_kz = 0d0  
  nex4th_kz = 0d0  


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
!$omp& lambda_3rd_v,lambda_3rd_c,zeig_vec_3rd_v,zeig_vec_3rd_c,&
!$omp& gamma_2nd_dot,deps_2nd_dot,eta_2nd_dot,yy_dot, &
!$omp& gamma_3rd, deps_3rd, eta_3rd, zz, &
!$omp& zeig_vec_4th_v, zeig_vec_4th_c) 
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
    nex1st_kz(ikz) = nex1st_kz(ikz) + abs(zCt(2,ikr,ikz))**2*kr(ikr)

! nex_2nd
    lambda_2nd_v = 0.5d0*(deps(ikr,ikz)-sqrt(deps(ikr,ikz)**2+4d0*gamma**2))
    lambda_2nd_c = 0.5d0*(deps(ikr,ikz)+sqrt(deps(ikr,ikz)**2+4d0*gamma**2))

    ss = 1d0/sqrt(1d0+xx**2)
    zeig_vec_2nd_v(1) = ss      ; zeig_vec_2nd_v(2) = zI*xx*ss
    zeig_vec_2nd_c(1) = zI*xx*ss; zeig_vec_2nd_c(2) = ss

    zx = sum(conjg(zeig_vec_2nd_v(:))*zCt(:,ikr,ikz))
    zy = sum(conjg(zeig_vec_2nd_c(:))*zCt(:,ikr,ikz))

    nex2nd_kz(ikz) = nex2nd_kz(ikz) + abs(zy)**2*kr(ikr)

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

    nex3rd_kz(ikz) = nex3rd_kz(ikz) + abs(zy)**2*kr(ikr)

! nex_4th
    gamma_2nd_dot = xx_dot2/(1d0+xx**2)-2d0*xx*gamma_2nd**2
    deps_2nd_dot  = (deps(ikr,ikz)*deps_dot(ikr,ikz)+4d0*gamma*gamma_dot)/deps_2nd
    eta_2nd_dot = 2d0*(gamma_2nd_dot*deps_2nd-gamma_2nd*deps_2nd_dot)/deps_2nd**2
    yy_dot = -eta_2nd_dot/((1d0+sqrt(1d0+eta_2nd**2))*sqrt(1d0+eta_2nd**2))
    gamma_3rd = yy_dot/(1d0+yy**2)
    deps_3rd  = lambda_3rd_c - lambda_3rd_v
    eta_3rd   = 2d0*gamma_3rd/deps_3rd
    zz = eta_3rd/(1d0+sqrt(1d0+eta_3rd**2))

    ss = 1d0/sqrt(1d0 + zz**2)

    zeig_vec_4th_v =    zeig_vec_3rd_v*ss       + zI*zeig_vec_3rd_c*ss*zz
    zeig_vec_4th_c = zI*zeig_vec_3rd_v*ss*zz    +    zeig_vec_3rd_c*ss

    zy = sum(conjg(zeig_vec_4th_c(:))*zCt(:,ikr,ikz))

    nex4th_kz(ikz) = nex4th_kz(ikz) + abs(zy)**2*kr(ikr)

  end do
  end do
!$omp end parallel

  nex0th_kz = nex0th_kz*2d0/((2d0*pi)**3)*(2d0*pi*dkr) 
  nex1st_kz = nex1st_kz*2d0/((2d0*pi)**3)*(2d0*pi*dkr) 
  nex2nd_kz = nex2nd_kz*2d0/((2d0*pi)**3)*(2d0*pi*dkr) 
  nex3rd_kz = nex3rd_kz*2d0/((2d0*pi)**3)*(2d0*pi*dkr) 
  nex4th_kz = nex4th_kz*2d0/((2d0*pi)**3)*(2d0*pi*dkr) 

end subroutine excited_electron_k_resolved

!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine excited_electron_k_resolved_old(nex0th_kz,nex1st_kz,nex2nd_kz,it)
  use global_variables
  implicit none
  real(8),intent(out) :: nex0th_kz(-NKz:NKz)
  real(8),intent(out) :: nex1st_kz(-NKz:NKz)
  real(8),intent(out) :: nex2nd_kz(-NKz:NKz)
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
  nex0th_kz = 0d0
  nex1st_kz = 0d0
  nex2nd_kz = 0d0

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
!$omp& eta,eta_dot,xx,xx_dot,lambda_acc_v,lambda_acc_c,zeig_vec_acc_v,zeig_vec_acc_c)
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

    nex1st_kz(ikz) = nex1st_kz(ikz) + abs(zy)**2*kr(ikr)
    nex0th_kz(ikz) = nex0th_kz(ikz) + abs(zCt(2,ikr,ikz))**2*kr(ikr)

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

    nex2nd_kz(ikz) = nex2nd_kz(ikz) + abs(zy)**2*kr(ikr)

  end do
  end do
!$omp end parallel

  nex0th_kz = nex0th_kz*2d0/((2d0*pi)**3)*(2d0*pi*dkr) 
  nex1st_kz = nex1st_kz*2d0/((2d0*pi)**3)*(2d0*pi*dkr) 
  nex2nd_kz = nex2nd_kz*2d0/((2d0*pi)**3)*(2d0*pi*dkr) 

end subroutine excited_electron_k_resolved_old
