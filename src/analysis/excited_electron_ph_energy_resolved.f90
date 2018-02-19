!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine excited_electron_ph_energy_resolved(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8) :: nex1st_krz(NKr,-NKz:NKz)
  real(8) :: nex2nd_krz(NKr,-NKz:NKz)
  real(8) :: nex3rd_krz(NKr,-NKz:NKz)
  real(8) :: nex4th_krz(NKr,-NKz:NKz)
  real(8),allocatable :: nph_dns(:)
  real(8),parameter :: dEph = 0.01d0/(2d0*Ry)
  real(8) :: Eph_max
  integer :: NEph, ieph
  character(len=512) :: cit, cfilename

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


  nex1st_krz = 0d0
  nex2nd_krz = 0d0
  nex3rd_krz = 0d0
  nex4th_krz = 0d0

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
    nex1st_krz(ikr,ikz) = abs(zCt(2,ikr,ikz))**2*kr(ikr)

! nex_2nd
    lambda_2nd_v = 0.5d0*(deps(ikr,ikz)-sqrt(deps(ikr,ikz)**2+4d0*gamma**2))
    lambda_2nd_c = 0.5d0*(deps(ikr,ikz)+sqrt(deps(ikr,ikz)**2+4d0*gamma**2))

    ss = 1d0/sqrt(1d0+xx**2)
    zeig_vec_2nd_v(1) = ss      ; zeig_vec_2nd_v(2) = zI*xx*ss
    zeig_vec_2nd_c(1) = zI*xx*ss; zeig_vec_2nd_c(2) = ss

    zx = sum(conjg(zeig_vec_2nd_v(:))*zCt(:,ikr,ikz))
    zy = sum(conjg(zeig_vec_2nd_c(:))*zCt(:,ikr,ikz))

    nex2nd_krz(ikr,ikz) = abs(zy)**2*kr(ikr)

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

    nex3rd_krz(ikr,ikz) = abs(zy)**2*kr(ikr)

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

    nex4th_krz(ikr,ikz) = abs(zy)**2*kr(ikr)

  end do
  end do
!$omp end parallel

  nex1st_krz = nex1st_krz*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
  nex2nd_krz = nex2nd_krz*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
  nex3rd_krz = nex3rd_krz*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 
  nex4th_krz = nex4th_krz*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 

  Eph_max = maxval(deps)
  NEph = aint(Eph_max/dEph) + 1
  allocate(nph_dns(0:NEph))
  nph_dns = 0d0
  do ikz = -NKz,NKz
     do ikr = 1,NKr
        ieph = aint(deps(ikr,ikz)/dEph)
        nph_dns(ieph) = nph_dns(ieph) + nex1st_krz(ikr,ikz)

     end do
  end do

  write(cit,"(I9.9)")it
  cfilename = "nph_dist_"//trim(cit)//".out"
  open(41,file=cfilename)
  do ieph = 0, NEph

     write(41,"(999e26.16e3)")ieph*dEph, nph_dns(ieph)

  end do
  close(41)

  return
end subroutine excited_electron_ph_energy_resolved
