!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine dt_evolve(it) ! Now coding
  use global_variables
  implicit none
  integer :: it
  real(8) :: lambda_v,lambda_c,theta_p,theta_m,eps_p,eps_m
  real(8) :: Etz0,Etz1,alpha,ss
  complex(8) :: zx,zy
  integer :: ikr,ikz
  complex(8) :: zeig_vec_v(2), zeig_vec_c(2)

  Etz0 = 0.5d0*(Act(it+1)-Act(it-1))/dt
  Etz1 = 0.5d0*(Act(it+2)-Act(it))/dt
  
!$omp parallel

!=== deps_int, deps ====
!$omp do private(ikz, ikr)
    do ikz = -NKz,NKz
      kz(ikz) = kz0(ikz) + Act(it)
      do ikr = 1,NKr
        deps(ikr,ikz) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
      end do
    end do
!=== deps_int, deps ====
  
  
!$omp do private(ikz,ikr,alpha,lambda_v,lambda_c,zx,zy,ss,zeig_vec_v,zeig_vec_c,Etz0)
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


    zx = zx * exp(-zI*lambda_v*0.5d0*dt)
    zy = zy * exp(-zI*lambda_c*0.5d0*dt)

    zCt(:,ikr,ikz) = zx*zeig_vec_v(:) + zy*zeig_vec_c(:)

  end do
  end do


!=== deps_int, deps ====
!$omp do private(ikz, ikr)
    do ikz = -NKz,NKz
      kz(ikz) = kz0(ikz) + Act(it+1)
      do ikr = 1,NKr
        deps(ikr,ikz) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
      end do
    end do
!=== deps_int, deps ====
  
  
!$omp do private(ikz,ikr,alpha,lambda_v,lambda_c,zx,zy,ss,zeig_vec_v,zeig_vec_c,Etz1)
  do ikz = -NKz,NKz
  do ikr = 1,NKr

    alpha = -piz_vc*Etz1/deps(ikr,ikz)
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


    zx = zx * exp(-zI*lambda_v*0.5d0*dt)
    zy = zy * exp(-zI*lambda_c*0.5d0*dt)

    zCt(:,ikr,ikz) = zx*zeig_vec_v(:) + zy*zeig_vec_c(:)

  end do
  end do

!$omp end parallel

  return
end subroutine dt_evolve
