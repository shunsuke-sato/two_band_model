!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine dt_evolve_2d(it)
  use global_variables_2d
  implicit none
  integer,intent(in) :: it
  real(8) :: Etx, Ety, Etz, Acx, Acy, Acz
  real(8) :: lambda_v,lambda_c,theta_p,theta_m,eps_p,eps_m
  real(8) :: alpha,ss,pi_dot_E
  complex(8) :: zx,zy
  integer :: ikx,iky
  complex(8) :: zeig_vec_v(2), zeig_vec_c(2)

  select case(npump_probe_type)
  case(N_COMBINED_PUMP_PROBE)
    Etx = -(Act_xyz(it+1,1)-Act_xyz(it,1))/dt
    Ety = -(Act_xyz(it+1,2)-Act_xyz(it,2))/dt
    Etz = -(Act_xyz(it+1,3)-Act_xyz(it,3))/dt
    Acx = 0.5d0*(Act_xyz(it+1,1) + Act_xyz(it,1))
    Acy = 0.5d0*(Act_xyz(it+1,2) + Act_xyz(it,2))
    Acz = 0.5d0*(Act_xyz(it+1,3) + Act_xyz(it,3))
  case default
    write(*,"(A,2x,A)")"Invalid npump_probe_type",npump_probe_type
  end select


!$omp parallel do private(ikx)
  do ikx = -NKx,NKx
    kx(ikx) = kx0(ikx) + Acx*fact_intra
  end do

!$omp parallel do private(iky)
  do iky = -NKy,NKy
    ky(iky) = ky0(iky) + Acy*fact_intra
  end do

  call set_deps_2d
  
  pi_dot_E = pix_vc*Etx + piy_vc*Ety + piz_vc*Etz

  if(mass_r < 0d0)then
    call dt_evolve_2d_nagative_mass(pi_dot_E)
    return
  end if

!$omp parallel do private(ikx,iky,alpha,lambda_v,lambda_c,zx,zy,ss,zeig_vec_v,zeig_vec_c) collapse(2)
  do ikx = -NKx,NKx
  do iky = -NKy,NKy

    alpha = pi_dot_E/deps(ikx,iky)
    lambda_v = 0.5d0*(deps(ikx,iky)-sqrt(deps(ikx,iky)**2+4d0*alpha**2))
    lambda_c = 0.5d0*(deps(ikx,iky)+sqrt(deps(ikx,iky)**2+4d0*alpha**2))
    zx = zi*alpha/(deps(ikx,iky)-lambda_v)
    zy = zi*alpha/lambda_c

    ss = 1d0/sqrt(1d0+abs(zx)**2)
    zeig_vec_v(1) = ss; zeig_vec_v(2) = zx*ss

    ss = 1d0/sqrt(abs(zy)**2+1d0)
    zeig_vec_c(1) = zy*ss; zeig_vec_c(2) = ss

    zx = sum(conjg(zeig_vec_v(:))*zCt(:,ikx,iky))
    zy = sum(conjg(zeig_vec_c(:))*zCt(:,ikx,iky))


    zx = zx * exp(-zI*lambda_v*dt)
    zy = zy * exp(-zI*lambda_c*dt)

    zCt(:,ikx,iky) = zx*zeig_vec_v(:) + zy*zeig_vec_c(:)

  end do
  end do


  return
end subroutine dt_evolve_2d
!-----------------------------------------------------------------------------------------
subroutine dt_evolve_2d_nagative_mass(pi_dot_E)
  use global_variables_2d
  implicit none
  real(8),intent(in) :: pi_dot_E
  real(8),parameter :: Ecut = 3d0/ev
  real(8) :: lambda_v,lambda_c,theta_p,theta_m,eps_p,eps_m
  real(8) :: alpha,ss
  complex(8) :: zx,zy
  integer :: ikx,iky
  complex(8) :: zeig_vec_v(2), zeig_vec_c(2)
  real(8)  :: d_dot_E(-NKx:NKx,-NKy:NKy)

  do ikx = -NKx,NKx
  do iky = -NKy,NKy
    if(deps(ikx,iky)>Ecut)then
      d_dot_E(ikx,iky) = pi_dot_E/deps(ikx,iky)
    else
      d_dot_E(ikx,iky) = 0d0
    end if
  end do
  end do


!$omp parallel do private(ikx,iky,alpha,lambda_v,lambda_c,zx,zy,ss,zeig_vec_v,zeig_vec_c) collapse(2)
  do ikx = -NKx,NKx
  do iky = -NKy,NKy

!    alpha = pi_dot_E/deps(ikx,iky)
    alpha = d_dot_E(ikx,iky)
    lambda_v = 0.5d0*(deps(ikx,iky)-sqrt(deps(ikx,iky)**2+4d0*alpha**2))
    lambda_c = 0.5d0*(deps(ikx,iky)+sqrt(deps(ikx,iky)**2+4d0*alpha**2))
    zx = zi*alpha/(deps(ikx,iky)-lambda_v)
    zy = zi*alpha/lambda_c

    ss = 1d0/sqrt(1d0+abs(zx)**2)
    zeig_vec_v(1) = ss; zeig_vec_v(2) = zx*ss

    ss = 1d0/sqrt(abs(zy)**2+1d0)
    zeig_vec_c(1) = zy*ss; zeig_vec_c(2) = ss

    zx = sum(conjg(zeig_vec_v(:))*zCt(:,ikx,iky))
    zy = sum(conjg(zeig_vec_c(:))*zCt(:,ikx,iky))


    zx = zx * exp(-zI*lambda_v*dt)
    zy = zy * exp(-zI*lambda_c*dt)

    zCt(:,ikx,iky) = zx*zeig_vec_v(:) + zy*zeig_vec_c(:)

  end do
  end do


end subroutine dt_evolve_2d_nagative_mass
