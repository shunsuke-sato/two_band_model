!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine dt_evolve_ms(Etz,im)
  use global_variables_ms
  implicit none
  integer, intent(in) :: im
  real(8), intent(in) :: Etz
  real(8) :: lambda,theta_p,theta_m,eps_p,eps_m
  real(8) :: ss
  complex(8) :: zsp,zsm
  integer :: ikr,ikz
  complex(8) :: zeig_vec_p(2), zeig_vec_m(2)

  do ikz = -NKz,NKz
  do ikr = 1,NKr

    lambda = deps_int(ikr,ikz,im)
    theta_p = lambda + 0.5d0*pi; theta_m = lambda - 0.5d0*pi
    zeig_vec_p(1) = 1d0/sqrt(2d0); zeig_vec_p(2) = exp(zi*theta_p)/sqrt(2d0)
    zeig_vec_m(1) = 1d0/sqrt(2d0); zeig_vec_m(2) = exp(zi*theta_m)/sqrt(2d0)

    eps_p = -piz_vc*Etz/deps(ikr,ikz,im); eps_m = -eps_p

    zsp = sum(conjg(zeig_vec_p(:))*zCt(:,ikr,ikz,im))
    zsm = sum(conjg(zeig_vec_m(:))*zCt(:,ikr,ikz,im))

    zsp = zsp * exp(-zI*eps_p*dt)
    zsm = zsm * exp(-zI*eps_m*dt)
    zCt(:,ikr,ikz,im) = zsp*zeig_vec_p(:) + zsm*zeig_vec_m(:)

  end do
  end do

  return
end subroutine dt_evolve_ms
