!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine prop_elec_ms
  use global_variables_ms
  implicit none
  integer :: im,ikz,ikr
  real(8) :: a,b,c,Act_t,Etz

  do im = Mx_s,Mx_e

    c = Az(im)
    b = (Az_new(im)-Az_old(im))/(2d0*dt)
    a = (Az_new(im) - 2d0*Az(im) + Az_old(im))/(dt**2)
    
!=== deps_int, deps ====
    Act_t = Az(im)
    do ikz = -NKz,NKz
      kz(ikz) = kz0(ikz) + Act_t
      do ikr = 1,NKr
        deps(ikr,ikz,im) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
      end do
    end do
    deps_int(:,:,im) = deps_int(:,:,im) + 0.5d0*(dt*0.5d0)*deps(:,:,im)

    Act_t = 0.5d0*a*(dt/2d0)**2 + b*dt/2d0 + c
    do ikz = -NKz,NKz
      kz(ikz) = kz0(ikz) + Act_t
      do ikr = 1,NKr
        deps(ikr,ikz,im) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
      end do
    end do
    deps_int(:,:,im) = deps_int(:,:,im) + 0.5d0*(dt*0.5d0)*deps(:,:,im)
!=== deps_int, deps ====

    Etz = -(Az_new(im)-Az(im))/dt
    call dt_evolve_ms(Etz,im)

!=== deps_int, deps ====
    Act_t = 0.5d0*a*(dt/2d0)**2 + b*dt/2d0 + c
    do ikz = -NKz,NKz
      kz(ikz) = kz0(ikz) + Act_t
      do ikr = 1,NKr
        deps(ikr,ikz,im) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
      end do
    end do
    deps_int(:,:,im) = deps_int(:,:,im) + 0.5d0*(dt*0.5d0)*deps(:,:,im)

    Act_t = 0.5d0*a*dt**2 + b*dt + c
    do ikz = -NKz,NKz
      kz(ikz) = kz0(ikz) + Act_t
      do ikr = 1,NKr
        deps(ikr,ikz,im) = eps_g + 0.5d0/mass_r*(kr(ikr)**2+kz(ikz)**2)
      end do
    end do
    deps_int(:,:,im) = deps_int(:,:,im) + 0.5d0*(dt*0.5d0)*deps(:,:,im)
!=== deps_int, deps ====
  end do

  return
end subroutine prop_elec_ms
