!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
! Two band model (3D)
module global_variables
  use mpi_module
! mathematical parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zI = (0d0,1d0)
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0

! Material parameter
  integer :: NKr,NKz
  real(8) :: kr_max,kz_max,dkr,dkz
  real(8),allocatable :: kr(:),kz0(:),kz(:)
  real(8),parameter :: piz_vc = 1.0d0, mass_r = 1.4d0, eps_g = 42.2d0/(2d0*Ry) 

! Time-propagation
  integer :: Nt
  real(8) :: dt
  real(8),allocatable :: deps_int(:,:),deps(:,:)
  complex(8),allocatable :: zCt(:,:,:)
  real(8),allocatable :: Act(:),Act_dt2(:),jtz(:),jtz_intra(:),jtz_inter(:)

! Laser parameters
  real(8) :: E0_1,omega_1,tpulse_1,omega_ev_1,tpulse_fs_1,Iwcm2_1
  real(8) :: E0_2,omega_2,tpulse_2,omega_ev_2,tpulse_fs_2,Iwcm2_2
  real(8) :: Tdelay_fs,Tdelay



end module global_variables
