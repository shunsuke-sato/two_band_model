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

! Control parameters
  integer,parameter :: N_PARABOLIC_BAND = 0
  integer,parameter :: N_NONPARABOLIC_BAND = 1

! Material parameters
  integer :: NKr,NKz
  real(8) :: kr_max,kz_max,dkr,dkz
  real(8),allocatable :: kr(:),kz0(:),kz(:)
  real(8),parameter :: eps_g = 9d0/(2d0*Ry)  !1.52d0/(2d0*Ry) 
  real(8),parameter :: mass_r = 0.4d0 !1d0/(1d0/0.57d0+1d0/0.067d0)
  real(8),parameter :: piz_vc = 0.5d0*sqrt(eps_g/mass_r)
  real(8),parameter :: fact_intra = 1d0
  integer,parameter :: nband_type = N_PARABOLIC_BAND !N_NONPARABOLIC_BAND



! Time-propagation
  integer :: Nt
  real(8) :: dt
  real(8),allocatable :: deps_int(:,:),deps(:,:),ddeps_dkz(:,:)
  complex(8),allocatable :: zCt(:,:,:)
  real(8),allocatable :: Act(:),Act_dt2(:),jtz(:),jtz_intra(:),jtz_inter(:)
  character(20) :: envelope_1,envelope_2

! Laser parameters
  real(8) :: E0_1,omega_1,tpulse_1,omega_ev_1,tpulse_fs_1,Iwcm2_1
  real(8) :: E0_2,omega_2,tpulse_2,omega_ev_2,tpulse_fs_2,Iwcm2_2
  real(8) :: Tdelay_fs,Tdelay



end module global_variables
