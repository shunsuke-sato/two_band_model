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
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0,ev=Ry*2d0

!! Control parameters
! Band structure parameters
  integer,parameter :: N_PARABOLIC_BAND = 0
  integer,parameter :: N_NONPARABOLIC_BAND = 1
  integer,parameter :: N_COS_BAND = 2
  integer,parameter :: N_COS4_BAND = 3
  integer,parameter :: N_HBN_BAND = 4
! Pump-probe parameters
  integer,parameter :: N_COMBINED_PUMP_PROBE = 0
  integer,parameter :: N_DECOMPOSED_PUMP_PROBE = 1


! Material parameters
  integer :: NKr,NKz, NKx, NKy
  real(8) :: kr_max,kz_max,dkr,dkz, kx_max, ky_max, dkx,dky
  real(8),allocatable :: kr(:),kz0(:),kz(:),kx0(:),ky0(:),kx(:),ky(:)
  real(8),parameter :: eps_g = 12.1d0/(2d0*Ry)  !1.52d0/(2d0*Ry) 
  real(8),parameter :: mass_r = 0.4d0 !1d0/(1d0/0.57d0+1d0/0.067d0)
  real(8),parameter :: piz_vc = 0d0 !0.5d0*sqrt(eps_g/mass_r)
  real(8),parameter :: pix_vc = 0.35d0
  real(8),parameter :: piy_vc = 0d0
  real(8),parameter :: fact_intra = 1d0
  integer,parameter :: nband_type = N_PARABOLIC_BAND !N_NONPARABOLIC_BAND
  real(8),parameter :: band_width = 2d0/ev !1d0/(1d0/0.57d0+1d0/0.067d0)



! Time-propagation
  integer,parameter :: npump_probe_type = N_COMBINED_PUMP_PROBE
  integer :: Nt
  real(8) :: dt
  real(8),allocatable :: deps_int(:,:),deps(:,:),ddeps_dkz(:,:)
  complex(8),allocatable :: zCt(:,:,:)
  real(8),allocatable :: Act(:),Act_dt2(:),jtz(:),jtz_intra(:),jtz_inter(:)
  real(8),allocatable :: Act_pump(:),Act_probe(:)
  real(8),allocatable :: Act_xyz(:,:),Act_pump_xyz(:,:),Act_probe_xyz(:,:),jt_xyz(:,:)
  character(20) :: envelope_1,envelope_2

! Laser parameters
  real(8) :: E0_1,omega_1,tpulse_1,omega_ev_1,tpulse_fs_1,Iwcm2_1,CEP_2pi_1
  real(8) :: E0_2,omega_2,tpulse_2,omega_ev_2,tpulse_fs_2,Iwcm2_2,CEP_2pi_2
  real(8) :: Tdelay_fs,Tdelay

!! pure intraband field
  logical :: if_pure_intraband_fields_exist = .false.
  real(8),allocatable :: Act_intra(:)
  real(8) :: E0_static, T_duration_static
  real(8) :: E0_static_V_AA



end module global_variables
