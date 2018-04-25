module global_variables
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0,1d0)
  real(8),parameter :: fs = 0.024189d0
  real(8),parameter :: ev = 27.2114d0
  real(8),parameter :: b_a = 0.529d0
  real(8),parameter :: clight = 137.035999139d0

! material parameters
  real(8),parameter :: tau_z = 1d0, delta_gap = 0d0
! Fermi velocity of graphene PRL 101, 226405 (2008).
  real(8),parameter :: velocity = clight*1.12d6/299792458d0


  real(8),parameter :: kx0=1d0, ky0=1d0
  real(8),parameter :: qx0=1d0, qy0=1d0
  real(8),parameter :: omega_phonon =1d0

  integer,parameter :: num_elec = 1
  integer,parameter :: nph_cutoff = 1

  integer,parameter :: ndim_elec, ndim_phonon, ndim_tot
  real(8),allocatable :: dHam(:,:)
  complex(8),allocatable :: zpsi(:,:)

  real(8),parameter :: gc_vv = 1d0, gc_cc = 1d0, gc_vc = 1d0
  

end module global_variables
!----------------------------------------------------------------------------------------!
program main
  use global_variables
  implicit none


  select case(num_elec)
  case(1)
    call one_electron_system
  case default
    write(*,*)"num_elec=",num_elec,"is invalid."
  end select


end program main
!----------------------------------------------------------------------------------------!
subroutine one_electron_system
  use global_variables
  implicit none

  call init_one_electron_system

end subroutine one_electron_system
!----------------------------------------------------------------------------------------!
subroutine init_one_electron_system
  use global_variables
  implicit none
  integer :: idim_elec, jdim_elec, idim_ph, jdim_ph
  integer :: idim, jdim
!LAPACK
  integer :: lwork
  real(8),allocatable :: work_lp(:),amat(:,:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info




  ndim_elec    = 4
  ndim_phonon  = 1 + nph_cutoff
  ndim_tot     = ndim_elec * ndim_phonon

  allocate(dHam(ndim_tot,ndim_tot))
  allocate(zpsi(ndim_elec,0:nph_cutoff))

! 1, |v_K>
! 2, |c_K>
! 3, |v_K+Q>
! 4, |c_K+Q>

  dHam = 0d0
  do idim_ph = 0,nph_cutoff
    do idim_elec = 1,4
      idim = idim_elec + 4*idim_ph


      do jdim_ph = 0,nph_cutoff
        do jdim_elec = 1,4
          jdim = jdim_elec + 4*jdim_ph

! diagonal part
          if(idim == jdim)then
            dHam(idim,jdim) = dHam(idim,jdim) + omega_phonon*idim_ph
            
            select case(idim_elec)
            case(1)
              dHam(idim,jdim) = dHam(idim,jdim) - velocity*sqrt(kx0**2+ky0**2)
            case(2)
              dHam(idim,jdim) = dHam(idim,jdim) + velocity*sqrt(kx0**2+ky0**2)
            case(3)
              dHam(idim,jdim) = dHam(idim,jdim) - velocity*sqrt((kx0+qx0)**2+(ky0+qy0)**2)
            case(4)
              dHam(idim,jdim) = dHam(idim,jdim) + velocity*sqrt((kx0+qx0)**2+(ky0+qy0)**2)
            case default
              stop "Error in construction of initial Hamiltonian!!"
            end select


          end if


! electron-phonon interactions          
          if(idim_elec*jdim_elec == 3)then ! 1*3
            if(abs(idim_ph - jdim_ph) == 1)then
              dHam(idim,jdim) = dHam(idim,jdim) + gc_vv*sqrt(dble(max(idim_ph,jdim_ph)))
            end if
          end if

          if(idim_elec*jdim_elec == 8)then ! 2*4
            if(abs(idim_ph - jdim_ph) == 1)then
              dHam(idim,jdim) = dHam(idim,jdim) + gc_cc*sqrt(dble(max(idim_ph,jdim_ph)))
            end if
          end if

        end do
      end do
    end do
  end do

Call dsyev('Vectors', 'Upper', n, a, lda, w, work, lwork, info)

end subroutine init_one_electron_system
