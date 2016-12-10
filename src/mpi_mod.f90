!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
! Two band model (3D)
module mpi_module

! MPI
  include 'mpif.h'
  integer :: Myrank,Nprocs,ierr

  contains
    subroutine err_finalize(err_message)
      implicit none
      character(*),intent(in) :: err_message
      if (Myrank == 0) then
        write(*,*) err_message
      endif
      call MPI_FINALIZE(ierr)
      stop
    end subroutine err_finalize
end module mpi_module
