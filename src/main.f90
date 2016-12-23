!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
program main
  use mpi_module
  implicit none
  character(50) :: calc_mode

  call MPI_init(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Myrank,ierr)
  Time_start=MPI_WTIME()

  if(myrank == 0)read(*,*)calc_mode
  call MPI_BCAST(calc_mode,50,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)  

  select case(calc_mode)
  case("single-cell","SC","sc")
    call single_cell
  case("multi-scale","MS","ms")
    call multi_scale
  case default
    write(*,"(A,2x,A)")"Invalid calc_mode",calc_mode
    stop
  end select
      
  Time_now=MPI_WTIME()
  if (Myrank == 0 ) write(*,"(A,2x,e16.6e3,A)") 'Total time =',(Time_now-Time_start),'sec'

end program main
