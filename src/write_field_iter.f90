!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine write_field_iter(it)
  use global_variables_ms
  implicit none
  integer :: it,ix
  real(8) :: Et
  character(50) :: citer,filename

  if(myrank == 0)then
    write(citer,"(I7.7)")it
    filename=trim(citer)//"_Et_Act.out"

    open(800,file=filename)
    do ix = Nx_L,Nx_R
      Et = -0.5d0*(Az_new(ix)-Az_old(ix))/dt
      write(800,"(999e26.16e3)")xn(ix),Et,Az(it)
    end do
    close(800)

  end if

end subroutine write_field_iter
