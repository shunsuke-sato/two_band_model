!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine write_results(it)
  use global_variables_ms
  implicit none
  integer :: it,im
  character(50) :: citer,filename

  do im = Mx_s,Mx_e
    write(citer,"(I7.7)")im
    filename="M"//trim(citer)//"_AcJ.out"
    open(800,file=filename)
    do it = 0,Nt+1
      write(800,"(999e26.16e3)")dt*dble(it),Az_store(it,im),Jz_store(it,im)
    end do
    close(800)
  end do


end subroutine write_results
