
! **************************************************************
subroutine initeloss
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m,only:elstopforce,ngrdel,ev2erg
  use var_pot,only:ntyp,cm
  use tab_imm_m,only:
!  use eam,only:
!  use eamerco,only:
!  use SMjuli,only:
!  use jqmod,only:
!  use neb_module,only:
!  use defcdp, ONLY :
!  use var_pot
#if(PARA)
  use mod_mpi
#endif 

  ! **************************************************************

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------

  integer::i,j,j1,j2,j0,npr,k
  real(double)::vel,vel2
!  real(double),parameter::clum=29979245800
  real(double), pointer, dimension(:):: veloc, stoppow
  real(double),pointer::vmaxel(:)

  allocate (elstopforce(ntyp,2,0:ngrdel))
  allocate(vmaxel(ntyp))
  open (unit=99,file='elstop.in')
  write(6,*)'electronic loss eV ; eV/Ang'
  do i=1,ntyp
     write(6,*)'TYPE ',i
     read(99,*)npr
     allocate(veloc(0:npr))
     allocate(stoppow(0:npr))
     veloc(0)=0; stoppow(0)=0
     do j=1,npr
        read(99,*)veloc(j),stoppow(j)
!        vel=clum*sqrt(1-1./((veloc(j)*ev2erg/(cm(i)*clum**2)+1)**2))
        vel2=dsqrt(2*veloc(j)*ev2erg/cm(i))
!        write(6,'(2G15.5)')vel,vel2
        veloc(j)=vel2 ! vitesse en cm.sec-1
        stoppow(j)=stoppow(j)*ev2erg*1e8 
!        write(62,*)j,veloc(j),stoppow(j)
     end do



     vmaxel(i)=veloc(npr)
     elstopforce(i,1,0)=0
     elstopforce(i,2,0)=0
     j0=0
     j1=0
     loopk: do k=1,ngrdel
        elstopforce(i,1,k)=vmaxel(i)*k/ngrdel
        loopj:   do j2=j0,npr
           if (veloc(j2).gt.elstopforce(i,1,k))then
              j0=j2-1
              j1=j2

              exit loopj
           end if
        end do loopj
        elstopforce(i,2,k)=stoppow(j0)+(stoppow(j1)-stoppow(j0))*(vmaxel(i)*k/ngrdel-veloc(j0))/(veloc(j1)-veloc(j0))
        if (mod(k,20)==0) write(6,*)i,elstopforce(i,1,k),elstopforce(i,2,k)
!        write(61,*)i,elstopforce(i,1,k),elstopforce(i,2,k)
     end do loopk
     deallocate(veloc)
     deallocate(stoppow)
!     stop
  end do

  return
end subroutine initeloss
