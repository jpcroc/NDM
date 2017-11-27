
! **************************************************************
subroutine initeloss
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m,only:elstopforce,ngrdel,ev2erg,rang,tstep,ecelec,ibrake
  use var_pot,only:ntyp,cm,gamlt
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

  integer::i,j,j1,j2,j0,npr,k,nv1,iti
  real(double)::vel,vel2,sp,vnlt,v1,f1
  !  real(double),parameter::clum=29979245800
  real(double), pointer, dimension(:):: veloc, stoppow
  real(double),pointer::vmaxel(:)

  allocate (elstopforce(ntyp,2,0:ngrdel))
  allocate(vmaxel(ntyp))
  open (unit=99,file='elstop.in')
  if(rang==0)  write(6,*)'electronic loss eV ; eV/Ang tstep',tstep
  do i=1,ntyp
     if(rang==0)     write(6,*)'TYPE ',i
     read(99,*)npr
     allocate(veloc(0:npr))
     allocate(stoppow(0:npr))
     veloc(0)=0; stoppow(0)=0
     do j=1,npr
        read(99,*)vel,sp
        !        vel=clum*sqrt(1-1./((veloc(j)*ev2erg/(cm(i)*clum**2)+1)**2))
        vel2=dsqrt(2*vel*ev2erg/cm(i))
        !        write(6,'(2G15.5)')vel,vel2
        veloc(j)=vel2 ! vitesse en cm.sec-1
        stoppow(j)=sp*ev2erg*1e8 
#if(CHECK)
        write(62,'(I3,7G15.7)')j,vel,sp,veloc(j),stoppow(j), stoppow(j)/veloc(j),tstep*stoppow(j)/veloc(j)/cm(i)
#endif

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
#if(CHECK)
        if((rang==0).and.(mod(k,20)==0)) write(6,*)i,elstopforce(i,1,k),elstopforce(i,2,k)
                write(61,*)i,elstopforce(i,1,k),elstopforce(i,2,k)
#endif

     end do loopk
     deallocate(veloc)
     deallocate(stoppow)
     !     stop

     if (ibrake==2) then

        vnlt=sqrt(2*Ecelec*ev2erg/cm(i))
        v1=elstopforce(i,1,1)
        nv1=1+INT(vnlt/v1)
        if (nv1.gt.ngrdel) then
           write(6,*)'elstop velocity > 49, rebuild elstop.in,nv1',nv1
           stop
        end if
        f1=elstopforce(i,2,nv1)-(elstopforce(i,2,nv1)-elstopforce(i,2,nv1-1))*(nv1-vnlt/v1)
        gamlt(i)=f1/(cm(i)*vnlt)
        if(rang==0) write(6,'(A,I3,2G15.7)')'typ gaml',i,gamlt(i),vnlt
     end if
  end do

  return
end subroutine initeloss
