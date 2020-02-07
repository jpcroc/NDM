module eloss
  USE T_kind_param_m, ONLY:  double
  use gen_com_m,only:ev2erg,rang,tstep,elosscel,tempc,l2T,erg2eV,im,iko,noxyz
  use var_pot,only:ntyp,cm,gamlt
  use tab_imm_m,only:fp,vp,ityp,ielat,num_at_glob


  !  use eam,only:
  !  use eamerco,only:
  !  use SMjuli,only:
  !  use jqmod,only:
  !  use neb_module,only:
  !  use defcdp, ONLY :
  !  use var_pot
#ifdef PARA
  use mod_para
#endif 

  ! **************************************************************

  implicit none

  real(double):: elosselec,elosselec1 ! electronic losses for all atoms ; the PKA
  real(double):: elosselectot,elosselectot1 ! electronic losses for all atoms ; the PKA
  real(double),pointer::elstopforce(:,:,:)
  real(double):: tcelec,Ecelec ! coupure pour les pertes 駘ectroniques
  integer::ibrake   ! electronic slowing in cascades : 0 none, 1 down to ecelec, tcelec , 2 connected to Langevin
  integer::ngrdel



contains

  ! **************************************************************
  subroutine initeloss
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
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
#ifdef PARA
	            call MPI_FINALIZE(ierr)
#endif 


             stop
          end if
          f1=elstopforce(i,2,nv1)-(elstopforce(i,2,nv1)-elstopforce(i,2,nv1-1))*(nv1-vnlt/v1)
          gamlt(i)=f1/(cm(i)*vnlt)
          if(rang==0) write(6,'(A,I3,2G15.7)')'typ gaml',i,gamlt(i),vnlt
       end if
    end do

    return
  end subroutine initeloss



  subroutine calceloss

    real(double)::ekin,vn,v1,f1,eta,etavc,f1vc,vc
    integer::koo,i,nv1,ic,iti
    integer, save:: icall=0

#ifdef PARA
    real(double), allocatable,dimension(:)::elosscel_tot
    if (allocated(elosscel)) allocate (elosscel_tot(noxyz))
#endif
    icall=icall+1
    if (icall==1) then
        elosselec=0
    	elosselec1=0
     endif

    if (L2T.eqv..true.)     elosscel(:)=0
    do i=1,im
       if(tcelec.gt.0) then
          koo = ielat(i)                          ! Numero de la cellule
          if (tempc(koo).le.tcelec) cycle
       end if

       vn= vp(1,i)**2+vp(2,i)**2+vp(3,i)**2
       ekin=0.5*erg2ev*vn*cm(ityp(i))
       if (ekin.gt.Ecelec) then


          !	write(6,*)'RG',rang,i,ekin
          vn=sqrt(vn)
          v1=elstopforce(ityp(i),1,1)
          !           write(6,*)v1,vn
          nv1=1+INT(vn/v1)
          if (nv1.gt.ngrdel) then
             write(6,*)'elstop velocity > 49, rebuild elstop.in'
#ifdef PARA
	            call MPI_FINALIZE(ierr)
#endif 
            stop
          end if
          f1=elstopforce(ityp(i),2,nv1)-(elstopforce(ityp(i),2,nv1)-elstopforce(ityp(i),2,nv1-1))*(nv1-vn/v1)
          !           write (6,'(A,4G15.7)')'felstop ',f1,vn, vn/v1,elstopforce(ityp(i),2,nv1)
          if (f1.le.0) then
             write(6,*)'f1<0 ?', f1
             stop
          end if
          if (ibrake==2) then
             iti=ityp(i)
!	     write(6,*)rang,i,iti
             eta=f1/vn
             vc=sqrt(2*Ecelec*ev2erg/cm(iti))
             v1=elstopforce(iti,1,1)
             nv1=1+INT(vc/v1)
             f1vc=elstopforce(iti,2,nv1)-(elstopforce(iti,2,nv1)-elstopforce(iti,2,nv1-1))*(nv1-vc/v1)
             etavc=f1vc/vc
             
             f1=f1-etavc*vn
          !           write(6,'(A,5G15.7)')'EL222',ekin,vn,f1,gamlt(ityp(i))*sqrt(cm(ityp(i))*2*Ecelec*ev2erg),f1/vn
          !           write(6,*)
          end if
          do ic=1,3
             fp(ic,i)=fp(ic,i)-vp(ic,i)*f1/vn
             !              Elosselec=Elosselec+(vp(ic,i)*f1/vn)*(vp(ic,i)*tstep)
             Elosselec=Elosselec+(vp(ic,i)*f1/vn)*(vp(ic,i)*tstep)*erg2ev
             if (L2T.eqv..true.) then
                elosscel(ielat(i))=elosscel(ielat(i))+(vp(ic,i)*f1/vn)*(vp(ic,i)*tstep)
             end if
             if (num_at_glob(i)==iko)then 

                !                write(6,*)'elfp',fp(ic,i)
                Elosselec1=Elosselec1+(vp(ic,i)*f1/vn)*(vp(ic,i)*tstep)*erg2ev
             end if
          end do
          !                 write (6,*)'felstop',f1,vn                



       end if

    end do
    !	write(6,*)'RG el',rang,elosselec,elosselec1
#ifdef PARA
    elosselectot=0
    elosselectot1=0
    call MPI_ALLREDUCE(elosselec,elosselectot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(elosselec1,elosselectot1,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    !if l2T
    if (allocated(elosscel)) then
       call MPI_ALLREDUCE(elosscel,elosscel_tot,noxyz,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
       elosscel=elosscel_tot
       deallocate (elosscel_tot)
    end if

#else
    elosselectot=elosselec
    elosselectot1=elosselec1

#endif
!  write(6,*)'TEST electronic losses ', elosselectot, elosselectot1


  end subroutine calceloss

end module eloss
