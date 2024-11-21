module eloss
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:ev2erg,rang,tstep,elosscel,l2T,erg2eV,iko,lspacendm
  USE var_pot, ONLY:ntyp,cm,gamlt
#ifdef PARA
  USE Tpara,only:COMM_space,myidsp,endmpi
#else
  use Tpara,only : nprocspace
  
#endif 
  use atomconfig,only:atom_config_d
   use newunit_mod,only:newunit

  use cellconfig,only:cell_config
  ! **************************************************************

  implicit none

  real(double):: elosselec,elosselec1 ! electronic losses for all atoms ; the PKA
  real(double):: elosselecstep,elosselecstep1 ! electronic losses for all atoms ; the PKA
  real(double),allocatable::elstopforce(:,:,:),gams(:)
  real(double):: tcelec,Ecelec ! coupure pour les pertes 駘ectroniques
  integer::ibrake   ! electronic slowing in cascades : 0 none, 1 down to ecelec, tcelec , 2 connected to Langevin
  integer::ngrdel,ikoloc



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


    integer::i,j,j1,j2,j0,npr,k,nv1
    real(double)::vel,vel2,sp,vnlt,v1,f1
    !  real(double),parameter::clum=29979245800
    real(double), allocatable, dimension(:):: veloc, stoppow
    real(double),allocatable::vmaxel(:)
    integer:: unit99
    allocate (elstopforce(ntyp,2,0:ngrdel))
    allocate(vmaxel(ntyp))
    allocate (gams(ntyp))
    call newunit(unit99)
    open (unit=unit99,file='elstop.in')
    if(rang==0)  write(6,*)'electronic loss eV ; eV/Ang tstep',tstep
    do i=1,ntyp
       select case(ibrake)
       case(3)
          if (rang==0) write (6,*)' ELSTOP READ : gams parameter in ps^-1'
          read(unit99,*) gams(i)
          gams(i)=gams(i)*1d12*cm(i)
       case default
          
          if(rang==0)     write(6,*)'TYPE ',i
       read(unit99,*)npr
       allocate(veloc(0:npr))
       allocate(stoppow(0:npr))
       veloc(0)=0; stoppow(0)=0
       do j=1,npr
          read(unit99,*)vel,sp
          !        vel=clum*sqrt(1-1./((veloc(j)*ev2erg/(cm(i)*clum**2)+1)**2))
          vel2=dsqrt(2*vel*ev2erg/cm(i))
          !        write(6,'(2G15.5)')vel,vel2
          veloc(j)=vel2 ! vitesse en cm.sec-1
          stoppow(j)=sp*ev2erg*1e8 
#ifdef CHECK
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
#ifdef CHECK
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
             write(6,*)'elstop velocity > ngrdel, rebuild elstop.in,nv1',nv1
#ifdef PARA
             call endMPI
#endif 
             stop
          end if
          f1=elstopforce(i,2,nv1)-(elstopforce(i,2,nv1)-elstopforce(i,2,nv1-1))*(nv1-vnlt/v1)
          gamlt(i)=f1/(cm(i)*vnlt)
!          if(rang==0) write(6,'(A,I3,2G15.7)')'typ gaml',i,gamlt(i),vnlt
       end if
    end select
 end do
 close (unit99)

    return
  end subroutine initeloss



  subroutine calceloss(celndm,atdml)
#ifdef PARA
    use mod_para,only : nprocspace
#endif

    class (cell_config),intent(in)::celndm
    class(atom_config_d)::atdml

    
    real(double)::ekin,vn,v1,f1,eta,etavc,f1vc,vc
    integer::koo,i,nv1,ic,iti
    integer, save:: icall=0

    icall=icall+1
!    if (icall==1) then
        elosselecstep=0
    	elosselecstep1=0
!     endif
#if PARA
  ikoloc=0
  do i=1,atdml%im
     if (atdml%num_at_glob(i)==iko) ikoloc=i
  end do
  
#else
ikoloc=iko
#endif

if (L2T.eqv..true.)     elosscel(:)=0
!write(6,*)gams
    do i=1,atdml%im
       if(tcelec.gt.0) then
          koo = atdml%ielat(i)                          ! Numero de la cellule
          if (celndm%tempc(koo).le.tcelec) cycle
       end if

       vn= atdml%vp(1,i)**2+atdml%vp(2,i)**2+atdml%vp(3,i)**2
       ekin=0.5*erg2ev*vn*cm(atdml%ityp(i))
!       if (ikoloc.gt.0)  then
!          if (i==ikoloc) write(6,*)ekin,ecelec
!       end if
       if (ekin.gt.Ecelec) then

          select case(ibrake)


       case(3)
          do ic=1,3
             atdml%fp(ic,i)=atdml%fp(ic,i)-atdml%vp(ic,i)*gams(atdml%ityp(i))
             Elosselecstep=Elosselecstep+(atdml%vp(ic,i)*gams(atdml%ityp(i)))*(atdml%vp(ic,i)*tstep)*erg2ev
             
             if (L2T.eqv..true.) then
                elosscel(atdml%ielat(i))=elosscel(atdml%ielat(i))+(atdml%vp(ic,i)*gams(atdml%ityp(i)))*(atdml%vp(ic,i)*tstep)
             end if
             if (i==ikoloc)then 

                !                write(6,*)'elfp',fp(ic,i)
                Elosselecstep1=Elosselecstep1+(atdml%vp(ic,i)*gams(atdml%ityp(i)))*(atdml%vp(ic,i)*tstep)*erg2ev
             end if
          end do

          case default
             vn=sqrt(vn)

          !	write(6,*)'RG',rang,i,ekin
          v1=elstopforce(atdml%ityp(i),1,1)
          !           write(6,*)v1,vn
          nv1=1+INT(vn/v1)
          if (nv1.gt.ngrdel) then
             write(6,*)'elstop velocity > 49, rebuild elstop.in'
#ifdef PARA
             call endMPI
#endif 
            stop
          end if
          f1=elstopforce(atdml%ityp(i),2,nv1)-(elstopforce(atdml%ityp(i),2,nv1)-&
               &elstopforce(atdml%ityp(i),2,nv1-1))*(nv1-vn/v1)
          !           write (6,'(A,4G15.7)')'felstop ',f1,vn, vn/v1,elstopforce(ityp(i),2,nv1)
          if (f1.le.0) then
             write(6,*)'f1<0 ?', f1
             stop
          end if
          if (ibrake==2) then
             iti=atdml%ityp(i)
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
             atdml%fp(ic,i)=atdml%fp(ic,i)-atdml%vp(ic,i)*f1/vn
             Elosselecstep=Elosselecstep+(atdml%vp(ic,i)*f1/vn)*(atdml%vp(ic,i)*tstep)*erg2ev
             if (L2T.eqv..true.) then
                elosscel(atdml%ielat(i))=elosscel(atdml%ielat(i))+(atdml%vp(ic,i)*f1/vn)*(atdml%vp(ic,i)*tstep)
             end if
             if (atdml%num_at_glob(i)==iko)then 

                !                write(6,*)'elfp',fp(ic,i)
                Elosselecstep1=Elosselecstep1+(atdml%vp(ic,i)*f1/vn)*(atdml%vp(ic,i)*tstep)*erg2ev
             end if
          end do
          !                 write (6,*)'felstop',f1,vn                

       end select
       end if

    end do
    !	write(6,*)'RG el',rang,elosselecstep,elosselecstep1
#ifdef PARA
if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
   call comm_space%sum(elosselecstep)
   call comm_space%sum(elosselecstep1)
   elosselec=elosselec+elosselecstep
   elosselec1=elosselec1+elosselecstep1
       !if l2T
       if (allocated(elosscel)) then
          call comm_space%sum(elosscel)
       end if
    end if
#else

#endif
!if (rang==0)   write(6,*)'TEST electronic losses ', elosselec, elosselec1


  end subroutine calceloss

end module eloss
