module elec_cell
  use T_kind_param_m
  use gen_com_m, only : nox,noy,noz, noxyz,nzl,bk,imm,nato,last,im_glob,tstep,erg2eV,pi,rang,elosscel,lenfnam,fnam,lrestart&
       &,joule2erg,erg2eV,it,timel,it,igen,lrestart
  use var_pot,only:cm
  use tab_imm_m, only : num_at_glob,ielat
  use eloss,only :Ecelec ,elstopforce,ngrdel
  !
  implicit none
  type :: ecelltype
     real(double)::temp
     real(double)::tempIon
     integer::ixb(3)
     logical::lionovlp
     real(double)::Qi2e
     integer::nIon
     integer::nIonS
     real(double)::En
  end type ecelltype

  !type :: voisceltype
  !  integer::ncelvois
  !  integer,allocatable::icelvois(:)
  !end type voisceltype


  !OR
  type(ecelltype), dimension(:,:,:),allocatable::ecell
  !type(voisceltye), 


  real(double)::cellside(3)  
  !  real(double)::e_atcell(3,3)

  real(double):: kappaE, Ce ! electronic thermal diff, thermal cond, Ce
  integer:: necycle,necyclemin ! number of electronic steps
  real(double):: etstep

  integer:: nex, ney, nez ,nexyz! number of cells in x,y,z directions
  integer:: nexov,neyov,nezov ! number of cells within the MD BOX in x,y,z directions
  integer:: nexmp,neymp,nezmp ! number of cells in x,y,z directions outside the MD box on each side
  real(double)::wex,wey,wez ! width of electronic system (wex=cellside(1)*nxe)
  real(double)::deltaxyz(3),Ax(3)
  real(double)::CeC,KeC,T0,k0T,GepC
  real(double)::Vecell
  integer::i2t ! type of 2T model   i2T=1 ! model Croc&Murphy (0= DD)model Croc&Murphy (0= DD)
  real(double):: t_cpl ! ep coupling after t_cpl iterations
  integer:: ibc ! borders conditions
  real(double)::Eelec,Teavg ! energy stored in cells
  real(double)::Tecmax
  integer::ietm(3)
  integer, parameter::nTmax=40000
  integer:: ncer
  integer:: itetec,integrTtype,igenelec
  real(double),dimension(:), allocatable:: CedT,Eedt,KedT,GepdT
  logical lalletemp 
contains

  subroutine readelec
    integer:: luelec=654
    integer::ic
    integer::ix,iy,iz
    character :: fnamedin*80
    namelist /inputelec/nexov,neyov,nezov,deltaxyz,Cec,KeC,ibc,T0,k0T,GepC,necyclemin,i2T&
         &,t_cpl,ncer,iteTec,lalletemp,integrTtype,igenelec
    i2T=1 ! model Croc&Murphy (0= DD)
    T0=300. !T0=borders temperature
    k0T=0.8
    necyclemin=10
    ibc=-1
    KeC=0
    CeC=0
    GepC=0
    t_cpl=-1
    itetec=-1
    lalletemp=.false. 
    integrTtype=1 ! type of temperature integrator 1= std , 2= for insulator
    deltaxyz(:)=0
    ncer=1000
    igenelec=igen

    if (rang.eq.0) write(6,*) 
    if (rang.eq.0) write(6,*) 
    if (rang.eq.0) write(6,*) '>>>>>>>>>>> entree readelec  input units are SI, internal units are cgs'

    fnamedin=  fnam(1:lenfnam)//'.edin'

!    open(unit=luelec, file='elec.in', status='unknown')
    open(unit=luelec, file=fnamedin, status='unknown')
    read (luelec, nml=inputelec)
    close(luelec)
    if (lrestart) igenelec=1
    necycle=necyclemin
    GepC=GepC*joule2erg*1d-6
    CeC=Cec*joule2erg*1d-6
    KeC=KeC*joule2erg*1d-2
    deltaxyz=deltaxyz*1d-8
    if (rang.eq.0)    write(6,*)nexov,neyov
    if(nexov==0)nex=nox
    if(neyov==0)ney=noy
    if(nezov==0)nez=noz
    if (ibc==-1) then
       if (rang.eq.0)       write(6,*)'ibc=-1, stop'
       stop
    end if
    !    call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
    !    do ic=1,3
    !       normat(ic)=sqrt(sum(bg(:,ic)**2))
    !       nzl(ic)=1.0/normat(ic)
    !    enddo

    cellside(1)=nzl(1)/nexov
    nexmp=1+int(deltaxyz(1)/cellside(1))
    nex=2*nexmp+nexov
    deltaxyz(1)=nexmp*cellside(1)

    cellside(2)=nzl(2)/neyov
    neymp=1+int(deltaxyz(2)/cellside(2))
    ney=2*neymp+neyov
    deltaxyz(2)=neymp*cellside(2)

    cellside(3)=nzl(3)/nezov
    nezmp=1+int(deltaxyz(3)/cellside(3))
    nez=2*nezmp+nezov
    deltaxyz(3)=nezmp*cellside(3)

    Ax(1)=cellside(2)*cellside(3)
    Ax(2)=cellside(3)*cellside(1)
    Ax(3)=cellside(1)*cellside(2)
    Vecell=cellside(1)*cellside(2)*cellside(3)
    !volume
    if (rang==0) then
       write(6,*)'nex ney nez',nex,ney,nez
       write(6,*)'nexmp neymp nezmp',nexmp,neymp,nezmp
       write(6,*)'deltayz= ',deltaxyz
       write(6,*)'cellside',cellside
    end if
    allocate(EedT(0:nTmax))
    if (CeC.lt.0) call prepCe
    if (GepC.lt.0) call prepGep
    if (KeC.lt.0) call prepKe

    allocate(ecell(nex,ney,nez))

    if (rang==0) then
       if (igenelec==1) then
          write(6,*)'ELECTRONIC TEMPERATURE READ FROM FILE, T0 NOT USED'
          if (igen.ne.1) write(6,*)'BUT NOT THE ATOMIC CONFIGURATION'
       else
          write(6,*)'ELECTRONIC TEMPERATURE INITIALZED AT T0=',T0
          if (igen==1)write(6,*)'BUT ATOMIC CONFIGURATION IS READ'
       end if
    end if
    if (igenelec==1) then
       call restartelec
    else
       ecell(:,:,:)%lionovlp=.false.
       ecell(:,:,:)%temp=T0
       !    ecell(5,5,5)%temp=1000
       ecell(:,:,:)%En=0
       do ix=1,nex
          do iy=1,ney
             do iz=1,nez
                ecell(ix,iy,iz)%ixb(1)=ix*cellside(1)-deltaxyz(1)
                ecell(ix,iy,iz)%ixb(2)=iy*cellside(2)-deltaxyz(2)
                ecell(ix,iy,iz)%ixb(3)=iz*cellside(3)-deltaxyz(3)
                !            if((ecell(ix,iy,iz)%ixb(1).ge.0).and.(ecell(ix,iy,iz)%ixb(1).le.zl(1))&
                !                 &.and.(ecell(ix,iy,iz)%ixb(2).ge.0).and.(ecell(ix,iy,iz)%ixb(2).le.zl(2))&
                !                 &.and.(ecell(ix,iy,iz)%ixb(3).ge.0).and.(ecell(ix,iy,iz)%ixb(3).le.zl(3)))&
                !                 & ecell%lionrecov=.true.
                if((ix.gt.nexmp).and.(ix.le.nexmp+nexov).and.&
                     &(iy.gt.neymp).and.(iy.le.neymp+neyov).and.&
                     (iz.gt.nezmp).and.(iz.le.nezmp+nezov))&
                     & ecell(ix,iy,iz)%lionovlp=.true.
             end do
          end do
       end do
    end if


    if (rang.eq.0) write(6,*) 
    if (rang.eq.0) write(6,*) 

    etstep=tstep/necyclemin
    !     if (etstep.gt.6d-17)then
    !        etstep=2d-16
    !        necycle=int(tstep/etstep)
    !        write(6,*)'chgt etstep',etstep,necycle
    !     end if


    return

  end subroutine readelec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine TTlangevin(xp, vp, fp,ityp,il,Gl)


#if(PARA)
    use mod_mpi
#endif

    real(double)  :: xp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: fp(3,imm)
    real(double)  :: Gl(3,imm)
    integer  :: ityp(imm)
    integer::il
    real(double)::rga
    integer :: i,ic,ko,i2,nv1
    real(double) :: u1,u2,gamlat,ekin,vpn2,v1,f1,vn,Gep
    !langevin codé à partir du poly de Gabriel Stolz page 84, dans une version avec expoentielle comme Manuel et Cosmin
    ! adapted to 2T model
    integer :: ixyze(3)

    if (i2t==0)elosscel(:)=0
    do ko = 1, noxyz
#if(PARA)
       if (proc_cell(ko).ne.myid) cycle
#endif
       if (nato(ko)==0) cycle
       call nox_2_nex(ko,ixyze)
       call GepT(Gep,ecell(ixyze(1),ixyze(2),ixyze(3))%temp)
       do i2 = 1, nato(ko)
          i = last(i2,ko)
          if (num_at_glob(i).gt.im_glob) cycle
          select case (i2t)
          case(1)
             !check for velcocity
             vpn2 = vp(1,i)**2+vp(2,i)**2+vp(3,i)**2
             ekin=0.5*erg2ev*vpn2*cm(ityp(i))
             if (ekin.gt.Ecelec) then
                gamlat=0
             else
                gamlat=VeCell*Gep/(3*bk*ecell(ixyze(1),ixyze(2),ixyze(3))%NionS)
                !if(i==1) write(6,*)'gamf',VeCell,Gep,bk,ecell(ixyze(1),ixyze(2),ixyze(3))%Nion
             end if
          case(0)
             vpn2 = vp(1,i)**2+vp(2,i)**2+vp(3,i)**2
             ekin=0.5*erg2ev*vpn2*cm(ityp(i))
             if (ekin.gt.Ecelec) then
                vn=sqrt(vpn2)
                v1=elstopforce(ityp(i),1,1)
                !           write(6,*)v1,vn
                nv1=1+INT(vn/v1)
                if (nv1.gt.ngrdel) then
                   if (rang.eq.0)  write(6,*)'elstop velocity > 49, rebuild elstop.in'
#if(PARA)
	            call MPI_FINALIZE(ierr)
#endif 

                   stop
                end if
                f1=elstopforce(ityp(i),2,nv1)-(elstopforce(ityp(i),2,nv1)-elstopforce(ityp(i),2,nv1-1))*(nv1-vn/v1)
                do ic=1,3
                   elosscel(ielat(i))=elosscel(ielat(i))+(vp(ic,i)*f1/vn)*(vp(ic,i)*tstep)
                end do
                gamlat=f1/(cm(ityp(i))*vn)
                if (timel.gt.t_cpl) then
                   gamlat=gamlat+VeCell*Gep/(3*bk*ecell(ixyze(1),ixyze(2),ixyze(3))%Nion)
                end if
             else
                if (timel.gt.t_cpl) then
                   gamlat=VeCell*Gep/(3*bk*ecell(ixyze(1),ixyze(2),ixyze(3))%Nion)
                else 
                   gamlat=0.
                end if
                !if(i==1) write(6,*)'gamf',VeCell,Gep,bk,ecell(ixyze(1),ixyze(2),ixyze(3))%Nion
             end if

          end select

          !"standard" Langevin algorythm         
          select case (il)
          case(1)
             rga=exp(-gamlat*tstep/2)

             !             if ((it.ge.1000).and.(i.lt.20)) write(6,'(A,3G15.7)')'gamstd 1 ',gamlat,rga
             !                if ((it.ge.1000).and.(i.lt.20))then
             !                write(6,'(A,2I5,3G15.7)')' vpa1',it,i,vp(:,i)
             !                write(6,'(A,2I5,3G15.7)')' xpa1',it,i,xp(:,i)
             !                write(6,'(A,2I5,3G15.7)')' fpa1',it,i,fp(:,i)
             !                write(6,'(A,3I2,G15.7)')'ixyze1',ixyze(1),ixyze(2),ixyze(3),ecell(ixyze(1),ixyze(2),ixyze(3))%Temp
             !             end if

             do ic=1,3
                !  write(6,*)'ct',cm(ityp(1)),tstep
                call random_number(u1)
                call random_number(u2)
                Gl(ic,i)=sqrt(-2.*log(u1))*cos(2.*pi*u2)   
                !if (i==1)write(6,*) vp(ic,i)

                vp(ic,i) = vp(ic,i)*rga+ fp(ic,i)*tstep/(cm(ityp(i))*2)+Gl(ic,i)*&
                     &sqrt(cm(ityp(i))*bk*ecell(ixyze(1),ixyze(2),ixyze(3))%temp*(1-rga))/cm(ityp(i))
                !if (i==1)write(6,*) vp(ic,i)

             end do
             !if (i==1)write(6,*) 
          case(2)
             rga=exp(-gamlat*tstep/2)
             !            if ((it.ge.1000).and.(i.lt.20)) write(6,'(A,3G15.7)')'gamstd 2 ',gamlat,rga
             !            if ((it.ge.1000).and.(i.lt.20))then
             !               write(6,'(A,2I5,3G15.7)')' vpa2',it,i,vp(:,i)
             !               write(6,'(A,2I5,3G15.7)')' xpa2',it,i,xp(:,i)
             !               write(6,'(A,2I5,3G15.7)')' fpa2',it,i,fp(:,i)
             !               write(6,'(A,3I2,G15.7)')'ixyze2',ixyze(1),ixyze(2),ixyze(3),ecell(ixyze(1),ixyze(2),ixyze(3))%Temp
             !            end if


             do ic=1,3
                !if (i==1)write(6,*) vp(ic,i)


                vp(ic,i) = vp(ic,i)*rga+ fp(ic,i)*tstep/(cm(ityp(i))*2)+Gl(ic,i)*&
                     &sqrt(cm(ityp(i))*bk*ecell(ixyze(1),ixyze(2),ixyze(3))%temp*(1-rga))/cm(ityp(i))
                !if (i==1)write(6,*) vp(ic,i)


             end do
             !            if ((it.ge.1000).and.(i.lt.20))write(6,'(A,2I5,3G15.7)')&
             !&' vpa2',it,i,vp(:,i)

          case default 
             if (rang.eq.0)             write(6,*)'check ilangevin'
             stop
          end select
       end do
    end do
  end subroutine TTlangevin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dynelec
    integer::ite,iex,iey,iez,iet
    real(double)nexttemp (nex,ney,nez)

    call calc_Qi2e
    do ite=1,necycle
       call Tevolv(ite)
    end do
    if ((itetec.gt.0).and.(mod(it,itetec)==0))then
       do iex=1,nex
          do iey=1,ney
             do iez=1,nez
                iet=iex+iey*nex+iez*nex*ney
!                write(6,*)iet,iex,iey,iez
                if (lalletemp.or.(ecell(iex,iey,iez)%lionovlp))then
                   !                   write(6,*)
                   !                   write(6,'(A,3I3,G15.7)')'ecell%temp',iex,iey,iez,ecell(iex,iey,iez)%lionovlp
                   !                   write(6,'(A,I5,2G15.7)')'ecell%temp',it,timel,ecell(iex,iey,iez)%temp
!                   if (rang.eq.0) write(iet,*)it,timel, ecell(iex,iey,iez)%temp
                end if
             end do
          end do
       end do
    end if

  end subroutine dynelec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_Qi2e

    integer::ko,ixe,iye,ize
    integer :: ixyze(3)
    real(double)::Gep
    ecell(:,:,:)%Qi2e=0
    !ES stopping
    do ko=1,noxyz
       call nox_2_nex(ko,ixyze)
       !       if (elosscel(ko).ne.0) then
       !          write(6,'(4I3)')ko,ixyze(:)
       !          write(6,*)'ko',elosscel(ko)
       !       end if
       ecell(ixyze(1),ixyze(2),ixyze(3))%Qi2e=ecell(ixyze(1),ixyze(2),ixyze(3))%Qi2e+elosscel(ko)/tstep
    end do
    !    write(6,*)ecell%Qi2e

    !EP coupling    
    do ixe=1,nex
       do iye=1,ney
          do ize=1,nez
             if(ecell(ixe,iye,ize)%lionovlp.eqv..true.)then
                call GepT(Gep,ecell(ixyze(1),ixyze(2),ixyze(3))%temp)
                if (timel.gt.t_cpl) then
                   ecell(ixe,iye,ize)%Qi2e=ecell(ixe,iye,ize)%Qi2e-&
                        &Gep*( ecell(ixe,iye,ize)%temp- ecell(ixe,iye,ize)%tempIon)*Vecell
                endif
             end if
             !             if (ecell(ixe,iye,ize)%Qi2e.ne.0) &
             !&  write(6,'(I6,A,3I4,G15.7)')it,' TRF ',ixe,iye,ize,ecell(ixe,iye,ize)%Qi2e

          end do
       end do
    end do

  end subroutine calc_Qi2e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine Tevolv(ite)

    real(double):: nexttemp(nex,ney,nez)
    real(double)::Ce, Ke,tfact
    integer:: iex,iey,iez,iet,ite,itn,itc
    real(double)Txm1,Txp1,Tym1,Typ1,Tzm1,Tzp1,TC,Tm,deltaE,Ebase,EcN



    call borders(ibc)
    nexttemp=ecell(:,:,:)%temp
    do iex=2,nex-1
       do iey=2,ney-1
          do iez=2,nez-1
             TC=ecell(iex,iey,iez)%temp
             Txm1=ecell(iex-1,iey,iez)%temp
             Txp1=ecell(iex+1,iey,iez)%temp
             Tym1=ecell(iex,iey-1,iez)%temp
             Typ1=ecell(iex,iey+1,iez)%temp
             Tzm1=ecell(iex,iey,iez-1)%temp
             Tzp1=ecell(iex,iey,iez+1)%temp 
             !             if (it.ge.900)             write(6,'(I5,A,3I3)')it,'it iex',iex,iey,iez
             !             if (it.ge.900)             write(6,'(I8,A,3I4,7G12.6)')it*necycle+ite,' TCpp ',iex,iey,iez,&
             !&TC,txm1,txp1,tym1,typ1,tzm1,tzp1
!             write(6,'(I5,A,3I3)')it,'it iex',iex,iey,iez
!             write(6,'(I8,A,3I4,7G15.6)')it*necycle+ite,' TCpp ',iex,iey,iez,&
!             &TC,txm1,txp1,tym1,typ1,tzm1,tzp1

             select case(integrTtype)
             case(1)
                call CeT(Ce, TC)


                !            write(6,'(2G15.7)')'ke ce alpha',iex,iey,iez,ke,ce,ke/ce
                !alpha=Ke/Ce
                !tfact=min(Ax(1)/cellside(1),Ax(2)/cellside(2),Ax(3)/cellside(3))/(6*alpha)
                !            write(6,*)'etstep,tfact',etstep,tfact


                !iex
                Tm=0.5*(Txm1+TC)
                call KeT(Ke, Tm)
                nexttemp(iex,iey,iez)=nexttemp(iex,iey,iez)+etstep*(Ax(1)/(cellside(1)*Vecell))*(Ke/Ce)&
                     &*(Txm1-TC)

                Tm=0.5*(Txp1+TC)
                call KeT(Ke, Tm)
                nexttemp(iex,iey,iez)=nexttemp(iex,iey,iez)+etstep*(Ax(1)/(cellside(1)*Vecell))*(Ke/Ce)&
                     &*(Txp1-TC)

                !iey
                Tm=0.5*(Tym1+TC)
                call KeT(Ke, Tm)
                nexttemp(iex,iey,iez)=nexttemp(iex,iey,iez)+etstep*(Ax(2)/(cellside(2)*Vecell))*(Ke/Ce)&
                     &*(Tym1-TC)

                Tm=0.5*(Typ1+TC)
                call KeT(Ke, Tm)
                nexttemp(iex,iey,iez)=nexttemp(iex,iey,iez)+etstep*(Ax(2)/(cellside(2)*Vecell))*(Ke/Ce)&
                     &*(Typ1-TC)


                !iez
                Tm=0.5*(Tzm1+TC)
                call KeT(Ke, Tm)
                nexttemp(iex,iey,iez)=nexttemp(iex,iey,iez)+etstep*(Ax(3)/(cellside(3)*Vecell))*(Ke/Ce)&
                     &*(Tzm1-TC)

                Tm=0.5*(Tzp1+TC)
                call KeT(Ke, Tm)
                nexttemp(iex,iey,iez)=nexttemp(iex,iey,iez)+etstep*(Ax(3)/(cellside(3)*Vecell))*(Ke/Ce)&
                     &*(Tzp1-TC)

                nexttemp(iex,iey,iez)=nexttemp(iex,iey,iez)+(etstep/Ce)*ecell(iex,iey,iez)%Qi2e/Vecell
                !                if (ecell(iex,iey,iez)%Qi2e.ne.0)write(6,'(I6,A,3I3,G15.6)')it,' qi2ets ',iex,iey,iez,&
                !&etstep*ecell(iex,iey,iez)%Qi2e
                ! integration for insulators (energy input then temperature change)
!                write(6,*)nexttemp(iex,iey,iez)
             case(2)
                deltaE=0
                Tm=0.5*(Txm1+TC)
                call KeT(Ke, Tm)
                deltaE=deltaE+etstep*(Ax(1)/(cellside(1)))*Ke*(Txm1-TC)

                Tm=0.5*(Txp1+TC)
                call KeT(Ke, Tm)
                deltaE=deltaE+etstep*(Ax(1)/(cellside(1)))*Ke*(Txp1-TC)

                !iey
                Tm=0.5*(Tym1+TC)
                call KeT(Ke, Tm)
                deltaE=deltaE+etstep*(Ax(2)/(cellside(2)))*Ke*(Tym1-TC)

                Tm=0.5*(Typ1+TC)
                call KeT(Ke, Tm)
                deltaE=deltaE+etstep*(Ax(2)/(cellside(2)))*Ke*(Typ1-TC)

                !iez
                Tm=0.5*(Tzm1+TC)
                call KeT(Ke, Tm)
                deltaE=deltaE+etstep*(Ax(3)/(cellside(3)))*Ke*(Tzm1-TC)

                Tm=0.5*(Tzp1+TC)
                call KeT(Ke, Tm)
                deltaE=deltaE+etstep*(Ax(3)/(cellside(3)))*Ke*(Tzp1-TC)

                !source term

                deltaE=deltaE+etstep*ecell(iex,iey,iez)%Qi2e

                ecell(iex,iey,iez)%En=deltaE+ecell(iex,iey,iez)%En
                EcN=ecell(iex,iey,iez)%En/Vecell
                !                write(6,*)iex,iey,iez,ecN,deltaE

                if (deltaE.gt.0) then
                   iTC=Int(TC)
                   itn=itC
                   !2                  continue
                   do while(EedT(itn).lt.EcN)
                      !                      write(6,*)itn,eedt(itn)
                      itn=itn+1
                   end do
                   !                   write(6,*)itn,eedt(itn)
                   Ebase=EedT(itn-1)
                   nexttemp(iex,iey,iez)=(itn-1)+(Ecn-Ebase)/(Eedt(itn)-EedT(itn-1))

                else if(deltaE.lt.0) then

                   iTC=1+Int(TC)
                   itn=itC

                   !2                  continue
                   do while(EedT(itn).gt.ecN)
                      !                   write(6,*)itn,eedt(itn)
                      itn=itn-1
                   end do
                   !                   write(6,*)itn,eedt(itn)
                   Ebase=EedT(itn+1)
                   nexttemp(iex,iey,iez)=(itn+1)+(Ecn-Ebase)/(Eedt(itn+1)-EedT(itn))

                else  
                   ! no change in energy or temperature
                   nexttemp(iex,iey,iez)=ecell(iex,iey,iez)%temp
                end if

             end select
          end do
       end do
    end do
    ecell(:,:,:)%temp=nexttemp

    !    write(6,*)'Elec222',ecell(2,2,2)%temp*Ce*Vecell
    Teavg=SUM(ecell(2:nex-1,2:ney-1,2:nez-1)%temp)/((nex-2)*(ney-2)*(nez-2))
    Tecmax=maxval(ecell(:,:,:)%temp)
    ietm(:)=maxloc(ecell(:,:,:)%temp)
    if (CeC.lt.0) then
       Eelec=0
       do iex=2,nex-1
          do iey=2,ney-1
             do iez=2,nez-1
                Eelec=Eelec+EedT(int(ecell(iex,iey,iez)%temp))*Vecell
             end do
          end do
       end do
    else
       Eelec=SUM(ecell(2:nex-1,2:ney-1,2:nez-1)%temp*Ce)*Vecell
    end if


    !    write(6,*)'Eelectr Teavg cell ',Eelec*erg2eV,Teavg
    return
  end subroutine Tevolv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine borders (ib)
    integer,intent(in):: ib
    integer:: iex,iey,iez
    select case(ib)
    case(1) ! T cst at borders
       do iey=1,ney
          do iez=1,nez
             ecell(1,iey,iez)%temp=T0
             ecell(nex,iey,iez)%temp=T0
          end do
       end do
       do iex=1,nex
          do iez=1,nez
             ecell(iex,1,iez)%temp=T0
             ecell(iex,ney,iez)%temp=T0
          end do
       end do
       do iex=1,nex
          do iey=1,ney
             ecell(iex,iey,1)%temp=T0
             ecell(iex,iey,nez)%temp=T0
          end do
       end do
    case(2) ! zero flux
       do iey=1,ney
          do iez=1,nez
             ecell(1,iey,iez)%temp=ecell(2,iey,iez)%temp
             ecell(nex,iey,iez)%temp=ecell(nex-1,iey,iez)%temp
          end do
       end do
       do iex=1,nex
          do iez=1,nez
             ecell(iex,1,iez)%temp=ecell(iex,2,iez)%temp
             ecell(iex,ney,iez)%temp=ecell(iex,ney-1,iez)%temp
          end do
       end do
       do iex=1,nex
          do iey=1,ney
             ecell(iex,iey,1)%temp=ecell(iex,iey,2)%temp
             ecell(iex,iey,nez)%temp=ecell(iex,iey,nez-1)%temp
          end do
       end do
    case(3) !partial T transfer
       do iey=1,ney
          do iez=1,nez
             ecell(1,iey,iez)%temp= T0+K0T*(ecell(2,iey,iez)%temp-T0)
             ecell(nex,iey,iez)%temp=T0+K0T*(ecell(nex-1,iey,iez)%temp-T0)
          end do
       end do
       do iex=1,nex
          do iez=1,nez
             ecell(iex,1,iez)%temp=T0+K0T*(ecell(iex,2,iez)%temp-T0)
             ecell(iex,ney,iez)%temp=T0+K0T*(ecell(iex,ney-1,iez)%temp-T0)
          end do
       end do
       do iex=1,nex
          do iey=1,ney
             ecell(iex,iey,1)%temp=T0+K0T*(ecell(iex,iey,2)%temp-T0)
             ecell(iex,iey,nez)%temp=T0+K0T*(ecell(iex,iey,nez-1)%temp-T0)
          end do
       end do
    end select
  end subroutine borders

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nox_2_nex (ko,ixyze)
    integer,intent(in)::ko
    integer,intent(out)::ixyze(3)
    integer::kx,ky,kz,koc
    koc=ko
    !    write(6,*)'ko',ko
    kx=mod(koc-1,nox)+1
    koc=(koc-kx)/nox
    !    write(6,*)'kx koc',kx,koc
    ky=mod(koc,noy)+1
    !    write(6,*)'ky,koc',ky,(koc-ky+1)/noy
    kz=(koc-ky+1)/noy+1
    !    write(6,*)'kz',kz

    !    write(6,*)'ko,kx,ky,kz'
    !    write(6,*)ko,kx,ky,kz

    ixyze(1)=1+(kx-1)*nexov/nox+nexmp
    ixyze(2)=1+(ky-1)*neyov/noy+neymp
    ixyze(3)=1+(kz-1)*nezov/noz+nezmp


    !    write(6,*)'           ',ixyze

  end subroutine nox_2_nex

  subroutine xp_2_ecell(pos,ixyze)
    real(double),intent(in)::pos(3)
    integer,intent(out) :: ixyze(3)

    ixyze(1)=int(pos(1)/nexov)+1+nexmp
    ixyze(2)=int(pos(2)/neyov)+1+neymp
    ixyze(3)=int(pos(3)/nezov)+1+nezmp
  end subroutine xp_2_ecell


  subroutine CeT(C,T)
    real(double),intent(in) :: T
    real(double),intent(out) ::C 
    if (Cec.lt.0) then
       C=CedT(int(T))
    else
       C=CeC
    end if
  end subroutine CeT
  subroutine KeT(K,T)
    real(double),intent(in) :: T
    real(double),intent(out) ::K
    if (KeC.lt.0) then
       K=KedT(int(T))
    else
       K=KeC
    endif
  end subroutine KeT

  subroutine GepT(G,T)
    real(double),intent(in) :: T
    real(double),intent(out) ::G
    if (GepC.lt.0) then
       G=GepdT(int(T))
    else
       G=GepC
    endif
!    write(6,*) 'g ',g
  end subroutine GepT

  subroutine prepCe

    integer::it,itr,ittc,iex,iey,iez
    real(double), allocatable:: cerf(:),temprf(:)
    real(double)::factT
    allocate(cerf(0:ncer))
    allocate(CedT(0:nTmax))
    allocate(temprf(0:ncer))
    if (rang==0) write(6,*) 'Ce read from Ce.in ! ATTENTION AUX UNITES 1O^5Jm-3K-1'
    open(unit=84, file='Ce.in',form='formatted')
    cerf(0)=0 ; temprf(0)=0
    do itr=1,ncer
       read(84,*)temprf(itr),cerf(itr)
    end do
    close(84)
    temprf(:)=temprf*1d4
    cerf(:)=cerf*1d5*joule2erg*1d-6
    itr=0
    do it=1,ntmax
       !1      continue
       !       write(6,*)'it itr temprf(itr)',it, itr, temprf(itr)
       do while (float(it).gt.temprf(itr))
          itr=itr+1
       end do
       factT=(it-temprf(itr-1))/(temprf(itr)-temprf(itr-1))
       CedT(it)=cerf(itr-1)+(cerf(itr)-cerf(itr-1))*factT

    end do
    !   write(6,*)'it,itr',it,itr

    Eedt(0)=0.
    do it=1,nTmax
       Eedt(it)=EedT(it-1)+CedT(it)
       !       write(6,'(A,I5,G15.7)')'T eedt',it,eedt(it)
    end do
    if (integrTtype==2) then
       do iex=2,nex-1
          do iey=2,ney-1
             do iez=2,nez-1
                ittc=INT(ecell(iex,iey,iez)%temp)
                ecell(iex,iey,iez)%En=(Eedt(ittc)+(ecell(iex,iey,iez)%temp-ittc)*(Eedt(ittc+1)-Eedt(ittc)))*Vecell
                !                write(6,*)'%En',iex,iey,iez,ecell(iex,iey,iez)%En
             end do
          end do
       end do
    end if
    deallocate(cerf)
    return
  end subroutine prepCe

  subroutine prepGep

    integer::it,itr
    real(double), allocatable:: Geprf(:),temprf(:)
    real(double)::factT
    allocate(Geprf(0:ncer))
    allocate(GepdT(0:nTmax))
    allocate(temprf(0:ncer))
    if (rang==0) write(6,*) 'Gep read from Gep.in! ATTENTION AUX UNITES 1O^17Jm-3K-1'
    open(unit=84, file='Gep.in',form='formatted')
    Geprf(0)=0 ; temprf(0)=0
    do itr=1,ncer
       read(84,*)temprf(itr),Geprf(itr)
    end do
    close(84)
    temprf(:)=temprf*1d4
    Geprf(:)=Geprf*1d17*joule2erg*1d-6
    itr=0
    do it=1,ntmax
1      continue
       if (float(it).gt.temprf(itr))then
          itr=itr+1
          goto 1
       else
          factT=(it-temprf(itr-1))/(temprf(itr)-temprf(itr-1))
          GepdT(it)=Geprf(itr-1)+(Geprf(itr)-Geprf(itr-1))*factT
       end if
    end do
    deallocate(Geprf)
    return
  end subroutine prepGep


  subroutine prepKe
    integer::it,itr
    real(double), allocatable:: Kerf(:),temprf(:)
    real(double)::factT
    allocate(Kerf(0:ncer))
    allocate(KedT(0:nTmax))
    allocate(temprf(0:ncer))
    if (rang==0) write(6,*) 'Ke read from Ke.in'
    open(unit=84, file='Ke.in',form='formatted')
    Kerf(0)=0 ; temprf(0)=0
    do itr=1,ncer
       read(84,*)temprf(itr),Kerf(itr)
    end do
    close(84)
    temprf(:)=temprf*1d4
    Kerf(:)=Kerf*1d5*joule2erg*1d-6
    itr=0
    do it=1,ntmax
1      continue
       if (float(it).gt.temprf(itr))then
          itr=itr+1
          goto 1
       else
          factT=(it-temprf(itr-1))/(temprf(itr)-temprf(itr-1))
          KedT(it)=Kerf(itr-1)+(Kerf(itr)-Kerf(itr-1))*factT
       end if
    end do
    deallocate(Kerf)
    return

  end subroutine prepKe

  subroutine sauveelec
    integer:: luecout
    character :: fnamecout*80

    integer::iex,iey,iez,ic
    luecout=65
!    write(6,*)'IN SVEL'
    if (rang.ne.0)then
       write(6,*) 'WTF sauvE rang <>0!'
       stop
    end if
    
    fnamecout = fnam(1:lenfnam)//'.ecout'
    open (unit=luecout,file=fnamecout,form='unformatted')


    write(luecout)ecell

    write(luecout)Eelec
    write(luecout)Teavg
    close(luecout)
  end subroutine sauveelec

  subroutine restartelec

    character ::  fnamecin*80
    integer:: luecin,iex,iey,iez,ic
    luecin=66
    write(6,*)'in Erestart',nez,nex,nez
    if (lrestart) then
       fnamecin = fnam(1:lenfnam)//'.ecout'
    else 
       fnamecin = fnam(1:lenfnam)//'.ecin'
    end if
    open (unit=luecin,file=fnamecin,form='unformatted')

    read(luecin)ecell
    read(luecin)Eelec
    read(luecin) Teavg
    close (luecin)
    write(6,*)'out restart'
  end subroutine restartelec

  subroutine fillrbuf (tabe,buftab,nx,ny,nz,ir)
    integer:: ir,nx,ny,nz,ix,iy,iz
    real(double):: tabe(nx,ny,nz),buftab(nx,ny,nz)

    if (ir==1) then
       do ix=1,nx
          do iy=1,ny
             do iz=1,nz
                buftab(ix,iy,iz)=tabe(ix,iy,iz)
             end do
          end do
       end do
    elseif (ir==-1) then
       do ix=1,nx
          do iy=1,ny
             do iz=1,nz
                tabe(ix,iy,iz)=buftab(ix,iy,iz)
             end do
          end do
       end do
       tabe=buftab
    end if
    return
  end subroutine fillrbuf

  subroutine fillibuf (tabe,buftab,nx,ny,nz,ir)
    integer:: ir,nx,ny,nz
    integer::tabe(nx,ny,nz),buftab(nx,ny,nz)

    if (ir==1) then
       buftab=tabe
    elseif (ir==-1) then
       tabe=buftab
    end if
    return
  end subroutine fillibuf

  subroutine filllbuf (tabe,buftab,nx,ny,nz,ir)
    integer:: ir,nx,ny,nz
    logical::tabe(nx,ny,nz),buftab(nx,ny,nz)

    if (ir==1) then
       buftab=tabe
    elseif (ir==-1) then
       tabe=buftab
    end if
    return
  end subroutine filllbuf


  subroutine eleccellmol
    integer ::  luvisue,lenfn2,iex,iey,iez
    real(double)::exmM(3)
    character :: extension*9

    luvisue=91
    if(rang==0) then
       if (it < 0  )   extension='iiiiiiiii' 
       if (it >= 0 )   write(extension,'(i9.9)') it
       open(luvisue, file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.mol', form='formatted', &
            status='unknown')
       write(luvisue,'(I9,A,I9,A, F12.6)')nex*ney*nez, ' IT =', it, ' Time = ', timel

       exmm(:)=ecell(nex,ney,nez)%ixb(:)+cellside(:)-ecell(1,1,1)%ixb(:)
       write (luvisue,'(9F12.6)')exmm(1),0,0,0,exmm(2),0,0,0,exmm(3)

       do iex=1,nex
          do iey=1,ney
             do iez=1,nez
                write(luvisue,'(5G15.7,2I5)')ecell(iex,iey,iez)%ixb(1),ecell(iex,iey,iez)%ixb(2),ecell(iex,iey,iez)%ixb(3),&
                     &ecell(iex,iey,iez)%temp,ecell(iex,iey,iez)%tempIon,ecell(iex,iey,iez)%nion,ecell(iex,iey,iez)%nion-ecell(iex,iey,iez)%nionS
             end do
          end do
       end do
       close(luvisue)
    end if

  end subroutine eleccellmol
end module elec_cell
