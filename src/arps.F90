module arps_mod
  USE T_kind_param_m,only:double,long
  USE arret_ndm_mod,only:arret_ndm
  USE analyseT_mod,only: analyseT
  USE controleT_mod,only: controleT
  USE boxconfig,only:box_config,periodbox
  USE atomconfig,only : atom_config_arps
  USE cellconfig, only:cell_config,caltabtc,cell_config_arps
  use Tpara,only:para_space_config
  USE var_pot, ONLY:ipotentiel,ntyp,potis1,potiseam,potisglue,potisrep,zz
  USE calfo2ccel_mod,only:calfo2ccel,sig2p
  USE calfoeamcel_mod,only:calfoeamcel
  USE calfoew_mod,only:calfozz
  USE gen_com_m, ONLY:potist,rang,sig,lspaceNDM,itmax,itloopmax,timemax,timeloopmax,latcomp,iteration,itesigma,timel,tstep,&
       &lperiod,sigkine,ltpcel,sigtot,potis2,bk,erg2ev,itetemp,ltberendsen,tautcon,text,lspacendm,dmtype,unitP,llangevin,pi,&
       &iterasmol,itetimestep
  USE deftimestep_mod,only: deftimestep  
  use calfocommon,only:sigcalfo,potistcalfo,test_sigma,sigcalfo,sigc,lcalcsigc
  use var_pot,only : cm,iewald,gamlt
  use vect_dist_mod,only:vect_dist
  USE calfoberend_mod,only: calfoberend
  use tempinstT_mod,only:tempinstT
  use calfoeamcel_mod,only:calfoeamcel
  USE arret_ndm_mod,only:arret_ndm
  USE calcfvp_mod
  USE eloss, ONLY : calceloss,ibrake 

#ifdef PARA
  use Tpara,only:nprocspace,para_space_config,comm_space,myidsp
  USE mod_para,only:maj_tabdensity_ftm,maj_atomes_frt_ftm,maj_atomes_frt_part


#else
  USE Tpara,only:nprocspace,para_space_config
#endif

  implicit none  


!!$  type, extends (atom_config_d)::atom_config_arps
!!$     real(double),allocatable::rho(:),fpr(:,:)
!!$     integer, allocatable::movxyz(:)
!!$  end type atom_config_arps

  real(double)::kmin,kmax
!  real(double),allocatable,dimension(:)::aspl,bspl,vmin,vmax,S0spl,xs
  real(double),allocatable,dimension(:) :: tabdensity
  !  logical::lxyz
  logical::lpartarps
  real(double)::sigem(3,3)  !sigép of calfo2ccel is ARTIFICIALLY used in the EAM case
  integer::noxyzkmin(3),noxyzkmax(3)
contains
  subroutine dmloop_arps(atdml,celndm,boxndm,psc)

    type(para_space_config)::psc
    class(box_config)::boxndm
    type(atom_config_arps)::atdml
    type(cell_config_arps),target:: celndm
    logical :: lreturn
    integer::imm,i,ilocal,ic,iti

    real(double), dimension(ntyp) :: aux
    real(double)::potisrep0,fvp,m,vn,xpar,u1,u2,u3,fvp2
    !    real(double)::tabtat(500,ntyp+1)
    if (rang==0) write (6, *) '***** FIRST ITERATION  ARPS****',itloopmax,timeloopmax,itesigma
    ! Appel de la routine generale des forces

    imm =atdml%imm


    atdml%lgul=.true.
    test_sigma=.true.
    if((test_sigma).and.(celndm%ltpcel))then
       lcalcsigc=.true.
       sigc=>celndm%sigc
    end if

    select case (ipotentiel)
    case(0,1,3,4,5,6,7,8,9)

       if ((iewald.gt.0).and.(iewald.ne.3))then
          write(6,*)'ewald not coded with arps stop'
          call arret_ndm
       end if
       atdml%fp=0
       potis1=0
       if(ltpcel)then
          sigc=>celndm%sigc
          sigc=0
       end if
       call calfo2ccel(atdml,celndm%cell_config,boxndm)
       sig=sig2P
       potist=potis1
       atdml%fpr=atdml%fp
    case(10,11,16)
       atdml%fp=0
       atdml%fpr=0
       potist=0
       if (ipotentiel==16) then
          write(6,*)'check algo ipotentiel=16, charge effects'
          call arret_ndm
       end if
       tabdensity(:)=0
       !       atdml%fpr=atdml%fpr-atdml%fpg ! fpr reduced to rep only
       !       potist=potist-potisglue
       sigcalfo=0
       select case (dmtype)
       case(41)
          call calfoglue_arps_1(atdml,celndm,boxndm)
          atdml%fpr=atdml%fpr+atdml%fp ! fpr =all active rep
          potist=potist+potisrep
          potisrep0=potisrep
          atdml%rho(1:atdml%im)=atdml%rho(1:atdml%im)+tabdensity(1:atdml%im) ! %rho +new active rho
          atdml%fp=0
          call calfoglue_arps_2(atdml,celndm,boxndm,psc)
          atdml%fpg=atdml%fp
          atdml%fpr=atdml%fpr+atdml%fpg ! total force =rep+glue
          potist=potist+potisglue
          sig=sig2p+sigem
       case(42)
          atdml%fp=0
          if(ltpcel)then
             sigc=>celndm%sigc
             sigc=0
          end if
          call calfoeamcel(atdml,celndm%cell_config,boxndm,psc)
          potist=potiseam
          atdml%fpr=atdml%fp
          sig=sigcalfo
       end select
!!$#endif


    end select
    if(ibrake.gt.0) call calceloss(celndm,atdml)
    atdml%mov=2
    call analyseT (atdml,celndm,boxndm,psc)
    call periodbox (boxndm,atdml)
    !A.2
#ifdef PARA    
    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call caltabtC(celndm%cell_config,atdml,lperiod,boxndm,psc=psc)
    else
       call caltabtC(celndm%cell_config,atdml,lperiod,boxndm)
    end if
#else
    call caltabtC(celndm%cell_config,atdml,lperiod,boxndm)
#endif

#ifdef PARA

    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
       !       select case (ipotentiel)
       !       case(0,1,3,4,5,6,7,8,9)
       call maj_atomes_frt_ftm(atdml,celndm%cell_config,boxndm,psc)
       !       write(6,*)'rang nat ',rang, iteration, atdml%im,atdml%imf
       !       end select
    end if
#endif
    !    call analyseT (atdml,celndm,boxndm,psc)
!!!!!!!!!!!!!!!LOOP START
    do while ((iteration.lt.itloopmax).and.(timel.lt.timeloopmax))

       iteration = iteration+1
       test_sigma=(mod(iteration,itesigma)==0)

       timel = timel+tstep
       aux(:ntyp) = tstep/cm(:ntyp)/2.d0       
       !B.1
       if(llangevin) then
          DO i=1, atdml%im
             do ic=1,3
                call random_number(u1)
                call random_number(u2)
                atdml%Glangv(ic,i)=sqrt(-2.*log(u1))*cos(2.*pi*u2)
             end do
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !             atdml%mov=2
             call calcfvp(fvp,atdml%ityp(i),atdml%vp(:,i),atdml%mov(i))
!!$             if (atdml%mov(i)==2) then
!!$                fvp=1
!!$             else
!!$                iti=atdml%ityp(i)
!!$                m=cm(iti)
!!$                vn=norm2(atdml%vp(:,i))
!!$                xpar=vn*m-vmin(iti)*m
!!$!                write(6,*)'xpar',xpar,xs(iti)
!!$                if (xpar.le.0)  then
!!$                   fvp=0
!!$                else if ((xpar.lt.xs(iti)).and.(xpar.gt.0)) then
!!$                   fvp=pol(xpar,aspl(iti),bspl(iti))/vn
!!$                else
!!$                   fvp=1
!!$!                   write(6,*)'POOOOO'
!!$!                   call arret_ndm
!!$                end if
!!$             end if
!!$             write(6,*)'factg', i, fvp,fvp2
             do ic=1,3
                !  write(6,*)'ct',cm(ityp(1)),tstep
!!$                write(6,*)'1-rga', fvp*gamlt(atdml%ityp(i))*tstep/2
!!$             write(6,*)'2',atdml%fpr(ic,i)*tstep/(cm(atdml%ityp(i))*2)
!!$             write(6,*)'3',atdml%Glangv(ic,i)*sqrt(cm(atdml%ityp(i))*bk*text*gamlt(atdml%ityp(i))*tstep*0.5)/cm(atdml%ityp(i))
             
                atdml%vp(ic,i) = atdml%vp(ic,i)*(1 -fvp*gamlt(atdml%ityp(i))*tstep/2) &
                     & + atdml%fpr(ic,i)*tstep/(cm(atdml%ityp(i))*2)&
                     &+atdml%Glangv(ic,i)*sqrt(cm(atdml%ityp(i))*bk*text*gamlt(atdml%ityp(i))*tstep*0.5)/cm(atdml%ityp(i))
             end do
          END DO
       else
          DO i=1, atdml%im
             atdml%vp(1:3,i) = atdml%vp(1:3,i) + aux(atdml%iTyp(i))*atdml%fpr(1:3,i)
          END DO
       end if
#ifdef PARA


       if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
          call maj_atomes_frt_part(atdml,celndm%cell_config,boxndm,psc,'v')          
       end if
#endif
       !B.2
       call setfree(atdml,celndm)
       ! TO CHECK    atdml%mov(:)=2 ;       atdml%lgul(:)=.true.

       !B.3
       select case (dmtype)
       case(41)
          select case (ipotentiel)
          case(0,1,3,4,5,6,7,8,9)
             if(ltpcel) then
                write(6,*)'LTPCEL pas programme pour pressions dmtype=41 et pot pair faire a l image de sigcalfo'
             end if
             atdml%fp=0 ; potis1=0;
             call calfo2ccel(atdml,celndm%cell_config,boxndm)
             sig=sig-sig2P
             atdml%fpr=atdml%fpr-atdml%fp
             potist=potist-potis1
          case(10,11)
             atdml%fp=0 
             potisrep=0
             tabdensity(:)=0
             atdml%fpr=atdml%fpr-atdml%fpg ! fpr reduced to rep only
             potist=potist-potisglue
             sig=sig-sigem
             call calfoglue_arps_1(atdml,celndm,boxndm)
             atdml%fpr=atdml%fpr-atdml%fp ! fpr (=rep) - old active rep
             potist=potist-potisrep
             potisrep0=potisrep0-potisrep
             sig=sig-sig2p
             atdml%rho(1:atdml%im)=atdml%rho(1:atdml%im)-tabdensity(1:atdml%im) ! %rho -old active rho
          end select
       case(42)
!!$          select case (ipotentiel)
!!$          case(0,1,3,4,5,6,7,8,9)
!!$             atdml%fp=0 ; potis1=0
!!$             call calfo2ccel(atdml,celndm,boxndm)
!!$             atdml%fpr=atdml%fp
!!$             potist=potis1
!!$          case(10,11)
!!$             atdml%fp=0
!!$             tabdensity=0
!!$             call calfoeamcel(atdml,celndm,boxndm,psc)
!!$             potist=potiseam
!!$             atdml%fpr=atdml%fp
!!$          end select
       end select
       if(ibrake.gt.0) call calceloss(celndm,atdml)



       !A.1
       call xpupdate(atdml)



       call periodbox (boxndm,atdml)
       !A.2
#ifdef PARA    
       if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
          call caltabtC(celndm,atdml,lperiod,boxndm,psc=psc)
       else
          call caltabtC(celndm,atdml,lperiod,boxndm)
       end if
#else
       call caltabtC(celndm,atdml,lperiod,boxndm)
#endif

#ifdef PARA
       if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
          ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
          !       select case (ipotentiel)
          !       case(0,1,3,4,5,6,7,8,9)
          call maj_atomes_frt_ftm(atdml,celndm%cell_config,boxndm,psc)
          !       end select
       end if
#endif
       call caltabtarps( celndm,atdml,lperiod,boxndm,psc)
       !!A.3
       select case(dmtype)
       case(41)
          select case (ipotentiel)
          case(0,1,3,4,5,6,7,8,9)
             atdml%fp=0; potis1=0
             call calfo2ccel(atdml,celndm%cell_config,boxndm)
             sig=sig+sig2P
             atdml%fpr=atdml%fpr+atdml%fp
             potist=potist+potis1
          case(10,11)
             atdml%fp=0 
             potisrep=0
             tabdensity(:)=0
             call calfoglue_arps_1(atdml,celndm,boxndm)
             atdml%fpr=atdml%fpr+atdml%fp ! fpr (=rep) + new active rep
             potist=potist+potisrep
             potisrep0=potisrep0+potisrep
             potisrep=potisrep0
             sig=sig+sig2p
             atdml%rho(1:atdml%im)=atdml%rho(1:atdml%im)+tabdensity(1:atdml%im) ! %rho +new active rho
             ! fpr contains the complete repulsion and %rho contains the complete density
             !A.4 calculation of glue force
             atdml%fp=0
             call calfoglue_arps_2(atdml,celndm,boxndm,psc)
             atdml%fpg=atdml%fp
             atdml%fpr=atdml%fpr+atdml%fpg ! total force =rep+glue
             sig=sig+sigem

             potist=potist+potisglue
          end select
       case(42)
          select case (ipotentiel)
          case(0,1,3,4,5,6,7,8,9)
             atdml%fp=0 ; potis1=0
             if(ltpcel)sigc=0
             call calfo2ccel(atdml,celndm%cell_config,boxndm)
             atdml%fpr=atdml%fp
             potist=potis1
          case(10,11)
             sigcalfo=0
             atdml%fp=0
             tabdensity=0
             if(ltpcel)sigc=0
             call calfoeamcel(atdml,celndm%cell_config,boxndm,psc)
             potist=potiseam
             atdml%fpr=atdml%fp
             sig=sigcalfo
!!$          write(6,*)'SIG',test_sigma,Sig(1,1)*unitP,Sig2p(1,1)*unitP,Sigem(1,1)*unitP
          end select
       end select
       if (lTberendsen) then
          block
            real(double)::tempm1,gamb,fact
            tempm1=tempinstT(atdml)

            !      write(6,*)'jy suis'
            gamb=1./(2.*tauTcon)
            !      write(6,*)gamb,text,tempm1

            do i=1,atdml%im
               fact=cm(atdml%ityp(i))*gamb*(Text/tempm1-1.0)
               do ic=1,3
                  !            write(6,*)fp(ic,i),fact*vp(ic,i)
                  atdml%fpr(ic,i)=atdml%fpr(ic,i)+fact*atdml%vp(ic,i)
               end do

            end do
          end block
       end if

       !O
       if(llangevin) then
          DO i=1, atdml%im
!             call random_number(u1)
!             call random_number(u2)
!             atdml%Glangv(ic,i)=sqrt(-2.*log(u1))*cos(2.*pi*u2)   
             call calcfvp(fvp,atdml%ityp(i),atdml%vp(:,i),atdml%mov(i))
             do ic=1,3
                !  write(6,*)'ct',cm(ityp(1)),tstep
!!$                write(6,*)'II-rga', fvp*gamlt(atdml%ityp(i))*tstep/2
!!$             write(6,*)'2',atdml%fpr(ic,i)*tstep/(cm(atdml%ityp(i))*2)
!!$             write(6,*)'3',atdml%Glangv(ic,i)*sqrt(cm(atdml%ityp(i))*bk*text*gamlt(atdml%ityp(i))*tstep*0.5)/cm(atdml%ityp(i))

                atdml%vp(ic,i) = atdml%vp(ic,i)*(1 -fvp*gamlt(atdml%ityp(i))*tstep/2) &
                     & + atdml%fpr(ic,i)*tstep/(cm(atdml%ityp(i))*2)&
                     &+atdml%Glangv(ic,i)*sqrt(cm(atdml%ityp(i))*bk*text*gamlt(atdml%ityp(i))*tstep*0.5)/cm(atdml%ityp(i))
             end do
          END DO
       else
          DO i=1, atdml%im
             atdml%vp(1:3,i) = atdml%vp(1:3,i) + aux(atdml%iTyp(i))*atdml%fpr(1:3,i)
          END DO
       end if

!       DO i=1, atdml%im
!          atdml%vp(1:3,i) = atdml%vp(1:3,i) + aux(atdml%iTyp(i))*atdml%fpr(1:3,i)
!!$          Tat=0.5*cm(atdml%ityp(i))*(atdml%vp(1,i)**2+atdml%vp(2,i)**2+atdml%vp(3,i)**2)*0.66666/bk
!!$          j=1+int(tat/10.)
!!$          tabtat(j,atdml%ityp(i))=tabtat(j,atdml%ityp(i))+1
!!$          tat=0.5*cm(atdml%ityp(i))*(atdml%vp(1,i)**2+atdml%vp(2,i)**2+atdml%vp(3,i)**2)*erg2eV
!!$          j=1+int(tat*10.)
!!$          tabtat(j,ntyp+1)=tabtat(j,ntyp+1)+1

!       END DO
       ! les positions et les vitesses sont synchrones en ce point ; les atomes sont bien r�partis en cellules

       if (test_sigma) then

          sigkine=0.
          do ilocal = 1, atdml%im
             sigkine(1:3,1) = sigkine(1:3,1) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(1,ilocal)
             sigkine(1:3,2) = sigkine(1:3,2) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(2,ilocal)
             sigkine(1:3,3) = sigkine(1:3,3) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(3,ilocal)
             if (lTPcel.EQV..true.) then
                celndm%sigc(1:3,1,atdml%ielat(ilocal)) = celndm%sigc(1:3,1,atdml%ielat(ilocal)) + &
                     cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(1,ilocal)*celndm%noxyz/boxndm%volu
                celndm%sigc(1:3,2,atdml%ielat(ilocal)) = celndm%sigc(1:3,2,atdml%ielat(ilocal)) + &
                     cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(2,ilocal)*celndm%noxyz/boxndm%volu
                celndm%sigc(1:3,3,atdml%ielat(ilocal)) = celndm%sigc(1:3,3,atdml%ielat(ilocal)) + &
                     cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(3,ilocal)*celndm%noxyz/boxndm%volu
             end if
          end do
          sigkine(1:3,1:3) = sigkine(1:3,1:3)/boxndm%volu

#ifdef PARA

          if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
             call comm_space%sum(sigkine)
             if (allocated(celndm%sigc)) then
                call comm_space%sum(celndm%sigc)
             end if
          end if
#endif
          sigtot = sigkine+sig
       end if

       call analysearps (atdml,celndm,boxndm,psc,'P2')

       call controleT(atdml,celndm%cell_config,boxndm,psc,lreturn)
       if (lreturn) return

    end do
!!$    do i=1,500
!!$       write(111,*)i*10,tabtat(i,1),tabtat(i,2),tabtat(i,3)
!!$    end do

  end subroutine dmloop_arps

  subroutine analysearps (atdml,celndm,boxndm,psc,chr)
    !    use gen_com_m,only:iterasmol
    USE rasmolT_mod,only: rasmolT
    type(para_space_config)::psc
    class(box_config)::boxndm
    type(atom_config_arps)::atdml
    type(cell_config_arps):: celndm
    character(len=*)::chr

    character*3,allocatable,target::tymov(:)

    real(double)::earps,kinps

    CHARACTER(len=89) :: namemov
    integer::i
    namemov='mov'

    if (itetemp>0) then
       if (mod(iteration,itetemp)==0) then
          kinps=kinarps(atdml)
          earps=potist+kinps
          if (rang==0) write(6,'(2A,I7,3E15.6)')'EARPS ',chr,iteration,timel, earps*erg2ev,kinps*erg2ev
       end if
    end if

    call analyseT (atdml,celndm,boxndm,psc)
    if (itetimestep>0) then
       if (mod(iteration,itetimestep)==0) call deftimestep(atdml,boxndm)
    endif

!!$    if (iterasmol>0) then     
!!$       if (mod(iteration,iterasmol)==0) then
!!$          allocate(tymov(atdml%im))
!!$          do i=1,atdml%im
!!$             select case (atdml%mov(i))
!!$             case(0)
!!$                tymov(i)=' Re'
!!$             case(1)
!!$                tymov(i)=' In'             
!!$             case(2)
!!$                tymov(i)=' Mo'
!!$             end select
!!$          end do
!!$          
!!$          call rasmolT(atdml,boxndm,iteration,latcomp=latcomp,rty=tymov,namefr='mov')
!!$!          if (l2T) call  eleccellmol
!!$
!!$       end if
!!$    endif
  end subroutine analysearps
  function kinarps(atdml)
    type(atom_config_arps)::atdml
    real(double)::kinarps
    real(double)::m,pmin,x,kin1,kin2,kin3
    integer::i,im,iti
    im =atdml%im
    kinarps=0.0
    kin1=0;kin2=0.;kin3=0
    do i=1,im
       iti=atdml%ityp(i)
       m=cm(iti)
       select case (atdml%mov(i))
       case(0)
          kin1=kin1+S0spl(iti)
       case(1)
          x=(norm2(atdml%vp(:,i))-vmin(iti))*m
          pmin=vmin(iti)*m
          kin2=kin2+uf(x,aspl(iti),bspl(iti),s0spl(iti),xs(iti),pmin,m)
       case(2)
          kin3=kin3+0.5*m*(norm2(atdml%vp(:,i)))**2
       end select

    end do
    kinarps=kiN1+kiN2+kiN3
#ifdef PARA

    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call comm_space%sum(kinarps)
    end if
#endif
    sigtot = sigkine+sig

    !    write(6,'(A,4G17.8)')'KINKIN ',kin1*erg2ev,kin2*erg2ev,kin3*erg2ev,kinarps*erg2ev
  end function kinarps
  !----------------------------------------------------------------------
  SUBROUTINE calfoglue_arps_1(atcf,celcf,boxcf)
    USE T_kind_param_m
    USE gen_com_m, ONLY:angst,low_limit,zero,pi,rang
    !  USE calfocommon

    !  USE cellconfig, only : cell_config
    !  use boxconfig,only: box_config

    USE var_pot, ONLY:ipotentiel,ngrid,potisrep,rue_pot,&
         &typ_and_pot,typ_pot_pair,ipotentiel,ngrid,potisrep,rhomax,rhomin,eamrho,ipo,eamrep,&
         &alpha,ntyp

#ifdef PARA
    use Tpara,only:nprocspace,para_space_config,comm_space

#else
    USE Tpara,only:nprocspace,para_space_config
#endif

    type(atom_config_arps),intent(inout)::atcf
    type(cell_config_arps),intent(in)::celcf
    type(box_config),intent(in)::boxcf

    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    ! eam variables
    !local variables
    integer :: i,j !atomes
    integer ::iti,itj !types
    integer :: l !paires
    integer:: koo,ko1,i2,i1 !cel.
    integer :: k ! aux pour splines

    REAL(double), dimension(1:3) :: dxp, gradij
    real(double) :: r!distance i-j
    real(double) :: dEembi ! potentiel et gradient de l'immersion
    real(double) :: rhoi,rhoj ! densite de i sur j et j sur i
    real(double) :: Erep,dErep ! potentiel et gradient de la repulsion de paire ij

    real(double) ::  drk,ktor, inv_ktor
    real(double),dimension(:),allocatable::ktorho(:), inv_ktorho(:)
    real(double) :: densityi !densite totale sur i
    real(double)::rue,alp,aux
    logical ::linter
    sig2p=0
    rue=rue_pot(ipotentiel)
    aux = 23.06134575D-20
    alp = alpha/sqrt(pi)*aux
    allocate(ktorho(ntyp))
    allocate(inv_ktorho(ntyp))
    ktor=rue/ngrid
    inv_ktor=1.d0/ktor
    ktorho(:)=(rhomax(:)-rhomin(:))/ngrid
    inv_ktorho(:) = 1.d0/ktorho(:)
    !    rue2=rue**2

    loop1at1: do i=1,atcf%im
       if (typ_and_pot(atcf%ityp(i),ipotentiel).eqv..false.)cycle

       ! --- Calcul de la densite sur i ---    

       !     nvi=0
       densityi=0.0 ; dEembi=0.0
       koo = atcf%ielat(i)                          ! Numero de la cellule
       iti = atcf%ityp(i)

       ! pour chaque cel. voisine
       loop1cel:   do i1 = 0, celcf%ncelvois(koo)
          ko1 = celcf%ncel(koo,i1)
          ! pour chaque atome ds la cel. voisine
          loop1at2: do i2 = 1, celcf%nato(ko1)
             j = celcf%atincel(i2,ko1)
             !             write(6,*)i,koo,i1, celcf%ncelvois(koo),i2,j
             if (typ_pot_pair(ipo(atcf%ityp(i),atcf%ityp(j))).ne.ipotentiel) cycle

             itj=atcf%ityp(j)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !CRC             if(i.eq.j) cycle
#ifdef PARA

             ! Methode pour ne prendre qu'une seule fois en compte
             ! le couple i,j en paralle :
             ! - i est necesairement local (boucle i<=im)
             ! - si j est local on ne retient que le couple i<j
             ! - si j n'est pas local, le couple n'est par definition
             !   pris qu'une fois puisque i est local
             if (nprocspace.gt.1) then

                if (j.le.atcf%im) then

                   ! les deux atomes sont locaux
                   if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme deja calcule
                else
                   ! j n'est pas local, on ne fait le calcul normal           
                endif
             else
                if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme dÃ£Â©ja calculÃ£Â©
             end if

#else
             if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme dÃ£Â©ja calculÃ£Â©
#endif
             if((.not.atcf%lgul(i)).and.(.not.atcf%lgul(j)))cycle ! both restrained, no change in density nor rep force

             call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp,indcv=i1, lperiod=boxcf%lperiod,rum=rue,linter=linter,dist=r)
             if(.not.linter) cycle
             k=Int(r*inv_ktor)
             gradij(1:3) = dxp(1:3)/r
             drk=r-k*ktor
             !             rhoj = eamrho(1,itj,k) + drk*( eamrho(2,itj,k) + drk*( eamrho(3,itj,k) + drk*eamrho(4,itj,k) ) )  !rho de j sur i
             !             densityi=densityi+rhoj

             rhoj = eamrho(1,itj,k) + drk*( eamrho(2,itj,k) + drk*( eamrho(3,itj,k) + drk*eamrho(4,itj,k) ) )  !rho de j sur i
             tabdensity(i)=tabdensity(i)+rhoj
             rhoi = eamrho(1,iti,k) + drk*( eamrho(2,iti,k) + drk*( eamrho(3,iti,k) + drk*eamrho(4,iti,k) ) )  !rho de i sur j
             tabdensity(j)=tabdensity(j)+rhoi
!!$             write(100,'(2I3,3G17.8)')i,j,r,tabdensity(i),tabdensity(j)
!!$             write(100,'(2I3,4G17.8)')i,j,eamrho(1,iti,k) , eamrho(2,iti,k),eamrho(3,iti,k),eamrho(4,iti,k)
!!$             write(100,'(2I3,4G17.8)')i,j,eamrho(1,itj,k) , eamrho(2,itj,k),eamrho(3,itj,k),eamrho(4,itj,k)
             !           nvi=nvi+1
             !           dxpij(1:3,nvi)=dxp(1:3)
             !           jvi(nvi)=j
             !           rij(nvi)=r            

             !terme de repulsion 
             l = ipo(iti,itj)
             Erep = eamrep(1,l,k) + drk*( eamrep(2,l,k) + drk*( eamrep(3,l,k) + drk*eamrep(4,l,k) ) )
!!$             if(lprteat.EQV..true.)then
!!$                select type (atcf)
!!$                class is (atom_config_e)
!!$                   atcf%eat(i)=atcf%eat(i)+Erep/2.d0
!!$                   if (j.le.atcf%im) atcf%eat(j)=atcf%eat(j)+Erep/2.d0
!!$                end select
!!$             end if
             dErep = eamrep(2,l,k) + drk*( 2.0*eamrep(3,l,k) + 3.0*drk*eamrep(4,l,k) )

             if (atcf%num_at_glob(i).lt.atcf%num_at_glob(j)) then
                potisrep = potisrep+Erep
             endif

             atcf%fp(1:3,i)=atcf%fp(1:3,i)-dErep*gradij(1:3)
             atcf%fp(1:3,j)=atcf%fp(1:3,j)+dErep*gradij(1:3)

             if (test_sigma) then        
                if (atcf%num_at_glob(i).lt.atcf%num_at_glob(j)) then          
                   sig2p(1:3,1) = sig2p(1:3,1)-dErep*gradij(1:3)*dxp(1)/boxcf%volu
                   sig2p(1:3,2) = sig2p(1:3,2)-dErep*gradij(1:3)*dxp(2)/boxcf%volu
                   sig2p(1:3,3) = sig2p(1:3,3)-dErep*gradij(1:3)*dxp(3)/boxcf%volu
!!$                   if (lTPcel.EQV..true.) then
!!$                      sigc(1:3,1,koo) =sigc(1:3,1,koo) -0.5*dErep*gradij(1:3)*dxp(1)*celcf%noxyz/boxcf%volu
!!$                      sigc(1:3,2,koo) =sigc(1:3,2,koo) -0.5*dErep*gradij(1:3)*dxp(2)*celcf%noxyz/boxcf%volu
!!$                      sigc(1:3,3,koo) =sigc(1:3,3,koo) -0.5*dErep*gradij(1:3)*dxp(3)*celcf%noxyz/boxcf%volu
!!$                      sigc(1:3,1,ko1) =sigc(1:3,1,ko1) -0.5*dErep*gradij(1:3)*dxp(1)*celcf%noxyz/boxcf%volu
!!$                      sigc(1:3,2,ko1) =sigc(1:3,2,ko1) -0.5*dErep*gradij(1:3)*dxp(2)*celcf%noxyz/boxcf%volu
!!$                      sigc(1:3,3,ko1) =sigc(1:3,3,ko1) -0.5*dErep*gradij(1:3)*dxp(3)*celcf%noxyz/boxcf%volu
!!$                   end if
!!$
                endif
             end if



          end do loop1at2
       end do loop1cel
    end do loop1at1
#ifdef PARA

    if (nprocspace.gt.1) then
       call comm_space%sum(potisrep)
       !       call comm_space%sum(potisglue)
       if (test_sigma) then 
          !          call comm_space%sum(sig)
          !          call comm_space%sum(sig2p)
          call comm_space%sum(sigem)

          !          if (associated(sigc)) then
          !             call comm_space%sum(sigc)
          !       endif

          !          call comm_space%sum(sig)
          call comm_space%sum(sig2p)


       endif


    endif


#endif


    return
  end SUBROUTINE calfoglue_arps_1


  SUBROUTINE calfoglue_arps_2(atcf,celcf,boxcf,psc)
    USE T_kind_param_m
    USE gen_com_m, ONLY:angst,low_limit,zero,pi
    !  USE calfocommon
    !  use vect_dist_mod,only:vect_dist
    !  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
    !  USE cellconfig, only : cell_config
    !  use boxconfig,only: box_config

    USE var_pot, ONLY:ipotentiel,ngrid,potisglue,rue_pot,&
         &typ_and_pot,typ_pot_pair,ipotentiel,ngrid,rhomax,rhomin,eamrho,ipo,&
         &eamglue,alpha,ntyp


    type(atom_config_arps),intent(inout)::atcf
    type(cell_config_arps),intent(in)::celcf
    type(box_config),intent(in)::boxcf

    type(para_space_config)::psc
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    ! eam variables
    !local variables
    integer :: i,j !atomes
    integer ::iti,itj !types

    integer:: koo,ko1,i2,i1 !cel.
    integer :: k ! aux pour splines

    !    real(double) :: rue2 !coupure**2
    REAL(double), dimension(1:3) :: dxp, gradij
    real(double) :: r!distance i-j
    real(double) :: Eembi ! potentiel et gradient de l'immersion


    REAL(double) :: Femb
    real(double) ::  drk,ktor, inv_ktor
    real(double),dimension(:),allocatable::ktorho(:), inv_ktorho(:)
    real(double) :: densityi !densite totale sur i
    real(double) :: tabdensity(atcf%imm)
    real(double)::rue,alp,aux
    logical ::linter

    sigem=0
    rue=rue_pot(ipotentiel)
    aux = 23.06134575D-20
    alp = alpha/sqrt(pi)*aux
    allocate(ktorho(ntyp))
    allocate(inv_ktorho(ntyp))
    ktor=rue/ngrid
    inv_ktor=1.d0/ktor
    ktorho(:)=(rhomax(:)-rhomin(:))/ngrid
    inv_ktorho(:) = 1.d0/ktorho(:)
    potisglue=0
    tabdensity(:)=0

    ! calcul et stockage de Eembi et dEembi
    loop2at1: do i=1,atcf%im
       if (typ_and_pot(atcf%ityp(i),ipotentiel).eqv..false.)cycle
       iti=atcf%ityp(i)
       k=Int((atcf%rho(i)-rhomin(iti))*inv_ktorho(iti))
       !       write(6,*)i,iti,k,tabdensity(i),rhomin(iti),inv_ktorho(iti)
       !       write(110,'(2I8,3G17.8)')i,k,tabdensity(i),rhomin(iti),inv_ktorho(iti)
       if(k.gt.ngrid) then
          write(6,*)k, ngrid, 'k> ngrid ; augmenter le facteur multiplicatif de rhomax dans calpo'
          write(6,*)'densityi',k,ngrid,densityi
          call arret_ndm
       end if
       drk=atcf%rho(i)-(rhomin(iti)+k*ktorho(iti))
       Eembi = eamglue(1,iti,k) + drk*( eamglue(2,iti,k) + drk*( eamglue(3,iti,k) + drk*eamglue(4,iti,k) ) )
       !       write(120,'(2I8,4G17.8)')i,k, eamglue(1,iti,k) , eamglue(2,iti,k),eamglue(3,iti,k),eamglue(4,iti,k)

       !       if(lprteat.EQV..true.)then
       !          select type (atcf)
       !          class is (atom_config_e)
       !             atcf%eat(i)=atcf%eat(i)+Eembi
       !          end select
       !       end if
       potisglue = potisglue+Eembi
       tabdensity(i) = eamglue(2,iti,k) + drk*( 2.0*eamglue(3,iti,k) + 3.0*drk*eamglue(4,iti,k) )
    end do loop2at1



#ifdef PARA

    if (nprocspace.gt.1) then
       call maj_tabdensity_ftm(tabdensity,atcf%imm,celcf%nato,atcf%num_at_glob,psc,atcf%im)
    end if

    !    write(3000+i,*)it
    !    do i=1,im
    !       write(6,*)i,num_at_glob(i),tabdensity(i)
    !    end do
#endif

    !    tabdensity=0
    !boucle des forces

    loop3at1: do i=1,atcf%im
       if (typ_and_pot(atcf%ityp(i),ipotentiel).eqv..false.)cycle

       ! --- Calcul de la densite sur i ---    

       koo = atcf%ielat(i)                          ! Numero de la cellule
       iti = atcf%ityp(i)
       ! pour chaque cel. voisine
       loop2cel:   do i1 = 0, celcf%ncelvois(koo)
          ko1 = celcf%ncel(koo,i1)
          ! pour chaque atome ds la cel. voisine
          loop2at2: do i2 = 1, celcf%nato(ko1)
             j = celcf%atincel(i2,ko1)
             if (typ_pot_pair(ipo(atcf%ityp(i),atcf%ityp(j))).ne.ipotentiel) cycle
             itj=atcf%ityp(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !CRC             if(i.eq.j) cycle
#ifdef PARA

             ! Methode pour ne prendre qu'une seule fois en compte
             ! le couple i,j en paralle :
             ! - i est necesairement local (boucle i<=im)
             ! - si j est local on ne retient que le couple i<j
             ! - si j n'est pas local, le couple n'est par definition
             !   pris qu'une fois puisque i est local
             if (nprocspace.gt.1) then

                if (j.le.atcf%im) then

                   ! les deux atomes sont locaux
                   if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme deja calcule
                else
                   ! j n'est pas local, on ne fait le calcul normal           
                endif
             else
                if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme dÃ£Â©ja calculÃ£Â©
             end if

#else
             if (atcf%num_at_glob(i).ge.atcf%num_at_glob(j)) cycle !terme dÃ£Â©ja calculÃ£Â©
#endif

             call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp,indcv=i1, lperiod=boxcf%lperiod,rum=rue,linter=linter,dist=r)
             if(.not.linter) cycle

             if (r.eq.zero) & 
                  write(*,*) '2. WARNING IN calfoeamcell TWO ATOMS VERY CLOSE i ,j , dist(angst)', i ,j , r*angst

             k=Int(r*inv_ktor)
             drk=r-k*ktor
             gradij(1:3) = dxp(1:3)/r
             Femb = ( eamrho(2,itj,k) + drk*( 2.0*eamrho(3,itj,k) + 3.0*drk*eamrho(4,itj,k) ) )*tabdensity(i) &
                  + ( eamrho(2,iti,k) + drk*( 2.0*eamrho(3,iti,k) + 3.0*drk*eamrho(4,iti,k) ) )*tabdensity(j)
             atcf%fp(1:3,i) = atcf%fp(1:3,i) - Femb*gradij(1:3)
             atcf%fp(1:3,j) = atcf%fp(1:3,j) + Femb*gradij(1:3)

             if (test_sigma) then                   
                if (atcf%num_at_glob(i).lt.atcf%num_at_glob(j)) then
                   sigem(1:3,1) = sigem(1:3,1) - Femb*gradij(1:3)*dxp(1)/boxcf%volu
                   sigem(1:3,2) = sigem(1:3,2) - Femb*gradij(1:3)*dxp(2)/boxcf%volu
                   sigem(1:3,3) = sigem(1:3,3) - Femb*gradij(1:3)*dxp(3)/boxcf%volu
!!$                   if (lTPcel.EQV..true.) then
!!$                      sigc(1:3,1,koo) =sigc(1:3,1,koo) - 0.5*Femb*gradij(1:3)*dxp(1)*celcf%noxyz/boxcf%volu
!!$                      sigc(1:3,2,koo) =sigc(1:3,2,koo) - 0.5*Femb*gradij(1:3)*dxp(2)*celcf%noxyz/boxcf%volu
!!$                      sigc(1:3,3,koo) =sigc(1:3,3,koo) - 0.5*Femb*gradij(1:3)*dxp(3)*celcf%noxyz/boxcf%volu
!!$                      sigc(1:3,1,ko1) =sigc(1:3,1,ko1) - 0.5*Femb*gradij(1:3)*dxp(1)*celcf%noxyz/boxcf%volu
!!$                      sigc(1:3,2,ko1) =sigc(1:3,2,ko1) - 0.5*Femb*gradij(1:3)*dxp(2)*celcf%noxyz/boxcf%volu
!!$                      sigc(1:3,3,ko1) =sigc(1:3,3,ko1) - 0.5*Femb*gradij(1:3)*dxp(3)*celcf%noxyz/boxcf%volu
!!$                   end if
!!$
                endif
             end if

          end do loop2at2
       end do loop2cel
    end do loop3at1



#ifdef PARA

    if (nprocspace.gt.1) then
       call comm_space%sum(potisglue)
       if (test_sigma) then 
          !          call comm_space%sum(sig)
          !          call comm_space%sum(sig2p)
          call comm_space%sum(sigem)

          !          if (associated(sigc)) then
          !             call comm_space%sum(sigc)
          !       endif
       endif
    end if
#endif
    !    if (test_sigma)sig=sig+sig2p+sigem
!!$    if (rang==0)    write(6,*)
!!$    if (rang==0)    write(6,*)sig2p
!!$    if (rang==0)    write(6,*)
!!$    if (rang==0)    write(6,*)sigem
    !   potiseam=potisglue+potisrep

    return
  end SUBROUTINE calfoglue_arps_2




  subroutine setfree(atcf,celcf)


    
    type(atom_config_arps)::atcf
    type(cell_config_arps):: celcf
    real(double)::vn
    !    real(double)::m,vmini,vmaxi,par1,par2,d1,d2,apar,bpar,pmin,pmax&
    !         &,xs,vn
    integer::i,iti,mov,nmov(0:2),imc
    real(double),save::nmovm(0:2)=0
    integer,save::icall=0
    integer::koxc(3),koo,ic,natvv
    logical::lok(3)
    nmov=0
#ifdef PARA
    if ((lspaceNDM.eqv..true.).and.(nprocspace.gt.1)) then
       imc=atcf%imf
    else
       imc=atcf%im
    end if
#else
    imc=atcf%im
#endif
    icall=icall+1
    natvv=0

    do i=1,imc
       lok=.false.
       koo=atcf%ielat(i)
       koxc(:)=celcf%koxyz(koo)
       !       write(6,*)koxc,noxyzkmin,noxyzkmax
       if (lpartarps) then 
          do ic=1,3
             if ((koxc(ic).ge.noxyzkmin(ic)).and.(koxc(ic).le.noxyzkmax(ic))) lok(ic)=.true.
          end do
          !      write(6,*)lok
       end if
          
       if(lok(1).and.lok(2).and.lok(3)) then
          natvv=natvv+1
          atcf%mov(i)=2
          if (i.le.atcf%im)          nmov(2)=nmov(2)+1
          mov=atcf%mov(i)
       else
          
          vn=norm2(atcf%vp(:,i))
          iti=atcf%ityp(i)
          if (vn.le.vmin(iti))  then
             mov=0
             if (i.le.atcf%im)   nmov(0)=nmov(0)+1
          else if ((vn.lt.vmax(iti)).and.(vn.gt.vmin(iti))) then
             mov=1
             if (i.le.atcf%im)          nmov(1)=nmov(1)+1
          else
             mov=2
             if (i.le.atcf%im)          nmov(2)=nmov(2)+1
          end if
          atcf%mov(i)=mov
       end if
       select case (mov)
       case(0)
          atcf%lgul(i)=.false.
       case(1,2)
          atcf%lgul(i)=.true.
       end select

    end do
#ifdef PARA
    call comm_space%sum(nmov)
    call comm_space%sum(natvv)
#endif
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             atcf%mov=2

    nmovm=(nmovm*(icall-1)+nmov)/icall
    if (itetemp>0) then
       if ((mod(iteration,itetemp)==0).and.(rang==0)) then
          write(6,'(A,I6,E15.6, 3I9)')'NMOV ',iteration, timel, nmov
          write(6,'(A,I6,E15.6, 3G15.5)')'NMOVM ',iteration, timel, nmovm
          if (lpartarps) write(6,*)'PARTARPS natvv nattot',natvv,atcf%im_glob
       end if
    end if


  end subroutine setfree


  subroutine initarps(atcf,celcf)
    type(atom_config_arps)::atcf
    type(cell_config_arps)::celcf
    integer::iti,ko,koxyz(3)
    real(double)::m,pmin,pmax,par1,par2,d1,d2
    !    allocate (atcf%fpr(3,atcf%imm))
    !    allocate (atcf%mov(atcf%imm))
    atcf%mov(atcf%imm)=2
    allocate (vmin(ntyp));    allocate(vmax(ntyp));allocate (xs(ntyp))
    allocate (aspl(ntyp));   allocate(bspl(ntyp));    allocate (s0spl(ntyp))
    do iti=1,ntyp
       vmin(iti)=sqrt(2*kmin/cm(iti))
       vmax(iti)=sqrt(2*kmax/cm(iti))
    end do
    if (ipotentiel.ge.10) then
       allocate (atcf%rho(atcf%imm))
       allocate (tabdensity(atcf%imm))
       allocate (atcf%fpg(3,atcf%imm))
    end if
    allocate(celcf%nmov(0:2,celcf%noxyz))

    do iti=1,ntyp

       m=cm(iti)
       !      vmini=vmin(iti);vmaxi=vmax(iti)
       pmin=vmin(iti)*m
       pmax=vmax(iti)*m
       xs(iti)=pmax-pmin


       par1=vmax(iti)
       par2=1/m
       d1=par1/(xs(iti)**2)
       d2=par2/xs(iti)
       aspl(iti)=d2/xs(iti)-2*d1/xs(iti)
       bspl(iti)=d1-aspl(iti)*xs(iti)
       S0spl(iti)=pmax**2/(2*m)-ppol(xs(iti),aspl(iti),bspl(iti))
!!$       do i=1,10000
!!$          vt=2*vmax(iti)*float(i)/10000
!!$          xpar=vt*m-vmin(iti)*m
!!$          if (xpar.le.0)  then
!!$             fact=0.!xp unchanged
!!$          else if ((xpar.lt.xs(iti)).and.(xpar.gt.0)) then
!!$             fact=pol(xpar,aspl(iti),bspl(iti))
!!$          else
!!$             fact = vt
!!$          end if
!!$          write(116,'(I7,3E17.6)')i,vt,xpar, fact
!!$       end do
!!$       stop
    end do
  end subroutine initarps


  subroutine xpupdate(atcf)
    type(atom_config_arps)::atcf
    integer::i,iti
    real(double)::vn,m,xpar,fact
    !    write(6,*)'in xpupdate'
    do i=1,atcf%im
          iti=atcf%ityp(i)
          m=cm(iti)
          vn=norm2(atcf%vp(:,i))
          xpar=vn*m-vmin(iti)*m

       if (atcf%mov(i)==2) then
          atcf%xp(:,i) = atcf%xp(:,i) + tstep*atcf%vp(:,i)
          fact=1
       else
          
!          write(6,*)'xpar2',xpar,xs(iti)
          if (xpar.le.0)  then
             !xp unchanged
                fact=0
          else if ((xpar.lt.xs(iti)).and.(xpar.gt.0)) then
             atcf%xp(:,i) = atcf%xp(:,i) + tstep*pol(xpar,aspl(iti),bspl(iti))*atcf%vp(:,i)/vn
                     fact=pol(xpar,aspl(iti),bspl(iti))/vn
          else
             atcf%xp(:,i) = atcf%xp(:,i) + tstep*atcf%vp(:,i)
!             fact=1
             write(6,*)'POOOOO'
             call arret_ndm
          end if
       end if
!       write(6,*)'xpar2',xpar,xs(iti),fact
       !       write(6,*)vn,fact

    end do
  end subroutine xpupdate


  subroutine caltabtarps( celcf,atcf,lperiod,boxcf,psc)
    type(cell_config_arps)::celcf
    type(atom_config_arps)::atcf
    logical::lperiod
    class(box_config)::boxcf
    type(para_space_config)::psc

    integer::ko,i,i1,imov
    celcf%nmov=0

    do ko=1,celcf%noxyz
#ifdef PARA
       if ( celcf%proc_cell(ko).ne.myidsp ) cycle
#endif
       do i1=1,celcf%nato(ko)
          i=celcf%atincel(i1,ko)
          imov=atcf%mov(i)
          celcf%nmov(imov,ko)=celcf%nmov(imov,ko)+1
       end do
    end do
#ifdef PARA
    call comm_space%sum(celcf%nmov)

#endif


  end subroutine caltabtarps

end module arps_mod
