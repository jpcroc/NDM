module specialinit_mod
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e
  USE cellconfig, only:cell_config,caltabtC
  use boxconfig,only: box_config
  USE arret_ndm_mod,only: arret_ndm
  USE initcasca_mod,only: initcasca,itko
  use newunit_mod,only:newunit
  USE T_kind_param_m, ONLY:  double
  
  use Tpara,only:para_space_config

  USE gen_com_m, ONLY:lcasca,xko,xx0,yko,yy0,zko,zz0,eko,iko,rang,bk,lperiod 


 use vect_dist_mod,only:distat
  USE calctemp_mod,only: calctemp
  implicit none
    logical::lpressinit,lheatinit
    integer::iusp
    integer::iheatinit,ipressinit
    real(double)::tempheatinit,Eheatinit,Rheatinit,Rpressinit,deltapressinit,centre(3)

contains
  ! **************************************************************
  subroutine specialinit(atdml,boxndm,celndm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
#ifdef PARA
    USE Tpara,only:COMM_space,nprocspace
#else
    use Tpara,only:nprocspace

#endif

    ! **************************************************************

    implicit none
    class(atom_config)::atdml
    type(cell_config)::celndm
    class(box_config)::boxndm
    type(para_space_config)::psc

    
    namelist /spinit/lpressinit,Lheatinit,eko,xko,yko,zko,xx0,yy0,zz0,iheatinit,tempheatinit,Eheatinit,Rheatinit,&
         &ipressinit,Rpressinit,deltapressinit,iko,itko

    iko=-1
    itko=-1
    xko=1.;    yko=1.;    zko=1.
    xx0=0.25;yy0=0.25;zz0=0.25
    iheatinit=1 ! 1 = from a sphere at 0.5 0.5 0.5 radius Rheatinit; 2 from a slice  at x=0.5 +:- Rheatinit
    Rheatinit=10.0 ! 
    tempheatinit=-1
    eheatinit=-1.
    centre(1:3)=0.5
    ipressinit =1
    Rpressinit=5.0
    deltapressinit=0.3
    eko=-1
    call newunit(iusp)
    open(unit=iusp,file='specialinit.in')
    read(iusp,nml=spinit)
    close(iusp)
    rheatinit=rheatinit*1d-8
    rpressinit=rpressinit*1d-8
    deltapressinit=deltapressinit*1d-8
    if (rang==0) write(6,*)'lpressinit,lheatinit',lpressinit,lheatinit
    if (lcasca) then
       select type(atdml)
       class is (atom_config_e)
          if (rang==0) write (6, *) '----CASCADE-----'
          if (rang==0) write (6, *) 'projectile=', iko, ' energie=', eko
          if (rang==0) write (6, *) 'direction=', xko, yko, zko
          if (rang==0) write (6, *) 'position de depart : ', xx0, yy0, zz0
          
          
          
          call initcasca(atdml,celndm,boxndm)
       class default
          write(6,*) 'lcasca and not atom_donfig_e ?'
          call arret_ndm
       end select
       
114    format(a3,1x,3(f10.4,1x),i5)
    end if

    if (lpressinit) call initsppress(atdml,boxndm,celndm)
    if (lheatinit) call initspheat(atdml,boxndm,celndm)
    

    return
  end subroutine specialinit

  subroutine initspheat(atdml,boxndm,celndm)
    class(atom_config)::atdml
    type(cell_config)::celndm
    class(box_config)::boxndm
    logical ::lgs(atdml%imm)
    integer::i,imsph
    type(atom_config_d)::atsph
    type(cell_config):: celsph
    real(double)::dist
    logical::linsph
    real(double)::tempfin,tempdec,tempsph,tempsph2
    real(double)::xproj,norma,kinesph

    lgs=atdml%lgul
    atdml%lgul=.false.
    celsph=celndm
    select case(iheatinit)
    case(1)
       do i=1,atdml%im

          call distat(atdml%xp(:,i),box=boxndm,dist=dist,x0red=centre,linter=linsph,rum=rheatinit)
          atdml%lgul(i)=linsph
       end do
    case(2)
       norma=norm2(boxndm%at(:,1))
       
       do i=1,atdml%im
          xproj=(dot_product(atdml%xp(:,i),boxndm%at(:,1))/norma)-norma/2.
          if(abs(xproj).le.rheatinit)atdml%lgul(i)=.true.
          
       end do
    end select
    

    call atdml%fab(atsph,lback=.true.)
    imsph=atsph%im_glob
!    if (rang==0)    write(6,*)'atsph',atsph%im_glob
!    call atsph%print
!    stop
    call caltabtC(celsph,atsph,lperiod,boxndm,lchktrav=.false.)
    call calctemp(tempsph,kinesph,atsph,celsph)
!    if (rang==0)write(6,*)'tempsph',tempsph
    if(tempheatinit.ge.0) then
       tempfin=tempheatinit
    else
       tempfin=tempsph+eheatinit*2/(3*imsph*bk)
    end if
    do i=1,atsph%im
       atsph%vp(:,i)=atsph%vp(:,i)*sqrt(tempfin/tempsph)
    end do
    call calctemp(tempsph2,kinesph,atsph,celsph)
    call atsph%backto(atdml)
    atdml%lgul=lgs
    if (rang==0) then
       write(6,*)'heating of ', imsph, 'atoms inside Rheatinit',rheatinit*1d8
       write(6,*)'from ', tempsph, ' to ', tempsph2
    end if
  end subroutine initspheat


  
  subroutine initsppress(atdml,boxndm,celndm)
    class(atom_config)::atdml
    type(cell_config)::celndm
    class(box_config)::boxndm
    logical ::lgs(atdml%imm)
    integer::i,impr,ic
!    type(atom_config_d)::atpr
!    type(cell_config):: celpr
!    real(double)::dist
    logical::linpr
!    real(double)::tempfin,tempdec,temppr,temppr2
    real(double)::xproj,norma,dx(3),x0(3),dxn

!    lgs=atdml%lgul
!    atdml%lgul=.false.
!    celpr=celndm
    x0=0
    do i=1,3
       do ic=1,3
          x0(i)=x0(i)+centre(ic)*boxndm%at(i,ic)
       end do
    end do
    impr=0
    select case(iheatinit)
    case(1)
       do i=1,atdml%im

          call distat(atdml%xp(:,i),box=boxndm,x0red=centre,linter=linpr,rum=rpressinit,dist=dxn)
          
          !         atdml%lgul(i)=linpr
          if (linpr) then
             impr=impr+1
             dx(:)=atdml%xp(:,i)-x0(:)

             atdml%xp(:,i)=x0(:)+dx(:)+deltapressinit*(dx(:)/Rpressinit)

          end if
          
       end do
    case(2)
       norma=norm2(boxndm%at(:,1))
       
       do i=1,atdml%im
          xproj=(dot_product(atdml%xp(:,i),boxndm%at(:,1))/norma)-norma/2.
!          write(6,*)i,xproj,rheatinit
          if(abs(xproj).le.rheatinit) then !atdml%lgul(i)=.true.
             impr=impr+1
             atdml%xp(1,i)=atdml%xp(1,i)+abs(xproj/rheatinit)*sign(deltapressinit,xproj)
          end if
       end do
    end select
    

!    call atdml%fab(atpr,lback=.true.)


    
!    atdml%lgul=lgs

    
  end subroutine initsppress
  
end module specialinit_mod
