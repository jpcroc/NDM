module boxconfig
  USE T_kind_param_m
  use recips_mod,only:recips,calcvol
  use atomconfig,only:atom_config,atom_config_d,atom_config_e
  use Tpara,only:mpi_communicator,endmpi
  implicit none
  type box_config
     real(double):: at(3,3)
     real(double):: bg(3,3)
     real(double):: zl(3),zls2(3),nzl(3),volu,normat(3),normbg(3)
     integer(long)::icaltabt
     logical::lperiod
     integer::ipbc(3) ! conditions périodiques sur les plan b-c,a-c,a-b
   contains
     procedure, pass::print=>boxprint
     procedure, pass::master2slave=>boxmaster2slave
  end type box_config
contains

  subroutine boxmaster2slave(box,rgem,mpic)

    type(mpi_communicator),intent(in)::mpic
    class(box_config)::box
    integer,intent(in)::rgem
    real(double)::atl(3,3)
    atl=box%at(:,:)
    call mpic%bcast(rgem,atl)
    if (mpic%rank.ne.rgem) then 
       call updatebox(box,atl)
    end if
    return
  end subroutine boxmaster2slave
   
  
  subroutine initbox(boxnew,at,ipbc,zl)
    type(box_config),intent(inout)::boxnew
    real(double),intent(in),optional::at(3,3)
    real(double),optional,intent(in)::zl(3)
    integer,intent(in) ::ipbc(3)
    integer::i,ic

    call updatebox(boxnew,at,zl,check=0)
    boxnew%ipbc(1:3)=ipbc(1:3)
    return
  end subroutine initbox

  subroutine updatebox(boxnew,at,zl,check)
    type(box_config),intent(inout)::boxnew
    real(double),intent(in),optional::at(3,3)
    real(double),optional,intent(in)::zl(3)
    integer,optional::check
    real(double)::nbg
    integer::i,ic,chk=0
    if (present(check))chk=check
    if (chk==1) then 
       nbg=boxnew%bg(1,1)**2+boxnew%bg(1,2)**2+boxnew%bg(1,3)**2
       if (nbg==0) then
          write(6,*) 'this is not an update as bg=0 stop'
          stop
       end if
    end if
    if((present(zl).eqv..false.).and.(present(at).eqv..false.)) then
       write(6,*)'box init at ET zl indéfinis : STOP'
       stop
    end if
    if(present(zl).and.(present(at))) then
       write(6,*)'box init at ET zl définis : STOP'
       stop
    end if
    if (present(at)) then
       boxnew%at=at
       !                    write(6,*)at
    else
       boxnew%at=0
       do ic = 1, 3
          boxnew%at(ic,ic)=boxnew%zl(ic)
       end do
    end if

    call recips (at(1,1), at(1,2), at(1,3), boxnew%bg(1,1), boxnew%bg(1,2), boxnew%bg(1,3))
    do ic = 1, 3
       boxnew%normat(ic) = 0
       boxnew%normat(ic) = boxnew%normat(ic)+sum(at(:,ic)**2)
       boxnew%normat(ic) = sqrt(boxnew%normat(ic))
       boxnew%zl(ic) = boxnew%normat(ic)
       boxnew%normbg(ic)=sqrt(sum(boxnew%bg(:,ic)**2))
        boxnew%nzl(ic)=1.0/ boxnew%normbg(ic)
    end do
    boxnew%zls2 = boxnew%zl/2.0
    boxnew%volu=calcvol(at(1:3,1),at(1:3,2),at(1:3,3))
    return
  end subroutine updatebox


  subroutine ndm2boxconfig(at,bg,zl,zls2,nzl,volu,normat,boxnew)
    real(double):: at(3,3)
    real(double):: bg(3,3)
    real(double):: zl(3),zls2(3),nzl(3),volu,normat(3)
    integer(long)::icaltabt 
    type(box_config)::boxnew
    boxnew%at=at
    boxnew%bg=bg
    boxnew%zl=zl
    boxnew%zls2=zls2
    boxnew%nzl=nzl
    boxnew%volu=volu
    boxnew%normat=normat
    boxnew%icaltabt=0
    return
  end subroutine ndm2boxconfig
  subroutine boxconfig2ndm(at,bg,zl,zls2,nzl,volu,normat,boxnew)
    real(double):: at(3,3)
    real(double):: bg(3,3)
    real(double):: zl(3),zls2(3),nzl(3),volu,normat(3)
    integer(long)::icaltabt 
    type(box_config)::boxnew
    at=boxnew%at
    bg=boxnew%bg
    zl=boxnew%zl
    zls2=boxnew%zls2
    nzl=boxnew%nzl
    volu=boxnew%volu
    normat=boxnew%normat
    return
  end subroutine boxconfig2ndm

  subroutine boxprint(boxprt,unit)
    class(box_config)::boxprt
    integer,optional::unit
    integer::unitw
    unitw=6
    if (present(unit))unitw=unit
     write(unitw,*)'boxprt at',boxprt%at(:,:)
     write(unitw,*)'boxprt bg',boxprt%bg(:,:)
     write(unitw,*)'boxprt volu',boxprt%volu
     write(unitw,*)'boxprt icaltabt',boxprt%icaltabt
   end subroutine boxprint

   subroutine periodbox(box,atcf)
     !-----------------------------------------------
     !   M o d u l e s
     !-----------------------------------------------
     USE cryst_to_cart_mod,only: cryst_to_cart
     USE T_kind_param_m, ONLY:  double
     USE gen_com_m, ONLY:lperiod,zero

    USE gen_com_m, ONLY:low_limit,zero

     ! *****************************************************************
     ! Cette routine applique les conditions périodiques par décalage
     !  ou +/-1 (triclinique).
     !Le décalage est fait sur xp xp ET ax. Ce décalage sur ax 
     !permet de mesurer correctement le déplacements des atomes
     ! depuis leur position de départ

     implicit none
     !-----------------------------------------------
     !   G l o b a l   P a r a m e t e r s
     !-----------------------------------------------
     !-----------------------------------------------
     !   D u m m y   A r g u m e n t s
     !-----------------------------------------------
     type(box_config),intent(in)::box
     class(atom_config)::atcf

     !-----------------------------------------------
     !   L o c a l   P a r a m e t e r s
     !-----------------------------------------------
     !-----------------------------------------------
     !   L o c a l   V a r i a b l e s
     !-----------------------------------------------
     integer :: i, ic,icp!,nbing,ibing(100)
     real(double)::dz,trav,ecav,ecap
     real(double):: cpp,xpici,cppzl,ctest
     !      integer,save  :: iperiod
     !  if (rang==0) write(6,*)'PARA-T entree period'


     !      iperiod=iperiod+1

     !      write(*,*) 'PBC PBC PBC capitala tarii e ....             period',iperiod


     call cryst_to_cart (atcf%imm, atcf%xp,  box%bg,  -1) !cart vers cryst
!        call atcf%print (unit=10)
     select type (atcf)
     type is (atom_config_d)
        call cryst_to_cart (atcf%imm, atcf%xpp, box%bg,  -1)
     type is (atom_config_e)

        call cryst_to_cart (atcf%imm, atcf%xpp, box%bg,  -1)
        if(atcf%lax)      call cryst_to_cart (atcf%imm, atcf%ax,  box%bg,  -1)
     end select
!     nbing=0
     loopdir:do ic=1,3
        select case(box%ipbc(ic))
        case(1) ! periodic boundary conditions if lperiod otherwise positions can become <0 or >1
           if(.not.lperiod) cycle loopdir
           do i=1,atcf%imm
              xpici=atcf%xp(ic,i)
              if ( (xpici < 0.d0 ).OR.( xpici >= 1.d0 ) ) then
                 if ( (xpici > -low_limit).and.(xpici<0.d0) ) then
                    atcf%xp(ic,i)=zero
                 else
                    cpp  = Dble(Floor(atcf%xp(ic,i)))
                    select type (atcf)
                    type is (atom_config_d)
                       atcf%xpp(ic,i) = atcf%xpp(ic,i) - cpp
                    type is (atom_config_e)
                       atcf%xpp(ic,i) = atcf%xpp(ic,i) - cpp
                       if(atcf%lax)   atcf%ax (ic,i) = atcf%ax(ic,i)  - cpp
                    end select
                    atcf%xp (ic,i) = xpici     - cpp
                 end if
              end if
           end do
        case(2) !wall conditions
           select type (atcf)
           class is (atom_config_d)
              call cryst_to_cart (atcf%imm, atcf%vp,  box%bg,  -1) !cart vers cryst
           end select
           do i=1,atcf%im
              xpici=atcf%xp(ic,i)
              if  (xpici < 0.d0 ) then
                 atcf%xp(ic,i)=-xpici
                 select type (atcf)
                 class is (atom_config_d)
                    atcf%vp(ic,i)=-atcf%vp(ic,i)
                 end select
              else if (xpici>1) then
                 atcf%xp(ic,i)=2-xpici
                 select type (atcf)
                 class is (atom_config_d)
                    atcf%vp(ic,i)=-atcf%vp(ic,i)
                 end select
              end if
           end do
           select type (atcf)
           class is (atom_config_d)
              call cryst_to_cart (atcf%imm, atcf%vp,  box%at,  1) !cryst vers cart
           end select
        end select
           
     end do loopdir
     call cryst_to_cart (atcf%imm, atcf%xp , box%at,  1)  !cryst vers cart
     select type (atcf)
     type is (atom_config_d)
        call cryst_to_cart (atcf%imm, atcf%xpp , box%at,  1)  !cryst vers cart
     type is (atom_config_e)
        call cryst_to_cart (atcf%imm, atcf%xpp , box%at,  1)  !cryst vers cart
        if(atcf%lax) call cryst_to_cart (atcf%imm, atcf%ax , box%at,  1)  !cryst vers cart
     end select
   end subroutine periodbox


 end module boxconfig


