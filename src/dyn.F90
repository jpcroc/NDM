module dyn_mod
  USE gen_com_m, ONLY:cunite,erg2ev,fnemd,iteration,itetconst,lcalcjq,leev,lnemd,lperiod,&
       &ltcon,text,timel,tstep,unite,usdh,bk
  use atomconfig,only:atom_config_d,atom_config_e
  implicit none
contains
  ! *************************************************************
  subroutine dyn (atdml)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    
    USE var_pot, ONLY:ntyp,cm
    USE jqmod

    implicit none

    class(atom_config_d)::atdml
    integer :: i,  ic,im
    real(double), dimension(ntyp) :: aux
    real(double) :: xprov
    real(double):: eatommoy
    im=atdml%im
    if(lEev) then
       unitE=erg2eV
       cunitE='  eV '
    else
       unitE=1.0
       cunitE=' erg '
    end if

    select type (atdml)
    class is (atom_config_e)
       if (atdml%lsigat) then
          atdml%eat(1:im)=atdml%eat(1:im)+&
               &0.5*cm(atdml%ityp(1:im))*(atdml%vp(1,1:im)**2+atdml%vp(2,1:im)**2+atdml%vp(3,1:im)**2)
          
          if (lnemd) then
             eatommoy=0.
             do i=1,im
          !        eatom(i)=eatom(i)+0.5*cm(ityp(i))*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
          eatommoy=eatommoy+atdml%eat(i)/float(im)
       end do
       do i=1,im

          atdml%fp(1,i)=atdml%fp(1,i)+(atdml%eat(i)-eatommoy)*Fnemd
       end do
       !     do i=1,im
       !        eatom(i)=eatom(i)-0.5*cm(ityp(i))*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
       !     end do
    end if
 end if
end select

    timel = timel+tstep
    aux(:ntyp) = tstep**2/cm(:ntyp)
    !      write(6,*)'aux ',aux


    !debug write(*,*) 'md_test1',  xp(1,1), xpp(1,1), vp(1,1) 
    do i = 1, im
       do ic = 1, 3
          xprov = (atdml%xp(ic,i)-atdml%xpp(ic,i))+atdml%xp(ic,i)+aux(atdml%ityp(i))*atdml%fp(ic,i)
          atdml%vp(ic,i) = (xprov-atdml%xpp(ic,i))*usdh
          atdml%xpp(ic,i) = atdml%xp(ic,i)
          atdml%xp(ic,i) = xprov
       end do
    end do
    !debug write(*,*) 'md_test2',  xp(1,1), xpp(1,1), vp(1,1)
    select type (atdml)
    class is (atom_config_e)
       if (atdml%lsigat) then
          if (lcalcjq) then
             eatommoy=0.
             jqp=jq ; jqk=0.0 !; expvect(:)=0.0
             do i=1,im
                !        eatom(i)=eatom(i)+0.5*cm(ityp(i))*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
                if(lnemd) eatommoy=eatommoy+atdml%eat(i)/float(atdml%im)
                expvect(:)=expvect(:)+atdml%eat(i)*atdml%xpp(:,i)
                do ic=1,3              
                   jqk(ic)=jqk(ic)+atdml%eat(i)*atdml%vp(ic,i)
                   !               jq(ic)=jq(ic)+eatom(i)*vp(ic,i)
                end do
             end do
             jq=jqp+jqk



       write(65,'(I8,3D16.8)')iteration-1,jq(1),jq(2),jq(3)
       !         write(66,'(I8,3D15.6)')it-1,expvect(1),expvect(2),expvect(3)
    end if
 end if
end select

    return
  end subroutine dyn

end module dyn_mod
