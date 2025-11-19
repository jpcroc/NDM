module scalebox_mod
  USE arret_ndm_mod,only:arret_ndm
  USE gen_com_m, ONLY:dmtype,itetabvois,lprahman,nvat,pi,iteration,rang,lperiod,lspacendm
  USE calpo_ew_mod,only: calpo_ew
  USE recips_mod,only: recips ,calcvol
  USE caltabi_mod,only: caltabi
  USE atomconfig,only : atom_config_d
  USE cellconfig, only:cell_config,caltabtC
  USE boxconfig, only:box_config,periodbox,box_config_lpr
#ifdef PARA
  USE mod_para,only:maj_atomes_frt_ftm
#endif
  use Tpara,only:para_space_config,nprocspace

  implicit none
contains
  ! ******************************************************************
  subroutine scalebox(atpr,celndm,boxndm,psc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE var_pot, ONLY:auxe,iewald,rumax

    ! ******************************************************************

    implicit none
    type(para_space_config)::psc
    class(atom_config_d)::atpr
    type(cell_config):: celndm
    class(box_config)::boxndm
    integer::noxn,noyn,nozn



    
    call periodbox (boxndm,atpr)

#ifndef PARA
    noxn = int(boxndm%zl(1)/rumax)
    noyn = int(boxndm%zl(2)/rumax)
    nozn = int(boxndm%zl(3)/rumax)
    
    if ((noxn==2).or.(noyn==2).or.(nozn==2))then
       noxn=1 ;noyn=1; nozn=1
    end if
    
    if ((celndm%nox(1).ne.noxn).or.(celndm%nox(2).ne.noyn).or.(celndm%nox(3).ne.nozn).or.((dmtype.eq.9).and.(iteration==1)))then
       write(6,*)'CHGT NOX'
       celndm%nox(1)=noxn; celndm%nox(2)=noyn; celndm%nox(3)=nozn

       if (dmtype.ne.9) then
          if (rang==0) write (6, *) 'IT =',ITeration,'chgt nox noy noz  = '&
               , celndm%nox(1),celndm%nox(2), celndm%nox(3)
       end if
       celndm%celsize(1) = boxndm%zl(1)/float(celndm%nox(1))
       celndm%celsize(2) = boxndm%zl(2)/float(celndm%nox(2))
       celndm%celsize(3) = boxndm%zl(3)/float(celndm%nox(3))
       celndm%noxyz = celndm%nox(1)*celndm%nox(2)*celndm%nox(3)
       celndm%natperc= INT(atpr%im/celndm%noxyz)
       celndm%natperc=max(3*celndm%natperc,10)
       nvat=3*celndm%natperc

       if (dmtype.ne.9) then
          if (rang==0)       write(6,*) ' natperc ', celndm%natperc
       end if
!       write(6,*)'BOUFFON!'
!       call arret_ndm
      call celndm%init(boxndm,celndm%nox(1),celndm%nox(2),celndm%nox(3),celndm%natperc) !contient dealloc

   end if
#endif
   if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
      call caltabtC(celndm,atpr,lperiod,boxndm,psc=psc,lchktrav=.true.)
    else
      call caltabtC(celndm,atpr,lperiod,boxndm,lchktrav=.false.)
    end if
    

#ifdef PARA
    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.))   call maj_atomes_frt_ftm(atpr,celndm,boxndm,psc)
#else
   if (atpr%ltabvois.and.(dmtype==9).and.((iteration==1).or.(mod(iteration,itetabvois)==0))) then
      call caltabi(atpr%atom_config,celndm,boxndm)
   end if
#endif
       if ((iewald.gt.0).and.(iewald.ne.3)) call calpo_ew(boxndm,atpr%imm)
    
    return
  end subroutine scalebox
end module scalebox_mod
