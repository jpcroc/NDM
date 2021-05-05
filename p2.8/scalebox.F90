module scalebox_mod
  USE gen_com_m, ONLY:dmtype,itetabvois,lpr,nvat,pi,it,rang,lperiod

  USE recips_mod,only: recips ,calcvol
  USE caltabi_mod,only: caltabi
  USE atomconfig,only : atom_config_d,ndm2config, config2ndm
  USE cellconfig, only:cell_config,ndm2cellconfig,cellconfig2ndm,caltabtC
  USE boxconfig, only:box_config,ndm2boxconfig,periodbox

  implicit none
contains
  ! ******************************************************************
  subroutine scalebox(atpr,celndm,boxndm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE var_pot, ONLY:alpha,auxe,rumax,tabv3,tabf3,ncoucx,ncoucy,ncoucz,iewald,q

    !           (version du 09 juin 2000)
    ! ******************************************************************

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
    integer :: i, nb1, nb2, nb3, i1, l,noxn,noyn,nozn
    real(double) :: zlx, zly, zlz, ux, uy, uz,  pi2, fact, fact1&
         , fact2, hk2, ex, ex1, ex2
    type(atom_config_d)::atpr
    type(cell_config):: celndm
    type(box_config)::boxndm
    !real(double), external :: calcvol
    !-----------------------------------------------
    !
    !
    !
    !
    ! reciprocal space ewald summation



    !debug      write (*,*) 'sub scalebox',it,xp(1,1)
    ! -------------------------------------------------------------
    ! Rescaling des positions
    ! -------------------------------------------------------------
    if (lperiod)    call periodbox (boxndm,atpr)

    !debug       write (*,*) 'sub scalebox',it,xp(1,1)


    ! ---------------------------------------------------------------
    ! Recalcul des quantites dependantes de la dimension
    ! ---------------------------------------------------------------
    if (lpr) then
       boxndm%zl(1) = Sqrt( Sum(boxndm%at(1:3,1)**2 ) )
       boxndm%zl(2) = Sqrt( Sum(boxndm%at(1:3,2)**2 ) )
       boxndm%zl(3) = Sqrt( Sum(boxndm%at(1:3,3)**2 ) )
    endif

    boxndm%volu=calcvol(boxndm%at(1:3,1),boxndm%at(1:3,2),boxndm%at(1:3,3))
    boxndm%zls2(1:3) = 0.5d0*boxndm%zl(1:3)

    noxn = int(boxndm%zl(1)/rumax)
    noyn = int(boxndm%zl(2)/rumax)
    nozn = int(boxndm%zl(3)/rumax)

#ifdef ML
    if (noxn==0) noxn=1
    if (noyn==0) noyn=1
    if (nozn==0) nozn=1

    if (celndm%nox==0) celndm%nox=1
    if (noy==0) celndm%noy=1
    if (noz==0) celndm%noz=1
#endif

    if ((noxn==2).or.(noyn==2).or.(nozn==2))then
       noxn=1 ;noyn=1; nozn=1
    end if

    if ((celndm%nox.ne.noxn).or.(celndm%noy.ne.noyn).or.(celndm%noz.ne.nozn).or.((dmtype.eq.9).and.(it==1)))then
!       call Deallocatecel !fait dans %init
       celndm%nox=noxn; celndm%noy=noyn; celndm%noz=nozn

       if (dmtype.ne.9) then
          if (rang==0) write (6, *) 'IT =',IT,'chgt nox noy noz  = '&
               , celndm%nox,celndm%noy, celndm%noz
       end if

       celndm%celsize(1) = boxndm%zl(1)/float(celndm%nox)
       celndm%celsize(2) = boxndm%zl(2)/float(celndm%noy)
       celndm%celsize(3) = boxndm%zl(3)/float(celndm%noz)

       celndm%noxyz = celndm%nox*celndm%noy*celndm%noz

       !write(*,*) 'inside scalebox1', nox, noy, noz
       !write(*,*) 'inside scalebox2', noxn, noyn, nozn
       !write(*,*) 'inside scalebox3', noxyz, zl(1), rumax,  im
       celndm%natperc= INT(atpr%im/celndm%noxyz)

       celndm%natperc=max(3*celndm%natperc,10)
       nvat=10*celndm%natperc

       if (dmtype.ne.9) then
          if (rang==0)       write(6,*) ' natperc ', celndm%natperc
       end if
!       write(6,*)'BOUFFON!'
!       stop
      call celndm%init(celndm%nox,celndm%noy,celndm%noz,celndm%natperc) !contient dealloc

    end if
             call caltabtC(celndm,atpr,lperiod,boxndm)
             if (atpr%ltabvois.and.(dmtype==9).and.((it==1).or.(mod(it,itetabvois)==0))) then
                call caltabi(atpr%atom_config,celndm,boxndm)
             end if

    if (iewald>0) then

       ! --- Tableaux des troisiemes termes de la sommation d'Ewald ---
       auxe = 23.06134575D-20                  ! en erg.cm (charge electron^2/4*pi*permitivite vide)
       pi2 = pi*pi
       !volu = zl(1)*zl(2)*zl(3)
       boxndm%volu=calcvol(boxndm%at(1:3,1),boxndm%at(1:3,2),boxndm%at(1:3,3))
       fact = pi2/alpha**2
       fact1 = auxe/2./pi/boxndm%volu
       fact2 = auxe*2./boxndm%volu
       do nb1 = -ncoucx, ncoucx
          do nb2 = -ncoucy, ncoucy
             do nb3 = -ncoucz, ncoucz
                if (nb1==0.and.nb2==0.and.nb3==0) cycle
                hk2 = nb1*nb1/boxndm%zl(1)**2+nb2*nb2/boxndm%zl(2)**2+nb3*nb3/boxndm%zl(3)**2
                ex = exp((-hk2*fact))/hk2
                ex1 = ex*fact1
                ex2 = ex*fact2
                tabv3(nb1,nb2,nb3) = ex1
                tabf3(:,nb1,nb2,nb3) = ex2*q(:)
             end do
          end do
       end do

    endif



    return
  end subroutine scalebox
end module scalebox_mod
