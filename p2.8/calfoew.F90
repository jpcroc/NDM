module calfoew_mod
  USE epme_mod,only: epme
  USE gen_com_m, ONLY:pi,potis3,zero,pi,potis2
  USE calfocommon
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only : cell_config
  use boxconfig,only: box_config
  USE var_pot, ONLY:alpha,iewald,nvecttot,ncoucx,ncoucy,ncoucz,q,tabv3,tabf3,auxe,ipo,zz
  USE recips_mod,only: calcvol
  implicit none
contains

    
  ! ***************************************************************
!  subroutine calfoew(im,imm,xp,fp,ityp,noxyz,at,bg,volu)
  subroutine calfozz(atcf)
#ifdef PARA
    use Tpara,only:nprocspace,comm_space
#endif
    
    class(atom_config),intent(in)::atcf

    integer::iti,l,i
    
    do i=1,atcf%im
       iti = atcf%ityp(i)
       ! --- Calcul du second potentiel de la somme d'Ewald ---
       l = ipo(iti,iti)
       potis2 = potis2-zz(l)*alpha/sqrt(pi)*23.06134575D-20
    end do
#ifdef PARA
    if (nprocspace.gt.1) then
       call comm_space%sum(potis2)
    end if
#endif       
  end subroutine calfozz
    
  subroutine calfoew(atcf,celcf,boxcf)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double


#ifdef PARA
    USE Tpara,only:COMM_space,nprocspace
#else
    USE Tpara,only:nprocspace
#endif
    ! ewald reciproque
    ! **************************************************************

    implicit none
    class(atom_config),intent(inout)::atcf
    type(cell_config),intent(in)::celcf
    type(box_config),intent(in)::boxcf
    integer :: Deb, Fin
    real(double), dimension(atcf%im) :: scalar
    integer nb1,nb2,nb3,i,i1,iti
    real(double) :: potisewg, hbn2, &
         hbv(3),phu
    real(double) :: scacos,scasin
    real(double), dimension(3,3) :: sige

    !parallelisation de ewald classique

    integer :: nv,debv,finv,ii,l
    real (double), dimension (3,3) :: sigep
    !  real(double), dimension (3,imm) :: fpewp

    ! EWALD RECIPROQUE METHODE CLASSIQUE (iewald=1)
    potisewg = zero
    sige = 0.0
    select case (iewald)
    case (1)

       debv=1
       finv=nvecttot

       !     fpewp(:,:)=0.0
       sigep(:,:)=0.0


       ! *** Somme sur les vecteurs NB1,NB2,NB3
       do nb1 = -ncoucx, ncoucx
          do nb2 = -ncoucy, ncoucy
             do nb3 = -ncoucz, ncoucz

                if (nb2==0.and.nb3==0.and.nb1==0) cycle

                do ii=1,size(hbv)
                   hbv(ii) = 2.d0*pi*(boxcf%bg(1,ii)*nb1+boxcf%bg(2,ii)*nb2+boxcf%bg(3,ii)*nb3)
                enddo

                do ii=1,atcf%im
                   scalar(ii)=atcf%xp(1,ii)*hbv(1)+atcf%xp(2,ii)*hbv(2)+atcf%xp(3,ii)*hbv(3)
                enddo

                ! Modification CCRT pour le compilo 6.0-4

                scacos = 0
                scasin = 0
                do ii=1,atcf%im
                   scacos = scacos + cos(scalar(ii))*q(atcf%ityp(ii))
                   scasin = scasin + sin(scalar(ii))*q(atcf%ityp(ii))
                enddo

#ifdef PARA
                ! Reduction MPI en interne de la boucle. Prefere au stockage dans des tableaux
                ! (pour scalar et hbv il faudrait ajouter des dimensions ncoux/y/z)
                if (nprocspace.gt.1) then
                   call comm_space%sum(scacos)
                   call comm_space%sum(scasin)
                end if
#endif

                do i = 1, atcf%im
                   iti = atcf%ityp(i)
                   phu = tabf3(iti,nb1,nb2,nb3)*(sin(scalar(i))*scacos-&
                        cos(scalar(i))*scasin)
                   atcf%fp(1,i) = atcf%fp(1,i)+phu*hbv(1)/(2.D0*pi)
!                   if (((nb1==0).or.(nb2==0).or.(nb3==0)).and.(i.lt.10))write(6,*)i,nb1,nb2,nb3,phu*hbv(1)/(2.D0*pi),phu,hbv(1)
                   atcf%fp(2,i) = atcf%fp(2,i)+phu*hbv(2)/(2.D0*pi)
                   atcf%fp(3,i) = atcf%fp(3,i)+phu*hbv(3)/(2.D0*pi)
                end do


                ! calcul de sig contrainte
                potisewg = tabv3(nb1,nb2,nb3)*(scacos**2+scasin**2)
                hbn2 = hbv(1)**2+hbv(2)**2+hbv(3)**2
                if (test_sigma) then
                   sige(1,1) = sige(1,1)+potisewg*hbv(1)*hbv(1)/hbn2*(hbn2/(4.0*&
                        &         alpha**2)+1)/boxcf%volu
                   sige(2,2) = sige(2,2)+potisewg*hbv(2)*hbv(2)/hbn2*(hbn2/(4.0*&
                        &         alpha**2)+1)/boxcf%volu
                   sige(3,3) = sige(3,3)+potisewg*hbv(3)*hbv(3)/hbn2*(hbn2/(4.0*&
                        &        alpha**2)+1)/boxcf%volu
                endif

                potis3 = potis3+potisewg
             end do
          end do
       end do

       ! --- Fin du calcul ---

       if (test_sigma) then
          do i1 = 1, 3
             sig(i1,i1) = sig(i1,i1)+sige(i1,i1)
             if (lTPcel.EQV..true.) then
                sigc(i1,i1,:celcf%noxyz) = sigc(i1,i1,:celcf%noxyz)+sige(i1,i1)
             end if
          end do
       endif



       ! FIN DU TERME EWALD dans l'espace reciproque



    case(2)  !Traitement par la methode PME


       ! Sequentiel

       Deb=1 !Test
       Fin=atcf%im !Test

       call epme (Deb,Fin,sige,atcf%im,atcf%xp,atcf%fp,atcf%ityp,boxcf%volu,boxcf%bg)

       if (test_sigma) then
          do i1 = 1, 3
             sig(i1,i1) = sig(i1,i1)+sige(i1,i1)
             if (lTPcel.EQV..true.) then
                sigc(i1,i1,:celcf%noxyz) = sigc(i1,i1,:celcf%noxyz)+sige(i1,i1)
             end if
          end do
       endif

    end select


    return

  end subroutine calfoew



end module calfoew_mod
