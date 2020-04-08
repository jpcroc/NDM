module calfoew_mod
  USE epme_mod
  USE gen_com_m, ONLY:imd,imm,ltpcel,tabf3,sig,sigc,noxyz
  implicit none
contains

  ! ***************************************************************
  subroutine calfoew
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE var_pot, ONLY:alpha,iewald,nvecttot,ncoucx,ncoucy,ncoucz,q
    USE tab_imm_m
#ifdef PARA
    USE mod_para
#endif
    ! ewald reciproque
    ! **************************************************************

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
    integer :: Deb, Fin
    real(double), dimension(imm) :: scalar
    integer nb1,nb2,nb3,i,i1,iti
    real(double) :: potisewg, hbn2, &
         hbv(3),phu
    real(double) :: scacos,scasin
#ifdef PARA
    real(double) :: scacos_glob, scasin_glob
#endif
    real(double), dimension(3,3) :: sige

    !parallelisation de ewald classique
    integer,dimension(nvecttot) :: nb1v,nb2v,nb3v ! tableaux des vecteurs du RRec
    integer :: nv,debv,finv,ii
    real (double)::potis3p
    real (double), dimension (3,3) :: sigep
    !  real(double), dimension (3,imm) :: fpewp

    ! EWALD RECIPROQUE METHODE CLASSIQUE (iewald=1)

    potisewg = zero
    sige = 0.0
    select case (iewald)
    case (1)
       nv=0
       ! repartition des vecteurs du RRec.
       do nb1 = -ncoucx, ncoucx
          do nb2 = -ncoucy, ncoucy
             do nb3 = -ncoucz, ncoucz
                if (nb2==0.and.nb3==0.and.nb1==0) cycle
                nv=nv+1
                nb1v(nv)=nb1; nb2v(nv)=nb2; nb3v(nv)=nb3
             enddo
          enddo
       enddo
       if (nv.ne.nvecttot) stop

       debv=1
       finv=nvecttot

       !     fpewp(:,:)=0.0
       sigep(:,:)=0.0
       potis3p=0.0


       ! *** Somme sur les vecteurs NB1,NB2,NB3
       do nb1 = -ncoucx, ncoucx
          do nb2 = -ncoucy, ncoucy
             do nb3 = -ncoucz, ncoucz

                if (nb2==0.and.nb3==0.and.nb1==0) cycle

                do ii=1,size(hbv)
                   hbv(ii) = 2.d0*pi*(bg(1,ii)*nb1+bg(2,ii)*nb2+bg(3,ii)*nb3)
                enddo

                do ii=1,imd
                   scalar(ii)=xp(1,ii)*hbv(1)+xp(2,ii)*hbv(2)+xp(3,ii)*hbv(3)
                enddo

                ! Modification CCRT pour le compilo 6.0-4

                scacos = 0
                scasin = 0
                do ii=1,imd
                   scacos = scacos + cos(scalar(ii))*q(ityp(ii))
                   scasin = scasin + sin(scalar(ii))*q(ityp(ii))
                enddo

#ifdef PARA
                ! Reduction MPI en interne de la boucle. Prefere au stockage dans des tableaux
                ! (pour scalar et hbv il faudrait ajouter des dimensions ncoux/y/z)
                call MPI_ALLREDUCE(scacos,scacos_glob,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
                scacos = scacos_glob
                call MPI_ALLREDUCE(scasin,scasin_glob,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
                scasin = scasin_glob
#endif
                do i = 1, imd
                   iti = ityp(i)
                   phu = tabf3(iti,nb1,nb2,nb3)*(sin(scalar(i))*scacos-&
                        cos(scalar(i))*scasin)
                   fp(1,i) = fp(1,i)+phu*hbv(1)/(2.D0*pi)
                   fp(2,i) = fp(2,i)+phu*hbv(2)/(2.D0*pi)
                   fp(3,i) = fp(3,i)+phu*hbv(3)/(2.D0*pi)
                end do


                ! calcul de sig contrainte
                potisewg = tabv3(nb1,nb2,nb3)*(scacos**2+scasin**2)
                hbn2 = hbv(1)**2+hbv(2)**2+hbv(3)**2
                if (itesigma>0) then
                   if (mod(it,itesigma)==0) then
                      sige(1,1) = sige(1,1)+potisewg*hbv(1)*hbv(1)/hbn2*(hbn2/(4.0*&
                           &         alpha**2)+1)/volu
                      sige(2,2) = sige(2,2)+potisewg*hbv(2)*hbv(2)/hbn2*(hbn2/(4.0*&
                           &         alpha**2)+1)/volu
                      sige(3,3) = sige(3,3)+potisewg*hbv(3)*hbv(3)/hbn2*(hbn2/(4.0*&
                           &        alpha**2)+1)/volu
                   endif
                endif

                potis3 = potis3+potisewg
             end do
          end do
       end do

       ! --- Fin du calcul ---

       if (itesigma>0) then
          if (mod(it,itesigma)==0) then
             do i1 = 1, 3
                sig(i1,i1) = sig(i1,i1)+sige(i1,i1)
                if (lTPcel.EQV..true.) then
                   sigc(i1,i1,:noxyz) = sigc(i1,i1,:noxyz)+sige(i1,i1)
                end if
             end do
          endif
       endif



       ! FIN DU TERME EWALD dans l'espace reciproque



    case(2)  !Traitement par la methode PME


       ! Sequentiel

       Deb=1 !Test
       Fin=imd !Test

       call epme (Deb,Fin,sige)

       if (itesigma>0) then
          if (mod(it,itesigma)==0) then
             do i1 = 1, 3
                sig(i1,i1) = sig(i1,i1)+sige(i1,i1)
                if (lTPcel.EQV..true.) then
                   sigc(i1,i1,:noxyz) = sigc(i1,i1,:noxyz)+sige(i1,i1)
                end if
             end do
          endif
       endif

    end select


    return

  end subroutine calfoew



end module calfoew_mod
