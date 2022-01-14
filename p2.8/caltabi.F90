!****************************************************************
module caltabi_mod
  USE arret_ndm_mod,only:arret_ndm
  USE gen_com_m, ONLY:it,lconstrtot,ldemitab,lperiod,rang
  use atomconfig,only: atom_config
  USE cellconfig,only:cell_config
  use boxconfig,only:box_config
  use vect_dist_mod,only:vect_dist
  implicit none
contains



  ! *****************************************************************
  subroutine caltabi(atvois,celvois,boxndm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double


    USE var_pot, ONLY:ipotentiel,npair,ipo
    !           version du 4 juin 2010, 14h38 - last chaged by MCM
    ! *****************************************************************

    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !----------------------------------------------1-
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    class(atom_config), intent(inout)::atvois
    type(cell_config), intent(in)::celvois
    type(box_config),intent(in)::boxndm
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: iw, iwph, i, ip, j, maxvoi, nvij,iwo
    integer :: itj,ll
    REAL(double) :: r2
    real(double), dimension(1:npair) :: rvois2,rvois

    real(double), dimension(3) :: xpi, dx, ds
    real(double),dimension(3,3)::at,bg
    integer :: iti, & !type de i
         koo, & !cel de i
         ko1, & !cel voisine de i
         i1,i2,itemp

    logical::linter
    
    at=boxndm%at ; bg=boxndm%bg
    !
    !-----------------------------------------------
    ! --------------------------
    !   OUVERTURE BOUCLE SUR I
    ! --------------------------

    if(celvois%icaltabt.ne.atvois%icaltabt) then
       write (6,*)'incoherence dans icaltabt caltabi'
       call arret_ndm
    end if

    iw = 0
    iwph = 0

    if (ipotentiel==12) then
       rvois2(1)=(7.0d-8)**2
       rvois2(2)=(3.5d-8)**2
       rvois2(3)=(2.8d-8)**2
       rvois(1)=(7.0d-8)
       rvois(2)=(3.5d-8)
       rvois(3)=(2.8d-8)
    else
       rvois(:)=atvois%rvois
    end if

    nvij=0


    !write(*,*) 'caltabi_inside  ', rvois, rvois2
    !*************construction par double boucle ****************
    if(lconstrtot) then  !construction par double boucle


       do i = 1, atvois%im
          iwo=iw
          iti=atvois%ityp(i)

          if(ldemitab)then
             ip = i+1
          else
             ip=1
          end if

          do j = ip, atvois%im
             if(i.eq.j) cycle
             itj=atvois%ityp(j)
             ll=ipo(iti,itj)
             
             call vect_dist(atvois,celvois,boxndm,i,j,lperiod=lperiod,linter=linter,rum=rvois(ll))
             if (.not.linter)cycle             
             iw = iw+1
             IF (iw.GT.atvois%nVois) THEN
                WRITE(0,'(a,i0)') 'Indice iw du tableau de voisin plus grand&
                     & que le max, nVois = ', atvois%nVois
                !                   WRITE(0,'(a)') 'Augmentez le nombre moyen de voisins par&
                !                        & atome dans le fichier *.din'
                WRITE(0,'(a,i0)') 'truc étrange dans setcellconf'
                STOP '< Caltabi >'
             END IF

             atvois%indi(iw) = j
             !           indi2(iw) = j
          end do
          atvois%iwmax(i) = iw
          nvij=iw-iwo

       end do   ! atvois%im
       maxvoi = iw
       !*************construction par celulle ****************
    else
       !     write(*,*) 'THE fist passage .........'
       do i = 1, atvois%im
          iwo=iw
          koo = atvois%ielat(i)                          ! Numero de la cellule
          iti=atvois%ityp(i)
          !        write(6,*)'atome i',i,iti
          ! pour chaque cel. voisine
          do i1 = 0, celvois%ncelvois(koo)
             ko1 = celvois%ncel(koo,i1)
             !           write(6,*)'i1 ko1 ',i1,ko1
             if (ko1==0) cycle
             loop_j: do i2 = 1, celvois%nato(ko1)
                j = celvois%atincel(i2,ko1)
                !                                write(6,*)'j ',j

                if(ldemitab) then
                   if(j.le.i) cycle !terme deja calcule
                else
                   if(j.eq.i) then
                      iwph=iwph+1
                      !                   indi2(iwph) = j
                      cycle
                   end if
                end if

                itj=atvois%ityp(j)
                ll=ipo(iti,itj)
                call vect_dist(atvois,celvois,boxndm,i,j,indcv=i1,lperiod=lperiod,linter=linter,rum=rvois(ll))

                if (.not.linter)cycle             
                iw = iw+1
                iwph = iwph+1
                !              write(6,*)i,koo,ko1,j,iw, at,bg,sqrt(r2)
                atvois%indi(iw) = j
                !              indi2(iwph) = j
             end do loop_j !i2
          end do !ncelvois
          atvois%iwmax(i) = iw
          nvij=iw-iwo
          !                 write(6,*)'NVIJ',i,nvij,iw
       end do ! fin i
       maxvoi=iw

    endif ! lconstrtot

    if ((rang==0).and.(it.le.100)) then
       !           write(6,*)'IT ',it,'  VOISINS ',maxvoi,' par atome ',float(maxvoi)/float(atvois%im)
    endif
    return
  end subroutine caltabi
end module caltabi_mod
