!****************************************************************
module caltabi_mod
  USE arret_ndm_mod,only:arret_ndm
  USE gen_com_m, only:uwrt,lwrt,lconstrtot,ldemitab,lperiod,rang
  use atomconfig,only: atom_config
  USE cellconfig,only:cell_config
  use boxconfig,only:box_config
  use vect_dist_mod,only:vect_dist
  use derived_types,only:system_state
  USE T_kind_param_m, ONLY:  double
  implicit none
contains
  subroutine caltabi(atvois,celvois,boxndm,lextr,lconstrtotR,cn2m)

    USE var_pot, ONLY:ipotentiel,npair,ipo
    implicit none
    class(atom_config), intent(inout)::atvois
    type(cell_config), intent(in)::celvois
    class(box_config),intent(in)::boxndm
    logical,optional::lextr,lconstrtotR
    logical::lextrait=.false.,lconstrtt

    type(system_state),optional::cn2m

    integer :: iw, iwph, i, ip, j, maxvoi, nvij,iwo,nvi
    integer :: itj,ll
    REAL(double) :: dij
    real(double), dimension(1:npair) :: rvois2,rvois

    real(double), dimension(3) :: VJI
    real(double),dimension(3,3)::at,bg
    integer :: iti, & !type de i
         koo, & !cel de i
         ko1, & !cel voisine de i
         i1,i2

    logical::linter
    integer::iml !last atom (=%im for standard; =%imm for extrait)


    atvois%iwmax(:)=0


    atvois%indi(:)=0
    
    if (present(lextr))lextrait=lextr
    if (present(lconstrtotR))then
       lconstrtt=lconstrtotR
    else
       lconstrtt=lconstrtot
    end if
       
    at=boxndm%at ; bg=boxndm%bg
    !
    !-----------------------------------------------
    ! --------------------------
    !   OUVERTURE BOUCLE SUR I
    ! --------------------------

    if(celvois%icaltabt.ne.atvois%icaltabt) then
       write (uwrt,*)'incoherence dans icaltabt caltabi'
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

    if (ipotentiel==20) then
      if (allocated(atvois%distance))deallocate(atvois%distance)
      allocate(atvois%distance(atvois%nVois, 4))
   end if
   ! write(*,*) 'testtttttt', ipotentiel


    !write(*,*) 'caltabi_inside  ', rvois, rvois2
    !*************construction par double boucle ****************
    if(lconstrtt) then  !construction par double boucle

       if (lextrait) then
          iml=atvois%imm
       else
          iml=atvois%im
       end if
       do i = 1, atvois%im
          iwo=iw
          iti=atvois%ityp(i)

          if(ldemitab)then
             ip = i+1
          else
             ip=1
          end if
          do j = ip, iml
             if(i.eq.j) cycle
             itj=atvois%ityp(j)
             ll=ipo(iti,itj)

            !  if (ipotentiel==20) then
            !    call vect_dist(atvois,celvois,boxndm,i,j,VJI,lperiod=lperiod,linter=linter,rum=rvois(ll))

             call vect_dist(atvois,celvois,boxndm,i,j,lperiod=lperiod,linter=linter,rum=rvois(ll))
             if (.not.linter)cycle             
             iw = iw+1
             IF (iw.GT.atvois%nVois) THEN
                WRITE(0,'(a,i0)') 'Indice iw du tableau de voisin plus grand&
                     & que le max, nVois = ', atvois%nVois
                !                   WRITE(0,'(a)') 'Augmentez le nombre moyen de voisins par&
                !                        & atome dans le fichier *.din'
                WRITE(0,'(a,i0)') 'truc étrange dans setcellconf'
                call arret_ndm(.true.)
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
          nvi=0
          iwo=iw
          koo = atvois%ielat(i)                          ! Numero de la cellule
          iti=atvois%ityp(i)
          !        write(uwrt,*)'atome i',i,iti
          ! pour chaque cel. voisine
          do i1 = 0, celvois%ncelvois(koo)
             ko1 = celvois%ncel(koo,i1)
             !           write(uwrt,*)'i1 ko1 ',i1,ko1
             if (ko1==0) cycle
             loop_j: do i2 = 1, celvois%nato(ko1)
                j = celvois%atincel(i2,ko1)
                !                                write(uwrt,*)'j ',j

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
                call vect_dist(atvois,celvois,boxndm,i,j,indcv=i1,lperiod=lperiod,linter=linter,rum=rvois(ll),VJI=VJI,dist=dij)
                if (.not.linter)cycle

                nvi=nvi+1
                iw = iw+1
                iwph = iwph+1
                !              write(uwrt,*)i,koo,ko1,j,iw, at,bg,sqrt(r2)
                if (present(cn2m)) then
                   call buildvoisext(cn2m,atvois%ityp(j),i,j,nvi,dij,VJI,koo,ko1)
                end if
                   
!                write(uwrt,*)rang,iw,size(atvois%indi)
                atvois%indi(iw) = j

                if (ipotentiel==20) then
                  atvois%distance(iw, 1) = dij
                  atvois%distance(iw, 2:4) = VJI(:)
                end if

                !              indi2(iwph) = j
             end do loop_j !i2
          end do !ncelvois
          atvois%iwmax(i) = iw
          nvij=iw-iwo
          !                 write(uwrt,*)'NVIJ',i,nvij,iw
       end do ! fin i
       maxvoi=iw

    endif ! lconstrtt

    return
  end subroutine caltabi

  subroutine buildvoisext(cn2m,itj,i,j,nvi,dij,VJI,koo,ko1)
    type(system_state)::cn2m
    integer::koo,ko1
    integer,intent(in)::i,j,nvi,itj
    real(double),intent(in)::dij,VJI(3)
    cn2m%n_neigh(i)=nvi ! ne devrait être fait qu'une fois à la fin de la boucle sur i , mais mis là pour éviter de polluer caltabi
    cn2m%r_ij(i,nvi)=dij
    cn2m%type_neigh(i,nvi)=itj
    cn2m%kind_neigh(i,nvi)=j
    cn2m%u_ij(i,nvi,1:3)=-VJI(1:3)
  end subroutine buildvoisext
end module caltabi_mod
