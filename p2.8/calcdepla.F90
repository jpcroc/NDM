module calcdepla_mod
  USE temp_com,only:zls2,at,bg,nad ! A EFFACER
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE var_pot, ONLY:ntyp,ty
  USE gen_com_m, ONLY:tdepla,lfilmext,iteration,timel,iko,lcasca,lfilm,rang
        implicit none 
        contains
! *******************************************************************
subroutine calcdepla(im,xp,ielat,ityp,ax)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double


#ifdef PARA
    USE mpi
    USE Tpara,only:MPI_COMM_space,status,ierr,nprocs,myidsp,NDM_MPI_REAl_DOUBLE
    use tab_imm_m,only:num_at_glob
#endif
  ! pas de conditions periodiques sur xp-ax
  !       version du 09 decembre 2003
  ! *******************************************************************

  implicit none
  integer,intent(in)::im
  real(double),allocatable::xp(:,:),ax(:,:)
  integer,allocatable::ityp(:),ielat(:)


  integer :: ndeplatot
  integer , dimension(ntyp) :: ndepla
  integer :: i, iti
  integer , dimension(im) :: indic
  real(double), dimension(im) :: distdepl
  integer :: lufilm, lufilmpaf, lufilmext, lutampon
  real(double), dimension(ntyp) :: dr2
  real(double) :: dri2, a1, a2, a3, c1, c2, c3, racdri2,depiko
  real(double), dimension(1,3) :: cv
  character :: fnamtampon*10, fnamfilmext*20, extension*10
  real(double), dimension(3) :: xp_iko
  integer :: ityp_iko
#ifdef PARA
  integer,      allocatable :: ityp_depla(:)
  integer,      allocatable :: indic_depla(:)
  real(double), allocatable :: xp_depla(:,:)
  real(double), allocatable :: dist_depla(:)
  real(double), dimension(ntyp) :: dr2_glob
  integer , dimension(ntyp) :: ndepla_glob
  integer :: ndeplatot_glob
  integer :: est_present
  integer :: ndeplatot_tmp
  integer :: proc_source
#endif
  !-----------------------------------------------
  !
  ! local variables
  !
  !
  !

  lufilm = 89                                ! index fichier film pour toutes les iterations
  lufilmpaf = 79                             ! index fichier film du paf pour toutes les iterations
  lufilmext = 69                             ! index fichiers positions pour une iterationterations



  !      write(6,*)'entree dans calcdepla'
  if (rang==0) then
     write (6, *)
     write (6, *) '----------- Deplacements ----------------'
  end if
!       write(6,*)'tdepla',tdepla
  ndeplatot = 0
  dr2(:ntyp) = 0.0
  ndepla(:ntyp) = 0
     call cryst_to_cart (im, xp, bg, -1)    !cart vers cryst
     call cryst_to_cart (im, ax, bg, -1)    !cart vers cryst

     do i = 1, im
        c1 = ax(1,i)-xp(1,i)
        c2 = ax(2,i)-xp(2,i)
        c3 = ax(3,i)-xp(3,i)
        ! pas de conditions périodiques sur les déplacements
        !            if (c1>0.5) c1 = c1-1.
        !            if (c1<(-0.5)) c1 = c1+1.
        !            if (c2>0.5) c2 = c2-1.
        !            if (c2<(-0.5)) c2 = c2+1.
        !            if (c3>0.5) c3 = c3-1.
        !            if (c3<(-0.5)) c3 = c3+1.
        cv(1,1) = c1
        cv(1,2) = c2
        cv(1,3) = c3
        call cryst_to_cart (1, cv, at, 1)    !cryst vers cart sur cv
        dri2 = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
        dr2(ityp(i)) = dr2(ityp(i))+dri2/nad(ityp(i))
        racdri2 = sqrt(dri2)
        if (racdri2<tdepla) cycle
        ndepla(ityp(i)) = ndepla(ityp(i))+1
        ndeplatot = ndeplatot+1
        indic(ndeplatot) = i
        distdepl(ndeplatot)=racdri2*1.0d8
     end do
     call cryst_to_cart (im, xp, at, 1)     !cryst vers cart
     call cryst_to_cart (im, ax, at, 1)     !cryst vers cart



#ifdef PARA
  call MPI_ALLREDUCE(dr2,dr2_glob,ntyp,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
  dr2 = dr2_glob
  call MPI_ALLREDUCE(ndepla,ndepla_glob,ntyp,MPI_INTEGER,MPI_SUM,MPI_COMM_space,ierr)
  ndepla = ndepla_glob
  call MPI_ALLREDUCE(ndeplatot,ndeplatot_glob,1,MPI_INTEGER,MPI_SUM,MPI_COMM_space,ierr)

  ! Allocation des tableaux d'emission/reception
  if (myidsp==0) then
     allocate(ityp_depla(ndeplatot_glob))
     allocate(xp_depla(1:3,ndeplatot_glob))
     allocate(indic_depla(ndeplatot_glob))
     allocate(dist_depla(ndeplatot_glob))
  else
     allocate(ityp_depla(ndeplatot))
     allocate(xp_depla(1:3,ndeplatot))
     allocate(indic_depla(ndeplatot))
     allocate(dist_depla(ndeplatot))
  endif
  ! Preparation des tableaux avec les donnees locales
  do i=1,ndeplatot
     ityp_depla(i)=ityp(indic(i))
     xp_depla(:,i)=xp(:,indic(i))
     indic_depla(i)=num_at_glob(indic(i))
     dist_depla(i)=distdepl(i)
  enddo

  ! Recuperation par l'ensemble des procs des différents deplacements

  if (myidsp==0) then
     do i=1,nprocs-1
        call MPI_RECV(ndeplatot_tmp,1,MPI_INTEGER,MPI_ANY_SOURCE,13001,MPI_COMM_space,status,ierr)
        if (ndeplatot_tmp.ne.0) then
           proc_source = status(MPI_SOURCE)
           call MPI_RECV(ityp_depla(ndeplatot+1),ndeplatot_tmp,MPI_INTEGER,proc_source,13002,MPI_COMM_space,status,ierr)
           call MPI_RECV(xp_depla(:,ndeplatot+1),3*ndeplatot_tmp,NDM_MPI_REAL_DOUBLE,proc_source,13003,MPI_COMM_space,status,ierr)
           call MPI_RECV(indic_depla(ndeplatot+1),ndeplatot_tmp,MPI_INTEGER,proc_source,13004,MPI_COMM_space,status,ierr)
           call MPI_RECV(dist_depla(ndeplatot+1),ndeplatot_tmp,NDM_MPI_REAL_DOUBLE,proc_source,13005,MPI_COMM_space,status,ierr)
           ndeplatot = ndeplatot + ndeplatot_tmp
        endif
     enddo
  else
     call MPI_SEND(ndeplatot,1,MPI_INTEGER,0,13001,MPI_COMM_space,ierr)
     if (ndeplatot.ne.0) then
        call MPI_SEND(ityp_depla,ndeplatot,MPI_INTEGER,0,13002,MPI_COMM_space,ierr)
        call MPI_SEND(xp_depla(:,1:ndeplatot),3*ndeplatot,NDM_MPI_REAL_DOUBLE,0,13003,MPI_COMM_space,ierr)
        call MPI_SEND(indic_depla,ndeplatot,MPI_INTEGER,0,13004,MPI_COMM_space,ierr)
        call MPI_SEND(dist_depla(1:ndeplatot),ndeplatot,NDM_MPI_REAL_DOUBLE,0,13005,MPI_COMM_space,ierr)
     endif
  endif
#endif


  if (rang==0) then
     write (6, '(A,I5,A,D10.3)') '*  ITERATION  = ', iteration, '  time = ', timel
        write (6, *) 'nombre total d-atomes deplaces = ', ndeplatot

     write (6, *)
     do iti = 1, ntyp
        if (nad(iti)==0) cycle
        write (6, *) 'deplacement moyen des atomes de type ', iti, ' = ', &
             sqrt(dr2(iti)*1D+16)
        write (6, *) 'nombre d-atomes de type ', iti, ' deplaces = ', ndepla(&
             iti)
     end do



  ! ***** Ecriture positions formattees dans un seul fichier *****
  ! ***** Fichier employe pour traitement avec xmakemol *****
  if (lfilm) then
     ! Ecriture des types et coordonnees des atomes deplaces de plus de
     ! tdepla angstroems dans le fichier film
     write (lufilm, '(I7,A,I7,A,2D10.3)')  ndeplatot+2, ' IT =', iteration, ' Time = ', timel
     at=at*1.d8
     write (lufilm,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
     at=at/1.d8 
     write (lufilm, 114)  zls2(1)*1d8,zls2(2)*1d8,zls2(3)*1d8
     write (lufilm, 114) -zls2(1)*1d8,-zls2(2)*1d8,-zls2(3)*1d8

#ifdef PARA
     do i = 1, ndeplatot_glob
        !       write(6,*)i,indic(i)
        write (lufilm, 113) ty(ityp_depla(i)), xp_depla(1,i)*1D+8, xp_depla(2,&
             i)*1D+8, xp_depla(3,i)*1D+8, indic_depla(i),dist_depla(i)
     end do
#else
     do i = 1, ndeplatot
        !       write(6,*)i,indic(i)
        write (lufilm, 113) ty(ityp(indic(i))), xp(1,indic(i))*1D+8, xp(2,&
             indic(i))*1D+8, xp(3,indic(i))*1D+8, indic(i),distdepl(i)
     end do
#endif
  endif
end if
     ! Ecriture coordonnees du premier atome frappe toutes les itedepla
     ! iterations
     if (lcasca.and.lfilm) then
#ifdef PARA
        ! Recherche du proc possedant iko
        est_present=0
        do i=1, im
           if (num_at_glob(i)==iko) then
              est_present=1
              exit
           endif
        enddo
        ! Emission/reception des infos vers le proc 0
        if (est_present==1.and.myidsp==0) then
           xp_iko(:) = xp(:,i)
           ityp_iko = ityp(i)
        else if (est_present==1) then
           call MPI_SEND(xp(:,i),3,NDM_MPI_REAL_DOUBLE,0,13006,MPI_COMM_space,ierr)
           call MPI_SEND(ityp(i),1,MPI_INTEGER,0,13007,MPI_COMM_space,ierr)
        else if (myidsp==0) then
           call MPI_RECV(xp_iko(:),3,NDM_MPI_REAL_DOUBLE,MPI_ANY_SOURCE,13006,MPI_COMM_space,status,ierr)
           call MPI_RECV(ityp_iko,1,MPI_INTEGER,MPI_ANY_SOURCE,13007,MPI_COMM_space,status,ierr)
        endif
#else
        xp_iko(:) = xp(:,iko)
        ityp_iko = ityp(iko)
#endif
        if(rang==0)then
           write (lufilmpaf, *) '   1'
           write (lufilmpaf, *) ' IT', iteration, ' time ', timel
           write (lufilmpaf, 113) ty(ityp_iko), xp_iko(1)*1D+8, xp_iko(2)*1D+8, &
                xp_iko(3)*1D+8, iko,distdepl(iko)
        end if
     endif


  ! ***** fin ecriture dans un seul fichier positions *****


  ! ***** Ecriture positions atomiques dans fichiers differents *****
  ! ***** 1 iteration -> 1 fichier pour traitement images animees gif *****
  if (lfilmext) then

     ! conversion entier-->alphanumerique par transfert du nombre
     ! de l'iteration vers fichier tampon relu sous format caractere.

     lutampon = 17
     fnamtampon = 'tampon'
     if(rang==0)then
        write(extension,'(i10.10)') iteration        


        !    ouverture d'un fichier filmext.(iteration) pour sauvegarde
        !    des positions toutes les itedepla iterations.

        fnamfilmext = 'film1it.'//extension
        open(unit=lufilmext, file=fnamfilmext, form='formatted', status=&
             'unknown')
        ! Ecriture des types et coordonnees des atomes deplaces de plus de
        ! tdepla angstroems dans le fichier film.extension
        write (lufilmext, *) ndeplatot+2
        write (lufilmext, *) ' IT', iteration, ' time ', timel
        write (lufilmext, 114)  zls2(1)*1d8,zls2(2)*1d8,zls2(3)*1d8
        write (lufilmext, 114) -zls2(1)*1d8,-zls2(2)*1d8,-zls2(3)*1d8
114     format('H ',1x,3(f10.4,1x))
#ifdef PARA
        do i = 1, ndeplatot_glob
           !     write(6,*)i,indic(i),iko
           if (indic_depla(i)==iko) then
              write (lufilmext, 112)  ty(ityp_depla(i)),xp_depla(1,i)*1d+8,xp_depla(2,i)*1d+8, xp_depla(&
                   3,i)*1d+8, indic_depla(i)
           else
              !              write (lufilmext, 113) ty(ityp(indic(i))), xp(1,indic(i))*1D+8, xp(&
              !                   2,indic(i))*1D+8, xp(3,indic(i))*1D+8, indic(i),distdepl(indic(i))
              ! Correction d'indice pour distdepl
              write (lufilmext, 113) ty(ityp_depla(i)), xp_depla(1,i)*1D+8, xp_depla(&
                   2,i)*1D+8, xp_depla(3,i)*1D+8, indic_depla(i),dist_depla(i)
           endif
        end do
#else
        do i = 1, ndeplatot
           !     write(6,*)i,indic(i),iko
           if (indic(i)==iko) then
              write (lufilmext, 112)  ty(ityp(indic(i))),xp(1,indic(i))*1d+8,xp(2,indic(i))*1d+8, xp(&
                   3,indic(i))*1d+8, indic(i)
           else
              write (lufilmext, 113) ty(ityp(indic(i))), xp(1,indic(i))*1D+8, xp(&
                   2,indic(i))*1D+8, xp(3,indic(i))*1D+8, indic(i),distdepl(indic(i))
           endif
        end do
#endif
        close(lufilmext)
     endif                                      ! lfilmext=TRUE

     ! ***** Fin ecriture positions dans plusieurs fichiers *****
  end if
112 format(a3,1x,3(f10.4,1x),1x,1x,i6,' PKA')
113 format(a3,1x,3(f10.4,1x),1x,1x,i6,G12.4)

#ifdef PARA
  ! liberation des tableaux
  deallocate(ityp_depla)
  deallocate(xp_depla)
  deallocate(indic_depla)
  deallocate(dist_depla)
#endif

end subroutine calcdepla
end module
