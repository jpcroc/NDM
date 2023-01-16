module param_det_mod
  USE arret_ndm_mod,only:arret_ndm
  USE arret_ndm_mod,only: arret_ndm
  USE gen_com_m, ONLY:lopt,zero,rang,pi,itab
  USE var_pot, ONLY:kpme,kpmex,kpmey,kpmez,n2max,ncouc3,ncoucx,ncoucy,ncoucz,npair,&
       &npotentiel,nvecttot,precisew,rue_pair,typ_pot_pair,ipotentiel,alpha,iewald,csive,&
       &ngrid,r3cm2,rumax,r3cm,tabv3,tabf3,ntyp,l3c
  use boxconfig,only:box_config
  use read_val,only:ltabvois,rvois
  implicit none
contains
  subroutine param_det(boxndm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    !determine rumax et rue si pas défini


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
    type(box_config)::boxndm
    integer :: nc1,l
    real(double) :: zlm, pparam, k00x, k00y, k00z,k001,rue
    real(double) :: ruex, ruey, ruez, pparax, pparay, pparaz,zl(3)
#ifdef ixia
    real(double) :: alpha_ixia
#endif


    if (npotentiel.ne.1)then
       if ((iewald.gt.0).and.(iewald.ne.3).and.(ncouc3==0)) then
          if (rang==0)  write(6,*)'npot>1 + ewald+ncouc3=0 : stop'
          call arret_ndm
       end if
       rue=0
       do l=1,npair
          if ((typ_pot_pair(l).lt.10).and.(typ_pot_pair(l).ne.2)) then
             if (rue_pair(l)==0) then
                if (rang==0)  write(6,*)'npot>1 + pot paire +ruepaire l =0 : stop',l
                call arret_ndm
             end if
          end if
          rue=max(rue,rue_pair(l))
       end do
       if ((iewald.gt.0).and.(iewald.ne.3))then
          if((alpha==0).or.(ncouc3==0)) then
             write(6,*)'npotentiel>1 and Ewald : specify ncouc3 and alpha'
             call arret_ndm
          end if
          ncoucx=ncouc3
          ncoucy=ncouc3
          ncoucz=ncouc3
          if (.not.allocated(tabv3)) then

             allocate (tabv3(-ncoucx:ncoucx,-ncoucy:ncoucy,-ncoucz:ncoucz))
             allocate (tabf3(ntyp,-ncoucx:ncoucx,-ncoucy:ncoucy,-ncoucz:ncoucz))
          end if

       end if

    else
       rue=0
       if (((ipotentiel.lt.10).and.(ipotentiel.ne.2)).or.(ipotentiel==16)) then
          do l=1,npair
             if ((typ_pot_pair(l).lt.10).or.(typ_pot_pair(l)==16)) then
                rue=max(rue,rue_pair(l))
             end if
          end do

          if  ((iewald.gt.0).and.(iewald.ne.3))then

             !  rue=rue*1.d-8
             !        alpha=alpha*1.d8
             !        write(6,*)'ZL',zl
             zl=boxndm%normat
             zlm=max(zl(1),zl(2),zl(3))
             k001=2.d0*pi/zlm

             k00x=2.d0*pi/zl(1)
             k00y=2.d0*pi/zl(2)
             k00z=2.d0*pi/zl(3)

             if (ncouc3 == 0 .and. ncoucx/=0 .and. (ncoucy==0 .or. ncoucz==0)) then
                write (6,*) rang, 'Parametres ncouc de la sommation d Exald mal definis'
                call arret_ndm
             endif

             if (ncouc3 == 0 .and. ncoucy/=0 .and. (ncoucx==0 .or. ncoucz==0)) then
                write (6,*) rang,'Parametres ncouc de la sommation d Exald mal definis'
                call arret_ndm
             endif

             if (ncouc3 == 0 .and. ncoucz/=0 .and. (ncoucx==0 .or. ncoucy==0)) then
                write (6,*) rang,'Parametres ncouc de la sommation d Exald mal definis'
                call arret_ndm
             endif

             if (ncoucx .lt. 0 .or. ncoucy .lt. 0 .or. ncoucz .lt. 0) then
                write (6,*) rang,'Parametres ncouc de la sommation d Exald mal definis'
                call arret_ndm
             endif

             if (ncoucx == 0 .and. ncoucy == 0 .and. ncoucz == 0) then
                nc1 = 0
             else
                nc1 = min(ncoucx, ncoucy, ncoucz)
             endif


             !        if (.not. lopt) then

             if (precisew == zero) then

                if (rue == zero .and. alpha==zero .and. ncouc3==0 .and. nc1 ==0) then
                   precisew=1.0d-3
                   pparam=-log(precisew)
                   rue=12.0d-8
                   alpha=dsqrt(pparam)/rue
                   if (iewald/=0) then
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   endif
                   goto 1
                endif

                if (rue == zero .and. alpha==zero .and. ncouc3/=0) then
                   if (iewald/=0) then
                      precisew=1.0d-3
                      pparam=-log(precisew)
                      ruex=2.d0*pparam/ncouc3/k00x
                      ruey=2.d0*pparam/ncouc3/k00y
                      ruez=2.d0*pparam/ncouc3/k00z
                      rue=max(ruex,ruey,ruez)
                      alpha=dsqrt(pparam)/rue
                      ncoucx=ncouc3
                      ncoucy=ncouc3
                      ncoucz=ncouc3
                   else
                      write(6,*) rang,'iewald=0 et ncouc3>0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue == zero .and. alpha==zero .and. ncouc3==0 .and. nc1/=0) then
                   if (iewald/=0) then
                      precisew=1.0d-3
                      pparam=-log(precisew)
                      ruex=2.d0*pparam/ncoucx/k00x
                      ruey=2.d0*pparam/ncoucy/k00y
                      ruez=2.d0*pparam/ncoucz/k00z
                      rue=max(ruex,ruey,ruez)
                      alpha=dsqrt(pparam)/rue
                   else
                      write(6,*) rang,'iewald=0 et ncouc/=0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue == zero .and. alpha/=zero .and. ncouc3==0 .and. nc1==0) then
                   precisew=1.0d-3
                   pparam=-log(precisew)
                   rue=dsqrt(pparam)/alpha
                   if (iewald/=0) then
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   endif
                   goto 1
                endif

                if (rue /= zero .and. alpha==zero .and. ncouc3==0 .and. nc1==0) then
                   precisew=1.0d-3
                   pparam=-log(precisew)
                   alpha=dsqrt(pparam)/rue
                   if (iewald/=0) then
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   endif
                   goto 1
                endif

                if (rue /= zero .and. alpha/=zero .and. ncouc3==0 .and. nc1==0) then
                   pparam=(rue*alpha)**2
                   if (iewald/=0) then
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   endif

                   goto 1
                endif

                if (rue /= zero .and. alpha==zero .and. ncouc3/=0) then
                   if (iewald/=0) then
                      pparax=rue*ncouc3*k00x/2.d0
                      pparay=rue*ncouc3*k00y/2.d0
                      pparaz=rue*ncouc3*k00z/2.d0
                      pparam=max(pparax,pparay,pparaz)
                      alpha=dsqrt(pparam)/rue
                      ncoucx=ncouc3
                      ncoucy=ncouc3
                      ncoucz=ncouc3
                   else
                      write(6,*) rang,'iewald=0 et ncouc3>0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue /= zero .and. alpha==zero .and. ncouc3==0 .and. nc1/=0) then
                   if (iewald/=0) then
                      pparax=rue*ncoucx*k00x/2.d0
                      pparay=rue*ncoucy*k00y/2.d0
                      pparaz=rue*ncoucz*k00z/2.d0
                      pparam=max(pparax,pparay,pparaz)
                      alpha=dsqrt(pparam)/rue
                   else
                      write(6,*) rang,'iewald=0 et ncouc3>0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue == zero .and. alpha/=zero .and. ncouc3/=0) then
                   if (iewald/=0) then
                      ruex=ncouc3*k00x/alpha**2/2.d0
                      ruey=ncouc3*k00y/alpha**2/2.d0
                      ruez=ncouc3*k00z/alpha**2/2.d0
                      rue=max(ruex,ruey,ruez)
                      ncoucx=ncouc3
                      ncoucy=ncouc3
                      ncoucz=ncouc3
                   else
                      write(6,*) rang,'iewald=0 et ncouc3>0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue == zero .and. alpha/=zero .and. ncouc3==0 .and. nc1/=0) then
                   if (iewald/=0) then
                      ruex=ncoucx*k00x/alpha**2/2.d0
                      ruey=ncoucy*k00y/alpha**2/2.d0
                      ruez=ncoucz*k00z/alpha**2/2.d0
                      rue=max(ruex,ruey,ruez)
                   else
                      write(6,*) rang,'iewald=0 et ncouc/=0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue /= zero .and. alpha/=zero .and. ncouc3/=0) then
                   ncoucx=ncouc3
                   ncoucy=ncouc3
                   ncoucz=ncouc3
                   if (iewald==0) then
                      write(6,*) rang,'iewald=0 et ncouc3>0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue /= zero .and. alpha/=zero .and. ncouc3==0 .and. nc1/=0) then
                   if (iewald==0) then
                      write(6,*) rang,'iewald=0 et ncouc/=0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                !         endif

             else !  if (precisew /= zero) then

                if (rue == zero .and. alpha==zero .and. ncouc3==0 .and. nc1==0) then
                   pparam=-log(precisew)
                   rue=12.0d-8
                   alpha=dsqrt(pparam)/rue
                   if (iewald/=0) then
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   endif

                   goto 1
                endif

                if (rue == zero .and. alpha==zero .and. ncouc3/=0) then
                   if (iewald/=0) then
                      pparam=-log(precisew)
                      ruex=2.d0*pparam/ncouc3/k00x
                      ruey=2.d0*pparam/ncouc3/k00y
                      ruez=2.d0*pparam/ncouc3/k00z
                      rue=max(ruex,ruey,ruez)
                      alpha=dsqrt(pparam)/rue
                      ncoucx=ncouc3
                      ncoucy=ncouc3
                      ncoucz=ncouc3
                   else
                      write(6,*) rang,'iewald=0 et ncouc3>0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue == zero .and. alpha==zero .and. ncouc3==0 .and. nc1/=0) then
                   if (iewald/=0) then
                      pparam=-log(precisew)
                      ruex=2.d0*pparam/ncoucx/k00x
                      ruey=2.d0*pparam/ncoucy/k00y
                      ruez=2.d0*pparam/ncoucz/k00z
                      rue=max(ruex,ruey,ruez)
                      alpha=dsqrt(pparam)/rue
                   else
                      write(6,*) rang,'iewald=0 et ncouc/=0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue == zero .and. alpha/=zero .and. ncouc3==0 .and. nc1==0) then
                   pparam=-log(precisew)
                   rue=dsqrt(pparam)/alpha
                   if (iewald/=0) then
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   endif

                   goto 1
                endif

                if (rue /= zero .and. alpha==zero .and. ncouc3==0 .and. nc1==0) then
                   pparam=-log(precisew)
                   alpha=dsqrt(pparam)/rue
                   if (iewald/=0) then
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   endif

                   goto 1
                endif

                if (rue /= zero .and. alpha/=zero .and. ncouc3==0 .and. nc1==0) then
                   if (rang==0) &
                        write(6,*) 'Seuls la precision et rue sont pris en compte'
                   pparam=-log(precisew)
                   alpha=dsqrt(pparam)/rue
                   if (iewald/=0) then
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   endif

                   goto 1
                endif

                if (rue /= zero .and. alpha==zero .and. ncouc3/=0) then
                   if (iewald/=0) then
                      if (rang==0) &
                           write(6,*) 'Seuls la precision et rue sont pris en compte'
                      pparam=-log(precisew)
                      alpha=dsqrt(pparam)/rue
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   else
                      write(6,*) rang,'iewald=0 et ncouc3>0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue /= zero .and. alpha==zero .and. ncouc3==0 .and. nc1/=0) then
                   if (iewald/=0) then
                      if (rang==0) &
                           write(6,*) 'Seuls la precision et rue sont pris en compte'
                      pparam=-log(precisew)
                      alpha=dsqrt(pparam)/rue
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   else
                      write(6,*) rang,'iewald=0 et ncouc/=0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue == zero .and. alpha/=zero .and. ncouc3/=0) then
                   if (iewald/=0) then
                      if (rang==0) &
                           write(6,*) 'Seuls la precision et alpha sont pris en compte'
                      pparam=-log(precisew)
                      rue=dsqrt(pparam)/alpha
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   else
                      write(6,*) rang,'iewald=0 et ncouc3>0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue == zero .and. alpha/=zero .and. ncouc3==0 .and. nc1/=0) then
                   if (iewald/=0) then
                      if (rang==0) &
                           write(6,*) 'Seuls la precision et alpha sont pris en compte'
                      pparam=-log(precisew)
                      rue=dsqrt(pparam)/alpha
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   else
                      write(6,*) rang,'iewald=0 et ncouc/=0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue /= zero .and. alpha/=zero .and. ncouc3/=0) then
                   if (iewald/=0) then
                      if (rang==0) &
                           write(6,*) 'Seuls la precision et rue sont pris en compte'
                      pparam=-log(precisew)
                      alpha=dsqrt(pparam)/rue
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   else
                      write(6,*) rang,'iewald=0 et ncouc3>0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

                if (rue /= zero .and. alpha/=zero .and. ncouc3==0 .and. nc1/=0) then
                   if (iewald/=0) then
                      if (rang==0)  &
                           write(6,*) 'Seuls la precision et rue sont pris en compte'
                      pparam=-log(precisew)
                      alpha=dsqrt(pparam)/rue
                      ncoucx=int(2.d0*pparam/rue/k00x)+1
                      ncoucy=int(2.d0*pparam/rue/k00y)+1
                      ncoucz=int(2.d0*pparam/rue/k00z)+1
                   else
                      write(6,*) rang,'iewald=0 et ncouc/=0 : incoherent !!!'
                      call arret_ndm
                   endif
                   goto 1
                endif

             endif


1            continue


             n2max=ncouc3*ncouc3

!!$           if (rang==0) &
!!$                write(6,*) 'Parametres utilises pour le traitement de EWALD :'

             !determination de kpme
             if (iewald == 2) then
                if (kpmex == 0) kpmex=int(2.5*(2*ncoucx+1))
                if (kpmey == 0) kpmey=int(2.5*(2*ncoucy+1))
                if (kpmez == 0) kpmez=int(2.5*(2*ncoucz+1))
                !          kpmex=100 ; kpmey=100 ; kpmez=100
                kpme=max(kpmex,kpmey,kpmez)
                if (rang==0) write(6,*) 'kpme x,y,x max ', kpmex,kpmey,kpmez,kpme
             endif

             !          if (kpmex <= 2*ncouc3) then
             !           if (rang==0) write(6,*) 'PME : Taille de grille trop faible STOP!'
             !          call arret_ndm


             if (rang==0) then
                if (iewald/=0) then
                   write(6,*) 'RUE=',rue,' ALPHA=',alpha,' NCOUC=',ncoucx,ncoucy,ncoucz,&
                        ' PRECISEW =',precisew
                   write(6,*) 'RUE_ANG=',rue*1d8,' ALPHA_ANGm1=',alpha*1d-8,' NCOUC=',ncoucx,ncoucy,ncoucz,&
                        ' PRECISEW =',precisew
                   write(6,*)'ncoucx_y_z',ncoucx,ncoucy,ncoucz


                else
                   write(6,*) 'RUE=',rue
                endif
             endif   ! rang = 0

             if ((iewald/=0).and.(.not.allocated(tabv3))) then

                allocate (tabv3(-ncoucx:ncoucx,-ncoucy:ncoucy,-ncoucz:ncoucz))
                allocate (tabf3(ntyp,-ncoucx:ncoucx,-ncoucy:ncoucy,-ncoucz:ncoucz))
             end if

             nvecttot=(2*ncoucx+1)*(2*ncoucy+1)*(2*ncoucz+1)-1
             !        write(6,*)'nvecttot',nvecttot

          end if

          rue_pair(:)=rue
       endif
    end if



    !C_debug
#if defined PHONDY || defined PARAPH || defined MAB || defined ML || defined PARAML
    !  write(6,*)rumax,rue_pair,maxval(rue_pair)
    rumax=0.0

    rumax = max(rumax,maxval(rue_pair))
    csive=rumax/float(ngrid)
    ! write(*,*) rumax, rvois
#else
    rumax = max(rumax,maxval(rue_pair))
    csive=rumax/float(ngrid)
    !    write(6,*)'RUMAX',rumax,csive,rue_pair
    !    call arret_ndm
#endif

    if(l3c) then
       if (r3cm.eq.0) r3cm=5.0d-8
       itab=1
       r3cm2=r3cm**2
    endif
    !if((rang==0).and.(appel==0))         write (6, *) ' rvois  ', rvois

    if (ltabvois) then
       if (rumax>rvois) then
          if((rang==0)) write (6, '(A,2F12.2)') ' rvois trop petit rvois rumax ', rvois*1d8, rumax*1d8
          call arret_ndm
       else
          if(rang==0) write (6,'(A,2F12.2)') ' rumax devient rvois&
               & pour le dimmensionnement en cel', rvois*1d8, rumax*1d8
          rumax=rvois
       end if


    endif

    return
  end subroutine param_det
end module param_det_mod
