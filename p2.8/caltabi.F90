!****************************************************************



! *****************************************************************
subroutine caltabi
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m
  !           version du 4 juin 2010, 14h38 - last chaged by MCM
  ! *****************************************************************

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !----------------------------------------------1-
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: iw, i, ip, j, maxvoi, nvij,iwo
  integer :: itj,ll
  REAL(double) :: r2
  real(double), dimension(1:npair) :: rvois2

  real(double), dimension(3) :: xpi, dx, ds  

  integer :: iti, & !type de i
       koo, & !cel de i
       ncelvois,ko1, & !cel voisine de i
       i1,i2


  real(double),dimension(:,:),allocatable :: xpnp  ! MODIF Cosmin                            

  !-----------------------------------------------
  ! --------------------------
  !   OUVERTURE BOUCLE SUR I
  ! --------------------------
  iw = 0

  if (ipotentiel==12) then
     rvois2(1)=(7.0d-8)**2
     rvois2(2)=(3.5d-8)**2
     rvois2(3)=(2.8d-8)**2
  else
     rvois2(:)=rvois**2 
  end if

  nvij=0
  
  ALLOCATE(xpnp(3,imm))
  if (lperiod) then
    xpnp(:,:)=xp(:,:)
   else
   call notperiod(xp,xpnp)
  end if  
  
  
  !*************construction par double boucle ****************
  if(lconstrtot) then  !construction par double boucle


     do i = 1, im
        iwo=iw
        xpi(:) = xpnp(:,i)
        iti=ityp(i)

        if(ldemitab)then
           ip = i+1
        else
           ip=1
        end if

        do j = ip, im
           if(i.eq.j) cycle
           dx(:) = xpi(:) - xpnp(:,j)
           ds(:) = MatMul( dx(:), bg(:,:) )
           
	   if ((ds(1)>0.5d0).or.(ds(1)<-0.5d0))   ds(1) = ds(1)-dble(Nint(ds(1)))
           if ((ds(2)>0.5d0).or.(ds(2)<-0.5d0))   ds(2) = ds(2)-dble(Nint(ds(2)))
           if ((ds(3)>0.5d0).or.(ds(3)<-0.5d0))   ds(3) = ds(3)-dble(Nint(ds(3)))

           dx(:) = MatMul( at(:,:), ds(:) )
           r2 = Sum( dx(:)**2 )

           itj=ityp(j) 
           ll=ipo(iti,itj)

           if (r2>rvois2(ll)) cycle

           iw = iw+1
           IF (iw.GT.nVois) THEN
                   WRITE(0,'(a,i0)') 'Indice iw du tableau de voisin plus grand&
                        & que le max, nVois = ', nVois
                   WRITE(0,'(a)') 'Augmentez le nombre moyen de voisins par&
                        & atome dans le fichier *.din'
                   WRITE(0,'(a,i0)') 'valeur actuelle: nvperat = ', nvperat
                   STOP '< Caltabi >'
           END IF

           indi(iw) = j
        end do
        iwmax(i) = iw
        nvij=iw-iwo

     end do   ! im 
     maxvoi = iw           

  !*************construction par celulle ****************
  else   


     do i = 1, im
        iwo=iw
        koo = ielat(i)                          ! Numero de la cellule
        iti=ityp(i) 
        xpi(:) = xpnp(:,i)

        ncelvois = min(noxyz,27)-1
        ! pour chaque cel. voisine
        do i1 = 0, ncelvois
           ko1 = ncel(koo,i1)
           !               write(6,*)'i1 ko1 ',i1,ko1
           if (ko1==0) cycle
           loop_j: do i2 = 1, nato(ko1)
              j = last(i2,ko1)
              !                  write(6,*)'j ',j

              if(ldemitab) then
                 if(j.le.i) cycle !terme deja calcule
              else
                 if(j.eq.i) cycle
              end if
              itj=ityp(j) ; ll=ipo(iti,itj)

              dx(:) = xpi(:) - xpnp(:,j)
              ds(:) = MatMul( dx(:), bg(:,:) )
              if ((ds(1)>0.5d0).or.(ds(1)<-0.5d0))   ds(1) = ds(1)-dble(Nint(ds(1)))
              if ((ds(2)>0.5d0).or.(ds(2)<-0.5d0))   ds(2) = ds(2)-dble(Nint(ds(2)))
              if ((ds(3)>0.5d0).or.(ds(3)<-0.5d0))   ds(3) = ds(3)-dble(Nint(ds(3)))
              dx(:) = MatMul( at(:,:), ds(:) )
              r2 = Sum( dx(:)**2 )

              if (r2>rvois2(ll)) cycle
              iw = iw+1
              !                  write(6,*)i,koo,ko1,j,iw
              indi(iw) = j
           end do loop_j !i2

        end do !ncelvois
        iwmax(i) = iw
        nvij=iw-iwo
        !         write(449,*)'NVIJ',i,nvij,iw       

     end do ! fin i
     maxvoi=iw





  endif ! lconstrtot 

   if ((rang==0).and.(it.le.100)) then
!           write(6,*)'IT ',it,'  VOISINS ',maxvoi,' par atome ',float(maxvoi)/float(im)
   endif
  DEALLOCATE (xpnp)  
  return
end subroutine caltabi

