module neb_module
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use contrainte

  !-----------------------------------------------
  implicit none


  integer, save                                  :: dragtest
  integer,dimension(:),allocatable,save          :: irelax,nebtest,icontrainte
  integer,dimension(:,:),allocatable,save        :: ielat_n, iwmax_n, ityp_n
  real(double),dimension(:,:,:),allocatable,save ::  xp_n, xpp_n, vp_n, ax_n, fp_n
  real(double),dimension(:,:),allocatable        :: fp_par,fp_perp
  real(double), dimension(:),allocatable, save   :: enePATH,enePATHev,norms,reaction_coord 
  real(double), dimension(:,:,:),allocatable,save:: s_path,force_neb,bruitneb 
  real(double)                                   :: forctot,formax,formaxperp, &
       formaxparl,masstot



contains 

  subroutine allocate_neb()
    implicit none
    allocate (ielat_n(imm,npath), iwmax_n(imm,npath), &
         ityp_n(imm,npath),&
         irelax(imm),                            &
	 icontrainte(imm),                       &
         xp_n(3,imm,npath),                      &
         xpp_n(3,imm,npath),                     &
         vp_n(3,imm,npath),                      &
         ax_n(3,imm,npath),                      &
         fp_n(3,imm,npath),                      &
         reaction_coord(npath))
    allocate  (enePATH(npath),enePATHev(npath),norms(npath),nebtest(npath))
    allocate  (fp_par(3,imm),fp_perp(3,imm))
    allocate  (s_path(3,imm,npath),force_neb(3,imm,npath))
    allocate   (bruitneb(3,imm,npath))

    return

  end subroutine allocate_neb

  ! **************************************************************
  subroutine into_path(iph, i_dir_path, &                         
       xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: ielat(imm)
    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    real(double)  :: fp(3,imm)
    !-----------------------------------------------
    !-----------------------------------------------
    integer   :: iph, i_dir_path
    !-----------------------------------------------


    if (i_dir_path==1) then
       ielat_n   (:,iph) = ielat   (:) 
       iwmax_n   (:,iph) = iwmax   (:)
       ityp_n    (:,iph) = ityp    (:)
       xp_n (:,:,iph)    = xp (:,:)
       xpp_n(:,:,iph)    = xpp(:,:)
       vp_n (:,:,iph)    = vp (:,:)
       !       vp_n (:,:,iph)    = 0.d0
       ax_n (:,:,iph)    = ax (:,:)
       fp_n (:,:,iph)    = fp (:,:)
    else 
       ielat   (:)      = ielat_n   (:,iph)
       iwmax   (:)      = iwmax_n   (:,iph)
       ityp    (:)      = ityp_n    (:,iph)
       xp (1:3,:)	      = xp_n (1:3,:,iph)   
       xpp(1:3,:)	      = xpp_n(1:3,:,iph)   
       vp (1:3,:)	      = vp_n (1:3,:,iph)   
       ax (1:3,:)	      = ax_n (1:3,:,iph)   
       fp (1:3,:)	      = fp_n (1:3,:,iph)   
    end if

    return

  end subroutine into_path



  subroutine init_neb(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: ielat(imm)
    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    real(double)  :: fp(3,imm)
    real(double)  :: dxx(3,imm)
    !-----------------------------------------------
    integer :: iph,ic,non_contr,i,idepmax
    real(double),dimension(:,:), allocatable   :: fp_buffer
    real(double)::deplamax,depla
    !-----------------------------------------------


    dxx(:,:)=xp_n(:,:,npath) - xp_n(:,:,1)
    open(unit=831,file='distimages')
    deplamax=0
    idepmax=0
    write(6,*)'im',im
    do i=1,im
       depla=1d8*sqrt(dxx(1,i)**2+dxx(2,i)**2+dxx(3,i)**2)
       write(831,*)i,depla
       if (depla.gt.deplamax) then 
          deplamax=depla
          idepmax=i
       end if
    end do
    write(6,*)'NEB: deplacement max entre configurations = ',deplamax,' pour l atome ',idepmax

    do iph=1,npath

       xp_n(:,:,iph) = xp_n(:,:,1) +                               &
            ( dxx(:,:) ) * dble(iph -1) / dble(npath-1)

       ielat_n   (:,iph) = ielat   (:) 
       iwmax_n   (:,iph) = iwmax   (:)
       ityp_n    (:,iph) = ityp    (:)
       xpp_n(:,:,iph)    = xp_n(:,:,iph)
       vp_n (:,:,iph)    = 0.d0
       ax_n (:,:,iph)    = xp_n(:,:,iph)
       fp_n (:,:,iph)    = 0.d0

    end do
   ! 	write(6,*)'toto'
!	write(6,*)cm
!	write(6,*)ityp
    masstot=SUM(cm(ityp(1:im)))
    if  (nebtype>=2) then
       write(*,'(" NEB: The kspring is in the eV/A^2                          :", f12.5)')  kspring
       kspring=kspring*angst**2/erg2eV 
       write(*,'(" NEB: The kspring is in the NDM internal units (copyright)  :", f12.5)')  kspring
    end if

    icontrainte(:)=1
    allocate (fp_buffer(3,imm))

    non_contr=0
    fp_buffer(:,:)=fp(:,:)
    fp(:,:)=1.d0
    if(lcontr) call contr (xp,vp,fp,ityp)
    do ic=1,im
       if(fp(1,ic).eq.0) then
          icontrainte(ic)=0
          non_contr=non_contr+1
       end if
    end do
    fp(:,:)=fp_buffer(:,:)

    write(*,'(" NEB: The number of atoms which are not included in the DRAG CONTRAINT :", i6)')  non_contr

    deallocate(fp_buffer)     


    return

  end subroutine init_neb

  subroutine build_s_path_neb(ityp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    integer  :: ityp(imm)
    integer  :: ip,ia
    real(double) :: e_i,e_i_p,e_i_m,dE_max, dE_min,norm,       &
         Rtemp_p,Rtemp_m, temp_m,temp_p
    real(double), dimension(3) :: tg_p(3),tg_m(3),Rtemp(3)
    real(double),dimension(3,npath)   :: rcm_loc
    real(double),dimension(3,imm)     :: dxx


    do ip=2,npath-1
       e_i  =enePATH(ip)
       e_i_p=enePATH(ip+1)
       e_i_m=enePATH(ip-1)
       temp_m=0.d0
       temp_p=0.d0
       do ia=1,im
          tg_p(:)=xp_n(:,ia,ip+1) - xp_n(:,ia,ip  ) 
          tg_m(:)=xp_n(:,ia,ip  ) - xp_n(:,ia,ip-1) 

          Rtemp_m= DOT_PRODUCT(tg_m,tg_m)
          Rtemp_p= DOT_PRODUCT(tg_p,tg_p)

          temp_m=temp_m+Rtemp_m
          temp_p=temp_p+Rtemp_p

          if  ((e_i_p > e_i ) .and. ( e_i > e_i_m)) then
	     ! 
             Rtemp(:)=tg_p(:)
             !
          else if ((e_i_p < e_i) .and. (e_i < e_i_m)) then
             !
             Rtemp(:)=tg_m(:)
             !
          else 
             !
             dE_max= max(dabs(e_i_p-e_i),dabs(e_i_m-e_i))
             dE_min= min(dabs(e_i_p-e_i),dabs(e_i_m-e_i))
             ! 
             if (e_i_p.ge.e_i_m) then
                !
	        Rtemp(:) = tg_p(:)*dE_max + tg_m(:)*dE_min
             else 
	        Rtemp(:) = tg_p(:)*dE_min + tg_m(:)*dE_max		
         !
             end if
             !
          end if

          s_path(:,ia,ip)=Rtemp(:)  
          !
       end do     !loop over ia
       !
       norms(ip)=SUM(s_path(:,:,ip)**2)
       force_neb(:,:,ip)=kspring*(dsqrt(temp_p)-dsqrt(temp_m))	
       !
    end do    ! loop over ip


    return



  end subroutine build_s_path_neb


  subroutine build_s_path_drag (ityp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: ityp(imm)
    integer  :: iph,ia,ic
    real(double)  :: dxx(3,imm),rcm_loc(3),masstot_contr


    dxx(:,:)=xp_n(:,:,npath) - xp_n(:,:,1)



    s_path(:,:,:)=0.d0
    rcm_loc(:)=0.d0
    masstot_contr=0.d0


    do ic=1,3
       do ia=1,im
          if (icontrainte(ia).eq.1) then
             rcm_loc(ic) = rcm_loc(ic) + dxx(ic,ia)
             masstot_contr=masstot_contr+cm(ityp(ia))
          end if
       end do
    end do

    do iph=2,npath-1
       do ia=1,im
          if (icontrainte(ia).eq.1) then
             s_path(1,ia,iph)= dxx(1,ia) -  rcm_loc(1)*cm(ityp(ia))/masstot_contr 
             s_path(2,ia,iph)= dxx(2,ia) -  rcm_loc(2)*cm(ityp(ia))/masstot_contr
             s_path(3,ia,iph)= dxx(3,ia) -  rcm_loc(3)*cm(ityp(ia))/masstot_contr
          end if
       end do
       norms(iph) = SUM(s_path(:,:,iph)**2)
    end do

    return

  end subroutine build_s_path_drag



  subroutine find_relax()

    implicit none
    integer   ::  ia
    real(double)   :: deltaR

    irelax(:)=0
    deltaRmax=deltaRmax**2
    do ia=1,im
       deltaR=DOT_PRODUCT(xp_n(:,ia,1)-xp_n(:,ia,npath),xp_n(:,ia,1)-xp_n(:,ia,npath))

       deltaR=deltaR*angst*angst

       if (deltaR>deltaRmax) then
          irelax(ia)=1
          write(*,*) 'NEB:    relaxation de l atome no ',ia,' deltaR= ',sqrt(deltaR)
       end if
    end do

    return
  end subroutine find_relax

  subroutine force_projection(ipath,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: ielat(imm)
    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    real(double)  :: fp(3,imm)      
    !-----------------------------------------------
    integer       :: ia,ipath
    real(double)  :: lbd,rcm_loc(3)
    !-----------------------------------------------

    lbd=0
    do ia=1,im

       if (irelax(ia)==1) then
          lbd = lbd + DOT_PRODUCT(s_path(:,ia,ipath),fp(:,ia))
       end if
    end do


    do ia=1,im
       if (irelax(ia)==1) then
          fp(:,ia) =fp(:,ia) - s_path(:,ia,ipath)*lbd/norms(ipath)
       end if
    end do

    if(lcontr) call contr (xp,vp,fp,ityp)      

    return

  end subroutine force_projection

  subroutine force_projection_neb(ipath,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: ielat(imm)
    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    real(double)  :: fp(3,imm)      
    !-----------------------------------------------
    integer       :: ia,ipath
    real(double)  :: lbd
    !-----------------------------------------------

    lbd=0
    do ia=1,im      
       if (irelax(ia)==1) then
          lbd = lbd + DOT_PRODUCT(s_path(:,ia,ipath),fp(:,ia))
       end if
    end do
    do ia=1,im

       if (irelax(ia)==1) then
          !
          fp(:,ia) =fp(:,ia) -  s_path(:,ia,ipath)*lbd/norms(ipath)  &
               + force_neb(:,ia,ipath)*s_path(:,ia,ipath) / dsqrt(norms(ipath))
          ! 	
       end if

    end do

    if(lcontr) call contr (xp,vp,fp,ityp)      

    return

  end subroutine force_projection_neb




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine configNEB(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
    USE T_kind_param_m, ONLY:  double
    use gen_com_m

    integer  :: ielat(imm)
    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    real(double)  :: fp(3,imm)

    integer :: ip,lucin,icintype,typmax,i,typmin
    character :: extension*2
    character :: fnamneb*80
    call allocate_neb()
    if (lrestart) then
       do ip=1,npath
          !          open(unit=17, file='tampon', form='formatted', status='unknown')
          !          if (ip==1) then
          !          fnamneb='deb.cout'
          !          goto 1
          !       end if
          !          if (ip==npath) then
          !          fnam='fin.cout'
          !          goto 1
          !       end if
         if (ip.ge.100) then
           write (6,*) 'ip >99 stop'
           stop 
         end if
         write(extension,'(i2.2)') ip

101       format(a1)
201       format(a2)
200       format(i2)

          fnamneb=fnam(1:lenfnam)//'.coutposition.'//extension
1         continue
          write(6,*)'image = ',fnam
          lucin = 93
          open(unit=lucin, file=fnamneb, form='unformatted', status='unknown')

          read (lucin) icintype
          if (icintype>=2) then
             read (lucin) at
             call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
             do ic = 1, 3
                normat(ic) = 0
                normat(ic) = normat(ic)+sum(at(:,ic)**2)
                normat(ic) = sqrt(normat(ic))
                zl(ic) = normat(ic)
             end do
             !                    write(6,*)at
             zls2 = zl/2.0

          else
             read (lucin) zl                      !size of the box
             if (rang==0) write (6, *) 'zl ', zl
             at(1,1)=zl(1)
             at(2,2)=zl(2)
             at(3,3)=zl(3)
             at(1,2)=zero
             at(1,3)=zero
             at(2,1)=zero
             at(2,3)=zero
             at(3,1)=zero
             at(3,2)=zero
             call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
             zls2 = zl/2.0

          endif

          read (lucin) im                         !number of atoms in the box
          if (im>imm) then
             if(rang==0)                    write (6, *) 'im > imM', im, imm
             stop
          endif


          read (lucin) ityp                       !types
          na(:ntyp) = 0
          typmax = -2
          typmin = 100

          typmax = max(maxval(ityp(:im)),typmax)

          i = 1
          if (im>0) then
             typmin = min(minval(ityp(:im)),typmin)
             i = im+1
          endif

          if (typmax>ntyp.or.typmin<1) then
             if(rang==0)            write (6, *) 'wrong ityp(', i, ')= ', ityp(i)
             stop
          endif

          do i=1,im
             na(ityp(i))=na(ityp(i))+1
          enddo
          read (lucin) xp
          write(6,*)'xp NEB image',ip

          !        xp_n (:,:,ip)    = xp (:,:)
          !        ityp_n    (:,ip) = ityp    (:)
          !        xpp_n(:,:,ip)    = 0.0
          !        vp_n (:,:,ip)    = 0.0
          !        ax_n (:,:,ip)    = 0.0
          !        fp_n (:,:,ip)    = 0.0

          ielat_n   (:,ip) = ielat   (:) 
          iwmax_n   (:,ip) = iwmax   (:)
          ityp_n    (:,ip) = ityp    (:)
          xp_n (:,:,ip)    = xp (:,:)
          xpp_n(:,:,ip)    = xpp(:,:)
          vp_n (:,:,ip)    = vp (:,:)
          ax_n (:,:,ip)    = ax (:,:)
          fp_n (:,:,ip)    = fp (:,:)
          close(lucin)
       end do
    else ! pas restart

       fnamneb=fnam
       fnam(1:lenfnam+4)='deb_'//fnamneb(1:lenfnam)
       write(6,*)'FNAM ',fnam
       lenfnam=lenfnam+4
       call config 
       fnam=fnamneb
       lenfnam=lenfnam-4
       call sauveposition(1)
       call into_path(1,1,xp, xpp, vp, ax, fp, ielat, iwmax,ityp)



       fnam(1:lenfnam+4)='fin_'//fnam(1:lenfnam)   
       lenfnam=lenfnam+4
       !       write(6,*)'FNAM ',fnam
       call config 
       fnam=fnamneb
       lenfnam=lenfnam-4
       !       write(6,*)'FNAM ',fnam

       call sauveposition(npath)
       call into_path(npath,1,xp, xpp, vp, ax, fp, ielat, iwmax,ityp)

       !       fnam=fnamneb
       !       lenfnam=lenfnam-4
       !       write(6,*)'FNAM ',fnam
       !---inNEB ATENTIE IN DETENTIE ASTA ESTE BUNA NUMAI PENTRU xp TOATE CELELALTE SUNT
       !---inNEB MINUNATE INSA NU EXISTA !!!!!!
       !---inNEB TREBUIE CA init SA CURGA PANA LIN PANA LA CAPAT  	  



    end if

  end subroutine configneb
  
  
  subroutine bruit_neb 
  implicit none
   integer    :: ia, ip
   integer, dimension(2) :: iseedt
   real(double)  :: zr1,zr2,zr3,zr4,totalbruit

  call system_clock (iseed)
   iseedt(1)=iseed
  call random_seed(iseedt(1))
 totalbruit=0.d0
 bruitneb(1:3,1:im,1:npath)=0.d0
 do ip=2,npath-1
  do ia=1,im
   call random_number(zr1)
   call random_number(zr2)
   call random_number(zr3)
   call random_number(zr4)
    if(zr1.eq.0.d0) zr1=0.000000001d0
    if(zr2.eq.0.d0) zr2=0.000000001d0
    if(zr3.eq.0.d0) zr3=0.000000001d0
    if(zr4.eq.0.d0) zr4=0.000000001d0
    
    bruitneb(1,ia,ip)=sqrt((-log(zr1)))*cos(2.0*pi*zr3)
    bruitneb(2,ia,ip)=sqrt((-log(zr1)))*sin(2.0*pi*zr3)
    bruitneb(3,ia,ip)=sqrt((-log(zr2)))*cos(2.0*pi*zr4)
    totalbruit=totalbruit + bruitneb(1,ia,ip)**2 + bruitneb(2,ia,ip)**2 + bruitneb(3,ia,ip)**2
    
   end do
  end do
  
  write(*,*) 'ISEED for neb, NORM of the noise ',iseed, neb_noise_scale, totalbruit
  bruitneb(1:3,1:im,2:npath-1) = bruitneb(1:3,1:im,2:npath-1) * neb_noise_scale * xp_n(1:3,1:im,2:npath-1) / (sqrt(totalbruit))
  xp_n(1:3,1:im,1:npath) = xp_n(1:3,1:im,1:npath) + bruitneb(1:3,1:im,1:npath)
  !debug write(*,*) xp_n(1,5,4), bruitneb(1,5,4)
  !debug stop
  
  end subroutine bruit_neb


end module neb_module
