module neb_module
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:iseed,neb_noise_scale,npath,lrestart,npath,deltarmax,kspring,lpathfromgin,&
       &lrestart,nebtype, fnam,pi,rang,im_glob,lenfnam,rang,zero,lcontr,&
       &angst,lenfnam,angst,erg2ev,normat,ltabvois,nvois,fnamcout,igen,lprteat
  USE constrconf_mod,only:constr_2gin,gin2ndm
    use cryst_to_cart_mod,only:cryst_to_cart
  USE contrainte,only:contr
!  USE config_mod,only: config
  USE recips_mod,only: recips
  USE sauveposition_mod,only: sauveposition
  USE rasmol_mod,only: rasmol
  use var_pot,only:ntyp,na,ipotentiel,cm,rumax
  USE dynalloccell,only:deallocateall
  !-----------------------------------------------
  USE atomconfig,only:atom_config,atom_config_d
  USE cellconfig, only:cell_config,init_cel
  USE constrconf_mod,only : config2data
  USE read_conf, only:read_cin,read_gin
  use boxconfig,only: box_config,ndm2boxconfig,boxconfig2ndm,initbox
  USE setcell,only:setcellconf,setnox
  USE sauvegardeT_mod,only:sauvegardeT
  implicit none


  integer, save                                  :: dragtest
  integer,dimension(:),allocatable,save          :: nebtest,icontrainte !irelax,
!  integer,dimension(:,:),allocatable,save        :: ielat_n, iwmax_n, ityp_n
!  real(double),dimension(:,:,:),allocatable,save ::  xp_n, xpp_n, vp_n, fp_n
!  real(double),dimension(:,:),allocatable        :: fp_par,fp_perp
  real(double), dimension(:),allocatable, save   :: enePATH,enePATHev,norms,reaction_coord
  real(double), dimension(:,:,:), allocatable, save :: sigPATH  ! Stress tensor
  real(double), dimension(:,:,:),allocatable,save:: s_path,force_neb,bruitneb 
  real(double)                                   :: forctot,formax,formaxperp,formaxparl,masstot
  logical:: lvzeroneb
  type(atom_config_d),allocatable,save::atneb(:)
  type(cell_config),allocatable,save:: cellneb(:)
  type(box_config)::boxneb
contains
  

  
  subroutine allocate_neb(im,imm)
    implicit none
    integer::im,imm
 !   integer ip2
 !   allocate(atneb(npath)) !   do ip2=1,npath
 !      call atneb(ip2)%init_atom_config(im,ltabvois,nvois)
    !   end do
    integer ipath
    allocate(atneb(npath))
    allocate(cellneb(npath))

    do ipath=1,npath
       call atneb(ipath)%init(im,imm,ltabvois,nvois,lprteat=lprteat)
    end do

         allocate (icontrainte(imm),reaction_coord(npath))
    allocate  (enePATH(npath),enePATHev(npath),norms(npath),nebtest(npath))
    allocate  (sigPATH(3,3,npath))    ! Stress tensor for each image

    allocate  (s_path(3,imm,npath),force_neb(3,imm,npath))
    allocate  (bruitneb(3,imm,npath))

    return

  end subroutine allocate_neb

  ! **************************************************************
  subroutine into_path(iph, i_dir_path, &                         
       xp, xpp, vp,  fp, ielat, iwmax, ityp,num_at_glob,imm)
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
    integer  :: num_at_glob(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: fp(3,imm)
    integer::imm
    !-----------------------------------------------
    !-----------------------------------------------
    integer   :: iph, i_dir_path
    !-----------------------------------------------


    if (i_dir_path==1) then
       atneb(iph)%ielat(:) = ielat(:)
       if(allocated(atneb(iph)%iwmax))      atneb(iph)%iwmax(:) = iwmax(:)
       atneb(iph)%ityp(:) = ityp(:)
       atneb(iph)%xp(:,:) = xp(:,:)
       atneb(iph)%vp(:,:) = 0
       atneb(iph)%xpp(:,:) = xpp(:,:)
       atneb(iph)%fp(:,:) = fp(:,:) 
       atneb(iph)%num_at_glob(:) =num_at_glob   (:) 
    else 
       ielat   (:)      = atneb(iph)%ielat(:)
       if(allocated(atneb(iph)%iwmax))  iwmax   (:)      = atneb(iph)%iwmax(:)
       num_at_glob   (:)      = atneb(iph)%num_at_glob(:)
       ityp    (:)      = atneb(iph)%ityp(:)
       xp (1:3,:)        = atneb(iph)%xp(:,:)
       xpp(1:3,:)       =atneb(iph)%xpp(:,:)
       vp(1:3,:) = 0.d0
       fp (1:3,:)      = atneb(iph)%fp(:,:)
end if

    return

  end subroutine into_path



  subroutine init_neb(im,imm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    integer::imm,im
    real(double)  :: dxx(3,imm)
    !-----------------------------------------------
    integer :: iph,ic,non_contr,i,idepmax,itread
    real(double),dimension(:,:), allocatable   :: fp_buffer
    real(double)::deplamax,depla
    CHARACTER(len=80) :: ginFile
    CHARACTER(len=9) :: extension
    CHARACTER(len=89) :: fnamneb
    logical ::ok
    !-----------------------------------------------

    dxx(:,:)=atneb(npath)%xp(:,:)-atneb(1)%xp(:,:) 
!    dxx(:,:)=xp_n(:,:,npath) - xp_n(:,:,1)
    open(unit=831,file='distimages')
    deplamax=0
    idepmax=0
    if(rang==0)write(6,*)'im',atneb(1)%im
    do i=1,im
       depla=1d8*sqrt(dxx(1,i)**2+dxx(2,i)**2+dxx(3,i)**2)
       write(831,*)i,depla
       if (depla.gt.deplamax) then 
          deplamax=depla
          idepmax=i
       end if
    end do
    if (rang==0) write(6,*)'NEB: deplacement max entre configurations = ',deplamax,' pour l atome ',idepmax

    do iph=1,npath

       IF (lrestart) THEN
          write(extension,'(i9.9)') iph
          fnamneb=fnam(1:lenfnam)//'.coutposition.'//extension
          itread=1
          call read_cin(boxneb,itread,atneb(iph),imm,fnamneb) !0=at seulement; 1=complet; 2 = at, xp et num_at_glob seulement , 3 trié par num_at_buff
       ELSE IF (lPathFromGin) THEN
               ! read initial path in gin files *.1.gin, *.2.gin, ...

          WRITE(ginFile, '(2a,i0,a)') Trim(fnam), '.', iph, '.gin'
          INQUIRE(file=ginFile, exist=ok)
          IF ((.NOT.ok).AND.(iph.LE.999999999)) THEN
             WRITE(ginFile, '(2a,i9.9,a)') Trim(fnam), '.', iph, '.gin'
             INQUIRE(file=ginFile, exist=ok)
          END IF
          IF ( (.NOT.ok).AND.(iph.LE.99999999)) THEN
             WRITE(ginFile, '(2a,i8.8,a)') Trim(fnam), '.', iph, '.gin'
             INQUIRE(file=ginFile, exist=ok)
          END IF
          IF ( (.NOT.ok).AND.(iph.LE.9999999)) THEN
             WRITE(ginFile, '(2a,i7.7,a)') Trim(fnam), '.', iph, '.gin'
             INQUIRE(file=ginFile, exist=ok)
          END IF
          IF ( (.NOT.ok).AND.(iph.LE.999999)) THEN
             WRITE(ginFile, '(2a,i6.6,a)') Trim(fnam), '.', iph, '.gin'
             INQUIRE(file=ginFile, exist=ok)
          END IF
          IF ( (.NOT.ok).AND.(iph.LE.99999)) THEN
             WRITE(ginFile, '(2a,i5.5,a)') Trim(fnam), '.', iph, '.gin'
             INQUIRE(file=ginFile, exist=ok)
          END IF
          IF ( (.NOT.ok).AND.(iph.LE.9999)) THEN
             WRITE(ginFile, '(2a,i4.4,a)') Trim(fnam), '.', iph, '.gin'
             INQUIRE(file=ginFile, exist=ok)
          END IF
          IF ( (.NOT.ok).AND.(iph.LE.999)) THEN
             WRITE(ginFile, '(2a,i3.3,a)') Trim(fnam), '.', iph, '.gin'
             INQUIRE(file=ginFile, exist=ok)
          END IF
          IF ( (.NOT.ok).AND.(iph.LE.99)) THEN
             WRITE(ginFile, '(2a,i2.2,a)') Trim(fnam), '.', iph, '.gin'
             INQUIRE(file=ginFile, exist=ok)
          END IF
          IF ( (.NOT.ok).AND. (iph.EQ.1) ) THEN
             WRITE(ginFile, '(3a)') 'deb_', Trim(fnam), '.gin'
             INQUIRE(file=ginFile, exist=ok)
          END IF
          IF ( (.NOT.ok).AND. (iph.EQ.nPath) ) THEN
             WRITE(ginFile, '(3a)') 'fin_', Trim(fnam), '.gin'
             INQUIRE(file=ginFile, exist=ok)
          END IF
          
          IF (ok ) THEN
             ! Load NEB image ip in file *.<ip>.gin
          if(rang==0)write(6,*)'FNAMneb  ',iph,ginfile
          call gin2ndm(atneb(iph),cellneb(iph),boxneb,ginfile,im_glob,rumax)
            do i=1,im
               atneb(iph)%num_at_glob(i)=i
            end do
         ELSE
            WRITE(0,'(a,i0)') 'Does not manage to find a backup file for image ', iph
            WRITE(0,'(3a)') 'File ', Trim(ginFile), ' does not exist'
            STOP '< Load_NEB_Image_Gin >'
         END IF
       ELSE

          atneb(iph)%xp(:,:)=atneb(1)%xp(:,:)+dxx(:,:)*dble(iph -1) / dble(npath-1)
          atneb(iph)%ityp(:)=atneb(1)%ityp(:)
          atneb(iph)%num_at_glob(:)=atneb(1)%num_at_glob(:)
          !CRC RESTE ITYP...
       END IF
       atneb(iph)%xpp(:,:)= atneb(iph)%xp(:,:)
    end do

    masstot=SUM(cm(atneb(1)%ityp(1:im)))
    if  (nebtype>=2) then
       if (rang==0) write(*,'(" NEB: The kspring is in the eV/A^2                          :", f12.5)')  kspring
       kspring=kspring*angst**2/erg2eV 
       if (rang==0) write(*,'(" NEB: The kspring is in the NDM internal units (cgs)  :", f12.5)')  kspring
    end if

    icontrainte(:)=1
    allocate (fp_buffer(3,imm))

    deallocate(fp_buffer)     

    return

  end subroutine init_neb


  !-------------------------------------------------------

  subroutine build_s_path_neb(im,imm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    integer  :: ip,ia,im,imm
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
          
          tg_p(:)=atneb(ip+1)%xp(:,ia) - atneb(ip)%xp(:,ia  )
          tg_m(:)=atneb(ip)%xp(:,ia) - atneb(ip-1)%xp(:,ia  ) 
!          tg_p(:)=xp_n(:,ia,ip+1) - xp_n(:,ia,ip  ) 
!          tg_m(:)=xp_n(:,ia,ip  ) - xp_n(:,ia,ip-1) 

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


  subroutine build_s_path_drag(im,imm)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: iph,ia,ic,imm,im
    real(double)  :: dxx(3,imm),rcm_loc(3)


    dxx(:,:)=atneb(npath)%xp(:,:) - atneb(1)%xp(:,:)
!    dxx(:,:)=xp_n(:,:,npath) - xp_n(:,:,1)


    s_path(:,:,:)=0.d0
    rcm_loc(:)=0.d0


    do ic=1,3
       do ia=1,im
          if (icontrainte(ia).eq.1) then
             rcm_loc(ic) = rcm_loc(ic) + dxx(ic,ia)
          end if
       end do
    end do

    do iph=2,npath-1
       do ia=1,im
          if (icontrainte(ia).eq.1) then
             s_path(1,ia,iph)= dxx(1,ia) -  rcm_loc(1)*cm(atneb(1)%ityp(ia))/masstot 
             s_path(2,ia,iph)= dxx(2,ia) -  rcm_loc(2)*cm(atneb(1)%ityp(ia))/masstot
             s_path(3,ia,iph)= dxx(3,ia) -  rcm_loc(3)*cm(atneb(1)%ityp(ia))/masstot
          end if
       end do
       norms(iph) = SUM(s_path(:,:,iph)**2)
    end do

    return

  end subroutine build_s_path_drag



  subroutine find_relax(im)

    implicit none
    integer   ::  ia,ipath,im
    real(double)   :: deltaR

!    atneb(1:)%lgul(:)=.false.
    !    irelax(:)=0
    deltaRmax=deltaRmax**2
    do ia=1,im
       deltaR=DOT_PRODUCT(atneb(1)%xp(:,ia)-atneb(npath)%xp(:,ia),atneb(1)%xp(:,ia)-atneb(npath)%xp(:,ia))


       deltaR=deltaR*angst*angst
       write(6,*)ia,deltar,deltarmax
       if (deltaR>deltaRmax) then
          atneb(1)%lgul(ia)=.true.
!          irelax(ia)=1
          if (rang==0) write(*,*) 'NEB:    relaxation de l atome no ',ia,' deltaR= ',sqrt(deltaR)
       end if
    end do
    do ipath=2,npath
       atneb(ipath)%lgul(:)=atneb(ipath)%lgul(:)
    end do
    
    return
  end subroutine find_relax
!CRC a gérer dans les appels
  subroutine force_projection(ipath,xp,  vp,  fp,  ityp,imm,im)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    integer::imm,im,ipath
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
!    integer  :: ielat(imm)
!    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
!    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: fp(3,imm)      
    !-----------------------------------------------
    integer       :: ia
    real(double)  :: lbd,rcm_loc(3)
    !-----------------------------------------------

    lbd=0
    do ia=1,im

       if (atneb(ipath)%lgul(ia)) then
          lbd = lbd + DOT_PRODUCT(s_path(:,ia,ipath),fp(:,ia))
       end if
    end do


    do ia=1,im
       if (atneb(ipath)%lgul(ia)) then
          fp(:,ia) =fp(:,ia) - s_path(:,ia,ipath)*lbd/norms(ipath)
       end if
    end do

    if(lcontr) call contr (xp,vp,fp,ityp)      

    return

  end subroutine force_projection

  subroutine force_projection_neb(ipath,xp, vp,  fp,  ityp,imm,im)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer::imm,im
    integer  :: ielat(imm)
    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: fp(3,imm)      
    !-----------------------------------------------
    integer       :: ia,ipath
    real(double)  :: lbd
    !-----------------------------------------------

    lbd=0
    do ia=1,im      
       if (atneb(ipath)%lgul(ia)) then
          lbd = lbd + DOT_PRODUCT(s_path(:,ia,ipath),fp(:,ia))
       end if
    end do
    do ia=1,im

       if (atneb(ipath)%lgul(ia)) then
          !
          fp(:,ia) =fp(:,ia) -  s_path(:,ia,ipath)*lbd/norms(ipath)  &
               + force_neb(:,ia,ipath)*s_path(:,ia,ipath) / dsqrt(norms(ipath))
          ! 	
       end if

    end do

    if(lcontr) call contr (xp,vp,fp,ityp)      

    return

  end subroutine force_projection_neb


  subroutine constrconfNEB !(xp, xpp, vp,  fp, ielat, iwmax, ityp)
    USE T_kind_param_m, ONLY:  double
    !    USE gen_com_m, ONLY:im,imm,nvois,
!    USE tab_imm_m,only:xp,xpp,vp,fp,ielat,iwmax,ityp,num_at_glob
    USE gen_com_m, ONLY:imm

    integer :: ic, ip,lucin,icintype,typmax,i,typmin,itread,fmt_cin,formatsauv,iti,lat(3)
    character :: extension*9
    character :: fnamneb*80
    real(double)::atg(3,3)
    type(atom_config)::atrgin
    type(box_config)::boxrgin
    call allocate_neb(0,imm)

    if (igen==1) then 
       itread=1;fmt_cin=2
       if (lrestart) then
          do ip=1, npath, npath-1 !CRC ne lit que deux images ??
             write(extension,'(i9.9)') ip
             fnamneb=fnam(1:lenfnam)//'.coutposition.'//extension
             if (rang==0) write(6,'(2a)')'image = ',fnamneb
             call read_cin(boxneb,itread,atneb(ip),imm,fnamneb,lrestart,fmt_cin)
             atneb(ip)%ielat(:)=0 !ielat(:)
             atneb(ip)%iwmax(:)=0 !iwmax(:)
             atneb(ip)%fp(:,:)= 0 !fp(:,:)
             atneb(ip)%vp(:,:)=0 !vp(:,:)
             atneb(ip)%xpp(:,:)=atneb(ip)%xp(:,:) !xpp(:,:)
             call setnox(boxneb,cellneb(ip),rumax)
             !          CALL fin allocation CELL et FIN DIVID
             do iti=1,ntyp
                na(iti)=count(atneb(ip)%ityp(1:atneb(ip)%im).eq.iti)
             end do
             close(lucin)
          end do
       else ! pas restart
          fnamneb='deb_'//fnam(1:lenfnam)//'.cin'
          if(rang==0)write(6,*)'FNAMneb 1 ',fnamneb
          call read_cin(boxneb,itread,atneb(1),imm,fnamneb,lrestart,fmt_cin)
          atneb(:)%im=atneb(1)%im
          atneb(1)%xpp=atneb(1)%xp
          atneb(1)%ielat=0
          atneb(1)%fp(:,:)= 0 !fp(:,:)
          atneb(1)%vp(:,:)=0 !vp(:,:)
          
          if (atneb(1)%ltabvois) then
             atneb(1)%iwmax=0 !iwmax(:)
             atneb(1)%indi=0
          end if
          call setnox(boxneb,cellneb(1),rumax)
          call setcellconf(cellneb(1),atneb(1),boxneb,im_glob,rumax)
          if (rang==0)then
             formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'neb.1.cout.'
             call sauvegardeT(atneb(1),cellneb(1),boxneb,formatsauv,fnamcout)
             call rasmol(atneb(1),boxneb,1)
          endif
          
          fnamneb='fin_'//fnam(1:lenfnam)//'.cin'
          if(rang==0)write(6,*)'FNAMneb npath ',fnamneb
          call read_cin(boxneb,itread,atneb(npath),imm,fnamneb,lrestart,fmt_cin)
          atneb(npath)%xpp=atneb(npath)%xp
          atneb(npath)%ielat=0
          atneb(npath)%fp(:,:)= 0 !fp(:,:)
          atneb(npath)%vp(:,:)=0 !vp(:,:)
          !       atneb(npath)%xpp(:,:)=0 !xpp(:,:)
          if (atneb(npath)%ltabvois) then
             atneb(npath)%iwmax=0 !iwmax(:)
             atneb(npath)%indi=0
          end if
          call setnox(boxneb,cellneb(npath),rumax)
          call setcellconf(cellneb(npath),atneb(npath),boxneb,im_glob,rumax)
          if (rang==0)then
             formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'neb.1.cout.'
             call sauvegardeT(atneb(npath),cellneb(npath),boxneb,formatsauv,fnamcout)
             call rasmol(atneb(npath),boxneb,npath)
          endif
          
       end if
    else
       fnamneb='deb_'//fnam(1:lenfnam)//'.gin'
       if(rang==0)write(6,*)'FNAMneb 1 ',fnamneb
       call gin2ndm(atneb(1),cellneb(1),boxneb,fnamneb,im_glob,rumax)
!!$       call read_gin(boxrgin,atrgin,fnamneb,lat)
!!$       do ic=1,3
!!$          atg(:,ic)=boxrgin%at(:,ic)*lat(ic)
!!$       end do
!!$       call initbox(boxneb,atg)
!!$
!!$       call setnox(boxneb,cellneb(1),rumax)
!!$       if (rang==0) then
!!$          write (6, '(A,D15.8,A,D15.8,A)') 'volume=', boxneb%volu,' cm3 ',boxneb%volu*1d24,' Ang3'
!!$       end if
!!$       
!!$       call constr_2gin (atneb(1),boxneb,cellneb(1),atrgin,boxrgin,lat,im_glob)
!!$       call cryst_to_cart (atneb(1)%imm, atneb(1)%xp, boxneb%at, 1)

       atneb(:)%im=atneb(1)%im
!       call atneb(1)%print(iwr=0)
          atneb(1)%xpp=atneb(1)%xp
          atneb(1)%ielat=0
          atneb(1)%fp(:,:)= 0 !fp(:,:)
          atneb(1)%vp(:,:)=0 !vp(:,:)
          if (atneb(1)%ltabvois) then
             atneb(1)%iwmax=0 !iwmax(:)
             atneb(1)%indi=0
          end if
!          call setcellconf(cellneb(1),atneb(1),boxneb,im_glob,rumax)
!          call cellneb(1)%print
!          call boxneb%print
          if (rang==0)then
             formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'neb.npath.cout.'
             call sauvegardeT(atneb(1),cellneb(1),boxneb,formatsauv,fnamcout)
             call rasmol(atneb(1),boxneb,1)
          endif

          
       fnamneb='fin_'//fnam(1:lenfnam)//'.gin'
       if(rang==0)write(6,*)'FNAMneb npath ',fnamneb
       call gin2ndm(atneb(npath),cellneb(npath),boxneb,fnamneb,im_glob,rumax)

!!$       call read_gin(boxrgin,atrgin,fnamneb,lat)
!!$       do ic=1,3
!!$          atg(:,ic)=boxrgin%at(:,ic)*lat(ic)
!!$       end do
!!$       call initbox(boxneb,atg)
!!$       call setnox(boxneb,cellneb(npath),rumax)
!!$!       if (rang==0) then
!!$!          write (6, '(A,D15.8,A,D15.8,A)') 'volume=', boxrcf%volu,' cm3 ',boxrcf%volu*1d24,' Ang3'
!!$!       end if
!!$       
!!$       call constr_2gin (atneb(npath),boxneb,cellneb(npath),atrgin,boxrgin,lat,im_glob)
!!$       call cryst_to_cart (atneb(npath)%imm, atneb(npath)%xp, boxneb%at, 1)
!       call atneb(npath)%print(iwr=0)
       atneb(:)%im=atneb(npath)%im
          atneb(npath)%xpp=atneb(npath)%xp
          atneb(npath)%ielat=0
          atneb(npath)%fp(:,:)= 0 !fp(:,:)
          atneb(npath)%vp(:,:)=0 !vp(:,:)
          if (atneb(npath)%ltabvois) then
             atneb(npath)%iwmax=0 !iwmax(:)
             atneb(npath)%indi=0
          end if
!          call setcellconf(cellneb(npath),atneb(npath),boxneb,im_glob,rumax)
          if (rang==0)then
             formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'neb.npath.cout.'
             call sauvegardeT(atneb(npath),cellneb(npath),boxneb,formatsauv,fnamcout)
             call rasmol(atneb(npath),boxneb,1)
          endif

       call atrgin%dealloc
       end if

    !#ifdef LAMMPS_VERSION

    if((ipotentiel==-10).or.(ipotentiel==-11)) then
       if(rang==0)write(6,*)'write configuration to conf.lmp'
       call config2data (atneb(1)%imm,atneb(1)%im,atneb(1)%xp,atneb(1)%ityp,boxneb%at,ntyp)
    end if
    !#endif     



  end subroutine constrconfNEB
  
  
  subroutine bruit_neb (im)
  implicit none
   integer    :: ia, ip,im
   integer, dimension(2) :: iseedt
   real(double)  :: zr1,zr2,zr3,zr4,totalbruit

  call system_clock (iseed)
   iseedt(1)=iseed
  call random_seed(iseedt(1))
 totalbruit=0.d0
 bruitneb(1:3,1:atneb(1)%im,1:npath)=0.d0
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
  
  if (rang==0) write(*,*) 'ISEED for neb, NORM of the noise ',iseed, neb_noise_scale, totalbruit
  do ip=2,npath-1
     bruitneb(1:3,1:im,ip) = bruitneb(1:3,1:im,ip) * neb_noise_scale * atneb(ip)%xp(1:3,1:im) / (sqrt(totalbruit))
  end do
  do ip=1,npath
     atneb(ip)%xp(1:3,1:im) = atneb(ip)%xp(1:3,1:im) + bruitneb(1:3,1:im,ip)
  end do
  !debug write(*,*) xp_n(1,5,4), bruitneb(1,5,4)
  !debug stop
  
  end subroutine bruit_neb

end module neb_module
