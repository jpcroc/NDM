subroutine neb(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  !-----------------------------------------------
  use neb_module
  use posana
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
  integer  :: iph,i_dir_path
  integer  :: ielat(imm)
  integer  :: iwmax(imm)
  integer  :: ityp(imm)
  real(double)  :: xp(3,imm)
  real(double)  :: xpp(3,imm)
  real(double)  :: vp(3,imm)
  real(double)  :: ax(3,imm)
  real(double)  :: fp(3,imm)
  integer :: ineb,ii,it_neb_inter
  real(double)  :: a_local,forneb

  open(unit=55,file='image_col_relax.out')
  open(unit=56,file='react_col_relax.out')


  call init_neb(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

  if (nebrelaxation==1) then
     write(6,*)'NEB: nebrelaxation == 1'
     write(6,*)'NEB: relaxation NEB que pour les atomes' 
     write(6,*)'TUTU',im
     call find_relax()
  else
     irelax(:)=1
     write(6,*)'NEB: nebrelaxation /= 1'
     write(6,*)'NEB: relaxation NEB pour TOUS les atomes'
  end if
  
  if (neb_noise.eq.1) then
     write(6,*)'NEB: We apply a random noise on the atoms '
    call bruit_neb
  end if



  do ii=1,npath
     it=1
     call into_path(ii,2,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
     call scalebox           (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)       
     call calfo 
     call analyse  
     call neb_controle    (ii,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)       
     enePATH(ii)=potist
     enePATHev(ii)=potist*erg2ev
     write(*,'(i5,3(g20.8,1x))') ii, enePATHev(ii),enePATHev(ii)-enePATHev(1)
     call into_path(ii,1,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

  end do



  select case (nebtype)


  case (1)
     write(6,*) 'NEB: !!!!-------this is DRAG----------!!!!!!'
     call build_s_path_drag(ityp)

     do ii=2,npath-1
        call into_path(ii,2,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

        it=0
        dragtest=0
        do while (dragtest==0)
           it=it+1
           call scalebox           (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)       
           call calfo
           call force_projection(ii,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
           call trempe             (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
           call analyse  
           call neb_controle    (ii,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
        end do   ! end do for a while

        enePATH(ii)=potist
        enePATHev(ii)=potist*erg2ev
        write(*,*) 'NEB: THIS IS THE DRAG IMAGE=====================', ii
        write(*,*) 'NEB: THE NUMBER OF ITERATIONS===================', it
        write(*,*) 'NEB: THE ENERGY OF THIS IMAGE===================', enePATHev(ii)
        call into_path(ii,1,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
     end do      !end ii,npath


  case(2)	
     write(6,*)'NEB: -------this is NEB--V2-------'
     write(6,*)'NEB: The MAX steps in NEB        :',maxneb
     nebtest(:)=0
     do ineb=1,maxneb     ! Main loop for NEB
        !    
        it=0
        !   
        call build_s_path_neb(ityp)       
        ! 
        print'("NEB:===============================================")'
        print'("NEB:pas-neb image    force      force_NEB          energie      statut    energie/stable")'
        do ii=2,npath-1
           !
           call into_path      (ii,2,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
           !
           it_neb_inter=0  
           !
           do while (it_neb_inter<=5)   ! drag-ize me that 5 steps while we keep NEB "attraction"
              !
              it_neb_inter=it_neb_inter+1
              it=it_neb_inter
              call scalebox               (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)       
              call calfo
              call force_projection_neb(ii,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
              call trempe                 (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
              call analyse 
              call neb_controle        (ii,xp, xpp, vp, ax, fp, ielat, iwmax, ityp) 
              !
           end do

           ! Backup image ii
           if (itesauv.GT.0) then
                   if (mod(ineb,itesauv)==0) call sauveposition(ii)
           else if (itesauvposition.GT.0) then
                   if (mod(ineb,itesauvposition)==0) call sauveposition ( ii)
           endif

           !debug          print'("NEB: ",2i5,3E14.5,E20.10,i3)', ineb, ii,  formax,       &
           !debug	            formaxperp, formaxparl, potist*erg2eV,nebtest(ii)
           forneb = SQRT(MAXVAL(force_neb(1,:,ii)**2 + force_neb(2,:,ii)**2         &
                + force_neb(3,:,ii)**2))*erg2eV/angst 
           print'("NEB: ",2i5, 2E14.5,E20.10,i3,F15.8)', ineb, ii,  formax, forneb,       &
                potist*erg2eV,nebtest(ii),potist*erg2eV-enepathev(1)
           ! 
           enePATH(ii)=potist
           enePATHev(ii)=potist*erg2eV       
           call into_path       (ii,1, xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
           !
        end do      ! end ii,path

        if (SUM(nebtest(:))==(npath-2)) then
           print'("NEB:=== bye,bye sweety is finished=================")'
           print'("NEB:===============================================")'
           exit
        end if

     end do        !maxneb 

  case default

     write(6,*) 'NO NEB-DYNAMICS FOR THIS  nebtype = ', nebtype
     write(6,*) 'CHANGE nebtype AND TRY AGAIN.'
     write(6,*) 'nebtype=1 for DRAG'
     write(6,*) 'nebtype=2 for  NEB'
     write(6,*) 'STOP in the neb.f90'
     stop

  end select


  reaction_coord(1)=0
  reaction_coord(npath)=1
  a_local=SUM((xp_n(:,:,npath)-xp_n(:,:,1))**2) 
  ! 
  do ii=1,npath      
     call into_path      (ii,2,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
     call sauveposition(ii)      
     call rasmol(ii)
     if (iteanaposneb.gt.0) call anapos(ii)
     !	 
     reaction_coord(ii) = SUM((xp_n(:,:,ii)-xp_n(:,:,1))*(xp_n(:,:,npath)-xp_n(:,:,1)))/a_local 
     !
  end do

  print'("NEB:--IMAGE-----REACT-COORD------ENERGY------ENERGY-ENERGY(1)")'
  do ii=1,npath
     !
     write(*,'(i5,3(g20.8,1x))')  ii, reaction_coord(ii), enePATHev(ii),enePATHev(ii)-enePATHev(1)
     write(55,'(3(g20.8,1x))') dble(ii-1)/dble(npath-1), enePATHev(ii),enePATHev(ii)-enePATHev(1)
     write(56,'(3(g20.8,1x))') reaction_coord(ii), enePATHev(ii),enePATHev(ii)-enePATHev(1)
     !       
  end do


  print*,'-------------------------------------------'
  write(*,'("          Min .. et .. Max")')
  print '("E:",2g14.7)',minval(enePATHev), maxval(enePATHev)
  write(*,'("          1 .. et .. NPATH")')
  print '("B:",2g14.7)',enePATHev(1), enePATHev(npath)
  write(*,'("les diffs")')
  print*,'1-NPATH    :',enePATHev(1)-enePATHev(npath)
  print*,'-------------------------------------------'
  print*,'MAX-1      :',maxval(enePATHev)-enePATHev(1)
  print*,'MAX-NPATH  :',maxval(enePATHev)-enePATHev(npath)



end subroutine neb
