module neb_mod
  USE calfo_mod
  USE analyse_mod
  USE trempe_mod
  USE neb_controle_mod
  USE scalebox_mod
  USE sauveforce_mod
  USE gen_com_m, ONLY:iteanaposneb,itesauvforce,itesauvposition,lfire,maxneb,neb_noise,nebrelaxation
  implicit none 
contains
  subroutine neb ! (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    !-----------------------------------------------
    USE neb_module
    USE posana
    USE FireModule
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
    !integer  :: iph,i_dir_path
    !integer  :: ielat(imm)
    !integer  :: iwmax(imm)
    !integer  :: ityp(imm)
    !real(double)  :: xp(3,imm)
    !real(double)  :: xpp(3,imm)
    !real(double)  :: vp(3,imm)
    !real(double)  :: ax(3,imm)
    !real(double)  :: fp(3,imm)
    integer :: ineb,ii,it_neb_inter
    real(double)  :: a_local,forneb

    ! Variables for Fire quench algorithm
    REAL(double), dimension(:), allocatable :: fire_dt, fire_alph
    INTEGER, dimension(:), allocatable :: fire_nstep

    if(lPkbar) then
       unitP=1.0d-9
       cunitP='kbar'
    else
       unitP=1.0
       cunitP='d/cm2'
    endif

    open(unit=55,file='image_col_relax.out')
    WRITE(55,'(a)')  '#  1: (i-1)/(nPath-1)'
    WRITE(55,'(a)')  '#  2: energy E(i) (eV)'
    WRITE(55,'(a)')  '#  3: energy difference E(i)-E(1) (eV)'
    WRITE(55,'(3a)') '#  4: stess s(1,1) (', cUnitP, ')'
    WRITE(55,'(a)')  '#  5:       s(2,2)'
    WRITE(55,'(a)')  '#  6:       s(3,3)'
    WRITE(55,'(a)')  '#  7:       s(2,3)'
    WRITE(55,'(a)')  '#  8:       s(1,3)'
    WRITE(55,'(a)')  '#  9:       s(1,2)'
    open(unit=56,file='react_col_relax.out')
    WRITE(56,'(a)')  '#  1: reaction coordinate z(i)'
    WRITE(56,'(a)')  '#  2: energy E(i) (eV)'
    WRITE(56,'(a)')  '#  3: energy difference E(i)-E(1) (eV)'
    WRITE(56,'(3a)') '#  4: stess s(1,1) (', cUnitP, ')'
    WRITE(56,'(a)')  '#  5:       s(2,2)'
    WRITE(56,'(a)')  '#  6:       s(3,3)'
    WRITE(56,'(a)')  '#  7:       s(2,3)'
    WRITE(56,'(a)')  '#  8:       s(1,3)'
    WRITE(56,'(a)')  '#  9:       s(1,2)'


    call init_neb(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

    ! Initialization of fire quench algorithm
    IF (lFire) THEN
       ALLOCATE(fire_dt(1:npath))
       ALLOCATE(fire_nstep(1:npath))
       ALLOCATE(fire_alph(1:npath))
       do ii=1,npath
          CALL init_trempe_fire(fire_dt(ii), fire_nstep(ii), fire_alph(ii))
       END DO
    END IF


    if (nebrelaxation==1) then
       if (rang==0) write(6,*)'NEB: nebrelaxation == 1'
       if (rang==0) write(6,*)'NEB: relaxation NEB que pour les atomes' 
       call find_relax()
    else
       irelax(:)=1
       if (rang==0) write(6,*)'NEB: nebrelaxation /= 1'
       if (rang==0) write(6,*)'NEB: relaxation NEB pour TOUS les atomes'
    end if

    if (neb_noise.eq.1) then
       if (rang==0) write(6,*)'NEB: We apply a random noise on the atoms '
       call bruit_neb
    end if


    !call calfo
    do ii=1,npath
       it=1
       call into_path(ii,2,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
       call scalebox           (xp, xpp, vp, ax, fp, ielat, iwmax, ityp)       
       !write(*,*) 'inside NEB debug1',ii, xp(1,1)
       call calfo 
       !write(*,*) 'inside NEB debug2',ii, xp(1,1)
       call analyse  
       call neb_controle(ii)  !  (ii,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)       
       enePATH(ii)=potist
       enePATHev(ii)=potist*erg2ev
       sigPATH(:,:,ii) = sigtot(:,:)      ! Contrainte
       if (rang==0) write(*,'(i5,3(g20.8,1x))') ii, enePATHev(ii),enePATHev(ii)-enePATHev(1)
       call into_path(ii,1,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

    end do


    select case (nebtype)


    case (1)
       if (rang==0) write(6,*) 'NEB: !!!!-------this is DRAG----------!!!!!!'
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
             IF (lFire) THEN
                call trempe_fire(xp, xpp, vp, ax, fp, ielat, iwmax, ityp, &
                     fire_dt(ii), fire_nstep(ii), fire_alph(ii))
             ELSE
                call trempe(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
             ENDIF
             call analyse  
             call neb_controle(ii) !   (ii,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
          end do   ! end do for a while

          enePATH(ii)=potist
          enePATHev(ii)=potist*erg2ev
          sigPATH(:,:,ii) = sigtot(:,:)      ! Contrainte
          if (rang==0) write(*,*) 'NEB: THIS IS THE DRAG IMAGE=====================', ii
          if (rang==0) write(*,*) 'NEB: THE NUMBER OF ITERATIONS===================', it
          if (rang==0) write(*,*) 'NEB: THE ENERGY OF THIS IMAGE===================', enePATHev(ii)
          if (rang==0) WRITE(6,'(3a)') 'NEB: stress tensor in Voigt notation (units: ', cunitP,' ):'
          if (rang==0) WRITE(6,'(a,6g20.8)') unitP*sigPATH(1,1,ii), unitP*sigPATH(2,2,ii), unitP*sigPATH(3,3,ii), &
               0.5*unitP*(sigPath(2,3,ii)+sigPath(3,2,ii)), &
               0.5*unitP*(sigPath(1,3,ii)+sigPath(3,1,ii)), &
               0.5*unitP*(sigPath(1,2,ii)+sigPath(2,1,ii))
          call into_path(ii,1,xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
       end do      !end ii,npath


    case(2)
       if (rang==0) write(6,*)'NEB: -------this is NEB--V2-------'
       if (rang==0) write(6,*)'NEB: The MAX steps in NEB        :',maxneb
       nebtest(:)=0
       do ineb=1,maxneb     ! Main loop for NEB
          !    
          it=0
          !   
          call build_s_path_neb(ityp)       
          ! 
          if (rang==0) print'("NEB:===============================================")'
          if (rang==0) print'("NEB:pas-neb image    force      force_NEB          energie      statut    energie/stable")'
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
                IF (lFire) THEN
                   call trempe_fire(xp, xpp, vp, ax, fp, ielat, iwmax, ityp, &
                        fire_dt(ii), fire_nstep(ii), fire_alph(ii))
                ELSE
                   call trempe(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
                ENDIF

                call analyse 
                call neb_controle(ii) !   (ii,xp, xpp, vp, ax, fp, ielat, iwmax, ityp) 
                !
             end do

             ! Backup image ii
             if (itesauv.GT.0) then
                if (mod(ineb,itesauv)==0) call sauveposition(ii)
             else if (itesauvposition.GT.0) then
                if (mod(ineb,itesauvposition)==0) call sauveposition ( ii)
             endif
             if (itesauvforce.GT.0) then
                if (mod(ineb,itesauvforce)==0) call sauveforce ( ii)
             endif

             !debug          print'("NEB: ",2i5,3E14.5,E20.10,i3)', ineb, ii,  formax,       &
             !debug	            formaxperp, formaxparl, potist*erg2eV,nebtest(ii)
             forneb = SQRT(MAXVAL(force_neb(1,:,ii)**2 + force_neb(2,:,ii)**2         &
                  + force_neb(3,:,ii)**2))*erg2eV/angst 
             if (rang==0) print'("NEB: ",2i5, 2g14.5,g20.10,i3,g15.8)', ineb, ii,  formax, forneb,       &
                  potist*erg2eV,nebtest(ii),potist*erg2eV-enepathev(1)
             ! 
             enePATH(ii)=potist
             enePATHev(ii)=potist*erg2eV       
             sigPATH(:,:,ii) = sigtot(:,:)      ! Contrainte
             call into_path       (ii,1, xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
             !
          end do      ! end ii,path

          if (SUM(nebtest(:))==(npath-2)) then
             if (rang==0) print'("NEB:=== bye,bye sweety is finished=================")'
             if (rang==0) print'("NEB:===============================================")'
             exit
          end if

       end do        !maxneb 

    case default

       if (rang==0) write(6,*) 'NO NEB-DYNAMICS FOR THIS  nebtype = ', nebtype
       if (rang==0) write(6,*) 'CHANGE nebtype AND TRY AGAIN.'
       if (rang==0) write(6,*) 'nebtype=1 for DRAG'
       if (rang==0) write(6,*) 'nebtype=2 for  NEB'
       if (rang==0) write(6,*) 'STOP in the neb.f90'
       stop

    end select

    ! Initialization of fire quench algorithm
    IF (lFire) THEN
       DEALLOCATE(fire_dt)
       DEALLOCATE(fire_nstep)
       DEALLOCATE(fire_alph)
    END IF

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

    if (rang==0) WRITE(*,'(3a)') "NEB:--IMAGE-----REACT-COORD------ENERGY------ENERGY-ENERGY(1)&
         &-------------STRESS-sVoigt(1:6)-(units:-", cunitP, ")"
    do ii=1,npath
       !
       if (rang==0) write(*,'(i5,9(g20.8,1x))')  ii, reaction_coord(ii), enePATHev(ii),enePATHev(ii)-enePATHev(1), &
            unitP*sigPATH(1,1,ii), unitP*sigPATH(2,2,ii), unitP*sigPATH(3,3,ii), &
            0.5*unitP*(sigPath(2,3,ii)+sigPath(3,2,ii)), &
            0.5*unitP*(sigPath(1,3,ii)+sigPath(3,1,ii)), &
            0.5*unitP*(sigPath(1,2,ii)+sigPath(2,1,ii))
       write(55,'(9(g20.8,1x))') dble(ii-1)/dble(npath-1), enePATHev(ii),enePATHev(ii)-enePATHev(1), &
            unitP*sigPATH(1,1,ii), unitP*sigPATH(2,2,ii), unitP*sigPATH(3,3,ii), &
            0.5*unitP*(sigPath(2,3,ii)+sigPath(3,2,ii)), &
            0.5*unitP*(sigPath(1,3,ii)+sigPath(3,1,ii)), &
            0.5*unitP*(sigPath(1,2,ii)+sigPath(2,1,ii))
       write(56,'(9(g20.8,1x))') reaction_coord(ii), enePATHev(ii),enePATHev(ii)-enePATHev(1), &
            unitP*sigPATH(1,1,ii), unitP*sigPATH(2,2,ii), unitP*sigPATH(3,3,ii), &
            0.5*unitP*(sigPath(2,3,ii)+sigPath(3,2,ii)), &
            0.5*unitP*(sigPath(1,3,ii)+sigPath(3,1,ii)), &
            0.5*unitP*(sigPath(1,2,ii)+sigPath(2,1,ii))
       !       
    end do


    if (rang==0) print*,'-------------------------------------------'
    if (rang==0) write(*,'("          Min .. et .. Max")')
    if (rang==0) print '("E:",2g14.7)',minval(enePATHev), maxval(enePATHev)
    if (rang==0) write(*,'("          1 .. et .. NPATH")')
    if (rang==0) print '("B:",2g14.7)',enePATHev(1), enePATHev(npath)
    if (rang==0) write(*,'("les diffs")')
    if (rang==0) print*,'1-NPATH    :',enePATHev(1)-enePATHev(npath)
    if (rang==0) print*,'-------------------------------------------'
    if (rang==0) print*,'MAX-1      :',maxval(enePATHev)-enePATHev(1)
    if (rang==0) print*,'MAX-NPATH  :',maxval(enePATHev)-enePATHev(npath)



  end subroutine neb
end module neb_mod
