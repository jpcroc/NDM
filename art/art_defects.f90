!****************************************************************
      subroutine first_neighbours_distance()
!     for the moments bulk_struncture can be only bcc or fcc 
      use defs
      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use tab_imm_m
      use random_art
      implicit none
       integer,parameter         :: ndr = 200
!       integer, parameter        :: art_ninst=1
           
      real(double)    ::lRmin,lRmax, volume
      integer :: n,nmin,i
      real(double),dimension(ndr)   :: r_distrib
      real(double),dimension(ndr+1) :: art_gdr
      integer, dimension(ndr+1)     :: norm_art_gdr
      integer, dimension(imm)       :: art_defect
      integer                       :: i_defect,id_defect
      real(8)                       :: ran3
      
    volume = volu*angst**3 / dble(im)

    ! Corresponding first nearest neighbour distance in a perfect fcc structure
    ! Radial distribution function will be calculated between lRmin and lRmax

      select case  (bulk_structure)
      !
        case('fcc')
        !
         art_ann = (4.d0*volume)**(1.d0/3.d0)/sqrt(2.d0)
         lRmin = 0.5d0*art_ann ; lRmax = art_ann*sqrt(2.d0)
        !
        case ('bcc')
        !
         art_ann = (2.d0*volume)**(1.d0/3.d0)/sqrt(2.d0)
         lRmin = 0.8d0*art_ann ; lRmax = art_ann*sqrt(2.d0)
	 write(*,*) 'ART:  in the first_neighbours_distance() we have'
	 write(*,*) 'ART:  Rmin',lRmin
	 write(*,*) 'ART:  Rmax',lRmax
	 write(*,*) 'ART:  The fist guess of art_ann:', art_ann 
        !
       case default
         write(*,*) 'ART: WARNING: the fisrt_neighbour_distance work ONLY for the:'
         write(*,*) 'ART: bcc AND fcc case - PLEASE RECPECT THE CASE'
         stop 
       !
       end select 	 
      
    ! Radial interval to calculate radial distribution function
    do n=1, ndr
       r_distrib(n) = lRmin + dble(n-1)/dble(ndr-1)*(lRmax-lRmin)
    end do
      
!      call art_radial_ditribution(ndr,r_distrib,art_gdr)
    call radial_distribuition_function(natoms,type,pos, &
                                      ndr,r_distrib,art_gdr)
      
      norm_art_gdr(1:ndr+1) = INT( NINT (dble(im) * art_gdr(1:ndr+1)  / 2.d0 ) ) 
 




      do n=1,ndr+1
         if (art_ninst.le.SUM(norm_art_gdr(1:n))) then
	 exit
	end if
      end do
      if (n.ge.ndr) then
        write(*,*) 'ART:  ERROR(1001) please increase lRmax in the first_neighbours_distance()'
        write(*,*) 'ART:  ERROR(1001) or icrease the grid in the first_neighbours_distance()'
        stop
      end if
      nmin=n

      art_Rinst=r_distrib(n) + 0.001d0
      


    ! Location of the maximum which should correspond to 1st nearest neighbour
    n = MaxLoc(art_gdr(1:ndr),1)
    if (n.eq.ndr) then
      write(*,*) 'ART: ERROR(1002) please increase lRmax in the first_neighbours_distance()'
      write(*,*) 'ART: ERROR(1002) or icrease the grid in the first_neighbours_distance()'
      stop
    end if
    ! Corresponding radius
    art_ann = 0.5d0*( (r_distrib(n)+r_distrib(n-1))*art_gdr(n) +        &
         (r_distrib(n-1)+r_distrib(n-2))*art_gdr(n-1) +             &
    (r_distrib(n+1)+r_distrib(n))*art_gdr(n+1) )/(art_gdr(n-1) +    &
       art_gdr(n) + art_gdr(n+1))

    WRITE(*,*)
    WRITE(*,'(a,g14.6,a,f0.5,a,i0,a)') 'ART:  Radial distribution maximum ', art_gdr(n), &
        'in R = ', r_distrib(n), ' (n=',n,')'
    WRITE(*,'(a,f0.5,a)') 'ART:  Found 1st nearest neighbour distance :    art_ann   = ', art_ann, ' A'
    WRITE(*,'(a,f0.5,a)') 'ART:  Found inst. dumbells dist. < than    :    Rinst = ', art_Rinst, ' A'

     call radial_distribuition_atoms(natoms,type,pos, &
                  nmin, ndr,r_distrib,art_defect,i_defect)
 

   write(*,*) ' I preferred_atom, preferred_atom', preferred_atom
   write(*,*) ' Number of defined        INST ', art_ninst
   write(*,*) ' Number of found atoms in INST ', i_defect
   
   id_defect = INT ( i_defect * ran3() + 1)    ! Between 1 and i_defect

   preferred_atom=art_defect(id_defect)
   write(*,*) ' II preferred_atom, preferred_atom', preferred_atom
      
   end subroutine first_neighbours_distance
      


!****************************************************************

      subroutine radial_distribuition_function(natoms,types,pos, &
                                      ndr,r_distrib,art_gdr)
 !****************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
   USE T_kind_param_m, ONLY:  double
   use tab_imm_m
   use gen_com_m
   use art_in_ndm_module


   real(double),dimension(ndr+1) :: art_gdr
   integer,dimension(ndr+1)      :: ngdr
   real(double),dimension(ndr)   :: r_distrib,r_distrib2
   real(double), dimension(3)    :: dxp(3)
   integer  :: i,j,iti,iw1,iw2

!input output  decalarations for the ART parameters
   integer, intent(in) :: natoms
   integer, intent(in), dimension(natoms):: types
   real(8), intent(in), dimension(3*natoms), target:: pos
   real(8), dimension(:), pointer :: x, y, z
  
    x => pos(1:NATOMS)
    y => pos(NATOMS+1:2*NATOMS)
    z => pos(2*NATOMS+1:3*NATOMS)
    do i=1,NATOMS 
      xp(1,i)=x(i)/angst
      xp(2,i)=y(i)/angst
      xp(3,i)=z(i)/angst
    enddo

!    write(*,*) 'at', at

!    write(*,*) 'bg', bg
    
    
    
    call caltabt 
    call caltabi 
    
    r_distrib2(:) = r_distrib(:)**2
    art_gdr(:)=0
    ngdr(:)=0 
    

      call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst

     
      
      iw2=0
      
      do i=1,im
       iti = ityp(i)
       iw1 = iw2+1
       iw2 = iwmax(i)
        do iw=iw1,iw2
          j = indi(iw)
!	  if (j.gt.i) then
          dxp(1:3) = xp(1:3,i) - xp(1:3,j)
             WHERE ( (dxp.GT.0.5d0).OR.(dxp.LT.-0.5d0) )
                      dxp(1:3) = dxp(1:3) - Dble(Nint(dxp(1:3)))
             END WHERE
!          if (ltriclin) then
                  dxp = MatMul(at,dxp)*angst
!              else
!	          dxp(:) = dxp(:)*dabs(zl(:))*angst
!	  end if  	  		  
	  ! Calcul du carré de la distance
          do izero=1,3
	   if (dabs(dxp(izero)).lt.low_limit) then
	    dxp(izero) = zero
	   end if
 	  end do 
          dxp2 = Sum( dxp(1:3)**2 )
            DO n=ndr, 0, -1
               IF (n.EQ.0) THEN
                    ngdr(1) = ngdr(1)+1
                 ELSE IF (dxp2.GE.r_distrib2(n)) THEN
                    ngdr(n+1) = ngdr(n+1) + 1
                    exit
                END IF
             END DO

!        end if
        end do
      end do 
      
    art_gdr(1:ndr+1) = 2.d0*dble(ngdr(1:ndr+1))/dble(im)
      
      end    subroutine radial_distribuition_function

!****************************************************************









!****************************************************************

      subroutine radial_distribuition_atoms(natoms,types,pos, &
                  nmin,ndr,r_distrib,art_defect,i_defect_max)
 !****************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
   USE T_kind_param_m, ONLY:  double
   use tab_imm_m
   use gen_com_m
   use art_in_ndm_module


   real(double),dimension(ndr)   :: r_distrib,r_distrib2
   real(double), dimension(3)    :: dxp(3)
   integer  :: i,j,iti,iw1,iw2,i_defect,i_defect_max

!input output  decalarations for the ART parameters
   integer, intent(in) :: natoms
   integer, intent(in), dimension(natoms):: types
   real(8), intent(in), dimension(3*natoms), target:: pos
   real(8), dimension(:), pointer :: x, y, z
   integer,dimension(imm)    :: art_defect 
  
    x => pos(1:NATOMS)
    y => pos(NATOMS+1:2*NATOMS)
    z => pos(2*NATOMS+1:3*NATOMS)
    do i=1,NATOMS 
      xp(1,i)=x(i)/angst
      xp(2,i)=y(i)/angst
      xp(3,i)=z(i)/angst
    enddo

!    write(*,*) 'at', at

!    write(*,*) 'bg', bg
    
    
    
    call caltabt 
    call caltabi 
    
    r_distrib2(:) = r_distrib(:)**2
    art_defect(:)=0
    
      call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst

     
      
      iw2=0
      i_defect=0
      
      do i=1,im
       iti = ityp(i)
       iw1 = iw2+1
       iw2 = iwmax(i)
        do iw=iw1,iw2
          j = indi(iw)
!	  if (j.gt.i) then
          dxp(1:3) = xp(1:3,i) - xp(1:3,j)
             WHERE ( (dxp.GT.0.5d0).OR.(dxp.LT.-0.5d0) )
                      dxp(1:3) = dxp(1:3) - Dble(Nint(dxp(1:3)))
             END WHERE
!          if (ltriclin) then
                  dxp = MatMul(at,dxp)*angst
!              else
!	          dxp(:) = dxp(:)*dabs(zl(:))*angst
!	  end if  	  		  
	  ! Calcul du carré de la distance
          do izero=1,3
	   if (dabs(dxp(izero)).lt.low_limit) then
	    dxp(izero) = zero
	   end if
 	  end do 
          dxp2 = Sum( dxp(1:3)**2 )
              IF (dxp2.LE.r_distrib2(nmin)) THEN
	       i_defect=i_defect+1
               art_defect(2*i_defect - 1)=i
               art_defect(2*i_defect)=j
	       write(*,*) i,j
               END IF

!        end if
        end do
      end do 

        i_defect_max=2*i_defect
      

      
      
      
      
      end    subroutine radial_distribuition_atoms

!****************************************************************
