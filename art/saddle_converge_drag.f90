! Subroutine saddle_convert
!
! This subroutine bring the configuration to a saddle point. It does that
! by first pushing the configuration outside of the harmonic well, using
! the initial direction selected in fi\left[ \left\lbrace n\left\langle d\left) \left] \left\rbrace _\left\rangle s\left. addle. Once outside the harmonic
! well, as defined by the appearance of a negative eigenvalue (or
! reasonnable size) the configuration follows the direction corresponding
! to this eigenvalue until the force components parallel and perpdendicular
! to the eigendirection become close to zero.
!
!  Normand Mousseau, June 2001

!old subroutine saddle_converge(ret, saddle_energy, fparl_max, fperp_max)
 subroutine saddle_converge(ret, saddle_energy, fparl_max, fperp_max)
  use random_art
  use defs
  use saddles
  use lanczos_defs
  use art_in_ndm_module
  implicit none

  integer, intent(out) :: ret
  real(8), intent(out) :: saddle_energy  !old , fpar, fperp
  real(8)               :: fpar,fperp
  logical :: new_projection

  integer :: i, j, k, kter, iter, iperp, k_rejected, eigen_rejected,ij
  integer :: preafor, maxvec,npart,itry, iter_init, kter_init, ierror
  real(8) :: fdotinit, fperp2, current_fperp, fnorm
  real(8) :: step, delr, diff, mu, nu, diff_eignorm
  real(8) :: one_art = 1.0d0
  real(8) :: boxl(3), current_energy,energy_before_md
  real(8), dimension(VECSIZE) :: posb, perp_force, forceb, perp_forceb, diff_eig
  real(8), dimension(VECSIZE) :: H_proj, test_pos_1, test_pos1, test_force_1, test_force1
  real(8), dimension(VECSIZE) :: forcei, force_i
!  real(8), dimension(VECSIZE) :: count_saddle   ! count_saddle(0) = nombre de points cols trouves
! local variables .... by me
  real(8)                              ::   posprov
  real(8), dimension(VECSIZE)          ::   vit_art,parl_force, posi, pos_i
  real(8), ALLOCATABLE, dimension(:,:) ::   hess
  real(8), dimension(VECSIZE/3)        ::   fparl_ia, fperp_ia
  real(8), intent(out)                 ::   fparl_max, fperp_max
  integer                              ::   ic_eigenvalue,i_before_inflex,i_after_inflex   
  character*1                          ::   coucou,moucou
  real(8)                              ::   deplasare
  real(8)                              ::   rcm_loc(3),masstot_contr,norm_loc,projection_md(VECSIZE),r_md(3)
  ! We compute at constant volume
  boxl(:) = box(:) * scala
  kter_init = 0

  if ( NEW_EVENT ) then 
write(*,*) 'New event..........'
     eigenvalue = 0.0d0
     call calcforce(NATOMS,type,pos,boxl,force,current_energy)
     
     ! We now project out the direction of the initial displacement from the
     ! minimum from the force vector so that we can minimize the energy in
       
      fdotinit=DOT_PRODUCT(force(:),initial_direction(:))
      perp_force(:)  = force(:) - fdotinit * initial_direction(:)  ! Vectorial force
     
!debug_new      step = INCREMENT

     ! On cherche a sortir du puits
     i_before_inflex=0
     do kter = kter_init, MAXKTER
        k = 0
        k_rejected = 0
        eigen_rejected = 0

        ! We relax perpendicularly using a simple variable-step steepest descent 
        vit_art(:) = 0.d0
        posb(:) = pos(:)
        call calcforce(NATOMS,type,pos,boxl,force,current_energy)
        
	fdotinit = DOT_PRODUCT(force,initial_direction)
        perp_force(:)  = force(:) - fdotinit * initial_direction(:)  ! Vectorial force

        ! on va minimiser dans l'hyperplan perpendiculaire

     
        energy_before_md=current_energy 
        do 
           do i=1,VECSIZE
              if (vit_art(i)*perp_force(i)>0) then
                 posprov = 2*posb(i) - pos(i) + tmass_art(i)*perp_force(i)
              else 
                 posprov = posb(i) + tmass_art(i) * perp_force(i)
              end if
              vit_art(i) = (posprov-pos(i))*usdh_art
              pos(i) = posb(i)
              posb(i) = posprov
           end do
           
           call calcforce(NATOMS,type,posb,boxl,forceb,total_energy)
           
	   evalf_number = evalf_number + 1
           i_before_inflex=i_before_inflex+1
	            
           fdotinit = DOT_PRODUCT(forceb, initial_direction)
           perp_forceb  = forceb - fdotinit * initial_direction  ! Vectorial force
           fperp2 = DOT_PRODUCT( perp_forceb(:),perp_forceb(:))
           fperp = sqrt(fperp2)
           force = forceb
           perp_force = perp_forceb
           current_energy = total_energy
           k = k + 1
           !old_art           if(fperp2 < FTHRESH2 .or. k > MAXKPERP .or. k_rejected > 5) then
           if(fperp2 < FTHRESH2 .or. k > MAXKPERP ) then
              pos(:)=posb(:)
              exit
           end if
        end do
        if ((total_energy-energy_before_md).gt.1.d-2) then
	  write(*, '(" WARNING FINDING INFLEXION POINT: After ",i9," MD steps the energy increased with: ",f18.2)')k, total_energy-energy_before_md
	end if 
        ! We now move the configuration along the initial direction and check the lowest
        ! eigendirection 

        fpar = DOT_PRODUCT(forceb, initial_direction)
        pos = pos - INCREMENT * sign(one_art,fpar) * initial_direction   !Vectorial operation        
        it_art=0

	! We start checking of negative eigenvalues only after a few steps
        if( kter>= KTER_MIN ) then
           if (kter==KTER_MIN) then
              new_projection = .true.    ! We do not  use previously computed lowest direction as seed
              maxvec = NVECTOR_LANCZOS
              call lanczos(maxvec,new_projection)
              i=0
              !	   write(*, '(i5,f18.8)') i, eigenvalue
              !	   write(26,'(i5,f18.8)') i, eigenvalue
              !crc         do i=1,300
              !         do i=1,5
              !           new_projection = .false.    
              !           maxvec = NVECTOR_LANCZOS
              !           call lanczos(maxvec,new_projection)
              !	   write(*, '(i5,f18.8)') i, eigenvalue
              !	   write(26,'(i5,f18.8)') i, eigenvalue
              !	 end do  
              !	   write(*,*) eigenvalue
              !	   stop
           else
              new_projection = .false.    
              maxvec = NVECTOR_LANCZOS
              call lanczos(maxvec,new_projection)
           end if
        endif
        current_energy = total_energy ! As computed in lanczos routine
        call displacement(posref, pos, delr,npart)
         
        if (print_details .and. mod(kter,kprint).eq.0) then 
           write(*,"(' ','kter: ',i4,' k min ', i4,' Energy: ',f16.6,'  e-val: ',d12.6,'  delr: ',f12.6,' evalf: ',i6)") &
                &  kter,k,current_energy, eigenvalue, delr, evalf_number
           write(FLOG,"(' ','kter: ',i4,'  Energy: ',f16.6,'  e-val: ',d12.6,'  delr: ',f12.6,' evalf: ',i6)") &
                &  kter,current_energy, eigenvalue, delr, evalf_number
        endif
        
        ! because we have updated the positions ....we compute the forces
        call calcforce(NATOMS,type,pos,boxl,force,total_energy)

        if (eigenvalue <  EIGEN_THRESH) exit
     end do
     

      ! The configuration is now out of the harmonic well, we can now bring
      ! it to the saddle point. Again, we split the move into a parallel and 
      ! perpendicular contribution.
     
      ! First, we must now orient the direction of the eigenvector corresponding to the
      ! negative eigendirection (called projection) such that it points away from minimum. 
     
      fpar =DOT_PRODUCT(force(:),projection(:))
!ooo      if(fpar > 0.0d0 ) projection = -1.0d0 * projection
      !WARNING - DANGEROUS FOR INCREMENT > 0.8 A - te joci cu focu' nepoate!      
      !      pos = pos +  INCREMENT * projection 
     
      fdotinit     =  DOT_PRODUCT(force(:), projection(:))
      perp_force(:)= force(:) - fdotinit * projection(:)  ! Vectorial force
     
      fperp2 =  DOT_PRODUCT(perp_force(:), perp_force(:))
      current_fperp = sqrt(fperp2)
      iter_init = 1
     
  else if (.not.NEW_EVENT) then       ! This is a convergence event
                   
     write(*,*) 'In NOT NEW EVENT.............................'                   
     ! We must compute the eigenvector with precision
     maxvec = NVECTOR_LANCZOS
     new_projection = .true. 
     call lanczos(maxvec,new_projection)

     new_projection = .false.
     call lanczos(maxvec,new_projection)

     new_projection = .false.
     call lanczos(maxvec,new_projection)

     fpar = DOT_PRODUCT(force(:), projection(:))
     if(fpar > 0.0d0 ) projection = -1.0d0 * projection

     fdotinit = DOT_PRODUCT(force(:),projection(:))
     perp_force  = force - fdotinit * projection  ! Vectorial force
    
     fperp2 = DOT_PRODUCT(perp_force(:),perp_force(:))
     fperp = sqrt(fperp2)
     current_fperp = fperp
     
  endif

! on va chercher a remonter vers le saddle point
  ic_eigenvalue=0
  i_after_inflex=0 
  do iter = iter_init, MAXITER 
    !    tstep_art=1.d0

! relaxation in the hyperplane        

       coucou=' '
 !      write(*,"('B','iter: ',i4,'  iper: ',i4,'  ener: ',f10.4,'  fpar: ',f10.4,'  fperp: ',f10.4,&
 !           & '  e-val: ', f10.4,'  delr: ',a1,f10.4,'  incr: ',a1,f10.4,'  evalf: ',i6)")  &
 !           & iter, iperp, total_energy, fpar, fperp_max, eigenvalue, coucou,(fpar / eigenvalue),coucou, INCREMENT, evalf_number
    
    posb(:)=pos(:)
    it_art=0
    rcm_loc(1-3)=0.d0
    do i=1,NATOMS
             rcm_loc(1) = pos(i)*tmass_art(i) + rcm_loc(1)
             rcm_loc(2) = pos(i+NATOMS)*tmass_art(i+NATOMS) + rcm_loc(2)
             rcm_loc(3) = pos(i+2*NATOMS)*tmass_art(i+2*NATOMS) + rcm_loc(3)
             masstot_contr=masstot_contr+ tmass_art(i)+tmass_art(i+NATOMS)+tmass_art(i+2*NATOMS)  
    end do
    rcm_loc(1)=rcm_loc(1)/masstot_contr
    rcm_loc(2)=rcm_loc(2)/masstot_contr
    rcm_loc(3)=rcm_loc(3)/masstot_contr
    norm_loc=SUM(rcm_loc(:)**2)
    rcm_loc(:)=rcm_loc(:)/dsqrt(norm_loc)
    vit_art(:)=0
    

      r_md(1)=SUM(projection(1:NATOMS))
      r_md(2)=SUM(projection(1+NATOMS:2*NATOMS))
      r_md(3)=SUM(projection(1+2*NATOMS:3*NATOMS))

    do i=1,NATOMS
!ooo      projection_md(i)=projection(i)-rcm_loc(1)
!ooo      projection_md(i+NATOMS)=projection(i+NATOMS)-rcm_loc(2)
!ooo      projection_md(i+2*NATOMS)=projection(i+2*NATOMS)-rcm_loc(3)
    projection_md(i)=r_md(1)
    projection_md(i+NATOMS)=r_md(2)
    projection_md(i+2*NATOMS)=r_md(3)
    end do 
    norm_loc=SUM(projection_md(:)**2) 
    fpar=DOT_PRODUCT(force(:),projection(:))
    parl_force(:)  = fpar * projection(:)
    
    perp_force(:)  = force(:) - parl_force(:)+tmass_art(i)*projection_md(:)/masstot_contr
     
    reject = .false.
    itry = 0
    iperp = 0
    
    energy_before_md=total_energy
    do  
        do i=1,VECSIZE
         
         
         if (vit_art(i)*perp_force(i)>0) then
             posprov = 2*posb(i) - pos(i) + tmass_art(i)*perp_force(i)!-0.0001*vit_art(i)
          else 
             posprov = posb(i) + tmass_art(i) * perp_force(i)!-0.0001*vit_art(i)
          end if
           vit_art(i) = (posprov-pos(i))*usdh_art
           pos(i) = posb(i)
          posb(i) = posprov
        end do

       call displacement(posb, pos, delr,npart)
       
       call calcforce(NATOMS,type,posb,boxl,forceb,total_energy)
       evalf_number = evalf_number + 1
       i_after_inflex=i_after_inflex+1
       
       fpar=DOT_PRODUCT(forceb(:),projection(:))
       parl_force(:)  = fpar * projection(:)
       perp_force(:)  = forceb(:) - parl_force(:)+tmass_art(i)*projection_md(:)/masstot_contr
       
       fperp2=DOT_PRODUCT(perp_force(:),perp_force(:))
       fperp = sqrt(fperp2)
       
       force = forceb
       current_energy = total_energy
       current_fperp = fperp
       iperp = iperp + 1
       fnorm = DOT_PRODUCT(force(:),force(:))
       
       if(fperp2 < FTHRESH2   .or.                &
!            iperp > (iter + 1)  .or.                &
            iperp > MAXIPERP    &  
	  !.or.   itry > 50             &
            ) then
          pos(:) = posb(:)
          ret=-30000
          exit
       end if
    end do      ! relaxation in the hyperplane   
    
        if ((total_energy-energy_before_md).gt.1.d-2) then
	  write(*, '(" WARNING CLIMBING: After ",i9," MD steps the energy increased: ",f18.2)') iperp,total_energy-energy_before_md
	end if 

    ! we check if the saddle point is reached
    do i=1,NATOMS
       fparl_ia(i)=dsqrt(parl_force(i)**2+        &
                        parl_force(i+NATOMS)**2+ &
			parl_force(i+2*NATOMS)**2)
      
      fperp_ia(i)=dsqrt(perp_force(i)**2+        &
                        perp_force(i+NATOMS)**2+ &
			perp_force(i+2*NATOMS)**2)
   end do
   fparl_max=MAXVAL(fparl_ia(:)) 
   fperp_max=MAXVAL(fperp_ia(:))
   
   call displacement(posref, pos, delr,npart)
   !write(*,*) 'new_projection = ', new_projection

   if (print_details .and. mod(iter,iprint) .eq. 0) then 
       write(FLOG,"(' ','iter: ',i4,'  iper: ',i4,'  ener: ',f10.4,'  fpar: ',f10.4,'  fperp: ',f10.4,&
            & '  e-val: ', d10.4,'  delr: ',f10.4,'  npart: ',i4,'  evalf: ',i6)")  &
            & iter,  iperp, total_energy, fparl_max, fperp_max, eigenvalue, delr, npart, evalf_number
   endif

!ooo       write(*,"(' ','iter: ',i4,'  iper: ',i4,'  ener: ',f10.4,'  fpar: ',f10.4,'  fperp: ',f10.4,&
!ooo            & '  e-val: ', f10.4,'  delr: ',f10.4,'  npart: ',i4,'  evalf: ',i6)")  &
!ooo            & iter, iperp, total_energy, DOT_PRODUCT(force(:), projection(:)), fperp_max, eigenvalue, delr, npart, evalf_number
!ooo    endif

    saddle_energy = current_energy
    
       if ( (fparl_max + fperp_max)< EXITTHRESH)  then
          ret = 20000 + iter
          write(*,'("debug purpose ..........:",i7,4f15.4)'), ret,fparl_max,fperp_max,fparl_max+fperp_max, EXITTHRESH
       else if ( (abs(fparl_max) < 0.1*EXITTHRESH) .and. (fperp_max<EXITTHRESH) ) then
          ret = 10000 + iter
          write(*,'("debug purpose ..........:",i7,4f15.4)'), ret,fparl_max,fperp_max,fparl_max+fperp_max, EXITTHRESH
       endif
!    end if
    
      nu = 1d-3
      if (((ret == 10000+iter) .or. (ret == 20000+iter)).and.(ic_eigenvalue == 0)) then
	count_saddle(1) = count_saddle(1) + 1
	
	write(*,'("kter  iter  hyper lanczos  avant  apres  tot ", 2i4, i5, i5, i5, 2i5, i6)') &
	 kter,iter, i_after_inflex + i_before_inflex,                       &
	 evalf_number-i_after_inflex -i_before_inflex,                      &
	 (kter-KTER_MIN)*NVECTOR_LANCZOS+1+i_before_inflex,                 &
	 evalf_number-(kter-KTER_MIN)*NVECTOR_LANCZOS-1-i_before_inflex,    &
	 evalf_number
       exit
      end if
     
      if (((ret == 10000+iter) .or. (ret == 20000+iter)).and.ic_eigenvalue >= 1) then  
       it_art=0
       new_projection = .false.    ! We start from the previuos direction each time
       maxvec = NVECTOR_LANCZOS
       call lanczos(maxvec,new_projection)
       if (eigenvalue < 0) then
	count_saddle(1) = count_saddle(1) + 1
	write(*,'("kter  iter  hyper lanczos  avant  apres  tot ", 2i4, i5, i5, i5, 2i5, i6)') &
	 kter,iter, i_after_inflex + i_before_inflex,                       &
	 evalf_number-i_after_inflex -i_before_inflex,                      &
	 (kter-KTER_MIN)*NVECTOR_LANCZOS+1+i_before_inflex,                 &
	 evalf_number-(kter-KTER_MIN)*NVECTOR_LANCZOS-1-i_before_inflex,    &
	 evalf_number
       exit
      end if
      end if


!? what is that ... releted to the old fashion to put the convergency.
    fpar = DOT_PRODUCT(force(:), projection(:))
!    if(fpar > 0.0d0 ) projection = -1.0d0 * projection
    fdotinit = fpar
    !write(*,*) 'fpar = ', fpar

     fnorm = DOT_PRODUCT(force(:),force(:))
     
    ! We now move the configuration along the eigendirection corresponding
    ! to the lowest eigenvalue
    
        
        deplasare=INCREMENT
	moucou=' '
	coucou=' '

	!if ((abs(fdotinit/eigenvalue) < 3 )) then
!ooo        if ((eigenvalue < -0.5).or.eigenvalue > 0  ) then ! .and. (fnorm < 20.0))  then
        if (dabs(eigenvalue) > 0.0 ) then ! .and. (fnorm < 20.0))  then
        !if (fnorm < 0.5d0) then 
           !if (fpar > 0) then
                !pos = pos - sign(one_art,fpar) * (fpar / eigenvalue) * projection       !KIMIYA!
           !write(*,*) 'We''re in my section now.....'
!ooo           if(fpar < 0.0d0 ) then
!ooo              projection = - 1.0d0 * projection
!ooo           end if

           ! KM6
	    deplasare =   sign(one_art,fpar) * dabs(fpar / eigenvalue)
	   ! deplasare =   fpar / eigenvalue
!Wales            deplasare=   sign(one_art,fpar) * 2.d0*fpar / eigenvalue**2*(1+dsqrt(1+4.D0*fpar**2/eigenvalue**4))
!	   deplasare =  (fpar / eigenvalue)
!debug	   write(*,*) 'deplasare: ', deplasare,  'normal: ', sign(one_art,fpar) * dabs(fpar / eigenvalue) , 'forced: ',INCREMENT * sign(one_art,fpar)
    	    moucou=' '
	    coucou='*'
!	   if ( dabs(deplasare) > INCREMENT ) deplasare =  INCREMENT * sign(one_art,fpar) !/ sqrt(1.0d0*iter)
	   if ( dabs(deplasare) > 0.001 )  then
	   
	    ! deplasare =  0.001 * sign(one_art,fpar) !/ sqrt(1.0d0*iter)
	    deplasare =  0.001*sign(one_art,fpar) 
	     coucou=' '
    	     moucou='*'
          end if
!debug	   write(*,*) '*deplasare: ', deplasare,  'normal: ', sign(one_art,fpar) * dabs(fpar / eigenvalue) , 'forced: ',INCREMENT * sign(one_art,fpar)
	   
           pos(:) = pos(:) - deplasare * projection(:)   !/ (1.0d0*iter**0.5)

        else

            if( fparl_max > EXITTHRESH*6.0 ) then 
               deplasare=  sign(one_art,fpar) * 0.8d0 * INCREMENT / sqrt(1.0d0*iter) 
            else
               deplasare=  sign(one_art,fpar) * 0.3d0 * INCREMENT / sqrt(1.0d0*iter) 
            endif
	    
            pos = pos - deplasare * projection / sqrt(1.0d0*iter) 
    	    moucou=' '
            coucou=' '

        endif

       write(*,"('D','iter: ',i4,'  iper: ',i4,'  ener: ',f10.4,'  fpar: ',f10.4,'  fperp: ',f10.4,&
            & '  e-val: ', d10.4,'  inc1: ',a1,f10.4,'  inc2: ',a1,f10.4,' overlap: ',d8.2,'  evalf: ',i6,'  delr: ',f10.4)")  &
            & iter, iperp, total_energy, fpar, fperp_max, eigenvalue, coucou,sign(one_art,fpar) * dabs(fpar / eigenvalue), &
	    moucou, deplasare , overlap, evalf_number,delr
 
!--------------LANCZOS
!--------------LANCZOS
    if ((abs(fparl_max)<EXITTHRESH).and.(iter>2)) then
         MAXIPERP=MAXIPERP_ORIG*4
	else 
	MAXIPERP=MAXIPERP_ORIG 
     end if	  
	 
         

!       it_art=0
!       write(*,*) 'NO LANCZOS'
    !
!    else
       !
       it_art=0
       new_projection = .false.    ! We start from the previuos direction each time
       maxvec = NVECTOR_LANCZOS
       call lanczos(maxvec,new_projection)
       if (reject) eigen_rejected = eigen_rejected + 1  
       !
!    end if
    
! if we have more than 2 consecutive positive eigenvalue 
! we leave this region

    if (eigenvalue < 0.0 ) then
       ic_eigenvalue=0 
    else 
       ic_eigenvalue=ic_eigenvalue+1
    end if      
    if (ic_eigenvalue> 2 ) then
      write(*,*)   'THIS IS NOT A SADDLE POINT ic_eigenvalue', ic_eigenvalue
      write(FLOG,*) 'THIS IS NOT A SADDLE POINT ic_eigenvalue', ic_eigenvalue
      count_saddle(2) = count_saddle(2) + 1
      ret = 60000 + iter
      exit
    end if
!--------------LANCZOS
!--------------LANCZOS

   it_art=0
   call calcforce(NATOMS,type,pos,boxl,force,total_energy)
   evalf_number = evalf_number + 1

   fpar = DOT_PRODUCT(force(:),projection(:))
!ooo   if(fpar > 0.0d0 ) projection = -1.0d0 * projection
   fdotinit = DOT_PRODUCT(force(:),projection(:))

   perp_force(:)  = force(:) - fdotinit * projection(:)  ! Vectorial force
   
   fperp2=DOT_PRODUCT(perp_force(:),perp_force(:))
   fperp = sqrt(fperp2)
   
   current_fperp = fperp
   
end do

!open(unit = CSAD,file = CSADDLE) !,status='unknown',action='write',position='append',iostat=ierror)
if (iter == MAXITER) then
   count_saddle(2) = count_saddle(2)+1
end if
write(CSAD,*) count_saddle(1),  count_saddle(2), evalf_number

write(*,*) '*** FOUND SADDLE TO NOT FOUND IS :::::  ', count_saddle(1), ' -- ', count_saddle(2) 
end subroutine saddle_converge
