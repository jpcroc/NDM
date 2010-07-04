! Subroutine saddle_convert
!
! This subroutine bring the configuration to a saddle point. It does that
! by first pushing the configuration outside of the harmonic well, using
! the initial direction selected in find_saddle. Once outside the harmonic
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

  integer :: i, j, k, kter, iter, iperp, k_rejected, eigen_rejected
  integer :: preafor, maxvec,npart,itry, iter_init, kter_init, ierror
  real(8) :: fdotinit, fperp2, current_fperp, Hd2, gHd, forceb2, force2, force2_, forceb2_
  real(8) :: step,delr, diff, tiny
  real(8) :: one_art = 1.0d0
  real(8) :: boxl(3), current_energy
  real(8), dimension(VECSIZE) :: posb, perp_force, forceb, perp_forceb, force_d, Hd, posb_d, d_prev
  real(8), dimension(VECSIZE) :: force_new, posprov_, force_k, perp_force_, force_, force_k_
! local variables .... by me
  real(8)                       ::   posprov
  real(8), dimension(VECSIZE)   ::   vit_art,parl_force
  real(8), dimension(VECSIZE/3) ::   fparl_ia,fperp_ia
  real(8),intent(out)           ::   fparl_max, fperp_max
  integer                       ::   ic_eigenvalue   

  ! We compute at constant volume
  boxl(:) = box(:) * scala
  kter_init = 0
 
  if ( NEW_EVENT ) then 

     eigenvalue = 0.0d0
     call calcforce(NATOMS,type,pos,boxl,force,current_energy)
     
     ! We now project out the direction of the initial displacement from the
     ! minimum from the force vector so that we can minimize the energy in
       
!old_art     fdotinit= 0.0d0
!old_art     do i=1, VECSIZE
!old_art        fdotinit = fdotinit + force(i) * initial_direction(i)
!old_art     end do
     fdotinit=DOT_PRODUCT(force(:),initial_direction(:))
     perp_force(:)  = force(:) - fdotinit * initial_direction(:)  ! Vectorial force
     
!debug_new      step = INCREMENT

     ! On cherche a sortir du puits
     do kter = kter_init, MAXKTER
        k = 0
        k_rejected = 0
        eigen_rejected = 0

        ! We relax perpendicularly using a simple variable-step steepest descent 
        vit_art(:)=0.d0
        posb(:)=pos(:)
        fdotinit=DOT_PRODUCT(force,initial_direction)
        perp_force(:)  = force(:) - fdotinit * initial_direction(:)  ! Vectorial force

        ! on va minimiser dans l'hyperplan perpendiculaire
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
           
           !old_art           fdotinit= 0.0d0  
           !old_art	   do i=1, VECSIZE
           !old_art              fdotinit = fdotinit + forceb(i) * initial_direction(i)
!old_art           end do
           fdotinit=DOT_PRODUCT(forceb,initial_direction)
           perp_forceb  = forceb - fdotinit * initial_direction  ! Vectorial force
           !old_art           fperp2 = 0.0d0
           !old_art           do i=1, VECSIZE
           !old_art              fperp2 = fperp2 + perp_forceb(i) * perp_forceb(i)
!old_art           end do
           fperp2=DOT_PRODUCT( perp_forceb(:),perp_forceb(:))
           fperp = sqrt(fperp2)
           
           !old_art           if(total_energy < current_energy ) then
           !old_art               tstep_art = 1.1 * tstep_art
           !old_art              k_rejected = 0
           !old_art           else
           !old_art              tstep_art = 0.8 * tstep_art
!old_art              k_rejected = k_rejected + 1
           !old_art           endif
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
        
        ! We now move the configuration along the initial direction and check the lowest
        ! eigendirection 
!original        pos = pos + INCREMENT * initial_direction   !Vectorial operation
        pos = pos + INCREMENT * initial_direction   !Vectorial operation
        
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
           write(*,"(' ','kter: ',i4,' k min ', i4,' Energy: ',f16.6,'  e-val: ',f12.6,'  delr: ',f12.6)") &
                &  kter,k,current_energy, eigenvalue, delr
           write(FLOG,"(' ','kter: ',i4,'  Energy: ',f16.6,'  e-val: ',f12.6,'  delr: ',f12.6)") &
                &  kter,current_energy, eigenvalue, delr
        endif
        
	
        ! because we have updated the positions ....we compute the forces
        call calcforce(NATOMS,type,pos,boxl,force,total_energy)

        if(eigenvalue <  EIGEN_THRESH) exit
     end do
     
     
     ! The configuration is now out of the harmonic well, we can now bring
     ! it to the saddle point. Again, we split the move into a parallel and 
     ! perpendicular contribution.
     
     ! First, we must now orient the direction of the eigenvector corresponding to the
     ! negative eigendirection (called projection) such that it points away from minimum. 
     
!     !WARNING - DANGEROUS FOR INCREMENT > 0.8 A - te joci cu focu' nepoate!      
!     !      pos = pos +  INCREMENT * projection 
!

     iter_init = 1
     
  else if (.not.NEW_EVENT) then       ! This is a convergence event
                                      
     ! We must compute the eigenvector with precision
     maxvec = NVECTOR_LANCZOS
     new_projection = .true. 
     call lanczos(maxvec,new_projection)
     
     new_projection = .false.
     call lanczos(maxvec,new_projection)
     
     new_projection = .false.
     call lanczos(maxvec,new_projection)
     
     
  endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! BEGIN KIMIYA SECTION !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       forceb2 = DOT_PRODUCT(forceb(:),forceb(:))
!        force2 = DOT_PRODUCT(force(:),force(:))
!       force_k = force + (forceb2 / force2) * d_prev;
!       fdotinit = DOT_PRODUCT(force_k(:),projection(:))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  END  KIMIYA SECTION !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     fpar= DOT_PRODUCT(force(:),projection(:))
     if(fpar > 0.0d0 ) projection = -1.0d0 * projection

     fdotinit     =  DOT_PRODUCT(force(:), projection(:))
     perp_force  = force - fdotinit * projection  ! Vectorial force
!     d_prev = perp_force
     fperp2 = DOT_PRODUCT(perp_force(:),perp_force(:));

     fperp = sqrt(fperp2)
     current_fperp = fperp


! on va chercher a remonter vers le saddle point
  ic_eigenvalue=0
  do iter = iter_init, MAXITER 
     
    reject = .false.
    itry = 0
    iperp = 0

    vit_art(:)=0.d0
    posb(:)=pos(:)


    force2 = DOT_PRODUCT(force(:),force(:))

    force_ = force
    perp_force_ = perp_force

!    do
       !    tstep_art=1.d0
       
       tiny = 0.001d0
       posb_d = posb + tiny * perp_force_
       call calcforce(NATOMS,type,posb_d,boxl,force_d,total_energy) 
       Hd = (-force_d + force_) / tiny
       Hd2 = DOT_PRODUCT(Hd(:),Hd(:))
       gHd = - DOT_PRODUCT(force_(:),Hd(:))

! relaxation in the hyperplane        
!    do  
!       do i=1,VECSIZE
!          if (vit_art(i)*perp_force(i)>0) then
              posprov_ = posb - (gHd / Hd2) * perp_force_                 !! USED BY KIMIYA  
!               posprov = 2*posb(i) - pos(i) + tmass_art(i)*perp_force(i)    !! USED BY COSMIN
             !	     posprov = 2*posb(i) - pos(i) + tstep_art*tmass_art(i)*perp_force(i)
!          else 
             !             posprov = posb(i) + tstep_art*tmass_art(i) * perp_force(i)
!             posprov = posb(i) + tmass_art(i) * perp_force(i)
!          end if

!          vit_art(i) = (posprov-pos(i))*usdh_art
          pos = posb
          posb = posprov_
!       end do

!debug       call displacement(posb, pos, delr,npart)
       
       call calcforce(NATOMS,type,posb,boxl,forceb,total_energy)
       evalf_number = evalf_number + 1

!       fpar=DOT_PRODUCT(forceb(:),projection(:))
!       if(fpar > 0.0d0 ) projection = -1.0d0 * projection
!       fdotinit =  DOT_PRODUCT(forceb(:), projection(:))


      forceb2_ = DOT_PRODUCT(forceb(:),forceb(:))
       force2_ = DOT_PRODUCT(force_(:),force_(:))
       force_k_ = force_ + (forceb2_ / force2_) * perp_force_;
      fdotinit = DOT_PRODUCT(force_k_(:),projection(:))


       
       parl_force(:)  = fdotinit * projection(:)
       perp_force_(:)  = force_k_(:) - parl_force(:)
       
       fperp2=DOT_PRODUCT(perp_force_(:),perp_force_(:))
       fperp = sqrt(fperp2)
       
       force_ = forceb

       current_energy = total_energy
       current_fperp = fperp
       iperp = iperp + 1
       
!       if(fperp2 < FTHRESH2   .or.                &
!            iperp > (iter - 2)  .or.                &
!            iperp > MAXIPERP    .or.                &
!            itry > 50) then
!         pos(:) = posb(:)
!         ret=-30000
!         exit
!     end if
!   end do      ! relaxation in the hyperplane   
    

    ! we check if the saddle point is reached

    do i=1,NATOMS
       fparl_ia(i)=dsqrt(parl_force(i)**2+        &
                        parl_force(i+NATOMS)**2+ &
			parl_force(i+2*NATOMS)**2)
      
      fperp_ia(i)=dsqrt(perp_force_(i)**2+        &
                        perp_force_(i+NATOMS)**2+ &
			perp_force_(i+2*NATOMS)**2)
   end do
   fparl_max=MAXVAL(fparl_ia(:)) 
   fperp_max=MAXVAL(fperp_ia(:))
   
   call displacement(posref, pos, delr,npart)

   if (print_details .and. mod(iter,iprint) .eq. 0) then 
       write(FLOG,"(' ','iter: ',i4,'  iper: ',i4,'  ener: ',f10.4,'  fpar: ',f10.4,'  fperp: ',f10.4,&
            & '  e-val: ', f10.4,'  delr: ',f10.4,'  npart: ',i4,'  evalf: ',i6)")  &
            & iter,  iperp, total_energy, fparl_max, fperp_max, eigenvalue, delr, npart, evalf_number

       write(*,"(' ','iter: ',i4,'  iper: ',i4,'  ener: ',f10.4,'  fpar: ',f10.4,'  fperp: ',f10.4,&
            & '  e-val: ', f10.4,'  delr: ',f10.4,'  npart: ',i4,'  evalf: ',i6)")  &
            & iter, iperp, total_energy, fparl_max, fperp_max, eigenvalue, delr, npart, evalf_number
    endif

    saddle_energy = current_energy
!    write(*,*) fparl_max+fperp_max, EXITTHRESH
    if (ic_eigenvalue==0) then
       if ( (fparl_max+fperp_max)< EXITTHRESH)  then
          ret = 20000 + iter
          exit
       else if ( (abs(fparl_max)<0.1*EXITTHRESH) .and. (fperp_max<EXITTHRESH) ) then
          ret = 10000 + iter
          exit
       endif
    end if
    


!? what is that ... releted to the old fashion to put the convergency.
!   force = force_
    fpar = DOT_PRODUCT(force(:),projection(:)) 	
    if(fpar > 0.0d0 ) projection = -1.0d0 * projection
    fdotinit =  DOT_PRODUCT(force(:), projection(:))

    ! We now move the configuration along the eigendirection corresponding
    ! to the lowest eigenvalue
    
!    if(abs(fpar > EXITTHRESH*0.4 ) then 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! BEGIN COSMIN SECTION !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    if( fparl_max > EXITTHRESH*6.0 ) then 
!       pos = pos - sign(one_art,fpar) * 0.8d0 * INCREMENT * projection / sqrt(1.0d0*iter) 
!    else
!       pos = pos - sign(one_art,fpar) * 0.3d0 * INCREMENT * projection / sqrt(1.0d0*iter) 
!    endif
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! END COSMIN SECTION !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! BEGIN KIMIYA SECTION !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    if( fparl_max > EXITTHRESH*6.0 ) then 
        pos = pos + (fdotinit/eigenvalue) * projection 
!    else
!        pos = pos - sign(one_art,fpar) * 0.3d0 * INCREMENT * projection / sqrt(1.0d0*iter) 
!    endif
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! END KIMIYA SECTION !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    if (abs(fparl_max) < 0.8*EXITTHRESH) then
       !
       it_art=0
    !
    else
       !
       it_art=0
       new_projection = .false.    ! We start from the previuos direction each time
       maxvec = NVECTOR_LANCZOS
       call lanczos(maxvec,new_projection)
       if (reject) eigen_rejected = eigen_rejected + 1  
       !
    end if
    
    ! if we have more than 2 consecutive positive eigenvalue 
    ! we leave this region
    
    if (eigenvalue < 0.0 ) then
       ic_eigenvalue=0 
    else 
       ic_eigenvalue=ic_eigenvalue+1
    end if
    if (ic_eigenvalue > 1 ) then
       write(*,*)   'THIS IS NOT A SADDLE POINT ic_eigenvalue', ic_eigenvalue
       write(FLOG,*) 'THIS IS NOT A SADDLE POINT ic_eigenvalue', ic_eigenvalue
      ret = 60000 + iter
      exit
   end if
   
   it_art=0
   call calcforce(NATOMS,type,pos,boxl,force,total_energy)
   evalf_number = evalf_number + 1
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! BEGIN KIMIYA SECTION !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       forceb2 = DOT_PRODUCT(force(:),force(:))
!        force2 = DOT_PRODUCT(force(:),force(:))
       force_k = force + (forceb2 / force2) * perp_force_
!       force_k = force + (forceb2 / forceb2_) * perp_force_
       fdotinit = DOT_PRODUCT(force_k(:),projection(:))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  END  KIMIYA SECTION !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   fpar=DOT_PRODUCT(force(:),projection(:))
   perp_force(:)  = force_k(:) - fdotinit * projection(:)  ! Vectorial force
 
!   perp_force = perp_force_;
  
   fperp2=DOT_PRODUCT(perp_force(:),perp_force(:))
   fperp = sqrt(fperp2)
   
   current_fperp = fperp
   
end do
end subroutine saddle_converge
