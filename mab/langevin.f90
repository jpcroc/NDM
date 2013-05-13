 subroutine prepare_langevin
 use gen_com_m, only:one,two,im,imm,tstep
 use var_pot
 use tab_imm_m
 use mab_in_ndm_module, only: sig_i,sig_ll, rga_i,m_i,temperature,damp_coef,langevin_type,dtlang
 implicit none
 real(double) :: gamma

  select case (langevin_type)
    case (1)
     sig_ll(1:3,1:im) = sqrt(2.d0*temperature*dtlang/(damp_coef*m_i(1:3,1:im)))
    case (2) 
     gamma=one/(tstep*1d2)
     rga_i(1:3,1:im) = exp(-gamma*tstep/two)
     sig_i(1:3,1:im) = sqrt(m_i(1:3,1:im)*temperature*(one-rga_i(1:3,1:im)**2))
  end select

 return  
 end subroutine prepare_langevin


 subroutine langevin()
    !On input ... the parameters 
    !dt -        integration step size (units, internal units ndm ?)
    !temperature (units ? )
    !rga         (?)

    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    use gen_com_m
    use tab_imm_m
    USE mab_in_ndm_module, only: sig_i,rga_i,m_i,it_mab,dtlang,Ecinetique,xbar,  &
                                 abf_type,block

    implicit none
    integer     :: ic,it_langevin
    real(double) :: pp(3,im)
    real(double) :: Ecin4
    real(double) :: vbar(3)
    real(double) :: xbari(3),xbarc(3)
    real(double), dimension(6,im+1)::gau
   
    it_langevin=it_mab 
   
    call genere_bruit2(sig_i,gau)

    Ecin4 = zero

 
    pp(1:3,1:im)=vp(1:3,1:im)*m_i(1:3,1:im)

    do ic=1,3
       vbar(ic)=sum(pp(ic,1:im))/dble(im) ! barycentre sur les particules
       xbar(ic)=sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules
       xbari(ic)=xbar(ic)
    enddo
! one force calculation ....

         if (itab/=0) then
          if (mod(it_langevin,itab)==0) then
           call caltabt
          endif
         endif
         if (ltabvois.and.mod(it_langevin,itetabvois)==0) call caltabi
        call calfo
        if (block) call calfoblock()

        call fill_histo()
        
        select case (abf_type)
           case (1)
                  continue
           case (2)  
                  call calfo_ABF_BIN
          ! case (3) 
          !        call calfo_ABF_GAUSSIAN ! ABF Gaussian
          ! case (4) 
          !        call calfo_ABF_EE
        end select 
        
    
    !step1: from p(1) -> p(1+1/4)
    pp(1:3,1:im)=pp(1:3,1:im)*rga_i(1:3,1:im) + gau(1:3,1:im)
    do ic=1,3
       vbar(ic)=sum(pp(ic,1:im))/dble(im) ! vit barycentre sur les particules
       pp(ic,1:im)=pp(ic,1:im)-vbar(ic)
    enddo
    !call control_angular_momenta(pp,xp)
    vp(1:3,1:im)=pp(1:3,1:im)/m_i(1:3,1:im)

    !step2: from p(1+1/4) -> p(1+1/2)
    pp(1:3,1:im) =  pp(1:3,1:im) + (fp(1:3,1:im))*dtlang/two
    
    !step3: from x(1) -> x(1+1)
    xp(1:3,1:im)=  xp(1:3,1:im) + pp(1:3,1:im)*dtlang/m_i(1:3,1:im)
    

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/dble(im)
       xbarc(ic)   = xbar(ic) - xbari(ic) ! déplacement du barycentre 
       xp(ic,1:im) = xp(ic,1:im) - xbarc(ic)  ! on recentre tout le systeme
    enddo

     !recompute the forces

      if (itab/=0) then
       if (mod(it_langevin,itab)==0) then
        call caltabt
       endif
      endif
      if (ltabvois.and.mod(it_langevin,itetabvois)==0) call caltabi
     call calfo
     if (block) call calfoblock()

     call fill_histo()
     
     select case (abf_type)
        case (1)
               continue
        case (2)  
               call calfo_ABF_BIN
       ! case (3) 
       !        call calfo_ABF_GAUSSIAN ! ABF Gaussian
       ! case (4) 
       !        call calfo_ABF_EE
     end select 


    !step4: p(1+1/2) -> p(1+3/4) 
    pp(1:3,1:im) = pp(1:3,1:im) + (fp(1:3,1:im))*dtlang/two
    vp(1:3,1:im) = pp(1:3,1:im) / m_i(1:3,1:im)
    
    
    !step5: p(1+3/4) -> p(1+1)
    pp(1:3,1:im) = pp(1:3,1:im)*rga_i(1:3,1:im) + gau(4:6,1:im)

    do ic=1,3
       vbar(ic)=sum(pp(ic,1:im))/dble(im) ! vit barycentre sur les particules
       pp(ic,1:im)=pp(ic,1:im)-vbar(ic)
    enddo
    !call control_angular_momenta(pp,xp)
    vp(1:3,1:im) = pp(1:3,1:im)/m_i(1:3,1:im)

    do ic=1,3
      Ecin4 = Ecin4 + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo

    Ecinetique = Ecin4
   
    return

 end subroutine langevin


