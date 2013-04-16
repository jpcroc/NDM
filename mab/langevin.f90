 subroutine prepare_langevin
 implicit none
 use gen_com_m, only:one,two
 use var_pot
 use tab_imm_m

  m_i(1:3,1:im) = cm(ityp(i))
  
  gamma=one/(tstep*1d2)

  rga_i(1:3,1:im) = exp(-gamma*tstep/two)
  sig_i(1:3,1:im) = sqrt(m_i(1:3,1:im)*temperature*(one-rga_i(1:3,1:im)**2))





 end subroutine prepare_lagevin


 subroutine langevin(dt,temperature,rga_i,sig_i)!,ielat,iwmax, ityp)!)!(xp, vp, fp,dt)!, ielat,iwmax, ityp)
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

    implicit none
    integer ic
    real(double),intent(in) :: dt,temperature
    real(double), dimension (3,N),intent(in) ::sig_i,rga_i
    real(double) :: xbar(3)
    real(double) :: pp(3,im)
    real(double) :: Ecin4
    real(double) :: vbar(3)
    real(double), dimension(6,N+1)::gau
    
    call genere_bruit2(sig_i,gau_i)
    !sig_i(:,:)=sig(:,:)
    !rga_i(:,:)=rga

    Ecin4 = zero

    write(*,*) im 
 
    pp(1:3,1:im)=vp(1:3,1:im)*m_i(1:3,1:im)

    do ic=1,3
       vbar(ic)=sum(vp(ic,1:im))/dble(im) ! barycentre sur les particules
    enddo

         if (itab/=0) then
          if (mod(it_langevin,itab)==0) then
           call caltabt
          endif
         endif
         if (ltabvois.and.mod(it_langevin,itetabvois)==0) call caltabi
        call calfo
    
    !step1: from p(1) -> p(1+1/4)
    pp(1:3,1:im)=pp(1:3,1:im)*rga_i(1:3,1:im) + gau(1:3,1:im)
    do ic=1,3
       vbar(ic)=sum(pp(ic,1:im))/dble(im) ! vit barycentre sur les particules
       pp(ic,1:im)=pp(ic,1:im)-vbar(ic)
    enddo
    !call control_angular_momenta(pp,xp)
    vp(1:3,1:im)=pp(1:3,1:im)/m_i(1:3,1:im)

    !step2: from p(1+1/4) -> p(1+1/2)
    pp(1:3,1:im) =  pp(1:3,1:im) + (fp(1:3,1:im))*dt/two
    
    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules
    enddo
    !step3: from x(1) -> x(1+1)
    xp(1:3,1:im)=  xp(1:3,1:im) + pp(1:3,1:im)*dt/m_i(1:3,1:im)
    

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/dble(im)-xbar(ic) ! déplacement du barycentre 
       xp(ic,1:im) = xp(ic,1:im) - xbar(ic)  ! on recentre tout le systeme
    enddo

     !recompute the forces
     if (itab/=0) then
       if (mod(it_langevin,itab)==0) then
        call caltabt
     endif
    endif
    if (ltabvois.and.mod(it_langevin,itetabvois)==0) call caltabi
    call calfo
  
    !step4: p(1+1/2) -> p(1+3/4) 
    pp(1:3,1:im) = pp(1:3,1:im) + (fp(1:3,1:im))*dt/two
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


