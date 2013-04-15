 subroutine langevin(dt,temperature,rga)!,ielat,iwmax, ityp)!)!(xp, vp, fp,dt)!, ielat,iwmax, ityp)
    !On input ... the parameters 
    !dt -        integration step size (units ?)
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
    real(double) :: pp(3,im)
    real(double), dimension (3,N)::sig
    real(double) :: xbar(3)
    real(double) :: Ecin0, Ecin1, Ecin3, Ecin4
    real(double) :: dt
    real(double) :: vbar(3)
    real(double), dimension(6,N+1)::gau
    real(double) ::temperature,rga
    sig(:,:) = sqrt(m_i(:,:)*temperature*(one-rga**2))
    call genere_bruit2(sig,gau)
    sig_i(:,:)=sig(:,:)
    rga_i(:,:)=rga

    Ecin0 = zero
    Ecin1 = zero
    Ecin3 = zero
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

     
     do ic =1,3
       Ecin0 = Ecin0   + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
     enddo
    pp(1:3,1:im)=pp(1:3,1:im)*rga_i(1:3,1:im) + gau(1:3,1:im)
    do ic=1,3
       vbar(ic)=sum(pp(ic,1:im))/dble(im) ! vit barycentre sur les particules
       pp(ic,1:im)=pp(ic,1:im)-vbar(ic)
    enddo
    !call control_angular_momenta(pp,xp)
    vp(1:3,1:im)=pp(1:3,1:im)/m_i(1:3,1:im)


    do ic =1,3
       Ecin1 = Ecin1   + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo
    
    pp(1:3,1:im) =  pp(1:3,1:im) + (fp(1:3,1:im))*dt/two
    
    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules
    enddo
    xp(1:3,1:im)=  xp(1:3,1:im) + pp(1:3,1:im)*dt/m_i(1:3,1:im)
    

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/dble(im)-xbar(ic) ! déplacement du barycentre 
       xp(ic,1:im) = xp(ic,1:im) - xbar(ic)  ! on recentre tout le systeme
    enddo

	  if (itab/=0) then
	    if (mod(it_langevin,itab)==0) then
	       call caltabt
	    endif
	 endif
	 if (ltabvois.and.mod(it_langevin,itetabvois)==0) call caltabi
	 call calfo
  
!debugC    write(*,*) 'Langevin.........:', itab, ltabvois, itetabvois, it_langevin
    
    pp(1:3,1:im) = pp(1:3,1:im) + (fp(1:3,1:im))*dt/two
    vp(1:3,1:im) = pp(1:3,1:im) / m_i(1:3,1:im)
    
    
    do ic=1,3
       Ecin3 = Ecin3 + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo

    pp(1:3,1:im) = pp(1:3,1:im)*rga_i(1:3,1:im) + gau(4:6,1:im)
    do ic=1,3
       vbar(ic)=sum(pp(ic,1:im))/dble(im) ! vit barycentre sur les particules
       pp(ic,1:im)=pp(ic,1:im)-vbar(ic)
    enddo
    !prot(1:3,1:im)=pp(1:3,1:im)
    call control_angular_momenta(pp,xp)!,prot)
    !pp(1:3,1:im) = prot(1:3,1:im)
!    call control_angular_momenta(pp,xp)
    vp(1:3,1:im) = pp(1:3,1:im)/m_i(1:3,1:im)

    do ic=1,3
      Ecin4 = Ecin4 + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
    enddo

    Ecinetique = Ecin4
   
    return

 end subroutine langevin


