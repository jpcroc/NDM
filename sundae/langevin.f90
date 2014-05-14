! Subroutines related to dynamics: 
!     - cal_hamilton, 
!     - langevin, 
!     - mapping_V_verlet 



Subroutine cal_hamilton(qx,px,N,ekin,epot) !!!!!!!!!!! MODIFIER AVEC BON POTENTIEL
 USE T_kind_param_m, ONLY:  double
 use gen_com_m
 !use jqmod
 use tab_imm_m
 use sundae_module, ONLY : it_trajectory,m_i

 implicit none
 integer   ::N
 real (double), dimension (1:3,1:N) :: qx   ! vector of position
 real (double), dimension (1:3,1:N) :: px   ! vector of impulsion
 real (double) :: pp(1:3,1:N),ekin,epot

 xp(1:3,1:N)=qx(1:3,1:N)
 pp(1:3,1:N)=px(1:3,1:N) 

 call calfo_teledyn(it_trajectory)

 ekin = SUM((pp(1:3,1:im))**2/m_i(1:3,1:im))/2.00
 epot = potist

! write(*,*) 'potist ' ,potist,ecin,potist+ecin

end subroutine cal_hamilton



subroutine langevin(dt,temperature,rga)!,ielat,iwmax, ityp)!)!(xp, vp, fp,dt)!, ielat,iwmax, ityp)
    USE T_kind_param_m, ONLY:  double
    use gen_com_m
    use tab_imm_m
    use sundae_module, ONLY:N, m_i, Ecinetique,sig_i,rga_i,gau, it_langevin

    implicit none
    integer ic
    real(double) :: pp(3,im)
    real(double), dimension (3,N)::siglocal
    real(double) :: xbar(3)
    real(double) :: Ecin0, Ecin1, Ecin3, Ecin4
    real(double) :: dt
    real(double) :: vbar(3)
    real(double) ::temperature,rga


    siglocal(:,:) = sqrt(m_i(:,:)*temperature*(one-rga**2))
    call genere_bruit2(siglocal,gau)
    sig_i(:,:)=siglocal(:,:)
    rga_i(:,:)=rga

    Ecin0 = zero
    Ecin1 = zero
    Ecin3 = zero
    Ecin4 = zero

    im=N
    
    pp(1:3,1:im)=vp(1:3,1:im)*m_i(1:3,1:im)

    do ic=1,3
       vbar(ic)=sum(pp(ic,1:im))/dble(im) ! barycentre sur les particules
    enddo

    call calfo_teledyn(it_langevin)
!deb	  if (itab/=0) then
!deb	    if (mod(it_langevin,itab)==0) then
!deb	       call caltabt
!deb	    endif
!deb	 endif
!deb	 if (ltabvois.and.mod(it_langevin,itetabvois)==0) call caltabi
!deb	 call calfo

     
     do ic =1,3
       Ecin0 = Ecin0   + DOT_PRODUCT(pp(ic,1:im),vp(ic,1:im))/two
     enddo
    pp(1:3,1:im)=pp(1:3,1:im)*rga_i(1:3,1:im) + gau(1:3,1:im)
    do ic=1,3
       vbar(ic)=sum(pp(ic,1:im))/dble(im) ! vit barycentre sur les particules
       pp(ic,1:im)=pp(ic,1:im)-vbar(ic)
    enddo
    !prot(1:3,1:im)=pp(1:3,1:im)
    call control_angular_momenta(pp,xp)
    !pp(1:3,1:im) = prot(1:3,1:im)
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

    call calfo_teledyn(it_langevin)

!deb	  if (itab/=0) then
!deb	    if (mod(it_langevin,itab)==0) then
!deb	       call caltabt
!deb	    endif
!deb	 endif
!deb	 if (ltabvois.and.mod(it_langevin,itetabvois)==0) call caltabi
!deb	 call calfo
  
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

Subroutine mapping_P_Verlet(qx,px,dt,N,q1s2) !!!!!!!!!!! MODIFIER AVEC BON POTENTIEL
 USE T_kind_param_m, ONLY:  double
    use gen_com_m
    !use jqmod
    use tab_imm_m
    use sundae_module, ONLY: m_i,it_trajectory
 implicit none
 integer :: ic
 integer   ::N
 real (double), dimension (3,1:N) :: qx   ! vector of position
 real (double), dimension (3,1:N) :: px   ! vector of impulsion
 real (double), dimension (3,1:N) :: q1s2 ! vector of intermediate position 
 real (double) ::dt
 real (double) :: pp(3,N),xbar(3),dxbar(3)

 xp(1:3,1:N)=qx(1:3,1:N)
 pp(1:3,1:N)=px(1:3,1:N)
 
 do ic=1,3
 pp(ic,1:im) = pp(ic,1:im) - sum(pp(ic,1:im))/dble(im)
 enddo
 call control_angular_momenta(pp,qx) 

 do ic=1,3
   xbar(ic)    = sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules de même masse
 enddo

 vp(1:3,1:im)  = pp(1:3,1:im)/m_i(1:3,1:im)

 xp(1:3,1:im)  = xp(1:3,1:im) + vp(1:3,1:im)*dt/two

 do ic=1,3
   dxbar(ic)    = sum(xp(ic,1:im))/dble(im) - xbar(ic) ! déplacement du barycentre 
   xp(ic,1:im) = xp(ic,1:im) - dxbar(ic)  ! on recentre tout le systeme
 enddo

 q1s2(1:3,1:N) = xp(1:3,1:N)

 call calfo_teledyn(it_trajectory)
 
!  ecin=  SUM((pp(1:3,1:im) + fp(1:3,1:im)*dt)**2/m_i(1:3,1:im))/2.00
!  write(*,*) 'potist ' ,potist,ecin,potist+ecin
    
 pp(1:3,1:im) =  pp(1:3,1:im) + fp(1:3,1:im)*dt

 do ic=1,3
 pp(ic,1:im) = pp(ic,1:im) - sum(pp(ic,1:im))/dble(im)
 enddo

 call control_angular_momenta(pp,qx) 


 vp(1:3,1:im) =  pp(1:3,1:im)/m_i(1:3,1:im)

 xp(1:3,1:im) =  xp(1:3,1:im) + vp(1:3,1:im)*dt/two

!debug write (*,*) 'lang1', xp (1:3,1),fp(1,1:3)

 do ic=1,3
   dxbar(ic)    = sum(xp(ic,1:im))/dble(im)-xbar(ic) ! déplacement du barycentre 
   xp(ic,1:im)  = xp(ic,1:im) - dxbar(ic)  ! on recentre tout le systeme
 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 qx(1:3,1:N)=xp(1:3,1:N)
 px(1:3,1:N)=pp(1:3,1:N)

! write(*,*) 'two dt ',two,dt  test : dt change bien de signe pour le backward

end subroutine mapping_P_Verlet

