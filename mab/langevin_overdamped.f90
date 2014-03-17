
 subroutine langevin_overdamped()
    !On input ... the parameters 
    !dt -        integration step size (units, internal units ndm ?)
    !temperature (units ? )
    !rga         (?)
    !schema Euler
 
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    use gen_com_m
    use tab_imm_m
    USE mab_in_ndm_module, only: sig_ll,m_i,dtlang,xbar,  &
                             Ecinetique,abf_type,block,gamma,fpeinstein,it_en 


    implicit none
    integer     :: ic
    real(double) :: xbari(3),xbarc(3)
    real(double), dimension(6,im+1)::gau
    real(double):: pp(3,im),vc(3,im)
    real(double) :: psum(3)
    real(double):: sig_mass(3,im) 
    real(double) :: Ecin4


    call genere_bruit(gau)
    sig_mass(1:3,1:im)=sig_ll(1:3,1:im)*gau(1:3,1:im)
    do ic=1,3
       xbar(ic)=sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules
       xbari(ic)=xbar(ic)
    enddo

 
 !one force calculation ....
   if (it_en > 0) then
     call calfo_einstein_solid ()
     fp(:,:) = fpeinstein (:,:)
    else 
      call calfo_mab()
    end if          
    
     ! xp(1:3,1:im)= xp(1:3,1:im)+fp(1:3,1:im)/(gamma*m_i(1:3,1:im))*dtlang + sig_ll(1:3,1:im)*gau(1:3,1:im)

      pp(1:3,1:im)= fp(1:3,1:im)/(gamma*m_i(1:3,1:im))*dtlang + sig_mass(1:3,1:im)

      do ic=1,3
      psum(ic)=sum(pp(ic,1:im))/dble(im)
      pp(ic,1:im) = pp(ic,1:im) -psum(ic)
      enddo

    !call control_angular_momenta(pp,xp)

    xp(1:3,1:im)=pp(1:3,1:im)+xp(1:3,1:im)

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/dble(im)
       xbarc(ic)   = xbar(ic) - xbari(ic) ! déplacement du barycentre 
       xp(ic,1:im) = xp(ic,1:im) - xbarc(ic)  ! on recentre tout le systeme
    enddo

    Ecin4=0.d0
    do ic=1,3
      !write(*,*) dtlang,m_i(1,1),gamma
      vc(ic,1:im)=pp(ic,1:im)*sqrt(m_i(ic,1:im))/dtlang  ! cm/dtlang
      !write(*,*)  DOT_PRODUCT(vc(ic,:),vc(ic,:))/two
      Ecin4 = Ecin4 + DOT_PRODUCT(vc(ic,1:im),vc(ic,1:im))/two
    enddo

    Ecinetique = Ecin4


     !write(35,*) it_mab,(2.d0*Ecinetique)/(KtoERG*3.d0*dble(im))
   
    return
 end subroutine langevin_overdamped




