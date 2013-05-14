
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
    USE mab_in_ndm_module, only: sig_ll,m_i,dtlang,Ecinetique,xbar,  &
                                 abf_type,block,gamma

    implicit none
    integer     :: ic
    real(double) :: xbari(3),xbarc(3)
    real(double), dimension(6,im+1)::gau
   
   
    call genere_bruit(gau)


    do ic=1,3
       xbar(ic)=sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules
       xbari(ic)=xbar(ic)
    enddo

 
 !one force calculation ....
      call calfo_mab()         
    
      xp(1:3,1:im)= xp(1:3,1:im)+fp(1:3,1:im)/(gamma*m_i(1:3,1:im))*dtlang + sig_ll(1:3,1:im)*gau(1:3,1:im)

    ! call control_angular_momenta(pp,xp)

    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/dble(im)
       xbarc(ic)   = xbar(ic) - xbari(ic) ! déplacement du barycentre 
       xp(ic,1:im) = xp(ic,1:im) - xbarc(ic)  ! on recentre tout le systeme
    enddo

 
   
    return
 end subroutine langevin_overdamped


