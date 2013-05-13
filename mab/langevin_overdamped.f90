
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
    USE mab_in_ndm_module, only: sig_ll,m_i,it_mab,dtlang,Ecinetique,xbar,  &
                                 abf_type,block,damp_coef

    implicit none
    integer     :: ic,it_langevin
    real(double) :: xbari(3),xbarc(3)
    real(double), dimension(6,im+1)::gau
   
    it_langevin=it_mab 
   
    call genere_bruit(gau)


    do ic=1,3
       xbar(ic)=sum(xp(ic,1:im))/dble(im) ! barycentre sur les particules
       xbari(ic)=xbar(ic)
    enddo

 
 !one force calculation ....

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
        
    
      xp(1:3,1:im)= xp(1:3,1:im)+fp(1:3,1:im)/(damp_coef*m_i(1:3,1:im))*dtlang + sig_ll(1:3,1:im)*gau(1:3,1:im)
    do ic=1,3
       xbar(ic)    = sum(xp(ic,1:im))/dble(im)
       xbarc(ic)   = xbar(ic) - xbari(ic) ! déplacement du barycentre 
       xp(ic,1:im) = xp(ic,1:im) - xbarc(ic)  ! on recentre tout le systeme
    enddo
  
   
    return
 end subroutine langevin_overdamped


