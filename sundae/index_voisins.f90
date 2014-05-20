  subroutine index_premier_voisin ()
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    use gen_com_m, ONLY: im,imm,at,bg
    use tab_imm_m, ONLY: xp
    USE sundae_module , ONLY: ipovois,xpvois,xtransla,xpvoisini
  implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer ic,ivois
    real(double)  :: xp_d(3,imm)
    real(double)  :: xpdist(1:im)
   !-----------------------------------------------
    
    !  calcul des distances par rapport à la lacune
    
    xp_d = xp
    call cryst_to_cart (imm, xp_d, bg, -1)    !cart vers cryst
    WHERE ( (xp_d.GT.0.5d0).OR.(xp_d.LT.-0.5d0) )
       xp_d(1:3,1:im) = xp_d(1:3,1:im) - Dble(Nint(xp_d(1:3,1:im)))
    END WHERE
    
    xp_d = MatMul(at,xp_d)
    
    xpdist(1:im) = xp_d(1,1:im)**2+xp_d(2,1:im)**2+xp_d(3,1:im)**2
    ivois=0
    do ic =1,im
       if (xpdist(ic).le.7e-16)  then 
          write(*,*) 'voisin',ic,' de la lacune ' , xpdist(ic)
          ivois=ivois+1
          ipovois(ivois)=ic
          xpvois(1:3,ivois) = xp(1:3,ic)
       endif
    enddo
    
    xpdist(1:im) = (xp_d(1,1:im)-xp_d(1,1))**2+(xp_d(2,1:im)-xp_d(2,1))**2+(xp_d(3,1:im)-xp_d(3,1))**2
    do ic =1,im
       if ((xpdist(ic).le.7e-16).and.(xpdist(ic).ne.0))  then 
          write(*,*) 'voisin',ic,' de la lacune ' , xpdist(ic)
          ivois=ivois+1
          ipovois(ivois)=ic
          xpvois(1:3,ivois) = xp(1:3,ic)
       endif
    enddo
    do ic =1,15
       write(*,*) 'voisin',ic,' de la lacune ' , xpvois(1:3,ic)
    enddo
    write(*,*)
    ! calcul de la translation du centre de masse relativement à la  lacune
    xtransla(1:3) = xp(1:3,1)/dble(im)
    xpvoisini = xpvois

  end subroutine index_premier_voisin

  subroutine distance_premier_voisin ()
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    use gen_com_m, ONLY: im,imm,at,bg,niteration
    use tab_imm_m, ONLY: xp
    USE sundae_module, ONLY : iteration,ipovois,xpvois,xpvoisini,d2vois,xtransla
   implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer ivois
    real(double)  :: xpvoiscentre(3,15)
   !-----------------------------------------------
    !-----------------------------------------------
    
    !  calcul des distances par rapport à la lacune    
    !imm=N
    do ivois=1,15
              write(*,*) ' ipovois(ivois) ',ipovois(ivois)
       xpvois(1:3,ivois) = xp(1:3,ipovois(ivois))
       !       write(*,*) ' xpvois(1:3,ivois) ',xpvois(1:3,ivois)
    enddo
    
    !    call cryst_to_cart (8, xpvois, bg, -1)
    !    WHERE ( (xpvois.GT.0.5d0).OR.(xpvois.LT.-0.5d0) )
    !       xpvois(1:3,1:8) = xpvois(1:3,1:8) - Dble(Nint(xpvois(1:3,1:8)))
    !    END WHERE
    !    xpvois = MatMul(at,xpvois)
    !    d2vois(1:8) = xpvois(1,1:8)**2+xpvois(2,1:8)**2+xpvois(3,1:8)**2
  
    write(*,*) 'text',niteration 
 
    xpvoiscentre(1:3,1)=xpvoisini(1:3,1) - xtransla(1:3)*dble(iteration)/dble(niteration)*dble(im-1)
    write(*,*) 'text' 
    xpvoiscentre(1,2:15)=xpvoisini(1,2:15) + xtransla(1)*dble(iteration)/dble(niteration)
    xpvoiscentre(2,2:15)=xpvoisini(2,2:15) + xtransla(2)*dble(iteration)/dble(niteration)
    xpvoiscentre(3,2:15)=xpvoisini(3,2:15) + xtransla(3)*dble(iteration)/dble(niteration)
    
    d2vois(1:15) = (xpvois(1,1:15)-xpvoiscentre(1,1:15))**2+(xpvois(2,1:15)-xpvoiscentre(2,1:15))**2+(xpvois(3,1:15)-xpvoiscentre(3,1:15))**2

    return

  end subroutine distance_premier_voisin

