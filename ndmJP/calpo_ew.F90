subroutine calpo_ew
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m
  use gen_com_m
  use var_pot
  implicit none
  real(double) :: pi2, fact, fact1, fact2, hk2, ex, ex1, ex2 ,hbv(3)
  integer ::nb1,nb2,nb3
     ! --- Tableaux des troisiemes termes de la sommation d'Ewald ---
     pi2 = pi*pi
     !         volu = zl(1)*zl(2)*zl(3)
     fact = pi2/alpha**2
     fact1 = auxe/2./pi/volu
     fact2 = auxe*2./volu

     do nb1 = -ncoucx, ncoucx
        do nb2 = -ncoucy, ncoucy
           do nb3 = -ncoucz, ncoucz
              if (nb1==0.and.nb2==0.and.nb3==0) cycle
              hbv(:) = bg(:,1)*nb1+bg(:,2)*nb2+bg(:,3)*nb3
              hk2 =hbv(1)**2+hbv(2)**2+hbv(3)**2
              ex = exp((-hk2*fact))/hk2
              ex1 = ex*fact1
              ex2 = ex*fact2
              tabv3(nb1,nb2,nb3) = ex1
              tabf3(:,nb1,nb2,nb3) = ex2*q(:)
           end do
        end do
     end do



  ! Traitement du cas PME

  if (iewald==2) then

     ntable=4*kpme+15
     call DynamicalAllocationPME   ! Allocation dynamique de memoire
     nfft1=kpmex  ! Nombre de points de la grille
     nfft2=kpmey
     nfft3=kpmez
     npoint=nfft1*nfft2*nfft3

     !                       ndim1=2*(nfft1/2)+1
     !                       ndim2=2*(nfft2/2)+1
     !                       ndim3=2*(nfft3/2)+1

     nff=nfft1*nfft2
     nf1=(nfft1+1)/2
     nf2=(nfft2+1)/2
     nf3=(nfft3+1)/2

     pterm=(pi/alpha)**2
     volterm=pi*volu

     call moduli !Initialisation des tableaux bsmod1, bsmod2, bsmod3

     table(:ntable,:3)=zero
     call fft_inter('init')


  endif !Fin du cas iewald=2
end subroutine calpo_ew
