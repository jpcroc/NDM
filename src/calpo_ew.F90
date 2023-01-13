module calpo_ew_mod
  USE moduli_mod,only: moduli
  USE gen_com_m, ONLY:pi,zero
  USE var_pot, ONLY:alpha,auxe,iewald,kpme,kpmex,kpmey,kpmez,ncoucx,ncoucy,ncoucz,nf1,&
       &nf2,nf3,nff,nfft1,nfft2,nfft3,npoint,ntable,pterm,volterm,table,q,&
       &nvecttot,tabf3,tabv3
  use boxconfig,only:box_config
    USE recips_mod,only: calcvol
  implicit none
contains
  subroutine calpo_ew(boxndm,immT)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m

    USE dynallocPME
    implicit none
    type(box_config)::boxndm
    integer,intent(in)::immT
    real(double) :: pi2, fact, fact1, fact2, hk2, ex, ex1, ex2 ,hbv(3)
    integer ::nb1,nb2,nb3
    ! --- Tableaux des troisiemes termes de la sommation d'Ewald ---
    pi2 = pi*pi
    !         volu = zl(1)*zl(2)*zl(3)
    fact = pi2/alpha**2
    boxndm%volu=calcvol(boxndm%at(1:3,1),boxndm%at(1:3,2),boxndm%at(1:3,3))
    fact1 = auxe/2./pi/boxndm%volu
    fact2 = auxe*2./boxndm%volu
!    write(6,*)'KPME',kpme
    do nb1 = -ncoucx, ncoucx
       do nb2 = -ncoucy, ncoucy
          do nb3 = -ncoucz, ncoucz
             if (nb1==0.and.nb2==0.and.nb3==0) cycle
             hbv(:) = boxndm%bg(:,1)*nb1+boxndm%bg(:,2)*nb2+boxndm%bg(:,3)*nb3
             hk2 =hbv(1)**2+hbv(2)**2+hbv(3)**2
             ex = exp((-hk2*fact))/hk2
             ex1 = ex*fact1
             ex2 = ex*fact2
             tabv3(nb1,nb2,nb3) = ex1
             tabf3(:,nb1,nb2,nb3) = ex2*q(:)
!             write(6,*)nb1,nb2,nb3,ex1,ex2,q
          end do
       end do
    end do

!!$    nv=0
!!$    if (.not.(allocated(nb1v))) then 
!!$       allocate (nb1v(nvecttot));  allocate (nb2v(nvecttot));allocate (nb3v(nvecttot))
!!$    end if
!!$       ! repartition des vecteurs du RRec.
!!$       do nb1 = -ncoucx, ncoucx
!!$          do nb2 = -ncoucy, ncoucy
!!$             do nb3 = -ncoucz, ncoucz
!!$                if (nb2==0.and.nb3==0.and.nb1==0) cycle
!!$                nv=nv+1
!!$                nb1v(nv)=nb1; nb2v(nv)=nb2; nb3v(nv)=nb3
!!$             enddo
!!$          enddo
!!$       enddo
!!$       if (nv.ne.nvecttot) stop
!!$

    ! Traitement du cas PME

    if (iewald==2) then

       ntable=4*kpme+15
       call DynamicalAllocationPME(immT)   ! Allocation dynamique de memoire
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
       volterm=pi*boxndm%volu

       call moduli !Initialisation des tableaux bsmod1, bsmod2, bsmod3

       table(:ntable,:3)=zero
       !DEBUG_WITH_INTEL
       !call fft_inter('init')


    endif !Fin du cas iewald=2
  end subroutine calpo_ew
end module calpo_ew_mod
