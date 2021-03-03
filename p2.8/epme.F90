module epme_mod
!  USE temp_com,only:volu,bg ! A EFFACER
  USE moduli_mod,only: moduli
    USE gen_com_m, ONLY:it,itesigma,pi,potis3,zero
        implicit none 
        contains
!                   Version du 10/12/2001
! ***********************************************************
subroutine epme (Deb,Fin,sige,im,xp,fp,ityp,volu,bg)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m
  use moduli_mod,only:bspline

  USE var_pot, ONLY:alpha,auxe,maxorder,kpmex,kpmey,kpmez,ncoucx,ncoucy,ncoucz,nf1,nf2,nf3,nff,nfft1,nfft2,nfft3,&
       &npoint,pterm,volterm,fr1,fr2,fr3,iiim,q,bsmod3,iiim,ijim,ikim,bsmod2,de3,bsmod1,de2,de1,tabv3
  USE fft_com_m

  implicit none

  !#ifdef para2c
  !     include '/usr/include/DXMLDEF.FOR' ! FFT para de ixia
  !#endif


  real(double), dimension(3,3) :: sige,bg
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !----------------------------------------------
  integer,intent(in)::im
  real(double), dimension(im) :: scalar
  integer :: deb, fin
  real(double),intent(inout),allocatable::fp(:,:)
  real(double),intent(in)::xp(:,:),volu
  integer,intent(in),allocatable::ityp(:)
    
  real(double) :: potisewg, hbn2  
  real(double)  :: theta1(maxorder,Deb:Fin),dtheta1(maxorder,Deb:Fin)
  real(double)  :: theta2(maxorder,Deb:Fin),dtheta2(maxorder,Deb:Fin)
  real(double)  :: theta3(maxorder,Deb:Fin),dtheta3(maxorder,Deb:Fin)
  real(double)  :: w1pme(Deb:Fin),w2pme(Deb:Fin),w3pme(Deb:Fin)

  real(double)  :: xi, yi, zi, w, hsq, term, expterm
  real(double)  :: h1, h2, h3, r1, r2, r3, t1, t2, t3
  integer :: i0i,i0,j0i,j0,k0i,k0,it1,it2,it3
  integer :: i, j, k, m, k1, k2, k3, m1, m2, m3
  integer :: mr1, mr2, mr3
  real(double)  :: e, denom, eterm, produc, struc2
  real(double)  :: dn1, dn2, dn3, dt1, dt2, dt3
  logical :: lvect
  integer :: nbatom


  allocate (qgrid(kpmex,kpmey,kpmez))

  ! *** Initialisation
  !          fp = zero !Test
  nbatom=Fin-Deb+1

  qgrid(:,:,:) = cmplx(0., 0.)

  theta1(:maxorder,Deb:Fin)=zero
  dtheta1(:maxorder,Deb:Fin)=zero
  theta2(:maxorder,Deb:Fin)=zero
  dtheta2(:maxorder,Deb:Fin)=zero
  theta3(:maxorder,Deb:Fin)=zero
  dtheta3(:maxorder,Deb:Fin)=zero

  !     write(6,*)'deb fin ',deb,fin !TestJM
  do i=Deb,Fin

     xi=xp(1,i)
     yi=xp(2,i)
     zi=xp(3,i)
     w=xi*bg(1,1)+yi*bg(1,2)+zi*bg(1,3)
     fr1(i)=kpmex*(w-dble(nint(w))+0.5d0)
     w=xi*bg(2,1)+yi*bg(2,2)+zi*bg(2,3)
     fr2(i)=kpmey*(w-dble(nint(w))+0.5d0)
     w=xi*bg(3,1)+yi*bg(3,2)+zi*bg(3,3)
     fr3(i)=kpmez*(w-dble(nint(w))+0.5d0)

     w1pme(i)=fr1(i)-dble(int(fr1(i)))
     w2pme(i)=fr2(i)-dble(int(fr2(i)))
     w3pme(i)=fr3(i)-dble(int(fr3(i)))

  enddo

  call bspline (nbatom,w1pme,maxorder,theta1,dtheta1)
  call bspline (nbatom,w2pme,maxorder,theta2,dtheta2)
  call bspline (nbatom,w3pme,maxorder,theta3,dtheta3)

  ! *** Remplissage du tableau Q(k1,k2,k3) de la reference Essmann ***
  do it3=1,maxorder
     do m=Deb,Fin
        k0i=int(fr3(m))-maxorder
        k0=k0i+it3
        ikim(it3,m)=k0+1+(kpmez-sign(kpmez,k0))/2
        j0i=int(fr2(m))-maxorder
        j0=j0i+it3
        ijim(it3,m)=j0+1+(kpmey-sign(kpmey,j0))/2
        i0i=int(fr1(m))-maxorder
        i0=i0i+it3
        iiim(it3,m)=i0+1+(kpmex-sign(kpmex,i0))/2
     enddo
  enddo


  do m=Deb,Fin
     do it3=1,maxorder
        t3=theta3(it3,m)
        k=ikim(it3,m)
        do it2=1,maxorder
           t2=theta2(it2,m)
           j=ijim(it2,m)
           do it1=1,maxorder
              t1=theta1(it1,m)
              i=iiim(it1,m)
              produc=t3*t2*t1*q(ityp(m))
              !            write(6,*)i,j,k
              !      write(6,*)kpmex,kpmey,kpmez
              qgrid(i,j,k)=qgrid(i,j,k)+produc

           enddo
        enddo
     enddo
  enddo


  ! *** Calcul de la transformee de Fourier discrete de qgrid() ***


  ! appel de la fft codees dans fft_inter

  !DEBUG_WITH_INTEL
  !>call fft_inter('cald') !Calcul avec la FFT sequentielle
  !     call fftfront (nfft1,nfft2,nfft3,kpmex,kpmey,kpmez,kpme,&
  !     ntable,table,qgrid)


  ! *** FFT vectorisee Fujitsu***
  !      call DFTCBM (ARPME,AIPME,NDIM,IDIM,WORK,TWORK,MODE,INIT,IFAIL)


  e=zero




  do i=1,npoint-1

     k3=i/nff+1
     j=i-(k3-1)*nff
     k2=j/nfft1+1
     k1=j-(k2-1)*nfft1+1
     m1=k1-1
     m2=k2-1
     m3=k3-1

     if (k1.gt.nf1) m1=m1-nfft1
     if (k2.gt.nf2) m2=m2-nfft2
     if (k3.gt.nf3) m3=m3-nfft3
     r1=dble(m1)
     r2=dble(m2)
     r3=dble(m3)
     mr1=abs(m1)
     mr2=abs(m2)
     mr3=abs(m3)
     lvect=.FALSE.
     if (mr1.le.ncoucx.and.mr2.le.ncoucy.and.mr3.le.ncoucz) lvect=.TRUE.
     if (mr1.eq.0.and.mr2.eq.0.and.mr3.eq.0) lvect=.FALSE.
     if (lvect) then
        h1=bg(1,1)*r1+bg(2,1)*r2+bg(3,1)*r3
        h2=bg(1,2)*r1+bg(2,2)*r2+bg(3,2)*r3
        h3=bg(1,3)*r1+bg(2,3)*r2+bg(3,3)*r3
        hsq=h1*h1+h2*h2+h3*h3
        term=-pterm*hsq
        !         if (term.gt.-50.0d0) then
        denom=volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
        expterm=exp(term)/denom
        struc2=real(qgrid(k1,k2,k3))**2+imag(qgrid(k1,k2,k3))**2
        eterm=0.5d0*expterm*struc2
        e=e+eterm
        !         endif
        qgrid(k1,k2,k3)=expterm*qgrid(k1,k2,k3)

        if (itesigma>0) then
           if (mod(it,itesigma)==0) then
              potisewg = tabv3(m1,m2,m3)*struc2/bsmod1(k1)/bsmod2(k2)/bsmod3(k3)
              hbn2 = (h1**2+h2**2+h3**2)*4*pi*pi

              sige(1,1) = sige(1,1)+potisewg*h1*h1/hbn2*(hbn2/(4.0*alpha&
                   **2)+1)/volu*4*pi*pi
              sige(2,2) = sige(2,2)+potisewg*h2*h2/hbn2*(hbn2/(4.0*alpha&
                   **2)+1)/volu*4*pi*pi
              sige(3,3) = sige(3,3)+potisewg*h3*h3/hbn2*(hbn2/(4.0*alpha&
                   **2)+1)/volu*4*pi*pi
           endif
        endif

     else
        qgrid(k1,k2,k3)=DCMPLX(zero,zero)
     endif
  enddo

  e=e*auxe ! Rajout Conversion en cgs et prise en compte de e*e


  potis3=e

  ! *** Debut du calcul des forces par la methode pme ***

  !DEBUG_WITH_INTEL
  !call fft_inter ('cali')  !FFT sequentielle inverse
  !     call fftback (nfft1,nfft2,nfft3,kpmex,kpmey,kpmez,kpme,ntable,table,qgrid)
  !      MODE='M'
  !        call DFTCBM (ARPME,AIPME,NDIM,IDIM,WORK,TWORK,MODE,INIT,IFAIL)
  !
  !     get first derivatives of the reciprocal space energy
  !
  dn1 = dble(nfft1)
  dn2 = dble(nfft2)
  dn3 = dble(nfft3)
  do m=Deb,Fin
     de1(m) = zero
     de2(m) = zero
     de3(m) = zero
  enddo
  do it3 = 1, maxorder
     do it2 = 1, maxorder
        do it1 = 1, maxorder
           do m = Deb,Fin
              i=iiim(it1,m)
              j=ijim(it2,m)
              k=ikim(it3,m)
              t3 = theta3(it3,m)
              dt3 = dn3 * dtheta3(it3,m)
              t2 = theta2(it2,m)
              dt2 = dn2 * dtheta2(it2,m)
              t1 = theta1(it1,m)
              dt1 = dn1 * dtheta1(it1,m)
              term = dreal(qgrid(i,j,k))
              de1(m) = de1(m) + term*dt1*t2*t3
              de2(m) = de2(m) + term*dt2*t1*t3
              de3(m) = de3(m) + term*dt3*t1*t2
           enddo
        enddo
     enddo
  enddo

  do m=Deb,Fin
     fp(1,m)=fp(1,m)-q(ityp(m))*(bg(1,1)*de1(m)+bg(2,1) &
          *de2(m)+bg(3,1)*de3(m))*auxe
     fp(2,m)=fp(2,m)-q(ityp(m))*(bg(1,2)*de1(m)+bg(2,2) &
          *de2(m)+bg(3,2)*de3(m))*auxe
     fp(3,m)=fp(3,m)-q(ityp(m))*(bg(1,3)*de1(m)+bg(2,3) &
          *de2(m)+bg(3,3)*de3(m))*auxe
  enddo

  deallocate (qgrid)
  return
end subroutine epme
end module
