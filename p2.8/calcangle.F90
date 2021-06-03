module calcangle_mod
  USE T_kind_param_m, ONLY:  double
  USE temp_com,only:ncel,noxyz,at,bg,deltadist,atincel,nato !A EFFACER
  USE notperiod_mod,only: notperiod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE var_pot, ONLY:ntyp,ty
  use atomconfig,only: atom_config
  use boxconfig,only:box_config
  USE cellconfig,only:cell_config, caltabtC
  USE gen_com_m, ONLY:lperiod,rang,it,pi,timel,lspacendm
#ifdef PARA
    USE mpi
    USE Tpara,only:COMM_space,nprocspace

#endif
 
  implicit none

  type adf_typ
     real(double)::rcangle,thetamax,thetamin
     real(double),allocatable::fda(:,:,:,:)
     integer::contmax,nfda
     logical::linstantfda
     character::nom(2)
     integer,allocatable::nad(:)
  end type adf_typ

  type(adf_typ)::adf0


contains

  subroutine initadf(adfc,contmax,rcangle,linstantfda,nom,thetamax,thetamin)
    type(adf_typ),intent(out)::adfc
    integer,intent(in)::contmax
    real(double),intent(in)::rcangle,thetamax,thetamin
    logical,intent(in)::linstantfda
    character,intent(in)::nom(2)
    allocate(adfc%nad(ntyp))

    allocate(adfc%fda(ntyp,ntyp,ntyp,contmax))
    adfc%nfda=0
    adfc%contmax=contmax
    adfc%rcangle=rcangle
    adfc%linstantfda=linstantfda
    adfc%nom=nom
    adfc%thetamax=thetamax
    adfc%thetamin=thetamin

  end subroutine initadf


  subroutine calcangle(atadf,celadf,boxadf,adfc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
#ifdef PARA
    USE mpi
    USE Tpara,only:COMM_space,nprocspace

#endif

    !******************************************************************
    implicit none
    integer :: i1, i2, i3, i4, koo, ko1, ko2, i, iti1,&
         j, k, m,m1,ka,ma,m2
    real(double) :: thetaijk, c11, c21, c31, c12, c22, c32, &
         incre,dij2,dik2,dij,dik
    real(double) :: costheta,invincre,rspace2,cv(1,3)
    real(double),allocatable::rc2(:,:),rc22(:,:)
    !-----------------------------------------------
    class(atom_config),intent(in)::atadf
    type(box_config),intent(in)::boxadf
    type(cell_config),intent(in):: celadf
    type(adf_typ)::adfc
    real(double),allocatable::fdatemp(:,:,:,:)
    real(double),allocatable :: xpnp(:,:)
    allocate (xpnp(3,atadf%imm))
    allocate(rc2(ntyp,ntyp),rc22(ntyp,ntyp))
    allocate(fdatemp(ntyp,ntyp,ntyp,adfc%contmax))
    fdatemp=0
    !repartition des atomes entre les petites cel.
    !  if (rang==0) write(6,*) 'PARA-T entree calcangle'
    rc2=adfc%rcangle
    rc22=rc2**2
    adfc%nfda=adfc%nfda+1
    incre = (adfc%thetamax-adfc%thetamin)/adfc%contmax
    invincre = 1/incre

    if (lperiod) then
       xpnp(:,:)=atadf%xp(:,:)
    else 
       call notperiod(atadf%imm,atadf%xp,xpnp,boxadf%at,boxadf%bg)
    end if
    adfc%nfda = adfc%nfda+1

    do i = 1, atadf%im-1
       koo = atadf%ielat(i)
       do i1 = 0, 26
          ko1 = celadf%ncel(koo,i1)
          do i2 = 1, celadf%nato(ko1) 
             j = celadf%atincel(i2,ko1)
             c11 = xpnp(1,i)-xpnp(1,j) 
             c21 = xpnp(2,i)-xpnp(2,j) 
             c31 = xpnp(3,i)-xpnp(3,j)
             c11 = c11+sum(boxadf%at(1,:)*celadf%deltadist(:,i1,koo))
             c21 = c21+sum(boxadf%at(2,:)*celadf%deltadist(:,i1,koo))
             c31 = c31+sum(boxadf%at(3,:)*celadf%deltadist(:,i1,koo))

             if (noxyz.ne.1) then
                if (abs(c11)>rc2(atadf%ityp(i),atadf%ityp(j))) cycle
                if (abs(c21)>rc2(atadf%ityp(i),atadf%ityp(j))) cycle
                if (abs(c31)>rc2(atadf%ityp(i),atadf%ityp(j))) cycle
             else
                cv(1,1) = c11
                cv(1,2) = c21
                cv(1,3) = c31
                call cryst_to_cart (1, cv, boxadf%bg, -1) !cryst vers cart sur cv
                WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                   cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
                END WHERE
                call cryst_to_cart (1, cv, boxadf%at, 1) !cryst vers cart sur cv
                c11=cv(1,1)
                c21=cv(1,2)
                c31=cv(1,3)
             end if


             dij2=c11**2+c21**2+c31**2
             if(dij2>rc22(atadf%ityp(i),atadf%ityp(j))) cycle
             do i3 = 0,26
                ko2 = ncel(koo,i3)
                do i4 = 1, celadf%nato(ko2)
                   k = celadf%atincel(i4,ko2)
                   if(j==i.or.k==i.or.j==k) cycle                 
                   c12 = xpnp(1,i)-xpnp(1,k) 
                   c22 = xpnp(2,i)-xpnp(2,k) 
                   c32 = xpnp(3,i)-xpnp(3,k)
                   c12 = c12+sum(at(1,:)*celadf%deltadist(:,i1,koo))
                   c22 = c22+sum(at(2,:)*celadf%deltadist(:,i1,koo))
                   c32 = c32+sum(at(3,:)*celadf%deltadist(:,i1,koo))
                   if (noxyz.ne.1) then
                      if (abs(c12)>rc2(atadf%ityp(i),atadf%ityp(k))) cycle
                      if (abs(c22)>rc2(atadf%ityp(i),atadf%ityp(k))) cycle
                      if (abs(c32)>rc2(atadf%ityp(i),atadf%ityp(k))) cycle
                   else
                      cv(1,1) = c12
                      cv(1,2) = c22
                      cv(1,3) = c32
                      call cryst_to_cart (1, cv, boxadf%bg, -1) !cryst vers cart sur cv
                      WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                         cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
                      END WHERE
                      call cryst_to_cart (1, cv, boxadf%at, 1) !cryst vers cart sur cv
                      c12=cv(1,1)
                      c22=cv(1,2)
                      c32=cv(1,3)
                   end if

                   !                    if (abs(c12)>rc2(ityp(i),ityp(k))) cycle
                   !                    if (abs(c22)>rc2(ityp(i),ityp(k))) cycle
                   !                    if (abs(c32)>rc2(ityp(i),ityp(k))) cycle
                   if(c12**2+c22**2+c32**2>rc22(atadf%ityp(i),atadf%ityp(k))) cycle

                   costheta = (c11*c12+c21*c22+c31*c32)/ &
                        (sqrt(c11*c11+c21*c21+c31*c31)*sqrt(c12*c12+ &
                        c22*c22+c32*c32))

                   thetaijk = acos(costheta)
                   ka = int(thetaijk*invincre)
                   m = ka+1
                   fdatemp(atadf%ityp(j),atadf%ityp(i),atadf%ityp(k),m) = fdatemp(atadf%ityp(j),atadf%ityp(i), &
                        atadf%ityp(k),m)+1

                end do
             end do
          end do
       end do
    end do
   do iti1=1,ntyp
       adfc%nad(iti1)=count(atadf%ityp(:)==iti1)
    end do
#ifdef PARA
    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call comm_space%sum(fdatemp)
       call comm_space%sum(adfc%nad)
    end if
#endif
    adfc%fda=adfc%fda+fdatemp
    deallocate (xpnp)
    return
  end subroutine calcangle

  !******************************************************************
  subroutine adfT(adfc)
    implicit none
    type(adf_typ)::adfc
    integer :: i1, i2, i3, i4, koo, ko1, ko2, i, &
         j, k, m,m1,lutriplet,ka,ma,m2

    real(double) :: thetaijk, c11, c21, c31, c12, c22, c32, &
         incre
    real(double) :: rspace2, invincre, &
         costheta
    real(double) :: aaa,bbb,ccc
    real(double), dimension(ntyp,ntyp,ntyp) :: intfda,angle
    character :: triplet1*20,triplet2*20,triplet3*20
    integer :: lenftriplet1,lenftriplet2,lenftriplet3
    character :: ftriplet1*80, ftriplet2*80, ftriplet3*80
    character :: ftripletangle*80,ftrip*80
    integer :: lusauvfda
    character :: charsauvfda*8
    !-----------------------------------------------

    incre = (adfc%thetamax-adfc%thetamin)/adfc%contmax
    invincre = 1/incre

    lutriplet = 13
    if (rang==0) then
       write (6, *)
       write(6,*) '-----------------------------------------'
       write (6, *) '------- Calcul des Distributions angulaires --------'
       write(6,*)

       write (6, '(A,I5,A,D10.3)') '*  ITERATION  = ', it, '  time = ', timel
    end if

    ! On calcule maintenant la distribution angulaire moyennee

    if(.not.adfc%linstantfda) then
       do i3=1,ntyp
          if(adfc%nad(i3)==0) cycle

          do i1=1,ntyp
             if(adfc%nad(i1)==0) cycle

             do i2=i3,ntyp
                if(adfc%nad(i2)==0) cycle
                triplet1= ty(i1)
                triplet2= ty(i2)
                triplet3= ty(i3)
                ftriplet1=triplet1
                ftriplet2=triplet2
                ftriplet3=triplet3
                lenftriplet1=index(ftriplet1,' ')-1
                lenftriplet2=index(ftriplet2,' ')-1
                lenftriplet3=index(ftriplet3,' ')-1
                ftripletangle = ftriplet2(1:lenftriplet2)// &
                     '_'//ftriplet1(1:lenftriplet1)// &
                     '_'//ftriplet3(1:lenftriplet3)//'.angle'

                ftrip = ftriplet2(1:lenftriplet2)// &
                     '_'//ftriplet1(1:lenftriplet1)// &
                     '_'//ftriplet3(1:lenftriplet3)


                do m1=int(adfc%thetamin*invincre),adfc%contmax
                   m2 = m1+1
                   rspace2 = (m2*incre)**2
                   aaa = adfc%fda(i2,i1,i3,m2)/adfc%nfda
                   intfda(i2,i1,i3)=intfda(i2,i1,i3)+aaa
                   angle(i2,i1,i3)=angle(i2,i1,i3)+aaa*m1*incre

                   if(rang==0) then
                      write (13,*) m2*incre*180/pi, aaa
                   endif
                end do


                if(angle(i2,i1,i3).gt.0.1) then
                   if (rang==0) open(lutriplet, file = ftripletangle, status = 'unknown')
                   angle(i2,i1,i3)=angle(i2,i1,i3)/intfda(i2,i1,i3)
                   if(rang.eq.0) write(6,'("Angle moyen du triplet",1X,A8,"=",1x,f7.3)')& 
                        ftrip,angle(i2,i1,i3)*180/pi
                   do m1=int(adfc%thetamin*invincre),adfc%contmax
                      m2 = m1+1
                      rspace2 = (m2*incre)**2
                      aaa = adfc%fda(i2,i1,i3,m2)/adfc%nfda
                      if(rang==0) then
                         write (13,*) m2*incre*180/pi, aaa
                      endif
                   end do

                endif
             end do    !i2
          end do
       end do

       if(rang==0)         close(lutriplet)
       if(rang==0) then
          write(6,*)'--------------------------------------'
          write(6,*) '--------------------------------------'
       end if
    else
       if (rang==0) then
          open(unit=32,file='tampfda',form='formatted',status='unknown')

          if (it<=9) write (32, '(I1)') it
          if (it<=99.and.it>9) write (32, 200) it
          if (it<=999.and.it>99) write (32, 300) it
          if (it<=9999.and.it>999) write (32, 400) it
          if (it<=99999.and.it>9999) write (32, 500) it
          if (it<=999999.and.it>99999) write(32, 600) it
          if (it<=9999999.and.it>999999) write(32, 700) it
          if (it<=99999999.and.it>9999999) write(32, 800) it
          if (it<=999999999.and.it>99999999) write(32, 900) it
          if  (it>999999999) then
             write (6, *) 'probleme de format dans calcangle.f90'
             stop
          endif
          rewind 32

          read (32, 1000) charsauvfda

          lusauvfda=31
       end if
       if (rang==0) then
          write(6,*) ' sauvegarde fda it=',it
          write(6,*)
       end if
       do i3=1,ntyp
          if(adfc%nad(i3)==0) cycle
          do i1=1,ntyp
             if(adfc%nad(i1)==0) cycle
             do i2=1,ntyp
                if(adfc%nad(i2)==0) cycle

                triplet1= ty(i1)
                triplet2= ty(i2)
                triplet3= ty(i3)
                ftriplet1=triplet1
                ftriplet2=triplet2
                ftriplet3=triplet3
                lenftriplet1=index(ftriplet1,' ')-1
                lenftriplet2=index(ftriplet2,' ')-1
                lenftriplet3=index(ftriplet3,' ')-1
                ftripletangle = ftriplet2(1:lenftriplet2)// &
                     '_'//ftriplet1(1:lenftriplet1)// &
                     '_'//ftriplet3(1:lenftriplet3)//'.'//charsauvfda
                ftrip = ftriplet2(1:lenftriplet2)// &
                     '_'//ftriplet1(1:lenftriplet1)// &
                     '_'//ftriplet3(1:lenftriplet3)

                if (rang==0) open(unit=lusauvfda, file = ftripletangle, status ='unknown')

                do m1=int(adfc%thetamin*invincre),adfc%contmax
                   m2 = m1+1
                   rspace2 = (m2*incre)**2
                   intfda(i2,i1,i3)=intfda(i2,i1,i3)+adfc%fda(i2,i1,i3,m2)
                   angle(i2,i1,i3)=angle(i2,i1,i3)+adfc%fda(i2,i1,i3,m2)*&
                        m1*incre
                   if(rang==0) then
                      write (lusauvfda,*) m2*incre*180/pi, adfc%fda(i2,i1,i3,m2)
                   endif
                end do

                angle(i2,i1,i3)=angle(i2,i1,i3)/intfda(i2,i1,i3)

                if(rang==0) then
                   write(6,'("Angle instantane moyen du triplet",1x,A8,"=",1x,f7.3)')&
                        ftrip,angle(i2,i1,i3)*180/pi
                endif
             enddo
          enddo
       enddo

200    format(i2)
300    format(i3)
400    format(i4)
500    format(i5)
600    format(i6)
700    format(i7)
800    format(i8)
900    format(i9)
1000   format(a6)
       if (rang==0)then
          write(6,*)'--------------------------------------'
          write(6,*) '--------------------------------------'
          close(lusauvfda)
          close(32,status='DELETE')
       endif
    endif

    adfc%fda(:,:,:,:)=0.D0
    intfda(:,:,:)=0.D0
    angle(:,:,:)=0.D0
    return
  end subroutine adfT

end module calcangle_mod


