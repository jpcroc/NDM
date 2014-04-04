
! *********************************************************************
subroutine bruit_xp 
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
 implicit none
   integer    :: ia, ip
   integer, dimension(2) :: iseedt
   real(double)  :: zr1,zr2,zr3,zr4,totalbruit
 
  call system_clock (iseed)
   iseedt(1)=iseed
 
  call random_seed(iseedt(1))
 totalbruit=0.d0
 bruitmd(1:3,1:im)=0.d0
  do ia=1,im
   call random_number(zr1)
   call random_number(zr2)
   call random_number(zr3)
   call random_number(zr4)
    if(zr1.eq.0.d0) zr1=0.000000001d0
    if(zr2.eq.0.d0) zr2=0.000000001d0
    if(zr3.eq.0.d0) zr3=0.000000001d0
    if(zr4.eq.0.d0) zr4=0.000000001d0
    
    bruitmd(1,ia)=sqrt((-log(zr1)))*cos(2.0*pi*zr3)
    bruitmd(2,ia)=sqrt((-log(zr1)))*sin(2.0*pi*zr3)
    bruitmd(3,ia)=sqrt((-log(zr2)))*cos(2.0*pi*zr4)
    totalbruit=totalbruit + bruitmd(1,ia)**2 + bruitmd(2,ia)**2 + bruitmd(3,ia)**2
   end do
  
  if (rang==0)  write(6,*) 'ISEED for MD, NORM of the noise ',iseed, neb_noise_scale, totalbruit
  bruitmd(1:3,1:im) = bruitmd(1:3,1:im) * mdcg_noise_scale * xp(1:3,1:im) / (sqrt(totalbruit))
  
  end subroutine bruit_xp




! *********************************************************************
subroutine initspeed
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif
  ! *********************************************************************


  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, ic, ia, ib
  integer, dimension(:), allocatable :: iseedt
  real(double), dimension(ntyp) :: temptyp
  ! ym      real(double), dimension(nce) :: tempc
!  real(double), dimension(noxyz) :: tempc
  real(double) :: vv, v0, v1, z1, z2, z3, z4
  real(double) :: tempsauv ,totmass
  real(double) :: rx, ry, rz, r2x, r2y, r2z, r2
  real(double) :: prx, pry, prz, px, py, pz, vrx, vry, vrz,ka
  real(double) :: omegax, omegay, omegaz
  real(double), dimension(3) :: vt1,  scom, pav,kinx
  real(double), dimension(3,3) :: ainer, aineri
  real(double), dimension(3,ntyp) :: vav
  ! ym      real(double), dimension(3,nce) :: sigkinec
  real(double), dimension(3,noxyz) :: sigkinec
  integer  :: i_glob
  integer  :: est_local
  integer :: seed_size
  integer::iti
  real(double)::sd,grnd,theta,fhi

#if(PARA)
  real(double) :: kinx_glob
  real(double), dimension(3)   :: scom_glob, pav_glob
  real(double), dimension(3,3) :: ainer_glob
  real(double) :: totmass_glob
  real(double) :: prx_glob
  real(double) :: pry_glob
  real(double) :: prz_glob
#endif
  !-----------------------------------------------
  !  external fucntions
  !-----------------------------------------------
  real(double) :: tempinst

  !  if (rang==0) write(6,*) 'PARA-T entree initspeed',iseed
  select case (dmtype)
  case(3,30,5,11,7)
     vp = 0.0
      return
  case(2)
     if (lvpread) then 
        return
     else
       if (mdcg_noise==0) then 
         vp=0.0
         return
       else
         vp=0.0 
         call bruit_xp
         xp(1:3,1:im) = xp(1:3,1:im) + bruitmd(1:3,1:im)
       return
       end if
     end if
  end select
1 continue
  !      write(6,*)'vp',vp(1,1)
  if (lvpread) then
     !       oldtstep=1.0d-15
     tempsauv=tempinst(vp,ityp)
     if (rang==0) write(6,*)'tempsauv ',tempsauv

     xpp(:,:im) = xp(:,:im)-(xp(:,:im)-xpp(:,:im))*tstep/oldtstep
!     vp(:,:im)=vp(:,:im)*tstep/oldtstep

     if (tinit<=0) then
        ! velocities are read from file and not modified
        if (rang==0) write (6, *) 'pas de chgt des vitesses= '
!        return
     else
        ! velocities are read from file and rescaled
        ! MPI
        if (rang==0) write (6, *) 'scaling read velocities at TINIT = ', &
             tinit, 'K'
        !   tempsauv=tempinst(vp,ityp)
        !   if (rang==0) write(6,*)'tempsauv ',tempsauv
        if (tempsauv.le.1.) then
           lvpread=.false. ; goto 1
        end if
        vv = sqrt(tinit/tempsauv)
        xpp(:,:im) = xp(:,:im)-(xp(:,:im)-xpp(:,:im))*vv
        vp(:,:im) = vp(:,:im)*vv
     endif

  else


     if (tinit<=0) then
        ! velocities are not read and no starting temperature is given
        if (dmtype==1 .or. dmtype==4) then                  !DM run
           write (6, *) rang,'no way to initiate the velocities stop'
           call arret_ndm
        else                                 !quench run
           vp(1,:im) = 0.0
           vp(2,:im) = 0.0
           vp(3,:im) = 0.0
           xpp(1,:im) = xp(1,:im)
           xpp(2,:im) = xp(2,:im)
           xpp(3,:im) = xp(3,:im)
        endif
     else
        !  a starting temperature is given
        ! MPI

        if (rang==0) write (6, *) 'random velocities at TINIT = ', tinit, &
             'K'

        if (iseed==0)  iseed=1
        
        call random_seed(size=seed_size)
        allocate(iseedt(seed_size))
        iseedt = 0
        

!        if (iseed==0)  call system_clock (iseed) 

!        write(6,*)''proc', myid, iseed pour tirage des vitesses',iseed

        iseedt(1)=iseed
        call    random_seed (put=iseedt)
        deallocate(iseedt)

        v0 = sqrt(2.D0*bk*tinit)
        vt1(:)=0.0
        !        do i_glob = 1, im_glob
        do i = 1, im
           !******************************************
           call random_number(z1)
           call random_number(z2)
           call random_number(z3)
           call random_number(z4)
           if(z1.eq.0.d0) z1=0.000000001d0
           if(z2.eq.0.d0) z2=0.000000001d0
           if(z3.eq.0.d0) z3=0.000000001d0
           if(z4.eq.0.d0) z4=0.000000001d0

 	   est_local=0
     !	   do i=1,im
     !	     if (num_at_glob(i)==i_glob) then
           est_local=1
           !		exit
           !             endif
           !           enddo

	   if (est_local==1) then
	      if(z1.eq.0.d0) z1=0.000000001d0
              if(z2.eq.0.d0) z2=0.000000001d0
              if(z3.eq.0.d0) z3=0.000000001d0
              if(z4.eq.0.d0) z4=0.000000001d0

              !**************************************
              v1 = one/sqrt(cm(ityp(i)))
!              vp(1,i) = v1*v0*sqrt((-log(z1)))*cos(2.0*pi*z3)
!              vp(2,i) = v1*v0*sqrt((-log(z1)))*sin(2.0*pi*z3)
!              vp(3,i) = v1*v0*sqrt((-log(z2)))*cos(2.0*pi*z4)
              theta=acos(1-2*z3)
              fhi=2*pi*z4
              vp(1,i) = v1*v0*sqrt((-log(z1)))*sin(theta)*cos(fhi)
              vp(2,i) = v1*v0*sqrt((-log(z1)))*sin(theta)*sin(fhi)
              vp(3,i) = v1*v0*sqrt((-log(z2)))*cos(theta)
              

	   endif
        end do

        do ic=1,3
           do i=1,im
              kinx(ic)=kinx(ic)+0.5*vp(ic,i)*vp(ic,i)*cm(ityp(i))/dfloat(im_glob)
           end do
           ka=0.5*bk*tinit 
#if(PARA)
           call MPI_ALLREDUCE(kinx(ic),kinx_glob,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
           kinx(ic)=kinx_glob
#endif          
           if (rang==0) write(6,*)'dir ',ic,' ka Ktinit ', kinx(ic),ka
        end do
        !***************************************************************
        !       Make total momentum zero
        !         Calculate the global average momemtum & center-of-mass


        totmass = 0.d0
        scom = 0.d0
        pav  = 0.d0

        do i = 1, im
           ic = ityp(i)
           totmass = totmass+cm(ic)
           do ia = 1, 3
              scom(ia) = scom(ia)+xp(ia,i)*cm(ic)
              pav (ia) = pav(ia) +vp(ia,i)*cm(ic)
           enddo
        enddo

#if(PARA)
        call MPI_ALLREDUCE(totmass,  totmass_glob,  1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(scom(1:3),scom_glob(1:3),3,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(pav(1:3), pav_glob(1:3), 3,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
	totmass = totmass_glob
        scom = scom_glob
        pav  = pav_glob
#endif          

        do ia = 1, 3
           scom(ia) = scom(ia)/totmass
           pav (ia) = pav (ia)/float(im_glob)
        enddo

        if (rang==0) then
           write (6,*) 'center of mass = ',scom(1)*1d8,&
                scom(2)*1d8,scom(3)*1d8
           write (6,*) 'Momentum/atom  = ',pav(1),pav(2),pav(3)
        endif

        !         Shift velocities to make the total momemtum zero
        do ic = 1, ntyp
           do ia = 1, 3
              vav(ia,ic) = pav(ia)/cm(ic)
           enddo
        enddo

        do i = 1, im
           do ia = 1, 3
              vp(ia,i) = vp(ia,i)-vav(ia,ityp(i))
           enddo
        enddo

        if (lcalcjq) then
           ka=0.5*bk*tinit
           kinx=0.
           do ic=1,3
              do i=1,im
                 kinx(ic)=kinx(ic)+0.5*vp(ic,i)*vp(ic,i)*cm(ityp(i))/dfloat(im_glob)
              end do

#if(PARA)
              call MPI_ALLREDUCE(kinx(ic),kinx_glob,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
              kinx(ic)=kinx_glob
#endif          

              if (rang==0)write(6,*)'dir ',ic,' ka Ktinit ', kinx(ic),ka
              do i=1,im
                 vp(ic,i)= vp(ic,i)*dsqrt(ka/kinx(ic))
              end do
           end do
           !          kinx=0.
           !          do ic=1,3
           !             do i=1,im
           !                kinx(ic)=kinx(ic)+0.5*vp(ic,i)*vp(ic,i)*cm(ityp(i))/dfloat(im)
           !             end do

           !            write(6,*)'dir ',ic,' ka Ktinit ', kinx(ic),ka

           !          end do


        else



           !      Make the angular momentum zero
           !         Calculate inertia tensor & angular momentum
           ainer = 0.d0
           prx = 0.d0
           pry = 0.d0
           prz = 0.d0

           do i = 1, im
              ic = ityp(i)
              rx = xp(1,i)-scom(1)
              ry = xp(2,i)-scom(2)
              rz = xp(3,i)-scom(3)
              r2x = rx*rx
              r2y = ry*ry
              r2z = rz*rz
              r2  = r2x+r2y+r2z
              ainer(1,1) = ainer(1,1)+cm(ic)*(r2-r2x)
              ainer(2,2) = ainer(2,2)+cm(ic)*(r2-r2y)
              ainer(3,3) = ainer(3,3)+cm(ic)*(r2-r2z)
              ainer(2,3) = ainer(2,3)-cm(ic)*ry*rz
              ainer(3,1) = ainer(3,1)-cm(ic)*rz*rx
              ainer(1,2) = ainer(1,2)-cm(ic)*rx*ry
              px  = cm(ic)*vp(1,i)
              py  = cm(ic)*vp(2,i)
              pz  = cm(ic)*vp(3,i)
              prx = prx+ry*pz-rz*py
              pry = pry+rz*px-rx*pz
              prz = prz+rx*py-ry*px
           enddo

           ainer(3,2) = ainer(2,3)
           ainer(1,3) = ainer(3,1)
           ainer(2,1) = ainer(1,2)

#if(PARA)
           call MPI_ALLREDUCE(ainer(1:3,1:3), ainer_glob(1:3,1:3), 9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
           ainer = ainer_glob
           call MPI_ALLREDUCE(prx, prx_glob, 1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
           prx = prx_glob
           call MPI_ALLREDUCE(pry, pry_glob, 1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
           pry = pry_glob
           call MPI_ALLREDUCE(prz, prz_glob, 1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
           prz = prz_glob
#endif          

           if (rang==0) then
              write(6,997) (ainer(1,ib),ib=1,3),prx
              write(6,997) (ainer(2,ib),ib=1,3),pry
              write(6,997) (ainer(3,ib),ib=1,3),prz
           endif
997        format('Inertia/anglm = ',3e12.4,5x,e12.4)

           !         calculate  angular velocity
           call matinv(ainer,aineri)
           omegax = aineri(1,1)*prx+aineri(1,2)*pry+aineri(1,3)*prz
           omegay = aineri(2,1)*prx+aineri(2,2)*pry+aineri(2,3)*prz
           omegaz = aineri(3,1)*prx+aineri(3,2)*pry+aineri(3,3)*prz

           !         shift velocities to make the angular momentum zero
           do i = 1, im
              rx = xp(1,i)-scom(1)
              ry = xp(2,i)-scom(2)
              rz = xp(3,i)-scom(3)
              vrx = omegay*rz-omegaz*ry
              vry = omegaz*rx-omegax*rz
              vrz = omegax*ry-omegay*rx
              vp(1,i) = vp(1,i)-vrx
              vp(2,i) = vp(2,i)-vry
              vp(3,i) = vp(3,i)-vrz
           enddo

           !      Old velocities translation
           !            vp(1,:im) = vp(1,:im)-vt1(1)/float(im)
           !            vp(2,:im) = vp(2,:im)-vt1(2)/float(im)
           !            vp(3,:im) = vp(3,:im)-vt1(3)/float(im)

           !      Initialisation of the previous position for Verlet

         end if

        xpp(1,:im) = xp(1,:im)-vp(1,:im)*tstep
        xpp(2,:im) = xp(2,:im)-vp(2,:im)*tstep
        xpp(3,:im) = xp(3,:im)-vp(3,:im)*tstep

      endif

   endif
  !     write(6,*)'sortie initspeed'

  tempsauv=tempinst(vp,ityp)
  if (rang==0) write(6,*)'temperature fin initspeed ',tempsauv

  if (lfrozen.EQV..true.) then
          WHERE (frozen(:,1:im))
                  vp(:,1:im) = 0.d0
                  xpp(:,1:im) = xp(:,1:im)
          END WHERE
  end if


  if (ldeplainit==.true.)then

     if (rang==0) then
        write(6,*)
        write(6,*)'-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
        write(6,*)'depla init Tempdeplainit',tempdeplainit,'debyetemp= ',debyetemp
        
        do iti=1,ntyp
           sd=sqrt((3*tempdeplainit*hbar**2)/(bk*cm(iti)*debyetemp**2))
           write(6,*)'sd de iti',sd,iti
        end do
     end if
     do i=1,im
        sd= sqrt((3*tempdeplainit*hbar**2)/(bk*cm(ityp(i))*debyetemp**2))
        
        do ic=1,3
           call gaussianrand(grnd)
!           write(6,*)grnd
           xp(ic,i)=xp(ic,i)+sd*grnd
        end do
     end do
     if (lperiod==.true.) call period
  end if


  return


    

end subroutine initspeed

subroutine gaussianrand(gr)
  USE T_kind_param_m
  use gen_com_m,only:pi
  implicit none
  real(double),intent(out)::gr
  
  real(double):: v1,v2,r,fac,z1,z2

1 continue
  call random_number(z1)  
  call random_number(z2)  
  gr=sqrt(-2*log(z1))*cos(2*pi*z2)

  
end subroutine  gaussianrand
