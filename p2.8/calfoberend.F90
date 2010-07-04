subroutine calfoberend(xp, vp, fp,ityp)
  use gen_com_m
  real(double)  :: xp(3,imm)
  real(double)  :: vp(3,imm)
  real(double)  :: fp(3,imm)
  integer  :: ityp(imm)
  integer :: i,ic
  real(double), external :: tempinst
  real(double) :: gamb,fact,tempm1

  tempm1=tempinst(vp,ityp)

  !      write(6,*)'jy suis'
  gamb=1./(2.*tauTcon)
  !      write(6,*)gamb,text,tempm1

  do i=1,imd
     fact=cm(ityp(i))*gamb*(Text/tempm1-1.0)
     do ic=1,3
        !            write(6,*)fp(ic,i),fact*vp(ic,i)
        fp(ic,i)=fp(ic,i)+fact*vp(ic,i)
     end do

  end do

end subroutine calfoberend

      subroutine calfolangevin(xp, vp, fp,ityp)
      use gen_com_m
#if(PARA)
  use mod_mpi
#endif

      real(double)  :: xp(3,imm)
      real(double)  :: vp(3,imm)
      real(double)  :: fp(3,imm)
      integer  :: ityp(imm)
      integer :: i,ic
      real(double) :: gamb,fact,kbtemp,u1,u2,b1,sombruit(3),bruit(3,imm),sombruitglob(3)

!      stop
      kbtemp=text*bk
!      call random_number(u1)
!      call random_number(u2)
!      b1=sqrt(-2.*log(u1))*cos(2.*pi*u2)   
!      write(6,*) 'fpvm11', fp(1,1),vp(1,1)*cm(ityp(1))*gamlang/tstep,+b1*sqrt(2.0*cm(ityp(1))*gamlang*kbtemp)/tstep

      sombruit=0.
      do i=1,im
         do ic=1,3
            call random_number(u1)
            bruit(ic,i)=2*( u1-0.5)
            sombruit(ic)=sombruit(ic)+bruit(ic,i)
         end do
      end do
#if(PARA)
  call MPI_ALLREDUCE(sombruit,sombruitglob,3,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  sombruit = sombruitglob
#endif

      do ic=1,3
         bruit(ic,:)=bruit(ic,:)-sombruit(ic)/float(im_glob)
      end do

      do i=1,im
         fact=cm(ityp(i))*gamlang
         do ic=1,3
            call random_number(u1)
            select case (ilangevin)
               case(1)
                  fp(ic,i)=fp(ic,i)-fact*vp(ic,i)/(2*tstep)+bruit(ic,i)*sqrt(6.0*fact*kbtemp)/tstep
               case(2)
                  fp(ic,i)=fp(ic,i)-fact*vp(ic,i)/(2*tstep)+bruit(ic,i)*sqrt(6.0*fact*kbtemp*(1-gamlang*0.5))/tstep
               end select
         end do
      end do


    end subroutine calfolangevin




