subroutine calfoberend(xp, vp, fp,ityp)
  use gen_com_m
  use var_pot
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

subroutine calfolangevin(xp, vp, fp,ityp,il,Gl)
  use gen_com_m
  use var_pot
#if(PARA)
  use mod_mpi
#endif

  real(double)  :: xp(3,imm)
  real(double)  :: vp(3,imm)
  real(double)  :: fp(3,imm)
  real(double)  :: Gl(3,imm)
  integer  :: ityp(imm)
  integer::il
  real(double)::rga
  integer :: i,ic
  real(double) :: u1,u2
!langevin codé à partir du poly de Gabriel Stolz page 84, dans une version avec expoentielle comme Manuel et Cosmin
  select case (il)
  case(1)
!     write(6,*)'rga',rga,Gl(ic,i)*sqrt(cm(ityp(1))*bk*text*(1-rga**2))/cm(ityp(1)),vp(1,1)
     do i=1,im
        rga=exp(-gamlt(ityp(i))*tstep/2)
        do ic=1,3
!  write(6,*)'ct',cm(ityp(1)),tstep
           call random_number(u1)
           call random_number(u2)
           Gl(ic,i)=sqrt(-2.*log(u1))*cos(2.*pi*u2)   
           vp(ic,i) = vp(ic,i)*rga+ fp(ic,i)*tstep/(cm(ityp(i))*2)+Gl(ic,i)*sqrt(cm(ityp(i))*bk*text*(1-rga))/cm(ityp(i))

        end do
     end do
  case(2)

     do i=1,im
        rga=exp(-gamlt(ityp(i))*tstep/2)
        do ic=1,3
           vp(ic,i) = vp(ic,i)*rga+ fp(ic,i)*tstep/(cm(ityp(i))*2)+Gl(ic,i)*sqrt(cm(ityp(i))*bk*text*(1-rga))/cm(ityp(i))
        end do
     end do
  case default 
     write(6,*)'check ilangevin'
     stop
  end select



end subroutine calfolangevin




