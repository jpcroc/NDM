module calfoberend_mod
!  USE tempinst_mod,only: tempinst
  USE T_kind_param_m, ONLY:  double
    USE var_pot, ONLY:gamlt,cm
  USE gen_com_m, ONLY:bk,pi,text,tstep,tautcon,text
  implicit none
contains
  subroutine calfoberend(im,xp, vp, fp,ityp)

    integer::im
    real(double)  :: xp(3,im)
    real(double)  :: vp(3,im)
    real(double)  :: fp(3,im)
    integer  :: ityp(im)
    integer :: i,ic
    !real(double), external :: tempinst
    real(double) :: gamb,fact,tempm1,mv2,v2

    do i = 1,im
       v2= vp(1,i)**2+ vp(2,i)**2+ vp(3,i)**2
       mv2= mv2 + cm(ityp(i))*v2
    enddo

#ifdef PARA
    call MPI_ALLREDUCE(mv2,mv2_glob,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    mv2 = mv2_glob
    tempm1=mv2/(3.d0*float(im_glob)*bk)

#else
    tempm1=mv2/(3.d0*float(im)*bk)


#endif


!    tempm1=tempinst(im,vp,ityp)

    !      write(6,*)'jy suis'
    gamb=1./(2.*tauTcon)
    !      write(6,*)gamb,text,tempm1

    do i=1,im
       fact=cm(ityp(i))*gamb*(Text/tempm1-1.0)
       do ic=1,3
          !            write(6,*)fp(ic,i),fact*vp(ic,i)
          fp(ic,i)=fp(ic,i)+fact*vp(ic,i)
       end do

    end do

  end subroutine calfoberend

  subroutine dynlangevin(im,xp, vp, fp,ityp,il,Gl)
    USE gen_com_m, ONLY:
    USE var_pot, ONLY:
#ifdef PARA
    USE mod_para
#endif
    integer::im
    real(double)  :: xp(3,im)
    real(double)  :: vp(3,im)
    real(double)  :: fp(3,im)
    real(double)  :: Gl(3,im)
    integer  :: ityp(im)
    integer::il
    real(double)::rga
    integer :: i,ic
    real(double) :: u1,u2
    !langevin code partir du poly de Gabriel Stolz page 84, dans une version avec expoentielle comme Manuel et Cosmin
    select case (il)
    case(1)
       !     write(6,*)'rga',rga,Gl(ic,i)*sqrt(cm(ityp(1))*bk*text*(1-rga**2))/cm(ityp(1)),vp(1,1)
       do i=1,im
          rga=exp(-gamlt(ityp(i))*tstep/2)
          !        write(6,'(A,2G15.7)')'gamstd ',gamlt(ityp(i)),rga
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



  end subroutine dynlangevin


end module calfoberend_mod

