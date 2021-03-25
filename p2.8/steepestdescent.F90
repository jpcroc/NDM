module steepestdescent_mod

  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:rang

  use WGC_mod,only:setV_F,nstep,ndir,betaguess,test_conv,ncalls,lvm

  !  real(double),allocatable,dimension (:,:)::X,R,G,H,F

implicit none

contains

  subroutine steepestdescent(N,R,V,F,lover)
    integer::N
    real(double),dimension(:)::R,F
    real(double)::V
    logical ::lover

    integer::idir
    real(double)::beta,forctot,formax
    real(double),dimension(N)::R0,F0,R1
    real(double)::V0,Vb,Vbs2
    logical ::lok
    beta =betaguess
!    write(6,*)'betainit',beta
    call setV_F (N,R,V,F,lover)
!    write(6,*)'post sVF0',V
!    call test_conv(N,F,lover,V,R)
!    call test_conv(N,F,lover,V,R)
!    write(6,*)'posttconf0',lover,V
    if (lover) then
       write(6,*)'NO NEED TO RELAX'
       return
    end if

    do idir=1,ndir
       write(6,*)
       R0(1:N)=R(1:N)
       F0(1:N)=F(1:N)
       V0=V
!       write(6,*)
!       write(6,*)'**************************'
!       write(6,*)'callmindir',idir,beta,V0
       call mindir(lover,beta,N,R0,V0,F0,R,V,F,lOK)

       write(6,*)'minimization idirection; lOVER;  beta ',idir,lover,beta
       if (lover) then
          write(6,*)
          write(6,*)'*************************************'
          if (lok) then
             write(6,*)'RELAXED AFTER ',idir,' DIRECTIONS and ', NCALLS,' force calculations'
             return
          else
             write(6,*)'NOT RELAXED!!!!!!!'
             call test_conv(N,F,lover,V,R)
          end if
       end if

    end do
    return
  end subroutine steepestdescent

  subroutine conjugategradient(N,R,V,F,lover)
    integer::N
    real(double),dimension(:)::R,F
    real(double)::V
    logical ::lover

    integer::idir,i
    real(double)::beta,forctot,formax,gamma
    real(double),dimension(N)::R0,F0,R1,G,H
    real(double)::V0,Vb,Vbs2,alpha,gigi,xixi
    logical ::lok
    alpha =betaguess
    write(6,*)'betainit',beta
    call setV_F (N,R,V,F,lover)
!    write(6,*)'post sVF0',V
!    call test_conv(N,F,lover,V,R)
!    call test_conv(N,F,lover,V,R)
    if (lover) then
       write(6,*)'NO NEED TO RELAX'
       return
    end if
  
    G=F
    H=F
    do idir=1,ndir
       write(6,*)
       R0(1:N)=R(1:N)
!       H0(1:N)=H(1:N)
       V0=V
!       write(6,*)
!       write(6,*)'**************************'
!       write(6,*)'callmindir',idir,beta,V0
       call mindir(lover,alpha,N,R0,V0,H,R,V,F,lOK)

       write(6,*)'minimization idirection; lOVER;  beta ',idir,lover,alpha
       if (lover) then
          write(6,*)
          write(6,*)'*************************************'
          if (lok) then
             write(6,*)'RELAXED AFTER ',idir,' DIRECTIONS and ', NCALLS,' force calculations'
             return
          else
             write(6,*)'NOT RELAXED!!!!!!!'
             call test_conv(N,F,lover,V,R)
          end if
       end if
       
       xixi=0
       gigi=0
       do i=1,N
          xixi=xixi+F(i)*F(i)
          gigi=gigi+G(i)*G(i)
       end do
       gamma=xixi/gigi
       G(:)=F(:)
       H(:)=G(:)+gamma*H(:)
       
    end do
    return
  end subroutine conjugategradient

  subroutine mindir(lover,beta,N,R0,V0,F0,R,V,F,lOK)

    logical,intent(out)::lover,lok
    real(double)::beta
    integer,intent(in)::N
    real(double),intent(in)::R0(N),V0,F0(N)
    real(double),intent(out)::R(N),V,F(N)

    real(double),dimension(N)::Fp,RBs2,Rbeta,Fa,Fb,Fc,Fbeta,Fbs2,Fbetatest,Rbetatest
    real(double)::normF02,Vbeta,Vbs2,va,vb,vc,Vmin,Vbetatest,betatest

    integer::i,istep
    real(double)::a,b,c,AA,BB
    logical ::ldir
    lOK=.true.
    normF02=SUM(F0(:)**2)

    beta=beta*10
    do i=1,nstep
!       write(6,*)
      beta=beta/10
!       write(6,*)'mindir1 betaS2',beta*0.5
       Rbs2(:)=R0(:)+0.5*beta*F0(:)
       call setV_F (N,Rbs2,Vbs2,Fbs2,lover)
!       write(6,*)'mindir1 betaS2',lover,beta*0.5, Vbs2
!       write(6,*)'DIFF BS2',VBS2-V0
       call checkline(lover,ldir,F0,normF02,N,R,V,F,Rbs2,Vbs2,Fbs2)
       if (lover) then
!          write(6,*)'betaS2 relaxed'
          if (VBS2.GT.V0) then
             V=V0;F=F0;R=R0
!             write(6,*)'STOP BETAS2 '
             lok=.false.
          end if
          return
       end if
       if (ldir) then
!          write(6,*)'line search over: betas2'
          beta=beta/2
          if (VBS2.GT.V0) then
             V=V0;F=F0;R=R0
             lover=.true.
!             write(6,*)'STOP BETAS2 line'
             lok=.false.
          end if

          return
       end if
       if  (VBS2.LT.V0) exit
    end do
!    write(6,*)
!    write(6,*)'BS2 OK-> BETA'
    Rbeta(:)=R0(:)+beta*F0(:)

    do i=1,nstep
!       write(6,*)'mindir2 beta',beta
       call setV_F (N,Rbeta,Vbeta,Fbeta,lover)
!       write(6,*)'mindir1 beta',lover, beta,Vbeta
       call checkline(lover,ldir,F0,normF02,N,R,V,F,Rbeta,Vbeta,Fbeta)
       if (lover) then
!          write(6,*)'beta relaxed'
          if ((VBS2.lt.Vbeta).or.(V0.lt.Vbeta)) then
             if (V0.LT.VBS2) then
                V=V0;R=R0;F=F0
                lover=.true.
!                write(6,*)'STOP BETA'
                lok=.false.
             else
                V=Vbs2;R=Rbs2;F=fbs2
             end if
          endif
          return
       end if
       if (ldir) then
!          write(6,*)'line search over: beta'
          if ((VBS2.lt.Vbeta).or.(V0.lt.Vbeta)) then
             if (V0.LT.VBS2) then
                V=V0;R=R0;F=F0
                lover=.true.
!                write(6,*)'STOP BETALINE'
                lok=.false.
             else
                V=Vbs2;R=Rbs2;F=fbs2
             end if
          endif
          return
       end if
!       write(6,*)'DIFF BETA',VBeta-VBS2
       if (Vbeta.gt.Vbs2) exit
 !   write(6,*)
       beta=beta*2
       Rbs2=Rbeta
       Vbs2=Vbeta
       Fbs2=Fbeta
       Rbeta(:)=R0(:)+beta*F0(:)
    end do
 !   write(6,*)
!   write(6,*)'BETA OK->PARA'

    a=0 ; b=beta ; c=beta/2
    vb=vbeta ; Vmin=Vbs2
    Va=V0
    Vc=Vbs2
 !   write(6,*)'V',Va,Vb,Vc
    do istep=1,nstep
       AA=(vc-va)/((c-a)*(c-b))-(VB-VA)/((b-a)*(c-b))
       BB=(VB-VA)/(b-a) -AA*(b+a)
       betatest=-0.5*BB/AA     
!       betatest=b-0.5*( ((b-a)**2)*(Vb-Vc)-(((b-c)**2)*(Vb-Va)))/ ((b-a)*(vb-vc)-(b-c)*(vb-va))
       Rbetatest(:)=R0(:)+betatest*F0(:)
!      write(6,*)'mindir3 beta',betatest
       call setV_F(N,Rbetatest,Vbetatest,Fbetatest,lover)
 !      write(6,*)'mindir3 beta',betatest,Vbetatest
       If (Vbetatest.gt.Vmin) then
!          write(6,*)'PARABOLIC SERACH FAILURE SWITHING TO BINARY'
          exit ! recherche parabolique en échec
       end If
       call checkline(lover,ldir,F0,normF02,N,R,V,F,Rbetatest,Vbetatest,Fbetatest)
       if (lover) then
!          write(6,*)'PARABOLIC relaxed'
          return
       end if
       if (ldir) then
 !         write(6,*)'PARABOLIC line search over'
          beta=betatest
          if (Vbetatest.GT.V0) then
             V=V0;R=R0;F=F0
             lover=.true.
!             write(6,*)'STOP PARABOLIC'
             lok=.false.
          end if
          return
       end if
       Vmin=Vbetatest
       if (betatest.LT.c) then
          b=c;Vb=Vc
          c=betatest;Vc=Vbetatest
       else
          a=c; Va=Vc
          c=betatest;Vc=Vbetatest
       end if
    end do
    ! debut de la recherche binaire      
 !      write(6,*)
 !      write(6,*)'BINARY START'

    a=0 ; b=beta ; c=beta/2
    vb=vbeta ; Vmin=Vbs2
    Va=V0
    Vc=Vbs2
    Vmin=Vbs2
 !   write(6,*)
 !   write(6,*)'BBIN',a,b,c
!    write(6,*)'VBIN',Va,Vb,Vc

    do istep=1,nstep
    if (Va.gt.Vb) then
          a=c; Va=Vc
       else
          b=c; Vb=Vc
       end if
 !      write(6,*)
 !      write(6,*)'BBIN',a,b
 !      write(6,*)'VBIN',Va,Vb
       beta=(a+b)/2
       c=beta
       Rbeta(:)=R0(:)+beta*F0(:)
       write(6,'(A,E15.5)')'mindir4 beta',beta
       call setV_F(N,Rbeta,Vbeta,Fbeta,lover)
 !      write(6,*)'mindir3 beta',beta,Vbeta
       call checkline(lover,ldir,F0,normF02,N,R,V,F,Rbeta,Vbeta,Fbeta)

       Vc=Vbeta
       if (lover) then
 !         write(6,*)'BINARY relaxed'
          return
       end if
       if (ldir) then
 !         write(6,*)'BINARY line search over'
          if (Vbeta.GT.V0) then
             V=V0;R=R0;F=F0
             lover=.true.
!             write(6,*)'STOP BINARY'
             lok=.false.
          end if

          return
       end if


    end do
    return
  end subroutine mindir


  subroutine checkline (lover,Ldir,F0,normF02,N,R,V,F,Rt,Vt,Ft)
    logical,intent(in)::lover
    logical,intent(out)::ldir
    integer::N
    real(double),dimension(N)::F0,R,Ft,Rt,F
    real(double)::V
    real(double),intent(in)::Vt,normF02
    real(double)::Fp(N)


    real(double)::scal

    integer::i
    ldir=.false.
    if (lover) then
       R=Rt
       V=Vt
       F=Ft
       return
    end if
    scal=0
    do i=1,N
       scal=scal+Ft(i)*F0(i)
    end do
    do i=1,N
       Fp(i)=scal*F0(i)/normF02
    end do
    call test_conv(N,Fp,ldir)
    if (ldir) then
       R=Rt
       V=Vt
       F=Ft
       return
    end if
  end subroutine checkline
end module steepestdescent_mod
