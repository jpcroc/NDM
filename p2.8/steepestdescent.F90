module steepestdescent_mod

  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:rang

  use WGC_mod,only:setV_F,nstep,ndir,beta,test_conv,ncalls,lvm,Rmin,unitgc

  !  real(double),allocatable,dimension (:,:)::X,R,G,H,F

  implicit none

contains

  subroutine steepestdescent(N,R,V,F,lover)
    integer::N
    real(double),dimension(:)::R,F
    real(double)::V
    logical ::lover

    integer::idir
    real(double)::forctot,formax
    real(double),dimension(N)::R0,F0,R1
    real(double)::V0,Vb,Vbs2
    logical ::lok,lvm
    
    call setV_F (N,R,V,F,lover,lvm)
        if (lvm)Rmin=R
    if (lover) then
       if (rang==0)write(unitgc,*)'NO NEED TO RELAX'
       if (rang==0)write(6,*)'NO NEED TO RELAX'
       return
    end if

    do idir=1,ndir
       if (rang==0)       write(unitgc,*)
       R0(1:N)=R(1:N)
       F0(1:N)=F(1:N)
       V0=V
       call mindir(lover,beta,N,R0,V0,F0,R,V,F,lOK)

       if (rang==0)       write(unitgc,*)'>>> minimization idirection; lOVER;  beta ',idir,lover,beta
       if (rang==0)       write(6,*)'>>> minimization idirection; lOVER ',idir,lover
       if (lover) then
!          write(6,*)
          if (rang==0)           write(unitgc,*)'*************************************'

          if (lok) then
             if (rang==0)           write(unitgc,*)
             if (rang==0)             write(unitgc,*)'RELAXED AFTER ',idir,' DIRECTIONS and ', NCALLS,' force calculations'
             if (rang==0)   write(6,*)'RELAXED AFTER ',idir,' DIRECTIONS and ', NCALLS,' force calculations'
             if (rang==0)   write(unitgc,*)'RELAXED AFTER ',idir,' DIRECTIONS and ', NCALLS,' force calculations'
             return
          else
             if (rang==0)  write(unitgc,*)'NOT RELAXED!!!!!!!'
             if (rang==0)  write(6,*)'NOT RELAXED!!!!!!!'
             call test_conv(N,F,lover,V,R,lvm)
             if (.not.lvm)R=Rmin
          end if
       end if

    end do
    return
  end subroutine steepestdescent

  subroutine conjugategradient(N,R,V,F,lover,lorig)
    integer::N
    real(double),dimension(:)::R,F
    real(double)::V
    logical ::lover
    logical,intent(in)::lorig

    integer::idir,i
    real(double)::forctot,formax,gamma
    real(double),dimension(N)::R0,F0,R1,G,H
    real(double)::V0,Vb,Vbs2,gigi,xixi
    logical ::lok

    call setV_F (N,R,V,F,lover,lvm)
    if (lvm)Rmin=R
!        write(unitgc,*)'post sVF0',V
    if (lover) then
       if (rang==0)       write(unitgc,*)'NO NEED TO RELAX'
       if (rang==0)       write(6,*)'NO NEED TO RELAX'
       return
    end if

    G=F
    H=F
    do idir=1,ndir
       R0(1:N)=R(1:N)
       V0=V
              write(unitgc,*)
              write(unitgc,*)'**************************'
              write(unitgc,*)'callmindir idir beta E0',idir,beta,V0
!              write(6,*)'callmindir idir  E0',idir,beta
       call mindir(lover,beta,N,R0,V0,H,R,V,F,lOK)

       if (rang==0) write(unitgc,'(A,I3,L2,E20.10)')' >>> minimization idirection; lOVER;  beta ',idir,lover,beta
       if (rang==0) write(6,'(A,I3,L2)')' >>> minimization idirection ',idir,lover
       if (lover) then
          if (lok) then
             if (rang==0)  write(unitgc,*)'RELAXED AFTER ',idir,' DIRECTIONS and ', NCALLS,' force calculations'
             if (rang==0)  write(6,*)'RELAXED AFTER ',idir,' DIRECTIONS and ', NCALLS,' force calculations'
             return
          else
             if (rang==0)  write(unitgc,*)'NOT RELAXED!!!!!!!'
             if (rang==0)  write(6,*)'NOT RELAXED!!!!!!!'
             call test_conv(N,F,lover,V,R,lvm)
             if (.not.lvm)R=Rmin
          end if
       end if

       xixi=0
       gigi=0
       if (lorig) then
          do i=1,N
             xixi=xixi+F(i)*F(i)
             gigi=gigi+G(i)*G(i)
          end do
       else
          do i=1,N
             xixi=xixi+(F(i)+G(i))*F(i)
             gigi=gigi+G(i)*G(i)
          end do
       end if
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
    real(double)::normF02,Vbeta,Vbs2,va,vb,vc,Vmin,Vbetatest,betatest,fhi

    integer::i,istep,id
    real(double)::a,b,c,AA,BB,betai,betaip1,Vbetai,Vbetaip1,Vd,d,ab
    logical ::ldir
    character*15::mic,mic2
    fhi=0.5*(1+sqrt(5.))
    lOK=.true.
    normF02=SUM(F0(:)**2)
    Rbeta(:)=R0(:)+beta*F0(:)

    call calcETcheck(N,Rbeta,Vbeta,Fbeta,lover,lvm,Rmin,R,V,F,R0,F0,V0,normF02,ldir)
    if (ldir) beta=beta/3
    if (lover.or.ldir) then
       lok=.true.
       return
    end if
    if  (VBeta.LT.V0) then
       write(unitgc,*)
       write(unitgc,*)'VBeta < V0',Vbeta,V0
       betai=beta;Vbetai=Vbeta
       loopG:          do i=1,nstep

          betaip1=betai*fhi
          write(unitgc,*)'betanew',betaip1
          Rbeta(:)=R0(:)+betaip1*F0(:)
          call calcETcheck(N,Rbeta,Vbeta,Fbeta,lover,lvm,Rmin,R,V,F,R0,F0,V0,normF02,ldir)
          write(unitgc,*)'betaold Vbeta',betaip1,Vbeta
          if (ldir) beta=beta/3
          if (lover.or.ldir) return
          Vbetaip1=Vbeta
          if (Vbeta.gt.Vbetai) then
             exit loopG
          end if
          betai=betaip1 ;Vbetai=Vbetaip1
       end do loopG
       b=betaip1 ; Vb=Vbetaip1
       c=betai; Vc=Vbetai


    else
       write(unitgc,*)
       write(unitgc,*)'VBeta > V0'
       betai=beta; Vbetai=Vbeta
       loopL:          do i=1,nstep
          betaip1=betai*(fhi-1)
          write(unitgc,*)'beta',betaip1
          Rbeta(:)=R0(:)+betaip1*F0(:)
          call calcETcheck(N,Rbeta,Vbeta,Fbeta,lover,lvm,Rmin,R,V,F,R0,F0,V0,normF02,ldir)
          Vbetaip1=Vbeta
          if (lover.or.ldir) return
          if (Vbeta.lt.Vbetai) then
             exit loopL
          end if
          betai=betaip1 ; Vbetai=Vbetaip1

       end do loopL
       b=betai ; Vb=Vbetai
       c=betaip1; Vc=Vbetaip1

    end if
    a=0 ; Va=V0

    write(unitgc,*)
    write(unitgc,'(A, 3E21.12)')"a0   b0   c0 ",a,b,c
    write(unitgc,'(A, 3E21.12)')"Va Vb Vc ",Va,Vb,Vc
    !    stop

    ab=(a+b)/2

    if (c.gt.ab) then
       d=a+(c-a)*(fhi-1)
       id=+1
    else
       d= b+(c-b)*(fhi-1)
       id=-1
    end if
    write(unitgc,'(A,E21.12)')'d=d',d

    do i=1,nstep
       beta=d
       !       write(unitgc,*)'mindir2 beta',beta
       Rbeta(:)=R0(:)+beta*F0(:)

       call calcETcheck(N,Rbeta,Vbeta,Fbeta,lover,lvm,Rmin,R,V,F,R0,F0,V0,normF02,ldir)
       !       write(unitgc,*)'vbeta', vbeta
       Vd=Vbeta
       if (lover.or.ldir) then
          lok=.true.
          return
       endif
       if (id==1) then
          if (Vd.lt.Vc) then
             b=c ;Vb=Vc
             c=d ; Vc=Vd
             d=a+(c-a)*(fhi-1)
             mic2='id=1 Vd<Vc'
             id=1
          else
             a=d ; Va=Vd
             id=-1
             d= b+(c-b)*(fhi-1)
             mic2='id=1 Vd>Vc'
          end if
       else
          if (Vd.lt.Vc) then
             a=c ; va=Vc
             c=d ; Vc=Vd
             d=b+(c-b)*(fhi-1)
             id=-1
             mic2='id=-1 Vd<Vc'
          else
             b=d ; Vb=vd
             d=a+(c-a)*(fhi-1)
             id=1
             mic2='id=-1 Vd>Vc'
          end if

       end if
       write(unitgc,*)
       write(unitgc,'(A,E21.12)')mic2,Vd
       write(unitgc,'(A, 4E21.12)')"a   b   c ",a,b,c,d
       write(unitgc,'(A, 3E21.12)')"Va Vb Vc ",Va,Vb,Vc
       write(unitgc,*)

    end do


    return
  end subroutine mindir


    subroutine calcETcheck (N,Rcalc,Vcalc,Fcalc,lover,lvm,Rmin,R,V,F,R0,F0,V0,normF02,ldir)
      logical,intent(out)::lover,lvm,ldir
      integer,intent(in)::N
      real(double)::Rcalc(N),Vcalc,Fcalc(N),Rmin(N)
      real(double),intent(in)::R0(N),V0,F0(N)
      real(double),intent(out)::R(N),V,F(N)
      real(double)::normF02
      lover=.false.
      ldir=.false.

      call setV_F (N,Rcalc,Vcalc,Fcalc,lover,lvm)
      if (lvm)Rmin=R
      call checkline(lover,ldir,F0,normF02,N,R,V,F,Rcalc,Vcalc,Fcalc)
      if (lover) then
         write(unitgc,*)'beta init relaxed'
         if (Vcalc.GT.V0) then
            V=V0;F=F0;R=R0
            write(unitgc,*)'STOP BETA '
         end if
         return
      end if
      if (ldir) then
         write(unitgc,*)'line search over: beta init'
         if (Vcalc.GT.V0) then
            V=V0;F=F0;R=R0
            lover=.true.
            write(unitgc,*)'STOP BETA init line'
         end if
         
         return
      end if
      return
    end subroutine calcETcheck

 
!!$  subroutine mindir0(lover,beta,N,R0,V0,F0,R,V,F,lOK)
!!$
!!$    logical,intent(out)::lover,lok
!!$    real(double)::beta
!!$    integer,intent(in)::N
!!$    real(double),intent(in)::R0(N),V0,F0(N)
!!$    real(double),intent(out)::R(N),V,F(N)
!!$
!!$    real(double),dimension(N)::Fp,RBs2,Rbeta,Fa,Fb,Fc,Fbeta,Fbs2,Fbetatest,Rbetatest
!!$    real(double)::normF02,Vbeta,Vbs2,va,vb,vc,Vmin,Vbetatest,betatest
!!$
!!$    integer::i,istep
!!$    real(double)::a,b,c,AA,BB
!!$    logical ::ldir
!!$    lOK=.true.
!!$    normF02=SUM(F0(:)**2)
!!$!    write(unitgc,*)
!!$!    write(unitgc,*)'normF02',normF02
!!$!        write(unitgc,*)
!!$    beta=beta*10
!!$    do i=1,nstep
!!$!       write(unitgc,*)
!!$      beta=beta/10
!!$       write(unitgc,*)'mindir1 betaS2',beta*0.5
!!$       Rbs2(:)=R0(:)+0.5*beta*F0(:)
!!$       call setV_F (N,Rbs2,Vbs2,Fbs2,lover,lvm)
!!$       if (lvm)Rmin=R
!!$!       write(unitgc,*)'mindir1 betaS2',lover,beta*0.5, Vbs2
!!$!       write(unitgc,*)'DIFF BS2',VBS2-V0
!!$       call checkline(lover,ldir,F0,normF02,N,R,V,F,Rbs2,Vbs2,Fbs2)
!!$       if (lover) then
!!$          write(unitgc,*)'betaS2 relaxed'
!!$          beta=beta/2
!!$          if (VBS2.GT.V0) then
!!$             V=V0;F=F0;R=R0
!!$             write(unitgc,*)'STOP BETAS2 '
!!$             lok=.false.
!!$          end if
!!$          return
!!$       end if
!!$       if (ldir) then
!!$          write(unitgc,*)'line search over: betas2'
!!$          beta=beta/3
!!$          if (VBS2.GT.V0) then
!!$             V=V0;F=F0;R=R0
!!$             lover=.true.
!!$             write(unitgc,*)'STOP BETAS2 line'
!!$             lok=.false.
!!$          end if
!!$
!!$          return
!!$       end if
!!$       if  (VBS2.LT.V0) exit
!!$    end do
!!$!    write(unitgc,*)
!!$    write(unitgc,*)'BS2 OK-> BETA'
!!$    Rbeta(:)=R0(:)+beta*F0(:)
!!$
!!$    do i=1,nstep
!!$!       write(unitgc,*)'mindir2 beta',beta
!!$       call setV_F (N,Rbeta,Vbeta,Fbeta,lover,lvm)
!!$           if (lvm)Rmin=R
!!$       write(unitgc,*)'mindir2 beta',lover, beta,Vbeta
!!$       call checkline(lover,ldir,F0,normF02,N,R,V,F,Rbeta,Vbeta,Fbeta)
!!$       if (lover) then
!!$          write(unitgc,*)'beta relaxed'
!!$          if ((VBS2.lt.Vbeta).or.(V0.lt.Vbeta)) then
!!$             if (V0.LT.VBS2) then
!!$                V=V0;R=R0;F=F0
!!$                lover=.true.
!!$                write(unitgc,*)'STOP BETA'
!!$                lok=.false.
!!$             else
!!$                V=Vbs2;R=Rbs2;F=fbs2
!!$             end if
!!$          endif
!!$          return
!!$       end if
!!$       if (ldir) then
!!$          write(unitgc,*)'line search over: beta'
!!$          if ((VBS2.lt.Vbeta).or.(V0.lt.Vbeta)) then
!!$             if (V0.LT.VBS2) then
!!$                V=V0;R=R0;F=F0
!!$                lover=.true.
!!$                write(unitgc,*)'STOP BETALINE'
!!$                lok=.false.
!!$             else
!!$                V=Vbs2;R=Rbs2;F=fbs2
!!$             end if
!!$          endif
!!$          return
!!$       end if
!!$!       write(unitgc,*)'DIFF BETA',VBeta-VBS2
!!$       if (Vbeta.gt.Vbs2) exit
!!$ !   write(unitgc,*)
!!$       beta=beta*2
!!$       Rbs2=Rbeta
!!$       Vbs2=Vbeta
!!$       Fbs2=Fbeta
!!$       Rbeta(:)=R0(:)+beta*F0(:)
!!$    end do
!!$ !   write(unitgc,*)
!!$   write(unitgc,*)'BETA OK->PARA'
!!$
!!$    a=0 ; b=beta ; c=beta/2
!!$    vb=vbeta ; Vmin=Vbs2
!!$    Va=V0
!!$    Vc=Vbs2
!!$   write(unitgc,*)'V',Va,Vb,Vc
!!$    do istep=1,nstep
!!$       AA=(vc-va)/((c-a)*(c-b))-(VB-VA)/((b-a)*(c-b))
!!$       BB=(VB-VA)/(b-a) -AA*(b+a)
!!$       betatest=-0.5*BB/AA     
!!$!       betatest=b-0.5*( ((b-a)**2)*(Vb-Vc)-(((b-c)**2)*(Vb-Va)))/ ((b-a)*(vb-vc)-(b-c)*(vb-va))
!!$       Rbetatest(:)=R0(:)+betatest*F0(:)
!!$      write(unitgc,*)'mindir3 beta IN',betatest
!!$       call setV_F(N,Rbetatest,Vbetatest,Fbetatest,lover,lvm)
!!$           if (lvm)Rmin=R
!!$      write(unitgc,*)'mindir3 beta OUT',betatest,Vbetatest
!!$       If (Vbetatest.gt.Vmin) then
!!$          write(unitgc,*)'PARABOLIC SERACH FAILURE SWITHING TO BINARY'
!!$          exit ! recherche parabolique en échec
!!$       end If
!!$       call checkline(lover,ldir,F0,normF02,N,R,V,F,Rbetatest,Vbetatest,Fbetatest)
!!$       if (lover) then
!!$          write(unitgc,*)'PARABOLIC relaxed'
!!$          beta=betatest
!!$          return
!!$       end if
!!$       if (ldir) then
!!$          write(unitgc,*)'PARABOLIC line search over'
!!$          beta=betatest
!!$          if (Vbetatest.GT.V0) then
!!$             V=V0;R=R0;F=F0
!!$             lover=.true.
!!$             write(unitgc,*)'STOP PARABOLIC'
!!$             lok=.false.
!!$          end if
!!$          return
!!$       end if
!!$       Vmin=Vbetatest
!!$       if (betatest.LT.c) then
!!$          b=c;Vb=Vc
!!$          c=betatest;Vc=Vbetatest
!!$       else
!!$          a=c; Va=Vc
!!$          c=betatest;Vc=Vbetatest
!!$       end if
!!$    end do
!!$    ! debut de la recherche binaire      
!!$    !      write(unitgc,*)
!!$!    return
!!$!      write(unitgc,*)'BINARY START'
!!$
!!$    a=0 ; b=beta ; c=beta/2
!!$    vb=vbeta ; Vmin=Vbs2
!!$    Va=V0
!!$    Vc=Vbs2
!!$    Vmin=Vbs2
!!$ !   write(unitgc,*)
!!$ !   write(unitgc,*)'BBIN',a,b,c
!!$!    write(unitgc,*)'VBIN',Va,Vb,Vc
!!$
!!$    do istep=1,nstep
!!$    if (Va.gt.Vb) then
!!$          a=c; Va=Vc
!!$       else
!!$          b=c; Vb=Vc
!!$       end if
!!$ !      write(unitgc,*)
!!$ !      write(unitgc,*)'BBIN',a,b
!!$ !      write(unitgc,*)'VBIN',Va,Vb
!!$       beta=(a+b)/2
!!$       c=beta
!!$       Rbeta(:)=R0(:)+beta*F0(:)
!!$       write(unitgc,'(A,E15.5)')'mindir4 beta',beta
!!$       call setV_F(N,Rbeta,Vbeta,Fbeta,lover,lvm)
!!$           if (lvm)Rmin=R
!!$      write(unitgc,*)'mindir4 beta',beta,Vbeta
!!$       call checkline(lover,ldir,F0,normF02,N,R,V,F,Rbeta,Vbeta,Fbeta)
!!$
!!$       Vc=Vbeta
!!$       if (lover) then
!!$          write(unitgc,*)'BINARY relaxed'
!!$          return
!!$       end if
!!$       if (ldir) then
!!$!          write(unitgc,*)'BINARY line search over'
!!$       
!!$          if (Vbeta.GT.V0) then
!!$             V=V0;R=R0;F=F0
!!$             lover=.true.
!!$!             write(unitgc,*)'STOP BINARY'
!!$             lok=.false.
!!$          end if
!!$
!!$          return
!!$       end if
!!$
!!$
!!$    end do
!!$    return
!!$  end subroutine mindir0


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
