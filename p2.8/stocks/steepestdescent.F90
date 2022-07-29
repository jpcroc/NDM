module steepestdescent_mod
   USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:rang,fpstop,fsumstop,sigstop

  use WGC_mod,only:setV_F,nstep,ndir,beta,test_conv,ncalls,lvm,Rmin,unitgc,Fmin,ityprel

  !  real(double),allocatable,dimension (:,:)::X,R,G,H,F

  implicit none

contains

  subroutine steepestdescent(N,R,V,F,lover,beta0)
    integer::N
    real(double),dimension(:)::R,F
    real(double)::V
    real(double)::beta0
    logical ::lover,ldirOK

    integer::idir,idesc
    real(double)::forctot,formax
    real(double),dimension(N)::R0,F0,R1
    real(double)::V0,Vb,Vbs2
    logical ::lok,lvm
    idesc=0
    call setV_F (N,R,V,F,lover,lvm,idesc)
    if (lvm)then
       Rmin=R
       Fmin=F
    end if
    if (lover) then
       if (rang==0)write(unitgc,*)'NO NEED TO RELAX'
       if (rang==0)write(6,*)'NO NEED TO RELAX'
       return
    end if
    R0(1:N)=R(1:N)
    F0(1:N)=F(1:N)
    V0=V

    do idir=1,ndir
       if (rang==0)       write(unitgc,*)

       R0(1:N)=R(1:N)
       F0(1:N)=F(1:N)
       V0=V
       idesc=0

       call mindir(lover,beta,N,R0,V0,F0,R,V,F,lOK,ldirOK,idesc)

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
          end if
       end if
       if (idesc==0) then
          write(6,*) 'no decrease of energy along this line '
          write(unitgc,*) 'no decrease of energy along this line '
          if (.not.lvm) then
             write(6,*) 'no force convergence along line ::RETRY'
             write(unitgc,*) 'no force convergence along line ::RETRY'
             fpstop=fpstop/3
             fsumstop=fsumstop/3
             write(unitgc,'(A,2G18.8)') 'fpstop, fsumstop divided by 3',fpstop,fsumstop
             write(6,'(A,2G18.8)') 'fpstop, fsumstop divided by 3',fpstop,fsumstop
!             call arret_ndm
          else
             R0(1:N)=R(1:N)
             F0(1:N)=F(1:N)
             V0=V
          end if
       else
          write(unitgc,*) 'position set at minimum energy found during line search'
          write(6,*) 'position set at minimum energy found during line search'
          R=Rmin
          F=Fmin
       end if
       if (ldirOK) then
          beta0=beta
       else
          beta=beta0
          write (unitgc,*)'RESET STEEP beta',beta
!          G=F
!          H=F
       endif

    end do
    return
  end subroutine steepestdescent

  subroutine conjugategradient(N,R,V,F,lover,lorig,beta0)
    integer::N
    real(double)::beta0
    real(double),dimension(:)::R,F
    real(double)::V
    logical ::lover
    logical,intent(in)::lorig

    integer::idir,i,idesc
    real(double)::forctot,formax,gamma
    real(double),dimension(N)::R0,F0,R1,G,H
    real(double)::V0,Vb,Vbs2,gigi,xixi
    logical ::lok,ldirOK
    idesc=0
    call setV_F (N,R,V,F,lover,lvm,idesc)
    idir=-1
    if (lvm)then
       Rmin=R
       Fmin=F
    end if
    !        write(unitgc,*)'post sVF0',V
    if (lover) then
       if (rang==0)       write(unitgc,*)'NO NEED TO RELAX'
       if (rang==0)       write(6,*)'NO NEED TO RELAX'
       return
    end if
    R0=R
    F0=F
    G=F
    H=F
    do idir=1,ndir
       !       R0(1:N)=R(1:N)
       V0=V
       idesc=0
       write(unitgc,*)
       write(unitgc,*)'**************************'
       write(unitgc,*)'callmindir idir beta E0',idir,beta,V0
       !              write(6,*)'callmindir idir  E0',idir,beta
       call mindir(lover,beta,N,R0,V0,H,R,V,F,lOK,ldirOK,idesc)


       if (rang==0) write(unitgc,'(A,I3,2L2,E20.10)')' >>> minimization idirection; lOVER; LDIROK; beta ',idir,lover,ldirOK,beta
       if (rang==0) write(6,'(A,I3,2L2,E20.10)')' >>> minimization idirection; lOVER; LDIROK; beta ',idir,lover,ldirOK,beta
       !       if (rang==0) write(6,'(A,I3,L2)')' >>> minimization idirection lover ldirOK ',idir,lover,ldirOK
       if (lover) then
          if (lok) then
             if (rang==0)  write(unitgc,*)'RELAXED AFTER ',idir,' DIRECTIONS and ', NCALLS,' force calculations'
             if (rang==0)  write(6,*)'RELAXED AFTER ',idir,' DIRECTIONS and ', NCALLS,' force calculations'
             return
          else
             if (rang==0)  write(unitgc,*)'NOT RELAXED!!!!!!!'
             if (rang==0)  write(6,*)'NOT RELAXED!!!!!!!'
             call test_conv(N,F,lover,V,R,lvm)
             if (.not.lvm)then
                R=Rmin
                F=Fmin
             end if
          end if
       end if
       if (ityprel==1) then 
          if (idesc.gt.0) then 
             if (ldirOK) then
                beta0=beta
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
                R0(1:N)=R(1:N)
             else
                write(unitgc,*) 'position set at minimum energy found during line search'
                write(6,*) 'position set at minimum energy found during line search'
                R0=Rmin
                G=Fmin
                H=Fmin
             endif
          else
             fpstop=fpstop/3
             fsumstop=fsumstop/3
             write(unitgc,'(A,2G18.8)') 'fpstop, fsumstop divided by 3',fpstop,fsumstop
             write(6,'(A,2G18.8)') 'fpstop, fsumstop divided by 3',fpstop,fsumstop
          end if
       else
          if (ldirOK) then
             beta0=beta
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
             R0(1:N)=R(1:N)
          else
             write(unitgc,*) 'position set at minimum energy found during line search'
             write(6,*) 'position set at minimum energy found during line search'
             R0=Rmin
             G=Fmin
             H=Fmin
          endif
       end if
          

    end do
    return

  end subroutine conjugategradient

  subroutine mindir(lover,beta,N,R0,V0,F0,R,V,F,lOK,ldir,idesc)
    logical,intent(out)::lover,lok
    real(double)::beta
    integer,intent(in)::N
    integer,intent(inout)::idesc
    real(double),intent(in)::R0(N),V0,F0(N)
    real(double),intent(out)::R(N),V,F(N)

    real(double),dimension(N)::Fp,RBs2,Rbeta,Fa,Fb,Fc,Fbeta,Fbs2,Fbetatest,Rbetatest
    real(double)::normF02,Vbeta,Vbs2,va,vb,vc,Vmin,Vbetatest,betatest,fhi

    integer::i,istep,id
    real(double)::a,b,c,AA,BB,betai,betaip1,Vbetai,Vbetaip1,Vd,d,ab
    logical ::ldir
    logical::linit
    character*15::mic,mic2
    ldir=.false.
    fhi=0.5*(1+sqrt(5.))
    lOK=.true.
    linit=.false.
    normF02=SUM(F0(:)**2)
    Rbeta(:)=R0(:)+beta*F0(:)
    write(6,*)'IN mindir betaIN',beta
    write(unitgc,*)'IN mindir betaIN',beta
    call calcETcheck(N,Rbeta,Vbeta,Fbeta,lover,lvm,Rmin,R,V,F,R0,F0,V0,normF02,ldir,idesc)
!    if (ldir) beta=beta/3
    if (lover.or.ldir) then
       write(unitgc,*)'retout direct de mindir',lover,ldir
       lok=.true.
       return
    end if
    if  (VBeta.LT.V0) then
       write(unitgc,*)
       write(unitgc,*)'VBeta < V0',Vbeta,V0
       betai=beta;Vbetai=Vbeta
       loopG:          do i=1,nstep

          betaip1=betai*fhi
          write(unitgc,*)'betanew step',betaip1,i,nstep
          Rbeta(:)=R0(:)+betaip1*F0(:)
          call calcETcheck(N,Rbeta,Vbeta,Fbeta,lover,lvm,Rmin,R,V,F,R0,F0,V0,normF02,ldir,idesc)
          write(unitgc,*)'betaold Vbeta',betaip1,Vbeta
          if (ldir) beta=betaip1/3
          if (lover.or.ldir) then
             write(unitgc,*)'retout de mindir dans etape1.1',lover,ldir
             lok=.true.
             return
          end if
          !if (lover.or.ldir) return
          Vbetaip1=Vbeta
          if (Vbeta.gt.Vbetai) then
             linit=.true.
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
!          write(unitgc,*)'beta',betaip1
          write(unitgc,*)'betanew step',betaip1,i,nstep
          Rbeta(:)=R0(:)+betaip1*F0(:)
          call calcETcheck(N,Rbeta,Vbeta,Fbeta,lover,lvm,Rmin,R,V,F,R0,F0,V0,normF02,ldir,idesc)
          Vbetaip1=Vbeta
          if (lover.or.ldir) then
             write(unitgc,*)'retout de mindir dans etape1.2',lover,ldir
             lok=.true.
             return
          end if
!          if (lover.or.ldir) return
          if (Vbeta.lt.Vbetai) then
             linit=.true.
             exit loopL
          end if
          betai=betaip1 ; Vbetai=Vbetaip1

       end do loopL
       b=betai ; Vb=Vbetai
       c=betaip1; Vc=Vbetaip1

    end if

    if (.not.linit)    then
       write(unitgc,*)'INIT mindir failed'
       ldir=.false.
       return
    end if
    a=0 ; Va=V0

    write(unitgc,*)
    write(unitgc,'(A, 3E21.12)')"a0   c0   b0 ",a,c,b
    write(unitgc,'(A, 3E21.12)')"Va Vc Vb ",Va,Vc,Vb
    !    call arret_ndm

    ab=(a+b)/2

    if (c.gt.ab) then
       d=a+(c-a)*(fhi-1)
       id=+1
    else
       d= b+(c-b)*(fhi-1)
       id=-1
    end if
    write(unitgc,'(A,E21.12)')'d=',d
    write(unitgc,'(A,I4)')'id=',id

    do i=1,nstep
       beta=d
       !       write(unitgc,*)'mindir2 beta',beta
       Rbeta(:)=R0(:)+beta*F0(:)
       write(unitgc,*)'betanew step',beta,i,nstep
       call calcETcheck(N,Rbeta,Vbeta,Fbeta,lover,lvm,Rmin,R,V,F,R0,F0,V0,normF02,ldir,idesc)
       !       write(unitgc,*)'vbeta', vbeta
       Vd=Vbeta
       if (lover.or.ldir) then
          write(unitgc,*)'retout de mindir dans etape2',lover,ldir
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
       write(unitgc,'(A, 4E21.12)')"a   c   d ",a,c,b,d
       write(unitgc,'(A, 3E21.12)')"Va Vc Vb ",Va,Vc,Vb
       write(unitgc,*)

    end do
    ldir=.false.
    write(unitgc,*)'ECHEC de mindir betaOUT',beta

    return
  end subroutine mindir


    subroutine calcETcheck (N,Rcalc,Vcalc,Fcalc,lover,lvm,Rmin,R,V,F,R0,F0,V0,normF02,ldir,idesc)
      logical,intent(out)::lover,lvm,ldir
      integer,intent(in)::N
      integer:: idesc
      real(double)::Rcalc(N),Vcalc,Fcalc(N),Rmin(N)
      real(double),intent(in)::R0(N),V0,F0(N)
      real(double),intent(out)::R(N),V,F(N)
      real(double)::normF02
      lover=.false.
      ldir=.false.
!     write(6,*)'DBG in calcetcheck'
      call setV_F (N,Rcalc,Vcalc,Fcalc,lover,lvm,idesc)
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
!            lover=.true.
            write(unitgc,*)'STOP BETA init line'
         end if
         
         return
      end if
      return
    end subroutine calcETcheck

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
    write(unitGC,*)'SCAL',scal
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
