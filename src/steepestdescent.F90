module steepestdescent_mod
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, only:uwrt,lwrt,rang,fpstop,fsumstop,dmtype,itmax

  use WGC_mod,only:setV_F,nstep,ndir,beta,test_conv,ncalls,lvm,Rmin,unitgc,Fmin,fpstopsig,ityprel,&
       &lvarstop,fstpdecr,ft2,fm2,fs2,beta35,gammas,gammav

  !  real(double),allocatable,dimension (:,:)::X,R,G,H,F

  implicit none

contains

  subroutine Adamrel(N,R,V,F,lover,beta0)
    integer::N
    real(double)::beta0
    real(double),dimension(:)::R,F
    real(double)::V
    logical,intent(out) ::lover


    integer::idesc,j,k

    real(double),dimension(N)::vel
    real(double)::sk,Fsq,eps,beta35

    beta35=1d-10

    lover=.false.
    
    call setV_F (N,R,V,F,lover,lvm,idesc)
    if (lvm)then
       Rmin=R
       Fmin=F
    end if
    !        write(unitgc,*)'post sVF0',V
    if (lover) then
       if (rang==0)       write(unitgc,*)'NO NEED TO RELAX'
       if (rang==0)       write(uwrt,*)'NO NEED TO RELAX'
       return
    end if

    vel(1:N)=0
    sk=0
    Fsq=0
    do j=1,N
       Fsq=Fsq+F(j)*F(j)
    end do
    eps=sqrt(Fsq)/1d8
    do k=1,itmax
       vel(1:N)=gammav*vel(1:N)+(1-gammav)*F(1:N)*beta35
       Fsq=0
       do j=1,N
          Fsq=Fsq+F(j)*F(j)
       end do
       sk=gammas*sk+(1-gammas)*(1-gammas)*fsq
       vel(:)=vel(:)/(1-gammav**k)
       sk=sk/(1-gammas**k)
       write(uwrt,*)'sk',sk,1/(eps+sqrt(sk))
       R(:)=R(:)+vel(:)/(eps+sqrt(sk))
       call setV_F (N,R,V,F,lover,lvm,idesc)
       if (lvm)then
          Rmin=R
          Fmin=F
       end if
       !        write(unitgc,*)'post sVF0',V
       if (lover) then
          if (rang==0)       write(unitgc,*)'NO NEED TO RELAX'
          if (rang==0)       write(uwrt,*)'NO NEED TO RELAX'
          return
       end if

    end do
    return
        
  end subroutine Adamrel
  
  subroutine conjugategradient(N,R,V,F,lover,lorig,beta0)
    integer::N
    real(double)::beta0
    real(double),dimension(:)::R,F
    real(double)::V
    logical ::lover
    logical,intent(in)::lorig
    integer::idesc

    logical:: ldecr=.false.,lcritI
    real(double)::fsumstopI,fpstopI,fpstopsigI,fmts,fsts,fsigts
    
    call setV_F (N,R,V,F,lover,lvm,idesc)
    if (lvm)then
       Rmin=R
       Fmin=F
    end if
    !        write(unitgc,*)'post sVF0',V
    if (lover) then
       if (rang==0)       write(unitgc,*)'NO NEED TO RELAX'
       if (rang==0)       write(uwrt,*)'NO NEED TO RELAX'
       return
    end if

    if (lvarstop) then
       write(uwrt,*)'VARIABLE FPSTOP/SIGSTOP/SUMSTOP'
       fpstopI=fpstop
       fsumstopI=fsumstop
       fpstopsigI=fpstopsig
       lcritI=.false.
       if (ityprel==1) then
          write(uwrt,*)'fpstop fsumstop at start',fpstop,fsumstop
       else
          write(uwrt,*)'fsigstop at start',fpstopsig
       end if
       do while (.not.lcritI)
!!$          if (ityprel==1) then
!!$             write(uwrt,*)'fpstop fsumstop at this poistart',fpstop,fsumstop
!!$          else
!!$             write(uwrt,*)'fsigstop at start',fpstopsig
!!$          end if

          if (ityprel==1) then 
             if (fpstop.gt.0) then
                fmts=fm2/fstpdecr
                if (fmts.gt.fpstopI) then
                   ldecr=.true.
                   fpstop=fmts
                else
                   fpstop=fpstopI
                end if
             end if
             if (fsumstop.gt.0) then
                fsts=ft2/fstpdecr
                if (fsts.gt.fsumstopI) then
                   ldecr=.true.
                   fsumstop=fsts
                else
                   fsumstop=fsumstopI
                end if
             end if
             write(uwrt,*)'fpstop and/or fsumstop set to',fpstop,fsumstop
          else
             fsigts=fs2/fstpdecr
             if (fsigts.gt.fpstopsig) then
                ldecr=.true.
                fpstopsig=fsigts
             else
                fpstopsig=fpstopsigI
             end if
             write(uwrt,*)'fpstopsig set to',fpstopsig
          end if
          
          call do_CG(N,R,V,F,lover,lorig,beta0)
          ldecr=.false.
          write(uwrt,*)'--- end of intermediate minimization---'
          if ((fpstop==fpstopI).and.(fsumstop==fsumstopI).and.(fpstopsig==fpstopsigI))lcritI=.true.
       end do
    else
       call do_CG(N,R,V,F,lover,lorig,beta0)
    end if
    return
    
  end subroutine conjugategradient

  subroutine do_CG(N,R,V,F,lover,lorig,beta0)
    integer::N
    real(double)::beta0
    real(double),dimension(:)::R,F
    real(double)::V
    logical ::lover
    logical,intent(in)::lorig

    integer::idir,i,idesc
    real(double)::gamma
    real(double),dimension(N)::R0,F0,G,H
    real(double)::V0,gigi,xixi
    logical ::lok,ldirOK

    
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
       !              write(uwrt,*)'callmindir idir  E0',idir,beta
       call mindir(lover,beta,N,R0,V0,H,R,V,F,lOK,ldirOK,idesc)


       if (rang==0) write(unitgc,'(A,I3,2L2,E20.10)')' >>> minimization idirection; lOVER; LDIROK; beta ',idir,lover,ldirOK,beta
       if (rang==0) write(uwrt,'(A,I3,2L2,E20.10)')' >>> minimization idirection; lOVER; LDIROK; beta ',idir,lover,ldirOK,beta
       !       if (rang==0) write(uwrt,'(A,I3,L2)')' >>> minimization idirection lover ldirOK ',idir,lover,ldirOK
       if (lover) then
          if (lok) then
             if (rang==0)  write(unitgc,*)'RELAXED AFTER ',idir,' DIRECTIONS and ', NCALLS,' force calculations'
             if (rang==0)  write(uwrt,*)'RELAXED AFTER ',idir,' DIRECTIONS and ', NCALLS,' force calculations'
             return
          else
             if (rang==0)  write(unitgc,*)'NOT RELAXED!!!!!!!'
             if (rang==0)  write(uwrt,*)'NOT RELAXED!!!!!!!'
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
                write(uwrt,*) 'position set at minimum energy found during line search'
                R0=Rmin
                G=Fmin
                H=Fmin
             endif
          else
             fpstop=fpstop/1.5
             fsumstop=fsumstop/1.5
             write(unitgc,'(A,2G18.8)') 'fpstop, fsumstop divided by 1.5',fpstop,fsumstop
             write(uwrt,'(A,2G18.8)') 'fpstop, fsumstop divided by 1.5',fpstop,fsumstop
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
             select case (dmtype)
             case(32)
                H=F  ! steepestdescent H is put back to the simple force direction
             case (33,34)
             end select
          else
             write(unitgc,*) 'position set at minimum energy found during line search'
             write(uwrt,*) 'position set at minimum energy found during line search'
             R0=Rmin
             G=Fmin
             H=Fmin
          endif
       end if
          

    end do
    return

  end subroutine do_CG

  subroutine mindir(lover,beta,N,R0,V0,F0,R,V,F,lOK,ldir,idesc)
    logical,intent(out)::lover,lok
    real(double)::beta
    integer,intent(in)::N
    integer,intent(inout)::idesc
    real(double),intent(in)::R0(N),V0,F0(N)
    real(double),intent(out)::R(N),V,F(N)

    real(double),dimension(N)::Rbeta,Fbeta
    real(double)::normF02,Vbeta,va,vb,vc,fhi

    integer::i,id
    real(double)::a,b,c,betai,betaip1,Vbetai,Vbetaip1,Vd,d,ab
    logical ::ldir
    logical::linit
    character*15::mic2
    ldir=.false.
    fhi=0.5*(1+sqrt(5.))
    lOK=.true.
    linit=.false.
    normF02=SUM(F0(:)**2)
    Rbeta(:)=R0(:)+beta*F0(:)
!!$    write(uwrt,*)'IN mindir betaIN',beta
!!$    write(unitgc,*)'IN mindir betaIN',beta
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
!     write(uwrt,*)'DBG in calcetcheck'
      call setV_F (N,Rcalc,Vcalc,Fcalc,lover,lvm,idesc)
      if (lvm)Rmin=R
      call checkline(lover,ldir,F0,normF02,N,R,V,F,Rcalc,Vcalc,Fcalc)
      if (lover) then
         write(unitgc,*)'beta init relaxed'
!         if (Vcalc.GT.V0) then
!            V=V0;F=F0;R=R0
            write(unitgc,*)'STOP BETA '
!         end if
         return
      end if
      if (ldir) then
         write(unitgc,*)'line search over: beta init'
!         if (Vcalc.GT.V0) then
!            V=V0;F=F0;R=R0
!            lover=.true.
            write(unitgc,*)'STOP BETA init line'
!         end if
         
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
