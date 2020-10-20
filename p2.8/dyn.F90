module dyn_mod
  USE gen_com_m, ONLY:cunite,erg2ev,fnemd,it,itetconst,lcalcjq,leev,lnemd,lperiod,&
       &ltcon,text,timel,tstep,unite,usdh,eatom,bk
  USE tempinst_mod,only: tempinst
!  USE period_mod,only: period
  implicit none
contains
  ! *************************************************************
  subroutine dyn (xp,xpp,vp,fp,ityp,im)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE var_pot, ONLY:ntyp,cm
    USE jqmod
!    USE tab_imm_m,only:ax,xp,vp,xpp,fp,ityp
    USE tempinst_mod,only: tempinst

!#ifdef PARA
!    USE mod_para,only:MPI_COMM_space
!#endif

    implicit none

    real(double),allocatable::xp(:,:),vp(:,:),xpp(:,:),fp(:,:)
    integer, allocatable::ityp(:)
    integer::im
    
    integer :: i, iti, ic
    real(double), dimension(ntyp) :: aux
    real(double) :: xprov ,vv
    real(double), save :: tmoyinst, imesureT
    !real(double), external :: tempinst
    real(double) :: tempavant,tmoy
    real(double) :: deltaE
    real(double):: eatommoy
    if(lEev) then
       unitE=erg2eV
       cunitE='  eV '
    else
       unitE=1.0
       cunitE=' erg '
    end if

    if (allocated(eatom))  eatom(1:im)=eatom(1:im)+0.5*cm(ityp(1:im))*(vp(1,1:im)**2+vp(2,1:im)**2+vp(3,1:im)**2)

    if (lnemd) then
       eatommoy=0.
       do i=1,im
          !        eatom(i)=eatom(i)+0.5*cm(ityp(i))*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
          eatommoy=eatommoy+eatom(i)/float(im)
       end do
       do i=1,im

          fp(1,i)=fp(1,i)+(eatom(i)-eatommoy)*Fnemd
       end do
       !     do i=1,im
       !        eatom(i)=eatom(i)-0.5*cm(ityp(i))*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
       !     end do
    end if



    timel = timel+tstep
    aux(:ntyp) = tstep**2/cm(:ntyp)
    !      write(6,*)'aux ',aux


    !debug write(*,*) 'md_test1',  xp(1,1), xpp(1,1), vp(1,1) 
    do i = 1, im
       do ic = 1, 3
          xprov = (xp(ic,i)-xpp(ic,i))+xp(ic,i)+aux(ityp(i))*fp(ic,i)
          vp(ic,i) = (xprov-xpp(ic,i))*usdh
          xpp(ic,i) = xp(ic,i)
          xp(ic,i) = xprov
       end do
    end do
    !debug write(*,*) 'md_test2',  xp(1,1), xpp(1,1), vp(1,1)

    if (lcalcjq) then
       eatommoy=0.
       jqp=jq ; jqk=0.0 !; expvect(:)=0.0
       do i=1,im
          !        eatom(i)=eatom(i)+0.5*cm(ityp(i))*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
          if(lnemd) eatommoy=eatommoy+eatom(i)/float(im)
          expvect(:)=expvect(:)+eatom(i)*xpp(:,i)
          do ic=1,3              
             jqk(ic)=jqk(ic)+eatom(i)*vp(ic,i)
             !               jq(ic)=jq(ic)+eatom(i)*vp(ic,i)
          end do
       end do
       jq=jqp+jqk



       write(65,'(I8,3D16.8)')it-1,jq(1),jq(2),jq(3)
       !         write(66,'(I8,3D15.6)')it-1,expvect(1),expvect(2),expvect(3)
    end if


    ! temp const moyenne depuis le dernier rescale
    if(lTcon) then
       tempavant=tempinst(vp,ityp,im, size(ityp))
       tmoyinst = tmoyinst +tempavant 
       imesureT = imesureT + 1
       if(mod(it,iteTconst).eq.0) then
          write(6,*)
          tmoy=tmoyinst/imesureT 
          write(6,'(A,3F12.2)')'> Temperature Constante = Text, Tmoy, Tinst ', Text,tmoy,tempavant


          deltaE=(tmoy-Text)*3*bk*unitE
          write(6,'(A,I0,D21.12,A)')'IT   modification d_energie par atome ',it,deltaE,&
               &cunitE
          !       write(6,'(A,I,D21.12,A)')'IT , modification d_energie par atome ',IT,deltaE, cunitE
          write(6,*)

          vv = sqrt(Text/(tmoy))
          xpp(:,:im) = xp(:,:im)-(xp(:,:im)-xpp(:,:im))*vv
          vp(:,:im) = vp(:,:im)*vv
          tmoyinst=0.
          imesureT=0
       endif
    endif

!    if (lperiod)       call period  (im,xp,xpp)





    return
  end subroutine dyn

end module dyn_mod
