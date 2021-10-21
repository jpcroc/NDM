module tempinst_mod
  USE gen_com_m, ONLY:bk,lspacendm
  USE atomconfig,only: atom_config_d
#ifdef PARA
    USE mpi
    USE Tpara,only:COMM_space,nprocspace,myidsp
#endif

  implicit none
contains
  !c******************************************************************
  function tempinst(atcf,latcomp)     !calcul de la T instant.
    !c******************************************************************

    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m
    USE var_pot, ONLY:cm
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    class(atom_config_d),intent(in)::atcf
    real(double)::  tempinst
    logical,optional::latcomp
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    real(double) ::  mv2,v2
    integer :: i,imtot
    logical ::latc=.false.

    if(present(latcomp))latc=latcomp
    
    mv2=0.0
    do i = 1,atcf%im
       v2= atcf%vp(1,i)**2+ atcf%vp(2,i)**2+ atcf%vp(3,i)**2
       mv2= mv2 + cm(atcf%ityp(i))*v2
    enddo
    if (latc) then
       tempinst=mv2/(3.d0*float(atcf%im)*bk)
    else

#ifdef PARA
       if ((lspaceNDM).and.(nprocspace.gt.1))then 
          imtot=atcf%im
          call comm_space%sum(imtot)
          call comm_space%sum(mv2)
          tempinst=mv2/(3.d0*imtot*bk)
       else
          tempinst=mv2/(3.d0*atcf%im*bk)
       end if
#else
       tempinst=mv2/(3.d0*float(atcf%im)*bk)
#endif
    end if



    return
  end function tempinst

  ! thermostat Andersen pour l'atome ia

  subroutine andersenth(vpa,ma,Text, nu,tstep,bk)
    USE T_kind_param_m
    implicit none
    real(double) :: vpa(3), ma   ! vitesse et masse de l'atome cible
    real(double) :: Text,bk         ! température cible et CBoltz
    real(double) :: nu, tstep    ! fréquence de chocs et pas de temps

    ! variables locales

    real(double) :: sigm         ! largeur de la gausssienne
    real(double) :: rnd1,v1,v2,tnu,rt,vt
    integer :: ic

    tnu=tstep*nu
    call random_number(rnd1)
    if(rnd1.lt. tnu) then
       sigm=sqrt(Text*bk/ma)
       rt=2.0

       do ic=1,3
          do while (rt.gt.1.)
             call random_number(v1)
             call random_number(v2)            
             v1=2*v1-1.
             v2=2*v2-1.
             rt=v1*v1+v2*v2
          end do
          vt=v1*sqrt(-2.*log(rt)/rt)
          vpa(ic)=vt*sigm
       end do
    end if
    return
  end subroutine andersenth



end module tempinst_mod
