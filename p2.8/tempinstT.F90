module tempinstT_mod
  USE gen_com_m, ONLY:bk,im_glob,lspaceNDM
#ifdef PARA
    USE Tpara,only:COMM_space,nprocspace

#endif
    use atomconfig,only:atom_config_d
  implicit none
contains
  !c******************************************************************
  function tempinstT(atcf,kine,latcomp) 
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
    real(double)::  tempinstT
    real(double),optional::kine
    logical,optional:: latcomp
    logical latc
    real(double)::kinetot
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    real(double) ::  mv2,v2
    integer :: i

    latc=.false.
    if (present(latcomp))latc=latcomp

    mv2=0.0
    do i = 1,atcf%im
       v2= atcf%vp(1,i)**2+ atcf%vp(2,i)**2+ atcf%vp(3,i)**2
       mv2= mv2 + cm(atcf%ityp(i))*v2
       if (present(kine))kine=kine+v2*0.5*cm(atcf%ityp(i))
    enddo

#ifdef PARA
    if ((lspaceNDM).and.(nprocspace.gt.1).and.(latc.eqv..false.))then 
    call comm_space%sum(mv2)
!    mv2 = mv2_glob
    tempinstT=mv2/(3.d0*float(im_glob)*bk)
       if(present(kine)) then
          call comm_space%sum(kine)
       end if
    else
       tempinstT=mv2/(3.d0*float(atcf%im)*bk)
    end if
#else
    tempinstT=mv2/(3.d0*float(atcf%im)*bk)


#endif




    return
  end function tempinstT

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



end module tempinstT_mod
