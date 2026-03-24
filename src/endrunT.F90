module endrunT_mod
  USE arret_ndm_mod,only:arret_ndm
  use calcangle_mod,only:adft,adf0
  USE arret_ndm_mod,only:arret_ndm
  USE sauvegardeT_mod,only:sauvegardeT!,cin2gin
  USE calcdigr_mod,only: rdfT,rdf0
  USE rasmolT_mod,only:rasmolT
  USE gen_com_m, only:uwrt,lwrt,itesauv,lprtfat,lwgin,angst,unitP,cunitP,erg2eV,&
       &iteanapos,iteangle,iterasmol,itesigma,itetemp,linstantfda,&
       &linstantrdf,lpkbar,lprteat,lprteattotm,lprtsigat,unitP,iterdf,&
       &lwgin, lposmoy,l2T,angst,dmtype,iteration,lenfnam,rang,timel,&
       &fnamcout,fnam,lspaceNDM
  use var_pot, only: eatref,eatref,eatref
  USE cellconfig,only:cell_config
  USE atomconfig,only:atom_config,atom_config_d,atom_config_e
  use boxconfig,only: box_config
  use posana,only:anapos
  USE Tpara,only:nprocspace
  implicit none
contains
  ! ****************************************************************
  subroutine endrunT(atdml,celndm,boxndm,latcomp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE posana,only:
    USE elec_cell, ONLY:  sauveelec
    
    implicit none
    class(box_config)::boxndm
    class(atom_config)::atdml
    type(cell_config):: celndm
    logical,intent(in)::latcomp


    integer :: formatsauv
    if (lPkbar) then
       unitP=1.0d-9
       cunitP='kbar'
    else
       unitP=1.0
       cunitP='d/cm2'
    endif


 if (rang==0) then

    write (uwrt, *)
    write (uwrt, *)

    write (uwrt, *) '####### END OF RUN  ######## = ', iteration, '  time = ', timel
 endif

 IF (iteSauv.GE.0) then
    fnamcout= fnam(1:lenfnam)//'.cout'
    select type(atdml)
    type is (atom_config)
       formatsauv=2
    class is (atom_config_d)
       formatsauv=5
    class is (atom_config_e)
       if (atdml%lxpp)then
          formatsauv=3
       else
          formatsauv=5
       end if
    end select
          
    call sauvegardeT(atdml,celndm,boxndm,formatsauv,fnamcout,latcomp)     ! Modif E. Clouet: sauvegarde seulement si voulu
    if (l2T.and.rang==0) call sauveelec
 end IF

 if (.not.rdf0%linstantrdf) then
    if (iterdf>=0) call rdfT (rdf0)
 endif

 if (.not.linstantfda) then
    if (iteangle>=0) call adfT(adf0)

 endif
 select case (dmtype)
 case(2,21,22,23,24,3,30,31,32,33,34)
    iteration=0
 end select



 if (iterasmol.GE.0) call rasmolT (atdml,boxndm,999999999,latcomp=latcomp)
 if (lwgin) call rasmolT (atdml,boxndm,999999999,latcomp=latcomp,ivisumol=5)
! a retravailler if (nprocspace==1.and.iteanapos>=0) call anapos (atdml,celndm,boxndm,iteration)

 call arret_ndm


 stop
 return
end subroutine endrunT
end module endrunT_mod
