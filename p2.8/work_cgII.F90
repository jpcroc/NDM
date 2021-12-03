module work_cgII
   USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:  inv_angst, lperiod, rang,itmax,leev,sig, &
       it, itesauv, itesauvposition, itesauvforce,itmax, fnam,lenfnam,fnamcout,&
       inv_angst, erg2ev, angst,fpstop,fsumstop,itetabvois, &
       dmtype, potist,mdcg_noise,formatsauv,lspaceNDM,latcomp
  USE sauvegardeT_mod,only: sauvegardeT
  USE endrunT_mod,only: endrunT
  USE arret_ndm_mod,only: arret_ndm
#ifdef PARA
use Tpara,only:COMM_space,myidsp,nprocspace,para_space_config
#else
use Tpara,only:nprocspace,para_space_config
#endif
  USE atomconfig,only : atom_config
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
  use paraconfig,only:para_config
  USE parautils,only:initcomp,pointer_caltabt_calfo


  implicit none

  type(box_config)::boxcg
  type(atom_config)::atcgcomp
  type(cell_config)::cellcgcomp
  class(atom_config),pointer::atcgloc
  type(cell_config),pointer::cellcgloc
  type(cell_config),target:: cellcible ! ne sert qu'à faire pointer cellnebloc sur quelquechose
  type(atom_config),target::atcible
  type(para_config)::gcpara
  
contains

  subroutine FUNCT(N,X,F,G,NCALLS,psc)
    type(para_space_config)::psc
    real(double),intent(in):: X(N)
    real(double),intent(out)::G(N),F
    integer,intent(in) ::N,NCALLS
    integer ::i
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    !-----------------------------------------------
    !-----------------------------------------------
    integer::iproc,proc_source,cellx,celly,cellz
    real(double)::aux,auy,auz
    character :: extension*2
    integer::lenfn2,ko,i1
    real(double) :: fpmax,fpn,forctot,formax,fpmax_glob
    logical::lover

    logical:: lchg

    latcomp=.true. ! GC ==> latcomp=.true.
    lover=.false.

    it=NCALLS-1

    !    IF (3*ims.NE.N) THEN
    !       WRITE(0,'(a,i0)') "3*imm = ", 3*ims
    !       WRITE(0,'(a,i0)') "N     = ", N
    !       STOP "< work_cgII >"
    !    END IF
    if (gcpara%mpi_image%rank==0) then
       IF (dmtype.EQ.30) THEN ! Variables = reduced coordinates
          do i=1,atcgcomp%im
             i1=atcgcomp%num_at_glob(i)             
             atcgcomp%xp(1:3,i) = MatMul( boxcg%at, X(3*i1-2:3*i1) )
          end do
       ELSE ! Variables = cartesian coordinates (in A)
          do i=1,atcgcomp%im
             i1=atcgcomp%num_at_glob(i)             
             atcgcomp%xp(1:3,i) = X(3*i1-2:3*i1)*inv_angst
          end do
       END IF
    end if


    lchg=.true.
    call pointer_caltabt_calfo(sig,potist,atcgcomp,cellcgcomp,boxcg,atcgloc,cellcgloc,gcpara,lperiod,&
         &atcgcomp%ltabvois,it,itetabvois,lchg,psc,'xft') 

    if (it==1) then
       if (lEev.EQV..true.) then 
          if (rang==0) write(6,*)'Resultats en eV, Ang'
       else
          if (rang==0) write(6,*)'Resultats en cgs'
       end if
       if (rang==0)      write(*,'(70("="))')
       if (rang==0)      write(*,'("CG:     ","iter",10(" "),"epsi",14(" "),"Fmax",14(" "), "Energy")')
       if (rang==0)      write(*,'(70("="))')
    end if
    !    end if
#ifdef PARA
if (nprocspace.gt.1) then
       call COMM_space%barrier
    end if
#endif    
    IF (it.GE.1) THEN
       forctot=sqrt( SUM(atcgcomp%fp(1:3,1:atcgcomp%im)**2) )
       formax = MaxVal( Abs(atcgcomp%fp(:,1:atcgcomp%im)) )
       lover=.false.
       if (gcpara%mpi_image%rank==0)then

          if (lEev.EQV..true.) then 
             forctot = forctot*erg2eV/angst
             formax  = formax*erg2eV/angst
             if (rang==0) write(*,'("GC: ",i6,3E20.10)') it,forctot, formax, potist*erg2eV
             !           if (rang==0) write(*,'("GC: ",3E20.10)') forctot, formax, potist*erg2eV
             if (fpstop>0) then   
                if (formax.le.fpstop) then
                   if (rang==0) write(6,*)'force par atome  max  ev/Ang ', formax
                   if (rang==0) write (6, *) 'energie ', potist*erg2eV
                   !                 if (it.le.1) xp(:,:)=ax(:,:)
                   lover=.true.
                end if
             end if
             if (fsumstop>0) then   
                if (forctot.le.fsumstop) then
                   if (rang==0) write(6,*)'  sqrt ( sum_f F_i^2 ):   ev/Ang ', forctot
                   if (rang==0) write(6, *) 'energie ', potist*erg2eV
                   !                 if (it.le.1) xp(:,:)=ax(:,:)
                   lover=.true.
                   !call endrunT(atcgcomp,cellcgcomp,boxcg,latcomp)
                end if
             end if

          else
             if (rang==0) write(*,'("GC: ",i6,3E20.10)') it,forctot, formax, potist
             if (fpstop>0) then   
                if (formax.le.fpstop) then
                   if (rang==0) write(6,*)'force par atome  max cgs ',formax
                   if (rang==0) write (6, *) 'energie ', potist

                   lover=.true.
                   !call endrunT(atcgcomp,cellcgcomp,boxcg,latcomp)

                end if
             end if

             if (fsumstop>0) then   
                if (forctot.le.fsumstop) then
                   if (rang==0) write(6,*)'  sqrt ( sum_f F_i^2 ):  cgs  ', forctot
                   if (rang==0) write (6, *) 'energie ', potist
                   lover=.true.
                   !call endrunT(atcgcomp,cellcgcomp,boxcg,latcomp)
                end if
             end if

          end if
       end if
        
    end IF! it .ge.1
#ifdef PARA
if (nprocspace.gt.1) then
   call gcpara%mpi_image%BCAST(0,lover)
end if
#endif
!    write(6,*)'LOVER',lover,rang,it
       if (it>=itmax) then
          if (rang==0) write (6, *) '*******Derniere iteration **** '
          lover=.true.

       endif




       if (lover) then

          call endrunT(atcgcomp,cellcgcomp,boxcg,latcomp)

          call arret_ndm
       end if
       
       if (rang==0)then

       fnamcout = fnam(1:lenfnam)//'.cout'

       if (it.ne.0) then
          !       if (rang==0) then
          !           write(6,*)'work_cg_II analyse -> sauvegarde',it
          if (itesauv.GT.0) then
             if (mod(it,itesauv)==0) call sauvegardeT(atcgcomp,cellcgcomp,boxcg,formatsauv,fnamcout,latcomp)
          endif
       end if
       !go to into eV, ang and GC world............................................      
          F=potist*erg2eV

          G(:)=0.d0
          IF (dmtype.EQ.30) THEN ! Variables = reduced coordinates
             do i=1,atcgcomp%im
                i1=atcgcomp%num_at_glob(i)
                G(3*i1-2:3*i1)=-MatMul(atcgcomp%fp(:,i),boxcg%at)*erg2eV
             end do
          ELSE ! Variables = cartesian coordinates (in A)
             do i=1,atcgcomp%im
                i1=atcgcomp%num_at_glob(i)
                G(3*i1-2:3*i1)=-1*atcgcomp%fp(:,i)*erg2eV/angst
             end do
          END IF

       end if

    return
  end subroutine FUNCT



end module work_cgII
