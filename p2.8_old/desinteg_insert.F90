! ***********************************************************
!           sous-programme controle.f
! ***********************************************************

subroutine desinteg_insert
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
#if(PARA)
  use mod_mpi
#endif

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------

  real(double),pointer,save:: Wchemin(:)
  real(double),save :: Wch0
  real(double)::wch1,testval,u1,taup,taum,mug
  real(double), parameter::bkev=8.617385d-5 
  real(double),save::taupm=0,taumm=0

  integer:: i,k,l,m,n,ic,iaccept
  integer,save:: nchemin=0,nchacc,nchup=0,nchdn=0

  logical:: laccept

  allocate(wchemin(1+itmax/nstepdes))

  nchemin=nchemin+1
	if(rang==0) then
	  write(6,*)

     write(6,'(A,I9,I4,D21.12)')'DES findes, it nchemin deltaF',it,nchemin,deltaF*erg2eV
endif

!     deltaF

     Wchemin(nchemin)=deltaF*erg2eV
     if (nchemin==1) then
        Wch0=deltaF*erg2eV
        nchacc=1
#if(PARA)
        itichdn=ityp
        imdesdn=im
        num_at_globdesdn=num_at_glob
#endif
        xpchdn=xp
        vpchdn=vp
     else
        wch1=deltaF*erg2eV
     if(rang==0)   write(6,*)'DES ',wch1,wch0
        if (wch1.lt.wch0) then
           laccept=.true.
 if(rang==0)           write(6,*)'DES accept direct'
           iaccept=2
           testval=exp(-0.5*(wch1-wch0)/(erg2eV*Tempdes*bk))
           u1=0
        else
           laccept=.false. ; u1=0
           iaccept=0
           testval=exp(-0.5*(wch1-wch0)/(erg2eV*Tempdes*bk))
#if(PARA)
if (rang==0)then
#endif
           call random_number(u1)
           write(6,*)'DES testval',testval,u1
           if (u1.le.testval)then
              laccept=.true.
              iaccept=1
           end if
           write(6,*)'DES laccept ',laccept
#if(PARA)
endif
!	write(6,*)'PO',rang,laccept,u1
  call MPI_BCAST(laccept,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(u1,1,NDM_MPI_REAL_DOUBLE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
!	write(6,*)'P1',rang,laccept,u1
#endif

        end if
        if (rang==0)write(6,'(A,2I5,I3,2D17.8, I3,2D15.6)')'DES ',nchemin,nchacc, pm1des,wch1,wch0,iaccept,testval,u1

        if (pm1des==1) then
           taup=1./(exp(wch1/(2*bkev*Tempdes))+exp(wch0/(2*bkev*Tempdes)))
           taupm=(taupm*nchup+taup)/float(nchup+1)
           nchup=nchup+1
        else
           taum=1./(exp(wch1/(-2*bkev*Tempdes))+exp(wch0/(-2*bkev*Tempdes)))
           taumm=(taumm*nchdn+taum)/float(nchdn+1)
           nchdn=nchdn+1
        end if
        if ((nchup.ge.1).and.(nchdn.ge.1)) then
           mug=-log(taupm/taumm)*bkev*Tempdes
           if (rang==0)write(6,*)'DES desin_mu',mug,nchup,nchdn
        end if
        if (laccept.EQV..true.) then
           if (rang==0)write(6,*)'DES ACCEPT',nchemin
           nchacc=nchemin
           wch0=wch1
           if (pm1des==1) then   
              xpchup=xp
              vpchup=vp

              xpchdn=xpchdeb
              vpchdn=vpchdeb

#if(PARA)
              itichdn=itichdeb
              itichup=ityp

              imdesup=im
              imdesdn=imdesdeb

              num_at_globdesdn=num_at_globdesdeb
              num_at_globdesup=num_at_glob
#endif

		if (typspr==1) then
			xpspr0=xpspr

#if(PARA)	
			do i=1,im
				if (num_at_glob(i)==1) then
					xpspr=xp(:,i)
					call MPI_BCAST(xpspr,3,NDM_MPI_REAL_DOUBLE,rang,MPI_COMM_WORLD,ierr)			
				endif
			enddo
#else
			xpspr(:)=xp(:,1)
#endif
	endif
	   else			!pm1=-1
              xpchup=xpchdeb
              vpchup=vpchdeb


              xpchdn=xp
              vpchdn=vp


#if(PARA)
              itichup=itichdeb
              itichdn=ityp

              imdesup=imdesdeb
              imdesdn=im

              num_at_globdesup=num_at_globdesdeb
              num_at_globdesdn=num_at_glob
#endif



           end if
           

        else ! rejet
           if (rang==0) write(6,*)'DES REJET retour a ', nchacc
           if (pm1des==-1)then
              xp=xpchdn
              vp=vpchdn
#if(PARA)
              ityp=itichdn
              im=imdesdn
              num_at_glob=num_at_globdesdn
#endif


           else			!pm1=1
              xp=xpchup
              vp=vpchup

#if(PARA)
              ityp=itichup
              im=imdesup
              num_at_glob=num_at_globdesup
#endif

	if (typspr==1) then
		xpspr=xpspr0
	endif


           end if
        end if
     end if

     xpchdeb=xp
     vpchdeb=vp

#if(PARA)
     itichdeb=ityp
     imdesdeb=im
     num_at_globdesdeb=num_at_glob
     itichdeb=ityp
#endif




     pm1des=-pm1des
     itdes=0

     deltaF=0.0
     deltaEspr=0.0
     call caltabt


	  write(6,*)'DES -------------------------------------'

  return
end subroutine desinteg_insert
