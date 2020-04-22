module test_pos_He_mod
        implicit none
        contains
! ***************************************************************
subroutine test_position_He (xp, at,ityp,rang,imm,im,it,ldesinteg,num_at_glob,nstepdes,itmax)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  ! USE 
  ! **************************************************************
  ! Programme par COUET Adrien, J-P Crocombette
  ! CEA-Saclay/DEN/DMN/SRMP
  ! Hiver 2008 version du 15 janvier 2009

  ! Test de la position des He dans une bulle de He par rapport Ã  son centre
  ! ***************************************************************


  USE T_kind_param_m, ONLY:  double
#ifdef PARA
  USE mod_para
#endif
  implicit none
  !variables passées
  INTEGER ityp(imm),rang,imm,im,nstepdes,itmax
  real(double) :: xp(3,imm), H
  real(double) :: at(3,3) ! valeur des 3 vecteurs definissant la taille de laboite
  logical::ldesinteg
  integer, dimension(:), allocatable       :: num_at_glob
  !variables locales

  real(double),save :: X, Y, Z, R, Rmax,a,b,c,d,e,f,g,q,l,m,n,o,p, Rmax_tot,rmax1,rmin ! position du centre et distance du centre   au plus proche voisin
  integer :: kount,kount_tot,tt,i,it,ikt!compteur d'He sorti 
  integer,parameter:: nkountmax=1000
  integer,save ::icall=0,nkount(nkountmax),itdesdeb,nprt
  logical::ldesdes
  logical,save::lsort
  integer,save::ntrm
  real(double),save,allocatable ::at_par_tr(:),natsortideplus(:),natsortidemoins(:)
  integer::ntr
  

  icall=icall+1
  tt=0
!  if (nstepdes.gt.0) then
!     nprt=nstepdes
!  else
!   nprt=1000
!end if
  if (icall==1)then
     open(unit=2121,file='bulle.in')
     read(2121,*) X,Y,Z,Rmin,nprt,ntrm,tt

     write(6,*) 'bulle.in', X,Y,Z,Rmin,nprt,ntrm
!     stop
     allocate(at_par_tr(ntrm))
     allocate(natsortideplus(ntrm))
     allocate(natsortidemoins(ntrm))
     Rmax=Rmin
	rmax1=rmin
     nkount(:)=0
     lsort=.false.
     itdesdeb=0
     at_par_tr(:)=0

  end if
  if (mod(it,nprt)==0) then
     open(unit=247,file='bulle.out')
     open(unit=248,file='b2.out')
endif
  H=0
  !write(*,*) X, Y, Z, R
  kount=0
  ldesdes=.false.
  DO I=1,im
     IF (ityp(I)==tt)then
        !write(*,*) xp(1,I), xp(2,I), xp(3,I) 
        H=sqrt((xp(1,I)-X)**2+ &
             &	(xp(2,I)-Y)**2+&
             &	(xp(3,I)-Z)**2)
        if ((ldesinteg.EQV..true.).and.(i==1))then

           IF (H.gt.Rmin)then
!              ldesdes=.true. ;lsort=.true.

!              if (itdesdeb==0)then
!                 itdesdeb=it
!                 write(6,*)'TPHE sortie',it,itdesdeb
!              end if
!           else
!              ldesdes=.false.
!              if (itdesdeb.ne.0)then
!                 write(6,*);write(6,*)'TPHE entree; at 1 sorti de  a ',itdesdeb,it-1
!                 itdesdeb=0
!              end if
              if(H.gt.rmax1) rmax1=h
              
           end IF
        endif
        IF (H.gt.Rmin)then
  if (mod(it,nprt)==0) then
              write(248,*)i,H*1d8
!              write(6,*)i,H*1d8
	endif
           !              write(247,*)it,i, H
           IF (H.ge.Rmax)  Rmax=H
           !		write(6,*) H-rmin, (H-rmin)*1.0d9
           ntr=min(1+int((H-rmin)*1.0d9),ntrm)
!           	write(6,*)'ntr ',ntr,nprt
           at_par_tr(ntr)=at_par_tr(ntr)+1
           !              kount=kount+1
           
        end if
        !	ENDIF
        
     ENDIF

  ENDDO
#ifdef PARA
  call MPI_ALLREDUCE(kount,kount_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  kount=kount_tot
  call MPI_ALLREDUCE(Rmax,Rmax_tot,1,NDM_MPI_REAL_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)
  Rmax=Rmax_tot



#endif
!  nkount(kount)=nkount(kount)+1
  if (mod(it,nprt)==0) then
     natsortideplus(:)=0.
     natsortidemoins(:)=0.
     if (rang==0)then
        write(6,*)
        write(6,*)'TPHE Rayon maximum', Rmax
        write(6,*)'TPHE Rayon maximum pour 1', Rmax1
        write(6,*)
        do ikt=1,ntrm
           if(at_par_tr(ikt).Gt.0.) then
              natsortideplus(1:ikt)=natsortideplus(1:ikt)+at_par_tr(ikt)/icall
              natsortidemoins(ikt:ntrm)=natsortidemoins(ikt:ntrm)+at_par_tr(ikt)/icall
           end if
        enddo
        do ikt=1,ntrm
           if(at_par_tr(ikt).Gt.0.) then
!              write(6,*)'TPHE tr',ikt,' avec ', at_par_tr(ikt)/icall,natsortideplus(ikt)
              write(247,'(I3,3F12.5)')ikt, at_par_tr(ikt)/icall,natsortideplus(ikt),natsortidemoins(ikt)
           endif
        enddo
     endif
 close (247)
 close (248)
    
   end if
!  if ((rang==0).and.(kount.ne.0)) then
!     if(ldesdes==.false.) write(6,*)'it ; nb he sorti ',it,kount
!  end if
  ! end if
  return
end subroutine test_position_He
		
end module
