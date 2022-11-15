! ***************************************************************
subroutine test_position_He (xp, at,ityp,rang,imm,im,it,ldesinteg,num_at_glob,nstepdes)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  ! USE 
  ! **************************************************************
  ! Programme par COUET Adrien, J-P Crocombette
  ! CEA-Saclay/DEN/DMN/SRMP
  ! Hiver 2008 version du 15 janvier 2009

  ! Test de la position des He dans une bulle de He par rapport à son centre
  ! ***************************************************************


  USE T_kind_param_m, ONLY:  double
#if(PARA)
  use mod_mpi
#endif
  implicit none
  !variables passes
  INTEGER ityp(imm),rang,imm,im,nstepdes
  real(double) :: xp(3,imm), H
  real(double) :: at(3,3) ! valeur des 3 vecteurs definissant la taille de laboite
  logical::ldesinteg
  integer, dimension(:), pointer       :: num_at_glob
  !variables locales

  real(double),save :: X, Y, Z, R, R_max,a,b,c,d,e,f,g,q,l,m,n,o,p, R_max_tot ! position du centre et distance du centre   au plus proche voisin
  integer :: kount,kount_tot,tt,i,it,ikt,nprt !compteur d'He sorti 
  integer,parameter:: nkountmax=1000
  integer,save ::icall=0,nkount(nkountmax),itdesdeb
  logical::ldesdes
  logical,save::lsort

  real(double),pointer,save:: natdt(:,:)
  integer::ntrrad,itrarad
  icall=icall+1
  tt=3
  if (nstepdes.gt.0) then
     nprt=nstepdes
  else
     nprt=1000
  end if
  if (icall==1)then
     open(unit=2121,file='bulle.in')
     read(2121,*) X,Y,Z,R
     R_max=R
     nkount(:)=0
     lsort=.false.
     itdesdeb=0
     ntrrad=Int(R/1e-9)+1
     allocate(ntrrad(ntrarad,ntyp)
     ntrrad(:,:)=0
  end if
  H=0
  !write(*,*) X, Y, Z, R
  kount=0
  ldesdes=.false.
  DO I=1,im

     !write(*,*) xp(1,I), xp(2,I), xp(3,I) 
     H=sqrt((xp(1,I)-X)**2+ &
          &	(xp(2,I)-Y)**2+&
          &	(xp(3,I)-Z)**2)


     itrrad=int(H/1e-9)+1
     ntrrad(itrrad,ityp(i))=ntrrad(itrrad,ityp(i))+1

     IF (ityp(I)==tt)then
        if ((ldesinteg==.true.).and.(i==1))then 	   
           IF (sqrt((xp(1,I)-X)**2+ &
                &	(xp(2,I)-Y)**2+&
                &	(xp(3,I)-Z)**2).gt.R)then
              ldesdes=.true. ;lsort=.true.

              if (itdesdeb==0)then
                 itdesdeb=it
                 write(6,*)'TPHE sortie',it,itdesdeb
              end if
           else
              !              ldesdes=.false.
              if (itdesdeb.ne.0)then
                 write(6,*);write(6,*)'TPHE entree; at 1 sorti de  a ',itdesdeb,it-1
                 itdesdeb=0
              end if

           end IF
        else
           IF (sqrt((xp(1,I)-X)**2+ &
                &	(xp(2,I)-Y)**2+&
                &	(xp(3,I)-Z)**2).gt.R)then
              !              write(6,*)'DESDES atome sorti',it
              ldesdes=.false.
              !              write(6,*)it,rang,I,'un He est sorti de la bulle!!!!!'

              IF (H.ge.R_max)then
                 R_max=H

              ENDIF

              kount=kount+1

           end if
	ENDIF

     ENDIF

  ENDDO
#if (PARA)
  call MPI_ALLREDUCE(kount,kount_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  kount=kount_tot
  call MPI_ALLREDUCE(R_max,R_max_tot,1,NDM_MPI_REAL_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)
  R_max=R_max_tot



#endif
  nkount(kount)=nkount(kount)+1
  if (mod(it,nprt)==0) then
     if (rang==0)then
        write(6,*)
        write(*,*)'TPHE Rayon maximum', R_max
        do ikt=0,nkountmax
           if(nkount(ikt).ne.0) write(6,*)'TPHE nbre de ',ikt,' : ',nkount(ikt),float(nkount(ikt))/it
        enddo

        write(914,*)
do
        write(914,*)

     endif
  end if
  !  if ((rang==0).and.(kount.ne.0)) then
  !     if(ldesdes==.false.) write(6,*)'it ; nb he sorti ',it,kount
  !  end if
  ! end if

  return
end subroutine test_position_He
		
