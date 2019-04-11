
! *************************************************************
subroutine calctemp(temptyp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use var_pot
  use gen_com_m
  use tab_imm_m
  use elec_cell, only: ecell,i2T,nex,ney,nez,nox_2_nex
  use eloss, only : tcelec,ecelec
#if(PARA)
  use mod_mpi
#endif

  ! *************************************************************

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(double) , intent(inout) :: temptyp(ntyp)
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: ic, i, iti, ko, i2,kx,ky,kz,koo,ixe,iye,ize
  real(double), dimension(ntyp) :: v2  ,tempmaxat(ntyp)
  real(double) :: vpn2,tat,ekin
  real(double), dimension(ntyp,3) :: vx2
  integer::ixyze(3),nats
  ! ym      real(double), dimension(ntyp,nce) :: v2c
  !  real(double), dimension (:),allocatable ::tempc,tempcm
#if(PARA)
  real(double), dimension(ntyp) :: v2_glob
  real(double), dimension(ntyp,3) :: vx2_glob
  real(double), dimension(noxyz)::tempc_tot
  real(double),allocatable::tempiontot(:,:,:)
  integer,allocatable::niontot(:,:,:)
  integer::natstot
  real(double):: tempEPtot
#endif

  !-----------------------------------------------
  !
  !
  ! local variables
  if ((ltpcel).or.(tempstopcel.gt.0).or.(tcelec.gt.0)) then
     !     allocate (tempc(noxyz))
     !     allocate (tempcm(noxyz))
     tempc(:)=0.
     tempcm(:)=0.
  endif
  temp = 0.0
  kine = 0.0
  temptyp(:ntyp) = 0.0
  v2(:ntyp) = 0.0
  vx2(:ntyp,:) = 0.0

  if(ltpcel.or.(tcelec.gt.0))tempmaxat(:)=0

  if (l2t)then
     ecell(:,:,:)%tempIon=0
     ecell(:,:,:)%nIon=0
     ecell(:,:,:)%nIonS=0
     nats=0
#if(PARA)
     allocate(tempiontot(nex,ney,nez))
     allocate(niontot(nex,ney,nez))
#endif
  end if

  tempEP=0

  do ko = 1, noxyz

     if (nato(ko)==0) cycle

#if(PARA)
     if (proc_cell(ko).ne.myid) cycle

#endif

     if (L2T)     call nox_2_nex(ko,ixyze)

     do i2 = 1, nato(ko)

        i = last(i2,ko)
        if (num_at_glob(i).gt.im_glob) cycle
        vpn2 = vp(1,i)**2+vp(2,i)**2+vp(3,i)**2
        !calculation of ionic temperature and number of ions in the electronic cell (only slow moving ions)  
        if (l2T.eqv..true.) then
           select case(i2t)
           case(1)
              ecell(ixyze(1),ixyze(2),ixyze(3))%nIon= ecell(ixyze(1),ixyze(2),ixyze(3))%nIon+1
              ekin=0.5*erg2ev*vpn2*cm(ityp(i))
              if (ekin.lt.Ecelec) then
                 ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon=ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon&
                      &+vpn2*cm(ityp(i))/(3.0*bk)
                 ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS= ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS+1
                 nats=nats+1
                 tempEP=tempEP+vpn2*cm(ityp(i))/(3.0*bk)
              end if
           case(0)
              ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon=ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon&
                   &+vpn2*cm(ityp(i))/(3.0*bk)
                 ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS= ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS+1
              ecell(ixyze(1),ixyze(2),ixyze(3))%nIon= ecell(ixyze(1),ixyze(2),ixyze(3))%nIon+1
              nats=nats+1
              tempEP=tempEP+vpn2*cm(ityp(i))/(3.0*bk)
           end select
        end if
        !       write(6,'(I4,D21.12)')i,vpn2
        v2(ityp(i)) = v2(ityp(i))+vpn2
        if (ltpcel) then
           if(tempmaxat(ityp(i)).lt.vpn2)tempmaxat(ityp(i))=vpn2
        end if
        vx2(ityp(i),:) = vx2(ityp(i),:)+vp(:,i)**2
        if ((ltpcel).or.(tempstopcel.gt.0).or.(tcelec.gt.0))then
           tat=vpn2*cm(ityp(i))/(3.0*bk)
           tempc(ko)=tempc(ko)+tat/nato(ko)
           if(tat.gt.tempcm(ko))tempcm(ko)=tat
        end if
     end do


  end do

#if(PARA)
  call MPI_ALLREDUCE(v2,v2_glob,ntyp,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  v2=v2_glob
  call MPI_ALLREDUCE(vx2(1:ntyp,1:3),vx2_glob(1:ntyp,1:3),ntyp*3,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  vx2=vx2_glob
  if (allocated(tempc)) then
     call MPI_ALLREDUCE(tempc,tempc_tot,noxyz,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     tempc=tempc_tot
  end if
  if (l2T.eqv..true.) then
     call MPI_ALLREDUCE(ecell%tempIon,tempiontot,nex*ney*nez,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     ecell(:,:,:)%tempIon=tempiontot(:,:,:)
     niontot=0
     call MPI_ALLREDUCE(ecell%nIon,niontot,nex*ney*nez,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     ecell(:,:,:)%nIon=niontot(:,:,:)
     niontot=0
     call MPI_ALLREDUCE(ecell%nIonS,niontot,nex*ney*nez,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     ecell(:,:,:)%nIonS=niontot(:,:,:)

     call MPI_ALLREDUCE(tempEP,tempEptot,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     tempEP=tempEPtot
     call MPI_ALLREDUCE(nats,natstot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     nats=natstot
     

     deallocate(tempiontot)
     deallocate(niontot)
  end if
#endif

  if (l2T.eqv..true.) then
     do ixe=1,nex
        do iye=1,ney
           do ize=1,nez
              if (ecell(ixe,iye,ize)%nIonS.gt.0) then
                 ecell(ixe,iye,ize)%tempIon=ecell(ixe,iye,ize)%tempIon/ecell(ixe,iye,ize)%nIonS
              end if
           end do
        end do
     end  do
        tempEP=tempEP/nats
  end if



     do iti = 1, ntyp
        !debug            write(6,*)'iti' ,iti, na(iti)
        if (na(iti)==0) cycle
        temptyp(iti) = v2(iti)*cm(iti)/(3.0*na(iti)*bk)
        !        if (ltpcel) then 
        !           tempmaxat(iti)=tempmaxat(iti)*cm(iti)/(3.0*bk)
        !           write(6,*)'tempmaxat(iti)',iti,tempmaxat(iti)
        !        end if
        kine = kine+v2(iti)*cm(iti)/2.0
        !debug write(*,*) 'KINE WAS HERE ....', kine, rangml
        temp = temp+temptyp(iti)*na(iti)
     end do

#if(PARA)
     if( associated(free)) then
        temp = temp/float(imfree)
     else
        temp = temp/float(im_glob)
     end if
#else
     temp = temp/float(imd)
#endif
!     write(6,*)'ZL',zl
     if ((ltpcel).or.(tempstopcel.gt.0)) then
        maxTcel=0.
        !     	if (ltpcel) then
        !           write (6, *)
        !           write (6, *) '----------valeurs par cellules------------'
        !        endif

        do kx=0,nox-1
           do ky=0,noy-1
              do kz=0,noz-1
                 ko=1+kx+nox*(ky+noy*kz)
                 !                 pmc=0.0
                 !                              write(6,*)'dans la celulle ',ko
         !     	if (ltpcel)  write(6,'(A,I7,I5,3I4,2F12.2)')'CEL-TEMP ', it,ko,kx,ky,kz,tempc(ko),tempcm(ko)
             	if (ltpcel)  write(743,'(3I4,3F15.5,F15.5)') kx,ky,kz,kx*Zl(1)*1d8/nox,ky*Zl(2)*1d8/noy,kz*Zl(3)*1d8/noz,tempc(ko)
                 maxTcel=max(maxTcel,tempc(ko))

                 !                              write (6, '(A11,I4,A15,F12.2)') 'Cellule: ', ko, &
                 !                                   'Temperature: ', tempc(ko)
              enddo
           end do
        end do
        !        write(6,*)'CEL-TEMPM',maxTcel
        !        deallocate (tempc)
     endif


     return
   end subroutine calctemp
