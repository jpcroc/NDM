
! *************************************************************
subroutine calctemp(temptyp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use var_pot
  use gen_com_m
  use tab_imm_m
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
  integer :: ic, i, iti, ko, i2,kx,ky,kz,koo
  real(double), dimension(ntyp) :: v2
  real(double) :: vpn2,pmc,tat
  real(double), dimension(ntyp,3) :: vx2
  ! ym      real(double), dimension(ntyp,nce) :: v2c
  real(double), dimension (:),allocatable ::tempc,tempcm
#if(PARA)
  real(double), dimension(ntyp) :: v2_glob
  real(double), dimension(ntyp,3) :: vx2_glob

#endif

  !-----------------------------------------------
  !
  !
  ! local variables
  if (ltpcel) then
     allocate (tempc(noxyz))
     allocate (tempcm(noxyz))
     tempc(:)=0.
     tempcm(:)=0.
  endif
  temp = 0.0
  kine = 0.0
  temptyp(:ntyp) = 0.0
  v2(:ntyp) = 0.0
  vx2(:ntyp,:) = 0.0


  do ko = 1, noxyz

     if (nato(ko)==0) cycle

#if(PARA)
     if (proc_cell(ko).ne.myid) cycle
#endif

     do i2 = 1, nato(ko)

        i = last(i2,ko)
        if (num_at_glob(i).gt.im_glob) cycle
        
        vpn2 = vp(1,i)**2+vp(2,i)**2+vp(3,i)**2
        !       write(6,'(I4,D21.12)')i,vpn2
        v2(ityp(i)) = v2(ityp(i))+vpn2
        
        vx2(ityp(i),:) = vx2(ityp(i),:)+vp(:,i)**2
        if (ltpcel==.true.) then
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
#endif


     !      enddo


     do iti = 1, ntyp
        !            write(6,*)'iti' ,iti
        if (na(iti)==0) cycle
        temptyp(iti) = v2(iti)*cm(iti)/(3.0*na(iti)*bk)
        kine = kine+v2(iti)*cm(iti)/2.0
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

     if (ltpcel) then

        write (6, *)
        write (6, *) '----------valeurs par cellules------------'

        do kx=0,nox-1
           do ky=0,noy-1
              do kz=0,noz-1
                 ko=1+kx+nox*(ky+noy*kz)
                 pmc=0.0
                 !                              write(6,*)'dans la celulle ',ko
                 write(6,'(A,I5,3I4,2F12.2)')'CEL-TEMP ', ko,kx,ky,kz,tempc(ko),tempcm(ko)
                 !                              write (6, '(A11,I4,A15,F12.2)') 'Cellule: ', ko, &
                 !                                   'Temperature: ', tempc(ko)
              enddo
           end do
        end do
        deallocate (tempc)
     endif


     return
   end subroutine calctemp
