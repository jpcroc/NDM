module calctemp_mod
   USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m, ONLY:  double
  USE var_pot, ONLY:ntyp,cm
  USE gen_com_m, ONLY:erg2ev,tempEP,bk,l2t,lspaceNDM,rang
  USE elec_cell, ONLY: ecell,i2T,nex,ney,nez,nox_2_nex
  USE eloss, ONLY : tcelec,ecelec
  USE atomconfig,only: atom_config_d
  USE cellconfig,only : cell_config
  USE atomconfig,only : atom_config_arps
  use calcfvp_mod,only:calcfvp
#ifdef PARA
  USE Tpara,only:myidsp,nprocspace,comm_space,nprocs
#else
    USE Tpara,only:nprocs
#endif

  ! *************************************************************

  implicit none
contains
  ! *************************************************************
  subroutine calctemp(temp,kine,atcf, cellcf,latcomp)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------

    class(atom_config_d),intent(in)::atcf
    type(cell_config),intent(inout)::cellcf
    real(double),intent(out)::temp,kine
    logical, optional::latcomp

    logical:: latc=.false.

    integer ::  i, ko, i2,kx,ixe,iye,ize
    real(double) :: sumtat2
    real(double) :: vpn2,tat,ekin,fvp
    !  real(double), dimension(ntyp,3) :: vx2
    integer::ixyze(3),nats,nat
    ! ym      real(double), dimension(ntyp,nce) :: v2c
    !  real(double), dimension (:),allocatable ::tempc,tempcm

    if (present(latcomp)) latc=latcomp
    if (nprocs==1) latc=.true.
    nat=0
    if (latc) then
       if(cellcf%icaltabt.ne.atcf%icaltabt) then
          write (6,*)'incoherence dans icaltabt calctemp'
          write(6,*)'cell atcf', cellcf%icaltabt,atcf%icaltabt
          call arret_ndm
       end if

       if ((cellcf%ltpcel).or.(tcelec.gt.0)) then
          !     allocate (tempc(noxyz))
          !     allocate (tempcm(noxyz))
          cellcf%tempc(:)=0.
       endif
       temp = 0.0
       kine = 0.0
       sumtat2 = 0.0
       !  vx2(:ntyp,:) = 0.0


       if (l2t)then
          ecell(:,:,:)%tempIon=0
          ecell(:,:,:)%nIon=0
          ecell(:,:,:)%nIonS=0
          nats=0
       end if

       tempEP=0
       do ko = 1, cellcf%noxyz
          if ((cellcf%nato(ko)==0).or.(cellcf%isghost(ko))) cycle

          if (L2T)     call nox_2_nex(ko,ixyze,cellcf)

          do i2 = 1, cellcf%nato(ko)

             i = cellcf%atincel(i2,ko)
             vpn2 = atcf%vp(1,i)**2+atcf%vp(2,i)**2+atcf%vp(3,i)**2
             select type(atcf)
             type is (atom_config_arps)
                call calcfvp(fvp,atcf%ityp(i),atcf%vp(:,i),atcf%mov(i))
                vpn2=vpn2*fvp
             end select
             kine=kine+vpn2*0.5*cm(atcf%ityp(i))
             tat=vpn2*cm(atcf%ityp(i))/(3.0*bk)
             sumtat2 = sumtat2+tat
             nat=nat+1
             !calculation of ionic temperature and number of ions in the electronic cell (only slow moving ions)  
             if (l2T.eqv..true.) then
                select case(i2t)
                case(1)
                   ecell(ixyze(1),ixyze(2),ixyze(3))%nIon= ecell(ixyze(1),ixyze(2),ixyze(3))%nIon+1
                   ekin=0.5*erg2ev*vpn2*cm(atcf%ityp(i))
                   if (ekin.lt.Ecelec) then
                      ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon=ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon&
                           &+vpn2*cm(atcf%ityp(i))/(3.0*bk)
                      ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS= ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS+1
                      nats=nats+1
                      tempEP=tempEP+vpn2*cm(atcf%ityp(i))/(3.0*bk)
                   end if
                case(0)
                   ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon=ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon&
                        &+vpn2*cm(atcf%ityp(i))/(3.0*bk)
                   ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS= ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS+1
                   ecell(ixyze(1),ixyze(2),ixyze(3))%nIon= ecell(ixyze(1),ixyze(2),ixyze(3))%nIon+1
                   nats=nats+1
                   tempEP=tempEP+vpn2*cm(atcf%ityp(i))/(3.0*bk)
                end select
             end if
             if ((cellcf%ltpcel).or.(tcelec.gt.0))then
                cellcf%tempc(ko)=cellcf%tempc(ko)+tat/cellcf%nato(ko)
             end if
          end do

       end do
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

       if (nat.ne.atcf%im) then
          write(6,*)'NAT NE atcf%im STOP',nat,atcf%im
          call arret_ndm
       end if
       temp = sumtat2/float(atcf%im)


    else  !LATC/LATCOMP=.false.

       if(cellcf%icaltabt.ne.atcf%icaltabt) then
          write (6,*)'incoherence dans icaltabt calctemp'
          write(6,*)'cell atcf', cellcf%icaltabt,atcf%icaltabt
          call arret_ndm
       end if

       if ((cellcf%ltpcel).or.(tcelec.gt.0)) then
          !     allocate (tempc(noxyz))
          !     allocate (tempcm(noxyz))
          cellcf%tempc(:)=0.
       endif
       temp = 0.0
       kine = 0.0
       sumtat2 = 0.0
       !  vx2(:ntyp,:) = 0.0
       if (l2t)then
          ecell(:,:,:)%tempIon=0
          ecell(:,:,:)%nIon=0
          ecell(:,:,:)%nIonS=0
          nats=0
       end if
       tempEP=0
       do ko = 1, cellcf%noxyz
          if ((cellcf%nato(ko)==0).or.(cellcf%isghost(ko))) cycle

#ifdef PARA

          if ((lspaceNDM).and.(cellcf%proc_cell(ko).ne.myidsp)) cycle
#endif

          if (L2T)  call nox_2_nex(ko,ixyze,cellcf)

          do i2 = 1, cellcf%nato(ko)
             nat=nat+1
             i = cellcf%atincel(i2,ko)
             vpn2 = atcf%vp(1,i)**2+atcf%vp(2,i)**2+atcf%vp(3,i)**2
             select type(atcf)
             type is (atom_config_arps)
                call calcfvp(fvp,atcf%ityp(i),atcf%vp(:,i),atcf%mov(i))
                vpn2=vpn2*fvp
             end select

             kine=kine+vpn2*0.5*cm(atcf%ityp(i))
             tat=vpn2*cm(atcf%ityp(i))/(3.0*bk)
             sumtat2 = sumtat2+tat

             !calculation of ionic temperature and number of ions in the electronic cell (only slow moving ions)  
             if (l2T.eqv..true.) then
                select case(i2t)
                case(1)
                   ecell(ixyze(1),ixyze(2),ixyze(3))%nIon= ecell(ixyze(1),ixyze(2),ixyze(3))%nIon+1
                   ekin=0.5*erg2ev*vpn2*cm(atcf%ityp(i))
                   if (ekin.lt.Ecelec) then
                      ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon=ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon&
                           &+vpn2*cm(atcf%ityp(i))/(3.0*bk)
                      ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS= ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS+1
                      nats=nats+1
                      tempEP=tempEP+vpn2*cm(atcf%ityp(i))/(3.0*bk)
                   end if
                case(0)
                   ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon=ecell(ixyze(1),ixyze(2),ixyze(3))%tempIon&
                        &+vpn2*cm(atcf%ityp(i))/(3.0*bk)
                   ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS= ecell(ixyze(1),ixyze(2),ixyze(3))%nIonS+1
                   ecell(ixyze(1),ixyze(2),ixyze(3))%nIon= ecell(ixyze(1),ixyze(2),ixyze(3))%nIon+1
                   nats=nats+1
                   tempEP=tempEP+vpn2*cm(atcf%ityp(i))/(3.0*bk)
                end select
             end if
             if ((cellcf%ltpcel).or.(tcelec.gt.0))then
                cellcf%tempc(ko)=cellcf%tempc(ko)+tat/cellcf%nato(ko)
             end if
          end do

       end do
#ifdef PARA
       if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
          call comm_space%sum(kine)
          call comm_space%sum(nat)
          call comm_space%sum(sumtat2)
          if (allocated(cellcf%tempc)) then
             call comm_space%sum(cellcf%tempc)
          end if
          if (l2T.eqv..true.) then
             call comm_space%sum(ecell%tempion)
             call comm_space%sum(ecell%nIon)
             call comm_space%sum(ecell%nIonS)
             call comm_space%sum(tempEP)
             call comm_space%sum(nats)
          end if
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

#ifdef PARA
       if ((nprocspace.gt.1).and.(lspacendm.eqv..true.))then 
          temp = sumtat2/nat
          !deallocate(tempc_tot)
       else
          if (nat.ne.atcf%im) then
             write(6,*)'NAT NE atcf%im STOP'
             call arret_ndm
          end if
          temp = sumtat2/float(atcf%im)
       end if
       !     end if

#else
       if (nat.ne.atcf%im) then
          write(6,*)'NAT NE atcf%im STOP'
          call arret_ndm
       end if
       temp = sumtat2/float(atcf%im)
#endif
    end if

    return
  end subroutine calctemp
end module calctemp_mod
