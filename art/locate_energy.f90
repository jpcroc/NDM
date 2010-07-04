!****************************************************************
      subroutine locate_energy ()

      use defs
      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use tab_imm_m
      use random_art
      implicit none






!input output  decalarations for the ART parameters
!   integer, intent(in) :: natoms
!   integer, intent(in), dimension(natoms):: types
!   real(8), intent(in), dimension(3*natoms), target:: pos

!   real(8), dimension(:), pointer :: x, y, z

   integer::i,J
   real(8)::eatommoy(ntyp),energyth,energythM
   real(8):: decat
   integer::ndef,idef,inddef(natoms)
      real(8)                       :: ran3
!    x => pos(1:NATOMS)
!    y => pos(NATOMS+1:2*NATOMS)
!    z => pos(2*NATOMS+1:3*NATOMS)
    do i=1,NATOMS 
      xp(1,i)=POS(i)/angst
      J=I+NATOMS
      xp(2,i)=POS(J)/angst
      J=I+2*NATOMS
      xp(3,i)=POS(J)/angst
    enddo

!    write(*,*) 'at', at

!    write(*,*) 'bg', bg
    
    
    
    call caltabt 
    call caltabi     
    call calfo 
    
    selected_atoms(:)=0
    
    eatommoy(:)=0
    do i=1,im
        eatommoy(ityp(i))=eatommoy(ityp(i))+eatom(i)/na(ityp(i))
     end do
     
    ndef=0
    energythM=0
    do i=1,natoms
       energyth=(eatom(i)-eatommoy(ityp(i)))*erg2ev
!       write(6,*) 'energyth, ENERGYDIFFERENCE'
!       write(6,*) energyth, ENERGYDIFFERENCE
       if (energyth.gt.ENERGYDIFFERENCE) then

          ndef=ndef+1
!          write(6,*)'i,ndef ',i,ndef,(eatom(i)-eatommoy(ityp(i)))*erg2eV
          inddef(ndef)=i
	  selected_atoms(ndef)=i
       end if
       if (energyth.gt.energythM) energythM=energyth
    end do

    if (ndef==0) then
       write(6,*)'energydifference est trop grand => ndef=0' 
       stop
    else
       write(6,*)'energie seuil ,max = ',energydifference,energythM
       write(6,*)'nombre de sites de depart possibles = ',ndef
       selected_atoms_max=ndef
   end if

   idef = INT ( ndef * ran3() + 1)    ! Between 1 and ndef
   write(6,*)idef,'eme site choisi aléatoirement = atome no ', inddef(idef)
   preferred_atom=inddef(idef)   

       end subroutine locate_energy

!****************************************************************
