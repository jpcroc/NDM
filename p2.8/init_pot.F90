module init_pot_mod
    USE input_pair_mod,only: input_pair
  USE inputtersoff_mod,only: inputtersoff
  USE calpoeam_mod,only: calpoeam
  USE calpo_mod,only: calpo
  USE tersoff_zbl_mod,only: tersoff_zbl
  USE gen_com_m, ONLY:firsttime_lammps,parallele,rang,umass,A2cm,rang
  USE var_pot, ONLY:npair,ntrip,r3cm,rumax,typ_and_pot,lpotentiel,l3c,npotmax,rue_pot,ipotentiel,ngrid,csive,npotentiel,&
       &typ_pot_pair,rue_pair,catom,cm,iewald,ipo,lu_roff_pair,lue_paire,lue_typ,ntyp,roff1,roff2,ty,typ_pot_pair,&
       &q,rue_lammps
  
  USE eam,only:inputeam
  USE eamerco,only:inputeamerco
  USE SMjuli,only:inputeamjl
  USE alloc_typ_mod,only: alloc_typ
  USE param_det_mod,only: param_det
  implicit none

contains
  ! **************************************************************
  subroutine init_pot
    integer::i,ipotcont
    
    if (npotentiel.gt.1)then
       ipotentiel=-1
       npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
       allocate(typ_and_pot(ntyp,npotmax))
       typ_and_pot(:,:)=.false.
       call  alloc_typ
       typ_pot_pair=0
    end if

    allocate(rue_pot(npotmax))
    rue_pot(:)=0.

    !>---------setting the potential------------------
    if (rang.eq.0) then
       write(6,*)
       write(6,*)'-*-*-*-*-*-*-*POTENTIELS*-*-*-*-*-'
       write(6,*)
    end if

    !#ifdef LAMMPS_VERSION
!    firsttime_lammps=.true.
    if ((ipotentiel==-10).or.(ipotentiel==-11))then
       call init_potential_simple(rue_lammps,rumax)
    else
       !#endif  

       do ipotcont=0,npotmax
          if(lpotentiel(ipotcont).EQV..true.) then
             ipotentiel=ipotcont
          else
             cycle
          end if

          if (ipotentiel.lt.10) then
             if (rang.eq.0) then
                write(6,*)
                write(6,*)'POTENTIEL DE PAIRES '
                write(6,*)
             end if

             call input_pair
          else
             select case(ipotentiel)
             case(10,16)
                if (rang.eq.0) then
                   write(6,*)
                   write(6,*)'POTENTIEL EAM'
                   if(ipotentiel==16)write(6,*)'+ CHARGE = CRG'
                   write(6,*)
                end if
                call inputeam(ntyp,npair,ntrip,cm,catom,ty,umass,rue_pot(ipotentiel),rumax,&
                     iewald,l3c,r3cm,roff1,roff2,typ_and_pot,npotmax,&
                     ipotentiel,typ_pot_pair,lue_typ,lue_paire,lu_roff_pair,&
                     npotentiel,ipo)
                do i=1,npair
                   if (typ_pot_pair(i)==ipotentiel) rue_pair(i)=rue_pot(ipotentiel)
                end do
             case(11)
                if (rang.eq.0) then
                   write(6,*)

                   write(6,*)'POTENTIEL EAM ERCOLESI'
                   write(6,*)
                end if
                call inputeamerco(ntyp,npair,ntrip,cm,catom,ty,umass,rue_pot(ipotentiel),rumax,&
                     iewald,l3c,rang,r3cm,roff1,roff2,typ_and_pot,npotmax,&
                     ipotentiel,typ_pot_pair,lue_typ,lue_paire,lu_roff_pair,&
                     npotentiel,ipo)
                do i=1,npair
                   if (typ_pot_pair(i)==ipotentiel) rue_pair(i)=rue_pot(ipotentiel)
                end do
             case(12)
                if (rang.eq.0) then
                   write(6,*)
                   write(6,*)'POTENTIEL Ju Li'
                   write(6,*)
                end if
                call inputeamjl(ntyp,npair,ntrip,cm,catom,ty,umass,rue_pot(ipotentiel),&
                     rumax,iewald,l3c,rang,r3cm,roff1,roff2,typ_and_pot,&
                     npotmax,ipotentiel,typ_pot_pair)
                rue_pair(:)=rue_pot(ipotentiel)
             case(13,14,15)
                !nguyen mettre input tersoff
                !        if (rang.eq.0) then
                !           if (rang==0)   write(6,*)'POTENTIEL tersoff.potin'
                !           if (ipotentiel==13) then
                !
                !           else
                !              write(6,*)'POTENTIEL tersoff.potin COUPURE MODIFIEE !!!!!!!!!!!!!!!!!!!!!!!!'
                !           end if
                !           write(6,*)
                !        end if
                if (rang.eq.0) then
                   write(6,*)
                   write(6,*)'POTENTIEL Tersoff-Brenner'
                   write(6,*)
                end if
                
                call inputtersoff

#ifdef ML
                ! MiLaDy
             case(20)
                if (rang.eq.0) then
                   write(6,*)
                   write(6,*)' ML ..... set-up MiLady potential'
                   write(6,*)
                end if
                !This comes with MiLaDy package
                call md_init_potential_ml
#endif
             end select

          endif
       end do

       !#ifdef LAMMPS_VERSION
    endif
    !#endif


    return
  end subroutine init_pot




  subroutine init_pot2(boxndm,immT)
    use boxconfig,only:box_config
    USE calpo_ew_mod,only: calpo_ew
    implicit none
#ifdef ML

#else
    
    type(box_config)::boxndm
    integer,intent(in)::immT
    integer::ipotcont

    if (rang==0)then
          write(6,*)
       write(6,*)' -------------------------------------------------------------------'
       write(6,*)'             definition des rayons de coupure'
    end if
    call param_det(boxndm)
       ! rumax défini en ce point
#endif 

    do ipotcont=0,npotmax
       if(lpotentiel(ipotcont).EQV..true.) then
          ipotentiel=ipotcont
       else
          cycle
       end if

       select case(ipotentiel)
       case(0:9)
          call calpo
          if (iewald==1.or.iewald==2) then
             call calpo_ew(boxndm,immT)
          end if

       case(10,11,12,16)
          call calpoeam
          if ((ipotentiel==16).and.(iewald.gt.0)) then
             call calpo_ew(boxndm,immT)
          end if

       case(13,14,15)
          if (.not.parallele) then
             if (maxval(roff1).gt.0) call tersoff_zbl
          end if
       end select
    end do
    if (npotentiel.gt.1) then
       !     write(6,*)
       !     write(6,*)'BILAN DES POTENTIELS'
       !     do i=1,ntyp
       !        do j=1,ntyp
       !           l=ipo(i,j)
       !           write(6,*)
       !           write(6,*)'paire i-j l',i,j,l
       !           write(6,*)'lue paire typ_pot_pair coupure'
       !           write(6,*)lue_paire(l),typ_pot_pair(l),rue_pair(l)*1d8
       !        end do
       !     end do
       if (rang==0) then
          write(6,*)
          write(6,*)'decoupage en cellule suivant'
          write(6,*)'rumax',rumax*1d8
          write(6,*)
       end if
    end if


    !<---------end setting the potential---------------

  end subroutine init_pot2
    
  subroutine init_potential_simple(rue,rum)
    USE T_kind_param_m, ONLY:  double
!    USE gen_com_m, ONLY: rang,A2cm,umass
!    USE var_pot, ONLY: ntyp, npair, ntrip,cm,catom, ty,rue_pair,ipotentiel,q
    implicit none
    integer :: i,error, beggin,  endding,lupotin
    character ::  fnampotin*80
    real(double)::rue,rum

    fnampotin = 'simple.potin'
    lupotin = 95
    open(unit=lupotin, file=fnampotin, status='old')
    read(lupotin,*)ntyp
    npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
    call  alloc_typ
    read(lupotin,*) rue
    rue=rue*A2cm
    rue_pair(:)=rue
    select case (ipotentiel)
    case(-10)
       do i = 1, ntyp
          read (lupotin,*) cm(i),catom(i),ty(i)
          if (rang/=0) cycle
          write (6, '(I4,2F9.3,A5)') i, cm(i),catom(i),ty(i)
       end do
    case(-11)
           do i = 1, ntyp
          read (lupotin,*) cm(i),catom(i),ty(i),q(i)
          if (rang/=0) cycle
          write (6, '(I4,2F9.3,A5,F9.3)') i, cm(i),catom(i),ty(i),q(i)
       end do

    end select
    cm(:ntyp) = cm(:ntyp)*umass
    rumax=max(rumax,rue)
    close (lupotin)

  end subroutine init_potential_simple
end module init_pot_mod
