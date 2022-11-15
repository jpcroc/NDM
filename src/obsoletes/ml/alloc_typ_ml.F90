module alloc_typ_ml_mod
        implicit none
        contains
subroutine alloc_typ_ml

  use var_pot
  implicit none
  integer::i,j,k

     npair=ntyp*(ntyp+1)/2
     ntrip=ntyp*ntyp*(ntyp+1)/2
     rumax=0.
     if (associated(na)) deallocate(na) ; allocate(na(ntyp))

     if (associated(ipo)) deallocate(ipo); allocate(ipo(ntyp,ntyp))
     ! initialisation de ipo
     k = 0
     do i = 1, ntyp
        do j = i, ntyp
           k = k+1
           ipo(i,j) = k
           ipo(j,i) = k
        end do
     end do
      if (associated(cm))    deallocate(cm)    ; allocate(cm(ntyp))
      if (associated(catom)) deallocate(catom) ; allocate(catom(ntyp))

      if (associated(ty)) deallocate(ty); allocate(ty(ntyp))
      if (associated(pot)) deallocate(pot); allocate(pot(4,npair,0:ngrid+1))
!      if (associated(pot_d)) deallocate(pot_d); allocate(pot_d(4,npair,0:ngrid+1))

      if (associated(q)) deallocate(q); allocate(q(ntyp))
       q(:)=0
      if (associated(rc)) deallocate(rc); allocate(rc(ntyp))
     rc(1:ntyp)=rclu(1:ntyp)*1.d-8
      if (associated(zz)) deallocate(zz); allocate(zz(npair))
      if (associated(lue_paire)) deallocate(lue_paire); allocate(lue_paire(npair))
      if (associated(lu_roff_pair)) deallocate(lu_roff_pair); allocate(lu_roff_pair(npair))
      lu_roff_pair(:)=.false.

      if (associated(rue_pair)) deallocate(rue_pair); allocate(rue_pair(npair))
     rue_pair(:)=0.


      if (associated(lue_typ)) deallocate(lue_typ)
      allocate(lue_typ(npair))
      if (associated(typ_pot_pair)) deallocate(typ_pot_pair)
      allocate(typ_pot_pair(npair))
      if (associated(lue_trip)) deallocate(lue_trip)
      allocate(lue_trip(ntrip))
     lue_trip(:)=.false.;lue_paire(:)=.false.; lue_typ(:)=.false.
      if (associated(ro)) deallocate(ro)
      allocate(ro(npair))
      if (associated(dip)) deallocate(dip)
      allocate(dip(npair))
      if (associated(pm)) deallocate(pm)
      allocate(pm(npair))
      if (associated(roff1)) deallocate(roff1)
      allocate(roff1(npair))
      if (associated(roff2)) deallocate(roff2)
      allocate(roff2(npair))
      if (associated(a_factor)) deallocate(a_factor)
      allocate(a_factor(npair))
      if (associated(r8p)) deallocate(r8p)
      allocate(r8p(npair))

!     if(iterdf.ge.0) then
      if (associated(coord)) deallocate(coord)
         allocate(coord(ntyp,ntyp,nkmax))
      if (associated(digr)) deallocate(digr)
         allocate(digr(ntyp,ntyp,nkmax))
        digr=0.d0
      if (associated(fda)) deallocate(fda)
         allocate(fda(ntyp,ntyp,ntyp,contmax))

     if (associated(nad)) deallocate(nad)
     allocate(nad(ntyp))
     if (associated(nas)) deallocate(nas)
     allocate(nas(ntyp))
     if (associated(nai)) deallocate(nai)
     allocate(nai(ntyp))

!  end if



  return
end subroutine alloc_typ_ml
end module
