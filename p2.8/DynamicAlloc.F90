subroutine DynamicalAllocationCell

  use gen_com_m
  use var_pot
  implicit none

  allocate(ncel(0:noxyz,0:26))
  allocate(nato(0:noxyz))
  allocate(last(natperc,0:noxyz))
  allocate(deltadist(3,0:26,noxyz))
  if (lTPcel.EQV..true.)then
         allocate(sigc(3,3,noxyz)); sigc(:,:,:noxyz)=0.
  endif
  if (iewald.ge.1) then
     allocate (tabv3(-ncoucx:ncoucx,-ncoucy:ncoucy,-ncoucz:ncoucz))
     allocate (tabf3(ntyp,-ncoucx:ncoucx,-ncoucy:ncoucy,-ncoucz:ncoucz))
  end if
 
  ncel(:noxyz,:26) = zero                 ! et petite initialisation

  if((ltpcel==.true.).or.(tempstopcel.gt.0)) then
     allocate (tempc(noxyz))
     allocate (tempcm(noxyz))
        allocate(celpm1(noxyz))
!        allocate(celpm2(noxyz))
        allocate(tm1(noxyz))
!        allocate(tm2(noxyz))
        allocate(celpp(noxyz))
!        allocate(celpp2(noxyz))
        allocate(tcp(noxyz))
!        allocate(tcp2(noxyz))
        allocate(pmc(noxyz))

        allocate (lprtcel(noxyz))

  end if
end subroutine DynamicalAllocationCell

subroutine Deallocatecel
  use gen_com_m
  use var_pot,only : iewald
  implicit none
  deallocate(ncel)
  deallocate(nato)
  deallocate(last)
  deallocate(deltadist)
  if(associated(sigc))deallocate(sigc)
  if (iewald.ge.1) then
     deallocate (tabv3)
     deallocate (tabf3)
  end if

end subroutine Deallocatecel

subroutine DeallocateAll

  use gen_com_m
  use var_pot
  implicit none

  deallocate(ncel)
  deallocate(nato)
  deallocate(last)
  deallocate(deltadist)
  if(associated(sigc))deallocate(sigc)
  if(associated(tabv3))deallocate(tabv3)
  if(associated(tabF3))deallocate(tabf3)
  deallocate(na)
  deallocate(ipo)
  deallocate(cm)
  deallocate(catom)
  deallocate(ty)
  deallocate(pot)
  deallocate(q)
  deallocate(rc)
  deallocate(lue_paire)
  deallocate(lue_typ)
  deallocate(lue_trip)
  deallocate(ro)
  deallocate(dip)
  deallocate(pm)
  deallocate(roff1)
  deallocate(roff2)
  deallocate(a_factor)
  deallocate(r8p)
  if(associated(ray))deallocate(ray)
  if(associated(bm))deallocate(bm)
  if(associated(shel))deallocate(shel)
  if(associated(Awat))deallocate(Awat)
  if(associated(Bwat))deallocate(Bwat)
  if(associated(qwat))deallocate(qwat)
  if(associated(rawat))deallocate(rawat)
  if(associated(potw))deallocate(potw)
  if(associated(bspw))deallocate(bspw)
  if(associated(cspw))deallocate(cspw)
  if(associated(dspw))deallocate(dspw)
  if(associated(eamrep))deallocate(eamrep)
  if(associated(eamrep_d))deallocate(eamrep_d)
  if(associated(eamglue))deallocate(eamglue)
  if(associated(eamglue_d))deallocate(eamglue_d)
  if(associated(eamrho))deallocate(eamrho)
  if(associated(eamrho_d))deallocate(eamrho_d)
  if(associated(lamb))deallocate(lamb)
  if(associated(gam))deallocate(gam)
  if(associated(cangle))deallocate(cangle)
  if(associated(coup3c))deallocate(coup3c)
  if(associated(ipo3c))deallocate(ipo3c)
  if(associated(coup3c2))deallocate(coup3c2)
  if(associated(l3ctyp))deallocate(l3ctyp)
  if(associated(l3cpair))deallocate(l3cpair)
  if(associated(coord))deallocate(coord)
  if(associated(digr))deallocate(digr)
  if(associated(fda))deallocate(fda)
  if(associated(nad))deallocate(nad)
  if(associated(nas))deallocate(nas)
  if(associated(nai))deallocate(nai)


  if(associated(lue_paire))deallocate (lue_paire)
  if(associated(lu_roff_pair))deallocate (lu_roff_pair)
  if(associated(lue_typ))deallocate (lue_typ)
  if(associated(typ_pot_pair))deallocate (typ_pot_pair)
  if(associated(lue_trip))deallocate (lue_trip)
  if(associated(rue_pair))deallocate (rue_pair)

  if(allocated(tempc))deallocate (tempc)
  if(allocated(tempcm))deallocate (tempcm)







end subroutine DeallocateAll
