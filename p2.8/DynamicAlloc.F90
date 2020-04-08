module dynalloccell
  USE gen_com_m, ONLY:l2t,lsigatcel,ltpcel,tempstopcel,tabv3,tabf3,sigc,tempcm,deltadist,nato,last,&
       &ncel,natchk,noxyz,pmc, tcp, celpp,tcp, lprtcel,tempc,tm1,patcelmax,zero,elossCel,&
       &patcelmax, celpm1,patcel,sigatcel,natperc
  USE var_pot, ONLY:iewald,ncoucx,ncoucy,ncoucz,ntyp,na,cm,ipo,catom,ty,pot,rc,lue_paire,lue_typ,lue_trip,dip,pm,roff1,&
       &roff2,a_factor,r8p,ray,bm,shel,awat,bwat,qwat,potw,bspw,cspw,bspw,eamrep,eamrep_d,eamglue,eamglue_d,eamrho,eamrho_d,&
       &lamb,gam,cangle, coup3c,ipo3c, coup3c2,l3ctyp,l3cpair,coord,digr,fda,nad,nas,nai,lu_roff_pair,lue_typ,&
       &typ_pot_pair,lue_trip,rue_pair,ipo,q,ro,rawat,dspw
  implicit none 
contains
  subroutine DynamicalAllocationCell



    USE eloss, ONLY :tcelec
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

    if((ltpcel.eqv..true.).or.(tempstopcel.gt.0).or.(tcelec.gt.0)) then
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
    if (l2T.eqv..true.)then
       if (.not.allocated(tempc))     allocate (tempc(noxyz))
       allocate (elossCel(noxyz))
    endif
    if (lsigatcel.eqv..true.) then
       allocate (patcel(noxyz))
       allocate (patcelmax(noxyz))
       allocate (sigatcel(3,3,noxyz))
       allocate (natchk(noxyz))
    end if


  end subroutine DynamicalAllocationCell

subroutine Deallocatecel
  USE gen_com_m, ONLY:

  implicit none
  if(associated(ncel))deallocate(ncel)
  if(associated(nato))deallocate(nato)
  if(associated(last))deallocate(last)
  if(associated(deltadist))deallocate(deltadist)
  if(associated(sigc))deallocate(sigc)
  if (iewald.ge.1) then
     deallocate (tabv3)
     deallocate (tabf3)
  end if

end subroutine Deallocatecel

subroutine DeallocateAll

  USE gen_com_m, ONLY:

  implicit none

  if(associated(ncel)) deallocate(ncel)
  if(associated(nato)) deallocate(nato)
  if(associated(last)) deallocate(last)
  if(associated(deltadist)) deallocate(deltadist)
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

  return


end subroutine DeallocateAll
end module
