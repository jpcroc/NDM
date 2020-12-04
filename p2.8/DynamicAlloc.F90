module dynalloccell
  use temp_com,only: ncel,nato,atincel,deltadist,sigc,tempc,tempcm,noxyz
  USE gen_com_m, ONLY:l2t,lsigatcel,ltpcel,tempstopcel,tabv3,tabf3,natchk, lprtcel,zero,elossCel
  USE var_pot, ONLY:iewald,ncoucx,ncoucy,ncoucz,ntyp,na,cm,ipo,catom,ty,pot,rc,lue_paire,lue_typ,lue_trip,dip,pm,roff1,&
       &roff2,a_factor,r8p,ray,bm,shel,awat,bwat,qwat,potw,bspw,cspw,bspw,eamrep,eamrep_d,eamglue,eamglue_d,eamrho,eamrho_d,&
       &lamb,gam,cangle, coup3c,ipo3c, coup3c2,l3ctyp,l3cpair,coord,digr,fda,nad,nas,nai,lu_roff_pair,lue_typ,&
       &typ_pot_pair,lue_trip,rue_pair,ipo,q,ro,rawat,dspw
  implicit none 
contains
  subroutine DynamicalAllocationCell



    USE eloss, ONLY :tcelec
    implicit none

!    allocate(ncel(0:noxyz,0:26))
!    allocate(nato(0:noxyz))
!    allocate(atincel(natperc,0:noxyz))
!    allocate(deltadist(3,0:26,noxyz))
    if (lTPcel.EQV..true.)then
!       allocate(sigc(3,3,noxyz)); sigc(:,:,:noxyz)=0.
    endif
    if (iewald.ge.1) then
       allocate (tabv3(-ncoucx:ncoucx,-ncoucy:ncoucy,-ncoucz:ncoucz))
       allocate (tabf3(ntyp,-ncoucx:ncoucx,-ncoucy:ncoucy,-ncoucz:ncoucz))
    end if

    ncel(:noxyz,:26) = zero                 ! et petite initialisation

    if((ltpcel.eqv..true.).or.(tempstopcel.gt.0).or.(tcelec.gt.0)) then
!       allocate (tempc(noxyz))
!       allocate (tempcm(noxyz))
!       allocate(celpm1(noxyz))
       !        allocate(celpm2(noxyz))
!       allocate(tm1(noxyz))
       !        allocate(tm2(noxyz))
!       allocate(celpp(noxyz))
       !        allocate(celpp2(noxyz))
!       allocate(tcp(noxyz))
       !        allocate(tcp2(noxyz))
!       allocate(pmc(noxyz))

!       allocate (lprtcel(noxyz))

    end if
    if (l2T.eqv..true.)then
 !      if (.not.allocated(tempc))     allocate (tempc(noxyz))
 !      allocate (elossCel(noxyz))
    endif
    if (lsigatcel.eqv..true.) then
 !      allocate (patcel(noxyz))
 !      allocate (patcelmax(noxyz))
 !      allocate (sigatcel(3,3,noxyz))
 !      allocate (natchk(noxyz))
    end if


  end subroutine DynamicalAllocationCell

subroutine Deallocatecel
  USE gen_com_m, ONLY:

  implicit none
  if(allocated(ncel))deallocate(ncel)
  if(allocated(nato))deallocate(nato)
  if(allocated(atincel))deallocate(atincel)
  if(allocated(deltadist))deallocate(deltadist)
  if(allocated(sigc))deallocate(sigc)
  if (iewald.ge.1) then
     deallocate (tabv3)
     deallocate (tabf3)
  end if

end subroutine Deallocatecel

subroutine DeallocateAll

  USE gen_com_m, ONLY:

  implicit none

  if(allocated(ncel)) deallocate(ncel)
  if(allocated(nato)) deallocate(nato)
  if(allocated(atincel)) deallocate(atincel)
  if(allocated(deltadist)) deallocate(deltadist)
  if(allocated(sigc))deallocate(sigc)
  if(allocated(tabv3))deallocate(tabv3)
  if(allocated(tabF3))deallocate(tabf3)
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
  if(allocated(ray))deallocate(ray)
  if(allocated(bm))deallocate(bm)
  if(allocated(shel))deallocate(shel)
  if(allocated(Awat))deallocate(Awat)
  if(allocated(Bwat))deallocate(Bwat)
  if(allocated(qwat))deallocate(qwat)
  if(allocated(rawat))deallocate(rawat)
  if(allocated(potw))deallocate(potw)
  if(allocated(bspw))deallocate(bspw)
  if(allocated(cspw))deallocate(cspw)
  if(allocated(dspw))deallocate(dspw)
  if(allocated(eamrep))deallocate(eamrep)
  if(allocated(eamrep_d))deallocate(eamrep_d)
  if(allocated(eamglue))deallocate(eamglue)
  if(allocated(eamglue_d))deallocate(eamglue_d)
  if(allocated(eamrho))deallocate(eamrho)
  if(allocated(eamrho_d))deallocate(eamrho_d)
  if(allocated(lamb))deallocate(lamb)
  if(allocated(gam))deallocate(gam)
  if(allocated(cangle))deallocate(cangle)
  if(allocated(coup3c))deallocate(coup3c)
  if(allocated(ipo3c))deallocate(ipo3c)
  if(allocated(coup3c2))deallocate(coup3c2)
  if(allocated(l3ctyp))deallocate(l3ctyp)
  if(allocated(l3cpair))deallocate(l3cpair)
  if(allocated(coord))deallocate(coord)
  if(allocated(digr))deallocate(digr)
  if(allocated(fda))deallocate(fda)
  if(allocated(nad))deallocate(nad)
  if(allocated(nas))deallocate(nas)
  if(allocated(nai))deallocate(nai)


  if(allocated(lue_paire))deallocate (lue_paire)
  if(allocated(lu_roff_pair))deallocate (lu_roff_pair)
  if(allocated(lue_typ))deallocate (lue_typ)
  if(allocated(typ_pot_pair))deallocate (typ_pot_pair)
  if(allocated(lue_trip))deallocate (lue_trip)
  if(allocated(rue_pair))deallocate (rue_pair)

  if(allocated(tempc))deallocate (tempc)
  if(allocated(tempcm))deallocate (tempcm)

  return


end subroutine DeallocateAll
end module
