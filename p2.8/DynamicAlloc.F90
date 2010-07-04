subroutine DynamicalAllocationCell

  use gen_com_m
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

end subroutine DynamicalAllocationCell

subroutine alloc_typ

  use gen_com_m
  implicit none
  integer::i,j,k ,ic
  integer,save ::ncall=0

  ncall=ncall+1
  if (ncall==1) then
     rumax=0.
     allocate(na(ntyp))
     allocate(ipo(ntyp,ntyp))
     ! initialisation de ipo
     k = 0
     do i = 1, ntyp
        do j = i, ntyp
           k = k+1
           ipo(i,j) = k
           ipo(j,i) = k
        end do
     end do
     allocate(cm(ntyp)) ; allocate(catom(ntyp))

     allocate(ty(ntyp))
     allocate(pot(4,npair,0:ngrid+1))

     allocate(q(ntyp))
     q(:)=0
     allocate(rc(ntyp))
     rc(1:ntyp)=rclu(1:ntyp)*1.d-8
     allocate(zz(npair))
     allocate(lue_paire(npair))
     allocate(lu_roff_pair(npair))
     lu_roff_pair(:)=.false.

     allocate(rue_pair(npair))
     rue_pair(:)=0.


     allocate(lue_typ(npair))
     allocate(typ_pot_pair(npair))
     allocate(lue_trip(ntrip))
     lue_trip(:)=.false.;lue_paire(:)=.false.; lue_typ(:)=.false.
     allocate(ro(npair)); allocate(dip(npair));allocate(pm(npair))
     allocate(roff1(npair));allocate(roff2(npair))
     allocate(a_factor(npair));allocate(r8p(npair))
     allocate(h2sm(ntyp))
     if(iterdf.ge.0) then
        allocate(coord(ntyp,ntyp,nkmax))
        allocate(digr(ntyp,ntyp,nkmax))
        digr=0.d0
        !     allocate(coorpart(contmax,ntyp))
     endif
     if(iteangle.ge.0) then
        allocate(fda(ntyp,ntyp,ntyp,contmax))
     endif

     allocate(nad(ntyp)) ; allocate(nas(ntyp)) ; allocate(nai(ntyp))

  end if



  if (ipotentiel==0) then
     allocate(ray(ntyp)) ;allocate(bm(ntyp)) ; allocate(shel(ntyp))
  endif

  if (ipotentiel==5) then
     allocate(Dmorse(npair)) ;allocate(amorse(npair)) ; allocate(Remorse(npair))
  endif

  if (ipotentiel==6) then
     allocate(ietaij(npair)) ;allocate(capHij(npair))
     allocate(capDij(npair));allocate(capWij(npair))
  endif


  if (ipotentiel==2) then
     allocate(Awat(npair));allocate(Bwat(npair));allocate(pwat(npair))
     allocate(qwat(npair));allocate(rawat(npair)) ; allocate(rawat2(npair))
     allocate(potw(npair,0:ngrid+1))
     allocate(bspw(npair,0:ngrid+1));allocate(cspw(npair,0:ngrid+1))
     allocate(dspw(npair,0:ngrid+1))
     allocate (gz(0:ngrid+1))
     allocate(fcr(0:ngrid+1))

     allocate(bspg(0:ngrid+1))
     allocate(cspg(0:ngrid+1))
     allocate(dspg(0:ngrid+1))
     allocate(cspf(0:ngrid+1))
     allocate(bspf(0:ngrid+1))
     allocate(dspf(0:ngrid+1))


  endif


  if (ipotentiel .ge.10) then
     allocate(eamrep(4,npair,0:ngrid+1))
     if (ipotentiel==12) then
        allocate(eamrho(4,npair,0:ngrid+1))
     else
        allocate(eamrho(4,ntyp,0:ngrid+1))
     end if
     allocate(eamglue(4,ntyp,0:ngrid+1))
  end if

!  if (l3c) then
     allocate(lamb(ntrip)) ; allocate(cangle(ntrip))
     allocate(gam(ntrip,npair)) ; allocate(coup3c(ntrip,npair))
     allocate(coup3c2(ntrip,npair))
     allocate(ipo3c(ntyp,ntyp,ntyp));allocate(l3ctyp(ntyp))
     allocate(l3cpair(ntyp))
     allocate(c3c(ntrip))
     k=0
     do ic=1,ntyp
        do i=1,ntyp
           do j=i,ntyp
              k=k+1
              ipo3c(ic,i,j)=k
              ipo3c(ic,j,i)=k
           enddo
        enddo
     enddo

!  endif




  return
end subroutine alloc_typ

subroutine Deallocatecel
  use gen_com_m
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
  if(associated(eamglue))deallocate(eamglue)
  if(associated(eamrho))deallocate(eamrho)
  if(associated(lamb))deallocate(lamb)
  if(associated(gam))deallocate(gam)
  if(associated(cangle))deallocate(cangle)
  if(associated(coup3c))deallocate(coup3c)
  if(associated(ipo3c))deallocate(ipo3c)
  if(associated(coup3c2))deallocate(coup3c2)
  if(associated(l3ctyp))deallocate(l3ctyp)
  if(associated(l3cpair))deallocate(l3cpair)
  deallocate(h2sm)
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







end subroutine DeallocateAll
