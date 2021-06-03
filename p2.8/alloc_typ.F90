module alloc_typ_mod
  USE var_pot, ONLY:ntyp,rumax,lpotentiel,ipo3c,l3cpair,l3ctyp,coup3c2,gam,lamb,&
       &cangle,eamglue,eamglue_d,eamrho,eamrho_d,eamrep,eamrep_d,dspf,bspf,cspf,dspg,cspg,&
       &cm,catom, ipo,ty,q,pot_d,rc,zz,lue_paire,lu_roff_pair,rue_pair,lue_typ,lue_trip,&
       &ro,dip,pm,r8p,roff1,roff2,pot,rclu,a_factor,fda,gamlt,shel,&
       &Dmorse,amorse,Remorse,ietaij,capHij,capwij,capDij,Awat,Bwat,pwat,qwat,rawat,&
       &rawat2,potw,bspw,cspw,dspw,gz,fcr,bspg,contmax,ngrid,nkmax,npair,ntrip,&
       &typ_pot_pair,ray,bm,coup3c,c3c
  implicit none 
contains

  subroutine alloc_typ


    USE gen_com_m, ONLY: llangevin,gamlg
    implicit none
    integer::i,j,k ,ic
    integer,save ::ncall=0
    integer :: ipot_loc

    ncall=ncall+1
    if (ncall==1) then
       rumax=0.
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
       allocate(pot_d(4,npair,0:ngrid+1))

       allocate(q(ntyp))
       q(:)=0
       allocate(rc(ntyp))
       rc(1:ntyp)=rclu(1:ntyp)*1.d-8
       allocate(zz(npair))
       zz=0
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

       !     if(iterdf.ge.0) then
       !     allocate(coorpart(contmax,ntyp))
       !     endif
       !     if(iteangle.ge.0) then
       allocate(fda(ntyp,ntyp,ntyp,contmax))


       if (llangevin) then 
          allocate (gamlt(ntyp))
          if (gamlg.ge.0) gamlt(:)=gamlg
       end if
    endif



    !  end if


    ipot_loc=0
    if (lpotentiel(ipot_loc).eqv..true.) then
       allocate(ray(ntyp)) ;allocate(bm(ntyp)) ; allocate(shel(ntyp))
    endif

    ipot_loc=5
    if (lpotentiel(ipot_loc).eqv..true.) then
       allocate(Dmorse(npair)) ;allocate(amorse(npair)) ; allocate(Remorse(npair))
    endif

    ipot_loc=6
    if (lpotentiel(ipot_loc).eqv..true.) then
       allocate(ietaij(npair)) ;allocate(capHij(npair))
       allocate(capDij(npair));allocate(capWij(npair))
    endif

    ipot_loc=2
    if (lpotentiel(ipot_loc).eqv..true.) then
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

    do ipot_loc=10,12
       if (lpotentiel(ipot_loc).eqv..true.) then
          allocate(eamrep(4,npair,0:ngrid+1))
          allocate(eamrep_d(4,npair,0:ngrid+1))
          if (ipot_loc==12) then
             allocate(eamrho(4,npair,0:ngrid+1))
             allocate(eamrho_d(4,npair,0:ngrid+1))
          else
             allocate(eamrho(4,ntyp,0:ngrid+1))
             allocate(eamrho_d(4,ntyp,0:ngrid+1))
          end if
          allocate(eamglue(4,ntyp,0:ngrid+1))
          allocate(eamglue_d(4,ntyp,0:ngrid+1))
       end if
    end do
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
end module alloc_typ_mod
