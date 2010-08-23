subroutine inputtersoff
  ! NE PAS FAIRE DE MOYENNE POUR LES n,c ... etc
  use gen_com_m
  USE T_kind_param_m, ONLY:  double
  use force_tersoff_facteurs
  implicit none
  integer :: lupotin=95
  character ::  fnampotin*80
  integer::iti,itj,ipr,i,l,nprns,ntypr,j,rue
  real(double)::rc1, psilu,rof1m,rof2m,cmr
  real(double),dimension (:), pointer :: lambda1lu,lambda2lu,lambda3lu,Aterlu,Bterlu
  real(double),dimension (:), pointer :: Rterlu,Sterlu !, betalu,nterlu,cterlu, dterlu, hterlu
  character :: tyr*3
  !real(double),dimension (:), pointer :: lambda1,lambda2,lambda3,Ater,Bter
  !real(double),dimension (:), pointer :: Rter,Ster, beta,nter,cter, dter, hter
  integer,pointer :: typtyp(:)

 
  fnampotin = 'tersoff.potin'

  open(unit=lupotin, file=fnampotin, status='old')

  if (npotentiel .gt.1)then
     read(lupotin,*) ntypr,psilu
  else
     iewald=0; l3c=.false.; r3cm=0.  
     read(lupotin,*) ntyp,psilu
     npair=  ntyp*(ntyp+1)/2 ; ntrip= ntyp*ntyp *(ntyp+1)/2
     call  alloc_typ
     ntypr=ntyp
  end if

  allocate(lambda1lu(ntyp));allocate(lambda2lu(ntyp));allocate(lambda3lu(ntyp))
  allocate(Aterlu(ntyp));allocate(Bterlu(ntyp));allocate(Rterlu(ntyp))
  allocate(Sterlu(ntyp))

  allocate(beta(ntyp));allocate(nter(ntyp))
  allocate(cter(ntyp));allocate(dter(ntyp));allocate(hter(ntyp))
  allocate(deltater(ntyp))

  allocate(lambda1(npair));allocate(lambda2(npair));allocate(lambda3(npair))
  allocate(Ater(npair));allocate(Bter(npair));allocate(Rter(npair))
  allocate(Ster(npair));allocate(psi(npair));
  allocate (typtyp(ntypr))
  
  do i=1,ntypr
     if (npotentiel .gt.1)then
        read(lupotin,*) cmr,tyr,iti
     else
        iti=i
        read(lupotin,*) cmr,tyr
     end if
     typtyp(i)=iti
     if(lue_typ(iti).EQV..true.)then
        if (rang==0)write(6,*) 'type',iti,'deja lu ; verification de la cohérence'
        if (cmr*umass.ne.cm(iti))then
           if (rang==0)write(6,*) 'pb avec cm'
           stop
        end if
        if (tyr.ne.ty(iti))then
           if (rang==0)write(6,*) 'pb avec ty'
           stop
        end if
     end if
     if (associated(typ_and_pot))typ_and_pot(iti,ipotentiel)=.true.

     cm(iti)=cmr*umass;ty(iti)=tyr
     read(lupotin,*) lambda1lu(iti),lambda2lu(iti),lambda3lu(iti),Aterlu(iti),Bterlu(iti)
     read(lupotin,*) Rterlu(iti),Sterlu(iti)
     read(lupotin,*) beta(iti),nter(iti),cter(iti),dter(iti),hter(iti),deltater(iti)
     if (deltater(iti)==0) deltater(iti)=1.0/(2*nter(iti))
     !     if (rang==0) write(6,*) iti,deltater(iti)
  end do
  read(lupotin,*,end=456)nprns
  if (npotentiel.gt.1) then
     if (rang==0) write(6,*)'lecture supplement potin pas possible avec npotentiel >1'
     call endrun
  end if
  do i=1,ntypr
     read(lupotin,*)catom(i)
  end do
  write(6,*)'lecture de roff1 et roff2 pour certaines paires ',nprns
  do i = 1, nprns
     read (lupotin, *) l, rof1m, rof2m
     roff1(l)=rof1m*1.d-8
     roff2(l)=rof2m*1.d-8
     lu_roff_pair(l)=.true.
  end do
456 continue

  close (lupotin)

  psi(:)=1.0
  rue =0
  do i=1,ntypr
     iti=typtyp(i)
     do j=i,ntypr       
        itj=typtyp(j)
        ipr=ipo(iti,itj)
        if (rang==0) write(6,*)'paire l active  ipotentiel: ',ipr, ipotentiel

        lambda1(ipr)=0.5*(lambda1lu(iti)+lambda1lu(itj))
        lambda2(ipr)=0.5*(lambda2lu(iti)+lambda2lu(itj))
        lambda3(ipr)=0.5*(lambda3lu(iti)+lambda3lu(itj))

        Ater(ipr)=sqrt(Aterlu(iti)*Aterlu(itj))
        Bter(ipr)=sqrt(Bterlu(iti)*Bterlu(itj))

        Rter(ipr)=sqrt(Rterlu(iti)*Rterlu(itj))
        Ster(ipr)=sqrt(Sterlu(iti)*Sterlu(itj))
        if (itj.ne.iti) then 
           psi(ipr)=psilu
        end if
        select case (ipotentiel)
        case(13)
           rc1= Ster(ipr)*1.0d-8  !Lisa Porter 89
           !rc1= Rter(ipr)+Ster(ipr)   !Tersoff 88
        case(14,15)
           rc1=(rter(ipr)+ster(ipr))*1.0d-8
        end select
        rue_pair(ipr)=rc1
        typ_pot_pair(ipr)=ipotentiel
     end do
  end do

  !CONVERSION en CGS
  lambda1(:)=lambda1(:)*1.0d+8
  lambda2(:)=lambda2(:)*1.0d+8
  lambda3(:)=lambda3(:)*1.0d+8
  Ater(:)=Ater(:)*1.60219d-12
  Bter(:)=Bter(:)*1.60219d-12

  Rter(:)=Rter(:)*1.0d-8
  Ster(:)=Ster(:)*1.0d-8


  rumax=max(rumax,maxval(rue_pair))
!  if(rang==0)write(6,*)'rumax',rumax

  ldemitab =.false.

!if(allocated (typ_and_pot).eqv..false.), i.e. si npotentiel==1 
  if(associated(typ_and_pot).eqv..false.) then
     allocate (typ_and_pot(ntyp,npotmax))
     typ_and_pot(:,:)=.false.
     typ_and_pot(1:ntyp,ipotentiel)=.true.
  end if


end subroutine inputtersoff
