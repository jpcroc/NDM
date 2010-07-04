
! **************************************************************
subroutine tr_fire(vp, fp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  implicit none

  real(double)  :: vp(3,imm)
  real(double)  :: fp(3,imm)

  real(double),parameter:: finc=1.1, fdec=0.5,alph_start=0.1,f_alph=0.99,tstep_MM=10
  integer,parameter:: Nmin=5


  integer,save:: nstep,ncall=0
  real(double):: norme_de_fp,norme_de_vp
  real(double),save::alph,pscal,tstep_max
  real(double):: fpn(3,imm)
  integer::i

  ncall=ncall+1
  if (ncall==1)then
     alph=alph_start
     nstep=0
     tstep_max=tstep_MM*tstep
  end if

  norme_de_fp=0.0 ; norme_de_vp=0.; pscal=0.
  do i=1,im
     norme_de_fp=norme_de_fp+fp(1,i)**2+fp(2,i)**2+fp(3,i)**2
     !         norme_de_vp=norme_de_vp+vp(1,i)**2+vp(2,i)**2+vp(3,i)**2
     pscal=pscal+vp(1,i)*fp(1,i)+vp(2,i)*fp(2,i)+vp(3,i)*fp(3,i)
  end do

  norme_de_fp=sqrt(norme_de_fp)
  !      norme_de_vp=sqrt(norme_de_vp)
  fpn(:,:)=fp(:,:)/norme_de_fp

  if (pscal.gt.0) then
     vp(:,:)=alph*vp(:,:)+(1.-alph)*fpn(:,:)*sqrt(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
     nstep=nstep+1
     if (nstep.gt.5) then
        tstep=max(tstep*finc,tstep_max)
        alph=alph*f_alph
     end if
  else
     vp(:,:)=0.
     tstep=tstep*fdec
     alph=alph_start
     nstep=0
  end if
  write(6,'(A,4E14.5)')'tr_f', pscal,tstep,alph,norme_de_vp
  return
end subroutine tr_fire

