module creadp_mod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE gen_com_m, ONLY:im_glob,it
  USE temp_com,only:imm,at,bg,im
  implicit none
contains
  subroutine creadp(xp, xpp, ityp,vp)

    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE defcdp
    USE var_pot, ONLY:ntyp,ty
    implicit none

    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: ic,j,ntry,iti,i
    integer :: itapp
    integer :: idep
    integer :: iposI,iint
    real(double) :: a1,a2,a3,c1,c2,c3,z1,r2,rd,z3,z2, edt
    real(double),dimension(3):: xdec, xavant,xapres,xposinttest
    real(double), dimension(1,3) :: cv
    real(double),dimension(:),allocatable:: edrat

    write(6,*)'creadp',it,ideftyp
    itapp=it-1
    select case (ideftyp)
    case(0)
       ntry=0
1      continue
       ntry=ntry+1
       ! tirer une position d'insertion



       select case (typint)
       case(0)

          call random_number(z1)
          iposI=1+Int(z1*nposI)
          !        write(6,*)iposI
          if (iposI.gt.nposI) idep=nposI
          !iposI est l'indice de la position int.
          !xposinttest est la position effective de l'int.

          xposinttest(:)=xposint(:,iposI)
          call cryst_to_cart (1, xposinttest(:), bg, -1) !cryst vers cart sur cv

       case(1)
          call random_number(z1)
          call random_number(z2)
          call random_number(z3)
          xposinttest(1)=z1
          xposinttest(2)=z2
          xposinttest(3)=z3
       end select

       ! verifier qu'elle est loin de tout atome

       call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst



       do j=1,im
          c1 = xposinttest(1)-xp(1,j)
          c2 = xposinttest(2)-xp(2,j)
          c3 = xposinttest(3)-xp(3,j)
          if (c1>0.5) c1 = c1-1.
          if (c1<(-0.5)) c1 = c1+1.
          if (c2>0.5) c2 = c2-1.
          if (c2<(-0.5)) c2 = c2+1.
          if (c3>0.5) c3 = c3-1.
          if (c3<(-0.5)) c3 = c3+1.

          cv(1,1) = c1
          cv(1,2) = c2
          cv(1,3) = c3
          call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
          r2 = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
          !           write(6,*)'toto',sqrt(r2),j


          rd=sqrt(r2)
          if (rd<dminins) then 
             !           write(6,*)'rate', rd, dminins
             call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
             goto 1                ! position proche d'un atome
          end if

       end do
       call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
       call cryst_to_cart (1, xposinttest(:), at, 1) !cryst vers cart sur cv

       ! choix du type de  l'atome à tirer si il y a des Ed non nuls
       edt=0.
       do iti=1,ntyp
          if (ed(iti).ne.0) Edt=edt+1./Ed(iti)
       end do
       !     write(6,*)
       if (Edt.ne.0) then
          if (allocated(edrat).EQV..false.) then
             allocate(edrat(ntyp))
             if (ed(1).ne.0) edrat(1)=1/(ed(1)*edt)
             !           write(6,*)'edrat (1)', edrat(1)
             do iti=2,ntyp
                if (ed(iti).ne.0) then
                   edrat(iti)=edrat(iti-1)+1/(ed(iti)*edt)
                else
                   edrat(iti)=edrat(iti-1)
                end if
                !              write(6,*)'edrat (iti)', edrat(iti)
             end do
          end if
          call random_number(z1)
          do iti=1,ntyp
             if (z1.le.edrat(iti)) then
                ioxdef=iti
                exit
             end if
          end do
          write(6,*)'Ed : INTRODUCTION PF de type ', ty(iti)
       end if

       ! deplacement acceptable
       ! tirer un atome
3      continue
       call random_number(z1)   
       idep=imin+Int(z1*(imax-imin-1))
       if (ioxdef.ne.0) then
          if (ioxdef.gt.0) then
             if (ityp(idep).ne.ioxdef) goto 3
          else
             if (ityp(idep).eq.ioxdef) goto 3
          end if
       end if
       !     select case (ioxdef)
       !     case(0)
       !     case(-2)
       !        if (ityp(idep)==2) goto 3
       !     case(2)
       !        if (ityp(idep).ne.2) goto 3
       !     case(-4)
       !        if (ityp(idep)==4) goto 3
       !     case(4)
       !        if (ityp(idep).ne.4) goto 3
       !     case(-5)
       !        if (ityp(idep)==5) goto 3
       !    case(5)
       !       if (ityp(idep).ne.5) goto 3




       xavant(:)=xp(:,idep)
       xdec(:)=xp(:,idep)-xpp(:,idep)
       xp(:,idep)=xposinttest(:)
       xapres(:)=xposinttest(:)
       xpp(:,idep)=xposinttest(:)-xdec(:)


       write(6,*)'INTRDUCTION PF de type ', ty(ityp(idep))
       write(6,*)'ntry',ntry
       !     write(6,*)
       write(6,*)'indice lac int', idep, iposI
       !     write(6,*)
       write(6,'(A,3F12.5)')'pos. lac.', xavant(1)*1.d8,xavant(2)*1.d8,xavant(3)*1.d8
       write(6,'(A,3F12.5)')'pos. int.', xapres(1)*1.d8,xapres(2)*1.d8,xapres(3)*1.d8
       write(6,*)


    case(1)
       !     write(6,*)'nintrodp',nintrodp
       do iint=1,nintrodp
          !     write(6,*)'nintrodp',iint,nintrodp
          ntry=0

11        continue
          ntry=ntry+1
          select case (typint)
          case(0)

             call random_number(z1)
             iposI=1+Int(z1*nposI)
             !        write(6,*)iposI
             if (iposI.gt.nposI) idep=nposI
             !iposI est l'indice de la position int.
             !xposinttest est la position effective de l'int.

             xposinttest(:)=xposint(:,iposI)
             call cryst_to_cart (1, xposinttest(:), bg, -1) !cryst vers cart sur cv

          case(1)
             call random_number(z1)
             call random_number(z2)
             call random_number(z3)
             xposinttest(1)=z1
             xposinttest(2)=z2
             xposinttest(3)=z3
          end select
          call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst
          do j=1,im
             c1 = xposinttest(1)-xp(1,j)
             c2 = xposinttest(2)-xp(2,j)
             c3 = xposinttest(3)-xp(3,j)
             if (c1>0.5) c1 = c1-1.
             if (c1<(-0.5)) c1 = c1+1.
             if (c2>0.5) c2 = c2-1.
             if (c2<(-0.5)) c2 = c2+1.
             if (c3>0.5) c3 = c3-1.
             if (c3<(-0.5)) c3 = c3+1.

             cv(1,1) = c1
             cv(1,2) = c2
             cv(1,3) = c3
             call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
             r2 = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
             !           write(6,*)'toto',sqrt(r2),j


             rd=sqrt(r2)
             if (rd<dminins) then 
                !           write(6,*)'rate', rd, dminins
                call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
                goto 11                ! position proche d'un atome
             end if

          end do


          if (rsphdef.ne.0) then
             c1 = xposinttest(1)-centresphdef(1)
             c2 = xposinttest(2)-centresphdef(2)
             c3 = xposinttest(3)-centresphdef(3)
             if (c1>0.5) c1 = c1-1.
             if (c1<(-0.5)) c1 = c1+1.
             if (c2>0.5) c2 = c2-1.
             if (c2<(-0.5)) c2 = c2+1.
             if (c3>0.5) c3 = c3-1.
             if (c3<(-0.5)) c3 = c3+1.

             cv(1,1) = c1
             cv(1,2) = c2
             cv(1,3) = c3
             call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
             r2 = cv(1,1)*cv(1,1)+cv(1,2)*cv(1,2)+cv(1,3)*cv(1,3)
             !           write(6,*)'toto',sqrt(r2),j


             rd=sqrt(r2)
             if (rd> rsphdef) then 
                !              write(6,*)'rate sphere', rd, rsphdef
                call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
                goto 11                ! position proche d'un atome
             end if

          end if
          call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart
          call cryst_to_cart (1, xposinttest(:), at, 1) !cryst vers cart sur cv

          im=im+1
          xp(:,im)= xposinttest(:)
          xpp(:,im)= xposinttest(:)
          vp(:,im)=0.
          ityp(im)=1
          write(6,*)'INT', im, xp(:,im)


          write(6,*)'nouvel im ', im
       end do
    end select


    itapp=it
    write(6,*)'outcdp'
    im_glob=im

    return
  end subroutine creadp
end module
