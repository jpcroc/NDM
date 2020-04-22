module contrainte

  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:imm,imm,im,dmtype,im,at,bg
  USE var_pot, ONLY:cm
  implicit none

  ! Constraint
  real(double),allocatable,save, private :: dg0(:,:)

  ! Tell if constraint is applied on reduced or cartesian coordinates
  LOGICAL, save, private :: reduced

  ! Variables USEd when reduced=.true.
  real(double), dimension(1:3,1:3), save, private :: inv_at, trans_at, trans_inv_at


contains


  ! **************************************************************
  subroutine initcontr(xp, xpp, vp, ax,ityp )
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    integer:: ityp(imm)
    integer :: i
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------

    real(double) :: masstot, sigg0(3)
    real(double) :: dg2, contrval
    real(double)  :: xc(3,imm)    
    CHARACTER(len=20) :: coordinate
    INTEGER :: io

    !-----------------------------------------------

    write(6,*)'entree initcontr'

    ! Lecture de la contrainte, i.e. la différence de déplacement
    ! entre l'état final et l'état initial
    allocate (dg0(3,imm))
    open (98,file='contrainte',action='read',status='old')
    do i=1,im
       read(98,*) dg0(1:3,i)
    end do
    READ(98,*,iostat=io) coordinate
    CLOSE(98)
    IF (io.NE.0) THEN
            reduced=.FALSE.
            WRITE(6,'(a)') 'Constraint is applied on cartesian coordinates (DEFAULT)'
    ELSE
            SELECT CASE(coordinate)
            CASE("CARTESIAN")
                    reduced=.FALSE.
                    WRITE(6,'(a)') 'Constraint is applied on cartesian coordinates'
            CASE("REDUCED")
                    reduced=.TRUE.
                    WRITE(6,'(a)') 'Constraint is applied on reduced coordinates'
            CASE DEFAULT
                    WRITE(6,'(3a)') 'read ', coordinate, " at the end of the constraint file 'contrainte'"
                    WRITE(6,'(a)')  '  allowed values: CARTESIAN or REDUCED'
                    STOP '< InitContr >'
            END SELECT
    END IF

    IF (reduced.AND.(dmtype.NE.30)) THEN
            WRITE(0,'(a)') "You have to USE &
                &minimization on reduced coordinate when applying &
                & a constraint on reduced coordinates"
            WRITE(0,'(a)') "Use dmtype=30"
            STOP '< InitContr >'
    END IF
    IF ((.NOT.reduced).AND.(dmtype.EQ.30)) THEN
            WRITE(0,'(a)') 'You cannot USE minimization on &
                &reduced coordinates when applying a constraint on &
                &cartesian coordinate'
            WRITE(0,'(a)') "Do not USE dmtype=30"
            STOP '< InitContr >'
    END IF


    ! Différence de déplacement du centre de masse
    sigg0(:) = 0.d0
    masstot = 0.d0
    do i=1,im
       sigg0(1:3) = sigg0(1:3) + dg0(1:3,i)*cm(ityp(i))
       masstot = masstot + cm(ityp(i))
    end do
    sigg0(:)=sigg0(:)/masstot

    ! On retranche de la contrainte le déplacement du centre de masse
    ! et calcul de la norm
    do i=1,im
       dg0(1:3,i) = dg0(1:3,i) -  sigg0(1:3)
    end do

    ! Calcul de la norme de la contrainte
    dg2 = Sum( dg0(:,1:im)**2 )

    ! Normalisation de la contrainte
    dg0(:,1:im) = dg0(:,1:im)/sqrt(dg2)

    ! Calcul de la valeur de la contrainte
    IF (reduced) THEN   ! Constraint is applied on reduced coordinates
            trans_at(:,:) = Transpose( at(:,:) )        ! Transposée des vecteurs de bases
            inv_at = Transpose( bg(:,:) )               ! Inverse
            trans_inv_at(:,:) = bg(:,:)                 ! Transposée de l'inverse
            xc(:,1:im) = MatMul( inv_at(:,:), xp(:,1:im) )
            contrval = Sum( xc(:,1:im)*dg0(:,1:im) )        
    ELSE                ! Constraint is applied on cartesian coordinates
            contrval = Sum( xp(:,1:im)*dg0(:,1:im) )
    END IF

    write(6,*)'contrval=',contrval

    vp=0.
    xpp=xp

    return
  end subroutine initcontr

  ! **************************************************************
  subroutine contr(xp, vp, fp,ityp )
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY : im, imm
    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    real(double)  :: xp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: fp(3,imm)
    integer :: ityp(imm)
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    real(double) :: lambda, contrval
    real(double)  :: xc(3,imm)    
    real(double)  :: fc(3,imm)    
    !-----------------------------------------------

    !-----------------------------------------------
    IF (reduced) THEN   ! Constraint is applied on reduced coordinates

            ! Transforme les positions des atomes et les forces en coordonnées réduites
            xc(:,1:im) = MatMul( inv_at(:,:), xp(:,1:im) )
            fc(:,1:im) = MatMul( trans_at(:,:), fp(:,1:im) )

            ! Valeur de la contrainte
            contrval = Sum( xc(:,1:im)*dg0(:,1:im) )
    
            ! Projection de la force dans la direction de la contrainte
            lambda = Sum( fc(:,1:im)*dg0(:,1:im) )

            ! Force perpendiculaire à la direction de la contrainte
            fc(:,1:im) = fc(:,1:im) - lambda*dg0(:,1:im)

            ! Transforme la force en coordonnées cartésiennes
            fp(:,1:im) = MatMul( trans_inv_at(:,:), fc(:,1:im) )
    
    !-----------------------------------------------
    ELSE                ! Constraint is applied on cartesian coordinates

            ! Valeur de la contrainte
            contrval = Sum( xp(:,1:im)*dg0(:,1:im) )
    
            ! Projection de la force dans la direction de la contrainte
            lambda = Sum( fp(:,1:im)*dg0(:,1:im) )

            ! Force perpendiculaire à la direction de la contrainte
            fp(:,1:im) = fp(:,1:im) - lambda*dg0(:,1:im)

    END IF

    write(6,*)'contrval=',contrval

    return

  end subroutine contr

end module contrainte
