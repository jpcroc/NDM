MODULE gin_module
  USE gen_com_m, ONLY : imm
CONTAINS
  SUBROUTINE ReadGin(xp, iTyp, im, at, inp)


    USE T_kind_param_m, ONLY:  double
    IMPLICIT NONE
    REAL(double), dimension(1:3,1:imm), intent(out) ::xp
    INTEGER, dimension(1:imm), intent(out) :: iTyp
    INTEGER, intent(out) :: im
    REAL(double), dimension(1:3,1:3), intent(out) :: at
    INTEGER, intent(in) :: inp

    INTEGER :: i, j,k,l, j0, n
    ! Number of atoms in unit cell
    INTEGER :: imCell, xLat, yLat, zLat
    REAL(kind(0.d0)), dimension(:,:), allocatable :: xc

    ! Number of cells in 3 directions
    READ(inp,*) xLat, yLat, zLat
    IF ( (xLat.LT.1).OR.(yLat.LT.1).OR.(zLat.LT.1) ) THEN
       WRITE(0,'(a)') 'Problem with number of times the unit cell have to be duplicated'
       WRITE(0,'(a,3(i0,1x),a)') ' lat(1:3) = ', xLat, yLat, zLat, ' (read in gin file)'
       STOP '< ReadGin >'
    END IF

    ! Lattice vector coordinates of the cell
    READ(inp,*) at(1:3,1)
    READ(inp,*) at(1:3,2)
    READ(inp,*) at(1:3,3)

    ! Lattice vectors in NDM units (copyright: JPC)
    at(:,:) = 1.e-8*at(:,:)

    ! Number of atoms in the cell
    READ(inp,*) imcell

    ! Number of atoms in the simulation box
    im=imcell*xLat*yLat*zLat

    ! Maximal number of atoms in simulation box (USEd to allocate tables)
    ! Maximal number of atoms in simulation box
    IF ( (im.GT.size(xp,2)).OR.(im.GT.size(iTyp,1)) ) THEN
       WRITE(0,'(a,i0)') '  dimension of xp(1:3,:): ', size(xp,2)
       WRITE(0,'(a,i0)') '  dimension of iTyp(:): ', size(iTyp,1)
       WRITE(0,'(a,i0)') '  number of atoms read in GIN file: ', im
       STOP '< ReadGin >'
    END IF

    IF (Allocated(xc)) Deallocate(xc)
    ALLOCATE(xc(1:3,1:imcell))

    DO i=1, imcell
       READ(inp,*) xc(1:3,i), ityp(i)
    END DO

    ! Atom real coordinates
    xp(1:3,1:imcell) = MatMul( at(1:3,1:3), xc(1:3,1:imcell) )
    Deallocate(xc)

    ! Duplicate unit cell
    DO l=1, zLat
       DO k=1, yLat
          IF (k*l.EQ.1) THEN
             j0=2
          ELSE
             j0=1
          END IF
          DO j=j0, xLat
             DO i=1, imcell
                n = i + imcell*( (l-1)*yLat*xLat + (k-1)*xLat +j-1 )
                xp(1:3,n ) = xp(1:3,i) &
                     + dble(j-1)*at(1:3,1) + dble(k-1)*at(1:3,2) + dble(l-1)*at(1:3,3)
                ityp(n) = ityp(i)
             END DO
          END DO
       END DO
    END DO

    ! Lattice vector coordinates of the simulation box
    at(1:3,1) = dble(xLat)*at(1:3,1)
    at(1:3,2) = dble(yLat)*at(1:3,2)
    at(1:3,3) = dble(zLat)*at(1:3,3)



  END SUBROUTINE ReadGin

END MODULE gin_module
