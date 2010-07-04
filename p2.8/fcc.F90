MODULE fcc_module


CONTAINS

  SUBROUTINE Count_Fcc_Neighbours(xp, iwmax, alat, fcc_nVoisins, fcc_cluster, out)
    ! For each atom i, count number of neighbours in fcc position
    !    nVoisins(i)=12 if atom i has 12 first nearest neighbours in fcc position
    !               =9  if atom i has 12 first nearest neighbours, 9 of them are
    !               in fcc position (hcp structure)
    !               =-1 otherwise
    ! alat is the lattice parameter of the fcc unit cell
    ! Defect atoms, ie atoms which do not have fcc stacking, are gathered into
    ! clusters and fcc_cluster(i) is the index of the corresponding cluster
    ! This subroutine assumes that the axes are x = [110]/2, y = [-112]/2, z = [1-11]
    USE gen_com_m
    IMPLICIT NONE

    REAL(double)  :: xp(3,imm)
    integer, intent(in)  :: iwmax(imm)
    REAL(double), intent(in) :: alat
    INTEGER, intent(in), optional :: out
    INTEGER, dimension(:), intent(out) :: fcc_nVoisins, fcc_cluster

    REAL(double) :: R1, Rmax2, Rmatch2, dxp2
    INTEGER :: i1, i2, i, n, nMin, nMax, nv, iter, iw, iw1, iw2
    REAL(double), dimension(1:3) :: xp1, dxc, dxp
    ! First nearest-neighbours coordinates for FCC lattice
    REAL(double), dimension(1:3,1:12) :: voisin 
    INTEGER, dimension(:), allocatable :: nVoisins
    INTEGER, parameter :: max_nVoisins=14
    INTEGER, dimension(:,:), allocatable :: iVoisins

    ! Variables used for clusters
    !  icl(i): label of the cluster which atom i belongs to
    !  ncl(i): number of atoms in cluster i
    INTEGER, dimension(:), allocatable :: icl, ncl
    INTEGER, dimension(1:max_nVoisins) :: ldif
    INTEGER :: newlab, maxlab, ndif, somme, ipr, nz
    INTEGER, parameter :: threshold_size=10     ! Minimal number of atoms in a cluster

    INTEGER :: MaxSize  ! Cluster maximal size
    INTEGER :: ncluster      ! Cluster size
    INTEGER, dimension(1:1) :: ncluster_vec      ! Cluster size
    INTEGER :: nSize
    ! cluster_size(i): number of atoms in cluster i
    INTEGER, dimension(:), allocatable :: cluster_size
    ! cluster_index(:,i): Indexes of atoms belonging to cluster i
    INTEGER, dimension(:,:), allocatable :: cluster_index
    ! Cluster mean position and center to calculate mean position
    REAL(double), dimension(1:3) :: Rcenter, Rmean


    !-------------------------------------------------------------
    ! Allocation and initialization
    R1=alat*0.5d0*sqrt(2.d0)    ! First nearest neighbour distance
    Rmax2 = ( 0.5d0*(R1+alat) )**2      ! Average of the first and second nearest neighbour distances
    ! Half distance between edges and center of 1st nearest-neighbour triangles
    Rmatch2 = ( 0.5d0*alat/sqrt(6.d0) )**2
    IF (Allocated(nVoisins)) Deallocate(nVoisins)
    Allocate(nVoisins(1:im))
    IF (Allocated(iVoisins)) Deallocate(iVoisins)
    Allocate(iVoisins(1:max_nVoisins,1:im))
    nVoisins(:) = 0
    fcc_nVoisins(:) = 0
    iVoisins(:,:) = 0

    !----------------------------------------------------------------
    ! First nearest neighbours postion for the axes: x = [110]/2, y = [-112]/2, z = [1-11]
    voisin(1:3,1)  = alat*(/  0.d0,  1.d0/sqrt(6.d0), -1.d0/sqrt(3.d0)/) 
    voisin(1:3,2)  = alat*(/  0.5d0/sqrt(2.d0), -0.5d0/sqrt(6.d0), -1.d0/sqrt(3.d0)/) 
    voisin(1:3,3)  = alat*(/  0.5d0/sqrt(2.d0),  0.5d0*sqrt(1.5d0),  0.d0/)
    voisin(1:3,4) = alat*(/  0.5d0/sqrt(2.d0), 0.5d0/sqrt(6.d0), 1.d0/sqrt(3.d0)/)
    voisin(1:3,5)  = alat*(/  1.d0/sqrt(2.d0), 0.d0,  0.d0/) 
    voisin(1:3,6)  = alat*(/  0.5d0/sqrt(2.d0), -0.5d0*sqrt(1.5d0),  0.d0/)
    voisin(1:3,7)  = alat*(/ -0.5d0/sqrt(2.d0),  0.5d0*sqrt(1.5d0),  0.d0/)
    voisin(1:3,8)  = alat*(/ -1.d0/sqrt(2.d0), 0.d0,  0.d0/)
    voisin(1:3,9)  = alat*(/ -0.5d0/sqrt(2.d0), -0.5d0/sqrt(6.d0), -1.d0/sqrt(3.d0)/) 
    voisin(1:3,10)  = alat*(/ -0.5d0/sqrt(2.d0), -0.5d0*sqrt(1.5d0),  0.d0/)
    voisin(1:3,11) = alat*(/ -0.5d0/sqrt(2.d0), 0.5d0/sqrt(6.d0), 1.d0/sqrt(3.d0)/)
    voisin(1:3,12) = alat*(/  0.d0,  -1.d0/sqrt(6.d0), 1.d0/sqrt(3.d0)/)


    !----------------------------------------------------------------
    ! Transform cartesians coordinates to reduced ones
    CALL Cryst_to_cart(imm, xp, bg, -1)


    !--------------------------------------------------------------------------
    ! Count number of 1st nearest neighbours for each atom i -> nVoisins(i)
    ! Keep index of each 1st nearest neighbours -> iVoisins(:,i)
    ! Count number of 1st nearest neighbours in fcc position -> fcc_nVoisins(i)
    iw2=0
    DO i1=1, im ! Loop on all atoms
       xp1(1:3) = xp(1:3,i1)      ! Reduced coordinates of atom i1
       iw1=iw2+1
       iw2=iwmax(i1)
       DO iw=iw1, iw2  ! Loop on all neighbours
          i2=indi(iw)
          IF (i2.LE.i1) Cycle ! Atom pairs has already been considered
          dxc(1:3) = xp1(1:3) - xp(1:3,i2)       ! Reduced coordinates 
          WHERE ((dxc(:).GT.0.5d0).OR.(dxc(:).LT.-0.5d0))
             dxc(1:3) = dxc(1:3) - ANInt(dxc(1:3))  ! Peridodic boundary conditions
          END WHERE
          dxp(1:3) = MatMul(at(1:3,1:3),dxc(1:3))! Real coordinates
          dxp2=Sum(dxp(1:3)**2)  ! Square of the distance between atoms i1 and i2
          IF ( dxp2.LE.Rmax2 ) THEN
             nVoisins(i1)=nVoisins(i1)+1
             iVoisins(nVoisins(i1),i1)=i2
             nVoisins(i2)=nVoisins(i2)+1
             iVoisins(nVoisins(i2),i2)=i1
             DO nv=1, 12  ! Loop on perfect FCC 1st nearest neighbours
                dxp2 = Sum( ( dxp(1:3)-Voisin(1:3,nV) )**2 )
                IF (dxp2.LE.Rmatch2) THEN
                   fcc_nVoisins(i1) = fcc_nVoisins(i1) + 1
                   fcc_nVoisins(i2) = fcc_nVoisins(i2) + 1
                   EXIT
                END IF
             ENDDO
          END IF
       END DO
    END DO


    !----------------------------------------------------------------------
    ! fcc_nVoisins(i) = 12 -> atom i neighbourhood agrees with fcc stacking
    ! ............... = 9  -> ................................ hcp
    ! ............... = -1 -> .................... disagrees with fcc and hcp stacking
    DO i1=1, im
       IF (nVoisins(i1).EQ.12) THEN
          IF ( (fcc_nVoisins(i1).NE.12).AND.(fcc_nVoisins(i1).NE.9) ) &
               fcc_nVoisins(i1) = -1
       ELSE
          fcc_nVoisins(i1) = -1
       END IF
    END DO

    !----------------------------------------------------------------
    ! Gather all atoms which stacking is not fcc in labelled clusters
    IF (Allocated(icl)) Deallocate(icl)
    Allocate(icl(1:im))
    IF (Allocated(ncl)) Deallocate(ncl)
    Allocate(ncl(1:im))
    icl(:)=0 ; ncl(:)=0
    ldif(:)=0
    maxlab=1
    DO i1=1, im ! Loop on all atoms
       IF (fcc_nVoisins(i1).EQ.12) Cycle
       newlab=maxlab
       ndif=0
       somme=0
       DO nz=1, nVoisins(i1)   ! Loop on all first nearest neighbours
          i2=iVoisins(nz,i1)
          IF (fcc_nVoisins(i2).EQ.12) Cycle
          ipr=icl(i2) ! Cluster which atom i belongs to
          IF (ipr.EQ.0) Cycle
          if(ncl(ipr).GT.0) goto 103
          ipr=-ncl(ipr)
          if(ncl(ipr).GT.0) goto 103
101       ipr=-ncl(ipr)
          if(ncl(ipr).GT.0) goto 102
          goto 101
102       ncl(icl(i2))=-ipr
103       newlab=min(newlab,ipr)
          do i=1,ndif
             if (ipr.EQ.ldif(i)) goto 104
          end do
          ndif=ndif+1
          ldif(ndif)=ipr
          somme=somme+ncl(ipr)
104       Continue
       END DO
       if(ndif.ne.0) then
          ! on fait pointer sur newlab tous les amas qu'il faut rattacher
          do i=1,ndif
             ncl(ldif(i))=-newlab
          end do
          ncl(newlab)=somme+1
          icl(i1)=newlab
       else
          ncl(newlab)=1
          icl(i1)=maxlab
          maxlab=maxlab+1
       end if
    END DO

    IF (Present(out)) THEN
       WRITE(out,*)
       WRITE(out,'(a)') 'COUNT 1st NEAREST NEIGHBOURS'
       WRITE(out,'((a,f0.5,a))') ' Rmax = ', sqrt(Rmax2)*1e8, ' A'
       WRITE(out,'((a,f0.5,a))') ' alat = ', alat*1e8, ' A'
       nMin = MinVal(nVoisins(1:im)) ; nMax = MaxVal(nVoisins(1:im))
       DO n=nMin, nMax
          i=Count(nVoisins(1:im).EQ.n)
          WRITE(out,'(1x,2(i0,a),f0.3,a)') i, ' atoms with ', &
               n, ' neighbours  (', 100.d0*dble(i)/dble(im) , '%)'
       END DO
    END IF

    DEALLOCATE(iVoisins,nVoisins)

    !------------------------------------------------------------------------
    ! Eliminate all atoms belonging to a cluster containing less than a critcal
    ! size and build tables:
    !   cluster_size(i): number of atoms in cluster i
    !   cluster_index(:,i): Indexes of atoms belonging to cluster i
    !   fcc_cluster(ia): index of cluster which atom i belongs to
    MaxSize=MaxVal(ncl(1:Maxlab))          ! Cluster maximal size
    IF (Allocated(cluster_size)) Deallocate(cluster_size)
    ALLOCATE(cluster_size(1:maxlab))    
    cluster_size(:)=0
    IF (Allocated(cluster_index)) Deallocate(cluster_index)
    ALLOCATE(cluster_index(1:maxsize,1:maxlab))   
    cluster_index(:,:)=0
    fcc_cluster(:)=0

    DO i1=1, im
       IF (fcc_nVoisins(i1).EQ.12) Cycle
       ipr=icl(i1)     ! Cluster index
       DO WHILE (ncl(ipr).LT.0)
          ipr=-ncl(ipr)
       END DO
       IF (ncl(ipr).LT.threshold_size) THEN
          fcc_nVoisins(i1)=12
          ncl(ipr)=0
       ELSE
          fcc_cluster(i1) = ipr
          cluster_size(ipr) = cluster_size(ipr)+1
          cluster_index(cluster_size(ipr),ipr) = i1
       END IF
    END DO

    DEALLOCATE(icl, ncl)

    !----------------------------------------------------------------
    IF (Present(out)) THEN
       WRITE(out,*)
       WRITE(out,'(a,i0,a)') 'DEFECT CLUSTERS (threshold size = ', threshold_size, ')'
       i=Count(fcc_nVoisins(1:im).EQ.12)
       WRITE(out,'(1x,i0,a,f0.3,a)') i, ' atoms with FCC neighbourhood &
            & (', 100.d0*dble(i)/dble(im) , '%)'
       i=Count(fcc_nVoisins(1:im).EQ.9)
       WRITE(out,'(1x,i0,a,f0.3,a)') i, ' atoms with HCP neighbourhood &
            & (', 100.d0*dble(i)/dble(im) , '%)'
       i=Count(fcc_nVoisins(1:im).EQ.-1)
       WRITE(out,'(1x,i0,a,f0.3,a)') i, ' atoms with neither FCC nor HCP neighbourhood &
            & (', 100.d0*dble(i)/dble(im) , '%)'

       !----------------------------------------------------------------
       ! Print information about the clusters (size, mean position)
       iter=0
       WRITE(out,*)
       DO     ! Loop on all clusters starting from the largest one
          iter=iter+1
          nCluster_vec(:)=MaxLoc(cluster_size)    ! Bigger cluster
          nCluster=nCluster_vec(1)
          nSize=cluster_size(nCluster)     ! Number of atoms in cluster
          IF (nSize.LT.threshold_size) EXIT     ! No more cluster to consider
          cluster_size(nCluster)=0

          ! Calculate mean position of the cluster
          i1=cluster_index(1,ncluster)     ! Random atom belonging to the cluster
          Rmean(1:3) = xp(1:3,i1)          ! This atom considered as the center of the cluster
          DO i=1, 2
             Rcenter(1:3) = Rmean(1:3)    ! Center to calculate mean position
             Rmean(1:3) = 0.d0
             DO n=1, nSize   ! Loop on all atoms belonging to cluster ipr
                i1=cluster_index(n,ncluster)      
                dxc(1:3) = xp(1:3,i1) - Rcenter(1:3) 
                WHERE ((dxc(:).GT.0.5d0).OR.(dxc(:).LT.-0.5d0))
                   dxc(1:3) = dxc(1:3) - ANInt(dxc(1:3))  ! Peridodic boundary conditions
                END WHERE
                Rmean(1:3) = Rmean(1:3) + dxc(1:3)
             END DO
             Rmean(:) = Rcenter(:) + Rmean(:)/dble(nSize) 
             WHERE ((Rmean(:).GT.0.5d0).OR.(Rmean(:).LT.-0.5d0))
                Rmean(1:3) = Rmean(1:3) - ANInt(Rmean(1:3))  ! Peridodic boundary conditions
             END WHERE
          END DO

          ! Output
          WRITE(out,'(3(a,i0),a, 3(f0.3,1x))') &
               ' cluster ', iter, ' labelled ', nCluster, ' containing ', nSize , &
               ' atoms located in X = ', 1.e8*MatMul(at(:,:), Rmean(1:3))
       END DO
    END IF
    !----------------------------------------------------------------

    DEALLOCATE(cluster_index, cluster_size)

    ! Transform back reduced coordinates to cartesian ones
    Call cryst_to_cart (imm, xp, at, 1)     

  END SUBROUTINE Count_Fcc_Neighbours

END MODULE fcc_module
