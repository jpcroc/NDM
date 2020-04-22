module tab_imm_m
  !
  USE T_kind_param_m
  USE gen_com_m, ONLY: lsuivinonpbc,lposmoy,mdcg_noise,llangevin,l2T
  !$ USE OMP_LIB
  ! 
  ! Module contenant les tableaux dimmensionnes sur le
  ! nombre d'atomes de la simulation

  integer, dimension(:), allocatable       :: ielat  ! numero de cellule
  integer, dimension(:), allocatable       :: iwmax  ! indice du dernier voisin
  integer, dimension(:), allocatable       :: iwmax2  ! indice du dernier voisin pour les constantes de force (only phondy)
  integer, dimension(:), allocatable       :: ityp   ! types 
  integer, dimension(:), allocatable       :: ityp_buffer   ! temp/iorary store the types buffer when we
                                                        ! perform temporary changes on ityp 
  real(double),dimension(:,:), allocatable :: xp     ! positions 
  real(double),dimension(:,:), allocatable :: xpp    ! positions precedentes
  real(double),dimension(:,:), allocatable :: vp     ! vitesses
  real(double),dimension(:,:), allocatable :: ax     ! positions d'origine
  real(double),dimension(:,:), allocatable :: posmoyx     ! positions moyennes
  real(double),dimension(:,:), allocatable :: fp     ! forces 
  real(double),dimension(:,:), allocatable :: bruitmd     ! bruitmd 
  real(double),dimension(:,:), allocatable :: xpnonpbc    ! only in the case, lsuivinonpbc  
  real(double),dimension(:,:), allocatable :: axnonpbc    ! only in the case, lsuivinonpbc  
  real(double),dimension(:,:), allocatable :: tmpsuivi    ! only in the case, lsuivinonpbc  
  real(double),dimension(:,:), allocatable :: Gl    ! random noise langevin

  integer, dimension(:), allocatable       :: num_at_glob ! numero global d'un atome
contains

  !--------------------------------------------------------------------------!
  !Routine d'allocation des tableaux dimensionnes sur le nombre max d'atomes

  subroutine alloc_all_tab_imm(nb_imm)
    implicit none

    integer :: nb_imm

    allocate(xp(3,nb_imm))
    xp = 0.0
    allocate(xpp(3,nb_imm))
    xpp = 0.0
    allocate(vp(3,nb_imm))
    vp = 0.0
    allocate(ax(3,nb_imm))
    ax = 0.0
    if (lposmoy.eqv..true.) then
       allocate(posmoyx(3,nb_imm))
       posmoyx = 0.0
    end if
     if ((llangevin.eqv..true.).or.(l2T.eqv..true.)) then
        allocate(Gl(3,nb_imm))
     end if

    if (mdcg_noise/=0) then
     allocate (bruitmd(3,nb_imm))
     bruitmd=0.0
    end if 

    allocate(fp(3,nb_imm))
    fp = 0.0
    allocate(ielat(nb_imm))
    ielat = 0
    allocate(iwmax(nb_imm))
    allocate(iwmax2(nb_imm))
    iwmax = 0
    allocate(ityp(nb_imm))
    allocate(ityp_buffer(nb_imm))
    ityp  = 0
    ityp_buffer = 0
    allocate(num_at_glob(nb_imm))
    num_at_glob  = 0

    if (lsuivinonpbc) then
      allocate(xpnonpbc(3,nb_imm))
      xpnonpbc = 0.0
      allocate(axnonpbc(3,nb_imm))
      axnonpbc = 0.0
      allocate(tmpsuivi(3,nb_imm))
      tmpsuivi = 0.0
    end if  
    

  end subroutine alloc_all_tab_imm


  !--------------------------------------------------------------------------!
  !Routine de reallocation des tableaux dimensionnes sur le nombre max 
  !d'atomes
  ! Comme il n'est pas possible de redimensionner un tableau existant, on 
  ! passe par une recopie dans un buffer, une destruction et une reallocation
  ! Cette routine est utilisee specifiquement lors de la reallocation des 
  ! tableaux d'atomes apres avoir eu connaissance du nombre d'atomes fantomes
  ! previsionnel -> il s'agit donc d'une augmentation de la taille

  subroutine realloc_all_tab_imm(new_nb_imm)
    implicit none

    integer :: new_nb_imm
    integer :: old_nb_imm
    integer, dimension(:), allocatable       :: ibuff
    real(double),dimension(:,:), allocatable :: rbuff

    old_nb_imm = size(xp,2)
    if ( old_nb_imm < new_nb_imm) then

       ! Tableaux de reels

       allocate(rbuff(3,old_nb_imm))

       rbuff = xp
       deallocate(xp)
       allocate(xp(3,new_nb_imm))
       xp = 0
       xp(:,1:old_nb_imm) = rbuff

       rbuff = xpp
       deallocate(xpp)
       allocate(xpp(3,new_nb_imm))
       xpp = 0
       xpp(:,1:old_nb_imm) = rbuff

       rbuff = vp
       deallocate(vp)
       allocate(vp(3,new_nb_imm))
       vp = 0
       vp(:,1:old_nb_imm) = rbuff

       rbuff = ax
       deallocate(ax)
       allocate(ax(3,new_nb_imm))
       ax = 0
       ax(:,1:old_nb_imm) = rbuff
       if (lsuivinonpbc) then
       ! 
	rbuff = xpnonpbc
        deallocate(xpnonpbc)
        allocate(xpnonpbc(3,new_nb_imm))
        xpnonpbc = 0
        xpnonpbc(:,1:old_nb_imm) = rbuff
       ! 
	rbuff = axnonpbc
        deallocate(axnonpbc)
        allocate(axnonpbc(3,new_nb_imm))
        axnonpbc = 0
        axnonpbc(:,1:old_nb_imm) = rbuff
        !
        rbuff = tmpsuivi
        deallocate(tmpsuivi)
        allocate(tmpsuivi(3,new_nb_imm))
        tmpsuivi = 0
        tmpsuivi(:,1:old_nb_imm) = rbuff
       !
       end if
       rbuff = fp
       deallocate(fp)
       allocate(fp(3,new_nb_imm))
       fp = 0
       fp(:,1:old_nb_imm) = rbuff
     
 
       if (mdcg_noise/=0) then
        rbuff = bruitmd
        deallocate(bruitmd)
        allocate(bruitmd(3,new_nb_imm))
        bruitmd = 0
        bruitmd(:,1:old_nb_imm) = rbuff
       end if

       deallocate(rbuff)
     ! Tableaux d'entiers

       allocate(ibuff(old_nb_imm))

       ibuff = ielat
       deallocate(ielat)
       allocate(ielat(new_nb_imm))
       ielat = 0
       ielat(1:old_nb_imm) = ibuff
       ibuff = iwmax
       deallocate(iwmax)
       deallocate(iwmax2)
       allocate(iwmax(new_nb_imm))
       allocate(iwmax2(new_nb_imm))
       iwmax = 0
       iwmax(1:old_nb_imm) = ibuff

       ibuff = ityp
       ibuff = ityp_buffer
       deallocate(ityp)
       deallocate(ityp_buffer)
       
       allocate(ityp(new_nb_imm))
       allocate(ityp_buffer(new_nb_imm))
       ityp = 0
       ityp_buffer = 0
       ityp(1:old_nb_imm) = ibuff
       ityp_buffer(1:old_nb_imm) = ibuff

       ibuff = num_at_glob
       deallocate(num_at_glob)
       allocate(num_at_glob(new_nb_imm))
       num_at_glob = 0
       num_at_glob(1:old_nb_imm) = ibuff

       deallocate(ibuff)

    endif

  end subroutine realloc_all_tab_imm


  !--------------------------------------------------------------------------!
  !Routine de deallocation des tableaux dimensionnes sur le nombre max 
  ! d'atomes

  subroutine dealloc_all_tab_imm
    implicit none

    if (allocated(xp))          deallocate(xp)
    if (allocated(xpp))         deallocate(xpp)
    if (allocated(vp))          deallocate(vp)
    if (allocated(ax))          deallocate(ax)
    if (allocated(fp))          deallocate(fp)
    if (allocated(bruitmd))      deallocate(bruitmd)
    if (allocated(ielat))       deallocate(ielat)
    if (allocated(iwmax))       deallocate(iwmax)
    if (allocated(iwmax2))      deallocate(iwmax2)
    if (allocated(ityp))        deallocate(ityp)
    if (allocated(ityp_buffer)) deallocate(ityp_buffer)
    if (allocated(num_at_glob)) deallocate(num_at_glob)
    if (lsuivinonpbc) then
      if (allocated(xpnonpbc))  deallocate(xpnonpbc)
    end if 

    if (lsuivinonpbc) then 
      if (allocated(axnonpbc))  deallocate(axnonpbc)
    end if 

    if (lsuivinonpbc) then 
            if (allocated(tmpsuivi)) deallocate(tmpsuivi)
    end if 

    if (llangevin) then 
       if(allocated(Gl)) deallocate(Gl)
    end if 

  end subroutine dealloc_all_tab_imm
end module tab_imm_m
