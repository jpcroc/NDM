module decoupage_mod
  USE arret_ndm_mod,only: arret_ndm
      USE T_kind_param_m, ONLY:  double
  implicit none
contains
  subroutine decoupage(nbr_cpuIN,ncore)

#ifdef PARA
    USE mpi
    USE mod_para,only:status,ierr,nprocs,myid,NDM_MPI_REAl_DOUBLE,res_cpu,coord_max,coord_min,proc_cell
    use tab_imm_m,only:realloc_all_tab_imm
#endif
    USE gen_com_m, ONLY:nox,noy,noz,rang,cell_debx,cell_deby,cell_debz,cell_finx,cell_finy,cell_finz,imm_glob,&
         &nb_cell_x,nb_cell_y,nb_cell_z,noxyz,imm

    USE var_pot, ONLY:

    implicit none

    !------------------
    integer :: nbr_cpuIN !Egal aussi au nombre de zone qu'on d�coupera dans la boite
    integer::ncore ! nb de coeur par noeud
    integer:: nnoeuds
    integer :: nb_sol  !nbr de decoupage possible (n+1)(n+2)/2
    integer :: num_sol !iteration du decoupage possible
    integer :: test

    integer, allocatable :: decoup(:,:) !tableau comprenant l'ensemble des decoupages 
    !possibles en fct des 3 dimensions
    real(double), allocatable :: specifs(:,:) 

#ifndef PARA
    integer, allocatable :: coord_min(:,:),coord_max(:,:),res_cpu(:,:)
#endif

    integer :: messages,messages_max,messages_min
    integer :: tailleminx,tailleminy,tailleminz   ! nbre de cellule min
    ! dans une direction
    integer :: taillemaxx,taillemaxy,taillemaxz   ! nbre de cellule max
    ! dans une direction
    integer :: restex,restey,restez
    integer :: interne_min,interne_max   ! nbre de cellules internes min et max
    integer :: local_min,local_max       ! nbre de cellules locales min et max
    ! (local = interne+ frontiere)
    integer :: fantome_min,fantome_max   ! nbre de cellules fantomes min et max
    real(double) :: equilibrage
    integer :: solution
    integer :: num_cpu
    integer :: kx,ky,kz,koo
    integer :: cellules_max
    integer :: imm_loc,nbr_cpu

    !------------------

    integer :: ii,jj,kk,ldecoup !indice de boucle


    !------------------
    !Allocation des tableaux



#ifdef DECOUP
    loop1:     do nbr_cpu=2,nbr_cpuIN

#else
       nbr_cpu=nbr_cpuIN
#endif



       nb_sol = 0

       if (rang==0) then

          print *
          print *
          print *,'-----------------------------------------------------------'
          print *,'Procedure de decoupage pour une boite de taille  :',nox,noy,noz
          print *,'sur ',nbr_cpu,' cpus',ncore
       endif

       do ii=1,nbr_cpu
          do jj=1,nbr_cpu
             do kk=1,nbr_cpu
                test=kk * jj * ii
                if ( test == nbr_cpu ) then
                   ! C'est une solution possible
                   ! on ne la conserve que si elle est plausible par rapport a la geometrie
                   if ( ii <= nox .and. jj <= noy .and. kk <= noz ) then
                      nb_sol = nb_sol + 1  
                   endif
                endif
             enddo
          enddo
       enddo
       allocate(decoup(nb_sol,3))
       nb_sol = 0
       do ii=1,nbr_cpu
          do jj=1,nbr_cpu
             do kk=1,nbr_cpu
                test=kk * jj * ii
                if ( test == nbr_cpu ) then
                   ! C'est une solution possible
                   ! on ne la conserve que si elle est plausible par rapport a la geometrie
                   if ( ii <= nox .and. jj <= noy .and. kk <= noz ) then
                      nb_sol = nb_sol + 1  
                      decoup(nb_sol,1)=ii
                      decoup(nb_sol,2)=jj
                      decoup(nb_sol,3)=kk
                   endif
                endif
             enddo
          enddo
       enddo
       allocate(specifs(nb_sol,7))

       allocate(res_cpu(0:nbr_cpu-1,3))
       allocate(coord_min(0:nbr_cpu-1,3))
       allocate(coord_max(0:nbr_cpu-1,3))

#ifdef PARA
       allocate(proc_cell(nox*noy*noz))
#endif

       !print *,'Voici les cas possibles :'
       !print *,'-------------------------'
       do ii=1,nb_sol
          !   print *,decoup(ii,:)
       enddo
       !print *

       if (nb_sol==0) then
          if (rang==0) then
             print *,'!!! Pas de possibilite de decoupage pour la configuration demandee !!!'
             print *,'!!! nx / ny / nz / nb_cpu :',nox,noy,noz,nbr_cpu
             print *,'!!! Arret du programme !!!'
          endif
#ifndef DECOUP
          call arret_ndm
#else
          deallocate(decoup)
          deallocate(specifs)

          deallocate(res_cpu)
          deallocate(coord_min)
          deallocate(coord_max)

          cycle loop1
#endif

       endif

       if (rang==0) print *,'Nbre de solutions possibles : ',nb_sol

       ! On scanne l'ensemble des solutions proposees pour en calculer 
       ! l'equilibrage de charge et le nombre de cellules fantomes
       messages_max = 0
       messages_min = nox * noy * noz
       do num_sol = 1,nb_sol
          tailleminx = int(nox/decoup(num_sol,1))
          tailleminy = int(noy/decoup(num_sol,2))
          tailleminz = int(noz/decoup(num_sol,3))
          taillemaxx = tailleminx
          if (tailleminx*decoup(num_sol,1) .ne. nox) taillemaxx = taillemaxx + 1
          taillemaxy = tailleminy
          if (tailleminy*decoup(num_sol,2) .ne. noy) taillemaxy = taillemaxy + 1
          taillemaxz = tailleminz
          if (tailleminz*decoup(num_sol,3) .ne. noz) taillemaxz = taillemaxz + 1
          interne_min = max(0,tailleminx-2) * max(0,tailleminy-2) * max(0,tailleminz-2)
          interne_max = max(0,taillemaxx-2) * max(0,taillemaxy-2) * max(0,taillemaxz-2)
          local_min = tailleminx * tailleminy * tailleminz
          local_max = taillemaxx * taillemaxy * taillemaxz
          fantome_min = (tailleminx+2)*(tailleminy+2)*(tailleminz+2)-local_min
          fantome_max = (taillemaxx+2)*(taillemaxy+2)*(taillemaxz+2)-local_max
          ! Les criteres retenus pour le decoupage sont :
          ! - l'equilibrage de charge
          ! - le minimum d'echange de messages
          equilibrage = real(local_min) / real(local_max)
          messages = fantome_max
          messages_min = min(messages,messages_min)
          messages_max = max(messages,messages_max)
          specifs(num_sol,1) = equilibrage
          specifs(num_sol,2) = messages
          specifs(num_sol,4) = local_min
          specifs(num_sol,5) = local_max
          specifs(num_sol,6) = fantome_min
          specifs(num_sol,7) = fantome_max

       enddo
       solution = 1
       do num_sol = 1,nb_sol
          ! La formule magique ! (a voir si il faut modifier les coefficients)
          specifs(num_sol,3) = 1.d-3*(specifs(num_sol,1)*1000+100*(-1+specifs(num_sol,4)/specifs(num_sol,7)))
          !     specifs(num_sol,3) = 1.d-3*(specifs(num_sol,1)*1000+2000/specifs(num_sol,2)+500*(-1+specifs(num_sol,4)/specifs(num_sol,7)))

          if (ncore.ne.0) then
             if (mod(nbr_cpu,ncore).ne.0)then
                nnoeuds=nbr_cpu/ncore+1
             else
                nnoeuds=nbr_cpu/ncore
             end if
             write(6,*)'GOGO',specifs(num_sol,3),ncore,nbr_cpu,nnoeuds
             specifs(num_sol,3)=specifs(num_sol,3)*nbr_cpu/(ncore*nnoeuds)
          end if
          if (specifs(num_sol,3)>specifs(solution,3)) solution = num_sol
       enddo

       if (rang==0) write(6,'(A,4I5,F10.4)')'LE MEILLEUR DECOUPAGE :',nbr_cpu, decoup(solution,1), & 
            decoup(solution,2), decoup(solution,3),specifs(solution,3)

       !Calcul des xmin, ymin, zmin pour chaque decoupage
       tailleminx = int(nox/decoup(solution,1))
       tailleminy = int(noy/decoup(solution,2))
       tailleminz = int(noz/decoup(solution,3))
       restex = nox - tailleminx*decoup(solution,1)
       restey = noy - tailleminy*decoup(solution,2)
       restez = noz - tailleminz*decoup(solution,3)

       do ii=1,decoup(solution,1)
          do jj=1,decoup(solution,2)
             do kk=1,decoup(solution,3)
                num_cpu = ii-1 + (jj-1)*decoup(solution,1) + (kk-1) * decoup(solution,1)*decoup(solution,2)

                coord_min(num_cpu,1) = (ii-1)*tailleminx +1
                if (ii>1) coord_min(num_cpu,1) = coord_min(num_cpu,1) + min(ii-1,restex)
                coord_max(num_cpu,1) = coord_min(num_cpu,1) + tailleminx -1
                if (ii<=restex) coord_max(num_cpu,1) = coord_max(num_cpu,1) + 1

                coord_min(num_cpu,2) = (jj-1)*tailleminy +1
                if (jj>1) coord_min(num_cpu,2) = coord_min(num_cpu,2) + min(jj-1,restey)
                coord_max(num_cpu,2) = coord_min(num_cpu,2) + tailleminy -1
                if (jj<=restey) coord_max(num_cpu,2) = coord_max(num_cpu,2) + 1

                coord_min(num_cpu,3) = (kk-1)*tailleminz +1
                if (kk>1) coord_min(num_cpu,3) = coord_min(num_cpu,3) + min(kk-1,restez)
                coord_max(num_cpu,3) = coord_min(num_cpu,3) + tailleminz -1
                if (kk<=restez) coord_max(num_cpu,3) = coord_max(num_cpu,3) + 1

                res_cpu(num_cpu,1) = coord_max(num_cpu,1) - coord_min(num_cpu,1) + 1
                res_cpu(num_cpu,2) = coord_max(num_cpu,2) - coord_min(num_cpu,2) + 1
                res_cpu(num_cpu,3) = coord_max(num_cpu,3) - coord_min(num_cpu,3) + 1

#ifdef PARA
                ! On affecte ces cellules au processeur concerne
                do kx = coord_min(num_cpu,1), coord_max(num_cpu,1)
                   do ky = coord_min(num_cpu,2), coord_max(num_cpu,2)
                      do kz = coord_min(num_cpu,3), coord_max(num_cpu,3)
                         koo = 1+(kx-1)+nox*((ky-1)+noy*(kz-1))
                         proc_cell(koo) = num_cpu
                      enddo
                   enddo
                enddo
#endif

             enddo
          enddo
       enddo


#ifdef PARA
       if (myid == 0) then
#endif
          Print *,'-----------------------------------------------'
          print *,'          FIN DU CALCUL DU DECOUPAGE :         '
          print *
          print *,'Nbre de cellule suivant x :',nox
          print *,'Nbre de cellule suivant y :',noy
          print *,'Nbre de cellule suivant z :',noz
          print *
          print *,'Le cas choisit :',decoup(solution,1),decoup(solution,2),decoup(solution,3)
          print *,'Nombre de cellules max echg  par proc :',nbr_cpu,int(specifs(solution,2))
          print *,'Taux d''equilibrage :', nbr_cpu,specifs(solution,1) 
          print *,'valeur ',nbr_cpu, specifs(solution,3)

          !     print *,'Nbre Min/Max de cel. echangees dans les conf. :', messages_min , messages_max
          print *,'Nbre Min/Max de cel. locales. :', int(specifs(solution,4)),int(specifs(solution,5))
          print *,'Nbre Min/Max de cel. fantomes. :',int(specifs(solution,6)),int(specifs(solution,7))
          print *,'locmin/fantomax :',nbr_cpu,specifs(solution,4)/specifs(solution,7)



#ifndef DECOUP
          ldecoup=1023
          open (unit=1023,file='decoup_out')
#else
          ldecoup=6
#endif
#ifdef PARA
#else
          write(ldecoup,*)'Taille des decoupages'
          do ii=0,nbr_cpu-1
             write(ldecoup,*)'Decoupage',ii,':',res_cpu(ii,1:3)
          enddo
          write(ldecoup,*)'----------------------------------------------'
          do ii = 0,nbr_cpu-1   
             write(ldecoup,*)'Debut/Fin en x pour ii',ii,'egal',coord_min(ii,1),coord_max(ii,1)
             write(ldecoup,*)'Debut/Fin en y pour ii',ii,'egal',coord_min(ii,2),coord_max(ii,2)
             write(ldecoup,*)'Debut/Fin en z pour ii',ii,'egal',coord_min(ii,3),coord_max(ii,3)
             write(ldecoup,*)
          enddo
          write(ldecoup,*)'-----------------------------------------------'

#endif
#ifdef PARA
       endif
#endif


#ifndef DECOUP
#ifdef PARA
       ! On est dans le code de calcul NDM, on realloue les tableaux sur le
       ! nombre d'atomes en tenant compte des cellules fantomes


       cellules_max=0
       do ii = 0,nbr_cpu-1
          cellules_max = max(cellules_max,(res_cpu(ii,1)+2) * (res_cpu(ii,2)+2)* (res_cpu(ii,3)+2))
       enddo
       cellules_max = min (cellules_max, noxyz)
       imm_loc = min( imm_glob, int(1.2 * imm_glob * cellules_max / noxyz) )
       ! Le processeur maitre recupere la valeur maximale des imm des
       ! differents processeurs afin de pouvoir receptionner les tableaux
       ! des autres processeurs lors d'I/O :
       imm = imm_loc
       call MPI_REDUCE(imm_loc,imm,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)

       !     print *,'test4' 
       call realloc_all_tab_imm(imm)

       !     print *,'test4' 
       ! Initialisation des donnees geometriques qui serviront pour le reste du code :
       cell_debx= coord_min(myid,1)
       cell_finx= coord_max(myid,1)
       cell_deby= coord_min(myid,2)
       cell_finy= coord_max(myid,2)
       cell_debz= coord_min(myid,3)
       cell_finz= coord_max(myid,3)
       nb_cell_x= cell_finx - cell_debx + 1
       nb_cell_y= cell_finy - cell_deby + 1
       nb_cell_z= cell_finz - cell_debz + 1
#endif
#endif

#ifdef DECOUP
       deallocate(decoup)
       deallocate(specifs)

       deallocate(res_cpu)
       deallocate(coord_min)
       deallocate(coord_max)

#ifdef PARA
       deallocate(proc_cell)
#endif

    enddo loop1



#else
    !        nbr_cpul=nbr_cpuIN
#endif


  end subroutine decoupage
end module decoupage_mod
