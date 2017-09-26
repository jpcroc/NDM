!  
subroutine read_mab_file()
 USE T_kind_param_m, ONLY:  double
 use gen_com_m, ONLY: lenfnam,fnam,angst,ev2erg,im,rang,A2cm,umass,erg2ev
 use tab_imm_m, ONLY: ityp
 use var_pot,   ONLY: cm
 USE mab_in_ndm_module, ONLY: dtlang,nlangevin,temperature,KtoERG,a0bcc,deltasph,  &
                              radiussph,nhisto,deltar1,deltar2,abf_type,block,     &
                              sim_mode,rtestlac,langevin_type,gamma,omega_abf,     &
                              omega_einstein,omega_veinstein, matfor, abf_mode,            & 
                              nwrite_histo,sigma_eta,ecart_eta, &
                              eta_mab,eta_ABFee,histo_equi,n_equilibre,           & 
                              maxforce,compute_mode,error_step,nom_deconvo,lang_factor, &
                              mode_zeta_potential, alpha_zeta,ntestvacancyjump,ha_mix,  &
                              temperature_zeta_min,temperature_zeta_max,          &
                              nsite_block, isite_block,atom_to_jump,              &
                              itype_reaction,itype_einstein, units_phondy,m_i,    &
                              nimage_neb,nimage_lambda,block_file,                          &
                              abf_mode_reaction, abf_mode_alchemical, abf_mode_temperature,rangmab, &
                              abf_restart, idebug

 implicit none
 namelist /input_mab/ dtlang,nlangevin,temperature,a0bcc,deltasph,radiussph,       &
                      nhisto,deltar1,deltar2,block,block_file, abf_type,sim_mode,rtestlac,     &
                      langevin_type,omega_einstein,gamma,omega_abf, nwrite_histo,  &
                      eta_mab,eta_ABFee,histo_equi,n_equilibre,            &
                      maxforce,compute_mode,abf_mode, error_step,nom_deconvo,lang_factor, &
                      mode_zeta_potential, alpha_zeta,ntestvacancyjump,ha_mix,            &
                      temperature_zeta_min,temperature_zeta_max,                          &
                      atom_to_jump, itype_reaction,itype_einstein,nimage_neb,nimage_lambda, &
                      abf_restart, idebug

 character(len=128) :: fnamtin, fnamt_lblock, fnamt_lfreq, fnamt_lm, fnamt_lu, fnamt_lv
 integer :: lumab,lublock,lcu,lcv,lcm
 integer :: ii,ia,ja,i,j,imax,iu,iv
 real(double) :: temp_read
 integer :: itemp_read


 idebug = 0
 rtestlac=0.1d0
 nwrite_histo=1000000
 omega_abf=1.d0
 maxforce=2.d0   ! in order to enhance the max force on the protective domains 
 nom_deconvo=10
 abf_mode = 1   !vacancy jump 
 abf_type = 1   !ABF BIN 
 omega_einstein=5.d0 ! einstein frequecy in THz
 lang_factor=1.d0
 mode_zeta_potential=0
 alpha_zeta=2.d0
 ntestvacancyjump=200
 ha_mix=0.d0
 temperature_zeta_min=150.d0
 temperature_zeta_max=1800.d0
 block=.false.
 block_file=.false.
! The atom which is desinged to jump ...
 atom_to_jump=7
 itype_reaction=0 ! 0 for vacancy, 1 for NEB 0 K reaction
 itype_einstein=0 ! 0 - einstein, 1 HA, 2 Morse (for rthe future)
      !neb part:

 nimage_neb=15
 nimage_lambda=70
 abf_restart=.false.

     
     do ia=1,3
      m_i(ia,1:im) = cm(ityp(1:im))
     end do


 fnamtin = fnam(1:lenfnam)//'.mab'
 write(*,*) 'file name', fnamtin
 lumab = 778
open(unit=lumab, file=fnamtin, status='unknown')
read (lumab, nml=input_mab)



 if (langevin_type==2) then
   !if (abf_mode==2) then
   !   write(6,*) 'Alchemical transition not yet implemented with the underdamped Langevin'
   !   write(6,*) 'put langevin_type = 1 and restart'
   !   write(6,*) 'stop in <read_mab_file>'
   !   stop 
   !end if 




   select case (abf_type)
      case (1)
             if (rangmab==0)   write(6,'("MAB:    Damped Langevin dynamics  ")')
      case (2)  
            if (rangmab==0) then
             if (abf_mode==1)  write(6,'("MAB:  Damped Langevin + ABF BIN dynamics + Reaction coordinate   ")')
             if (abf_mode==2)  write(6,'("MAB:  Damped Langevin + ABF BIN dynamics + Alchemical coordinate ")')
             if (abf_mode==22) write(6,'("MAB:  Damped Langevin + ABF BIN dynamics + Temperature coordinate")')
            end if  
      case (3) 
            if (rangmab==0) then
             if (abf_mode==1)  write(6,'("MAB:  Damped Langevin + ABF BIN dynamics with Omega + Reaction coordinate   ")')
             if (abf_mode==2)  write(6,'("MAB:  Damped Langevin + ABF BIN dynamics with Omega + Alchemical coordinate ")')
             if (abf_mode==22) write(6,'("MAB:  Damped Langevin + ABF BIN dynamics with Omega + Temperature coordinate")')
            end if  
      case (4) 
             if (rangmab==0) write(6,*) ' Damped Langevin + ABF GAUSSIAN dynamics'
      case (5) 
             if (rangmab==0) write(6,*) ' Damped Langevin + ABF EE dynamics'
      case (6) 
             if (rangmab==0) write(6,*) ' Damped Langevin + ABF EE dynamics + constant biais'
      case (7) 
             if (rangmab==0) write(6,*) ' Damped Langevin + ABF BIN dynamics + constant biais'
      case (8) 
             if (rangmab==0) write(6,*) ' Damped Langevin + ABF EE dynamics + iterative'
   end select 
 end if

 if (langevin_type==1) then
   select case (abf_type)
      case (1)
             if (rangmab==0)   write(6,'("MAB:    OverDamped Langevin dynamics  ")')
      case (2)
             if (abf_mode==1) write(6,*) ' Damped Langevin + ABF BIN dynamics'
             if (abf_mode==2) write(6,*) ' Damped Langevin + ABF BIN dynamics + External parameter'
      case (3) 
             if (abf_mode==1) write(6,*) ' Overdamped Langevin + ABF BIN dynamics with Omega'
             if (abf_mode==2) write(6,*) ' Overdamped Langevin + ABF BIN dynamics with Omega + External parameter'
      case (4) 
              if (rangmab==0) write(6,*) ' Overdamped Langevin + ABF GAUSSIAN dynamics'
      case (5) 
              if (rangmab==0) write(6,*) ' Overdamped Langevin + ABF EE dynamics'
      case (6) 
              if (rangmab==0) write(6,*) ' Overdamped Langevin + ABF EE dynamics + constant biais'
      case (7) 
              if (rangmab==0) write(6,*) ' Overdamped Langevin + ABF BIN dynamics + constant biais'
      case (8) 
              if (rangmab==0) write(6,*) ' Overdamped Langevin + ABF EE dynamics + iterative'
   end select 
 end if


if (block)  then 
  if (rangmab==0) write(6,'("MAB: WARNING -->>> Some atoms are in protective domains! <<-- WARNING")')
   if (block_file) then
     if (rangmab==0) write(6,'("MAB:  Atoms in protective domains are readed in on *.mab.lblock file!")')
     fnamt_lblock = fnam(1:lenfnam)//'.mab.lblock'
     !debug  write(*,*) fnamtin
     lublock=775
     open(unit=lublock, file=fnamt_lblock, status='unknown')
     read(lublock,*) nsite_block
     allocate (isite_block(nsite_block))
     do ii=1,nsite_block
     !
      read(lublock,*)  isite_block(ii)
      if (isite_block(ii) > im) then
        if (rang==0) then
          write(*,*) 'MAB: this atom cannot be blocked', isite_block(ii)
          write(*,*) 'MAB: the value exceeds the number of atoms im ', im
          write(*,*) 'MAB: stop in read_mab_file.f90'
        end if 
        stop
      end if  
      if (abf_mode==1) then
        if (isite_block(ii)==atom_to_jump) then
          if (rang==0) then
            write(*,*) 'MAB: The atom which is designed to jump is BLOCKED by the list *.mab.lblock ', isite_block(ii)
            write(*,*) 'MAB: The atom to jump is set by atom_to_jump, currently set to ', atom_to_jump
            write(*,*) 'MAB: stop in read_mab_file.f90'
          end if 
          stop
        end if
      end if
      !
     end do   
     close(lublock)
  else ! block_file  == .false.
     if (rangmab==0) write(6,'("MAB:  ALL the atoms are in protective domains *.mab.lblock file is IGNORATED!")')
     nsite_block=im
     allocate (isite_block(im))
     do ii=1,nsite_block
       isite_block(ii)=ii
     end do

     if (abf_mode==abf_mode_reaction) then
      if (rangmab==0) then
       write(6,'("MAB: block, block_file, abf_mode are not compatible")')
       write(6,'("MAB: Change the type of abf_mode others than 1")')
       write(6,'("MAB: ... or switch block_file to T and provide an *.mab.lblock file which is compatible with")') 
       write(6,'("MAB: your abf_mode reaction coordinate 1")')
       write(6,*) 'MAB: stop in read_mab_file.f90'
       end if 
       stop
     end if 
  end if !block 

end if 

     if  (.not.( (abf_mode/=abf_mode_alchemical).or.  &
           (abf_mode/=abf_mode_temperature).or. &
           (abf_mode/=abf_mode_reaction))  ) then
         if (rangmab==0) then
         write(6,*) 'MAB: No implementation for this ABF reaction coordinate which can be:'
         write(6,*) 'MAB: abf_mode = 1 for  geometric  reaction coordiante'
         write(6,*) 'MAB: abf_mode = 2 for  alchemical reaction  coordiante'
         write(6,*) 'MAB: abf_mode = 22 for temperature reaction  coordiante'
         write(6,*) 'MAB: abf_mode =  ', abf_mode
         write(6,*) 'MAB: stop in read_mab_file'
         end if 
         stop
     end if 

     if  ((abf_mode==abf_mode_alchemical) .and. ((abf_type==1).or.(abf_type==4).or. &
                             (abf_type==6).or.(abf_type==7).or.  &
                             (abf_type==9))) then
         write(6,*) 'MAB: The is no ABF implementation for this mode'
         write(6,*) 'MAB: abf_type .....',abf_type
         write(6,*) 'MAB: abf_mode .....',abf_mode
         if (abf_type==1) write(6,*) 'MAB: You cannot use Langevin dynamics having abf_mode in the input file.' 
         if (abf_type==1) write(6,*) 'MAB: Put abf_mode=1 and restart the calcultations' 
         write(6,*) '<stop in read_mab_file>'
         stop
    end if 


     if ((abf_mode==abf_mode_temperature) .and. ((abf_type==1).or.(abf_type==4).or. &
                             (abf_type==6).or.(abf_type==7).or.  &
                             (abf_type==9))) then
         write(6,*) 'MAB: The is no ABF implementation for this mode'
         write(6,*) 'MAB: abf_type .....',abf_type
         write(6,*) 'MAB: abf_mode .....',abf_mode
         if (abf_type==1) write(6,*) 'MAB: You cannot use Langevin dynamics having abf_mode in the input file.' 
         if (abf_type==1) write(6,*) 'MAB: Put abf_mode=1 and restart the calcultations' 
         write(6,*) '<stop in read_mab_file>'
         stop
    end if 

    if ((abf_mode==abf_mode_alchemical).or.(abf_mode==abf_mode_temperature)) then
       allocate (omega_veinstein(3,im))

      if (itype_einstein==0 ) then
       if (rangmab==0) write(6,'("MAB: Einstein frequency (omega_einstein)...................:",D15.4)') omega_einstein
       omega_veinstein(:,:)=omega_einstein
      end if 

      if (itype_einstein==1) then

       units_phondy = erg2eV*A2cm*A2cm*umass
       allocate (matfor(im,im,3,3))
       fnamt_lfreq = fnam(1:lenfnam)//'.mab.lfreq'
       fnamt_lm = fnam(1:lenfnam)//'.mab.m'
       fnamt_lv = fnam(1:lenfnam)//'.mab.v'
       fnamt_lu = fnam(1:lenfnam)//'.mab.u'
       lcm=71
       lcu=72
       lcv=73
       open(unit=lcm, file=fnamt_lm, status='unknown', form='unformatted')
       open(unit=lcu, file=fnamt_lu, status='unknown', form='unformatted')
       open(unit=lcv, file=fnamt_lv, status='unknown', form='unformatted')
       open(unit=lublock, file=fnamt_lfreq, status='unknown')
       do ii=1,im
          do ia=1,3
           read(lublock,*)  itemp_read, omega_veinstein(ia,ii)
          end do
       end do
       close(lublock)

       read (lcm) imax
       do ii=1,imax
        read(lcu) iu
        read(lcv) iv
        ia=mod(iu,3)+1
        ja=mod(iv,3)+1
         i=(iu-ia)/3+1
         j=(iv-ja)/3+1
        read(lcm) temp_read
          matfor (i,j,ia,ja) = temp_read *dsqrt(m_i(ia,i)*m_i(ja,j))/units_phondy
       end do

      close(lcu)
      close(lcm)
      close(lcv)

     end if !  itype_einstein==1

    end if ! abf_mode == 2 and 22


        !instead that I will a file with all the einstein  frequencies 

     select case (sim_mode)
      case (1) 
           if (rangmab==0) write(6,'("MAB:The simulation check the first passage time and   ")')
           if (rangmab==0) write(6,'("MAB: MAB will stop once the vacancy goes into the final postion ")')
      case (2) 
           if (rangmab==0) write(6,'("MAB: The simulation stops after nlangevin steps.....................:",i9)')  nlangevin
     end select 




if (rangmab==0) write(6,'("MAB: a0 of the cubic unit cell............................:",D15.4)') a0bcc 
if (rangmab==0) write(6,'("MAB: Langevin time step in s..............................:",D15.4)') dtlang 
if (rangmab==0) write(6,'("MAB: Total number of steps ...............................:",I9)')  nlangevin
if (rangmab==0) write(6,'("MAB: Langevin temperature in K............................:",F8.1)') temperature
if (rangmab==0) write(6,'("MAB: Langevin dumping coefficient ........................:",D15.4)') gamma
if (rangmab==0) write(6,'("MAB: Factor for Langevin coefficient .....................:",D15.4)') lang_factor
if (rangmab==0) write(6,'("MAB: Debug verbosity idebug ..............................:",I9)') idebug


if (abf_type==3) then
 if (rangmab==0) write(6,'("MAB: Omega ABF BIN.........................................:",E25.12E3)') omega_abf
end if

if (abf_type==4) then
 if (rangmab==0) write(6,'("MAB: eta_mab the width of the Gaussian in bins............:",E25.12E3)') eta_mab 
end if

if (block) then
  if (rangmab==0) write(6,'("MAB: Radius of the blocking spheres (1nn units) ..........:",f15.4)') radiussph
  if (rangmab==0) write(6,'("MAB: Width of the FD function in A........................:",f15.4)') deltasph
end if

if (rangmab==0) write(6,'("MAB: Number of the bins of histo.........................:",i7)') nhisto
if (rangmab==0) write(6,'("MAB: The frequency of writing histo......................:",i7)') nwrite_histo
if (rangmab==0) write(6,'("MAB: The first shell of the histo (1nn units)............:",f15.4)') deltar1
if (rangmab==0) write(6,'("MAB: The second shell of the histo (1nn units)...........:",f15.4)') deltar2
if ((abf_mode==abf_mode_temperature).or.(abf_mode==abf_mode_alchemical)) then
 if (mode_zeta_potential==1) then
  if (rangmab==0) write(6,'("MAB: ===============THERE IS AN EXTRA POTENTIAL FOR ZETA===========")') 
  if (rangmab==0) write(6,'("MAB: The prefactor of zeta potential.......................:",E25.12E3)') alpha_zeta
  alpha_zeta=alpha_zeta*ev2erg
 end if 
end if

if (abf_mode==22) then
if (rangmab==0) write(6,'("MAB: ha_mix, U(\zeta,q)=\zeta*[U(q)+ha_mix*U_HA(q))].......:",f15.4)') ha_mix
if (rangmab==0) write(6,'("MAB: Temperature min \zeta.................................:",f15.4)') temperature_zeta_min
if (rangmab==0) write(6,'("MAB Temperature max \zeta..................................:",f15.4)') temperature_zeta_max
end if 
 
if (rangmab==0) write(6,'("MAB: The cutoff radius for ending sim (1nn unit)B.........:",D15.4)') rtestlac

 temperature=temperature*KtoERG
 temperature_zeta_min=temperature_zeta_min*KtoERG
 temperature_zeta_max=temperature_zeta_max*KtoERG
 deltasph=deltasph/angst
 radiussph=dsqrt(3.d0)*a0bcc*radiussph/(angst*2.d0)

    if (deltar2 < deltar1) then
     write(6,*) 'MAB: Increase deltar2. deltar2 should be greater than deltar1'
     write(6,*) 'MAB: deltar2,deltar1',deltar2,deltar1
     stop
    end if 



close (lumab)

 if (abf_restart) then 
  if (.not.((abf_mode==abf_mode_reaction).or.(abf_mode==abf_mode_alchemical))) then
   if (rangmab==0) write(6,'("MAB: The restart mode is not implemented for this abf_mode", i6)') abf_mode
   if (rangmab==0) write(6,'("MAB stop in read_mab_file.f90")')
   stop
  end if 
 end if 




end subroutine read_mab_file


subroutine print_mab(rangloc)
implicit none 
integer, intent(in) :: rangloc
character(len=1) :: quote,dquote


 quote=char(39)
dquote=char(34)

if (rangloc==0) then

write(6,'("/---------------------------------------------------------\")') 
write(6,'("                                                           ")')
write(6,'("               __  __   ___             /|                 ")')
write(6,'("              |  |/  `.",a,"   `.           ||                 ")')quote
write(6,'("              |   .-.  .-.   ",a,"          ||                 ")')quote
write(6,'("              |  |  |  |  |  |    __    ||  __             ")')
write(6,'("              |  |  |  |  |  | .:--.",a,".  ||/",a,"__ ",a,".          ")')quote,quote,quote
write(6,'("              |  |  |  |  |  |/ |   \ | |:/`  ",a,". ",a,"         ")')quote,quote
write(6,'("              |  |  |  |  |  |`",a," __ | | ||     | |         ")')dquote
write(6,'("              |__|  |__|  |__| .",a,".",a,a"| | ||\    / ",a,"         ")')quote,quote,quote,quote
write(6,'("                              / /   | |_|/\",a,"..",a," /          ")')quote,quote
write(6,'("                              \ \._,\ ",a,"/",a,"  `",a,"-",a,"`           ")')quote,quote,quote,quote
write(6,'("                               `--",a,"  `",a,"                    ")')quote,dquote
write(6,'("copyleft CEA by ...                                        ")')
!write(6,'("... M.-C. Marinica, M. Athenes                             ")')
!write(6,'("others contributions                                       ")')
!write(6,'("L. Cao, G. Stoltz, T. Lelievre                             ")')
write(6,'("email:mihai-cosmin.marinica@cea.fr                         ")')
write(6,'("\---------------------------------------------------------/")') 

end if 


end subroutine print_mab



subroutine read_gin_file
 USE T_kind_param_m, ONLY:  double
 use gen_com_m, ONLY: lenfnam,fnam,im,imm,erg2ev,at
 use tab_imm_m, ONLY: xp,ityp
 use var_pot,   ONLY: cm
 use gin_module 
 USE mab_in_ndm_module, ONLY: a0bcc, xp0, block, block_file, abf_restart, rangmab, &
                              atom_to_jump
 implicit none
 character(len=89) :: ginFile_block
 integer :: im_local, i , j
 real(double), dimension(3,3) :: at_local
 integer     :: ityp_local(imm)
 real(double) ::  xp_local(3,imm)
 real(double) :: atmp
 logical :: ok
 
 
! ginFile_block=fnam(1:lenfnam)//'.gin'
! if (abf_restart) 
 ginFile_block=fnam(1:lenfnam)//'.gin.block'
 inquire(file=ginFile_block, exist=ok)

    IF (ok ) THEN
            write(6,'("MAB: Read gin file for blocking spheres.....",a)') TRIM(ginFile_block)
            OPEN(unit=93, file=ginFile_block, status='old', action='read')
            CALL ReadGin(xp_local(:,:), ityp_local(:), im_local, at_local, 93)
            CLOSE(93)
    ELSE
            WRITE(0,'("MAB: Does not manage to find ...", a)') TRIM(ginFile_block)
            WRITE(0,'(3a)') 'File ', Trim(ginFile_block), ' does not exist'
            STOP '< read_gin_file >'
    END IF

atmp=0.d0

  do i=1,3
    do j=1,3
      atmp = atmp + (at_local(i,j)*1.d8-at(i,j)*1.d8)**2
    end do
  end do

  if (atmp >= 1.d-8) then
    if (rangmab==0) write(6,'("MAB: the at box of the *.gin file and .gin.block file is not the same")')
    stop ' < read_gin_file >'
  end if 


  if ((im - im_local)**2/=0) then
    if (rangmab==0) write(6,'("MAB: the im no of atoms of the *.gin file and .gin.block file is not the same")')
    stop ' < read_gin_file >'
  end if 

  do i=1,im
   atmp = dsqrt(SUM(xp(:,i)-xp_local(:,i))**2)
   if (i/=atom_to_jump) then
     if (atmp*1.d8>=(a0bcc*dsqrt(3.d0)/2.d0)) then
       if (rangmab==0) write(6,'("MAB: the position of atom no",i0," of the *.gin file and .gin.block file is not related")') i
       stop ' < read_gin_file >'
     end if 
   end if 
  end do

xp0(:,:) = xp_local(:,:)

return


end subroutine read_gin_file
