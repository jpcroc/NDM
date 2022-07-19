implicit none 
integer :: nb_procsph, i_start, i_final,im
integer :: nratio1,nratio2,rangph,rest

write(*,*) 'The number of atoms'
read(*,*) im
write(*,*) 'The number of possibles procs'
read(*,*) nb_procsph



rest=mod(im,nb_procsph)

write(*,*) rest

do rangph=0,nb_procsph-1
 !
  nratio1=im/nb_procsph+1
  nratio2=im/nb_procsph
 if (rangph<rest) then
   i_start=rangph*nratio1+1
   if (rangph /= (nb_procsph-1) ) i_final=(rangph+1)*nratio1
 end if

 if (rangph>=rest) then
 i_start = rest*nratio1 + (rangph-rest)*nratio2 + 1
 if (rangph /= (nb_procsph-1) ) i_final=rest*nratio1 + (rangph-rest+1)*nratio2
 endif 
  


 if (rangph == (nb_procsph-1) ) i_final= im 
  write(*,*) rangph, i_start, i_final, i_final-i_start
 !
end do


end


