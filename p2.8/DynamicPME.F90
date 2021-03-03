module dynallocPME
!  USE temp_com,only:imm ! A EFFACER
        implicit none
        contains
subroutine DynamicalAllocationPME(immT)

  !  USE gen_com_m, ONLY:imm
  
  USE var_pot, ONLY:bsmod1,bsmod2,bsmod3,table,iiim,ijim,ikim,fr1,fr2,fr3,de1,de2,de3,kpmex,kpmey,kpmez,maxorder,ntable
  implicit none
  integer,intent(in)::immT
  allocate(bsmod1(kpmex))
  allocate(bsmod2(kpmey))
  allocate(bsmod3(kpmez))
  allocate(table(ntable,3))
  allocate(iiim(maxorder,immT))
  allocate(ijim(maxorder,immT))
  allocate(ikim(maxorder,immT))
  allocate(fr1(immT),fr2(immT),fr3(immT))
  allocate(de1(immT),de2(immT),de3(immT))
  !  allocate(w1pme(immT),w2pme(immT),w3pme(immT))

end subroutine DynamicalAllocationPME
end module
