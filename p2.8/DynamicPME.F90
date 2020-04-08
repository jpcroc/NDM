module dynallocPME
        implicit none
        contains
subroutine DynamicalAllocationPME

  USE gen_com_m, ONLY:imm
  USE var_pot, ONLY:bsmod1,bsmod2,bsmod3,table,iiim,ijim,ikim,fr1,fr2,fr3,de1,de2,de3,kpmex,kpmey,kpmez,maxorder,ntable
  implicit none

  allocate(bsmod1(kpmex))
  allocate(bsmod2(kpmey))
  allocate(bsmod3(kpmez))
  allocate(table(ntable,3))
  allocate(iiim(maxorder,imm))
  allocate(ijim(maxorder,imm))
  allocate(ikim(maxorder,imm))
  allocate(fr1(imm),fr2(imm),fr3(imm))
  allocate(de1(imm),de2(imm),de3(imm))
  !  allocate(w1pme(imm),w2pme(imm),w3pme(imm))

end subroutine DynamicalAllocationPME
end module
