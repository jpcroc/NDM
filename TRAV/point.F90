program pointertest
  real*8,pointer::p
  real*8,target::t

  t=4
  !  p=>t
  allocate(p)
  p=3
  write(6,*)'p',p,t
end program pointertest
