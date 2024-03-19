module calfocommon
  USE T_kind_param_m
  implicit none
  real(double)::potistcalfo,sigcalfo(3,3),sig2p(3,3)
  real(double),pointer:: sigc(:,:,:)
  logical ::lsigat,lprteat
  logical::test_sigma,lcalcsigc
end module calfocommon
