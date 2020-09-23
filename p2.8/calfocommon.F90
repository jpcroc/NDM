module calfocommon
  USE T_kind_param_m
  implicit none
  real(double)::potist,sig(3,3)
  real(double),pointer:: eat(:), sigat(:,:,:),sigc(:,:,:)
  logical ::lsigat,lprteat
  logical::test_sigma,ltpcel
end module calfocommon
