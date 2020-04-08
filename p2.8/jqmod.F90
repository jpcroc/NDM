module jqmod
  USE T_kind_param_m
  USE gen_com_m, ONLY:
  integer,save ::icall
  real(double) :: jqf,jqfk,jqfp
  real(double) ::jq(3),jqk(3),jqp(3),expvect(3)
  real(double),save ::jq0(3),jqk0(3),jqp0(3)
  character ::fnamjq*80,fnamjqbis*80

end module jqmod
