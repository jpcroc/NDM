! nox,noy,noz, natperc doivent être connus pour l'instant

module boxconfig
  USE T_kind_param_m


  implicit none
  type box_config
     real(double):: at(3,3)
     real(double):: bg(3,3)
     real(double):: zl,zls2,nzl,volu,normat(3)
     integer(long)::icaltabt 
  end type box_config

end module boxconfig


