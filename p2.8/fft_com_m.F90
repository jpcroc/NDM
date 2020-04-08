module fft_com_m
  USE T_kind_param_m

  integer*4 :: STATUS
#ifdef para2c
  complex(double), dimension(:,:,:), pointer :: qgrid1   !tampon de qgrid
#endif
  complex(double), dimension(:,:,:), pointer :: qgrid   !tampon de qgrid

end module fft_com_m
