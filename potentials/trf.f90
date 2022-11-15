program trf
real*8:: xr,yr1,yr2,yr3
read(5,*)nl
do i=1,nl
read(10,*)xr,yr1,yr2,yr3
write(11,'(4G20.12)')xr,yr1,yr2,yr3
end do
end program trf
