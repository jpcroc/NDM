module neigcel_mod
  USE gen_com_m, ONLY:deltadist,nox,noxyz,noy,noz,ncel
  implicit none
contains
  !**************************************************************
  !                                                             *
  !                       subroutine cells                      *
  !                        forwards cells                       *
  !                    version du 28 septembre 2000             *
  !**************************************************************
  subroutine neigcel
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double



    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: kx, ky, kz, koo, l, lz, mz, ly, my, lx, mx, kxy

    !-----------------------------------------------

    !     write(6,*)'entree neigcel'

    if (noxyz==1) then
       ncel(1,0)=1
       deltadist=0
    else

       do kz = 1, noz
          do ky = 1, noy
             do kx = 1, nox
                koo = 1+(kx-1)+nox*((ky-1)+noy*(kz-1))

                ncel(koo,0) = koo
                deltadist(:,0,koo) = 0

                l = 1
                do lz = -1, 1
                   do ly = -1, 1
                      do lx = -1, 1
                         deltadist(:,l,koo) = 0

                         mz = kz+lz
                         if (mz<1) then
                            mz = mz+noz
                            deltadist(3,l,koo) = 1
                         endif
                         if (mz>noz) then
                            mz = mz-noz
                            deltadist(3,l,koo) = -1
                         endif

                         my = ky+ly
                         if (my<1) then
                            my = my+noy
                            deltadist(2,l,koo) = 1
                         endif
                         if (my>noy) then
                            my = my-noy
                            deltadist(2,l,koo) = -1
                         endif

                         mx = kx-lx
                         if (mx<1) then
                            mx = mx+nox
                            deltadist(1,l,koo) = 1
                         endif
                         if (mx>nox) then
                            mx = mx-nox
                            deltadist(1,l,koo) = -1
                         endif

                         kxy = 1+(mx-1)+nox*((my-1)+noy*(mz-1))
                         if (kxy==koo) cycle
                         ncel(koo,l) = kxy
                         !                        write(6,*)koo,lz,ly,lx,l,kxy
                         !                        if ((kz==noz).and.(lz==1))write(6,*)koo,lz,l,kxy
                         !                        if ((kz==1).and.(lz==-1))write(6,*)koo,lz,l,kxy
                         l = l+1
                      end do
                   end do
                end do
             end do
          end do
       end do
    end if
    !      if (ltranche) then
    !         do kz = 1, noz
    !            do ky = 1, noy
    !               do kx = 1, nox
    !                  if (kz==noz) then
    !                     koo = 1+(kx-1)+nox*((ky-1)+noy*(kz-1))
    !                     ncel(koo,1) = 0
    !                     ncel(koo,2) = 0
    !                     ncel(koo,3) = 0
    !                     ncel(koo,4) = 0
    !                     ncel(koo,5) = 0
    !                     ncel(koo,6) = 0
    !                     ncel(koo,7) = 0
    !                     ncel(koo,8) = 0
    !                     ncel(koo,9) = 0
    !                  endif
    !                  if (kz/=1) cycle
    !                  koo = 1+(kx-1)+nox*((ky-1)+noy*(kz-1))
    !                  ncel(koo,18) = 0
    !                  ncel(koo,19) = 0
    !                  ncel(koo,20) = 0
    !                  ncel(koo,21) = 0
    !                  ncel(koo,22) = 0
    !                  ncel(koo,23) = 0
    !                  ncel(koo,24) = 0
    !                  ncel(koo,25) = 0
    !                  ncel(koo,26) = 0
    ! 
    !               end do
    !            end do
    !         end do
    !     endif

    return
  end subroutine neigcel
end module neigcel_mod
