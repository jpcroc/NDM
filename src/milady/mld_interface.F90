module mld_interface_mod
    use var_pot
    USE atomconfig,only:atom_config
    USE cellconfig,only:cell_config
    USE boxconfig, only: box_config
    use mod_test, only: md_calfo_ml, ndm2mld_var_pot, ndm2mld_gen_comm, mld2ndm_var_pot, copy_fnam, md_config_wrap
    USE NDM_ML,only :rue_ml
    use gen_com_m, only:dmtype

    contains
    
    subroutine mld_calfo(atcf,boxcf,potistcalfo,sigcalfo,celcf)
    
        implicit none
        type(atom_config), INTENT(INOUT):: atcf
        type(cell_config), INTENT(IN):: celcf
        type(box_config), INTENT(INOUT):: boxcf
        real(double), INTENT(INOUT)::potistcalfo, sigcalfo(3,3)
    
        call ndm2mld_var_pot(ipotentiel, rue_pot, typ_pot_pair, npair, &
                        npotentiel, lprtpot, ntyp, &
                        ntyp_buffer, ntrip, eatref, lpotentiel, ipotrep, ngr, &
                        l3c, iewald,ngrid, eta, rumax, &
                        csive, ncouc3, n2max, ncoucx, ncoucy, ncoucz, nvecttot, &
                        precisew, epswat, sigmawat, gm1, gm2, gm3, gm4, gm5, &
                        gR, gd, csive_g, r3cm, r3cm2, rbp5, rp5p3, rp3c, &
                        alpha, lambda, xsi, potisrep, potisglue, potiseam, &
                        kpmex, kpmey, kpmez, kpme, maxorder, iorder, npoint, nfft1, &
                        nfft2, nfft3, nff, nf1, nf2, nf3, ntable, pterm, volterm)
        call md_calfo_ml(atcf%im,atcf%imm,atcf%ityp,atcf%xp,atcf%fp,boxcf%volu,boxcf%at,boxcf%bg, potistcalfo,sigcalfo,celcf%ncel)
    
    end subroutine mld_calfo

    subroutine mld_init_potential
        use gen_com_m, only: rang
        use read_val,only:rvois

        implicit none
        call ndm2mld_var_pot(ipotentiel, rue_pot, typ_pot_pair, npair, &
                        npotentiel, lprtpot, ntyp, &
                        ntyp_buffer, ntrip, eatref, lpotentiel, ipotrep, ngr, &
                        l3c, iewald,ngrid, eta, rumax, &
                        csive, ncouc3, n2max, ncoucx, ncoucy, ncoucz, nvecttot, &
                        precisew, epswat, sigmawat, gm1, gm2, gm3, gm4, gm5, &
                        gR, gd, csive_g, r3cm, r3cm2, rbp5, rp5p3, rp3c, &
                        alpha, lambda, xsi, potisrep, potisglue, potiseam, &
                        kpmex, kpmey, kpmez, kpme, maxorder, iorder, npoint, nfft1, &
                        nfft2, nfft3, nff, nf1, nf2, nf3, ntable, pterm, volterm)
        call ndm2mld_gen_comm(dmtype, rvois,rang)
        call md_init_potential_ml
        call mld2ndm_var_pot(ipo,cm, typ_pot_pair, npair)
    end subroutine mld_init_potential


    subroutine mld_init_config(atdml)
        implicit none
        type(atom_config), INTENT(INOUT):: atdml

        call md_config_wrap(atdml%im, atdml%imm)
        call md_init_config_ml
    
    end subroutine mld_init_config

    subroutine mld_copy_fnam(fnam)
        character, intent(in) :: fnam*80
        call copy_fnam(fnam)
    end subroutine mld_copy_fnam


end module mld_interface_mod