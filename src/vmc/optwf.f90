module optwf_mod
    implicit None
    public :: optwf
    PRIVATE
    save

contains

    subroutine optwf()

        !> Main switch for optimization
        use method_opt, only: method
        use optwf_contrl, only: ioptwf, idl_flag, ilbfgs_flag
        use optwf_dl_mod, only: optwf_dl
        use optwf_sr_mod, only: optwf_sr
        ! debug lines below  ravindra
        use dets, only: ndet
        use coefs, only: norb
        use atom, only: nctype_tot, ncent_tot
        use csfs, only: nstates
        use const, only: nelec
        use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE
        implicit None

        if (ioptwf .gt. 0) then
            if (idl_flag .gt. 0) then
                call optwf_dl()
            elseif (ilbfgs_flag .gt. 0) then
                call optwf_olbfgs
            elseif (method .eq. 'sr_n') then
                call optwf_sr
            elseif (method .eq. 'lin_d') then
                call optwf_lin_d
            elseif (method .eq. 'mix_n') then
                call optwf_mix
            else
                call optwf_matrix_corsamp
            endif
        else
            print*, "printing nctype_tot, ncent_tot ", nctype_tot, ncent_tot
            print*, "printing ndet ", ndet
            print*, "printing norb ", norb
            print*, "printing nstates ", nstates
            print*, "printing nelec ", nelec
            print*, "MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE"
            print*, MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE
            call vmc
        endif

    end subroutine optwf

end module optwf_mod
