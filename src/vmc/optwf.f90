module optwf_mod
    implicit None
    public :: optwf
    PRIVATE
    save

contains

    subroutine optwf()

        !> Main switch for optimization
      use optwf_control, only: idl_flag,ilbfgs_flag,ioptwf,method
      use optwf_dl_mod, only: optwf_dl
      use optwf_lin_dav, only: optwf_lin_d
      use optwf_matrix_corsamp_mod, only: optwf_matrix_corsamp
      use optwf_mix_sa, only: optwf_mix
      use optwf_olbfgs_mod, only: optwf_olbfgs
      use optwf_sr_mod, only: optwf_sr
      use vmc_f_mod, only: vmc

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
            call vmc
        endif

    end subroutine optwf

end module optwf_mod
