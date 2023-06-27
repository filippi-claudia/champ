module qmmm_pot
contains
        subroutine qmmm_extpot_read
        implicit none
        return
        end

        subroutine qmmm_extpot_ene(coord,nelec,ext_pot)
      use precision_kinds, only: dp
        implicit none
        integer :: nelec
        real(dp) :: ext_pot
        real(dp), dimension(*) :: coord
        return
        end

        subroutine qmmm_extpot_final(nelec)
        implicit none
        integer :: nelec
        return
        end
end module
