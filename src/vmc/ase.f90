module ase_mod
    implicit None
    public :: export_forces
    PRIVATE
    save

contains

    subroutine export_forces()

        use force_fin, only: da_energy_ave
        use atom,      only: ncent
        use forcewt,   only: wcum
        use estcum,    only: ecum

        integer :: i

        open(75, FILE = 'write_forces',form='formatted',status='replace')

        write(75,'(1p6e14.5)') (ecum(1,1)/wcum(1,1)) 

        do i = 1, ncent
            write(75,'(1p6e14.5)') (da_energy_ave(:,i)) 
        enddo

        close(75)

    end subroutine export_forces

end module ase_mod