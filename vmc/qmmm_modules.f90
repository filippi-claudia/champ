        MODULE qmmm_density
        character (len=80) :: filename_dens,title_dens(2)
        integer :: n_xd,n_yd,n_zd
        integer :: n_atomsd,cdipole
        double precision :: x0d(3),deltad(3),Ang
        double precision outofbox,inofbox,outofboxs,inofboxs
        double precision :: cc_nuc(3),cc_ele(3),cc_ele2(3),dipole(3),dipole2(3)
        integer , dimension(:), allocatable :: id_atomd
        double precision, dimension(:,:), allocatable :: x_atomd
        double precision, dimension(:), allocatable :: chrg_atomd
        double precision, dimension(:,:,:), allocatable :: dens,sme
        parameter ( Ang = 0.529177249 )
        end

        MODULE qmmm_extpot
        character (len=80) :: filename_cube,title_cube(2)
        integer :: n_x,n_y,n_z
        integer :: n_atoms
        double precision :: x0(3),delta(3)
        integer , dimension(:), allocatable :: id_atom
        double precision, dimension(:,:), allocatable :: x_atom
        double precision, dimension(:), allocatable :: chrg_atom
        double precision, dimension(:,:,:), allocatable :: pot
        double precision, dimension(:), allocatable :: xdata,ydata,zdata
        integer :: nout,ncount
        double precision :: ave, ave2,err
        end


        MODULE qmmm_splines
        integer :: kxord, kyord, kzord, nxknot, nyknot,nzknot
        integer ::  nxcoef, nycoef, nzcoef
        double precision :: x,y,z
        double precision :: deltavec(3)
        double precision,dimension(:,:,:), allocatable ::  bscoef
        double precision, dimension(:),allocatable :: &
     &        xknot,yknot,zknot
        end

        MODULE qmmm_vector
        integer ::  nxvec,nyvec,nzvec
        logical :: first_call_vec
        double precision, dimension(:),allocatable :: xvec,yvec,zvec
        double precision,dimension(:,:,:), allocatable :: value
        end
