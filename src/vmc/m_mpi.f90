
 module mpi_qmc
  !> Arguments NPROCX
  integer, parameter :: NPROCX=1024
  private
  public :: NPROCX
  save
 end module mpi_qmc

 module mpiconf
   !> Arguments: idtask, nproc, wid
    integer  :: idtask
    integer  :: nproc
    logical  :: wid 

    private 
    public :: idtask, nproc, wid 
    save
 end module mpiconf