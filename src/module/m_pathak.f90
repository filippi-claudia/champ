module pathak_mod
    use precision_kinds, only: dp
    use dmc_mod, only: mwalk
    use force_pth, only: PTH
  
      implicit none
  
      integer :: ipathak
      real(dp) :: eps_max, deps
      real(dp), dimension(:), allocatable :: eps_pathak, pnew !(PTH)
      real(dp), dimension(:, :), allocatable :: pold !(MWALK, PTH)
  
  
      private
      public :: ipathak, eps_max, deps, pnew
      public :: eps_pathak, pold
      public :: init_pathak, init_eps_pathak
      public :: allocate_pathak, deallocate_pathak, pathak
  
    contains
  
      subroutine init_pathak()
        if (ipathak .gt. 0) then
           PTH = ipathak
        else
           PTH = 1
        endif
      end subroutine init_pathak
  
      subroutine init_eps_pathak()
        use contrl_file, only: ounit
        use mpiconf, only: wid
  
        integer :: iph
  
        if (ipathak .gt. 0) then
           do iph=1,PTH
              eps_pathak(iph) = eps_max - deps * (iph - 1)
           enddo ![iph=1,PTH]
  !         eps_pathak(PTH) = 0.d0
           if (wid) write(ounit,*) 'eps ', eps_pathak
        endif      
      end subroutine init_eps_pathak
  
      subroutine allocate_pathak()      
        if (.not. allocated(eps_pathak)) allocate (eps_pathak(PTH))
        if (.not. allocated(pnew)) allocate (pnew(PTH))
        if (.not. allocated(pold)) allocate (pold(mwalk, PTH))
      end subroutine allocate_pathak
  
      subroutine deallocate_pathak()
        if (allocated(eps_pathak)) deallocate (eps_pathak)
        if (allocated(pnew)) deallocate (pnew)
        if (allocated(pold)) deallocate (pold)
      end subroutine deallocate_pathak
  
      subroutine pathak(distance, f, epsilon)
        use precision_kinds, only: dp
  
        implicit none
  
        real(dp) :: x1, x2, x4, x6
        real(dp) :: f, distance, epsilon
  
        x1 = distance / epsilon
  
        if (x1 .gt. 1.d0) then
           f = 1.d0
        else
           x2 = x1 * x1
           x4 = x2 * x2
           x6 = x2 * x4
           f = 7.d0*x6 - 15.d0*x4 + 9.d0*x2
        endif
  
        return
      end subroutine pathak
  
  
  end module pathak_mod