      subroutine qmc
      use vmc_f_mod, only: vmc

      implicit none

      call vmc

      return
      end
!---------------------------------------------------
      subroutine reset_configs_start
      use mc_configs, only: mc_configs_start

      implicit none

      call mc_configs_start

      return
      end
