!> Subroutines to allocate/deallocate memory of specific VMC-related
!> variables.
!> The allocation will occurr only in via read_input.f
subroutine allocate_vmc()

      call allocate_m_common()
      call allocate_m_basis
      call allocate_m_control
      call allocate_m_deriv
      call allocate_m_efield
      call allocate_m_estimators
      call allocate_m_ewald
      call allocate_m_force
      call allocate_m_gradhess
      call allocate_m_grdnt
      call allocate_m_grid
      call allocate_m_jastrow
      call allocate_m_mixderiv
      call allocate_m_mmpol
      call allocate_m_mstates
      call allocate_m_optci
      call allocate_m_optorb
      call allocate_m_optwf
      call allocate_m_pcm
      call allocate_m_prop
      call allocate_m_pseudo
      call allocate_m_sampling
      call allocate_m_sr
      call allocate_m_state_avrg

end subroutine allocate_vmc

!subroutine deallocate_vmc()
!
!      call deallocate_m_common
!      call deallocate_m_basis
!      call deallocate_m_control
!      call deallocate_m_deriv
!      call deallocate_m_efield
!      call deallocate_m_estimators
!      call deallocate_m_ewald
!      call deallocate_m_force
!      call deallocate_m_gradhess
!      call deallocate_m_grdnt
!      call deallocate_m_grid
!      call deallocate_m_jastrow
!      call deallocate_m_mixderiv
!      call deallocate_m_mmpol
!      call deallocate_m_mstates
!      call deallocate_m_optci
!      call deallocate_m_optorb
!      call deallocate_m_optwf
!      call deallocate_m_pcm
!      call deallocate_m_prop
!      call deallocate_m_pseudo
!      call deallocate_m_sampling
!      call deallocate_m_sr
!      call deallocate_m_state_avrg
!
!end subroutine deallocate_vmc