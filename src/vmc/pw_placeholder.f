      subroutine do_read_lattice(iu)
      implicit none

      integer :: iu

      return
      end

      subroutine pw_setup_input
      implicit none

      return
      end

      subroutine read_orb_pw
      implicit none

      return
      end

      subroutine read_orb_pw_tm
      implicit none

      return
      end

      subroutine set_ewald
      implicit none

      return
      end

      subroutine find_image3(r,rnorm)
      use precision_kinds, only: dp
      implicit none


      real(dp) :: r, rnorm
      return
      end

      subroutine find_image4(rs,r,rnorm)
      use precision_kinds, only: dp
      implicit none


      real(dp) :: r, rnorm, rs
      return
      end

      subroutine orbitals_pw(x,orb,dorb,ddorb)
      use precision_kinds, only: dp
      implicit none


      real(dp) :: ddorb, dorb, orb, x
      return
      end

      subroutine orbitals_pwe(iel,x,orb)
      use precision_kinds, only: dp
      implicit none

      integer :: iel
      real(dp) :: orb, x
      return
      end

      subroutine orbitals_pw_grade(iel,x,orb,dorb,ddorb)
      use precision_kinds, only: dp
      implicit none

      integer :: iel
      real(dp) :: ddorb, dorb, orb, x
      return
      end

      subroutine pot_nn_ewald
      implicit none

      return
      end

      subroutine pot_en_ewald(x,pe_en)
      use precision_kinds, only: dp
      implicit none


      real(dp) :: pe_en, x
      return
      end

      subroutine pot_ee_ewald(x,pe_ee)
      use precision_kinds, only: dp
      implicit none


      real(dp) :: pe_ee, x
      return
      end
