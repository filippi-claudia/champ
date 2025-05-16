module branching
    use branch, only: eest
    use contrldmc, only: tau, ibranching_c
    use error, only: fatal_error
    use fragments, only: ibranching_cfrag, sqrt_nelecfrag, eestfrag, nfrag, ifragelec
    use precision_kinds, only: dp  
    use system, only: nelec
    
    implicit none

    real(dp), parameter :: sqrt_pi_o2 = 0.88622692545d0
    private

    public :: calculate_fratio, calculate_reweight
contains
    !> Calculates the ratio |v̄(R)|/|v(R)| used to cutoff the energy. The calculation method depends on the choice of icut_e. \\
    !> icut_e = 0 (default) uses the method of UNR93 to calculate v̄(R)
    !> icut_e = 1 does not exist (dont ask me why)
    !> icut_e = 2 
    subroutine calculate_fratio(icut_e, adrift, tratio, taunow, sqrt_nelec, v_dmc, eest, e, fratio, eest_i, eloc_i, fratio_i, eestfrag, elocfrag, fratiofrag)
        !> selects the method of cutting off
        integer, intent(in) :: icut_e 
        real(dp), intent(in) :: adrift, tratio, sqrt_nelec, taunow
        real(dp), dimension(3,nelec) :: v_dmc 
        real(dp), intent(in) :: eest 
        real(dp), intent(in) :: e 
        real(dp), dimension(nelec), intent(in) :: eest_i, eloc_i
        real(dp), dimension(nfrag), intent(in) :: eestfrag, elocfrag

        real(dp), intent(out) :: fratio 
        real(dp), dimension(nelec), intent(out) :: fratio_i
        real(dp), dimension(nfrag), intent(out) :: fratiofrag

        integer :: i
        real(dp) :: v2, v2sum, vav2sum
        real(dp) :: vavfac
        real(dp), dimension(nelec) :: v2_i, vav2_i
        real(dp), dimension(nfrag) :: v2frag, vav2frag
        
        select case (icut_e)

        case (0)
            call dmc_eloc_cutoff(v_dmc(:,:), adrift, tratio, vav2sum, v2sum)
            fratio = dsqrt(vav2sum/v2sum)
        
        case (1)


        case (2)
            v2sum = 0
            do i=1,nelec
               v2sum = v2sum + v_dmc(1,i)**2 + v_dmc(2,i)**2 + v_dmc(3,i)**2
            enddo
            fratio = ibranching_c * sqrt(v2sum) / sqrt_nelec

        case(3)
            fratio = ibranching_c * dabs(eest-e) / (0.2d0 * sqrt_nelec)

        ! Electron fragments
        case(-1) 
            do i = 1, nelec
              v2_i(i) = v_dmc(1,i)**2 + v_dmc(2,i)**2 + v_dmc(3,i)**2
              vavfac = (-1.d0 + dsqrt(1.d0 + 2.d0 * adrift * v2_i(i) * tau)) / (adrift * v2_i(i) * tau)
              vav2_i(i) = vavfac**2 * v2_i(i)
              fratio_i(i) = dsqrt(vav2_i(i) / v2_i(i))
            end do

        case(-2) 
            do i = 1, nelec
              fratio_i(i) = ibranching_c * norm2(v_dmc(:,i)) / sqrt_nelec
            end do

        case(-3)
            fratio_i(:) = ibranching_c * dabs(eest_i(:) - eloc_i(:)) / (0.2d0 * sqrt_nelec)

        ! Fragments
        case(-4)
            v2frag = 0.d0
            vav2frag = 0.0d0
            do i = 1, nelec
              v2 = v_dmc(1,i)**2 + v_dmc(2,i)**2 + v_dmc(3,i)**2
              vavfac = (-1.d0 + dsqrt(1.d0 + 2.d0 * adrift * v2 * taunow)) / (adrift * v2 * taunow)
              v2frag(ifragelec(i)) = v2frag(ifragelec(i)) + v2
              vav2frag(ifragelec(i)) = vav2frag(ifragelec(i)) + vavfac**2 * v2
            end do
            fratiofrag(:) = dsqrt(vav2frag / v2frag)

        case(-5)
            v2frag = 0.d0
            do i = 1, nelec
              v2frag(ifragelec(i)) = v2frag(ifragelec(i)) + v_dmc(1,i)**2 + v_dmc(2,i)**2 + v_dmc(3,i)**2
            enddo
            fratiofrag(:) = ibranching_cfrag(:) * sqrt(v2frag(:)) / sqrt_nelecfrag(:)

        case(-6) 
            fratiofrag(:) = ibranching_cfrag(:) * dabs(eestfrag(:) - elocfrag(:)) / (0.2d0 * sqrt_nelecfrag(:))

        case default
            call fatal_error("CALCULATE_FRATIO: The chosen icut_e is not implemented")
        end select
    end subroutine calculate_fratio

    subroutine calculate_reweight(idmc, icut_e, icut_br, taunow, e_cutoff, &
                                  etrial, eest, eold, enew, fratio, fration, &
                                  eest_i, eloco_i, eloc_i, fratio_i, fration_i, &
                                  etrialfrag, eestfrag, elocofrag, elocfrag, fratiofrag, frationfrag, tauefffrag, &
                                  dwt)
                                  
        integer, intent(in) :: idmc, icut_e, icut_br
        real(dp), intent(in) :: taunow, e_cutoff
        real(dp), intent(in) :: etrial, eest, eold, enew, fratio, fration
        real(dp), dimension(nelec), intent(in) :: eest_i, eloco_i, eloc_i, fratio_i, fration_i
        real(dp), dimension(nfrag), intent(in) :: etrialfrag, eestfrag, elocofrag, elocfrag, fratiofrag, frationfrag, tauefffrag
        real(dp), intent(out) :: dwt
        
        integer :: i
        real(dp) :: deo, den, ecuto, ecutn, ewto, ewtn, fratio_aux, expon
        real(dp), dimension(nelec) :: deo_i, den_i, fratio_aux_i
        real(dp), dimension(nfrag) :: deofrag, denfrag, fratio_auxfrag, ewtofrag, ewtnfrag

        real(dp), parameter :: half = 0.5d0

        deo = eest - eold
        den = eest - enew
        ecuto=min(e_cutoff, dabs(deo))
        ecutn=min(e_cutoff, dabs(den))

        select case (icut_e)

        case (0)
            ewto = eest - deo * fratio
            ewtn = eest - den * fration
        case (1)
            ewto =eest - sign(1.d0, deo) * ecuto
            ewtn =eest - sign(1.d0, den) * ecutn
        case (2,3)
            fratio_aux = max(fratio * taunow, 1e-9)
            ewto = eest - deo * sqrt_pi_o2 * derf(fratio_aux) / (fratio_aux)
            fratio_aux = max(fration * taunow, 1e-9)
            ewtn = eest - den * sqrt_pi_o2 * derf(fratio_aux) / (fratio_aux)
        
        ! electron fragments
        case (-1)
            deo_i(:) = eest_i(:) - eloco_i(:)
            den_i(:) = eest_i(:) - eloc_i(:) 
            ewto=eest-sum(deo_i(:) * fratio_i(:))
            ewtn=eest-sum(den_i(:) * fration_i(:))
        case (-2,-3)
            deo_i(:) = eest_i(:) - eloco_i(:)
            den_i(:) = eest_i(:) - eloc_i(:) 
            do i = 1, nelec
                fratio_aux_i(i) = max(fratio_i(i) * taunow, 1e-9)
            end do
            ewto = eest - sum(deo_i(:) * sqrt_pi_o2 * derf(fratio_aux_i(:)) / (fratio_aux_i(:)))
            do i = 1, nelec
                fratio_aux_i(i) = max(fration_i(i) * taunow, 1e-9)
            end do
            ewtn = eest - sum(den_i(:) * sqrt_pi_o2 * derf(fratio_aux_i(:)) / (fratio_aux_i(:)))
        
        ! fragments
        case (-4)
            deofrag(:) = eestfrag(:) - elocofrag(:)
            denfrag(:) = eestfrag(:) - elocfrag(:) 
            ewto = eest - sum(deofrag(:) * fratiofrag(:))
            ewtn = eest - sum(denfrag(:) * frationfrag(:))

        case (-5,-6)
            deofrag(:) = eestfrag(:) - elocofrag(:)
            denfrag(:) = eestfrag(:) - elocfrag(:) 
            do i = 1, nelec
                fratio_auxfrag(ifragelec(i)) = max(fratiofrag(ifragelec(i)) * tauefffrag(ifragelec(i)), 1e-9)
            end do
            ewtofrag = eestfrag - deofrag * sqrt_pi_o2 * derf(fratio_auxfrag(:)) / (fratio_auxfrag(:))
            do i = 1, nelec
                fratio_auxfrag(ifragelec(i)) = max(frationfrag(ifragelec(i)) * tauefffrag(ifragelec(i)), 1e-9)
            end do
            ewtnfrag = eestfrag - denfrag(:) * sqrt_pi_o2 * derf(fratio_auxfrag(:)) / (fratio_auxfrag(:))
        
        case default 
            call fatal_error("CALCULATE_REWEIGHT: The chosen icut_e is not implemented")
        end select

        if(idmc.gt.0) then
            if (nfrag.gt.1) then
                expon = sum((etrialfrag - half * (ewtofrag + ewtnfrag)) * tauefffrag) 
            else
                expon=(etrial - half * (ewto + ewtn)) * taunow
            endif

            if(icut_br.le.0) then
                dwt = dexp(expon)
            else
                dwt = 0.5d0 + 1 / (1 + exp(-4 * expon))
            endif
        endif
    end subroutine calculate_reweight

    subroutine dmc_eloc_cutoff(v, adrift, tratio, vav2sum, v2sum)  
        integer  :: i
        real(dp) :: adrift, tratio
        real(dp) :: v2, vavvt, vavvn
        real(dp) :: vav2sum, v2sum
        real(dp), dimension(3, nelec) :: v
  
        vav2sum = 0.d0
        v2sum = 0.d0
        do i=1,nelec
          v2    = v(1,i)**2 + v(2,i)**2 + v(3,i)**2
          vavvt = (dsqrt(1.d0+2.d0*adrift*v2*tau*tratio)-1.d0)/(adrift*v2)
          vavvn = vavvt/(tau*tratio)
  
          vav2sum = vav2sum + vavvn**2 * v2
          v2sum = v2sum + v2  
        enddo
  
    endsubroutine

end module branching
