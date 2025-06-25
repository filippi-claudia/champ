module jastrow_qmckl_mod
    
contains

subroutine jastrow_init_qmckl(ictx)

    use jastrow, only: norda,nordb,nordc,a4,b,c,scalek
    use jastrow4_mod, only: nterms4
    use precision_kinds, only: dp
    use scale_dist_mod, only: scale_dist2,switch_scale2
    use system,  only: iwctype,ncent,nelec,nup,nctype
    use error,   only: fatal_error
    use contrl_file,    only: ounit


    use qmckl_data

    implicit none
    integer :: i
    integer, intent(in) :: ictx
    integer(qmckl_exit_code) :: rc, rc2
    double precision :: scalek_en(nctype)
    integer*8 :: itypes(ncent)
    integer*8 :: dimc, norda_l, nordb_l
    real(dp) :: x(3, nelec)

    do i=1,nelec
        x(1,i) = 0.0d0
        x(2,i) = 0.0d0
        x(3,i) = 0.0d0
    end do

    rc = qmckl_set_electron_coord(qmckl_ctx(ictx), 'N', 1_8, x, nelec*3_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow x-coords.')

    norda_l = max(norda, 1)
    nordb_l = max(nordb, 1)

    rc = qmckl_set_jastrow_champ_aord_num (qmckl_ctx(ictx), norda_l*1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow aord num.')

    rc = qmckl_set_jastrow_champ_bord_num (qmckl_ctx(ictx), nordb_l*1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow bord num.')

    rc = qmckl_set_jastrow_champ_cord_num (qmckl_ctx(ictx), nordc*1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow cord num.')

    rc = qmckl_set_jastrow_champ_type_nucl_num (qmckl_ctx(ictx), nctype*1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow nucl num.')

    do i = 1, ncent
        itypes(i) = int(iwctype(i),8)-1
    end do
    rc = qmckl_set_jastrow_champ_type_nucl_vector (qmckl_ctx(ictx), itypes, ncent*1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow nucl vec.')

    rc = qmckl_set_jastrow_champ_a_vector (qmckl_ctx(ictx), a4(1:norda_l+1,1:nctype,1), (norda_l+1)*nctype)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow a vec.')

    rc = qmckl_set_jastrow_champ_b_vector (qmckl_ctx(ictx), b(1:(nordb_l+1),1,1), (nordb_l+1)*2_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow b vec.')
    
    if (nordc.gt.1) then
        dimc = nterms4(nordc)
        rc = qmckl_set_jastrow_champ_c_vector (qmckl_ctx(ictx), c(1:dimc,1:nctype,1), dimc*nctype*1_8)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow c vec.')
    endif

    rc = qmckl_set_jastrow_champ_rescale_factor_ee (qmckl_ctx(ictx), scalek(1))
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow rescale ee.')

    scalek_en(:) = scalek(1)
    rc = qmckl_set_jastrow_champ_rescale_factor_en (qmckl_ctx(ictx),scalek_en, nctype*1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow rescale en.')

    rc = qmckl_set_jastrow_champ_spin_independent(qmckl_ctx(ictx), 0)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow spin.')

    rc = qmckl_set_numprec_precision(qmckl_ctx(ictx), 53) ! 24
    if (rc .ne. QMCKL_SUCCESS) call fatal_error('INPUT: QMCkl error: Unable to set precision')

    write(ounit, *) " QMCkl precision set to 53 bits"
end subroutine



!---------------------------------------------------------------------------------

    subroutine jastrow_qmckl(x,fjo,d2o,fsumo)

    use qmckl_data
    use multiple_geo, only: iwf
    use jastrow4_mod, only: jastrow_factor4
    use precision_kinds, only: dp
    use jastrow, only: asymp_jasa,asymp_jasb
    use system,  only: iwctype,ncent,nelec,nup,nctype
    use error,   only: fatal_error
    use m_force_analytic, only: iforce_analy
    use contrl_file, only: ounit
    use da_jastrow, only: da_d2j, da_j, da_vj

    implicit none
    real(dp)  :: x(3,*)
    real(dp), dimension(3, *) :: fjo
    real(dp) :: old_x(3,nelec)
    real(dp) :: fsumo, d2o
    real(dp) :: value

    real(dp) :: jen(1), jee(1), jeen(1)
    real(dp) :: jen_gl(4* nelec), jee_gl(4* nelec), jeen_gl(4* nelec)

    real(dp) :: da_j_en(3,ncent), da_j_een(3,ncent)
    real(dp) :: da_vj_en(3,nelec,3,ncent), da_vj_een(nelec,3,ncent,3)
    real(dp) :: da_d2j_en(3,ncent), da_d2j_een(ncent,3)
    integer(qmckl_exit_code) :: rc

    integer :: it,i ,j, same, k, ic, l
    double precision :: xx(10000)

    do i = 1, nelec
        fjo(1,i) = 0.0d0
        fjo(2,i) = 0.0d0
        fjo(3,i) = 0.0d0
    enddo
    d2o = 0.0d0

    rc = qmckl_get_point(qmckl_ctx(qmckl_no_ctx), 'N', old_x, 3_8*nelec)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl Jastrow x-coords.')
   
    same = 1
    do i = 1, nelec
        do k = 1, 3
            if (abs(old_x(k,i) - x(k,i)) > 1e-6) then
                same = 0
            end if
        end do
        !if (same == 0) exit
    end do

    ! if (same == 0) then
        rc = qmckl_set_point(qmckl_ctx(qmckl_no_ctx), 'N', 1_8*nelec, x, 3_8*nelec)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow x-coords.')
    ! endif
    
    rc = qmckl_get_jastrow_champ_factor_ee(qmckl_ctx(qmckl_no_ctx), jee, 1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e Jastrow.')

    rc = qmckl_get_jastrow_champ_factor_en(qmckl_ctx(qmckl_no_ctx), jen, 1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow.')

    rc = qmckl_get_jastrow_champ_factor_een(qmckl_ctx(qmckl_no_ctx), jeen, 1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow.')

    fsumo = jee(1) + jen(1) + jeen(1)

    rc = qmckl_get_jastrow_champ_factor_en_gl(qmckl_ctx(qmckl_no_ctx), jen_gl, 1_8*4*nelec)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow.')
    do i = 1, nelec
        fjo(1,i) = fjo(1,i) + jen_gl(i+nelec*0)
        fjo(2,i) = fjo(2,i) + jen_gl(i+nelec*1)
        fjo(3,i) = fjo(3,i) + jen_gl(i+nelec*2)
        d2o = d2o + jen_gl(i+nelec*3)
    enddo

    rc = qmckl_get_jastrow_champ_factor_ee_gl(qmckl_ctx(qmckl_no_ctx), jee_gl, 1_8*4*nelec)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e Jastrow.')
    do i = 1, nelec
        fjo(1,i) = fjo(1,i) +  jee_gl(i+nelec*0)
        fjo(2,i) = fjo(2,i) +  jee_gl(i+nelec*1)
        fjo(3,i) = fjo(3,i) +  jee_gl(i+nelec*2)
        d2o = d2o + jee_gl(i+nelec*3)
    enddo

   rc = qmckl_get_jastrow_champ_factor_een_gl(qmckl_ctx(qmckl_no_ctx), jeen_gl, 1_8*4*nelec)
   if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow.')
   do i = 1, nelec
       fjo(1,i) = fjo(1,i) +  jeen_gl(i+nelec*0)
       fjo(2,i) = fjo(2,i) +  jeen_gl(i+nelec*1)
       fjo(3,i) = fjo(3,i) +  jeen_gl(i+nelec*2)
       d2o = d2o + jeen_gl(i+nelec*3)
   enddo


   if (iforce_analy.eq.1) then

      rc = qmckl_get_forces_jastrow_en(qmckl_ctx(qmckl_no_ctx), da_j_en, ncent*3_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow forces.')

      rc = qmckl_get_forces_jastrow_een(qmckl_ctx(qmckl_no_ctx), da_j_een, ncent*3_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow forces.')

      rc = qmckl_get_forces_jastrow_en_g(qmckl_ctx(qmckl_no_ctx), da_vj_en, 3_8*nelec*ncent*3_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow gradient forces.')

      rc = qmckl_get_forces_jastrow_een_g(qmckl_ctx(qmckl_no_ctx), da_vj_een, 3_8*nelec*ncent*3_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow gradient forces.')

      rc = qmckl_get_forces_jastrow_en_l(qmckl_ctx(qmckl_no_ctx), da_d2j_en, 3_8*ncent)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow Laplacian forces.')

      rc = qmckl_get_forces_jastrow_een_l(qmckl_ctx(qmckl_no_ctx), da_d2j_een, 3_8*ncent)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow Laplacian forces.')

      da_j = 0.d0
      do ic=1,ncent
        do k=1,3
            da_j(k,1,1,ic) = da_j_en(k,ic)+da_j_een(k,ic)
            da_d2j(k,ic) = da_d2j_en(k,ic)+da_d2j_een(ic,k)
            do i=1,nelec
                do l=1,3
                    da_vj(k,l,i,ic)=da_vj_en(l,i,k,ic)+da_vj_een(i,l,ic,k)
                enddo
            enddo
        enddo
      enddo 

   endif

    end subroutine

!---------------------------------------------------------------------------------

    subroutine jastrowe_qmckl(iel,x,fjn,d2n,fsumn,iflag)

    use qmckl_data
    use jastrow4e_mod, only: jastrow4e
    use precision_kinds, only: dp
    use system,  only: nelec
    use error,   only: fatal_error
    use jastrow_update, only: d2ijo, d2o, fijo, fjo, fso, fsumo
    use contrl_file, only: ounit

    implicit none
    integer            :: iel, iflag,i
    real(dp)  :: x(3)
    real(dp), dimension(3, *) :: fjn
    real(dp) :: fsumn, d2n
    real(dp) :: value

    real(dp) :: jen(1), jee(1), jeen(1)
    real(dp) :: jen_gl(4,nelec), jee_gl(4,nelec), jeen_gl(4*nelec)

    integer(qmckl_exit_code) :: rc
    character*(1024) :: err_message = ''

    do i = 1, nelec
        fjn(1,i) = 0.0d0
        fjn(2,i) = 0.0d0
        fjn(3,i) = 0.0d0
    enddo
    d2n = 0.0d0

    ! iflag
    ! 0 = value only
    ! 1 = value and gradient 
    ! 2 = value, gradient and laplacian ! not implemented

 
    rc = qmckl_set_single_point(qmckl_ctx(qmckl_no_ctx), 'N', iel*1_8, x, 3_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting single point.')

    rc = qmckl_get_jastrow_champ_single_en(qmckl_ctx(qmckl_no_ctx), jen, 1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow.')

    rc = qmckl_get_jastrow_champ_single_ee(qmckl_ctx(qmckl_no_ctx), jee, 1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e Jastrow.')

    rc = qmckl_get_jastrow_champ_single_een(qmckl_ctx(qmckl_no_ctx), jeen, 1_8)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow.')

    fsumn =  jen(1)  + jee(1) + jeen(1)

    if (iflag .ge. 1) then

        rc = qmckl_get_jastrow_champ_single_ee_gl(qmckl_ctx(qmckl_no_ctx), jee_gl, 1_8*4*nelec)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e Jastrow gl.')
        do i = 1, nelec
            fjn(1,i) = fjn(1,i) +  jee_gl(1,i)
            fjn(2,i) = fjn(2,i) +  jee_gl(2,i)
            fjn(3,i) = fjn(3,i) +  jee_gl(3,i)
        enddo

        rc = qmckl_get_jastrow_champ_single_en_gl(qmckl_ctx(qmckl_no_ctx), jen_gl, 1_8*4*nelec)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow gl.')
        do i = 1, nelec
            fjn(1,i) = fjn(1,i) +  jen_gl(1,i)
            fjn(2,i) = fjn(2,i) +  jen_gl(2,i)
            fjn(3,i) = fjn(3,i) +  jen_gl(3,i)
        enddo
        
        if (iflag.eq.1) then
            rc = qmckl_get_jastrow_champ_single_een_g(qmckl_ctx(qmckl_no_ctx), jeen_gl, 1_8*4*nelec)
            if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow gl.')
        else
            rc = qmckl_get_jastrow_champ_single_een_gl(qmckl_ctx(qmckl_no_ctx), jeen_gl, 1_8*4*nelec)
            if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow gl.')

            do i = 1, nelec
                d2n = d2n + jee_gl(4,i) + jen_gl(4,i) + jeen_gl(i+3*nelec)
            end do
        endif
        
        do i = 1, nelec
            fjn(1,i) = fjn(1,i) +  jeen_gl(i+nelec*0)
            fjn(2,i) = fjn(2,i) +  jeen_gl(i+nelec*1)
            fjn(3,i) = fjn(3,i) +  jeen_gl(i+nelec*2)
        enddo

    endif


    end subroutine

!    subroutine jastrow_quad_qmckl(nquad,indices,x,ratio_jn,vjn,da_psij_ratio,iflag)
!
!        use qmckl_data
!        use jastrow4e_mod, only: jastrow4e
!        use precision_kinds, only: dp
!        use system,  only: nelec, ncent, ncent_tot
!        use error,   only: fatal_error
!        use jastrow_update, only: d2ijo, d2o, fijo, fjo, fso, fsumo
!        use contrl_file, only: ounit
!    
!        implicit none
!        integer            :: iflag,i, nq, nquad, indx,k, ic
!        integer*8 :: indices(nquad)
!        real(dp)  :: x(3,nquad)
!        real(dp), dimension(3, nquad) :: vjn
!        real(dp), dimension(3,ncent_tot,*) :: da_psij_ratio
!        real(dp), dimension(nquad) :: ratio_jn
!        real(dp) :: value
!    
!        real(dp) :: jen(nquad), jee(nquad), jeen(nquad)
!        real(dp) :: jen_gl(3,nelec,nquad), jee_gl(3,nelec,nquad), jeen_gl(3,nquad)
!        real(dp), dimension(3,ncent,nquad) :: da_single_een, da_single_en
!    
!        integer(qmckl_exit_code) :: rc
!        character*(1024) :: err_message = ''
!    
!
!        ! iflag
!        ! 0 = value only
!        ! 1 = value and gradient 
!
!        do nq = 1, nquad
!            indices(nq) = indices(nq)-1
!        end do
!     
!        rc = qmckl_set_quad_points(qmckl_ctx(qmckl_no_ctx), 'N', nquad, indices, x, 3_8*nquad)
!        if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting single point.')
!    
!        rc = qmckl_get_jastrow_champ_quad_en(qmckl_ctx(qmckl_no_ctx), jen, 1_8*nquad)
!        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow.')
!    
!        rc = qmckl_get_jastrow_champ_quad_ee(qmckl_ctx(qmckl_no_ctx), jee, 1_8*nquad)
!        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e Jastrow.')
!    
!        rc = qmckl_get_jastrow_champ_quad_een(qmckl_ctx(qmckl_no_ctx), jeen, 1_8*nquad)
!        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow.')
!    
!        do nq = 1, nquad
!            ratio_jn(nq) =  jen(nq)  + jee(nq) + jeen(nq)
!        enddo
!    
!        if (iflag .ge. 1) then
!    
!        rc = qmckl_get_jastrow_champ_quad_ee_gl(qmckl_ctx(qmckl_no_ctx), jee_gl, 1_8*3*nelec*nquad)
!        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e Jastrow gl.')
!    
!        rc = qmckl_get_jastrow_champ_quad_en_gl(qmckl_ctx(qmckl_no_ctx), jen_gl,  1_8*3*nelec*nquad)
!        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow gl.')
!    
!        rc = qmckl_get_jastrow_champ_quad_een_gl(qmckl_ctx(qmckl_no_ctx), jeen_gl,  1_8*3*nquad)
!        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow gl.')
!        do nq = 1, nquad
!            !fjn(1,i,nq) = jeen_gl(i+nelec*0,nq) + jen_gl(1,i,nq) + jee_gl(1,i,nq)
!            !fjn(2,i,nq) = jeen_gl(i+nelec*1,nq) + jen_gl(2,i,nq) + jee_gl(2,i,nq)
!            !fjn(3,i,nq) = jeen_gl(i+nelec*2,nq) + jen_gl(3,i,nq) + jee_gl(3,i,nq)
!            indx = indices(nq)+1
!            !print *, 'jeen_gl', jeen_gl(1,1,nq), jeen_gl(2,1,nq), jeen_gl(3,1,nq)
!            vjn(1,nq) = jeen_gl(1,nq)  + jen_gl(1,indx,nq) + jee_gl(1,indx,nq) + fjo(1,indx,1)
!            vjn(2,nq) = jeen_gl(2,nq)  + jen_gl(2,indx,nq) + jee_gl(2,indx,nq) + fjo(2,indx,1)
!            vjn(3,nq) = jeen_gl(3,nq)  + jen_gl(3,indx,nq) + jee_gl(3,indx,nq) + fjo(3,indx,1)
!            ! vjn(1,nq) = jeen_gl(indx,1,nq) 
!            ! vjn(2,nq) = jeen_gl(indx,2,nq) 
!            ! vjn(3,nq) = jeen_gl(indx,3,nq)
!      
!        end do
!
!        rc = qmckl_get_forces_jastrow_quad_en(qmckl_ctx(qmckl_no_ctx), da_single_en, 3*ncent*nquad)
!        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl Jastrow single en force.')
!        rc = qmckl_get_forces_jastrow_quad_een(qmckl_ctx(qmckl_no_ctx), da_single_een, 3*ncent*nquad)
!        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl Jastrow single een force.')
!        do nq = 1, nquad
!            do ic=1,ncent
!                do k=1,3
!                    da_psij_ratio(k,ic,nq)=da_single_en(k,ic,nq)+da_single_een(k,ic,nq)
!                enddo
!            enddo
!        enddo
!    
!       endif
!
!    
!        end subroutine

    end module
