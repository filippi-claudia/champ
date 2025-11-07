module jastrow_qmckl_mod
    
contains

subroutine jastrow_init_qmckl(ictx)

    use jastrow, only: norda,nordb,nordc,a4,b,c,scalek
    use jastrow4_mod, only: nterms4
    use precision_kinds, only: dp
    use scale_dist_mod, only: scale_dist2,switch_scale2
    use system,  only: iwctype,ncent,nelec,nctype
    use error,   only: fatal_error
    use contrl_file,    only: ounit
    use qmckl_data

    implicit none
    integer :: i
    integer, intent(in) :: ictx
    integer(qmckl_exit_code) :: rc
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
    
    if (nordc.gt.0) then
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
    use jastrow4_mod, only: jastrow_factor4
    use precision_kinds, only: dp
    use system,  only: ncent,nelec
    use error,   only: fatal_error
    use m_force_analytic, only: iforce_analy
    use da_jastrow, only: da_d2j, da_j, da_vj

    implicit none
    real(dp)  :: x(3,nelec)
    real(dp), dimension(3, nelec) :: fjo
    real(dp) :: fsumo, d2o
    real(dp) :: jen(1), jee(1), jeen(1)
    real(dp) :: jen_gl(4* nelec), jee_gl(4* nelec), jeen_gl(4* nelec)
    real(dp) :: da_j_en(3,ncent), da_j_een(3,ncent)
    real(dp) :: da_vj_en(3,nelec,3,ncent), da_vj_een(nelec,3,ncent,3)
    real(dp) :: da_d2j_en(3,ncent), da_d2j_een(ncent,3)
    integer(qmckl_exit_code) :: rc

    integer :: i, k, ic, l

    do i = 1, nelec
        fjo(1,i) = 0.0d0
        fjo(2,i) = 0.0d0
        fjo(3,i) = 0.0d0
    enddo
    d2o = 0.0d0

    rc = qmckl_set_point(qmckl_ctx(qmckl_no_ctx), 'N', 1_8*nelec, x, 3_8*nelec)
    if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow x-coords.')
    
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

    implicit none
    integer            :: iel, iflag,i
    real(dp)  :: x(3)
    real(dp), dimension(3, *) :: fjn
    real(dp) :: fsumn, d2n
    real(dp) :: jen(1), jee(1), jeen(1)
    real(dp) :: jen_gl(4,nelec), jee_gl(4,nelec), jeen_gl(4*nelec)
    integer(qmckl_exit_code) :: rc

    do i = 1, nelec
        fjn(1,i) = 0.0d0
        fjn(2,i) = 0.0d0
        fjn(3,i) = 0.0d0
    enddo
    d2n = 0.0d0

    ! iflag
    ! 0 = value only
    ! 1 = value, gradient and laplacian
    ! 2 = value and gradient 

 
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
            rc = qmckl_get_jastrow_champ_single_een_gl(qmckl_ctx(qmckl_no_ctx), jeen_gl, 1_8*4*nelec)
        else
            rc = qmckl_get_jastrow_champ_single_een_g(qmckl_ctx(qmckl_no_ctx), jeen_gl, 1_8*4*nelec)
        endif
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow gl.')
        do i = 1, nelec
            fjn(1,i) = fjn(1,i) +  jeen_gl(i+nelec*0)
            fjn(2,i) = fjn(2,i) +  jeen_gl(i+nelec*1)
            fjn(3,i) = fjn(3,i) +  jeen_gl(i+nelec*2)
        enddo

    endif


    end subroutine

end module
