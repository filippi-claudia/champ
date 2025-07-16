module deriv_jastrow_qmckl_mod
  contains
    !> Calculates the Jastrow factor AND its derivative.
    !> One shoul consider calculating the value and derivative in separate subroutines
    subroutine deriv_jastrow4_qmckl(x,fjo,d2o,fsumo,g,d2g,gvalue)

      use contrl_file, only: ounit
      use cuspmat4, only: iwc4, nterms
      use da_jastrow, only: da_d2j, da_j, da_vj
      use error, only: fatal_error
      use jastrow, only: norda, nordb, nordc
      use jaspointer, only: npoint
      use m_force_analytic, only: iforce_analy
      use optwf_nparmj, only: nparma ,nparmb, nparmc
      use optwf_parms, only: nparmj
      use optwf_wjas, only: iwjasa, iwjasb, iwjasc
      use precision_kinds, only: dp
      use qmckl
      use qmckl_data
      use system,  only: iwctype, nctype, nelec, ncent
      use vardep,  only: cdep, iwdepend, nvdepend
      implicit none

      real(dp), intent(in) , dimension(3, nelec) :: x
      real(dp), intent(out), dimension(3, nelec) :: fjo
      real(dp), intent(out) :: d2o
      real(dp), intent(out) :: fsumo
      real(dp), intent(out), dimension(nparmj) :: gvalue
      real(dp), intent(out), dimension(nparmj) :: d2g
      real(dp), intent(out), dimension(3, nelec, nparmj) :: g

      integer(qmckl_exit_code) :: rc
      
      integer :: i, ic, id, ideriv, iparm, it 
      integer :: j, jj, jparm
      integer :: k, l, l_hi, ll, m, n
   
      real(dp) :: cd

      real(dp), dimension(1) :: jen, jee, jeen
      real(dp), dimension(4 * nelec) :: jen_gl, jee_gl, jeen_gl

      real(dp), dimension(norda+1,nctype) :: deriv_a
      real(dp), dimension(nordb+1) :: deriv_b
      real(dp), dimension(nterms, nctype) :: deriv_c

      real(dp), dimension(4, nelec, norda+1, nctype) :: gl_deriv_a
      real(dp), dimension(4, nelec, nordb+1) :: gl_deriv_b
      real(dp), dimension(4, nelec, nterms, nctype) :: gl_deriv_c   

      real(dp) :: da_j_en(3,ncent), da_j_een(3,ncent)
      real(dp) :: da_vj_en(3,nelec,3,ncent), da_vj_een(nelec,3,ncent,3)
      real(dp) :: da_d2j_en(3,ncent), da_d2j_een(ncent,3)

      ! Initialize output arrys to 0
      do i = 1, nelec
        fjo(1,i) = 0.0d0
        fjo(2,i) = 0.0d0
        fjo(3,i) = 0.0d0
      enddo

      d2o = 0.0d0

      do iparm=1, nparmj
        gvalue(iparm) = 0.d0
        g(:,:,iparm) = 0.d0
        d2g(iparm) = 0.d0  
      enddo
    
      ! Set the electron coordinates
      rc = qmckl_set_point(qmckl_ctx(qmckl_no_ctx), 'N', 1_8*nelec, x, 3_8*nelec)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting QMCkl Jastrow x-coords.')

      ! Calculate the contributions to the Jastrow
      rc = qmckl_get_jastrow_champ_factor_ee(qmckl_ctx(qmckl_no_ctx), jee, 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e Jastrow.')

      rc = qmckl_get_jastrow_champ_factor_en(qmckl_ctx(qmckl_no_ctx), jen, 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow.')

      rc = qmckl_get_jastrow_champ_factor_een(qmckl_ctx(qmckl_no_ctx), jeen, 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow.')

      ! Calculate Jastrow gradients and laplacians
      rc = qmckl_get_jastrow_champ_factor_en_gl(qmckl_ctx(qmckl_no_ctx), jen_gl, 1_8*4*nelec)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow.')
      
      rc = qmckl_get_jastrow_champ_factor_ee_gl(qmckl_ctx(qmckl_no_ctx), jee_gl, 1_8*4*nelec)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e Jastrow.')

      rc = qmckl_get_jastrow_champ_factor_een_gl(qmckl_ctx(qmckl_no_ctx), jeen_gl, 1_8*4*nelec)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow.')

      ! Actual calculations of the derivatives to the paramaters
      rc = qmckl_get_jastrow_champ_factor_en_pderiv(qmckl_ctx(qmckl_no_ctx), deriv_a, nctype*(norda + 1) * 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting e-n Jastrow Parameter derivatives.')

      rc = qmckl_get_jastrow_champ_factor_ee_pderiv(qmckl_ctx(qmckl_no_ctx), deriv_b, (nordb + 1) * 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting e-e Jastrow Parameter derivatives.')
      
      if (nterms.gt.0) then
        rc = qmckl_get_jastrow_champ_factor_een_pderiv(qmckl_ctx(qmckl_no_ctx), deriv_c, nterms * nctype * 1_8)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting e-e-n Jastrow Parameter derivatives.')
      endif

      ! Calculation of the parameter derivatives of the gradient and Laplacian
      rc = qmckl_get_jastrow_champ_factor_en_gl_pderiv(qmckl_ctx(qmckl_no_ctx), gl_deriv_a, 4 * nelec * (norda + 1) * nctype * 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting e-n Jastrow Parameter derivatives of the gradient and Laplacian.')

      rc = qmckl_get_jastrow_champ_factor_ee_gl_pderiv(qmckl_ctx(qmckl_no_ctx), gl_deriv_b, 4 * nelec * (nordb + 1) * 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting e-e Jastrow Parameter derivatives of the gradient and Laplacian.')

      if (nterms .gt. 0) then
        rc = qmckl_get_jastrow_champ_factor_een_gl_pderiv(qmckl_ctx(qmckl_no_ctx), gl_deriv_c, 4 * nelec * nterms * nctype * 1_8)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting e-e-n Jastrow Parameter derivatives of the gradient and Laplacian.')
      endif 

      ! Collect Jastrow value
      fsumo = jee(1) + jen(1) + jeen(1)

      ! Collect Jastrow gradients and Laplacians
      do i = 1, nelec
        fjo(1,i) = fjo(1,i) + jen_gl(i+nelec*0) + jee_gl(i+nelec*0) + jeen_gl(i+nelec*0)
        fjo(2,i) = fjo(2,i) + jen_gl(i+nelec*1) + jee_gl(i+nelec*1) + jeen_gl(i+nelec*1)
        fjo(3,i) = fjo(3,i) + jen_gl(i+nelec*2) + jee_gl(i+nelec*2) + jeen_gl(i+nelec*2)
        d2o = d2o + jen_gl(i+nelec*3) + jee_gl(i+nelec*3) + jeen_gl(i+nelec*3)
      enddo

      iparm = 0

      ! Collect relevant a parameters
      do it = 1, nctype
        do i = 1, nparma(it)
          iparm = iparm + 1
          gvalue(iparm) = deriv_a(iwjasa(i,it), it)
          g(1,:,iparm) = gl_deriv_a(1,:,iwjasa(i,it),it)
          g(2,:,iparm) = gl_deriv_a(2,:,iwjasa(i,it),it)
          g(3,:,iparm) = gl_deriv_a(3,:,iwjasa(i,it),it)
          do j = 1, nelec
            d2g(iparm) = d2g(iparm) + gl_deriv_a(4,j,iwjasa(i,it),it)
          enddo
        enddo
      enddo
      
      ! Collect relevant b parameters
      do i = 1, nparmb(1)
        iparm = iparm + 1
        gvalue(iparm) = deriv_b(iwjasb(i,1))
        g(1,:,iparm) = gl_deriv_b(1,:,iwjasb(i,1))
        g(2,:,iparm) = gl_deriv_b(2,:,iwjasb(i,1))
        g(3,:,iparm) = gl_deriv_b(3,:,iwjasb(i,1))
        do j = 1, nelec
          d2g(iparm) = d2g(iparm) + gl_deriv_b(4,j,iwjasb(i,1))
        enddo
      enddo
      
      ! Resolving the dependencies in the c parameters due to e-n and e-e cusp conditions
      if(nordc.gt.1) then
        do it = 1, nctype  
          ll=0
          jj=1
          jparm=1

          ! Loop over all allowed combinations of indices k, l, n, m
          do n=2,nordc
            do k=n-1,0,-1
              if(k.eq.0) then
                l_hi=n-k-2
              else
                l_hi=n-k
              endif
              do l=l_hi,0,-1
                m=(n-k-l)/2
                if(2*m.eq.n-k-l) then
                  ll=ll+1
                  
                  ! Check if the current c is dependent or independent
                  ! ideriv = 0 -> Parameter is not varied and not dependent (Not sure if or when this happens)
                  ! ideriv = 1 -> Dependent variable
                  ! ideriv = 2 -> Independent variable
                  ideriv=0
                  if(ll.eq.iwjasc(jparm,it)) then
                    ideriv=2
                  else
                    do id=1,2*(nordc-1)
                      if(ll.eq.iwc4(id)) then
                        jj=id
                        if(nvdepend(jj,it).gt.0) ideriv=1
                      endif
                    enddo
                  endif
                  
                  if(ideriv.gt.0) then
                    
                    ! Dependent c parameters
                    if(ideriv.eq.1) then
                      do id=1,nvdepend(jj,it)
                        iparm=npoint(it)+iwdepend(jj,id,it)
                        cd=cdep(jj,id,it)

                        gvalue(iparm)=gvalue(iparm) + cd * deriv_c(ll,it)
                        g(1,:,iparm) = g(1,:,iparm) + cd * gl_deriv_c(1,:,ll,it)
                        g(2,:,iparm) = g(2,:,iparm) + cd * gl_deriv_c(2,:,ll,it)
                        g(3,:,iparm) = g(3,:,iparm) + cd * gl_deriv_c(3,:,ll,it)
                        do j = 1, nelec
                          d2g(iparm) = d2g(iparm) + cd * gl_deriv_c(4,j,ll,it)
                        enddo
                      enddo
                    
                    ! Independent c parameters
                    elseif(ideriv.eq.2) then
                      iparm=npoint(it)+jparm

                      gvalue(iparm)=gvalue(iparm) + deriv_c(ll,it)
                      g(1,:,iparm) = g(1,:,iparm) + gl_deriv_c(1,:,ll,it)
                      g(2,:,iparm) = g(2,:,iparm) + gl_deriv_c(2,:,ll,it)
                      g(3,:,iparm) = g(3,:,iparm) + gl_deriv_c(3,:,ll,it)
                      do j = 1, nelec
                        d2g(iparm) = d2g(iparm) + gl_deriv_c(4,j,ll,it)
                      enddo

                      jparm=jparm+1
                    endif
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
      endif

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
        do ic = 1, ncent
          do k = 1, 3
            da_j(k,1,1,ic) = da_j_en(k,ic) + da_j_een(k,ic)
            da_d2j(k,ic) = da_d2j_en(k,ic) + da_d2j_een(ic,k)
            do i = 1, nelec
              do l = 1, 3
                da_vj(k,l,i,ic) = da_vj_en(l,i,k,ic) + da_vj_een(i,l,ic,k)
              enddo
            enddo
          enddo
        enddo 

      endif
    return
    end

!--------------------------------------------------------------------------------------
    !> Calculates the single electron move one body Jastrow factor and its derivative
    subroutine deriv_jastrowe_qmckl(iel,x,fjn,d2n,fsumn,dpsij_ratio,iflag)
      use cuspmat4, only: iwc4, nterms
      use error, only: fatal_error
      use jaspointer, only: npoint
      use jastrow, only: norda, nordb, nordc
      use optwf_nparmj, only: nparma ,nparmb, nparmc
      use optwf_parms, only: nparmj
      use optwf_wjas, only: iwjasa, iwjasb, iwjasc
      use precision_kinds, only: dp
      use qmckl
      use qmckl_data
      use system,  only: nctype, nelec, ncent_tot
      use vardep,  only: cdep, iwdepend, nvdepend
      implicit none
      integer                          , intent(in) :: iel
      real(dp), dimension(3)           , intent(in) :: x
      real(dp), dimension(3, nelec)    , intent(out) :: fjn
      real(dp),                          intent(out) :: d2n
      real(dp),                          intent(out) :: fsumn
      real(dp), dimension(nparmj)      , intent(out) :: dpsij_ratio

      !> Sets if the Jastrow parameters should be computed \ 
      !> iflag = 0 -> Do not calculate derivatives         \ 
      !> iflag = 1 -> Calculate the gradient and the Laplacian
      !> iflag = 2 -> Only calculate the gradient          \ 
      integer,                           intent(in) :: iflag

      !for forces (not yet implemented)
      !real(dp), dimension(3)           , intent(out) :: vjn
      !real(dp), dimension(3, ncent_tot), intent(out) :: da_psij_ratio

      integer :: i, id, ideriv, iparm, it 
      integer :: j, jj, jparm
      integer :: k, l, l_hi, ll, m, n
      integer :: rc
      
      real(dp) :: cd
      
      real(dp) :: jee(1), jen(1), jeen(1)
      real(dp), dimension(4,nelec) :: jee_gl, jen_gl
      real(dp), dimension(4*nelec) :: jeen_gl

      real(dp), dimension(norda+1,nctype) :: sderiv_a
      real(dp), dimension(nordb+1) :: sderiv_b
      real(dp), dimension(nterms, nctype) :: sderiv_c

      ! Set arrays to zero
      do i = 1, nelec
        fjn(1,i) = 0.0d0
        fjn(2,i) = 0.0d0
        fjn(3,i) = 0.0d0
      enddo
      d2n = 0.0d0
      dpsij_ratio = 0.d0

      ! Set the position after the single electron move
      rc = qmckl_set_single_point(qmckl_ctx(qmckl_no_ctx), 'N', iel*1_8, x, 3_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error setting single point.')

      ! Calcculate the single electron Jastrow value
      rc = qmckl_get_jastrow_champ_single_en(qmckl_ctx(qmckl_no_ctx), jen, 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow.')

      rc = qmckl_get_jastrow_champ_single_ee(qmckl_ctx(qmckl_no_ctx), jee, 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e Jastrow.')

      rc = qmckl_get_jastrow_champ_single_een(qmckl_ctx(qmckl_no_ctx), jeen, 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow.')

      ! Actual calculations of the derivatives to the paramaters
      rc = qmckl_get_jastrow_champ_single_en_pderiv(qmckl_ctx(qmckl_no_ctx), sderiv_a, (norda + 1) * nctype * 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow single move parameter derivative.')

      rc = qmckl_get_jastrow_champ_single_ee_pderiv(qmckl_ctx(qmckl_no_ctx), sderiv_b, (nordb + 1) * 1_8)
      if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e Jastrow single move parameter derivative.')

      if (nterms .gt. 0) then
        rc = qmckl_get_jastrow_champ_single_een_pderiv(qmckl_ctx(qmckl_no_ctx), sderiv_c, nterms * nctype * 1_8)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow single move parameter derivative.')
      endif

      iparm = 0
      
      ! Collect Jastrow value
      fsumn = jen(1) + jee(1) + jeen(1)

      ! Collect relevant a parameters
      do it = 1, nctype
        do i = 1, nparma(it)
          iparm = iparm + 1
          dpsij_ratio(iparm) = sderiv_a(iwjasa(i,it), it)
        enddo
      enddo

      ! Collect relevant b parameters
      do i = 1, nparmb(1)
        iparm = iparm + 1
        dpsij_ratio(iparm) = sderiv_b(iwjasb(i,1))
      enddo

      ! Resolving the dependencies in the c parameters due to e-n and e-e cusp conditions
      if(nordc.gt.1) then
        do it = 1, nctype  
          ll=0
          jj=1
          jparm=1
        
          ! Loop over all allowed combinations of indices k, l, n, m
          do n=2,nordc
            do k=n-1,0,-1
              if (k.eq.0) then
                l_hi=n-k-2
              else
                l_hi=n-k
              endif
              do l=l_hi,0,-1
                m=(n-k-l)/2
                if (2*m.eq.n-k-l) then
                  ll=ll+1
                  
                  ! Check if the current c is dependent or independent
                  ! ideriv = 0 -> Parameter is not varied and not dependent (Not sure if or when this happens)
                  ! ideriv = 1 -> Dependent variable
                  ! ideriv = 2 -> Independent variable
                  ideriv=0
                  if(ll.eq.iwjasc(jparm,it)) then
                    ideriv=2
                  else
                    do id=1,2*(nordc-1)
                      if(ll.eq.iwc4(id)) then
                        jj=id
                        if(nvdepend(jj,it).gt.0) ideriv=1
                      endif
                    enddo
                  endif
                  
                  if(ideriv.gt.0) then
                    
                    ! Dependent c parameters
                    if(ideriv.eq.1) then
                      do id=1,nvdepend(jj,it)
                        iparm=npoint(it)+iwdepend(jj,id,it)
                        cd=cdep(jj,id,it)
                        dpsij_ratio(iparm)=dpsij_ratio(iparm) + cd * sderiv_c(ll,it)
                      enddo
                    
                    ! Independent c parameters
                    elseif(ideriv.eq.2) then
                      iparm=npoint(it)+jparm
                      dpsij_ratio(iparm)=dpsij_ratio(iparm) + sderiv_c(ll,it)
                      jparm=jparm+1
                    endif
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
      endif
      
      ! Calculate derivatives (if necessary)
      if (iflag .ge. 1) then

        rc = qmckl_get_jastrow_champ_single_ee_gl(qmckl_ctx(qmckl_no_ctx), jee_gl, 1_8*4*nelec)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e Jastrow gl.')
        
        rc = qmckl_get_jastrow_champ_single_en_gl(qmckl_ctx(qmckl_no_ctx), jen_gl, 1_8*4*nelec)
        if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-n Jastrow gl.')
        
        if (iflag .eq. 1) then
          rc = qmckl_get_jastrow_champ_single_een_gl(qmckl_ctx(qmckl_no_ctx), jeen_gl, 1_8*4*nelec)
          if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow gl.')

          do i = 1, nelec
            d2n = d2n + jeen_gl(i+nelec*3)
          enddo
        elseif (iflag .eq. 2) then
          rc = qmckl_get_jastrow_champ_single_een_g(qmckl_ctx(qmckl_no_ctx), jeen_gl, 1_8*4*nelec)
          if (rc /= QMCKL_SUCCESS) call fatal_error('Error getting QMCkl e-e-n Jastrow g.')
        else 
          call fatal_error('Error in deriv_jastrowe_qmckl: iflag greater than 2 does not existe')
        endif
        
        do i = 1, nelec
          fjn(1,i) = fjn(1,i) + jee_gl(1,i) + jen_gl(1,i) + jeen_gl(i+nelec*0)
          fjn(2,i) = fjn(2,i) + jee_gl(2,i) + jen_gl(2,i) + jeen_gl(i+nelec*1)
          fjn(3,i) = fjn(3,i) + jee_gl(3,i) + jen_gl(3,i) + jeen_gl(i+nelec*2)
        enddo

      endif
    end 

end module

