module acuest_mod
! Written by Cyrus Umrigar, modified by Claudia Filippi
! routine to accumulate estimators for energy etc.

      use acuest_reduce_mod, only: acues1_reduce,acuest_reduce
      use config,  only: eold,nearesto,psi2o,psido,psijo,rmino,rvmino
      use config,  only: vold,xold
      use contrl_file, only: ounit
      use control, only: ipr, mode
      use csfs,    only: nstates
      use determinant_psig_mod, only: determinant_psig
      use determinante_mod, only: compute_determinante_grad
      use distance_mod, only: r_en,rvec_en
      use distances_mod, only: distances
      use est2cm,  only: ecm2,ecm21,pecm2,tpbcm2
      use estcum,  only: ecum,ecum1,iblk,pecum,tpbcum
      use estpsi,  only: apsi,aref,detref
      use estsig,  only: ecm21s,ecum1s
      use estsum,  only: acc,esum,esum1,pesum,tpbsum
      use ewald, only: cos_n_sum,sin_n_sum
      use force_analytic, only: force_analy_cum,force_analy_init
      use force_analytic, only: force_analy_save
      use forcewt, only: wcum,wsum
      use hpsi_mod, only: hpsi
      use inputflags, only: eps_node_cutoff,node_cutoff
      use mmpol,   only: mmpol_init
      use mmpol_vmc, only: mmpol_cum,mmpol_save
      use mpi
      use mpiconf, only: nproc,wid
      use mstates_ctrl, only: iguiding
      use mstates_mod, only: MSTATES
      use multiple_geo, only: MFORCE,fcm2,fcum,nforce,pecent
      use multiple_states, only: efficiency_init
      use multislater, only: detiab
      use nodes_distance_mod, only: nodes_distance,rnorm_nodes_num
      use optci_mod, only: optci_cum,optci_init,optci_save
      use optjas_mod, only: optjas_cum,optjas_init,optjas_save
      use optorb_cblock, only: ns_current
      use optorb_f_mod, only: optorb_cum,optorb_init,optorb_save
      use optwf_control, only: ioptorb
      use optx_jas_ci, only: optx_jas_ci_init
      use optx_jas_orb, only: optx_jas_orb_init
      use optx_orb_ci, only: optx_orb_ci_init
      use pcm_mod, only: pcm_init
      use pcm_vmc, only: pcm_cum,pcm_save
      use pot,     only: pot_nn
      use precision_kinds, only: dp
      use prop_vmc, only: prop_save
      use properties_mod, only: prop_cum,prop_init
      use pseudo,  only: nloc
      use qua,     only: nquad,wq,xq,yq,zq
      use rotqua_mod, only: gesqua
      use slater,  only: kref
      use strech_mod, only: strech
      use system,  only: cent,iwctype,ncent,nelec,znuc,nup
      use vmc_mod, only: nwftypeorb, stoj
      use pathak_mod, only: init_eps_pathak, ipathak

      implicit none

contains
      subroutine acuest

      implicit none

      integer :: i, ic, ierr, ifr, istate, jel, k
      real(dp) :: ajacob, distance_node, eave
      real(dp) :: psidg, rnorm_nodes
      real(dp) :: penow, tpbnow
      real(dp), dimension(3,nelec) :: xstrech

      real(dp), dimension(:,:), allocatable :: enow

      allocate(enow(MSTATES, MFORCE))

! xsum = sum of values of x from metrop
! xnow = average of values of x from metrop
! xcum = accumulated sums of xnow
! xcm2 = accumulated sums of xnow**2

! collect cumulative averages
      do ifr=1,nforce
        do istate=1,nstates

        enow(istate,ifr)=esum(istate,ifr)/wsum(istate,ifr)
        wcum(istate,ifr)=wcum(istate,ifr)+wsum(istate,ifr)
        ecum(istate,ifr)=ecum(istate,ifr)+esum(istate,ifr)
        ecm2(istate,ifr)=ecm2(istate,ifr)+esum(istate,ifr)*enow(istate,ifr)

        if(ifr.eq.1) then
          penow=pesum(istate)/wsum(istate,ifr)
          tpbnow=tpbsum(istate)/wsum(istate,ifr)

          pecm2(istate)=pecm2(istate)+pesum(istate)*penow
          tpbcm2(istate)=tpbcm2(istate)+tpbsum(istate)*tpbnow

          pecum(istate)=pecum(istate)+pesum(istate)
          tpbcum(istate)=tpbcum(istate)+tpbsum(istate)

         else
          fcum(istate,ifr)=fcum(istate,ifr)+wsum(istate,1)*(enow(istate,ifr)-esum(istate,1)/wsum(istate,1))
          fcm2(istate,ifr)=fcm2(istate,ifr)+wsum(istate,1)*(enow(istate,ifr)-esum(istate,1)/wsum(istate,1))**2
        endif
        enddo
      enddo

! only called for ifr=1
      call optjas_cum(wsum(1,1),enow(1,1))
      call optorb_cum(wsum(1,1),esum(1,1))
      call optci_cum(wsum(1,1))
      call prop_cum(wsum(1,1))
      call pcm_cum(wsum(1,1))
      call mmpol_cum(wsum(1,1))

      if(wid) eave=ecum(1,1)/wcum(1,1)
      call MPI_BCAST(eave,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      call force_analy_cum(wsum(1,1),eave)

! zero out xsum variables for metrop

      do istate=1,nstates
        do ifr=1,nforce
          esum(istate,ifr)=0
          wsum(istate,ifr)=0
        enddo
        pesum(istate)=0
        tpbsum(istate)=0
      enddo

      call prop_init(1)
      call optorb_init(1)
      call optci_init(1)
      call pcm_init(1)
      call mmpol_init(1)
      call force_analy_init(1)

      call acuest_reduce(enow)

      if(allocated(enow)) deallocate(enow)

      end subroutine
!-----------------------------------------------------------------------
      subroutine acues1(wtg)
      implicit none
      real(dp), dimension(*) :: wtg
      integer :: istate, k
! statistical fluctuations without blocking
      do istate=1,nstates
        ecum1(istate)=ecum1(istate)+esum1(istate)*wtg(istate)
        ecm21(istate)=ecm21(istate)+esum1(istate)**2*wtg(istate)
        esum1(istate)=0

        apsi(istate)=apsi(istate)+dabs(psido(istate))
      enddo
      do k=1,nwftypeorb

        aref(k)=aref(k)+dabs(detiab(kref,1,k)*detiab(kref,2,k))

        detref(1,k)=detref(1,k)+dlog10(dabs(detiab(kref,1,k)))
        detref(2,k)=detref(2,k)+dlog10(dabs(detiab(kref,2,k)))
      enddo

      call acues1_reduce

      end subroutine
!-----------------------------------------------------------------------
      subroutine acusig(wtg)
      implicit none
      real(dp), dimension(*) :: wtg
      integer :: istate
! sigma evaluation
      do istate=1,nstates
        ecum1s(istate)=ecum1s(istate)+esum1(istate)*wtg(istate)
        ecm21s(istate)=ecm21s(istate)+esum1(istate)**2*wtg(istate)
        esum1(istate)=0
      enddo
      end subroutine
!-----------------------------------------------------------------------
      subroutine zerest
      use slater,  only: d2dx2, ddx
      use backflow_mod, only: single_rios_backflow, backflow
      use m_backflow, only: quasi_x, dquasi_dx, d2quasi_dx2
      use system, only: nelec
      use hpsie,   only: psie
      implicit none
      integer :: i, ic, ifr, istate, jel, k, iel, j, kk, l
      real(dp) :: psidg, rnorm_nodes, tmp, tmp2
      real(dp) :: ajacob, distance_node
      real(dp), dimension(3,nelec) :: xstrech
      real(dp), dimension(3) :: ddx_l
      real(dp), dimension(MSTATES) :: ekino, psidoo, psidoo2
      real(dp), dimension(3,nelec) :: quasi_x_new
      real(dp), dimension(3,nelec,3,nelec) :: dquasi_dx_new
      real(dp), dimension(3,nelec,nelec) :: d2quasi_dx2_new
      integer, dimension(nelec) :: indices
      real(dp), dimension(3) :: xnew
      real(dp), dimension(3,nelec) :: xnew2
      real(dp), dimension(MSTATES) :: psidnn, psidn
      real(dp) :: delta, dist, val_j_1, val_d_1, val_j_2, val_d_2, div_j, div_d, epsilon
      real(dp), dimension(3) :: rvec
      real(dp), dimension(3,nelec) :: x_orig

! entry point to zero out all averages etc.
! the initial values of energy psi etc. is also calculated here
! although that really only needs to be done before the equil. blocks.

      iblk=0

! set quadrature points
      if(nloc.gt.0) call gesqua(nquad,xq,yq,zq,wq)

! zero out estimators
      acc=0
      do istate=1,nstates
        pecum(istate)=0
        tpbcum(istate)=0
        ecum1(istate)=0
        ecum1s(istate)=0

        pecm2(istate)=0
        tpbcm2(istate)=0
        ecm21(istate)=0
        ecm21s(istate)=0

        pesum(istate)=0
        tpbsum(istate)=0

        apsi(istate)=0
      enddo
      do k=1,nwftypeorb
        detref(1,k)=0
        detref(2,k)=0

        aref(k)=0
      enddo

      call optjas_init
      call optci_init(0)
      call optorb_init(0)
      call optx_jas_orb_init
      call optx_jas_ci_init
      call optx_orb_ci_init

      call prop_init(0)
      call pcm_init(0)
      call mmpol_init(0)
      call force_analy_init(0)
      call efficiency_init
      if(ipathak.gt.0) call init_eps_pathak()

      do ifr=1,nforce
        do istate=1,nstates
          ecum(istate,ifr)=0
          ecm2(istate,ifr)=0
          wcum(istate,ifr)=0
          esum(istate,ifr)=0
          wsum(istate,ifr)=0
          fcum(istate,ifr)=0
          fcm2(istate,ifr)=0
        enddo
      enddo

! get nuclear potential energy
      call pot_nn(cent,znuc,iwctype,ncent,pecent,cos_n_sum,sin_n_sum)

! get wavefunction etc. at initial point

! secondary configs
! set n- and e-coords and n-n potentials before getting wavefn. etc.
      do ifr=2,nforce
        call strech(xold,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psido,psijo,ekino,eold(1,ifr),0,ifr)
        do istate=1,nstates
          psi2o(istate,ifr)=2*(dlog(dabs(psido(istate)))+psijo(stoj(istate)))+dlog(ajacob)
        enddo
      enddo

! primary config
! set n- and e-coords and n-n potentials before getting wavefn. etc.
      if(nforce.gt.1) call strech(xold,xstrech,ajacob,1,0)
      call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)

      ! Cusp check block
      if (.false.) then
         write(ounit,*) '--- Cusp Check (Spherical Average) ---'
         delta = 0.00001d0
         epsilon = 0.00000d0 ! Avoid exactly 0
         
         x_orig = xold
         
         ! Electron-Nucleus
         do ic = 1, ncent
            do i = 1, nelec
               div_j = 0.0d0
               div_d = 0.0d0
               
               ! Integrate over 6 directions (+x,-x,+y,-y,+z,-z)
               do k = 1, 3
                  do l = -1, 1, 2 ! -1 and +1
                     rvec = 0.0d0
                     rvec(k) = dble(l)
                     
                     ! 1. Point at epsilon
                     xold(:,i) = cent(:,ic) + epsilon * rvec
                     call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)
                     val_j_1 = psijo(1)
                     val_d_1 = psido(1)
                     
                     ! 2. Point at epsilon + delta
                     xold(:,i) = cent(:,ic) + (epsilon + delta) * rvec
                     call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)
                     val_j_2 = psijo(1)
                     val_d_2 = psido(1)
                     
                     div_j = div_j + (val_j_2 - val_j_1)/delta
                     div_d = div_d + (val_d_2 - val_d_1)/delta/val_d_1
                  end do
               end do
               
               ! Restore coordinate
               xold(:,i) = x_orig(:,i)
               
               div_j = div_j / 6.0d0
               div_d = div_d / 6.0d0
               
               write(ounit,'(A,I3,I3,A,F12.6,A,F12.6,A,F12.6)', advance='no') &
                    'e-n: ', i, ic, ' r: ', epsilon, ' dJ: ', div_j, ' dD: ', div_d
               
               if (abs(div_j) > 0.001) write(ounit,'(A)', advance='yes') ' WARN: J!=0'
               if (abs(div_d) > 0.001) write(ounit,'(A)', advance='yes') ' WARN: D!=0'
               write(ounit,*) ''

            end do
         end do
         
         ! Electron-Electron (Spherical Average)
         epsilon = 0.00001d0 
         do i = 1, nelec
            do j = i+1, nelec
               div_j = 0.0d0
               div_d = 0.0d0
               
               do k = 1, 3
                  do l = -1, 1, 2
                     rvec = 0.0d0
                     rvec(k) = dble(l)
                     
                     ! Keep j fixed, move i relative to j
                     
                     ! 1. Point at epsilon
                     xold(:,i) = xold(:,j) + epsilon * rvec
                     call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)
                     val_j_1 = psijo(1)
                     
                     ! 2. Point at epsilon + delta
                     xold(:,i) = xold(:,j) + (epsilon + delta) * rvec
                     call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)
                     val_j_2 = psijo(1)
                     
                     div_j = div_j + (val_j_2 - val_j_1)/delta
                  end do
               end do
               
               ! Restore coordinate
               xold(:,i) = x_orig(:,i)
               
               div_j = div_j / 6.0d0
               
               if ( (i .le. nup .and. j .le. nup) .or. (i .gt. nup .and. j .gt. nup) ) then
                   ! Parallel
                   write(ounit,'(A,I3,I3,A,F12.6,A,F12.6)', advance='no') &
                        'e-e(P): ', i, j, ' r: ', epsilon, ' dJ: ', div_j
                   if (abs(div_j - 0.25d0) > 0.001) write(ounit,'(A)', advance='yes') ' WARN: J!=0.25'
               else
                   ! Anti-Parallel
                   write(ounit,'(A,I3,I3,A,F12.6,A,F12.6)', advance='no') &
                        'e-e(A): ', i, j, ' r: ', epsilon, ' dJ: ', div_j
                   if (abs(div_j - 0.5d0) > 0.001) write(ounit,'(A)', advance='yes') ' WARN: J!=0.5'
               end if
               write(ounit,*) ''
            end do
         end do
         
         ! Final restore of all coordinates to be safe
         xold = x_orig
         call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)
         
         stop
      end if



      ! call backflow(xold)
      ! do iel=1,nelec
      !    do k=1,3
      !      xnew2 = xold
      !      xnew2(k,iel) = xnew2(k,iel) + 1.0d0
      !      call psie(iel,xnew2,psidn,psijo,1,0)
      !      xold(k,iel) = xold(k,iel) + 1.0d0
      !      call hpsi(xold,psidnn,psijo,ekino,eold(1,1),0,1)
      !      if (abs(psidnn(1) - psidn(1)) > 1d-10) then
      !        print *, ' psidnn=', psidnn(1) , ' psidn=', psidn(1)
      !        stop
      !      end if
      !      xold(k,iel) = xold(k,iel) - 1.0d0
      !      call hpsi(xold,psidnn,psijo,ekino,eold(1,1),0,1)
      !    enddo
      ! enddo
      ! stop

      ! call backflow(xold)
      ! do iel=1,nelec
      !   do k=1,3
      !     xnew(:) = xold(:,iel)
      !     xnew(k) = xnew(k) + 1.0d0
      !     call single_rios_backflow(iel, xold, xnew, quasi_x_new, dquasi_dx_new, d2quasi_dx2_new, indices)
      !     xold(k,iel) = xold(k,iel) + 1.0d0
      !     call backflow(xold)
      !     xold(k,iel) = xold(k,iel) - 1.0d0
      !     do i = 1, nelec
      !       do kk=1,3
      !         if (abs(quasi_x_new(kk,i) - quasi_x(kk,i)) > 1d-10) then
      !           print *, 'Error in quasi_x for iel=', iel, ' k=', k, ' i=', i, ' kk=', kk
      !           print *, ' quasi_x_new=', quasi_x_new(kk,i) , ' quasi_x=', quasi_x(kk,i)
      !           stop
      !         end if
      !         do j = 1, nelec
      !           if (abs(d2quasi_dx2_new(kk,i,j) - d2quasi_dx2(kk,i,j)) > 1d-10) then
      !             print *, 'Error in d2quasi_dx2 for iel=', iel, ' k=', k, ' i=', i, ' j=', j, ' kk=', kk
      !             print *, ' d2quasi_dx2_new=', d2quasi_dx2_new(kk,i,j) , ' d2quasi_dx2=', d2quasi_dx2(kk,i,j)
      !             stop
      !           end if
      !           do l = 1, 3
      !             if (abs(dquasi_dx_new(kk,i,l,j) - dquasi_dx(kk,i,l,j)) > 1d-10) then
      !               print *, 'Error in dquasi_dx for iel=', iel, ' k=', k, ' i=', i, ' j=', j, ' l=', l, ' kk=', kk
      !               print *, ' dquasi_dx_new=', dquasi_dx_new(kk,i,l,j) , ' dquasi_dx=', dquasi_dx(kk,i,l,j)
      !               stop
      !             end if
      !           end do
      !         end do
      !       end do
      !     end do
      !     !print *, quasi_x_new(k,iel) , quasi_x(k,iel)
      !     !print *, dquasi_dx_new(k,iel,k,iel) , dquasi_dx(k,iel,k,iel)
      !     !print *, d2quasi_dx2_new(k,iel,iel) , d2quasi_dx2(k,iel, iel)
      !     call backflow(xold)
      !   enddo
      !   print *, '-----'
      ! enddo 
      ! stop

      ! do iel=1,nelec
      !   do k=1,3
      !     xold(k,iel) = xold(k,iel) + 0.0001
      !     call hpsi(xold,psidoo,psijo,ekino,eold(1,1),0,1)
      !     tmp = quasi_x(1,iel)
      !     xold(k,iel) = xold(k,iel) - 0.0001
      !     call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)
      !     tmp2 = quasi_x(1,iel)
      !     print *, (tmp-tmp2)/0.0001, dquasi_dx(1,iel,k,iel)
      !     xold(k,iel) = xold(k,iel) - 0.0001
      !     call hpsi(xold,psidoo2,psijo,ekino,eold(1,1),0,1)
      !     xold(k,iel) = xold(k,iel) + 0.0001
      !     ddx_l(k) = (tmp + quasi_x(1,iel) -2*tmp2)/0.0001/0.0001
      !   enddo
      !   print *, ddx_l(1) + ddx_l(2) + ddx_l(3)
      !   print *, d2quasi_dx2(1,iel,iel)
      !   print *, '-----'
      ! enddo 
      ! stop



      ! do iel=1,nelec
      !   do k=1,3
      !     xold(k,iel) = xold(k,iel) + 0.00001
      !     call hpsi(xold,psidoo,psijo,ekino,eold(1,1),0,1)
      !     tmp = r_en(iel,4,0,2)
      !     xold(k,iel) = xold(k,iel) - 0.00001
      !     call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)
      !     tmp2 = r_en(iel,4,0,2)
      !     print *, (tmp-tmp2)/0.00001, r_en_gl(iel,k,4,0,2)
      !     xold(k,iel) = xold(k,iel) - 0.00001
      !     call hpsi(xold,psidoo2,psijo,ekino,eold(1,1),0,1)
      !     xold(k,iel) = xold(k,iel) + 0.00001
      !     ddx_l(k) = (tmp + r_en(iel,4,0,2) -2*tmp2)/0.00001/0.00001
      !   enddo
      !   print *, ddx_l(1) + ddx_l(2) + ddx_l(3)
      !   print *, r_en_gl(iel,4,4,0,2)
      !   print *, '-----'
      ! enddo 
      ! stop

      ! do iel=1,nelec
      !   do k=1,3
      !     xold(k,iel) = xold(k,iel) + 0.00001
      !     call hpsi(xold,psidoo,psijo,ekino,eold(1,1),0,1)
      !     xold(k,iel) = xold(k,iel) - 0.00001
      !     call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)
      !     print *, ddx(k,iel,1) , (psidoo(1) - psido(1))/0.00001/psido(1)
      !     xold(k,iel) = xold(k,iel) - 0.00001
      !     call hpsi(xold,psidoo2,psijo,ekino,eold(1,1),0,1)
      !     xold(k,iel) = xold(k,iel) + 0.00001
      !     ddx_l(k) = (psidoo2(1) + psidoo(1) -2*psido(1))/0.00001/0.00001/psido(1)
      !   enddo
      !   call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)
      !   print *, ddx_l(1) + ddx_l(2) + ddx_l(3)
      !   print *, d2dx2(iel,1)
      !   print *, '-----'
      ! enddo 
      ! stop

      ! do iel=1,nelec
      !   do k=1,3
      !     call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)
      !     xold(k,iel) = xold(k,iel) + 0.00001
      !     call hpsi(xold,psidoo,psijo,ekino,eold(1,1),0,1)
      !     call compute_determinante_grad(iel,psidoo(1),psidoo,psijo,vold(1,iel),1)
      !     xold(k,iel) = xold(k,iel) - 0.00001
      !     print *, ddx(k,iel,1) , (psidoo(1) - psido(1))/0.00001/psido(1), vold(k,iel)
      !   enddo
      !   print *, '-----'
      ! enddo 
      ! stop

      ! do iel=1,nelec
      !   do k=1,3
      !     call hpsi(xold,psido,psijo,ekino,eold(1,1),0,1)
      !     xnew2 = xold
      !     xnew2(k,iel) = xnew2(k,iel) + 0.00001d0
      !     call psie(iel,xnew2,psidn,psijo,1,0)
      !     call compute_determinante_grad(iel,psidn(1),psidn,psijo,vold(1,iel),0)
      !     print *, ddx(k,iel,1) , (psidn(1) - psido(1))/0.00001/psido(1), vold(k,iel)
      !   enddo
      !   print *, '-----'
      ! enddo 
      ! stop
      
      do istate=1,nstates
        psi2o(istate,1)=2*(dlog(dabs(psido(istate)))+psijo(stoj(istate)))
        if(ipr.gt.1) write(ounit,'(''zerest STATE,psido,psijo,psi2o='',i4,3d12.4)') &
                      istate,psido(istate),psijo(stoj(istate)),psi2o(istate,1)
      enddo

      if(iguiding.gt.0) then
        call determinant_psig(psido,psijo,psidg)
! rewrite psi2o if you are sampling guiding
        psi2o(1,1)=2*(dlog(dabs(psidg)))
!        psi2o(1,1)=2*(dlog(dabs(psidg))+psijo(1))
        if(ipr.gt.1) write(ounit,'(''zerest after guiding: psig,psi2o='',2d12.4)') psidg/exp(psijo(1)),psi2o(1,1)
      endif

      if(index(mode,'all').ne.0.or.node_cutoff.gt.0) then
        do jel=1,nelec
          call compute_determinante_grad(jel,psido(1),psido,psijo,vold(1,jel),1)
        enddo
      endif

      if(node_cutoff.gt.0) then
        call nodes_distance(vold,distance_node,1)
        rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node

        psi2o(1,1)=psi2o(1,1)+2*dlog(rnorm_nodes)

        if(ipr.gt.1) then
          write(ounit,'(''distance_node='',d12.4)') distance_node
          write(ounit,'(''rnorm_nodes='',d12.4)') rnorm_nodes
          write(ounit,'(''psid2o_ncut='',f9.4)') psi2o(1,1)
          do i=1,nelec
              write(ounit,'(''vd'',3e20.10)') (vold(k,i),k=1,3)
          enddo
        endif
      endif

      if(ioptorb.gt.0) ns_current=0

      call prop_save
      call pcm_save
      call mmpol_save
      call optjas_save
      call optci_save
      call optorb_save
      call force_analy_save

! get interparticle distances
      call distances(0,xold)

      do i=1,nelec
        rmino(i)=99.d9
        do ic=1,ncent
          if(r_en(i,ic).lt.rmino(i)) then
            rmino(i)=r_en(i,ic)
            nearesto(i)=ic
          endif
        enddo
        do  k=1,3
          rvmino(k,i)=rvec_en(k,i,nearesto(i))
        enddo
      enddo

      return
      end
end module
