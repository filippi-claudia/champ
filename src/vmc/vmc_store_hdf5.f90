module vmc_store_hdf5_mod
        contains
        subroutine vmc_store_hdf5(restart_filename)
        !> @brief Store the VMC data to a HDF5 file for later restart
        !> @details This subroutine stores the VMC data to a HDF5 file for restarting purposes
        !> @author Ravindra Shinde
        !> @date 2022-11-03
        !> @email r.l.shinde@utwente.nl

        use hdf5, only: hid_t
        use custom_broadcast, only: bcast
        use hdf5_utils, only: hdf5_file_create, hdf5_file_close, hdf5_file_open
        use hdf5_utils, only: hdf5_group_create, hdf5_group_close, hdf5_group_open
        use hdf5_utils, only: hdf5_write, hdf5_read
        use mpi

        ! union of all the required arrays
        use basis,   only: ns, np, nd, nf, ng, zex
        use coefs,   only: nbasis
        use config,  only: eold,nearesto,psi2o,psido,psijo,rmino,rvmino
        use config,  only: tjfo,vold,xnew,xold
        use constants, only: hb
        use contrl_file, only: errunit,ounit
        use control, only: mode
        use control_vmc, only: vmc_idump, vmc_irstar, vmc_isite, vmc_nconf, vmc_nblk, vmc_nblk_max
        use control_vmc, only: vmc_nblkeq, vmc_nconf_new, vmc_nstep, vmc_icharged_atom, vmc_nblk_ci
        use csfs,    only: ncsf, nstates, ccsf
        use denupdn, only: rprobdn,rprobup

        use determinante_mod, only: compute_determinante_grad
        use error,   only: fatal_error
        use est2cm,  only: ecm2,ecm21,pecm2,r2cm2,tjfcm2,tpbcm2
        use estcum,  only: ecum,ecum1,iblk,pecum,r2cum,tjfcum,tpbcum
        use estsig,  only: ecm21s,ecum1s
        use estsum,  only: acc,esum,pesum,r2sum,tjfsum,tpbsum
        use force_analytic, only: force_analy_dump,force_analy_rstrt
        use forcewt, only: wcum,wsum
        use hpsi_mod, only: hpsi
        use inputflags, only: eps_node_cutoff,node_cutoff
        use metropolis, only: delta,deltar,deltat
        use mstates_ctrl, only: iguiding
        use mstates_mod, only: MSTATES
        use multiple_geo, only: fcm2,fcum,fgcm2,fgcum,iwftype,nforce,nwftype,pecent

        use mpiconf, only: idtask,nproc,wid
        use nodes_distance_mod, only: nodes_distance,rnorm_nodes_num
        use optci_mod, only: optci_dump,optci_rstrt,optci_save
        use optjas_mod, only: optjas_dump,optjas_rstrt,optjas_save
        use optorb_cblock, only: ns_current
        use optorb_f_mod, only: optorb_dump,optorb_rstrt,optorb_save
        use optwf_control, only: ioptorb, ioptci, ioptjas, ioptwf
        use optx_jas_ci, only: optx_jas_ci_dump,optx_jas_ci_rstrt
        use optx_jas_orb, only: optx_jas_orb_dump,optx_jas_orb_rstrt
        use optx_orb_ci, only: optx_orb_ci_dump,optx_orb_ci_rstrt
        use pcm_mod, only: pcm_dump,pcm_rstrt
        use precision_kinds, only: dp
        ! use prop_vmc, only: prop_save
        use properties_mod, only: prop_dump,prop_rstrt
        use pseudo,  only: nloc
        use qua,     only: nquad,wq,xq,yq,zq
        use random_mod, only: savern
        use slater,  only: cdet,coef,ndet,norb
        use stats,   only: rejmax
        use step,    only: ekin,ekin2,rprob,suc,trunfb,try
        use system,  only: cent,iwctype,ncent,ncent_tot,nctype,nctype_tot,ndn,nelec,newghostype
        use system,  only: nghostcent,nup,znuc
        use strech_mod, only: setup_force,strech
        use vmc_mod, only: norb_tot,nrad

        ! optorb
        use optorb_cblock, only: nefp_blocks,norb_f_bcum,norbprim,norbterm
        use optorb_cblock, only: nreduced
        use optwf_control, only: iapprox,ioptorb
        use orb_mat_003, only: orb_o_cum
        use orb_mat_004, only: orb_oe_cum
        use orb_mat_005, only: orb_ho_cum
        use orb_mat_006, only: orb_oo_cum
        use orb_mat_007, only: orb_oho_cum
        use orb_mat_024, only: orb_f_bcm2,orb_f_bcum
        use orb_mat_030, only: orb_ecum,orb_wcum

        ! optci
        use ci000,   only: nciprim,nciterm
        use ci005_blk, only: ci_o_cum
        use ci008_blk, only: ci_oe_cm2,ci_oe_cum
        use ci009_blk, only: ci_oo_cm2,ci_oo_cum
        use ci010_blk, only: ci_ooe_cum
        use optwf_control, only: ioptci,method

        ! properties
        use prp000,  only: iprop,nprop
        use prp003,  only: vprop_cm2,vprop_cum

        ! efficiency
        use mstates2, only: effcm2,effcum
        use mstates_ctrl, only: iefficiency,nstates_psig

        ! force analytic
        use da_energy_sumcum, only: da_energy_cm2,da_energy_cum,da_psi_cum
        use m_force_analytic, only: iforce_analy

        ! Jastrow optimization
        use gradhessj, only: d2j,d2j_e,de,de_de,de_e,dj,dj_de,dj_dj
        use gradhessj, only: dj_dj_e,dj_e,dj_e2,e2
        use gradjerr, only: grad_jas_bcm2,grad_jas_bcum
        use gradjerrb, only: ngrad_jas_bcum,ngrad_jas_blocks
        use optwf_control, only: ioptjas
        use optwf_parms, only: nparmj

        ! optx_jas_orb
        use mix_jas_orb, only: de_o,dj_ho,dj_o,dj_oe

        ! optx_orb_ci
        use mix_jas_ci, only: de_o_ci,dj_de_ci,dj_o_ci,dj_oe_ci

        ! optx_orb_ci
        use mix_orb_ci, only: ci_de_o,ci_o_ho,ci_o_o,ci_o_oe

        implicit none

        ! HDF5 related variables
        character(len=*), intent(in)  ::  restart_filename
        integer(hid_t)                 ::  file_id
        integer(hid_t)                 ::  group_id
        character(len=20)              ::  author = "CHAMP"

        integer :: i, ib, ic, id, idfrom, idget, ierr
        integer :: ifr, istate, j, k
        integer :: nelecx, nforcex, nlocx, nproco
        integer :: nq_id, nqd_id, nqx, nscounts
        integer, dimension(4,0:nproc) :: irn
        integer, dimension(MPI_STATUS_SIZE) :: istatus
        integer, dimension(4,0:nproc) :: irn_tmp
        integer, dimension(0:nproc) :: ircounts
        integer, dimension(0:nproc) :: idispls
        integer :: irequest, iw
        real(dp) :: rnd, wq_id, x_id, xq_id, yq_id, zq_id

        ! optorb
        integer :: matdim


        do i=0,nproc-1
            ircounts(i)=4
            idispls(i)=i*4
        enddo

        idispls(nproc)=4*nproc
        nscounts=ircounts(idtask)

        call savern(irn(1,idtask))

        call mpi_gatherv(irn(1,idtask),nscounts,mpi_integer,irn_tmp,ircounts,idispls,mpi_integer,0,MPI_COMM_WORLD,ierr)

        ! Bradcast the data which is spread across the processors
        call bcast(irn)
        call bcast(xold)
        call bcast(xq)
        call bcast(yq)
        call bcast(zq)

        ! Only the master process will write the data to the HDF5 file
        if (wid) then
        ! Open the HDF5 file
        write(ounit, *) " HDF5 Group saved :: HDF5 Restart file name:: ", restart_filename
        call hdf5_file_create(restart_filename, file_id)

        call hdf5_group_create(file_id, "Metadata", group_id)
        call hdf5_group_open(file_id, "Metadata", group_id)

        call get_environment_variable ("USER", author)
        call hdf5_write(file_id, group_id, "Author", author)
        call hdf5_write(file_id, group_id, " Code Compilation Date ", __DATE__)
        call hdf5_write(file_id, group_id, " Code Compilation Time ", __TIME__)

#if defined(GIT_HEAD_BRANCH)
        call hdf5_write(file_id, group_id, " Git Branch ", GIT_HEAD_BRANCH)
#endif

#if defined(GIT_REVISION_HASH)
        call hdf5_write(file_id, group_id, " Git Commit Hash ", GIT_REVISION_HASH)
#endif

#if defined(CMAKE_Fortran_COMPILER)
        call hdf5_write(file_id, group_id, " Compiler ", CMAKE_Fortran_COMPILER)
#endif

#if defined(CMAKE_Fortran_COMPILER_VERSION)
        call hdf5_write(file_id, group_id, " Compiler Version ", CMAKE_Fortran_COMPILER_VERSION)
#endif

#if defined(TARGET_ARCHITECTURE)
        call hdf5_write(file_id, group_id, " Vectorization Instructions ", TARGET_ARCHITECTURE)
#endif

#if defined(HDF5_VERSION)
        call hdf5_write(file_id, group_id, " HDF5 Version ", HDF5_VERSION)
#endif
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: Metadata "

        call hdf5_group_create(file_id, "Electrons", group_id)
        call hdf5_group_open(file_id, "Electrons", group_id)
        call hdf5_write(file_id, group_id, "Number of Up-Spin Electrons", nup)
        call hdf5_write(file_id, group_id, "Number of Down-Spin Electrons", ndn)
        call hdf5_write(file_id, group_id, "Total Number of Electrons", nelec)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: Electrons "

        call hdf5_group_create(file_id, "System", group_id)
        call hdf5_group_open(file_id, "System", group_id)
        call hdf5_write(file_id, group_id, "Number of Center Types", nctype)
        call hdf5_write(file_id, group_id, "Number of Centers", ncent)
        call hdf5_write(file_id, group_id, "Center Coordinates", cent)
        call hdf5_write(file_id, group_id, "Nuclear Charge Znuc", znuc)
        call hdf5_write(file_id, group_id, "Number of Ghost Center Types", newghostype)
        call hdf5_write(file_id, group_id, "Number of Ghost Centers", nghostcent)
        call hdf5_write(file_id, group_id, "Index of Which Center Type", iwctype)
        call hdf5_write(file_id, group_id, "PE Centers", pecent)
        call hdf5_write(file_id, group_id, "Nforce", nforce)
        call hdf5_write(file_id, group_id, "Nloc", nloc)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: System "

        call hdf5_group_create(file_id, "ECP", group_id)
        call hdf5_group_open(file_id, "ECP", group_id)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: ECP "

        call hdf5_group_create(file_id, "Basis", group_id)
        call hdf5_group_open(file_id, "Basis", group_id)
        call hdf5_write(file_id, group_id, "Zex", zex(1:nbasis,1))
        if (nloc .gt. 0) then
            call hdf5_write(file_id, group_id, "nquad", nquad)
            call hdf5_write(file_id, group_id, "xq", xq(1:nquad))
            call hdf5_write(file_id, group_id, "yq", yq(1:nquad))
            call hdf5_write(file_id, group_id, "zq", zq(1:nquad))
            call hdf5_write(file_id, group_id, "wq", wq(1:nquad))
        endif
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: Basis "

        call hdf5_group_create(file_id, "AO", group_id)
        call hdf5_group_open(file_id, "AO", group_id)
        call hdf5_write(file_id, group_id, "Number of Basis", nbasis)
        call hdf5_write(file_id, group_id, "Number of S Type AOs", ns)
        call hdf5_write(file_id, group_id, "Number of P Type AOs", np)
        call hdf5_write(file_id, group_id, "Number of D Type AOs", nd)
        call hdf5_write(file_id, group_id, "Number of F Type AOs", nf)
        call hdf5_write(file_id, group_id, "Number of G Type AOs", ng)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: AO "

        call hdf5_group_create(file_id, "MO", group_id)
        call hdf5_group_open(file_id, "MO", group_id)
        call hdf5_write(file_id, group_id, "Number of Orbitals", norb)
        call hdf5_write(file_id, group_id, "Number of Orbitals Total", norb_tot)
        call hdf5_write(file_id, group_id, "MO Coefficients", coef(1:nbasis,1:norb,1))
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: MO "

        call hdf5_group_create(file_id, "Determinants", group_id)
        call hdf5_group_open(file_id, "Determinants", group_id)
        call hdf5_write(file_id, group_id, "Number of Determinants", ndet)
        call hdf5_write(file_id, group_id, "Determinant Coefficients", cdet(1:ndet,1,1))
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: Determinants "

        call hdf5_group_create(file_id, "CSFs", group_id)
        call hdf5_group_open(file_id, "CSFs", group_id)
        call hdf5_write(file_id, group_id, "Number of CSFs", ncsf)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: CSFs "

        call hdf5_group_create(file_id, "States", group_id)
        call hdf5_group_open(file_id, "States", group_id)
        call hdf5_write(file_id, group_id, "Number of States", nstates)
        call hdf5_write(file_id, group_id, "iguiding", iguiding)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: States "

        call hdf5_group_create(file_id, "UnitCell", group_id)
        call hdf5_group_open(file_id, "UnitCell", group_id)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: UnitCell "

        call hdf5_group_create(file_id, "Periodic", group_id)
        call hdf5_group_open(file_id, "Periodic", group_id)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: Periodic "

        call hdf5_group_create(file_id, "QMC", group_id)
        call hdf5_group_open(file_id, "QMC", group_id)
        call hdf5_write(file_id, group_id, "Mode", mode)
        call hdf5_write(file_id, group_id, "Number of Processors", nproc)
        call hdf5_write(file_id, group_id, "Random Numbers Each Processor", irn_tmp(1:4,0:nproc-1))
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: QMC "

        call hdf5_group_create(file_id, "VMC", group_id)
        call hdf5_group_open(file_id, "VMC", group_id)
        call hdf5_write(file_id, group_id, "Number of VMC Blocks", vmc_nblk)
        call hdf5_write(file_id, group_id, "Number of Equilibration VMC Blocks", vmc_nblkeq)
        call hdf5_write(file_id, group_id, "Maximum Number of VMC Blocks", vmc_nblk_max)
        call hdf5_write(file_id, group_id, "Number of VMC Steps per Block", vmc_nstep)
        call hdf5_write(file_id, group_id, "Number of VMC Configurations", vmc_nconf)
        call hdf5_write(file_id, group_id, "Number of New VMC Configurations", vmc_nconf_new)
        call hdf5_write(file_id, group_id, "iblk", iblk)
        call hdf5_write(file_id, group_id, "xold", xold)
        call hdf5_write(file_id, group_id, "Delta", delta)
        call hdf5_write(file_id, group_id, "Deltar", deltar)
        call hdf5_write(file_id, group_id, "Deltat", deltat)


        call hdf5_write(file_id, group_id, "ecum1", ecum1(1:nstates))
        call hdf5_write(file_id, group_id, "ecum", ecum(1:nstates,1:nforce))
        call hdf5_write(file_id, group_id, "pecum", pecum(1:nstates))
        call hdf5_write(file_id, group_id, "tpbcum", tpbcum(1:nstates))
        call hdf5_write(file_id, group_id, "r2cum", r2cum)
        call hdf5_write(file_id, group_id, "acc", acc)
        call hdf5_write(file_id, group_id, "ecm21", ecm21(1:nstates))
        call hdf5_write(file_id, group_id, "ecm2", ecm2(1:nstates,1:nforce))
        call hdf5_write(file_id, group_id, "pecm2", pecm2(1:nstates))
        call hdf5_write(file_id, group_id, "tpbcm2", tpbcm2(1:nstates))
        call hdf5_write(file_id, group_id, "r2cm2", r2cm2)

        if (nforce .gt. 1) then
            call hdf5_write(file_id, group_id, "wcum", wcum(1:nstates,1:nforce))
            call hdf5_write(file_id, group_id, "fcum", fcum(1:nstates,1:nforce))
            call hdf5_write(file_id, group_id, "fcm2", fcm2(1:nstates,1:nforce))
        else
            call hdf5_write(file_id, group_id, "wcum", wcum(1:nstates,1))
        end if

        call hdf5_write(file_id, group_id, "ecum1s", ecum1s(1:nstates))
        call hdf5_write(file_id, group_id, "ecm21s", ecm21s(1:nstates))

        call hdf5_write(file_id, group_id, "try", try(1:nrad))
        call hdf5_write(file_id, group_id, "suc", suc(1:nrad))
        call hdf5_write(file_id, group_id, "trunfb", trunfb(1:nrad))
        call hdf5_write(file_id, group_id, "rprob", rprob(1:nrad))
        call hdf5_write(file_id, group_id, "rprobup", rprobup(1:nrad))
        call hdf5_write(file_id, group_id, "rprobdn", rprobdn(1:nrad))
        call hdf5_write(file_id, group_id, "ekin", ekin(1:nrad))
        call hdf5_write(file_id, group_id, "ekin2", ekin2(1:nrad))
        call hdf5_write(file_id, group_id, "rejmax", rejmax)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: VMC "

        ! Orbital Optimization
        if(ioptwf.ne.0 .and. ioptorb.ne.0) then
                matdim=nreduced*(nreduced+1)/2
                if(iapprox.gt.0) matdim=nreduced
                call hdf5_group_create(file_id, "Orbital Optimization", group_id)
                call hdf5_group_open(file_id, "Orbital Optimization", group_id)
                call hdf5_write(file_id, group_id, "norbprim", norbprim)
                call hdf5_write(file_id, group_id, "norbterm", norbterm)
                call hdf5_write(file_id, group_id, "nreduced", nreduced)
                call hdf5_write(file_id, group_id, "nefp_blocks", nefp_blocks)
                call hdf5_write(file_id, group_id, "norb_f_bcum", norb_f_bcum)
                call hdf5_write(file_id, group_id, "orb_o_cum", orb_o_cum(1:norbterm,1:nstates))
                call hdf5_write(file_id, group_id, "orb_oe_cum", orb_oe_cum(1:norbterm,1:nstates))
                call hdf5_write(file_id, group_id, "orb_ho_cum", orb_ho_cum(1:norbterm,1:nstates))
                call hdf5_write(file_id, group_id, "orb_f_bcum", orb_f_bcum(1:norbterm,1:nstates))
                call hdf5_write(file_id, group_id, "orb_f_bcm2", orb_f_bcm2(1:norbterm,1:nstates))
                call hdf5_write(file_id, group_id, "orb_oo_cum", orb_oo_cum(1:matdim,1:nstates))
                call hdf5_write(file_id, group_id, "orb_oho_cum", orb_oho_cum(1:nreduced*nreduced,1:nstates))
                call hdf5_write(file_id, group_id, "orb_wcum", orb_wcum(1:nstates))
                call hdf5_write(file_id, group_id, "orb_ecum", orb_ecum(1:nstates))
                call hdf5_group_close(group_id)
                write(ounit, *) " HDF5 Group saved :: Orbital Optimization "
        endif


        ! CI optimization
        if (ioptwf.ne.0 .and.  .not.(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d')) then
                call hdf5_group_create(file_id, "CI Optimization", group_id)
                call hdf5_group_open(file_id, "CI Optimization", group_id)
                matdim=nciterm*(nciterm+1)/2

                call hdf5_write(file_id, group_id, "nciprim", nciprim)
                call hdf5_write(file_id, group_id, "nciterm", nciterm)
                call hdf5_write(file_id, group_id, "ci_o_cum", ci_o_cum(1:nciterm,1:nstates))
                call hdf5_write(file_id, group_id, "ci_oe_cum", ci_oe_cum(1:nciterm,1:nciterm,1:nstates))
                call hdf5_write(file_id, group_id, "ci_oe_cm2", ci_oe_cm2(1:nciterm,1:nciterm,1:nstates))
                call hdf5_write(file_id, group_id, "ci_oo_cum", ci_oo_cum(1:matdim,1:nstates))
                call hdf5_write(file_id, group_id, "ci_oo_cm2", ci_oo_cm2(1:matdim,1:nstates))
                call hdf5_write(file_id, group_id, "ci_ooe_cum", ci_ooe_cum(1:matdim,1:nstates))
                call hdf5_group_close(group_id)
                write(ounit, *) " HDF5 Group saved :: CI Optimization "
        endif


        ! properties
        if (iprop.ne.0) then
                call hdf5_group_create(file_id, "Properties", group_id)
                call hdf5_group_open(file_id, "Properties", group_id)
                call hdf5_write(file_id, group_id, "iprop", iprop)
                call hdf5_write(file_id, group_id, "nprop", nprop)
                call hdf5_write(file_id, group_id, "vprop_cum", vprop_cum(1:nprop))
                call hdf5_write(file_id, group_id, "vprop_cm2", vprop_cm2(1:nprop))
                call hdf5_group_close(group_id)
                write(ounit, *) " HDF5 Group saved :: Properties "
        endif


        ! Efficiency
        if (iefficiency.ne.0) then
                call hdf5_group_create(file_id, "Efficiency", group_id)
                call hdf5_group_open(file_id, "Efficiency", group_id)
                call hdf5_write(file_id, group_id, "nstates_psig", nstates_psig)
                call hdf5_write(file_id, group_id, "effcum", effcum(1:nstates_psig))
                call hdf5_write(file_id, group_id, "effcm2", effcm2(1:nstates_psig))
                call hdf5_group_close(group_id)
                write(ounit, *) " HDF5 Group saved :: Efficiency "
        endif


        ! Force Analytical
        if (iforce_analy.ne.0) then
                call hdf5_group_create(file_id, "Force Analytical", group_id)
                call hdf5_group_open(file_id, "Force Analytical", group_id)
                call hdf5_write(file_id, group_id, "da_energy_cum", da_energy_cum(1:3,1:ncent))
                call hdf5_write(file_id, group_id, "da_psi_cum", da_psi_cum(1:3,1:ncent))
                call hdf5_write(file_id, group_id, "da_energy_cm2", da_energy_cm2(1:3,1:ncent))
                call hdf5_group_close(group_id)
                write(ounit, *) " HDF5 Group saved :: Force Analytical "
        endif


        ! Jastrow Optimization
        if (ioptwf.ne.0 .and. ioptjas.ne.0) then
                call hdf5_group_create(file_id, "Jastrow Optimization", group_id)
                call hdf5_group_open(file_id, "Jastrow Optimization", group_id)
                call hdf5_write(file_id, group_id, "nparmj", nparmj)
                call hdf5_write(file_id, group_id, "dj", dj(1:nparmj,1:nstates))
                call hdf5_write(file_id, group_id, "de", de(1:nparmj,1:nstates))
                call hdf5_write(file_id, group_id, "dj_e", dj_e(1:nparmj,1:nstates))
                call hdf5_write(file_id, group_id, "de_e", de_e(1:nparmj,1:nstates))
                call hdf5_write(file_id, group_id, "dj_e2", dj_e2(1:nparmj,1:nstates))
                call hdf5_write(file_id, group_id, "e2", e2(1:nparmj,1:nstates))
                call hdf5_write(file_id, group_id, "dj_de", dj_de(1:nparmj,1:nparmj,1:nstates))
                call hdf5_write(file_id, group_id, "dj_dj", dj_dj(1:nparmj,1:nparmj,1:nstates))
                call hdf5_write(file_id, group_id, "dj_dj_e", dj_dj_e(1:nparmj,1:nparmj,1:nstates))
                call hdf5_write(file_id, group_id, "d2j", d2j(1:nparmj,1:nparmj,1:nstates))
                call hdf5_write(file_id, group_id, "d2j_e", d2j_e(1:nparmj,1:nparmj,1:nstates))
                call hdf5_write(file_id, group_id, "de_de", de_de(1:nparmj,1:nparmj,1:nstates))
                if(ngrad_jas_blocks.gt.0) then
                        call hdf5_write(file_id, group_id, "grad_jas_bcum", grad_jas_bcum(1:nparmj,1:nstates))
                        call hdf5_write(file_id, group_id, "grad_jas_bcm2", grad_jas_bcm2(1:nparmj,1:nstates))
                        call hdf5_write(file_id, group_id, "ngrad_jas_bcum", ngrad_jas_bcum)
                endif
                call hdf5_group_close(group_id)
                write(ounit, *) " HDF5 Group saved :: Jastrow Optimization "
        endif



        ! Opt Jas Orb
        if (.not. (ioptjas.eq.0.or.ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d')) then
                call hdf5_group_create(file_id, "Optx_Jas_Orb", group_id)
                call hdf5_group_open(file_id, "Optx_Jas_Orb", group_id)
                call hdf5_write(file_id, group_id, "dj_o", dj_o(1:nparmj,1:nreduced,1:nstates))
                call hdf5_write(file_id, group_id, "dj_oe", dj_oe(1:nparmj,1:nreduced,1:nstates))
                call hdf5_write(file_id, group_id, "dj_ho", dj_ho(1:nparmj,1:nreduced,1:nstates))
                call hdf5_write(file_id, group_id, "de_o", de_o(1:nparmj,1:nreduced,1:nstates))
                call hdf5_group_close(group_id)
                write(ounit, *) " HDF5 Group saved :: Opt Jas Orb "
        endif



        ! Opt Jas CI
        if (.not. (ioptjas.eq.0.or.ioptci.eq.0)) then
                call hdf5_group_create(file_id, "Optx_Jas_CI", group_id)
                call hdf5_group_open(file_id, "Optx_Jas_CI", group_id)
                call hdf5_write(file_id, group_id, "dj_o_ci", dj_o_ci(1:nparmj,1:nciterm,1:nstates))
                call hdf5_write(file_id, group_id, "dj_oe_ci", dj_oe_ci(1:nparmj,1:nciterm,1:nstates))
                call hdf5_write(file_id, group_id, "dj_de_ci", dj_de_ci(1:nparmj,1:nciterm,1:nstates))
                call hdf5_write(file_id, group_id, "de_o_ci", de_o_ci(1:nparmj,1:nciterm,1:nstates))
                call hdf5_group_close(group_id)
                write(ounit, *) " HDF5 Group saved :: Opt Jas CI "
        endif




        ! Opt Orb CI
        if (.not. (ioptorb.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d')) then
                call hdf5_group_create(file_id, "Optx_Orb_CI", group_id)
                call hdf5_group_open(file_id, "Optx_Orb_CI", group_id)
                call hdf5_write(file_id, group_id, "ci_o_o", ci_o_o(1:nciterm,1:nreduced))
                call hdf5_write(file_id, group_id, "ci_o_oe", ci_o_oe(1:nciterm,1:nreduced))
                call hdf5_write(file_id, group_id, "ci_o_ho", ci_o_ho(1:nciterm,1:nreduced))
                call hdf5_write(file_id, group_id, "ci_de_o", ci_de_o(1:nciterm,1:nreduced))
                call hdf5_group_close(group_id)
                write(ounit, *) " HDF5 Group saved :: Opt Orb CI "
        endif



        call hdf5_file_close(file_id)
        ! Close the HDF5 file

        endif   ! master thread

        end subroutine vmc_store_hdf5
end module vmc_store_hdf5_mod
