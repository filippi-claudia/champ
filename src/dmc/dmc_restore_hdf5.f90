module dmc_restore_hdf5_mod
        contains
        subroutine dmc_restore_hdf5(restart_filename)
        !> @brief Restore the VMC data to a HDF5 file for later restart
        !> @details This subroutine restores the VMC data to a HDF5 file for restarting purposes
        !> @author Ravindra Shinde
        !> @date 2023-06-15
        !> @email r.l.shinde@utwente.nl

        ! From restart.f
        use age,     only: iage,ioldest,ioldestmx
        use basis,   only: zex, ns, np, nd, nf, ng
        use branch,  only: eest,eigv,eold,ff,fprod,nwalk,wdsumo,wgdsumo,wt
        use branch,  only: wtgen
        use casula,  only: i_vpsp,icasula
        use coefs,   only: nbasis
        use config,  only: psido_dmc,psijo_dmc,vold_dmc,xold_dmc
        use constants, only: hb
        use contrl_file, only: ounit
        use contrldmc, only: idmc,nfprod,rttau,tau
        use control, only: ipr,mode
        use control_dmc, only: dmc_nconf
        use denupdn, only: rprobdn,rprobup
        use derivest, only: derivcm2,derivcum,derivsum,derivtotave_num_old
        use da_energy_sumcum, only: da_energy_cm2,da_energy_cum,da_psi_cum
        use force_pth, only: PTH
        use m_force_analytic, only: iforce_analy
        use multiple_geo, only: nwprod
        use pathak_mod, only: ipathak,pold
        use vd_mod, only: da_branch_cum,deriv_eold,dmc_ivd,ehist,esnake
        use determinante_mod, only: compute_determinante_grad
        use error,   only: fatal_error
        use est2cm,  only: ecm21_dmc,ecm2_dmc,efcm2,efcm21,egcm2,egcm21
        use est2cm,  only: ei1cm2,ei2cm2,ei3cm2,pecm2_dmc,r2cm2_dmc,ricm2
        use est2cm,  only: tpbcm2_dmc,wcm2,wcm21,wdcm2,wdcm21
        use est2cm,  only: wfcm2,wfcm21,wgcm2,wgcm21,wgdcm2
        use estcum,  only: ecum1_dmc,ecum_dmc,efcum,efcum1,egcum,egcum1
        use estcum,  only: ei1cum,ei2cum,ei3cum,iblk,ipass,pecum_dmc
        use estcum,  only: r2cum_dmc,ricum,taucum,tpbcum_dmc
        use estcum,  only: wcum1,wcum_dmc,wdcum,wdcum1,wfcum,wfcum1,wgcum
        use estcum,  only: wgcum1,wgdcum
        use estsum,  only: efsum,egsum,ei1sum,ei2sum,esum_dmc,pesum_dmc
        use estsum,  only: wsum1
        use estsum,  only: r2sum,risum,tausum,tpbsum_dmc,wdsum
        use estsum,  only: wfsum,wgdsum,wgsum,wsum_dmc
        use general, only: write_walkalize
        use hpsi_mod, only: hpsi
        use jacobsave, only: ajacob,ajacold
        use mmpol,   only: mmpol_init,mmpol_rstrt
        use mmpol_dmc, only: mmpol_save
        use mpi
        use mpiblk,  only: iblk_proc
        use mpiconf, only: idtask,nproc,wid
        use multiple_geo, only: fgcm2,fgcum,istrech,nforce,pecent
        use nonloc_grid_mod, only: t_vpsp_sav
        use pcm_dmc, only: pcm_save
        use pcm_mod, only: pcm_init,pcm_rstrt
        use precision_kinds, only: dp
        use prop_dmc, only: prop_save_dmc
        use properties_mod, only: prop_init,prop_rstrt
        use pseudo,  only: nloc
        use qua,     only: nquad,wq,xq,yq,zq
        use random_mod, only: setrn
        use restart_gpop, only: startr_gpop
        use slater,  only: cdet,coef,ndet,norb
        use stats,   only: acc,dfus2ac,dfus2un,dr2ac,dr2un,nacc,nbrnch
        use stats,   only: nodecr,trymove
        use step,    only: rprob
        use strech_mod, only: strech
        use system,  only: cent,iwctype,ncent,ncent_tot,nctype,ndn,nelec
        use system,  only: nghostcent,nup,znuc
        use velratio, only: fratio
        use vmc_mod, only: norb_tot,nrad
        use walksav_det_mod, only: walksav_det
        use walksav_jas_mod, only: walksav_jas        


        ! union of all the required arrays
        use age,     only: iage,ioldest,ioldestmx
        use basis,   only: ns, np, nd, nf, ng, zex
        use branch,  only: eest,eigv,ff,fprod,nwalk,wdsumo,wgdsumo,wt
        use branch,  only: wtgen
        use coefs,   only: nbasis
        use config,  only: xold_dmc
        use constants, only: hb
        use contrl_file, only: ounit
        use contrldmc, only: idmc,nfprod,rttau,tau
        use control, only: mode
        use control_dmc, only: dmc_idump,dmc_irstar,dmc_isite,dmc_nblk
        use control_dmc, only: dmc_nblkeq,dmc_nconf,dmc_nconf_new
        use control_dmc, only: dmc_nstep
        use csfs,    only: ncsf, nstates, ccsf
        use denupdn, only: rprobdn,rprobup
        use derivest, only: derivcm2,derivcum,derivtotave_num_old
        use dmc_mod, only: MWALK
        use dumper_gpop_mod, only: dumper_gpop
        use est2cm,  only: ecm21_dmc,ecm2_dmc,efcm2,efcm21,egcm2,egcm21
        use est2cm,  only: ei1cm2,ei2cm2,ei3cm2,pecm2_dmc,r2cm2_dmc,ricm2
        use est2cm,  only: tpbcm2_dmc,wcm2,wcm21,wdcm2,wdcm21
        use est2cm,  only: wfcm2,wfcm21,wgcm2,wgcm21,wgdcm2
        use estcum,  only: ecum1_dmc,ecum_dmc,efcum,efcum1,egcum,egcum1
        use estcum,  only: ei1cum,ei2cum,ei3cum,iblk,ipass,pecum_dmc
        use estcum,  only: r2cum_dmc,ricum,taucum,tpbcum_dmc
        use estcum,  only: wcum1,wcum_dmc,wdcum,wdcum1,wfcum,wfcum1,wgcum
        use estcum,  only: wgcum1,wgdcum
        use jacobsave, only: ajacob
        use mmpol,   only: mmpol_dump
        use mpi
        use mpiblk,  only: iblk_proc
        use mpiconf, only: idtask,nproc,wid
        use multiple_geo, only: fgcm2,fgcum,nforce,pecent
        use mstates_ctrl, only: iguiding
        use pcm_mod, only: pcm_dump
        use precision_kinds, only: dp
        use properties_mod, only: prop_dump
        use pseudo,  only: nloc
        use qua,     only: nquad,wq,xq,yq,zq
        use random_mod, only: savern
        use slater,  only: cdet,coef,ndet,norb
        use stats,   only: acc,dfus2ac,dfus2un,dr2ac,dr2un,nacc,nbrnch
        use stats,   only: nodecr,trymove
        use step,    only: rprob
        use strech_mod, only: strech
        use system,  only: cent,iwctype,ncent,nctype,ndn,nelec,newghostype
        use system,  only: nghostcent,nup,znuc
        use velratio, only: fratio
        use vmc_mod, only: nrad

        ! properties
        use prp000,  only: iprop,nprop
        use prp003,  only: vprop_cm2,vprop_cum




        use hdf5, only: hid_t
        use custom_broadcast, only: bcast
        use hdf5_utils, only: hdf5_file_create, hdf5_file_close, hdf5_file_open
        use hdf5_utils, only: hdf5_group_create, hdf5_group_close, hdf5_group_open
        use hdf5_utils, only: hdf5_write, hdf5_read
        use mpi

        implicit none

        ! HDF5 related variables
        character(len=*), intent(in)    ::  restart_filename
        ! character(len=*)       ::  date, time, read_author
        ! character(len=*)       ::  git_branch, git_commit
        ! character(len=*)       ::  compiler, compiler_version
        ! character(len=*)       ::  arch, hdf5_version
        integer(hid_t)                  ::  file_id
        integer(hid_t)                  ::  group_id


        integer :: i, iage_id, ib, ic, id
        integer :: ie, ierr, ifr, ioldest_id, ioldestmx_id
        integer :: iw, j, k, n1_id
        integer :: n2_id, nbasx, ncentx, nctypex
        integer :: ndetx, ndnx, nelecx, newghostypex
        integer :: nghostcentx, nprock, nproco, nq_id, num
        integer :: nupx, nwalk_id
        character(len=20) :: s
        integer, dimension(8, 0:nproc-1) :: irn_tmp
        integer, dimension(MPI_STATUS_SIZE) :: istatus
        character(len=32) :: cnum
        character(len=12) :: mode_stored
        integer :: iph, nwprod_stored, ipathak_stored
        integer :: nblk_stored, nstep_stored, nconf_stored
        integer :: nelec_input, nforce_input, nquad_input
        real(dp) :: tau_input, hb_input
        integer, dimension(0:nproc-1) :: nwalk_all
        real(dp), allocatable :: deriv_eold_send(:,:,:), esnake_send(:,:,:,:)
        real(dp), allocatable :: ehist_send(:,:,:,:,:), pold_send(:,:)
        integer :: nwalk_send, ioldest_send, ioldestmx_send
        real(dp) :: fprod_send, eigv_send, eest_send, wdsumo_send
        real(dp), allocatable :: xold_dmc_send(:,:,:,:)
        real(dp), allocatable :: wt_send(:), ff_send(:), fratio_send(:,:)
        integer, allocatable :: iage_send(:)
        integer, dimension(nctype)      :: nsx,npx,ndx,nfx,ngx
        real(dp) :: different, eest_id
        real(dp) :: eigv_id, ff_id, fmt, fprod_id
        real(dp) :: fratio_id, hbx, taux, wdsumo_id
        real(dp) :: wq_id, wt_id, xold_dmc_id, xq_id
        real(dp) :: yq_id, zq_id
        real(dp) :: ekino(1)
        real(dp), dimension(nbasis, norb_tot) :: coefx
        real(dp), dimension(nbasis) :: zexx
        real(dp), dimension(3, ncent_tot) :: centx
        real(dp), dimension(ncent_tot) :: znucx
        real(dp), dimension(ndet) :: cdetx
        real(dp), parameter :: zero = 0.d0
        real(dp), parameter :: one = 1.d0
        real(dp), parameter :: small = 1.e-6
        real(dp) :: egave_rstrt, peave_rstrt, tpbave_rstrt
        real(dp) :: egerr_rstrt, peerr_rstrt, tpberr_rstrt, rn_eff

        character*13 filename

        ! Remember the values coming from the restart input before any of them is
        ! overwritten by a dataset, so the two sets can be reported side by side.
        nelec_input  = nelec
        nforce_input = nforce
        nquad_input  = nquad
        tau_input    = tau
        hb_input     = hb

        ! Only the master process will read the data to the HDF5 file
        if (wid) then

        ! Open the HDF5 file
        write(ounit, *) " HDF5 Restart file name:: ", restart_filename
        call hdf5_file_open(restart_filename, file_id)

        ! call hdf5_group_open(file_id, "Metadata", group_id)
        ! call hdf5_read(file_id, group_id, "Author", read_author)
        ! write(ounit, *) " Author:: ", read_author
        ! call hdf5_read(file_id, group_id, " Code Compilation Date ",date)
        ! write(ounit, *) " Code Compilation Date:: ", date
        ! call hdf5_read(file_id, group_id, " Code Compilation Time ", time)
        ! write(ounit, *) " Code Compilation Time:: ", time
        ! call hdf5_read(file_id, group_id, " Git Branch ", git_branch)
        ! write(ounit, *) " Git Branch:: ", git_branch
        ! call hdf5_read(file_id, group_id, " Git Commit Hash ", git_commit)
        ! write(ounit, *) " Git Commit Hash:: ", git_commit
        ! call hdf5_read(file_id, group_id, " Compiler ", compiler)
        ! write(ounit, *) " Compiler:: ", compiler
        ! call hdf5_read(file_id, group_id, " Compiler Version ", compiler_version)
        ! write(ounit, *) " Compiler Version:: ", compiler_version
        ! call hdf5_read(file_id, group_id, " Vectorization Instructions ", arch)
        ! write(ounit, *) " Vectorization Instructions:: ", arch
        ! call hdf5_read(file_id, group_id, " HDF5 Version ", hdf5_version)
        ! write(ounit, *) " HDF5 Version:: ", hdf5_version
        ! call hdf5_group_close(group_id)
        ! write(ounit, *) " HDF5 Group read :: Metadata "

        call hdf5_group_open(file_id, "Electrons", group_id)
        call hdf5_read(file_id, group_id, "Number of Up-Spin Electrons", nup)
        call hdf5_read(file_id, group_id, "Number of Down-Spin Electrons", ndn)
        call hdf5_read(file_id, group_id, "Total Number of Electrons", nelec)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: Electrons "

        call hdf5_group_open(file_id, "System", group_id)
        call hdf5_read(file_id, group_id, "Number of Center Types", nctype)
        call hdf5_read(file_id, group_id, "Number of Centers", ncent)
        call hdf5_read(file_id, group_id, "Center Coordinates", cent)
        call hdf5_read(file_id, group_id, "Nuclear Charge Znuc", znuc)
        call hdf5_read(file_id, group_id, "Number of Ghost Center Types", newghostype)
        call hdf5_read(file_id, group_id, "Number of Ghost Centers", nghostcent)
        call hdf5_read(file_id, group_id, "Index of Which Center Type", iwctype)
        call hdf5_read(file_id, group_id, "PE Centers", pecent)
        call hdf5_read(file_id, group_id, "Nuclear Charge Znuc", znuc)
        call hdf5_read(file_id, group_id, "Nforce", nforce)
        call hdf5_read(file_id, group_id, "Nloc", nloc)
        call hdf5_read(file_id, group_id, "hb", hbx)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: System "

        call hdf5_group_open(file_id, "ECP", group_id)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: ECP "

        call hdf5_group_open(file_id, "Basis", group_id)
        call hdf5_read(file_id, group_id, "Zex", zex)
        if (nloc .gt. 0) then
            call hdf5_read(file_id, group_id, "nquad", nquad)
            call hdf5_read(file_id, group_id, "xq", xq(1:nquad))
            call hdf5_read(file_id, group_id, "yq", yq(1:nquad))
            call hdf5_read(file_id, group_id, "zq", zq(1:nquad))
            call hdf5_read(file_id, group_id, "wq", wq(1:nquad))
        endif
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: Basis "

        call hdf5_group_open(file_id, "AO", group_id)
        call hdf5_read(file_id, group_id, "Number of Basis", nbasis)
        call hdf5_read(file_id, group_id, "Number of S Type AOs", ns)
        call hdf5_read(file_id, group_id, "Number of P Type AOs", np)
        call hdf5_read(file_id, group_id, "Number of D Type AOs", nd)
        call hdf5_read(file_id, group_id, "Number of F Type AOs", nf)
        call hdf5_read(file_id, group_id, "Number of G Type AOs", ng)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: AO "

        call hdf5_group_open(file_id, "MO", group_id)
        call hdf5_read(file_id, group_id, "Number of Orbitals", norb)
        call hdf5_read(file_id, group_id, "Number of Orbitals Total", norb_tot)
        call hdf5_read(file_id, group_id, "MO Coefficients", coef(1:nbasis,1:norb,1))
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: MO "

        call hdf5_group_open(file_id, "Determinants", group_id)
        call hdf5_read(file_id, group_id, "Number of Determinants", ndet)
        call hdf5_read(file_id, group_id, "Determinant Coefficients", cdet(1:ndet,1,1))
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: Determinants "

        call hdf5_group_open(file_id, "CSFs", group_id)
        call hdf5_read(file_id, group_id, "Number of CSFs", ncsf)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: CSFs "

        call hdf5_group_open(file_id, "States", group_id)
        call hdf5_read(file_id, group_id, "Number of States", nstates)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: States "

        call hdf5_group_open(file_id, "UnitCell", group_id)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: UnitCell "

        call hdf5_group_open(file_id, "Periodic", group_id)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: Periodic "

        call hdf5_group_open(file_id, "QMC", group_id)
        call hdf5_read(file_id, group_id, "Mode", mode_stored)
        if (mode_stored .ne. mode) &
          write(ounit, '(a,a,a,a,a)') "Warning: HDF5 restart mode '", trim(mode_stored), &
              "' differs from current mode '", trim(mode), "' -- using current mode"
        call hdf5_read(file_id, group_id, "Number of Processors", nproco)
        if(nproco.ne.nproc) then
          write(ounit, '(a)' )          "Error: different number of processors in the restart file"
          write(ounit, '(a,i4,a,i4,a)') "Number of processors from restart file = ", nproco, ", Current number of processors = ", nproc
          call fatal_error('DMC_RESTORE_HDF5: different num procs')
        endif
        call hdf5_read(file_id, group_id, "Random Numbers Each Processor", irn_tmp(1:8,0:nproc-1))
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: QMC "

        call hdf5_group_open(file_id, "DMC", group_id)
        call hdf5_read(file_id, group_id, "Number of DMC Blocks", nblk_stored)
        call hdf5_read(file_id, group_id, "Number of DMC Steps per Block", nstep_stored)
        call hdf5_read(file_id, group_id, "Number of DMC Configurations ", nconf_stored)
        call hdf5_read(file_id, group_id, "Number of Processors", nproco)
        call hdf5_read(file_id, group_id, "tau", taux)

        ! Report the control parameters found in the checkpoint next to the ones
        ! parsed from the restart input, and say which of the two is used.
        write(ounit,'(/,a)') " DMC restart :: control parameters"
        write(ounit,'(a)')   " parameter        checkpoint        input       used"
        write(ounit,'(a,i12,i12,4x,a)') " dmc_nblk    ", nblk_stored,  dmc_nblk,  "input"
        write(ounit,'(a,i12,i12,4x,a)') " dmc_nstep   ", nstep_stored, dmc_nstep, "input"
        write(ounit,'(a,i12,i12,4x,a)') " dmc_nconf   ", nconf_stored, dmc_nconf, "checkpoint"
        write(ounit,'(a,i12,i12,4x,a)') " nproc       ", nproco,       nproc,     "input"
        write(ounit,'(a,i12,i12,4x,a)') " nelec       ", nelec,  nelec_input,  "checkpoint"
        write(ounit,'(a,i12,i12,4x,a)') " nforce      ", nforce, nforce_input, "checkpoint"
        write(ounit,'(a,i12,i12,4x,a)') " nquad       ", nquad,  nquad_input,  "checkpoint"
        write(ounit,'(a,f12.6,f12.6,4x,a)')   " hb          ", hbx,  hb_input,  "input"
        write(ounit,'(a,f12.6,f12.6,4x,a,/)') " tau         ", taux, tau_input, "checkpoint"

        ! dmc_nconf is part of the branching state (it is the target population that
        ! produced the dumped weights), so it is taken from the checkpoint, exactly as
        ! the legacy binary startr does. The block/step counts govern how long *this*
        ! run lasts and are therefore taken from the restart input, as in VMC.
        dmc_nconf = nconf_stored

        allocate(xold_dmc_send(3, nelec, MWALK, nforce))
        allocate(wt_send(MWALK))
        allocate(ff_send(0:nfprod))
        allocate(fratio_send(MWALK, nforce))
        allocate(iage_send(MWALK))

        do id=0, nproco-1
            write (unit=s,fmt="(i0)") id
            if (id .eq. 0) then
                call hdf5_read(file_id, group_id, "Number of Walkers proc_"//trim(s), nwalk)
                nwalk_all(id) = nwalk
                call hdf5_read(file_id, group_id, "xold_dmc_proc_"//trim(s), xold_dmc)
                call hdf5_read(file_id, group_id, "nfprod_proc_"//trim(s), nfprod)
                call hdf5_read(file_id, group_id, "ff_proc_"//trim(s), ff)
                call hdf5_read(file_id, group_id, "wt_proc_"//trim(s), wt)
                call hdf5_read(file_id, group_id, "fprod_proc_"//trim(s), fprod)
                call hdf5_read(file_id, group_id, "eigv_proc_"//trim(s), eigv)
                call hdf5_read(file_id, group_id, "eest_proc_"//trim(s), eest)
                call hdf5_read(file_id, group_id, "wdsumo_proc_"//trim(s), wdsumo)
                call hdf5_read(file_id, group_id, "iage_proc_"//trim(s), iage)
                call hdf5_read(file_id, group_id, "ioldest_proc_"//trim(s), ioldest)
                call hdf5_read(file_id, group_id, "ioldestmx_proc_"//trim(s), ioldestmx)
                call hdf5_read(file_id, group_id, "fratio_proc_"//trim(s), fratio)
            else
                call hdf5_read(file_id, group_id, "Number of Walkers proc_"//trim(s), nwalk_send)
                nwalk_all(id) = nwalk_send
                call hdf5_read(file_id, group_id, "xold_dmc_proc_"//trim(s), xold_dmc_send)
                call hdf5_read(file_id, group_id, "nfprod_proc_"//trim(s), nfprod)
                call hdf5_read(file_id, group_id, "ff_proc_"//trim(s), ff_send)
                call hdf5_read(file_id, group_id, "wt_proc_"//trim(s), wt_send)
                call hdf5_read(file_id, group_id, "fprod_proc_"//trim(s), fprod_send)
                call hdf5_read(file_id, group_id, "eigv_proc_"//trim(s), eigv_send)
                call hdf5_read(file_id, group_id, "eest_proc_"//trim(s), eest_send)
                call hdf5_read(file_id, group_id, "wdsumo_proc_"//trim(s), wdsumo_send)
                call hdf5_read(file_id, group_id, "iage_proc_"//trim(s), iage_send)
                call hdf5_read(file_id, group_id, "ioldest_proc_"//trim(s), ioldest_send)
                call hdf5_read(file_id, group_id, "ioldestmx_proc_"//trim(s), ioldestmx_send)
                call hdf5_read(file_id, group_id, "fratio_proc_"//trim(s), fratio_send)
                if (id .lt. nproc) then
                    call mpi_send(nwalk_send,1,mpi_integer,id,1,MPI_COMM_WORLD,ierr)
                    call mpi_send(xold_dmc_send,3*nelec*MWALK*nforce,mpi_double_precision,id,2,MPI_COMM_WORLD,ierr)
                    call mpi_send(wt_send,MWALK,mpi_double_precision,id,3,MPI_COMM_WORLD,ierr)
                    call mpi_send(ff_send(0),1+nfprod,mpi_double_precision,id,4,MPI_COMM_WORLD,ierr)
                    call mpi_send(fprod_send,1,mpi_double_precision,id,5,MPI_COMM_WORLD,ierr)
                    call mpi_send(fratio_send,MWALK*nforce,mpi_double_precision,id,6,MPI_COMM_WORLD,ierr)
                    call mpi_send(eigv_send,1,mpi_double_precision,id,7,MPI_COMM_WORLD,ierr)
                    call mpi_send(eest_send,1,mpi_double_precision,id,8,MPI_COMM_WORLD,ierr)
                    call mpi_send(wdsumo_send,1,mpi_double_precision,id,9,MPI_COMM_WORLD,ierr)
                    call mpi_send(iage_send,MWALK,mpi_integer,id,10,MPI_COMM_WORLD,ierr)
                    call mpi_send(ioldest_send,1,mpi_integer,id,11,MPI_COMM_WORLD,ierr)
                    call mpi_send(ioldestmx_send,1,mpi_integer,id,12,MPI_COMM_WORLD,ierr)
                endif
            endif
        enddo
        deallocate(xold_dmc_send, wt_send, ff_send, fratio_send, iage_send)
        call hdf5_read(file_id, group_id, "nforce", nforce)

        call hdf5_read(file_id, group_id, "wgcum", wgcum(1:nforce))
        call hdf5_read(file_id, group_id, "egcum", egcum(1:nforce))
        call hdf5_read(file_id, group_id, "pecum_dmc", pecum_dmc(1:nforce))
        call hdf5_read(file_id, group_id, "tpbcum_dmc", tpbcum_dmc(1:nforce))
        call hdf5_read(file_id, group_id, "taucum", taucum(1:nforce))
        call hdf5_read(file_id, group_id, "wgcm2", wgcm2(1:nforce))
        call hdf5_read(file_id, group_id, "egcm2", egcm2(1:nforce))
        call hdf5_read(file_id, group_id, "pecm2_dmc", pecm2_dmc(1:nforce))
        call hdf5_read(file_id, group_id, "tpbcm2_dmc", tpbcm2_dmc(1:nforce))

        call hdf5_read(file_id, group_id, "ipass", ipass)
        call hdf5_read(file_id, group_id, "iblk", iblk)
        call hdf5_read(file_id, group_id, "iblk_proc", iblk_proc)

        tau = taux
        call hdf5_read(file_id, group_id, "rttau", rttau)
        call hdf5_read(file_id, group_id, "idmc", idmc)
        call hdf5_read(file_id, group_id, "wtgen", wtgen(0:nfprod))
        call hdf5_read(file_id, group_id, "wgdsumo", wgdsumo)

        call hdf5_read(file_id, group_id, "wcum_dmc", wcum_dmc)
        call hdf5_read(file_id, group_id, "wfcum", wfcum)
        call hdf5_read(file_id, group_id, "wdcum", wdcum)
        call hdf5_read(file_id, group_id, "wgdcum", wgdcum)
        call hdf5_read(file_id, group_id, "wcum1bynproc", wcum1)
        call hdf5_read(file_id, group_id, "wfcum1bynproc", wfcum1)
        call hdf5_read(file_id, group_id, "wdcum1", wdcum1)
        call hdf5_read(file_id, group_id, "ecum_dmc", ecum_dmc)
        call hdf5_read(file_id, group_id, "efcum", efcum)
        call hdf5_read(file_id, group_id, "ecum1_dmcbynproc", ecum1_dmc)
        call hdf5_read(file_id, group_id, "efcum1bynproc", efcum1)
        call hdf5_read(file_id, group_id, "ei1cum", ei1cum)
        call hdf5_read(file_id, group_id, "ei2cum", ei2cum)
        call hdf5_read(file_id, group_id, "ei3cum", ei3cum)
        call hdf5_read(file_id, group_id, "r2cum_dmc", r2cum_dmc)
        call hdf5_read(file_id, group_id, "ricum", ricum)
        call hdf5_read(file_id, group_id, "wgcum1bynproc", wgcum1(1:nforce))
        call hdf5_read(file_id, group_id, "egcum1bynproc", egcum1(1:nforce))

        call hdf5_read(file_id, group_id, "wcm2", wcm2)
        call hdf5_read(file_id, group_id, "wfcm2", wfcm2)
        call hdf5_read(file_id, group_id, "wdcm2", wdcm2)
        call hdf5_read(file_id, group_id, "wgdcm2", wgdcm2)
        call hdf5_read(file_id, group_id, "wdcm21", wdcm21)
        call hdf5_read(file_id, group_id, "ecm2_dmc", ecm2_dmc)
        call hdf5_read(file_id, group_id, "efcm2", efcm2)

        call hdf5_read(file_id, group_id, "wcm21bynproc", wcm21)
        call hdf5_read(file_id, group_id, "wfcm21bynproc", wfcm21)
        call hdf5_read(file_id, group_id, "wgcm21bynproc", wgcm21(1:nforce))
        call hdf5_read(file_id, group_id, "ecm21_dmcbynproc", ecm21_dmc)
        call hdf5_read(file_id, group_id, "efcm21bynproc", efcm21)
        call hdf5_read(file_id, group_id, "egcm21bynproc", egcm21(1:nforce))
        call hdf5_read(file_id, group_id, "ei1cm2", ei1cm2)
        call hdf5_read(file_id, group_id, "ei2cm2", ei2cm2)
        call hdf5_read(file_id, group_id, "ei3cm2", ei3cm2)
        call hdf5_read(file_id, group_id, "r2cm2_dmc", r2cm2_dmc)
        call hdf5_read(file_id, group_id, "ricm2", ricm2)

        call hdf5_read(file_id, group_id, "fgcum", fgcum(1:nforce))
        call hdf5_read(file_id, group_id, "fgcm2", fgcm2(1:nforce))
        call hdf5_read(file_id, group_id, "derivcum", derivcum)
        call hdf5_read(file_id, group_id, "derivcm2", derivcm2)
        if (allocated(derivtotave_num_old)) &
            call hdf5_read(file_id, group_id, "derivtotave_num_old", derivtotave_num_old(1:nforce))

        if (nrad > 0 .and. allocated(rprob)) then
            call hdf5_read(file_id, group_id, "rprobbynproc", rprob(1:nrad))
            call hdf5_read(file_id, group_id, "rprobup", rprobup(1:nrad))
            call hdf5_read(file_id, group_id, "rprobdn", rprobdn(1:nrad))
        end if
        call hdf5_read(file_id, group_id, "dfus2ac", dfus2ac)
        call hdf5_read(file_id, group_id, "dfus2un", dfus2un)
        call hdf5_read(file_id, group_id, "dr2ac", dr2ac)
        call hdf5_read(file_id, group_id, "dr2un", dr2un)
        call hdf5_read(file_id, group_id, "acc", acc)
        call hdf5_read(file_id, group_id, "trymove", trymove)
        call hdf5_read(file_id, group_id, "nacc", nacc)
        call hdf5_read(file_id, group_id, "nbrnch", nbrnch)
        call hdf5_read(file_id, group_id, "nodecr", nodecr)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: DMC "

        ! properties
        if (iprop.ne.0) then
                call hdf5_group_open(file_id, "Properties", group_id)
                call hdf5_read(file_id, group_id, "iprop", iprop)
                call hdf5_read(file_id, group_id, "nprop", nprop)
                call hdf5_read(file_id, group_id, "vprop_cum", vprop_cum(1:nprop))
                call hdf5_read(file_id, group_id, "vprop_cm2", vprop_cm2(1:nprop))
                call hdf5_group_close(group_id)
                write(ounit, *) " HDF5 Group read :: Properties "
        endif


        ! analytical forces (same content as force_analy_rstrt)
        if (iforce_analy.ne.0) then
            call hdf5_group_open(file_id, "Force Analytical", group_id)
            call hdf5_read(file_id, group_id, "da_energy_cum", da_energy_cum)
            call hdf5_read(file_id, group_id, "da_psi_cum", da_psi_cum)
            call hdf5_read(file_id, group_id, "da_energy_cm2", da_energy_cm2)
            if (dmc_ivd.gt.0) then
                call hdf5_read(file_id, group_id, "da_branch_cum", da_branch_cum)
                call hdf5_read(file_id, group_id, "nwprod", nwprod_stored)
                call hdf5_read(file_id, group_id, "ipathak", ipathak_stored)
                if (nwprod_stored.ne.nwprod) &
                    call fatal_error('DMC_RESTORE_HDF5: nwprod differs from the checkpoint')
                if (ipathak_stored.ne.ipathak) &
                    call fatal_error('DMC_RESTORE_HDF5: ipathak differs from the checkpoint')

                write (unit=s,fmt="(i0)") 0
                call hdf5_read(file_id, group_id, "deriv_eold_proc_"//trim(s), deriv_eold(1:3,1:ncent,1:nwalk))
                call hdf5_read(file_id, group_id, "esnake_proc_"//trim(s), esnake(1:3,1:ncent,1:nwalk,1:PTH))
                if (ipathak.gt.0) &
                    call hdf5_read(file_id, group_id, "pold_proc_"//trim(s), pold(1:nwalk,1:PTH))
                do iph=1,PTH
                    write (unit=s,fmt="(i0,a,i0)") 0, "_", iph
                    call hdf5_read(file_id, group_id, "ehist_proc_"//trim(s), &
                                   ehist(1:3,1:ncent,1:nwalk,0:nwprod-1,iph))
                enddo

                if (nproc .gt. 1) then
                    allocate(deriv_eold_send(3, ncent, MWALK))
                    allocate(esnake_send(3, ncent, MWALK, PTH))
                    allocate(ehist_send(3, ncent, MWALK, 0:nwprod-1, PTH))
                    allocate(pold_send(MWALK, PTH))
                    do id=1, nproc-1
                        nwalk_send = nwalk_all(id)
                        write (unit=s,fmt="(i0)") id
                        call hdf5_read(file_id, group_id, "deriv_eold_proc_"//trim(s), &
                                       deriv_eold_send(1:3,1:ncent,1:nwalk_send))
                        call hdf5_read(file_id, group_id, "esnake_proc_"//trim(s), &
                                       esnake_send(1:3,1:ncent,1:nwalk_send,1:PTH))
                        if (ipathak.gt.0) &
                            call hdf5_read(file_id, group_id, "pold_proc_"//trim(s), pold_send(1:nwalk_send,1:PTH))
                        do iph=1,PTH
                            write (unit=s,fmt="(i0,a,i0)") id, "_", iph
                            call hdf5_read(file_id, group_id, "ehist_proc_"//trim(s), &
                                           ehist_send(1:3,1:ncent,1:nwalk_send,0:nwprod-1,iph))
                        enddo

                        call mpi_send(deriv_eold_send(1:3,1:ncent,1:nwalk_send),3*ncent*nwalk_send, &
                                      mpi_double_precision,id,13,MPI_COMM_WORLD,ierr)
                        call mpi_send(esnake_send(1:3,1:ncent,1:nwalk_send,1:PTH),3*ncent*nwalk_send*PTH, &
                                      mpi_double_precision,id,14,MPI_COMM_WORLD,ierr)
                        call mpi_send(ehist_send(1:3,1:ncent,1:nwalk_send,0:nwprod-1,1:PTH), &
                                      3*ncent*nwalk_send*nwprod*PTH, &
                                      mpi_double_precision,id,15,MPI_COMM_WORLD,ierr)
                        if (ipathak.gt.0) &
                            call mpi_send(pold_send(1:nwalk_send,1:PTH),nwalk_send*PTH, &
                                          mpi_double_precision,id,16,MPI_COMM_WORLD,ierr)
                    enddo
                    deallocate(deriv_eold_send, esnake_send, ehist_send, pold_send)
                endif
            endif
            call hdf5_group_close(group_id)
            write(ounit, *) " HDF5 Group read :: Force Analytical "
        endif

        call hdf5_file_close(file_id)
        ! Close the HDF5 file

        write(ounit, *) ' HDF5 file read successfully :: ', restart_filename
        if (nforce.gt.1) then
          write(ounit,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32 &
            &,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66 &
            &,''fgave'',t79,''(fgerr)'',t93,''npass'',t102,''wgsum'',t112   &
            &,''ioldest'')')
        else
          write(ounit,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32&
            &,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t67&
            &,''npass'',t77,''wgsum'',t85,''ioldest'')')
        endif

        egave_rstrt = egcum(1)/wgcum(1)
        peave_rstrt = pecum_dmc(1)/wgcum(1)
        tpbave_rstrt = tpbcum_dmc(1)/wgcum(1)
        rn_eff = wgcum(1)**2 / wgcm2(1)
        if (rn_eff .gt. 1.d0) then
            egerr_rstrt = dsqrt(max((egcm2(1)/wgcum(1) - egave_rstrt**2)/(rn_eff-1.d0), 0.d0))
            peerr_rstrt = dsqrt(max((pecm2_dmc(1)/wgcum(1) - peave_rstrt**2)/(rn_eff-1.d0), 0.d0))
            tpberr_rstrt = dsqrt(max((tpbcm2_dmc(1)/wgcum(1) - tpbave_rstrt**2)/(rn_eff-1.d0), 0.d0))
        else
            egerr_rstrt = 0.d0
            peerr_rstrt = 0.d0
            tpberr_rstrt = 0.d0
        endif
        if (nforce.gt.1) then
          write(ounit,'(f10.5,3(f10.5,''('',i5,'')''),62x,3i10)') &
              egave_rstrt, egave_rstrt, nint(100000*egerr_rstrt), peave_rstrt, nint(100000*peerr_rstrt), &
              tpbave_rstrt, nint(100000*tpberr_rstrt), iblk_proc*dmc_nstep, nint(wgcum(1)/nproc), ioldest
        else
          write(ounit,'(f10.5,3(f10.5,''('',i5,'')''),3i10)') &
              egave_rstrt, egave_rstrt, nint(100000*egerr_rstrt), peave_rstrt, nint(100000*peerr_rstrt), &
              tpbave_rstrt, nint(100000*tpberr_rstrt), iblk_proc*dmc_nstep, nint(wgcum(1)/nproc), ioldest
        endif

        endif ! master thread (wid)

        if (.not. wid) then
          if (idtask .lt. nproc) then
            call mpi_recv(nwalk,1,mpi_integer,0,1,MPI_COMM_WORLD,istatus,ierr)
            call mpi_recv(xold_dmc,3*nelec*MWALK*nforce,mpi_double_precision,0,2,MPI_COMM_WORLD,istatus,ierr)
            call mpi_recv(wt,MWALK,mpi_double_precision,0,3,MPI_COMM_WORLD,istatus,ierr)
            call mpi_recv(ff(0),1+nfprod,mpi_double_precision,0,4,MPI_COMM_WORLD,istatus,ierr)
            call mpi_recv(fprod,1,mpi_double_precision,0,5,MPI_COMM_WORLD,istatus,ierr)
            call mpi_recv(fratio,MWALK*nforce,mpi_double_precision,0,6,MPI_COMM_WORLD,istatus,ierr)
            call mpi_recv(eigv,1,mpi_double_precision,0,7,MPI_COMM_WORLD,istatus,ierr)
            call mpi_recv(eest,1,mpi_double_precision,0,8,MPI_COMM_WORLD,istatus,ierr)
            call mpi_recv(wdsumo,1,mpi_double_precision,0,9,MPI_COMM_WORLD,istatus,ierr)
            call mpi_recv(iage,MWALK,mpi_integer,0,10,MPI_COMM_WORLD,istatus,ierr)
            call mpi_recv(ioldest,1,mpi_integer,0,11,MPI_COMM_WORLD,istatus,ierr)
            call mpi_recv(ioldestmx,1,mpi_integer,0,12,MPI_COMM_WORLD,istatus,ierr)
            if (iforce_analy.ne.0 .and. dmc_ivd.gt.0) then
              call mpi_recv(deriv_eold(1:3,1:ncent,1:nwalk),3*ncent*nwalk, &
                            mpi_double_precision,0,13,MPI_COMM_WORLD,istatus,ierr)
              call mpi_recv(esnake(1:3,1:ncent,1:nwalk,1:PTH),3*ncent*nwalk*PTH, &
                            mpi_double_precision,0,14,MPI_COMM_WORLD,istatus,ierr)
              call mpi_recv(ehist(1:3,1:ncent,1:nwalk,0:nwprod-1,1:PTH),3*ncent*nwalk*nwprod*PTH, &
                            mpi_double_precision,0,15,MPI_COMM_WORLD,istatus,ierr)
              if (ipathak.gt.0) &
                call mpi_recv(pold(1:nwalk,1:PTH),nwalk*PTH,mpi_double_precision,0,16,MPI_COMM_WORLD,istatus,ierr)
            endif
          endif
        endif

        call bcast(nforce)
        call bcast(nproco)
        call bcast(dmc_nblk)
        call bcast(dmc_nstep)
        call bcast(dmc_nconf)

        call bcast(wgcum(1:nforce))
        call bcast(egcum(1:nforce))
        call bcast(pecum_dmc(1:nforce))
        call bcast(tpbcum_dmc(1:nforce))
        call bcast(taucum(1:nforce))
        call bcast(wgcm2(1:nforce))
        call bcast(egcm2(1:nforce))
        call bcast(pecm2_dmc(1:nforce))
        call bcast(tpbcm2_dmc(1:nforce))
        call bcast(wgcum1(1:nforce))
        call bcast(egcum1(1:nforce))
        call bcast(wgcm21(1:nforce))
        call bcast(egcm21(1:nforce))
        call bcast(fgcum(1:nforce))
        call bcast(fgcm2(1:nforce))

        call bcast(ipass)
        call bcast(iblk)
        call bcast(iblk_proc)

        call bcast(tau)
        call bcast(rttau)
        call bcast(idmc)
        call bcast(wtgen(0:nfprod))
        call bcast(wgdsumo)

        call bcast(wcum_dmc)
        call bcast(wfcum)
        call bcast(wdcum)
        call bcast(wgdcum)
        call bcast(wcum1)
        call bcast(wfcum1)
        call bcast(wdcum1)
        call bcast(ecum_dmc)
        call bcast(efcum)
        call bcast(ecum1_dmc)
        call bcast(efcum1)
        call bcast(ei1cum)
        call bcast(ei2cum)
        call bcast(ei3cum)
        call bcast(r2cum_dmc)
        call bcast(ricum)

        call bcast(wcm2)
        call bcast(wfcm2)
        call bcast(wdcm2)
        call bcast(wgdcm2)
        call bcast(wdcm21)
        call bcast(ecm2_dmc)
        call bcast(efcm2)
        call bcast(wcm21)
        call bcast(wfcm21)
        call bcast(ecm21_dmc)
        call bcast(efcm21)
        call bcast(ei1cm2)
        call bcast(ei2cm2)
        call bcast(ei3cm2)
        call bcast(r2cm2_dmc)
        call bcast(ricm2)

        call bcast(derivcum)
        call bcast(derivcm2)

        call bcast(dfus2ac)
        call bcast(dfus2un)
        call bcast(dr2ac)
        call bcast(dr2un)
        call bcast(acc)
        call bcast(trymove)
        call bcast(nacc)
        call bcast(nbrnch)
        call bcast(nodecr)
        if (.not. wid) then
          acc = zero
          nacc = 0
          trymove = 0
          nodecr = 0
        endif

        call bcast(pecent)
        if (nloc .gt. 0) then
          call bcast(xq(1:nquad))
          call bcast(yq(1:nquad))
          call bcast(zq(1:nquad))
          call bcast(wq(1:nquad))
        endif

        if (iprop .ne. 0) then
          call bcast(vprop_cum(1:nprop))
          call bcast(vprop_cm2(1:nprop))
        endif

        call MPI_Bcast(irn_tmp, 8*nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (idtask .le. nproc-1) call setrn(irn_tmp(1, idtask))

      do iw=1,nwalk
        if(istrech.eq.0) then
          do ifr=2,nforce
            do ie=1,nelec
              do k=1,3
                xold_dmc(k,ie,iw,ifr)=xold_dmc(k,ie,iw,1)
              enddo
            enddo
          enddo
        endif
        do ifr=1,nforce
          if(nforce.gt.1) then
            if(ifr.eq.1.or.istrech.eq.0) then
              call strech(xold_dmc(1,1,iw,1),xold_dmc(1,1,iw,ifr),ajacob,ifr,0)
               else
              call strech(xold_dmc(1,1,iw,1),xold_dmc(1,1,iw,ifr),ajacob,ifr,1)
            endif
           else
            ajacob=one
          endif
          ajacold(iw,ifr)=ajacob
          if(icasula.lt.0) i_vpsp=icasula
          call hpsi(xold_dmc(1,1,iw,ifr),psido_dmc(iw,ifr),psijo_dmc(iw,ifr),ekino,eold(iw,ifr),0,ifr)
          i_vpsp=0
          do i=1,nelec
            call compute_determinante_grad(i,psido_dmc(iw,ifr),psido_dmc(iw,ifr),psijo_dmc(iw,ifr),vold_dmc(1,i,iw,ifr),1)
          enddo
          if(ifr.eq.1) then
            call walksav_det(iw)
            call walksav_jas(iw)
            call t_vpsp_sav
            call prop_save_dmc(iw)
            call pcm_save(iw)
            call mmpol_save(iw)
          endif
        enddo
      enddo

! zero out xsum variables for metrop

      wsum_dmc=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      esum_dmc=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      r2sum=zero
      risum=zero

      do ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum_dmc(ifr)=zero
        tpbsum_dmc(ifr)=zero
        tausum(ifr)=zero
        wsum1(ifr)=zero
        derivsum(:,:,:,ifr)=zero
      enddo

      call prop_init(1)
      call pcm_init(1)
      call mmpol_init(1)

      if(write_walkalize) then
        if(idtask.le.9) then
          write(filename,'(''walkalize.'',i1)') idtask
         elseif(idtask.le.99) then
          write(filename,'(''walkalize.'',i2)') idtask
         elseif(idtask.le.999) then
          write(filename,'(''walkalize.'',i3)') idtask
         else
          call fatal_error('STARTR: idtask > 999')
        endif
        open(unit=11,file=filename,status='old')
        do i=1,2000000000
          read(11,fmt=*,end=100)
        enddo
  100   backspace 11
        backspace 11
      endif

        end subroutine dmc_restore_hdf5
end module dmc_restore_hdf5_mod
