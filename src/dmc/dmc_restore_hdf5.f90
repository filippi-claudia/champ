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
        integer :: ie, ifr, ioldest_id, ioldestmx_id
        integer :: iw, j, k, n1_id
        integer :: n2_id, nbasx, ncentx, nctypex
        integer :: ndetx, ndnx, nelecx, newghostypex
        integer :: nghostcentx, nprock, nq_id, num
        integer :: nupx, nwalk_id
        integer, dimension(4, 0:nproc) :: irn_tmp
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

        character*13 filename

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
        call hdf5_read(file_id, group_id, "iguiding", iguiding)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: States "

        call hdf5_group_open(file_id, "UnitCell", group_id)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: UnitCell "

        call hdf5_group_open(file_id, "Periodic", group_id)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: Periodic "

        call hdf5_group_open(file_id, "QMC", group_id)
        call hdf5_read(file_id, group_id, "Mode", mode)
        call hdf5_read(file_id, group_id, "Number of Processors", nproc)
        call hdf5_read(file_id, group_id, "Random Numbers Each Processor", irn_tmp(1:4,0:nproc-1))
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group read :: QMC "

        call hdf5_group_open(file_id, "DMC", group_id)
        call hdf5_read(file_id, group_id, "Number of DMC Blocks", dmc_nblk)
        call hdf5_read(file_id, group_id, "Number of DMC Steps per Block", dmc_nstep)
        call hdf5_read(file_id, group_id, "Number of DMC Configurations ", dmc_nconf)
        call hdf5_read(file_id, group_id, "Number of Walkers", nwalk)
        call hdf5_read(file_id, group_id, "Number of Processors", nproc)
        call hdf5_read(file_id, group_id, "xold_dmc", xold_dmc)
        call hdf5_read(file_id, group_id, "nfprod", nfprod)
        call hdf5_read(file_id, group_id, "ff", ff)
        call hdf5_read(file_id, group_id, "wt", wt)
        call hdf5_read(file_id, group_id, "fprod", fprod)
        call hdf5_read(file_id, group_id, "eigv", eigv)
        call hdf5_read(file_id, group_id, "eest", eest)
        call hdf5_read(file_id, group_id, "wdsumo", wdsumo)
        call hdf5_read(file_id, group_id, "iage", iage)
        call hdf5_read(file_id, group_id, "ioldest", ioldest)
        call hdf5_read(file_id, group_id, "ioldestmx", ioldestmx)
        call hdf5_read(file_id, group_id, "nforce", nforce)
        call hdf5_read(file_id, group_id, "fratio", fratio)
        call hdf5_read(file_id, group_id, "wgcum", wgcum(1:nforce))
        call hdf5_read(file_id, group_id, "egcum", egcum(1:nforce))
        call hdf5_read(file_id, group_id, "pecum_dmc", pecum_dmc(1:nforce))
        call hdf5_read(file_id, group_id, "tpbcum_dmc", tpbcum_dmc(1:nforce))
        call hdf5_read(file_id, group_id, "taucum", taucum(1:nforce))

        call hdf5_read(file_id, group_id, "ipass", ipass)
        call hdf5_read(file_id, group_id, "iblk", iblk)
        call hdf5_read(file_id, group_id, "iblk_proc", iblk_proc)

        call hdf5_read(file_id, group_id, "tau", tau)
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
        call hdf5_read(file_id, group_id, "derivcum", derivcum(1:3,1:nforce))
        call hdf5_read(file_id, group_id, "derivcm2", derivcm2(1:nforce))
        call hdf5_read(file_id, group_id, "derivtotave_num_old", derivtotave_num_old(1:nforce))

        call hdf5_read(file_id, group_id, "rprobbynproc", rprob(1:nrad))
        call hdf5_read(file_id, group_id, "rprobup", rprobup(1:nrad))
        call hdf5_read(file_id, group_id, "rprobdn", rprobdn(1:nrad))
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


        call hdf5_file_close(file_id)
        ! Close the HDF5 file

      write(ounit,'(t5,''egnow'',t15,''egave'',t21 ,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55 &
                        ,''(tpberr)'' ,t66,''npass'',t77,''wgsum'',t88 ,''ioldest'')')

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
        do k=1,3
          derivsum(k,ifr)=zero
        enddo
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
      endif
  100 backspace 11
      backspace 11

        endif   ! master thread
        end subroutine dmc_restore_hdf5
end module dmc_restore_hdf5_mod
