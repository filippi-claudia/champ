        module dmc_store_hdf5_mod
        contains
        subroutine dmc_store_hdf5(restart_filename)
        !> @brief Store the DMC data to a HDF5 file for later restart
        !> @details This subroutine stores DMC the data to a HDF5 file for restarting purposes
        !> @author Ravindra Shinde
        !> @date 2022-11-03
        !> @email r.l.shinde@utwente.nl

        use age,     only: iage,ioldest,ioldestmx
        use basis,   only: ns, np, nd, nf, ng, zex
        use branch,  only: eest,eigv,ff,fprod,nwalk,wdsumo,wgdsumo,wt
        use branch,  only: wtgen
        use coefs,   only: nbasis
        use csfs,    only: ncsf, nstates, ccsf
        use config,  only: xold_dmc
        use constants, only: hb
        use contrl_file, only: ounit
        use contrldmc, only: idmc,nfprod,rttau,tau
        use control, only: mode
        use control_vmc, only: vmc_isite,vmc_nblk, vmc_nstep
        use control_dmc, only: dmc_idump,dmc_irstar,dmc_isite,dmc_nblk, dmc_nstep
        use control_dmc, only: dmc_nblkeq,dmc_nconf,dmc_nconf_new
        use custom_broadcast, only: bcast
        use denupdn, only: rprobdn,rprobup
        use derivest, only: derivcm2,derivcum,derivtotave_num_old
        use dmc_mod, only: MWALK
        use est2cm,  only: ecm21_dmc,ecm2_dmc,efcm2,efcm21,egcm2,egcm21
        use est2cm,  only: ei1cm2,ei2cm2,ei3cm2,pecm2_dmc,r2cm2_dmc,ricm2
        use est2cm,  only: tjfcm_dmc,tpbcm2_dmc,wcm2,wcm21,wdcm2,wdcm21
        use est2cm,  only: wfcm2,wfcm21,wgcm2,wgcm21,wgdcm2
        use estcum,  only: ecum1_dmc,ecum_dmc,efcum,efcum1,egcum,egcum1
        use estcum,  only: ei1cum,ei2cum,ei3cum,iblk,ipass,pecum_dmc
        use estcum,  only: r2cum_dmc,ricum,taucum,tjfcum_dmc,tpbcum_dmc
        use estcum,  only: wcum1,wcum_dmc,wdcum,wdcum1,wfcum,wfcum1,wgcum
        use estcum,  only: wgcum1,wgdcum
        use jacobsave, only: ajacob
        use mmpol,   only: mmpol_dump
        use mpi
        use mpiblk,  only: iblk_proc
        use mpiconf, only: idtask,nproc,wid
        use multiple_geo, only: fgcm2,fgcum,nforce,pecent
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
        use vmc_mod, only: nrad, norb_tot
        use hdf5
        use hdf5_utils, only: hdf5_file_create, hdf5_file_close, hdf5_file_open
        use hdf5_utils, only: hdf5_group_create, hdf5_group_close, hdf5_group_open
        use hdf5_utils, only: hdf5_write!, hdf5_read

        ! properties
        use prp000,  only: iprop,nprop
        use prp003,  only: vprop_cm2,vprop_cum



        implicit none

        ! HDF5 related variables
        character(len=*), intent(in)   ::  restart_filename
        integer(hid_t)                 ::  file_id
        integer(hid_t)                 ::  group_id
        character(len=20)              ::  author = "CHAMP"

        integer :: i, ib, ic, id, ierr
        integer :: ifr, irequest, iw, j
        integer :: k, nscounts
        integer, dimension(4, 0:nproc) :: irn
        integer, dimension(MPI_STATUS_SIZE) :: istatus
        integer, dimension(4, 0:nproc) :: irn_tmp

        real(dp), parameter             :: zero = 0.0d0
        real(dp), parameter             :: one  = 1.0d0
        real(dp), parameter             :: small = 1.0d-6

        if(nforce.gt.1) call strech(xold_dmc,xold_dmc,ajacob,1,0)

        call savern(irn(1,idtask))

        nscounts=4
        call mpi_gather(irn(1,idtask),nscounts,mpi_integer,irn_tmp,nscounts,mpi_integer,0,MPI_COMM_WORLD,ierr)

        ! Bradcast the data which is spread across the processors
        call bcast(nwalk)
        call bcast(xold_dmc)
        call bcast(wt)
        call bcast(ff(0))
        call bcast(fprod)
        call bcast(fratio)
        call bcast(eigv)
        call bcast(eest)
        call bcast(wdsumo)
        call bcast(iage)
        call bcast(ioldest)
        call bcast(ioldestmx)
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
        call hdf5_write(file_id, group_id, "Author", trim(author))
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
        call hdf5_write(file_id, group_id, "Number of Ghost Center Types", newghostype)
        call hdf5_write(file_id, group_id, "Number of Ghost Centers", nghostcent)
        call hdf5_write(file_id, group_id, "Index of Which Center Type", iwctype)
        call hdf5_write(file_id, group_id, "PE Centers", pecent)
        call hdf5_write(file_id, group_id, "Nuclear Charge Znuc", znuc)
        call hdf5_write(file_id, group_id, "Nforce", nforce)
        call hdf5_write(file_id, group_id, "Nloc", nloc)
        call hdf5_write(file_id, group_id, "hb", hb)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: System "

        call hdf5_group_create(file_id, "ECP", group_id)
        call hdf5_group_open(file_id, "ECP", group_id)
        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: ECP "

        call hdf5_group_create(file_id, "Basis", group_id)
        call hdf5_group_open(file_id, "Basis", group_id)
        call hdf5_write(file_id, group_id, "Zex", zex)
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

        call hdf5_group_create(file_id, "DMC", group_id)
        call hdf5_group_open(file_id, "DMC", group_id)
        call hdf5_write(file_id, group_id, "Number of DMC Blocks", dmc_nblk)
        call hdf5_write(file_id, group_id, "Number of DMC Steps per Block", dmc_nstep)
        call hdf5_write(file_id, group_id, "Number of DMC Configurations ", dmc_nconf)
        call hdf5_write(file_id, group_id, "Number of Walkers", nwalk)
        call hdf5_write(file_id, group_id, "Number of Processors", nproc)
        call hdf5_write(file_id, group_id, "xold_dmc", xold_dmc)
        call hdf5_write(file_id, group_id, "nfprod", nfprod)
        call hdf5_write(file_id, group_id, "ff", ff)
        call hdf5_write(file_id, group_id, "wt", wt)
        call hdf5_write(file_id, group_id, "fprod", fprod)
        call hdf5_write(file_id, group_id, "eigv", eigv)
        call hdf5_write(file_id, group_id, "eest", eest)
        call hdf5_write(file_id, group_id, "wdsumo", wdsumo)
        call hdf5_write(file_id, group_id, "iage", iage)
        call hdf5_write(file_id, group_id, "ioldest", ioldest)
        call hdf5_write(file_id, group_id, "ioldestmx", ioldestmx)
        call hdf5_write(file_id, group_id, "nforce", nforce)
        call hdf5_write(file_id, group_id, "fratio", fratio)
        call hdf5_write(file_id, group_id, "wgcum", wgcum(1:nforce))
        call hdf5_write(file_id, group_id, "egcum", egcum(1:nforce))
        call hdf5_write(file_id, group_id, "pecum_dmc", pecum_dmc(1:nforce))
        call hdf5_write(file_id, group_id, "tpbcum_dmc", tpbcum_dmc(1:nforce))
        call hdf5_write(file_id, group_id, "taucum", taucum(1:nforce))

        call hdf5_write(file_id, group_id, "ipass", ipass)
        call hdf5_write(file_id, group_id, "iblk", iblk)
        call hdf5_write(file_id, group_id, "iblk_proc", iblk_proc)

        call hdf5_write(file_id, group_id, "tau", tau)
        call hdf5_write(file_id, group_id, "rttau", rttau)
        call hdf5_write(file_id, group_id, "idmc", idmc)
        call hdf5_write(file_id, group_id, "wtgen", wtgen(0:nfprod))
        call hdf5_write(file_id, group_id, "wgdsumo", wgdsumo)

        call hdf5_write(file_id, group_id, "wcum_dmc", wcum_dmc)
        call hdf5_write(file_id, group_id, "wfcum", wfcum)
        call hdf5_write(file_id, group_id, "wdcum", wdcum)
        call hdf5_write(file_id, group_id, "wgdcum", wgdcum)
        call hdf5_write(file_id, group_id, "wcum1bynproc", wcum1)
        call hdf5_write(file_id, group_id, "wfcum1bynproc", wfcum1)
        call hdf5_write(file_id, group_id, "wdcum1", wdcum1)
        call hdf5_write(file_id, group_id, "ecum_dmc", ecum_dmc)
        call hdf5_write(file_id, group_id, "efcum", efcum)
        call hdf5_write(file_id, group_id, "ecum1_dmcbynproc", ecum1_dmc)
        call hdf5_write(file_id, group_id, "efcum1bynproc", efcum1)
        call hdf5_write(file_id, group_id, "ei1cum", ei1cum)
        call hdf5_write(file_id, group_id, "ei2cum", ei2cum)
        call hdf5_write(file_id, group_id, "ei3cum", ei3cum)
        call hdf5_write(file_id, group_id, "r2cum_dmc", r2cum_dmc)
        call hdf5_write(file_id, group_id, "ricum", ricum)
        call hdf5_write(file_id, group_id, "wgcum1bynproc", wgcum1(1:nforce))
        call hdf5_write(file_id, group_id, "egcum1bynproc", egcum1(1:nforce))

        call hdf5_write(file_id, group_id, "wcm2", wcm2)
        call hdf5_write(file_id, group_id, "wfcm2", wfcm2)
        call hdf5_write(file_id, group_id, "wdcm2", wdcm2)
        call hdf5_write(file_id, group_id, "wgdcm2", wgdcm2)
        call hdf5_write(file_id, group_id, "wdcm21", wdcm21)
        call hdf5_write(file_id, group_id, "ecm2_dmc", ecm2_dmc)
        call hdf5_write(file_id, group_id, "efcm2", efcm2)

        call hdf5_write(file_id, group_id, "wcm21bynproc", wcm21)
        call hdf5_write(file_id, group_id, "wfcm21bynproc", wfcm21)
        call hdf5_write(file_id, group_id, "wgcm21bynproc", wgcm21(1:nforce))
        call hdf5_write(file_id, group_id, "ecm21_dmcbynproc", ecm21_dmc)
        call hdf5_write(file_id, group_id, "efcm21bynproc", efcm21)
        call hdf5_write(file_id, group_id, "egcm21bynproc", egcm21(1:nforce))
        call hdf5_write(file_id, group_id, "ei1cm2", ei1cm2)
        call hdf5_write(file_id, group_id, "ei2cm2", ei2cm2)
        call hdf5_write(file_id, group_id, "ei3cm2", ei3cm2)
        call hdf5_write(file_id, group_id, "r2cm2_dmc", r2cm2_dmc)
        call hdf5_write(file_id, group_id, "ricm2", ricm2)

        call hdf5_write(file_id, group_id, "fgcum", fgcum(1:nforce))
        call hdf5_write(file_id, group_id, "fgcm2", fgcm2(1:nforce))
        call hdf5_write(file_id, group_id, "derivcum", derivcum(1:3,1:nforce))
        call hdf5_write(file_id, group_id, "derivcm2", derivcm2(1:nforce))
        call hdf5_write(file_id, group_id, "derivtotave_num_old", derivtotave_num_old(1:nforce))

        call hdf5_write(file_id, group_id, "rprobbynproc", rprob(1:nrad))
        call hdf5_write(file_id, group_id, "rprobup", rprobup(1:nrad))
        call hdf5_write(file_id, group_id, "rprobdn", rprobdn(1:nrad))
        call hdf5_write(file_id, group_id, "dfus2ac", dfus2ac)
        call hdf5_write(file_id, group_id, "dfus2un", dfus2un)
        call hdf5_write(file_id, group_id, "dr2ac", dr2ac)
        call hdf5_write(file_id, group_id, "dr2un", dr2un)
        call hdf5_write(file_id, group_id, "acc", acc)
        call hdf5_write(file_id, group_id, "trymove", trymove)
        call hdf5_write(file_id, group_id, "nacc", nacc)
        call hdf5_write(file_id, group_id, "nbrnch", nbrnch)
        call hdf5_write(file_id, group_id, "nodecr", nodecr)

        call hdf5_group_close(group_id)
        write(ounit, *) " HDF5 Group saved :: DMC "

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

        call hdf5_file_close(file_id)
        ! Close the HDF5 file
        endif   ! master thread

        end subroutine dmc_store_hdf5
        end module dmc_store_hdf5_mod
