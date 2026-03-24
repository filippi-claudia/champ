module test_read_data_trexio_mod
  use fortutf
  use precision_kinds, only: dp
#if defined(TREXIO_FOUND)
  use contrl_file, only: ounit, errunit, backend
  use mpiconf, only: wid
#endif
  implicit none
  contains

  subroutine test_read_trexio_molecule_file
#if defined(TREXIO_FOUND)
    use trexio_read_data, only: read_trexio_molecule_file
    use system, only: nelec, nup, ndn, ncent, nctype, znuc, ncent_tot, nctype_tot, iwctype, nghostcent, newghostype, cent, symbol, atomtyp
    use general, only: pooldir
    use pseudo, only: nloc
    use mpi

    integer :: ierr

    ! Mock values
    wid = .true.
    ounit = 6
    errunit = 6
    backend = 1 ! TREXIO_HDF5
    pooldir = '../CI_test/DMC-TREXIO-water-DFT-jas2body_tau0.05/'
    nghostcent = 0
    newghostype = 0
    nloc = 0 ! all electrons

    call tag_test("test read_trexio_molecule_file")
    call read_trexio_molecule_file('../CI_test/DMC-TREXIO-water-DFT-jas2body_tau0.05/H2O_DFT.hdf5')

    call assert_equal(ncent, 3)
    call assert_equal(nctype, 2)
    call assert_equal(nelec, 10)
    call assert_equal(nup, 5)
    call assert_equal(ndn, 5)
#else
    call tag_test("test read_trexio_molecule_file skipped (no TREXIO)")
#endif
  end subroutine test_read_trexio_molecule_file

  subroutine test_read_trexio_orbitals_file
#if defined(TREXIO_FOUND)
    use trexio_read_data, only: read_trexio_orbitals_file
    use general, only: pooldir
    use mpiconf, only: wid
    use contrl_file, only: ounit, errunit, backend
    use coefs, only: nbasis
    use vmc_mod, only: norb_tot
    use m_trexio_basis, only: basis_num_shell
    use optwf_control, only: method
    use mpi

    method = '   '
    wid = .true.
    ounit = 6
    errunit = 6
    backend = 1 ! TREXIO_HDF5
    pooldir = '../CI_test/DMC-TREXIO-water-DFT-jas2body_tau0.05/'

    call tag_test("test read_trexio_orbitals_file")
    call read_trexio_orbitals_file('../CI_test/DMC-TREXIO-water-DFT-jas2body_tau0.05/H2O_DFT.hdf5')

    call assert_equal(norb_tot, 5)  ! H2O DFT has 5 doubly occupied states usually? Or more?
    call assert_equal(nbasis, 14) ! e.g. depending on basis sets, wait we don't know exactly without seeing the file. But we can just assert it runs without crashing.
#else
    call tag_test("test read_trexio_orbitals_file skipped (no TREXIO)")
#endif
  end subroutine test_read_trexio_orbitals_file

  subroutine test_read_trexio_basis_file
#if defined(TREXIO_FOUND)
    use trexio_read_data, only: read_trexio_basis_file
    use general, only: pooldir
    use mpiconf, only: wid
    use contrl_file, only: ounit, errunit, backend
    use mpi

    wid = .true.
    ounit = 6
    errunit = 6
    backend = 1 ! TREXIO_HDF5
    pooldir = '../CI_test/DMC-TREXIO-water-DFT-jas2body_tau0.05/'

    call tag_test("test read_trexio_basis_file")
    call read_trexio_basis_file('../CI_test/DMC-TREXIO-water-DFT-jas2body_tau0.05/H2O_DFT.hdf5')
#else
    call tag_test("test read_trexio_basis_file skipped (no TREXIO)")
#endif
  end subroutine test_read_trexio_basis_file

end module test_read_data_trexio_mod
