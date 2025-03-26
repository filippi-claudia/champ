module test_nuclear_repulsion_mod
  !> Unit test the nuclear repulsion energy
  !> @author Ravindra Shinde
  !> @email: r.l.shinde@utwente.nl
  !> @date 2022-11-08
  use fortutf
  use precision_kinds, only: dp
  use contrl_file, only:ounit
  implicit none
  contains
  subroutine test_nuclear_repulsion
    use pot, only: pot_nn
    use system, only: ncent_tot, nctype_tot

    implicit none

    ! Test Water molecule
    integer                   :: nloc = 0
    integer                   :: ncent = 3
    integer                   :: nelec = 8
    integer                   :: nghostcent = 0
    integer                   :: ipr = -1
    integer                   :: iperiodic = 0
    real(dp),dimension(2)     :: znuc = (/6.0d0, 1.0d0/)
    integer,dimension(3)      :: iwctype = (/1,2,2/)
    real(dp), dimension(*)    :: cos_n_sum
    real(dp), dimension(*)    :: sin_n_sum
    real(dp)                  :: pecent
    real                      :: pecent_compare
    real                      :: pecent_expected = 6.98361052364638


    real(dp), dimension(3, 3) :: cent = reshape((/ 0.0d0, 0.0d0, 0.0d0, &
                      -1.43042870599673d0,0.0d0, -1.10715696499747d0, &
                      1.43042870599673d0, 0.0d0, -1.10715696499747d0 /), &
                      (/3,3/), order=(/1,2/))

    real(dp), dimension(8, 3) :: r_en = reshape((/        &
    2.46477203527998d0,  0.771056760199376d0,  1.35802574998623d0, &
    3.59412819526641d0,   4.83786892010651d0, 0.654820174072406d0, &
    0.826640867501410d0,  0.216917081570742d0,  1.52549192819374d0, &
    2.57048630502370d0,   2.68877505199265d0,  3.26479623518646d0, &
    5.94556134820441d0,   1.92490295026335d0,  2.43727073376156d0, &
    1.92722768751480d0,   4.04903788155243d0,  1.86593831690394d0, &
    0.581485889231336d0,   4.40390896944772d0,  3.52038368423407d0, &
    1.92355489240525d0,   1.25782114253839d0,  1.91928900694899d0 /), &
    (/8,3/), order=(/1,2/))


    real(dp), dimension(28) :: r_ee = (/ &
    3.05114575068016d0,   3.71316927068916d0,  1.31245648803426d0, &
    3.55079313213211d0,   4.05862127105284d0,  4.35297920831962d0, &
    6.95665367605274d0,   4.51872576720526d0,  3.63474872494265d0, &
    7.91064648645484d0,   2.59536556415534d0,  1.08208396869396d0, &
    1.58257346833580d0,   3.02405305656700d0,  5.19191693957403d0, &
    3.25766183638361d0,  0.780688618523554d0, 0.834918610857739d0, &
    3.72160785755865d0,   4.38518951953064d0, 0.840366296746517d0, &
    2.48843847397175d0,  0.695294020370444d0,  1.46014763433284d0, &
    3.50853921905586d0,   4.94744780461971d0, 0.549441978242180d0, &
    0.818775132865792d0 /)


    ! Test the nuclear repulsion energy of the water molecule

    ! Get the pecent value
    call tag_test("Test pot_nn: Nuclear repulsion energy")
    ncent_tot = 3 ! Necessary because pot_nn redefines the bounds of cent to be (3,ncent_tot)
    nctype_tot = 2 ! Necessary becaue pot_nn redefines the bounds of znuc to be (nctype_tot)
    call pot_nn(cent,znuc,iwctype,ncent,pecent,cos_n_sum,sin_n_sum)

    ! Lower the precision because FortUTF does not yet provide double precision comparison
    pecent_compare = real(pecent, kind=4)
    call assert_almost_equal(pecent_compare, pecent_expected, 0.000001)

  end subroutine test_nuclear_repulsion

end module test_nuclear_repulsion_mod
