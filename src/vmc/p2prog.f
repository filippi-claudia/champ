C -------------------------- gen2p -------------------------------------
C keywords from file all.p2
C ----------------------------------------------------------------------
      subroutine p2init
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      call p2ini0
      call p2inid
      sfile='all.p2'
      nkey=49
      keys(1)='znuc'
      keylen(1)=4
      nargs(1)=0
      fmts(1)='*'
      keys(2)='lcao'
      keylen(2)=4
      nargs(2)=4
      fmts(2)='iiia'
      keys(3)='geometry'
      keylen(3)=8
      nargs(3)=0
      fmts(3)='*'
      keys(4)='exponents'
      keylen(4)=9
      nargs(4)=1
      fmts(4)='i'
      keys(5)='determinants'
      keylen(5)=12
      nargs(5)=2
      fmts(5)='ii'
      keys(6)='multideterminants'
      keylen(6)=17
      nargs(6)=1
      fmts(6)='i'
      keys(7)='jastrow_parameter'
      keylen(7)=17
      nargs(7)=1
      fmts(7)='i'
      keys(8)='basis'
      keylen(8)=5
      nargs(8)=1
      fmts(8)='i'
      keys(9)='qmc_bf_info'
      keylen(9)=11
      nargs(9)=1
      fmts(9)='i'
      keys(10)='lattice'
      keylen(10)=7
      nargs(10)=0
      fmts(10)='*'
      keys(11)='forces_displace'
      keylen(11)=15
      nargs(11)=0
      fmts(11)='*'
      keys(12)='csf'
      keylen(12)=3
      nargs(12)=3
      fmts(12)='iia'
      keys(13)='csfmap'
      keylen(13)=6
      nargs(13)=1
      fmts(13)='a'
      keys(14)='jasderiv'
      keylen(14)=8
      nargs(14)=0
      fmts(14)='*'
      keys(15)='sym_labels'
      keylen(15)=10
      nargs(15)=3
      fmts(15)='iia'
      keys(16)='optorb_mixvirt'
      keylen(16)=14
      nargs(16)=3
      fmts(16)='iia'
      keys(17)='energies'
      keylen(17)=8
      nargs(17)=2
      fmts(17)='ia'
      keys(18)='eigenvalues'
      keylen(18)=11
      nargs(18)=2
      fmts(18)='ia'
      keys(19)='dmatrix'
      keylen(19)=7
      nargs(19)=3
      fmts(19)='iia'
      keys(20)='cavity_spheres'
      keylen(20)=14
      nargs(20)=1
      fmts(20)='i'
      keys(21)='gradients_cartesian'
      keylen(21)=19
      nargs(21)=0
      fmts(21)='*'
      keys(22)='gradients_zmatrix'
      keylen(22)=17
      nargs(22)=0
      fmts(22)='*'
      keys(23)='modify_zmatrix'
      keylen(23)=14
      nargs(23)=0
      fmts(23)='*'
      keys(24)='hessian_zmatrix'
      keylen(24)=15
      nargs(24)=0
      fmts(24)='*'
      keys(25)='zmatrix_connectionmatrix'
      keylen(25)=24
      nargs(25)=0
      fmts(25)='*'
      keys(26)='efield'
      keylen(26)=6
      nargs(26)=3
      fmts(26)='iia'
      keys(27)='quit'
      keylen(27)=4
      nargs(27)=0
      fmts(27)='*'
      keys(28)='fit_input'
      keylen(28)=9
      nargs(28)=0
      fmts(28)='*'
      keys(29)='array'
      keylen(29)=5
      nargs(29)=2
      fmts(29)='ai'
      keys(30)='vector'
      keylen(30)=6
      nargs(30)=2
      fmts(30)='ai'
      keys(31)='table'
      keylen(31)=5
      nargs(31)=2
      fmts(31)='ai'
      keys(32)='printmacros'
      keylen(32)=11
      nargs(32)=1
      fmts(32)='a'
      keys(33)='savemacros'
      keylen(33)=10
      nargs(33)=1
      fmts(33)='a'
      keys(34)='skipto'
      keylen(34)=6
      nargs(34)=2
      fmts(34)='ai'
      keys(35)='gotol'
      keylen(35)=5
      nargs(35)=1
      fmts(35)='i'
      keys(36)='loop'
      keylen(36)=4
      nargs(36)=4
      fmts(36)='aiii'
      keys(37)='end'
      keylen(37)=3
      nargs(37)=0
      fmts(37)='*'
      keys(38)='?'
      keylen(38)=1
      nargs(38)=0
      fmts(38)='*'
      keys(39)='??'
      keylen(39)=2
      nargs(39)=0
      fmts(39)='*'
      keys(40)='load'
      keylen(40)=4
      nargs(40)=1
      fmts(40)='a'
      keys(41)='fop'
      keylen(41)=3
      nargs(41)=5
      fmts(41)='aadad'
      keys(42)='@@'
      keylen(42)=2
      nargs(42)=5
      fmts(42)='aadad'
      keys(43)='iop'
      keylen(43)=3
      nargs(43)=5
      fmts(43)='aaiai'
      keys(44)='@'
      keylen(44)=1
      nargs(44)=5
      fmts(44)='aaiai'
      keys(45)='info'
      keylen(45)=4
      nargs(45)=2
      fmts(45)='ii'
      keys(46)='finfo'
      keylen(46)=5
      nargs(46)=1
      fmts(46)='i'
      keys(47)='rewind'
      keylen(47)=6
      nargs(47)=1
      fmts(47)='a'
      keys(48)='autorewind'
      keylen(48)=10
      nargs(48)=2
      fmts(48)='ai'
      keys(49)='open_file'
      keylen(49)=9
      nargs(49)=3
      fmts(49)='aaa'
      end
      subroutine p2inid
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      include 'inc/p2defv.inc'
      do i=1,MXKEY
       ideflt(i)=0
       do j=1,MXIDL
        idefpp(j,i)=0
       enddo
      enddo
      ideflt(2)=3
      idefpp(3,2)=1
      idefvv(1)=1
      idefpp(4,2)=1
      adefvv(1)='<input>'
      ideflt(4)=1
      idefpp(1,4)=2
      idefvv(2)=1
      ideflt(5)=2
      idefpp(2,5)=3
      idefvv(3)=1
      ideflt(7)=1
      idefpp(1,7)=4
      idefvv(4)=1
      ideflt(12)=2
      idefpp(2,12)=5
      idefvv(5)=1
      idefpp(3,12)=2
      adefvv(2)='<input>'
      ideflt(13)=1
      idefpp(1,13)=3
      adefvv(3)='<input>'
      ideflt(15)=3
      idefpp(3,15)=4
      adefvv(4)='<input>'
      ideflt(16)=3
      idefpp(3,16)=5
      adefvv(5)='<input>'
      ideflt(17)=2
      idefpp(2,17)=6
      adefvv(6)='<input>'
      ideflt(18)=2
      idefpp(2,18)=7
      adefvv(7)='<input>'
      ideflt(19)=3
      idefpp(3,19)=8
      adefvv(8)='<input>'
      ideflt(26)=3
      idefpp(3,26)=9
      adefvv(9)='<input>'
      ideflt(29)=2
      idefpp(2,29)=6
      idefvv(6)=1
      ideflt(30)=2
      idefpp(2,30)=7
      idefvv(7)=1
      ideflt(31)=2
      idefpp(2,31)=8
      idefvv(8)=1
      ideflt(32)=1
      idefpp(1,32)=10
      adefvv(10)='stdout'
      ideflt(33)=1
      idefpp(1,33)=11
      adefvv(11)='stdout'
      ideflt(34)=2
      idefpp(2,34)=9
      idefvv(9)=1
      ideflt(36)=4
      idefpp(4,36)=10
      idefvv(10)=1
      ideflt(41)=5
      idefpp(5,41)=1
      ddefvv(1)=0
      ideflt(42)=5
      idefpp(5,42)=2
      ddefvv(2)=0
      ideflt(43)=3
      idefpp(3,43)=11
      idefvv(11)=0
      idefpp(4,43)=12
      adefvv(12)='x'
      idefpp(5,43)=12
      idefvv(12)=0
      ideflt(44)=3
      idefpp(3,44)=13
      idefvv(13)=0
      idefpp(4,44)=13
      adefvv(13)='x'
      idefpp(5,44)=14
      idefvv(14)=0
      ideflt(45)=2
      idefpp(2,45)=15
      idefvv(15)=-1
      ideflt(46)=1
      idefpp(1,46)=16
      idefvv(16)=0
      ideflt(48)=1
      idefpp(1,48)=14
      adefvv(14)='on'
      idefpp(2,48)=17
      idefvv(17)=0
      ideflt(49)=2
      idefpp(2,49)=15
      adefvv(15)='f'
      idefpp(3,49)=16
      adefvv(16)='app'
      ip2dfl=1
      end
      subroutine p2call(ikw,itmp,ftmp,is1,is2,lne,iend,MXF,iu)
      implicit double precision (a-h,o-z)
      include 'inc/p2etc.inc'
      dimension itmp(MXF)
      dimension ftmp(MXF)
      character lne*(*)
      dimension is1(MXF)
      dimension is2(MXF)
      iend=0
      goto(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
     $    ,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43
     $    ,44,45,46,47,48,49) ikw
      call fatal('p2call: bad keyword-ID')
 1    continue
       call read_znuc(iu)
      goto 9999
 2    continue
       call read_lcao(itmp(1),itmp(2),itmp(3),lne(is1(1):is2(1)))
      goto 9999
 3    continue
       call read_geometry(iu)
      goto 9999
 4    continue
       call read_exponents(iu,itmp(1))
      goto 9999
 5    continue
       call read_determinants(iu,itmp(1),itmp(2))
      goto 9999
 6    continue
       call read_multideterminants(iu,itmp(1))
      goto 9999
 7    continue
       call read_jastrow_parameter(iu,itmp(1))
      goto 9999
 8    continue
       call read_bas_num_info(iu,itmp(1))
      goto 9999
 9    continue
       call read_bas_num_info(iu,itmp(1))
      goto 9999
 10   continue
       call read_lattice(iu)
      goto 9999
 11   continue
       call read_forces(iu)
      goto 9999
 12   continue
       call read_csf(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 13   continue
       call read_csfmap(lne(is1(1):is2(1)))
      goto 9999
 14   continue
       call read_jasderiv(iu)
      goto 9999
 15   continue
       call read_sym(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 16   continue
       call read_optorb_mixvirt(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 17   continue
       call read_energies(itmp(1),lne(is1(1):is2(1)))
      goto 9999
 18   continue
       call read_energies(itmp(1),lne(is1(1):is2(1)))
      goto 9999
 19   continue
       call read_dmatrix(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 20   continue
       call read_cavity_spheres(iu,itmp(1))
      goto 9999
 21   continue
       call read_gradnts_cart(iu)
      goto 9999
 22   continue
       call read_gradnts_zmat(iu)
      goto 9999
 23   continue
       call read_modify_zmat(iu)
      goto 9999
 24   continue
       call read_hessian_zmat(iu)
      goto 9999
 25   continue
       call read_zmat_conn(iu)
      goto 9999
 26   continue
       call read_efield(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 27   continue
       iend=1
      goto 9999
 28   continue
       iend=1
      goto 9999
 29   continue
       call p2arry(iu,lne(is1(1):is2(1)),0,itmp(1))
      goto 9999
 30   continue
       call p2arry(iu,lne(is1(1):is2(1)),1,itmp(1))
      goto 9999
 31   continue
       call p2arry(iu,lne(is1(1):is2(1)),2,itmp(1))
      goto 9999
 32   continue
       call p2vin(lne(is1(1):is2(1)),0)
      goto 9999
 33   continue
       call p2vin(lne(is1(1):is2(1)),1)
      goto 9999
 34   continue
       call skipto(lne(is1(1):is2(1)),itmp(1))
      goto 9999
 35   continue
       call gotol(itmp(1))
      goto 9999
 36   continue
       call loop(lne(is1(1):is2(1)),itmp(1),itmp(2),itmp(3))
      goto 9999
 37   continue
       call ectrl
      goto 9999
 38   continue
       call cmdlst(3)
      goto 9999
 39   continue
       call cmdlst(0)
      goto 9999
 40   continue
       call ldcmdf(lne(is1(1):is2(1)))
      goto 9999
 41   continue
       call p2fop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),ftmp(1),lne(is
     $    1(3):is2(3)),ftmp(2))
      goto 9999
 42   continue
       call p2fop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),ftmp(1),lne(is
     $    1(3):is2(3)),ftmp(2))
      goto 9999
 43   continue
       call p2iop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),itmp(1),lne(is
     $    1(3):is2(3)),itmp(2))
      goto 9999
 44   continue
       call p2iop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),itmp(1),lne(is
     $    1(3):is2(3)),itmp(2))
      goto 9999
 45   continue
       call infox(itmp(1),itmp(2))
      goto 9999
 46   continue
       call finfo(itmp(1))
      goto 9999
 47   continue
       call rwf(lne(is1(1):is2(1)))
      goto 9999
 48   continue
       call arwnd(lne(is1(1):is2(1)),itmp(1))
      goto 9999
 49   continue
       call fn_open(lne(is1(1):is2(1)),lne(is1(2):is2(2)),lne(is1(3):i
     $    s2(3)))
      goto 9999
 9999 continue
      end 
