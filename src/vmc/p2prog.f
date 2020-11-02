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
      nkey=50
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
      keys(4)='masses'
      keylen(4)=6
      nargs(4)=0
      fmts(4)='*'
      keys(5)='exponents'
      keylen(5)=9
      nargs(5)=1
      fmts(5)='i'
      keys(6)='determinants'
      keylen(6)=12
      nargs(6)=2
      fmts(6)='ii'
      keys(7)='multideterminants'
      keylen(7)=17
      nargs(7)=1
      fmts(7)='i'
      keys(8)='jastrow_parameter'
      keylen(8)=17
      nargs(8)=1
      fmts(8)='i'
      keys(9)='basis'
      keylen(9)=5
      nargs(9)=1
      fmts(9)='i'
      keys(10)='qmc_bf_info'
      keylen(10)=11
      nargs(10)=1
      fmts(10)='i'
      keys(11)='lattice'
      keylen(11)=7
      nargs(11)=0
      fmts(11)='*'
      keys(12)='forces_displace'
      keylen(12)=15
      nargs(12)=0
      fmts(12)='*'
      keys(13)='csf'
      keylen(13)=3
      nargs(13)=3
      fmts(13)='iia'
      keys(14)='csfmap'
      keylen(14)=6
      nargs(14)=1
      fmts(14)='a'
      keys(15)='jasderiv'
      keylen(15)=8
      nargs(15)=0
      fmts(15)='*'
      keys(16)='sym_labels'
      keylen(16)=10
      nargs(16)=3
      fmts(16)='iia'
      keys(17)='optorb_mixvirt'
      keylen(17)=14
      nargs(17)=3
      fmts(17)='iia'
      keys(18)='energies'
      keylen(18)=8
      nargs(18)=2
      fmts(18)='ia'
      keys(19)='eigenvalues'
      keylen(19)=11
      nargs(19)=2
      fmts(19)='ia'
      keys(20)='dmatrix'
      keylen(20)=7
      nargs(20)=3
      fmts(20)='iia'
      keys(21)='cavity_spheres'
      keylen(21)=14
      nargs(21)=1
      fmts(21)='i'
      keys(22)='gradients_cartesian'
      keylen(22)=19
      nargs(22)=0
      fmts(22)='*'
      keys(23)='gradients_zmatrix'
      keylen(23)=17
      nargs(23)=0
      fmts(23)='*'
      keys(24)='modify_zmatrix'
      keylen(24)=14
      nargs(24)=0
      fmts(24)='*'
      keys(25)='hessian_zmatrix'
      keylen(25)=15
      nargs(25)=0
      fmts(25)='*'
      keys(26)='zmatrix_connectionmatrix'
      keylen(26)=24
      nargs(26)=0
      fmts(26)='*'
      keys(27)='efield'
      keylen(27)=6
      nargs(27)=3
      fmts(27)='iia'
      keys(28)='quit'
      keylen(28)=4
      nargs(28)=0
      fmts(28)='*'
      keys(29)='fit_input'
      keylen(29)=9
      nargs(29)=0
      fmts(29)='*'
      keys(30)='array'
      keylen(30)=5
      nargs(30)=2
      fmts(30)='ai'
      keys(31)='vector'
      keylen(31)=6
      nargs(31)=2
      fmts(31)='ai'
      keys(32)='table'
      keylen(32)=5
      nargs(32)=2
      fmts(32)='ai'
      keys(33)='printmacros'
      keylen(33)=11
      nargs(33)=1
      fmts(33)='a'
      keys(34)='savemacros'
      keylen(34)=10
      nargs(34)=1
      fmts(34)='a'
      keys(35)='skipto'
      keylen(35)=6
      nargs(35)=2
      fmts(35)='ai'
      keys(36)='gotol'
      keylen(36)=5
      nargs(36)=1
      fmts(36)='i'
      keys(37)='loop'
      keylen(37)=4
      nargs(37)=4
      fmts(37)='aiii'
      keys(38)='end'
      keylen(38)=3
      nargs(38)=0
      fmts(38)='*'
      keys(39)='?'
      keylen(39)=1
      nargs(39)=0
      fmts(39)='*'
      keys(40)='??'
      keylen(40)=2
      nargs(40)=0
      fmts(40)='*'
      keys(41)='load'
      keylen(41)=4
      nargs(41)=1
      fmts(41)='a'
      keys(42)='fop'
      keylen(42)=3
      nargs(42)=5
      fmts(42)='aadad'
      keys(43)='@@'
      keylen(43)=2
      nargs(43)=5
      fmts(43)='aadad'
      keys(44)='iop'
      keylen(44)=3
      nargs(44)=5
      fmts(44)='aaiai'
      keys(45)='@'
      keylen(45)=1
      nargs(45)=5
      fmts(45)='aaiai'
      keys(46)='info'
      keylen(46)=4
      nargs(46)=2
      fmts(46)='ii'
      keys(47)='finfo'
      keylen(47)=5
      nargs(47)=1
      fmts(47)='i'
      keys(48)='rewind'
      keylen(48)=6
      nargs(48)=1
      fmts(48)='a'
      keys(49)='autorewind'
      keylen(49)=10
      nargs(49)=2
      fmts(49)='ai'
      keys(50)='open_file'
      keylen(50)=9
      nargs(50)=3
      fmts(50)='aaa'
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
      ideflt(5)=1
      idefpp(1,5)=2
      idefvv(2)=1
      ideflt(6)=2
      idefpp(2,6)=3
      idefvv(3)=1
      ideflt(8)=1
      idefpp(1,8)=4
      idefvv(4)=1
      ideflt(13)=2
      idefpp(2,13)=5
      idefvv(5)=1
      idefpp(3,13)=2
      adefvv(2)='<input>'
      ideflt(14)=1
      idefpp(1,14)=3
      adefvv(3)='<input>'
      ideflt(16)=3
      idefpp(3,16)=4
      adefvv(4)='<input>'
      ideflt(17)=3
      idefpp(3,17)=5
      adefvv(5)='<input>'
      ideflt(18)=2
      idefpp(2,18)=6
      adefvv(6)='<input>'
      ideflt(19)=2
      idefpp(2,19)=7
      adefvv(7)='<input>'
      ideflt(20)=3
      idefpp(3,20)=8
      adefvv(8)='<input>'
      ideflt(27)=3
      idefpp(3,27)=9
      adefvv(9)='<input>'
      ideflt(30)=2
      idefpp(2,30)=6
      idefvv(6)=1
      ideflt(31)=2
      idefpp(2,31)=7
      idefvv(7)=1
      ideflt(32)=2
      idefpp(2,32)=8
      idefvv(8)=1
      ideflt(33)=1
      idefpp(1,33)=10
      adefvv(10)='stdout'
      ideflt(34)=1
      idefpp(1,34)=11
      adefvv(11)='stdout'
      ideflt(35)=2
      idefpp(2,35)=9
      idefvv(9)=1
      ideflt(37)=4
      idefpp(4,37)=10
      idefvv(10)=1
      ideflt(42)=5
      idefpp(5,42)=1
      ddefvv(1)=0
      ideflt(43)=5
      idefpp(5,43)=2
      ddefvv(2)=0
      ideflt(44)=3
      idefpp(3,44)=11
      idefvv(11)=0
      idefpp(4,44)=12
      adefvv(12)='x'
      idefpp(5,44)=12
      idefvv(12)=0
      ideflt(45)=3
      idefpp(3,45)=13
      idefvv(13)=0
      idefpp(4,45)=13
      adefvv(13)='x'
      idefpp(5,45)=14
      idefvv(14)=0
      ideflt(46)=2
      idefpp(2,46)=15
      idefvv(15)=-1
      ideflt(47)=1
      idefpp(1,47)=16
      idefvv(16)=0
      ideflt(49)=1
      idefpp(1,49)=14
      adefvv(14)='on'
      idefpp(2,49)=17
      idefvv(17)=0
      ideflt(50)=2
      idefpp(2,50)=15
      adefvv(15)='f'
      idefpp(3,50)=16
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
     $    ,44,45,46,47,48,49,50) ikw
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
       call read_masses(iu)
      goto 9999
 5    continue
       call read_exponents(iu,itmp(1))
      goto 9999
 6    continue
       call read_determinants(iu,itmp(1),itmp(2))
      goto 9999
 7    continue
       call read_multideterminants(iu,itmp(1))
      goto 9999
 8    continue
       call read_jastrow_parameter(iu,itmp(1))
      goto 9999
 9    continue
       call read_bas_num_info(iu,itmp(1))
      goto 9999
 10   continue
       call read_bas_num_info(iu,itmp(1))
      goto 9999
 11   continue
       call read_lattice(iu)
      goto 9999
 12   continue
       call read_forces(iu)
      goto 9999
 13   continue
       call read_csf(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 14   continue
       call read_csfmap(lne(is1(1):is2(1)))
      goto 9999
 15   continue
       call read_jasderiv(iu)
      goto 9999
 16   continue
       call read_sym(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 17   continue
       call read_optorb_mixvirt(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 18   continue
       call read_energies(itmp(1),lne(is1(1):is2(1)))
      goto 9999
 19   continue
       call read_energies(itmp(1),lne(is1(1):is2(1)))
      goto 9999
 20   continue
       call read_dmatrix(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 21   continue
       call read_cavity_spheres(iu,itmp(1))
      goto 9999
 22   continue
       call read_gradnts_cart(iu)
      goto 9999
 23   continue
       call read_gradnts_zmat(iu)
      goto 9999
 24   continue
       call read_modify_zmat(iu)
      goto 9999
 25   continue
       call read_hessian_zmat(iu)
      goto 9999
 26   continue
       call read_zmat_conn(iu)
      goto 9999
 27   continue
       call read_efield(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 28   continue
       iend=1
      goto 9999
 29   continue
       iend=1
      goto 9999
 30   continue
       call p2arry(iu,lne(is1(1):is2(1)),0,itmp(1))
      goto 9999
 31   continue
       call p2arry(iu,lne(is1(1):is2(1)),1,itmp(1))
      goto 9999
 32   continue
       call p2arry(iu,lne(is1(1):is2(1)),2,itmp(1))
      goto 9999
 33   continue
       call p2vin(lne(is1(1):is2(1)),0)
      goto 9999
 34   continue
       call p2vin(lne(is1(1):is2(1)),1)
      goto 9999
 35   continue
       call skipto(lne(is1(1):is2(1)),itmp(1))
      goto 9999
 36   continue
       call gotol(itmp(1))
      goto 9999
 37   continue
       call loop(lne(is1(1):is2(1)),itmp(1),itmp(2),itmp(3))
      goto 9999
 38   continue
       call ectrl
      goto 9999
 39   continue
       call cmdlst(3)
      goto 9999
 40   continue
       call cmdlst(0)
      goto 9999
 41   continue
       call ldcmdf(lne(is1(1):is2(1)))
      goto 9999
 42   continue
       call p2fop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),ftmp(1),lne(is
     $    1(3):is2(3)),ftmp(2))
      goto 9999
 43   continue
       call p2fop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),ftmp(1),lne(is
     $    1(3):is2(3)),ftmp(2))
      goto 9999
 44   continue
       call p2iop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),itmp(1),lne(is
     $    1(3):is2(3)),itmp(2))
      goto 9999
 45   continue
       call p2iop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),itmp(1),lne(is
     $    1(3):is2(3)),itmp(2))
      goto 9999
 46   continue
       call infox(itmp(1),itmp(2))
      goto 9999
 47   continue
       call finfo(itmp(1))
      goto 9999
 48   continue
       call rwf(lne(is1(1):is2(1)))
      goto 9999
 49   continue
       call arwnd(lne(is1(1):is2(1)),itmp(1))
      goto 9999
 50   continue
       call fn_open(lne(is1(1):is2(1)),lne(is1(2):is2(2)),lne(is1(3):i
     $    s2(3)))
      goto 9999
 9999 continue
      end 
