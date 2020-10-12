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
      nkey=47
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
      keys(8)='lattice'
      keylen(8)=7
      nargs(8)=0
      fmts(8)='*'
      keys(9)='forces_displace'
      keylen(9)=15
      nargs(9)=0
      fmts(9)='*'
      keys(10)='csf'
      keylen(10)=3
      nargs(10)=3
      fmts(10)='iia'
      keys(11)='csfmap'
      keylen(11)=6
      nargs(11)=1
      fmts(11)='a'
      keys(12)='jasderiv'
      keylen(12)=8
      nargs(12)=0
      fmts(12)='*'
      keys(13)='sym_labels'
      keylen(13)=10
      nargs(13)=3
      fmts(13)='iia'
      keys(14)='optorb_mixvirt'
      keylen(14)=14
      nargs(14)=3
      fmts(14)='iia'
      keys(15)='energies'
      keylen(15)=8
      nargs(15)=2
      fmts(15)='ia'
      keys(16)='eigenvalues'
      keylen(16)=11
      nargs(16)=2
      fmts(16)='ia'
      keys(17)='dmatrix'
      keylen(17)=7
      nargs(17)=3
      fmts(17)='iia'
      keys(18)='cavity_spheres'
      keylen(18)=14
      nargs(18)=1
      fmts(18)='i'
      keys(19)='gradients_cartesian'
      keylen(19)=19
      nargs(19)=0
      fmts(19)='*'
      keys(20)='gradients_zmatrix'
      keylen(20)=17
      nargs(20)=0
      fmts(20)='*'
      keys(21)='modify_zmatrix'
      keylen(21)=14
      nargs(21)=0
      fmts(21)='*'
      keys(22)='hessian_zmatrix'
      keylen(22)=15
      nargs(22)=0
      fmts(22)='*'
      keys(23)='zmatrix_connectionmatrix'
      keylen(23)=24
      nargs(23)=0
      fmts(23)='*'
      keys(24)='efield'
      keylen(24)=6
      nargs(24)=3
      fmts(24)='iia'
      keys(25)='quit'
      keylen(25)=4
      nargs(25)=0
      fmts(25)='*'
      keys(26)='fit_input'
      keylen(26)=9
      nargs(26)=0
      fmts(26)='*'
      keys(27)='array'
      keylen(27)=5
      nargs(27)=2
      fmts(27)='ai'
      keys(28)='vector'
      keylen(28)=6
      nargs(28)=2
      fmts(28)='ai'
      keys(29)='table'
      keylen(29)=5
      nargs(29)=2
      fmts(29)='ai'
      keys(30)='printmacros'
      keylen(30)=11
      nargs(30)=1
      fmts(30)='a'
      keys(31)='savemacros'
      keylen(31)=10
      nargs(31)=1
      fmts(31)='a'
      keys(32)='skipto'
      keylen(32)=6
      nargs(32)=2
      fmts(32)='ai'
      keys(33)='gotol'
      keylen(33)=5
      nargs(33)=1
      fmts(33)='i'
      keys(34)='loop'
      keylen(34)=4
      nargs(34)=4
      fmts(34)='aiii'
      keys(35)='end'
      keylen(35)=3
      nargs(35)=0
      fmts(35)='*'
      keys(36)='?'
      keylen(36)=1
      nargs(36)=0
      fmts(36)='*'
      keys(37)='??'
      keylen(37)=2
      nargs(37)=0
      fmts(37)='*'
      keys(38)='load'
      keylen(38)=4
      nargs(38)=1
      fmts(38)='a'
      keys(39)='fop'
      keylen(39)=3
      nargs(39)=5
      fmts(39)='aadad'
      keys(40)='@@'
      keylen(40)=2
      nargs(40)=5
      fmts(40)='aadad'
      keys(41)='iop'
      keylen(41)=3
      nargs(41)=5
      fmts(41)='aaiai'
      keys(42)='@'
      keylen(42)=1
      nargs(42)=5
      fmts(42)='aaiai'
      keys(43)='info'
      keylen(43)=4
      nargs(43)=2
      fmts(43)='ii'
      keys(44)='finfo'
      keylen(44)=5
      nargs(44)=1
      fmts(44)='i'
      keys(45)='rewind'
      keylen(45)=6
      nargs(45)=1
      fmts(45)='a'
      keys(46)='autorewind'
      keylen(46)=10
      nargs(46)=2
      fmts(46)='ai'
      keys(47)='open_file'
      keylen(47)=9
      nargs(47)=3
      fmts(47)='aaa'
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
      ideflt(10)=2
      idefpp(2,10)=5
      idefvv(5)=1
      idefpp(3,10)=2
      adefvv(2)='<input>'
      ideflt(11)=1
      idefpp(1,11)=3
      adefvv(3)='<input>'
      ideflt(13)=3
      idefpp(3,13)=4
      adefvv(4)='<input>'
      ideflt(14)=3
      idefpp(3,14)=5
      adefvv(5)='<input>'
      ideflt(15)=2
      idefpp(2,15)=6
      adefvv(6)='<input>'
      ideflt(16)=2
      idefpp(2,16)=7
      adefvv(7)='<input>'
      ideflt(17)=3
      idefpp(3,17)=8
      adefvv(8)='<input>'
      ideflt(24)=3
      idefpp(3,24)=9
      adefvv(9)='<input>'
      ideflt(27)=2
      idefpp(2,27)=6
      idefvv(6)=1
      ideflt(28)=2
      idefpp(2,28)=7
      idefvv(7)=1
      ideflt(29)=2
      idefpp(2,29)=8
      idefvv(8)=1
      ideflt(30)=1
      idefpp(1,30)=10
      adefvv(10)='stdout'
      ideflt(31)=1
      idefpp(1,31)=11
      adefvv(11)='stdout'
      ideflt(32)=2
      idefpp(2,32)=9
      idefvv(9)=1
      ideflt(34)=4
      idefpp(4,34)=10
      idefvv(10)=1
      ideflt(39)=5
      idefpp(5,39)=1
      ddefvv(1)=0
      ideflt(40)=5
      idefpp(5,40)=2
      ddefvv(2)=0
      ideflt(41)=3
      idefpp(3,41)=11
      idefvv(11)=0
      idefpp(4,41)=12
      adefvv(12)='x'
      idefpp(5,41)=12
      idefvv(12)=0
      ideflt(42)=3
      idefpp(3,42)=13
      idefvv(13)=0
      idefpp(4,42)=13
      adefvv(13)='x'
      idefpp(5,42)=14
      idefvv(14)=0
      ideflt(43)=2
      idefpp(2,43)=15
      idefvv(15)=-1
      ideflt(44)=1
      idefpp(1,44)=16
      idefvv(16)=0
      ideflt(46)=1
      idefpp(1,46)=14
      adefvv(14)='on'
      idefpp(2,46)=17
      idefvv(17)=0
      ideflt(47)=2
      idefpp(2,47)=15
      adefvv(15)='f'
      idefpp(3,47)=16
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
     $    ,44,45,46,47) ikw
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
       call read_lattice(iu)
      goto 9999
 9    continue
       call read_forces(iu)
      goto 9999
 10   continue
       call read_csf(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 11   continue
       call read_csfmap(lne(is1(1):is2(1)))
      goto 9999
 12   continue
       call read_jasderiv(iu)
      goto 9999
 13   continue
       call read_sym(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 14   continue
       call read_optorb_mixvirt(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 15   continue
       call read_energies(itmp(1),lne(is1(1):is2(1)))
      goto 9999
 16   continue
       call read_energies(itmp(1),lne(is1(1):is2(1)))
      goto 9999
 17   continue
       call read_dmatrix(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 18   continue
       call read_cavity_spheres(iu,itmp(1))
      goto 9999
 19   continue
       call read_gradnts_cart(iu)
      goto 9999
 20   continue
       call read_gradnts_zmat(iu)
      goto 9999
 21   continue
       call read_modify_zmat(iu)
      goto 9999
 22   continue
       call read_hessian_zmat(iu)
      goto 9999
 23   continue
       call read_zmat_conn(iu)
      goto 9999
 24   continue
       call read_efield(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 25   continue
       iend=1
      goto 9999
 26   continue
       iend=1
      goto 9999
 27   continue
       call p2arry(iu,lne(is1(1):is2(1)),0,itmp(1))
      goto 9999
 28   continue
       call p2arry(iu,lne(is1(1):is2(1)),1,itmp(1))
      goto 9999
 29   continue
       call p2arry(iu,lne(is1(1):is2(1)),2,itmp(1))
      goto 9999
 30   continue
       call p2vin(lne(is1(1):is2(1)),0)
      goto 9999
 31   continue
       call p2vin(lne(is1(1):is2(1)),1)
      goto 9999
 32   continue
       call skipto(lne(is1(1):is2(1)),itmp(1))
      goto 9999
 33   continue
       call gotol(itmp(1))
      goto 9999
 34   continue
       call loop(lne(is1(1):is2(1)),itmp(1),itmp(2),itmp(3))
      goto 9999
 35   continue
       call ectrl
      goto 9999
 36   continue
       call cmdlst(3)
      goto 9999
 37   continue
       call cmdlst(0)
      goto 9999
 38   continue
       call ldcmdf(lne(is1(1):is2(1)))
      goto 9999
 39   continue
       call p2fop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),ftmp(1),lne(is
     $    1(3):is2(3)),ftmp(2))
      goto 9999
 40   continue
       call p2fop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),ftmp(1),lne(is
     $    1(3):is2(3)),ftmp(2))
      goto 9999
 41   continue
       call p2iop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),itmp(1),lne(is
     $    1(3):is2(3)),itmp(2))
      goto 9999
 42   continue
       call p2iop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),itmp(1),lne(is
     $    1(3):is2(3)),itmp(2))
      goto 9999
 43   continue
       call infox(itmp(1),itmp(2))
      goto 9999
 44   continue
       call finfo(itmp(1))
      goto 9999
 45   continue
       call rwf(lne(is1(1):is2(1)))
      goto 9999
 46   continue
       call arwnd(lne(is1(1):is2(1)),itmp(1))
      goto 9999
 47   continue
       call fn_open(lne(is1(1):is2(1)),lne(is1(2):is2(2)),lne(is1(3):i
     $    s2(3)))
      goto 9999
 9999 continue
      end 
