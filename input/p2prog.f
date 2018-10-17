C -------------------------- gen2p -------------------------------------
C keywords from file ../input/all.p2
C ----------------------------------------------------------------------
      subroutine p2init
      implicit double precision (a-h,o-z)
      include 'inc/p2_dim.inc'
      include 'inc/p2.inc'
      call p2ini0
      call p2inid
      sfile='../input/all.p2'
      nkey=29
      keys(1)='znuc'
      keylen(1)=4
      nargs(1)=0
      fmts(1)='*'
      keys(2)='lcao'
      keylen(2)=4
      nargs(2)=3
      fmts(2)='iia'
      keys(3)='geometry'
      keylen(3)=8
      nargs(3)=0
      fmts(3)='*'
      keys(4)='exponents'
      keylen(4)=9
      nargs(4)=0
      fmts(4)='*'
      keys(5)='determinants'
      keylen(5)=12
      nargs(5)=1
      fmts(5)='i'
      keys(6)='jastrow_parameter'
      keylen(6)=17
      nargs(6)=0
      fmts(6)='*'
      keys(7)='basis'
      keylen(7)=5
      nargs(7)=1
      fmts(7)='i'
      keys(8)='fit_input'
      keylen(8)=9
      nargs(8)=0
      fmts(8)='*'
      keys(9)='array'
      keylen(9)=5
      nargs(9)=2
      fmts(9)='ai'
      keys(10)='vector'
      keylen(10)=6
      nargs(10)=2
      fmts(10)='ai'
      keys(11)='table'
      keylen(11)=5
      nargs(11)=2
      fmts(11)='ai'
      keys(12)='printmacros'
      keylen(12)=11
      nargs(12)=1
      fmts(12)='a'
      keys(13)='savemacros'
      keylen(13)=10
      nargs(13)=1
      fmts(13)='a'
      keys(14)='skipto'
      keylen(14)=6
      nargs(14)=2
      fmts(14)='ai'
      keys(15)='gotol'
      keylen(15)=5
      nargs(15)=1
      fmts(15)='i'
      keys(16)='loop'
      keylen(16)=4
      nargs(16)=4
      fmts(16)='aiii'
      keys(17)='end'
      keylen(17)=3
      nargs(17)=0
      fmts(17)='*'
      keys(18)='?'
      keylen(18)=1
      nargs(18)=0
      fmts(18)='*'
      keys(19)='??'
      keylen(19)=2
      nargs(19)=0
      fmts(19)='*'
      keys(20)='load'
      keylen(20)=4
      nargs(20)=1
      fmts(20)='a'
      keys(21)='fop'
      keylen(21)=3
      nargs(21)=5
      fmts(21)='aadad'
      keys(22)='@@'
      keylen(22)=2
      nargs(22)=5
      fmts(22)='aadad'
      keys(23)='iop'
      keylen(23)=3
      nargs(23)=5
      fmts(23)='aaiai'
      keys(24)='@'
      keylen(24)=1
      nargs(24)=5
      fmts(24)='aaiai'
      keys(25)='info'
      keylen(25)=4
      nargs(25)=2
      fmts(25)='ii'
      keys(26)='finfo'
      keylen(26)=5
      nargs(26)=1
      fmts(26)='i'
      keys(27)='rewind'
      keylen(27)=6
      nargs(27)=1
      fmts(27)='a'
      keys(28)='autorewind'
      keylen(28)=10
      nargs(28)=2
      fmts(28)='ai'
      keys(29)='open_file'
      keylen(29)=9
      nargs(29)=3
      fmts(29)='aaa'
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
      adefvv(1)='<input>'
      ideflt(9)=2
      idefpp(2,9)=1
      idefvv(1)=1
      ideflt(10)=2
      idefpp(2,10)=2
      idefvv(2)=1
      ideflt(11)=2
      idefpp(2,11)=3
      idefvv(3)=1
      ideflt(12)=1
      idefpp(1,12)=2
      adefvv(2)='stdout'
      ideflt(13)=1
      idefpp(1,13)=3
      adefvv(3)='stdout'
      ideflt(14)=2
      idefpp(2,14)=4
      idefvv(4)=1
      ideflt(16)=4
      idefpp(4,16)=5
      idefvv(5)=1
      ideflt(21)=5
      idefpp(5,21)=1
      ddefvv(1)=0
      ideflt(22)=5
      idefpp(5,22)=2
      ddefvv(2)=0
      ideflt(23)=3
      idefpp(3,23)=6
      idefvv(6)=0
      idefpp(4,23)=4
      adefvv(4)='x'
      idefpp(5,23)=7
      idefvv(7)=0
      ideflt(24)=3
      idefpp(3,24)=8
      idefvv(8)=0
      idefpp(4,24)=5
      adefvv(5)='x'
      idefpp(5,24)=9
      idefvv(9)=0
      ideflt(25)=2
      idefpp(2,25)=10
      idefvv(10)=-1
      ideflt(26)=1
      idefpp(1,26)=11
      idefvv(11)=0
      ideflt(28)=1
      idefpp(1,28)=6
      adefvv(6)='on'
      idefpp(2,28)=12
      idefvv(12)=0
      ideflt(29)=2
      idefpp(2,29)=7
      adefvv(7)='f'
      idefpp(3,29)=8
      adefvv(8)='app'
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
     $    ,24,25,26,27,28,29) ikw
      call fatal('p2call: bad keyword-ID')
 1    continue
       call getznuc(iu)
      goto 9999
 2    continue
       call getlcao(itmp(1),itmp(2),lne(is1(1):is2(1)))
      goto 9999
 3    continue
       call getgeometry(iu)
      goto 9999
 4    continue
       call getexponents(iu)
      goto 9999
 5    continue
       call getdeterminants(iu,itmp(1))
      goto 9999
 6    continue
       call getjastrow_parameter(iu)
      goto 9999
 7    continue
       call getbasis(iu,itmp(1))
      goto 9999
 8    continue
       iend=1
      goto 9999
 9    continue
       call p2arry(iu,lne(is1(1):is2(1)),0,itmp(1))
      goto 9999
 10   continue
       call p2arry(iu,lne(is1(1):is2(1)),1,itmp(1))
      goto 9999
 11   continue
       call p2arry(iu,lne(is1(1):is2(1)),2,itmp(1))
      goto 9999
 12   continue
       call p2vin(lne(is1(1):is2(1)),0)
      goto 9999
 13   continue
       call p2vin(lne(is1(1):is2(1)),1)
      goto 9999
 14   continue
       call skipto(lne(is1(1):is2(1)),itmp(1))
      goto 9999
 15   continue
       call gotol(itmp(1))
      goto 9999
 16   continue
       call loop(lne(is1(1):is2(1)),itmp(1),itmp(2),itmp(3))
      goto 9999
 17   continue
       call ectrl
      goto 9999
 18   continue
       call cmdlst(3)
      goto 9999
 19   continue
       call cmdlst(0)
      goto 9999
 20   continue
       call ldcmdf(lne(is1(1):is2(1)))
      goto 9999
 21   continue
       call p2fop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),ftmp(1),lne(is
     $    1(3):is2(3)),ftmp(2))
      goto 9999
 22   continue
       call p2fop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),ftmp(1),lne(is
     $    1(3):is2(3)),ftmp(2))
      goto 9999
 23   continue
       call p2iop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),itmp(1),lne(is
     $    1(3):is2(3)),itmp(2))
      goto 9999
 24   continue
       call p2iop(lne(is1(1):is2(1)),lne(is1(2):is2(2)),itmp(1),lne(is
     $    1(3):is2(3)),itmp(2))
      goto 9999
 25   continue
       call infox(itmp(1),itmp(2))
      goto 9999
 26   continue
       call finfo(itmp(1))
      goto 9999
 27   continue
       call rwf(lne(is1(1):is2(1)))
      goto 9999
 28   continue
       call arwnd(lne(is1(1):is2(1)),itmp(1))
      goto 9999
 29   continue
       call fn_open(lne(is1(1):is2(1)),lne(is1(2):is2(2)),lne(is1(3):i
     $    s2(3)))
      goto 9999
 9999 continue
      end 
