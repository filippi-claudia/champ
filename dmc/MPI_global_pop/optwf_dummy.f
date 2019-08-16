cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     placeholders for wave function optimization related subroutines
c     Modified by C. Filippi
c     dummy versions for use in DMC code
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine optorb_deriv(psid,denergy)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_compute(eloc)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optjas_deloc(psid,energy,dvpsp_dj,vj)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optci_deloc(eloc_det,e_other,psid,energy)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optjas_sum(wtg_new,wtg_old,enew,eold,iflag)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_define
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_init(iflg)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_save
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_restore
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_cum(wsum,enow)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_sum(p,q,enew,eold,iflag)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_dump(iu)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_rstrt(iu)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_fin(iblk,wcum,ecum)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_prt(w,iblk,etot,iu)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optorb_setup(ndetorb)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optci_setup(ndetorb)
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optci_define
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optci_sum       
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optx_jas_orb_sum       
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optx_jas_ci_sum       
      implicit real*8(a-h,o-z)
      return
      end

      subroutine optx_orb_ci_sum
      implicit real*8(a-h,o-z)
      return
      end

      subroutine do_read_multiple_cistates(iu,ns)
      implicit real*8(a-h,o-z)
      character line*800

      do j=1,ns
       read(iu,*) line
       call incpos(iu,itmp,1)
      enddo
      call p2chkend(iu, 'multiple_cistates')

      end
