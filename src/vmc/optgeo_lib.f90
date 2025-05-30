module optgeo_lib
contains
      subroutine write_geometry(iter)

      use system, only: cent, iwctype, nctype, ncent
      use mpiconf, only: wid

      implicit none

      integer :: i, index, iter, k
      character(len=40) filename,itn

      if (.not.wid) return

      if(iter.lt.0) then
        filename='geo_optimal_final'
       else
        if(iter.lt.10) then
          write(itn,'(i1)') iter
         elseif(iter.lt.100) then
          write(itn,'(i2)') iter
         elseif(iter.lt.1000) then
          write(itn,'(i3)') iter
        endif
        filename='geo_optimal.iter'//itn(1:index(itn,' ')-1)
      endif

      open(2,file=filename,status='unknown')

      write(2,'(''# geometry iter '',i5)') iter
      write(2,'(''&atoms nctype '',i5,'' natom '',i5)') nctype,ncent
      write(2,'(''geometry'')')
      do i=1,ncent
        write(2,'(3f14.6,i5)') (cent(k,i),k=1,3),iwctype(i)
      enddo
      write(2,'(''end'')')

      close(2)

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine compute_positions

        use coords_int
        use system, only: cent, ncent
        use zmatrix, only: czint, izcmat
        use m_force_analytic, only: iforce_analy, iuse_zmat, alfgeo, da_energy_ave
        use contrl_file,    only: ounit
      implicit none

      integer :: ic, k

        if (iforce_analy.eq.0) return

        call compute_position_bcast

        if(iuse_zmat.eq.1) then
          call coords_init (ncent, cent, izcmat)
          call coords_compute_wilson (cent, izcmat)
          call coords_transform_gradients (da_energy_ave(:,:,1))
          call coords_compute_step (alfgeo)
          call coords_transform_step (czint, cent, izcmat)

          write (ounit,*) 'INTERNAL'
          do ic=1,ncent
            write (ounit,'(1x, 3f10.5)') czint(1:3, ic)
          enddo

          write (ounit,*) 'CENT'
          do ic=1,ncent
            write(ounit,'(3f10.5)') (cent(k,ic),k=1,3)
          enddo

        else
          do ic=1,ncent
            do k=1,3
              cent(k,ic)=cent(k,ic)-alfgeo*da_energy_ave(k,ic,1)
            enddo
            write(ounit,*)'CENT ',(cent(k,ic),k=1,3)
          enddo
        endif

        return
      end

      subroutine compute_position_bcast

      use system, only: ncent
      use m_force_analytic, only: iforce_analy, da_energy_ave
      use mpi
      use system,  only: ncent

      implicit none

      integer :: i

      if(iforce_analy.eq.0)return

      call MPI_BCAST(da_energy_ave,3*ncent,MPI_REAL8,0,MPI_COMM_WORLD,i)

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine force_store(l)

      use system, only: ncent
      use da_energy_now, only: da_energy, da_psi
      use m_force_analytic, only: force_o

      implicit none

      integer :: i, ii, k, l

      ii=0
      do i=1,ncent
        do k=1,3
          ii=ii+1
          force_o(ii,l)=da_psi(k,i)
        enddo
      enddo

      do i=1,ncent
        do k=1,3
          ii=ii+1
          force_o(ii,l)=da_energy(k,i)
        enddo
      enddo

      return
      end
end module 
