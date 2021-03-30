      subroutine optci_reduce
      use optci, only: MXCITERM, MXCIREDUCED, MXCIMATDIM
      use optwf_contrl, only: ioptci
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum
      use ci000, only: nciterm
      use ci005_blk, only: ci_o_cum
      use ci006_blk, only: ci_de_cum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum
      use ci010_blk, only: ci_ooe_cum
      use method_opt, only: method
      use csfs, only: nstates
      use mpi

      implicit real*8(a-h,o-z)

      parameter(MXTMP=MXCITERM+MXCIMATDIM)
      dimension collect(MXTMP),collect2(MXCITERM,MXCIREDUCED)

      if (iefficiency.gt.0) then
         call mpi_reduce(effcum,collect,nstates_psig
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         call mpi_bcast(collect,nstates_psig
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

         do j=1,nstates_psig
            effcum(j)=collect(j)
         enddo

         call mpi_reduce(effcm2,collect,nstates_psig
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         call mpi_bcast(collect,nstates_psig
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

         do j=1,nstates_psig
            effcm2(j)=collect(j)
         enddo
      endif

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates
         call mpi_reduce(ci_o_cum(1,istate),collect(1),nciterm
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         call mpi_bcast(collect,nciterm
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

         do i=1,nciterm
            ci_o_cum(i,istate)=collect(i)
         enddo

         call mpi_reduce(ci_de_cum(1,istate),collect(1),nciterm
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         call mpi_bcast(collect,nciterm
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

         do i=1,nciterm
            ci_de_cum(i,istate)=collect(i)
         enddo

         call mpi_reduce(ci_oe_cum(1,1,istate),collect2,MXCITERM*nciterm
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         call mpi_bcast(collect2,MXCITERM*nciterm
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

         do i=1,nciterm
            do j=1,nciterm
               ci_oe_cum(i,j,istate) = collect2(i,j)
            enddo
         enddo

         call mpi_reduce(ci_oe_cm2(1,1,istate),collect2,MXCITERM*nciterm
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         call mpi_bcast(collect2,MXCITERM*nciterm
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

         do i=1,nciterm
            do j=1,nciterm
               ci_oe_cm2(i,j,istate)=collect2(i,j)
            enddo
         enddo

         matdim=nciterm*(nciterm+1)/2

         call mpi_reduce(ci_oo_cum(1,istate),collect(1),matdim
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         call mpi_bcast(collect,matdim
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

         do i=1,matdim
            ci_oo_cum(i,istate)=collect(i)
         enddo

         call mpi_reduce(ci_oo_cm2(1,istate),collect(1),matdim
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         call mpi_bcast(collect,matdim
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

         do i=1,matdim
            ci_oo_cm2(i,istate)=collect(i)
         enddo

         call mpi_reduce(ci_ooe_cum(1,istate),collect(1),matdim
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         call mpi_bcast(collect,matdim
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

         do i=1,matdim
            ci_ooe_cum(i,istate)=collect(i)
         enddo
      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      end subroutine
