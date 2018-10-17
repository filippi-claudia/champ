      subroutine dumper
c Subroutine separate from dumper_more.f since this is the only
c part different in the MPI run
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'pseudo.h'

      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /config/ xold(3,MELEC),xnew(3,MELEC),vold(3,MELEC)
     &,vnew(3,MELEC),psi2o(MSTATES,MFORCE),psi2n(MFORCE),eold(MSTATES,MFORCE),enew(MFORCE)
     &,peo(MSTATES),pen,tjfn,tjfo(MSTATES),psido(MSTATES),psijo
     &,rmino(MELEC),rminn(MELEC),rvmino(3,MELEC),rvminn(3,MELEC)
     &,rminon(MELEC),rminno(MELEC),rvminon(3,MELEC),rvminno(3,MELEC)
     &,nearesto(MELEC),nearestn(MELEC),delttn(MELEC)
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc
      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad

      dimension irn(4)

      call savern(irn)
      rewind 10
c write out number of processors = 1 (serial run)
c it allows to do a parallel restart with nproc > 1 from a serial restart_vmc
      write(10) 1
      write(10) irn
      write(10) nelec,nforce,nloc
      write(10) ((xold(ic,i),ic=1,3),i=1,nelec)
      if(nloc.gt.0) write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)

      call dumper_more
      return
c-----------------------------------------------------------------------
      entry startr

      rewind 10
      read(10) nproc
      read(10) irn
      call setrn(irn)
      read(10) nelecx,nforcex,nlocx
      if (nelecx.ne.nelec) call fatal_error('STARTR: nelec')
      if (nforcex.ne.nforce) call fatal_error('STARTR: nforce')
      if (nlocx.ne.nloc) call fatal_error('STARTR: nloc')
      read(10) ((xold(ic,i),ic=1,3),i=1,nelec)
      if(nloc.gt.0) read(10) nqx,(xq(i),yq(i),zq(i),wq(i),i=1,nqx)
      if(nqx.ne.nquad) call fatal_error('STARTR: nquad')

      call startr_more
      return
      end
