      integer MCHARGES
      parameter(MCHARGES=100)

      common /efield/ iefield,ncharges,iscreen
      common /efield_blk/ xcharge(MCHARGES),ycharge(MCHARGES),zcharge(MCHARGES)
      common /efield_blk/ qcharge(MCHARGES),ascreen(MCHARGES),bscreen(MCHARGES)
