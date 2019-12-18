      integer MCHARGES
      parameter(MCHARGES=100)

      common /efield/ iefield,ncharges,iscreen
      common /efield_blk/ xcharge(MCHARGES),ycharge(MCHARGES),zcharge(MCHARGES)
     &                   ,qcharge(MCHARGES),ascreen(MCHARGES),bscreen(MCHARGES)
