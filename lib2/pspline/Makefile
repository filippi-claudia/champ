# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.
# This Makefile should not be used on its own but should instead be called
# from the Makefile in lib2 (to set FC and FFLAGS).

# OBJS = f2test.o f3test.o ibc_ck.o lookup_test.o pspline_calls.o pspltest.o r8akherm1.o r8akherm2.o r8akherm3.o r8bcspeval.o r8bcspgrid.o r8bcspline.o r8bcspvec.o r8bpsplinb.o r8bpspline.o r8cspeval.o r8cspline.o r8cubsplb.o r8dnherm1.o r8dnherm2.o r8dnherm3.o r8evbicub.o r8evspline.o r8evtricub.o r8genxpkg.o r8gridbicub.o r8gridherm1.o r8gridherm2.o r8gridherm3.o r8gridpc1.o r8gridpc2.o r8gridpc3.o r8gridspline.o r8gridtricub.o r8herm1ev.o r8herm2ev.o r8herm3ev.o r8lookup_test.o r8mkbicub.o r8mkbicubw.o r8mkherm1.o r8mkherm2.o r8mkherm3.o r8mkspl2p.o r8mkspl2pb.o r8mkspl2z.o r8mkspl2zb.o r8mkspl3pb.o r8mkspl3zb.o r8mkspline.o r8mktricub.o r8mktricubw.o r8nspline.o r8pc1ev.o r8pc2ev.o r8pc3ev.o r8pspline.o r8seval.o r8seval2.o r8seval3.o r8speval.o r8spgrid.o r8splaan.o r8splbrk.o r8spleen.o r8splinck.o r8spline.o r8spvec.o r8tcspeval.o r8tcspgrid.o r8tcspline.o r8tcspvec.o r8tpsplinb.o r8tpspline.o r8v_spline.o r8vecbicub.o r8vecherm1.o r8vecherm2.o r8vecherm3.o r8vecpc1.o r8vecpc2.o r8vecpc3.o r8vecspline.o r8vectricub.o r8xlookup.o r8zonfind.o v_spline.o

# OBJS = ibc_ck.o r8bcspeval.o r8bcspline.o r8bcspvec.o r8cspeval.o r8cspline.o r8cubsplb.o r8dnherm1.o r8dnherm2.o r8dnherm3.o r8evbicub.o r8evspline.o r8evtricub.o r8genxpkg.o r8herm1ev.o r8herm2ev.o r8herm3ev.o r8mkbicubw.o r8mkherm1.o r8mkherm2.o r8mkherm3.o r8mkspl2p.o r8mkspl2pb.o r8mkspl2z.o r8mkspl2zb.o r8mkspl3pb.o r8mkspl3zb.o r8mkspline.o r8mktricubw.o r8nspline.o r8pc1ev.o r8pc2ev.o r8pc3ev.o r8speval.o r8splaan.o r8splbrk.o r8spleen.o r8splinck.o r8spline.o r8spvec.o r8tcspeval.o r8tcspline.o r8tcspvec.o r8v_spline.o r8vecbicub.o r8vecherm1.o r8vecherm2.o r8vecherm3.o r8vecpc1.o r8vecpc2.o r8vecpc3.o r8vecspline.o r8vectricub.o r8zonfind.o

OBJS = ibc_ck.o r8tcspline.o r8mktricubw.o r8tcspeval.o r8splinck.o \
r8v_spline.o r8bcspline.o r8zonfind.o r8bcspeval.o r8fvtricub.o \
tcspline.o mktricubw.o tcspeval.o splinck.o \
v_spline.o bcspline.o zonfind.o bcspeval.o fvtricub.o 

libpspline.a: $(OBJS)
	ar rs $@ $(OBJS)

clean:
	- rm -f $(OBJS)
	- rm -f libpspline.a
