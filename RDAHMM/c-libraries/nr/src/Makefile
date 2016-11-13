# Makefile for Recipes: BSD UNIX distribution
#
# $Id: Makefile,v 1.11 1999/07/20 15:30:10 granat Exp granat $
# $Log: Makefile,v $
# Revision 1.11  1999/07/20 15:30:10  granat
# added nr_dludcmp.c
# ,
#
# Revision 1.10  1998/04/22 17:00:51  roden
# Added test to all target dependencies so test progs will be made automatically.
#
# Revision 1.9  1998/04/20 21:01:14  roden
# - Simplified $(LIBRARY) rule so that Gnu make doesn't call ar unnecessarily.
# - Removed %.o rule all together since implicit rule was being used mostly,
#   and explicit rule caused two compilations under certain circumstances.
#
# Revision 1.8  1998/04/20 18:52:43  roden
# Removed $(OUTPUT_OPTION) from %.o compile rule, because on Linux
# this caused extra "-o libnr.a" to be added to compile lines which
# is just wrong.  $(OUTPUT_OPTION) did nothing under Solaris and it
# did the wrong thing under Linux, so removing it is a good thing.
#
# Revision 1.7  1997/07/31 22:25:42  agray
# added makedepend facility.
#
# Revision 1.6  1997/07/31 22:13:11  granat
# revised Makefile to adhere to new "NR_" standard
#
# Revision 1.5  1997/05/15 17:54:15  agray
# change made by Alex Fukunaga, unsanctioned
#
# Revision 1.4  1996/09/24 17:16:51  granat
# Changed object list and dependencies list to use da_ hacked files rather
# than the originals
#
# Revision 1.3  1996/07/25 17:14:52  granat
# added svdcmp3 and nr_util2 to objects list
#
# Revision 1.2  1996/07/25 17:09:29  agray
# added comments.
#
# Note: all ranlib calls were taken out in the port from SunOS to Solaris.
# 7/25/96 AG

# Library, source code, and objects.

LIBRARY		= libnr.a

DESTDIR		= ../../lib

INCDIR		= ../../include

SOURCES.c	= \
		nr_addint.c nr_airy.c nr_amebsa.c nr_amoeba.c nr_amotry.c nr_amotsa.c nr_anneal.c \
		nr_anorm2.c nr_arcmak.c nr_arcode.c nr_arcsum.c nr_asolve.c nr_atimes.c nr_avevar.c \
		nr_balanc.c nr_banbks.c nr_bandec.c nr_banmul.c nr_bcucof.c nr_bcuint.c nr_beschb.c \
		nr_bessi.c nr_bessi0.c nr_bessi1.c nr_bessik.c nr_bessj.c nr_bessj0.c nr_bessj1.c \
		nr_bessjy.c nr_bessk.c nr_bessk0.c nr_bessk1.c nr_bessy.c nr_bessy0.c nr_bessy1.c \
		nr_beta.c nr_betacf.c nr_betai.c nr_bico.c nr_bksub.c nr_bnldev.c nr_brent.c nr_broydn.c \
		nr_bsstep.c nr_caldat.c nr_chder.c nr_chebev.c nr_chebft.c nr_chebpc.c nr_chint.c \
		nr_chixy.c nr_choldc.c nr_cholsl.c nr_chsone.c nr_chstwo.c nr_cisi.c nr_complex.c nr_cntab1.c \
		nr_cntab2.c nr_convlv.c nr_copy.c nr_correl.c nr_cosft1.c nr_cosft2.c nr_covsrt.c \
		nr_crank.c nr_cyclic.c nr_daub4.c nr_dawson.c nr_dbrent.c nr_dcholdc.c nr_dcholsl.c nr_ddpoly.c nr_decchk.c \
		nr_df1dim.c nr_dfour1.c nr_dfpmin.c nr_dfridr.c nr_dftcor.c nr_dftint.c nr_dgammln.c nr_dgammp.c nr_dgasdev.c nr_dgcf.c nr_dgser.c nr_difeq.c \
		nr_djacobi.c nr_dlinmin.c nr_dlubksb.c nr_dludcmp.c nr_dpythag.c nr_dqrdcmp.c nr_drealft.c nr_dsprsax.c nr_dsprstx.c nr_dsvbksb.c \
		nr_dsvdcmp.c nr_eclass.c nr_eclazz.c nr_ei.c nr_eigsrt.c nr_elle.c nr_ellf.c nr_ellpi.c \
		nr_elmhes.c nr_erfcc.c nr_erff.c nr_erffc.c nr_eulsum.c nr_evlmem.c nr_expdev.c \
		nr_expint.c nr_f1dim.c nr_factln.c nr_factrl.c nr_fasper.c nr_fdjac.c nr_fgauss.c \
		nr_fill0.c nr_fit.c nr_fitexy.c nr_fixrts.c nr_fleg.c nr_flmoon.c nr_fmin.c nr_four1.c \
		nr_fourew.c nr_fourfs.c nr_fourn.c nr_fpoly.c nr_fred2.c nr_fredin.c nr_frenel.c \
		nr_frprmn.c nr_ftest.c nr_gamdev.c nr_gammln.c nr_gammp.c nr_gammq.c nr_gasdev.c \
		nr_gaucof.c nr_gauher.c nr_gaujac.c nr_gaulag.c nr_gauleg.c nr_gaussj.c nr_gcf.c \
		nr_golden.c nr_gser.c nr_hpsel.c nr_hpsort.c nr_hqr.c nr_hufapp.c nr_hufdec.c \
		nr_hufenc.c nr_hufmak.c nr_hunt.c nr_hypdrv.c nr_hypgeo.c nr_hypser.c nr_icrc.c \
		nr_icrc1.c nr_igray.c nr_iindexx.c nr_indexx.c nr_interp.c nr_irbit1.c nr_irbit2.c \
		nr_jacobi.c nr_jacobn.c nr_julday.c nr_kendl1.c nr_kendl2.c nr_kermom.c nr_ks2d1s.c \
		nr_ks2d2s.c nr_ksone.c nr_kstwo.c nr_laguer.c nr_lfit.c nr_linbcg.c nr_linmin.c \
		nr_lnsrch.c nr_locate.c nr_lop.c nr_lubksb.c nr_ludcmp.c nr_machar.c nr_matadd.c \
		nr_matsub.c nr_medfit.c nr_memcof.c nr_metrop.c nr_mgfas.c nr_mglin.c nr_midexp.c \
		nr_midinf.c nr_midpnt.c nr_midsql.c nr_midsqu.c nr_miser.c nr_mmid.c nr_mnbrak.c \
		nr_mnewt.c nr_moment.c nr_mp2dfr.c nr_mpdiv.c nr_mpinv.c nr_mpmul.c nr_mpops.c nr_mppi.c \
		nr_mprove.c nr_mpsqrt.c nr_mrqcof.c nr_mrqmin.c nr_newt.c nr_util.c nr_odeint.c nr_orthog.c \
		nr_pade.c nr_pccheb.c nr_pcshft.c nr_pearsn.c nr_period.c nr_piksr2.c nr_piksrt.c \
		nr_pinvs.c nr_plgndr.c nr_poidev.c nr_polcoe.c nr_polcof.c nr_poldiv.c nr_polin2.c \
		nr_polint.c nr_powell.c nr_predic.c nr_probks.c nr_psdes.c nr_pwt.c nr_pwtset.c \
		nr_pythag.c nr_pzextr.c nr_qgaus.c nr_qrdcmp.c nr_qromb.c nr_qromo.c nr_qroot.c \
		nr_qrsolv.c nr_qrupdt.c nr_qsimp.c nr_qtrap.c nr_quad3d.c nr_quadct.c nr_quadmx.c \
		nr_quadvl.c nr_ran0.c nr_ran1.c nr_ran2.c nr_ran3.c nr_ran4.c nr_rank.c nr_ranpt.c \
		nr_ratint.c nr_ratlsq.c nr_ratval.c nr_rc.c nr_rd.c nr_realft.c nr_rebin.c nr_red.c \
		nr_relax.c nr_relax2.c nr_resid.c nr_revcst.c nr_reverse.c nr_rf.c nr_rj.c nr_rk4.c \
		nr_rkck.c nr_rkdumb.c nr_rkqs.c nr_rlft3.c nr_rofunc.c nr_rotate.c nr_rsolv.c \
		nr_rstrct.c nr_rtbis.c nr_rtflsp.c nr_rtnewt.c nr_rtsafe.c nr_rtsec.c nr_rzextr.c \
		nr_savgol.c nr_scrsho.c nr_select.c nr_selip.c nr_shell.c nr_shoot.c nr_shootf.c \
		nr_simp1.c nr_simp2.c nr_simp3.c nr_simplx.c nr_simpr.c nr_sinft.c nr_slvsm2.c \
		nr_slvsml.c nr_sncndn.c nr_snrm.c nr_sobseq.c nr_solvde.c nr_sor.c nr_sort.c nr_sort2.c \
		nr_sort3.c nr_spctrm.c nr_spear.c nr_sphbes.c nr_splie2.c nr_splin2.c nr_spline.c \
		nr_splint.c nr_spread.c nr_sprsax.c nr_sprsin.c nr_sprspm.c nr_sprstm.c nr_sprstp.c \
		nr_sprstx.c nr_stifbs.c nr_stiff.c nr_stoerm.c nr_svbksb.c nr_svdcmp.c \
		nr_svdfit.c nr_svdvar.c nr_toeplz.c nr_tptest.c nr_tqli.c nr_trapzd.c nr_tred2.c nr_tridag.c \
		nr_trncst.c nr_trnspt.c nr_ttest.c nr_tutest.c nr_twofft.c nr_vander.c nr_vegas.c \
		nr_voltra.c nr_wt1.c nr_wtn.c nr_wwghts.c nr_zbrac.c nr_zbrak.c nr_zbrent.c nr_zrhqr.c \
		nr_zriddr.c nr_zroots.c

SOURCES.h	= \
			nr.h nr_util.h nr_complex.h

SOURCES		= $(SOURCES.h) $(SOURCES.c)

OBJECTS		= $(SOURCES.c:%.c=%.o)

# Compiler options.

CC 			= gcc
#CC 			= purify $(NR_CC)
CFLAGS 			= -m64
#CFLAGS 		= -pg $(NR_CFLAGS)
CPPFLAGS	= -I../../include

#If any of the above have problems with an optimization option or other
#compiler flags, a specific rule will be needed for them. This is very
#simple to write. For example, for svdcmp (which sometimes has problems
#with optimizers), a very plain rule is:
#svdcmp.o: nr_svdcmp.c ; $(CC) -I../include -c $<
#Uncomment or write analogously for other recipes as needed.

# Special make commands.

.PRECIOUS:	$(LIBRARY)

# Standard targets.

#all:		depend $(SOURCES) $(OBJECTS) $(LIBRARY) test install
all:		$(SOURCES) $(OBJECTS) $(LIBRARY) test install
objects:	$(SOURCES) $(OBJECTS)
sources:	$(SOURCES)
test:		badluk fredex sfroid sphfpt sphoot

install: 
			cp $(LIBRARY) $(DESTDIR)
			cp libnr.dylib $(DESTDIR)/../dylib
			cp $(SOURCES.h) $(INCDIR)
			cp badluk fredex sfroid sphfpt sphoot $(DESTDIR)/../bin

clean:	
			/bin/rm -f core $(OBJECTS) $(LIBRARY)
			/bin/rm -f badluk fredex sfroid sphfpt sphoot

depend:
#			makedepend -I$(C_CODE_HOME)/include *.c
#			make -I$(C_CODE_HOME)/include *.c

# Targets' dependencies.

$(LIBRARY):	$(OBJECTS)
#			ar rv $@ $?
			libtool -c -static -o $@ $?
			libtool -v -arch_only x86_64 -dynamic -o libnr.dylib $? -lc

# Using implicit .c.o compile rule

badluk: nr_badluk.c $(LIBRARY)
	$(CC) $(CFLAGS) -o badluk nr_badluk.c $(LIBRARY) -lm

fredex: nr_fredex.c $(LIBRARY)
	$(CC) $(CFLAGS) -o fredex nr_fredex.c $(LIBRARY) -lm

sfroid: nr_sfroid.c $(LIBRARY)
	$(CC) $(CFLAGS) -o sfroid nr_sfroid.c $(LIBRARY) -lm

sphfpt: nr_sphfpt.c $(LIBRARY)
	$(CC) $(CFLAGS) -o sphfpt nr_sphfpt.c $(LIBRARY) -lm

sphoot: nr_sphoot.c $(LIBRARY)
	$(CC) $(CFLAGS) -o sphoot nr_sphoot.c $(LIBRARY) -lm

# Specific objects' dependencies.  Use "makedepend *nr_.c".
# DO NOT DELETE THIS LINE -- make depend depends on it.
