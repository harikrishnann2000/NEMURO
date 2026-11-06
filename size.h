ccccccccccccccccccccccccccc include size.h cccccccccccccccccccccccccccccccc
c
c	dimension definitions for offline model
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	parameter (imt=1,jmt=1,km=51,lm=13, NPP=2, nn=2+15, nproc=96,
     &     noutput=1, nvar3d=5)


!!!!!!!!!!!!!!!! Double precision Start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!5-D	
	double precision  t , buffer_tr

!4-D	
	double precision u,v,w

!3-D
	double precision jpo4_3d, jdop_3d, jo2_3d, jdic_3d 
	double precision jalk_3d, jrrcv_3d
	double precision rbgc_Sb_po4
	double precision rbgc_Sb_dic, rbgc_Sb_dop, rbgc_Sb_Alk
	double precision buffer_mpp_top, buffer_mpp_bot

!2-D	
	double precision rimpl_vmix, rCFC_save
	double precision dzfact, dzfact_rev
	double precision rkpp_nonlocal
	double precision ptb
	double precision co2_flux, alk_flux
	
!1-D
	double precision rt_mix, rlamda 
	double precision adv_vn_iso, adv_vs_iso, adv_ve_iso
	double precision adv_vw_iso, adv_vt_iso, adv_fe_iso, adv_vb_iso
        double precision adv_fw_iso, adv_fn_iso, adv_fs_iso
	double precision rK13_east, rK13_west, rK23_north, rK23_south
	double precision adv_ft_iso, adv_fb_iso
	double precision rD, rD1, rD2, rDDu, rDDl, biotic_source_sink
	double precision rsum_volint, tr_mix_coeff
	double precision adv_vt_bup, adv_vb_bup
        double precision adv_ve_bup, adv_vw_bup
        double precision adv_vn_bup, adv_vs_bup
	double precision dz, dzb
	double precision jpo4, jdop, jo2, jdic, jalk, ppo4star, jrrcv
	double precision rbgc_Sb, rbgc_Sb_Ca
	double precision rbgc_B, rbgc_B_Ca

!5-D	
	common /size/ t(imt,km,jmt,nn,-1:1)
	common /size/ buffer_tr(imt, km, nproc, 4, 2)

!4-D
	common /size/ u(imt,km,jmt,-1:1)
	common /size/ v(imt,km,jmt,-1:1)
	common /size/ w(imt,km,jmt,-1:1)
	
!3-D
	common /size/ jpo4_3d(imt,km,jmt), jdop_3d(imt,km,jmt)
	common /size/ jo2_3d(imt,km,jmt), jdic_3d(imt,km,jmt)
	common /size/ jalk_3d(imt,km,jmt), jrrcv_3d(imt,km,jmt)
	common /szie/ rbgc_Sb_po4(imt,km,jmt), rbgc_Sb_dic(imt,km,jmt)
	common /szie/ rbgc_Sb_dop(imt,km,jmt), rbgc_Sb_Alk(imt,km,jmt)
	common /size/ buffer_mpp_top(imt,km,2)
        common /size/ buffer_mpp_bot(imt,km,2)

!2-D
	common /size/ rimpl_vmix(km, km)
	common /size/ rCFC_save(imt,jmt)
	common /size/ dzfact(km, km), dzfact_rev(km, km)
	common /size/ rkpp_nonlocal(km,jmt)
	common /size/ ptb(km,5)
	common /size/ co2_flux(imt, jmt), alk_flux(imt, jmt)
!1-D
	
	common /size/ rt_mix(km), rlamda(km)
	common /size/ adv_vn_iso(km), adv_vs_iso(km)
        common /size/ adv_ve_iso(km), adv_vw_iso(km)
        common /size/ adv_vt_iso(km), adv_vb_iso(km)
        common /size/ adv_ve_bup(km), adv_vw_bup(km)
        common /size/ adv_vn_bup(km), adv_vs_bup(km)
        common /size/ adv_vt_bup(km), adv_vb_bup(km)
        common /size/ adv_fe_iso(km), adv_fw_iso(km)
        common /size/ adv_fn_iso(km), adv_fs_iso(km)
        common /size/ adv_ft_iso(km), adv_fb_iso(km)
	common /size/ rK13_east(0:km),rk13_west(0:km)
        common /size/ rK23_north(0:km),rk23_south(0:km)
	common /size/ rD(km), rD1(km), rD2(km)
        common /size/ rDDu(km), rDDl(km), biotic_source_sink(km)
	common /size/ rsum_volint(nproc)
	common /size/ tr_mix_coeff(0:km)
	common /size/ dz(0:km), dzb(0:km)
	common /size/ jpo4(km), jdop(km)
        common /size/ jo2(km), jdic(km), jalk(km), jrrcv(km)
        common /size/ ppo4star(km), rbgc_Sb(km), rbgc_Sb_Ca(km)
	common /size/ rbgc_B(km), rbgc_B_Ca(km)

!!!!!!!!!!!!!!!! Double precision end !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!! Real as well as undeclared but alligned
!3-D
	real umask, vmask, tmask, lmask,bmask
	real muv_read, mld_read, dic_read, alk_read, tau_read, eta
	real mask
!2-D
	real mld, tau,  muv
	real kmask, max_mask

!5-D
	common /size/ data_offline(imt,km,jmt,5,2)
!4-D
	common /size/ var3d_recv(imt,km,jmt,nvar3d)
	common /size/ var3d_send(imt,km,jmt,nvar3d)
	common /size/ var3d_mpp(imt, km, jmt,noutput)
	common /size/ phosp_read(imt, km, jmt,2)
	common /size/ phospmean_read(imt, km, jmt, 2)

!3-D
	common /size/ umask(imt,km,jmt)
	common /size/ vmask(imt,km,jmt)
	common /size/ tmask(imt,km+1,jmt), lmask(imt,km,jmt)
	common /size/ bmask(imt, km, jmt), muv_read(imt,jmt,2)
	common /size/ mld_read(imt, jmt, 2), eta(imt, jmt, -1:1)
	common /size/ dic_read(imt, km, jmt), alk_read(imt,km,jmt)
	common /size/ mask(imt,km,jmt), tau_read(imt, jmt, 2)

	common /size/ rho(imt,km,jmt)
	common /size/ DIFF_4_X_BH(imt, km, jmt)
	common /size/ DIFF_4_Y_BH(imt, km, jmt)
	common /size/ rvf_read(imt,jmt,2)
	common /size/ pme_read(imt, jmt, 2)
	common /size/ pmr_read(imt, jmt, 2)
	common /size/ eta_read(imt, jmt, 2)
	common /size/ sflux_read(imt, jmt, 2)
	common /size/ hflux_read(imt, jmt, 2)
	common /size/ u10m_read(imt, jmt, 2)
	common /size/ fice_read(imt,jmt,2), xkw_read(imt,jmt,2)
	common /size/ cfc_mix_ratio(imt,km,jmt)
	common /size/ p_read(imt,jmt,2), sst_read(imt, jmt, 2)
	common /size/ ssn_read(imt,jmt,2)
	common /size/ true_mask(imt,km,jmt)
	common /size/ uvel_read(imt,km,jmt)
	common /size/ vvel_read(imt,km,jmt)
	common /size/ trestart_read(imt, km, jmt)
	common /size/ srestart_read(imt, km, jmt)
	common /size/ error_corr(imt,km,jmt)
	common /size/ vmix_coeff(imt,km,jmt)
	common /size/ correction_kh(19,lm-1,km)
	common /size/ optimal_kh(19,lm-1,km)
	common /size/ oxyg_read(imt, km, jmt), dop_read(imt,km,jmt)
	common /size/ ralkalinity_read(imt, jmt,2)
	common /size/ salt_read(imt,jmt,2)
	common /size/ chl_read(imt, jmt,2), par_read(imt, jmt,2)
	common /size/ phosp(imt, km, jmt), phospmean(imt,km,jmt)

!2-D
	common /size/ mld(imt,jmt), tau(imt,jmt)
	common /size/ muv(imt,jmt), kmask(imt,jmt)
	common /size/ max_mask(imt,jmt)

	common /size/ sflux(imt,jmt), pmf(imt, jmt)
!	common /size/ pme(imt, jmt), pmr(imt, jmt)
	common /size/ rvf(imt,jmt), hflux(imt, jmt)
	common /size/ u10m(imt, jmt)
	common /size/ rimpl_vmix_send(km,km)
	common /size/ fice(imt,jmt), xkw(imt,jmt), p(imt,jmt)
	common /size/ sst(imt,jmt), ssn(imt,jmt)
	common /size/ rCFC_print (imt,jmt), restore(imt,jmt)
	common /size/ dxdr(imt,jmt), dydr(imt,jmt), dxdy(imt, jmt)
	common /size/ pCFC_11_read(2,2), pCFC_12_read(2,2)
	common /size/ ocanic_zone(imt, jmt)
	common /size/ co2_mask(imt, jmt)
        common /size/ pH(imt,jmt)
	common /size/ salt(imt,jmt), chl(imt, jmt)
	common /size/ ralkalinity(imt, jmt)
	common /size/ par(imt, jmt)
	common /size/ rVA (2,6), rZA(2,6)
	common /size/ p_region_mask(imt, jmt)
	common /size/ p_export_mask(imt, jmt)
	common /size/ ro2ex_2d(imt,jmt), rco2ex_2d(imt,jmt)
	common /size/ ph_2d(imt, jmt), pco2_2d(imt, jmt)


!1-D

	common /size/ fft_send_1(imt+2), fft_filtered_1(imt+2)
	common /size/ pCFC_11_North(jmt), pCFC_11_South(jmt)
	common /size/ pCFC_12_North(jmt), pCFC_12_South(jmt)
	common /size/ rK31(km), rK32(km), rK33(km)
	common /size/ rmaskave(imt)
	common /size/ jrow(jmt), phi(imt), theta(jmt)
	common /size/ Sh(km-1), Sm(km-1),rkpp_nlocal (0:km)
	common /size/ sigma(km)
        common /size/ Rf(km-1)
        common /size/ zdepth(0:km), ddz(0:km)
        common /size/ q(km-1), q_backup(km-1)
        common /size/ rl(km-1)
	common /size/ dpm(lm)
	common /size/ xmatric(jmt)
	common /size/ Am(jmt), Am_bak(jmt), Am_bh(jmt)
	common /size/ racc(km)
	common /size/ co2_anthropogenic_read(2)
	common /size/ rsum_dic(nproc), rsum_alk(nproc), rsum_pme(nproc)
	common /size/ saltglobal_ave(lm)
	common /size/ saltglobal_ave_interannual (0:1000)
	common /size/ dic_global_ave_init_km(km)
	common /size/ rI_par(km)
	common /size/ rbgc_F(km)

!!!!!!!!!!!!!! Integer declarations !!!!!!!!!!!!!!!!!!!!!!!

	integer jstart, jend, i_co2_injection_region, j_co2_injection_region
	
	common /size/ jstart(nproc), jend(nproc)
        common /size/ i_co2_injection_region(2)
	common /size/ j_co2_injection_region(2)
