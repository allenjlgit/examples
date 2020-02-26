() = evalfile("/nfs/grs1915/h6/allenjl/isis_start/.isisrc");
() = evalfile("readascii"); 
require("xspec");
xspec_abund("wilm");
xspec_xsect("vern");

%======================================================================
%Load all observations of SAX J1808
%Compute hardness ratio, energy flux, and perform rough spectral fit
%Save spectral plot, fit parameters (in txt doc and fit files)
%Create time plot of count rate, hardness ratio, flux, and spectral index 
%======================================================================

%Delete all previously used/loaded data
	delete_data (all_data);
	delete_arf (all_arfs);
	delete_rmf (all_rmfs);
	define marquardt(){ set_fit_method("marquardt"); }
	marquardt;

%Plotting choices
	%Plot time in days since April 9, 2015
	variable mjd_offset = 57121 ;
	variable xtime_label = "Days since April 9, 2015";

variable obsid, obs_cal_date, obs_mj_date, obs_mode, arf_corr, evt_fl, src_fl, bkg_fl, num_obs;
num_obs = readascii( "dates_obsid_files.txt",  &obsid, &obs_cal_date, &obs_mj_date, 
	&obs_mode, &arf_corr, &evt_fl, &src_fl, &bkg_fl;format="%s %s %i %s %f %s %s %s",comment="#");

variable ctrt_05_10 = Double_Type[num_obs,2];
variable ctrt_2_10 = Double_Type[num_obs,2];
variable ctrt_05_2 = Double_Type[num_obs,2];

variable cts_01_15 = Double_Type[num_obs];
variable cts_05_10 = Double_Type[num_obs];

variable hr = Double_Type[num_obs,2];
variable gamma = Double_Type[num_obs,3];
variable eflux = Double_Type[num_obs,3];
variable uneflux = Double_Type[num_obs,3];
variable lum = Double_Type[num_obs,3];
variable exp_arr = Double_Type[num_obs];

%Distance to SAX J1808 in kpc
	  variable dist = 3.5 ;
	  variable kpc_cm = 3.08e21;

variable spec_plots = open_plot("saxj1808_spec_7_20_2015.ps/vcps", 2, int(ceil(num_obs/2)));
variable flux_tab = fopen("saxj1808_flux_7_20_2015_full.dat", "a");
	 
variable i;
for (i = 0; i<num_obs; i++){
%for (i = 0; i <3; i++){

    %Delete all previously used/loaded data
	delete_data (all_data);
	delete_arf (all_arfs);
	delete_rmf (all_rmfs);

    () = load_data(obsid[i]+"/src_"+obs_mode[i]+".pi");
    () = load_arf(obsid[i]+"/src_"+obs_mode[i]+".arf");
    assign_arf(1,1);

    if (obs_mode[i] == "PC") {
       variable rmf_file = "swxpc0to12s6_20130101v014.rmf";
    }
    if (obs_mode[i] =="WT") {
       variable rmf_file = "swxwt0to2s6_20131212v015.rmf";
    }    

    () = load_rmf(rmf_file+";strict=0");
    assign_rmf(1,1);

    variable exp_sec = get_data_exposure(1);
    set_arf_exposure(1, exp_sec);
    exp_arr[i] = exp_sec;

    define_back(1, obsid[i]+"/bkg_"+obs_mode[i]+".pi");
    set_back_exposure(1, exp_sec);
    
    list_data ;

    notice_values( 1, 0.3, 10;unit="keV");
    %rebin_data(1, 1);
    %notice_values( 1, 0.3, 10;unit="keV");

%Total counts in 0.5-10 keV range   

       %Get counts with background
       variable data_counts = get_data_counts(1);
 
	variable gd_05_10 = where( data_counts.bin_lo >= _A(10.) and data_counts.bin_hi <= _A(0.5) );
	variable gd_05_2 = where( data_counts.bin_lo >= _A(2.) and data_counts.bin_hi <= _A(0.5) );
	variable gd_2_10 = where( data_counts.bin_lo >= _A(10.) and data_counts.bin_hi <= _A(2) );

	%Subtract background from data counts
	variable data_counts_bkg = data_counts.value - get_back(1);
	
	%Total background subtracted counts in each band
	variable cts_05_10 = sum( data_counts_bkg[gd_05_10]);
	variable cts_05_2 = sum( data_counts_bkg[gd_05_2]);
	variable cts_2_10 = sum( data_counts_bkg[gd_2_10]);
	
	%Old
    %variable cts_05_10 = region_counts(1, _A(10.), _A(0.5));
    %variable cts_05_2 = region_counts(1, _A(2.), _A(0.5));
    %variable cts_2_10 = region_counts(1, _A(10.), _A(2.));

    ctrt_05_10[i,0] = cts_05_10 / exp_sec * 100. / arf_corr[i];
    ctrt_05_10[i,1] = sqrt(cts_05_10) / exp_sec ; %* 100. / arf_corr[i];

    ctrt_05_2[i,0] = cts_05_2 / exp_sec * 100. / arf_corr[i];
    ctrt_05_2[i,1] = sqrt(cts_05_2) / exp_sec ; %* 100. / arf_corr[i];
    
    ctrt_2_10[i,0] = cts_2_10 / exp_sec * 100. / arf_corr[i];
    ctrt_2_10[i,1] = sqrt(cts_2_10) / exp_sec ; %* 100. / arf_corr[i];

    	%Old
       %ctrt_05_10[i, 0] = cts_05_10.sum / exp_sec * 100. / arf_corr[i];
       %ctrt_05_10[i, 1] = cts_05_10.sum_err / exp_sec ;
       %ctrt_05_2[i, 0] = cts_05_2.sum / exp_sec * 100. / arf_corr[i];
       %ctrt_05_2[i, 1] = cts_05_2.sum_err / exp_sec ;
       %ctrt_2_10[i, 0] = cts_2_10.sum / exp_sec * 100. / arf_corr[i] ;
       %ctrt_2_10[i, 1] = cts_2_10.sum_err / exp_sec;

       hr[i,0] = ctrt_2_10[i,0] /  ctrt_05_2[i,0];
       hr[i,1] = hr[i,0] * sqrt( (ctrt_2_10[i,1]/ctrt_2_10[i,0])^2 +  (ctrt_05_2[i,1]/ctrt_05_2[i,0])^2);

       if (hr[i,0] < 0){
       	  hr[i,0] = 0.0 ;
	  hr[i,1] = 0.0 ;
       }

}

%Choose binning based on number of counts
	if (cts_05_10< 100.){
	   variable bin_cts = 5;
	}
	if (cts_05_10 >100.){
	   variable bin_cts = 20;
	}

      %group(1; min_sn=4.5, bounds=0.3, unit="keV");
      	rebin_data(1,bin_cts);
	notice_values( 1, 0.3, 10;unit="keV");
      
      fit_fun("TBnew(1) * (bbodyrad(1) + powerlaw(1)) ");
      load_par("start_bbpl.p");
      set_par("TBnew(1).nH", 0.13, 1);
      renorm_counts;
      fit_counts;

      conf_loop(,1,0.1);
      gamma[i,0] = get_par("powerlaw(1).PhoIndex");
      variable amin, amax;
      (amin, amax) = conf("powerlaw(1).PhoIndex", 1, 0.1);
      gamma[i,1] = gamma[i,0] - amin;
      gamma[i,2] = amax - gamma[i,0];

print("Finished BB+PL fit");
      variable fit_stat;
     fit_counts(&fit_stat);
     variable red_chi = fit_stat.statistic / (fit_stat.num_bins - fit_stat.num_variable_params);
                 
print("Plotting Spec");
      xlog; ylog;
      charsize(2.5);
      title(obs_cal_date[i]+" ObsID "+obsid[i]+" "+obs_mode[i]+" Mode, X^2="+sprintf("%.2f", red_chi)+"");
      plot_data(1;xrange={0.3,10},xlabel="");
      

print("Computing unabsorbed flux");
      fit_fun("cflux(1, TBnew(1) * ( bbodyrad(1) + powerlaw(1) ) )");
      freeze("bbodyrad(1).norm", "powerlaw(1).norm");
      fit_counts ;
      eflux[i,0] = get_par("cflux(1).lg10Flux");
      (amin,amax) = conf("cflux(1).lg10Flux");
      eflux[i,1] = amin ;
      eflux[i,2] = amax ;

      fit_fun( "TBnew(1) * cflux(1, bbodyrad(1) + powerlaw(1) ) " );   
      fit_counts ;
      uneflux[i,0] = get_par("cflux(1).lg10Flux");
      (amin,amax) = conf("cflux(1).lg10Flux");
      uneflux[i,1] = amin ;
      uneflux[i,2] = amax ;   

      lum[i,*] = 4. * 3.14 * dist^2 * kpc_cm * kpc_cm * 10.^uneflux[i,*];
      
      thaw("bbodyrad(1).norm", "powerlaw(1).norm");

      %Full flux, luminosity numbers with errors
      %() = fprintf(flux_tab, "%s %s %s %.1f %6.4f %6.4f %6.4f %6.3e %6.3e %6.3e \n", 
      	  %obsid[i], obs_mode[i], obs_cal_date[i], obs_mj_date[i],
     	  %eflux[i,0], eflux[i,1], eflux[i,2],
      	  %lum[i,0], lum[i,1], lum[i,2] );

      %Count rate, hardness ratio, photon index, abs flux, unabs luminosity
      () = fprintf(flux_tab, "%s %s %s %i %.1f %6.3f %6.3f %6.3f %6.3f %6.2f %6.2f %6.4f %6.3e \n", 
      	 obsid[i], obs_mode[i], obs_cal_date[i], obs_mj_date[i], exp_arr[i],
      	 %Count rate plus hardness ratio (one-sided error)
      	 ctrt_05_10[i,0], ctrt_05_10[i,1], hr[i,0], hr[i,1],
      	 %Photon index (one-sided error)
      	 gamma[i,0], gamma[i,2],
      	 %flux and luminosity (no error)
      	 eflux[i,0], lum[i,0],);

      %Flux & luminosity with no errors
      %() = fprintf(flux_tab, "%s %s %s %i %6.4f %6.3e \n", 
      	  %obsid[i], obs_mode[i], obs_cal_date[i], obs_mj_date[i], eflux[i,0], lum[i,0]); 

print("Moving onto next ObsID");
}

fclose(flux_tab);

close_plot(spec_plots);

%Calculate energy flux and errors
eflux[*,*] = 10.^ (eflux[*,*]);
eflux[*,1] = eflux[*,0] - eflux[*,1];
eflux[*,2] = eflux[*,2] - eflux[*,0];

%Plot count rate vs time
variable 
plots = open_plot("saxj1808_ctrt_7_20_2015.ps/vcps", 1, 2);
      resize(18, 0.6);
      ylog;
      title("SAX J1808.4-3658 Outburst");
      plotxy( obs_mj_date - mjd_offset, , , ctrt_05_10[*,0], ctrt_05_10[*,1], ctrt_05_10[*,1];
      	      xlabel=xtime_label, ylabel="Swift XRT 0.5-10 keV Count Rate (c/s)", yrange={0.0005,50.}, xrange={0,NULL} );
close_plot(plots);

	variable good = where( hr[*,0] != 0.);

plots = open_plot("saxj1808_param_7_20_2015.ps/vcps", 1, 3);
      resize(18, 0.7);
      charsize(2);
      ylog;
      title("SAX J1808.4-3658 Outburst");
      plotxy( obs_mj_date[good] - mjd_offset, , , hr[good,0], hr[good,1], hr[good,1];
      	      xlabel=xtime_label, ylabel="Hardness Ratio (2-10keV / 0.5-2 keV) ", yrange={0.01,1.}, xrange={0,NULL} );
	      
	      xlin; ylin;title("");
	plotxy( obs_mj_date[good] - mjd_offset, , , gamma[good,0], gamma[good,1], gamma[good,2];
      	      xlabel=xtime_label, ylabel="Powerlaw Photon Index", yrange={0.5,5}, xrange={0,NULL} );

	      xlin; ylog;title("");
	plotxy( obs_mj_date[good] - mjd_offset, , , eflux[good,0], eflux[good,1], eflux[goode,2];
      	      xlabel=xtime_label, ylabel="Energy Flux (erg/scm^2/s)", yrange={NULL,NULL}, xrange={0,NULL} );

close_plot(plots);
