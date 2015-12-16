%Finds the best slope and intercept values for the S1a/S1b -> Field dependence mapping lines
%Based on a chi^2 minimization of the offset between ER band slices compared to simulation

%This assumes the relation between the two lines is correct, but the measurement of the first line is wrong


%% Pick random slope and intercept values over a loop: 
%s1ab_s2_xyz_fit_map.p1 from 0 to -4, s1ab_s2_xyz_fit_map.p2 from -6 to 6, 
%s1ab_s1_xyz_fit_map.p1 from 0 to 4, s1ab_s1_xyz_fit_map.p1 from -6 to 6




param_counter=1;
    best_chi2_midtobot=10^100;
    best_chi2_midtotop=10^100;
    best_chi2_combined=10^100;
    
    
    %initialize matrices
a2_list=zeros(1,2225158);
b2_list=zeros(1,2225158);
a1_list=zeros(1,2225158);
b1_list=zeros(1,2225158);
g1_list=zeros(1,2225158);
g2_list=zeros(1,2225158);
total_h3_chi2_list=zeros(1,2225158);
total_kr_chi2_list=zeros(1,2225158);
total_chi2_list=zeros(1,2225158);
alternate_total_chi2_list=zeros(1,2225158);
h3_chi2_2014_list=zeros(1,2225158);
h3_chi2_2015_list=zeros(1,2225158);
kr_chi2_2014_list=zeros(1,2225158);
kr_chi2_2015_list=zeros(1,2225158);
electron_lifetime_2014_list=zeros(1,2225158);
electron_lifetime_2015_list=zeros(1,2225158);
kr_staterr_2014_list=zeros(1,2225158);
kr_staterr_2015_list=zeros(1,2225158);
kr_sigma_2014_list=zeros(1,2225158);
kr_sigma_2015_list=zeros(1,2225158);
kr_2014_num_events=zeros(1,2225158);
kr_2015_num_events=zeros(1,2225158);
kr_means_2014_list=zeros(1,2225158);
kr_means_2015_list=zeros(1,2225158);
full_kr_chi2_2014_list=zeros(1,2225158);
full_kr_chi2_2015_list=zeros(1,2225158);
full_kr_2014_resolution_list=zeros(1,2225158);
full_kr_2015_resolution_list=zeros(1,2225158);
full_kr_2014_resolution_err_list=zeros(1,2225158);
full_kr_2015_resolution_err_list=zeros(1,2225158);

kr_chi2_2014_zslice_list=zeros(2225158,length(40+20/2:20:300-20/2));
kr_chi2_2015_zslice_list=zeros(2225158,length(40+20/2:20:300-20/2));
kr_chi2_zslice_zcenter_list=zeros(2225158,length(40+20/2:20:300-20/2));
        
total_h3_chi2_map=zeros(61,61,13,46);
total_kr_chi2_map=zeros(61,61,13,46);
total_chi2_map=zeros(61,61,13,46);
alternate_total_chi2_map=zeros(61,61,13,46);
h3_chi2_2014_map=zeros(61,61,13,46);
kr_chi2_2014_map=zeros(61,61,13,46);
h3_chi2_2015_map=zeros(61,61,13,46);
kr_chi2_2015_map=zeros(61,61,13,46);
electron_lifetime_2014_map=zeros(61,61,13,46);
electron_lifetime_2015_map=zeros(61,61,13,46);
full_kr_2014_resolution_map=zeros(61,61,13,46);
full_kr_2015_resolution_map=zeros(61,61,13,46);
full_kr_2014_resolution_err_map=zeros(61,61,13,46);
full_kr_2015_resolution_err_map=zeros(61,61,13,46);

best_total_chi2=1000000000000000000;
best_total_alternate_chi2=1000000000000000000;

inpaint_on=1;
% Measurement is a1=-0.765; b1=2.939; a2=0.1934; b2=-0.744;
    
load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p18\NEST_CH3T_Spectrum.mat')

a1_counter=1;
for a1=-3:0.1:3; %p1*x^2 coefficient for s1ab_s2_xyz_fit.p1
    a2_counter=1;
    for a2=-3:0.1:3; %p1*x^2 coefficient for s1ab_s2_xyz_fit.p1
           
             
      
g1_counter=1;
g2_counter=1;

 %%UPDATE THIS
             
s1ab_s2_xyz_fit_map.p1=a1;
b1=1-2.7322.*a1;
s1ab_s2_xyz_fit_map.p2=b1;

s1ab_s1_xyz_fit_map.p1=a2;
b2=1-2.7322.*a2;
s1ab_s1_xyz_fit_map.p2=b2;
%% Load the Uncorrected Sep 2015 Kr data


%HAVE TO MAKE THIS .MAT FILE STILL - USE class 1 and class 9 for S1
load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p18\MappingParametersFit_Start_2015_CutDown.mat')


%% Calculate 2015 Corrections

try
 %If we can use 3D map, do so.  If there aren't enough events do 1D then 2D
%If we can use 3D map, do so.  If there aren't enough events do 1D then 2D
 if length(s1ab_z)>100000; %require 100,000 events to do 3D binning
%Removing zeros from map with nearest neighbor method
s1ab_xyz_mean_old=s1ab_xyz_mean;
temp_map=s1ab_xyz_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s1ab_xyz_mean=temp_map; temp_mapQ=temp_map(q); s1ab_xyz_mean(~q)=temp_mapQ(r);s1ab_xyz_mean=reshape(s1ab_xyz_mean,size(s1ab_xyz_mean_old,1),size(s1ab_xyz_mean_old,2),size(s1ab_xyz_mean_old,3));

%interpolating with cubic so that we get NaN when extrapolating
 s1as1b_values=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1ab_xyz_mean,s2x,s2y,drift_time,'cubic');
 s1as1b_at_center=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1ab_xyz_mean,x_center,y_center,z_center,'cubic');
 s2_3Dfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_values+s1ab_s2_xyz_fit_map.p2)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_at_center+s1ab_s2_xyz_fit_map.p2); 
 s1_3Dfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_values+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_at_center+s1ab_s1_xyz_fit_map.p2); 
 
%Filling in NaN values with the s1a/s1b Z map polynomial 
  s1as1b_zvalues= polyval(s1ab_P,drift_time);
  s1as1b_z_at_center=polyval(s1ab_P,z_center);   
  s2_zfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zvalues+s1ab_s2_xyz_fit_map.p2)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_z_at_center+s1ab_s2_xyz_fit_map.p2);
  s1_zfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zvalues+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_z_at_center+s1ab_s1_xyz_fit_map.p2); 
  s2_3Dfieldremoval(isnan(s2_3Dfieldremoval))=s2_zfieldremoval(isnan(s2_3Dfieldremoval));
  s1_3Dfieldremoval(isnan(s1_3Dfieldremoval))=s1_zfieldremoval(isnan(s1_3Dfieldremoval));
 
s2_phe_bottom=s2_phe_bottom./s2_3Dfieldremoval;
s2_phe_both=s2_phe_both./s2_3Dfieldremoval;
s1_phe_bottom=s1_phe_bottom./s1_3Dfieldremoval;
s1_phe_both=s1_phe_both./s1_3Dfieldremoval;


 else 
%Removing zeros from map with nearest neighbor method
s1ab_xy_mean_old=s1ab_xy_mean;
temp_map=s1ab_xy_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s1ab_xy_mean=temp_map; temp_mapQ=temp_map(q); s1ab_xy_mean(~q)=temp_mapQ(r);s1ab_xy_mean=reshape(s1ab_xy_mean,size(s1ab_xy_mean_old,1),size(s1ab_xy_mean_old,2));

%interpolating with cubic so that we get NaN when extrapolating
 s1as1b_zvalues= polyval(s1ab_P,drift_time);
 s1as1b_zcorr_xyvalues=interp2(s1ab_xbins,s1ab_ybins,s1ab_xy_mean,s2x,s2y,'cubic');
 s1as1b_z_at_center=polyval(s1ab_P,z_center); 
 s1as1b_xy_at_center=interp2(s1ab_xbins,s1ab_ybins,s1ab_xy_mean,x_center,y_center,'cubic');
 
 s2_zfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zvalues+s1ab_s2_xyz_fit_map.p2)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_z_at_center+s1ab_s2_xyz_fit_map.p2);
 s2_xyfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zcorr_xyvalues+s1ab_s2_xyz_fit_map.p2)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_xy_at_center+s1ab_s2_xyz_fit_map.p2); 
 s1_zfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zvalues+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_z_at_center+s1ab_s1_xyz_fit_map.p2); 
 s1_xyfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zcorr_xyvalues+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_xy_at_center+s1ab_s1_xyz_fit_map.p2); 

 %Filling in NaN values with no XY correction.  (Z dependence should never be NaN)
 s2_xyfieldremoval(isnan(s2_xyfieldremoval))=1;
 s2_zfieldremoval(isnan(s2_zfieldremoval))=1;
 s1_xyfieldremoval(isnan(s1_xyfieldremoval))=1;
 s1_zfieldremoval(isnan(s1_zfieldremoval))=1;

 
s2_phe_bottom=s2_phe_bottom./(s2_zfieldremoval.*s2_xyfieldremoval);
s2_phe_both=s2_phe_both./(s2_zfieldremoval.*s2_xyfieldremoval);
s1_phe_bottom=s1_phe_bottom./(s1_zfieldremoval.*s1_xyfieldremoval);
s1_phe_both=s1_phe_both./(s1_zfieldremoval.*s1_xyfieldremoval);



 end 
    clear y_s1_both y_s1_bottom Fit_s1_both Fit_s1_bottom mean_s1_both mean_s1_both
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both x xfit S
    
    %set up energy bins for the gaussian fit
    bin_s1=8;%5
    s1_bin_min=100;
    s1_bin_max=1000;
    x=s1_bin_min:bin_s1:s1_bin_max; %150:bin_s1:500;
    x=x';
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s1_phe_both,[s1_bin_min,s1_bin_max]) & inrange(s2_width,[100,1500]) & inrange(s1_width,[30,400]);   
    
    A = sort(size(s1_phe_both(cut))); %%sort, incase 1xn or nx1 vector is defined.
    bins = ceil(A(2)./min_evts_per_bin); %get the number of points, adjust binning if stats are low.
    if bins>max_num_bins
       bins = max_num_bins; % usually we will have 65 bin. 0:5:325
    end
%     bin_size = (488/bins);%span the entire detector. ~49.0 cm. Bin centers every 4\mus
     bin_size = (floor(det_edge/bins*4)/4);%span the entire detector. 325us. typically = 5us. Round to nearest 0.25.

   
    hist_s1_both = zeros(length(x),bins);
    Fit_s1_both = cell(1,bins);
    means = zeros(bins,1);
    means_error = zeros(bins,1);
    mean_index = zeros(bins,1);
    
      
    for bin=1:bins

     binStart = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   

%      time_cut = inrange(drift_distance_\mus,[binStart,binEnd]) & cut==1;
       time_cut = inrange(drift_time,[binStart,binEnd]) & cut==1;
     
        hist_s1_both(:,bin)= hist(s1_phe_both(time_cut),x)'/bin_s1;
        temp_hist=hist_s1_both(:,bin);   
            amp_start=max(temp_hist);
            mean_start=sum(temp_hist.*x)/sum(temp_hist);
            sigma_start=std(temp_hist.*x);
        Fit_s1_both{bin}=fit(x(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges

     means(bin)=Fit_s1_both{bin}.b1;
     means_error(bin) = Fit_s1_both{bin}.c1/sqrt(2)/sqrt(sum(hist_s1_both(:,bin))*bin_s1); % == sigma/sqrt(N);
     mean_index(bin) = (binStart+binEnd)/2;     

    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
    dT_range= inrange( mean_index,[zcut_min, zcut_max] );
    dT_fit=0:1:det_edge;    
    [P_s1_both S] = polyfit(mean_index(dT_range),means(dT_range),2);
    sigma_P_s1_both = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';    

    s1_both_means=means;
    s1_both_means_sigma=means_error;
    s1_both_means_bin=mean_index;
    
    
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Calculating the S2 Z-dependence (i.e. electron lifetime)
%   Using Both PMT arrays
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
    s1_z_correction=polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);
    
    %s1_both_z_correction=interp1(s1_both_means_bin,s1_both_means,160)./interp1(s1_both_means_bin,s1_both_means,drift_time);
    s1_phe_both_z=s1_phe_both.*s1_z_correction;

    clear y_s1_both y_s1_both Fit_s1_both Fit_EL xEL yfitEL xfit
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both dT_range dT_fit S x2
    
    %set up energy bins for the gaussian fit. Using both S2 only
   bin_s2=500;
   s2_bin_min=1000;
   s2_bin_max=40000;
   x2=s2_bin_min:bin_s2:s2_bin_max;
   x2=x2';
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s2_phe_both,[s2_bin_min,s2_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]) ;
    A = sort(size(s2_phe_both(cut))); %%sort, incase 1xn or nx1 vector is defined.
    bins = ceil(A(2)./min_evts_per_bin); %get the number of points, adjust binning if stats are low.
    if bins>max_num_bins
       bins = max_num_bins; % usually we will have 65 bin. 0:5:325
    end
%     bin_size = (488/bins);%span the entire detector. 325us. typically = 5us
    bin_size = (floor(det_edge/bins*4)/4);%span the entire detector. 325us. typically = 5us. Round to nearest 0.25.
   
        hist_s2_both = zeros(length(x2),bins);
        Fit_s2_both = cell(1,bins);
        means=zeros(bins,1);
        means_error = zeros(bins,1);
        mean_index= zeros(bins,1);
        
    for bin=1:bins

     binStart = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   

%      time_cut = inrange(drift_distance_\mus,[binStart,binEnd]) & cut==1;
       time_cut = inrange(drift_time,[binStart,binEnd]) & cut==1; 
       
        hist_s2_both(:,bin)= hist(s2_phe_both(time_cut),x2)'/bin_s2;
        temp_hist=hist_s2_both(:,bin);
            amp_start=max(temp_hist);
            mean_start=sum(temp_hist.*x2)/sum(temp_hist);
            sigma_start=std(temp_hist.*x2);
        Fit_s2_both{bin}=fit(x2(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges
          

         means(bin)=Fit_s2_both{bin}.b1;
         means_error(bin) = Fit_s2_both{bin}.c1/sqrt(2)/sqrt(sum(hist_s2_both(:,bin))*bin_s2); % == sigma/sqrt(N);
         mean_index(bin) = (binStart+binEnd)/2;     
     
    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
    dT_range= inrange( mean_index,[5+bin_size, 500] );%from 50 to 300us for the fit, based on when s1as1b values get strange

    s2_both_norm_z_index=mean_index(dT_range);
    s2_both_norm_z_means=means(dT_range);
    s2_at_top=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4); %normalize to right below the gate
    s2_both_norm_z=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4)./s2_both_norm_z_means;

    s2_both_means=means;
    s2_both_means_sigma=means_error;
    s2_both_means_bin=mean_index;      

    %%Tracker for pseudo-lifetime
    s2_at_bottom=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,320);
    pseudo_lifetime_2015=-320/log(s2_at_bottom/s2_at_top);
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S2 in XY (after correcting for Z-dep). Both PMT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s2_phe_both_z=s2_phe_both.*RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);

    clear x2 hist_bin

    %set up energy bins for the gaussian fit. Using both S2 only
   bin_s2=500;
   s2_bin_min=1000;
   s2_bin_max=40000;
   x2=s2_bin_min:bin_s2:s2_bin_max;
   x2=x2';

%set up matricies to be filled
mean_S2_both_xy = zeros(floor(s2_xbins),floor(s2_ybins));
Sigma_S2_both = zeros(floor(s2_xbins),floor(s2_ybins));
Count_S2_both = zeros(floor(s2_xbins),floor(s2_ybins));

time_energy_cut = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[100,s1area_bound_max]) & inrange(s2_phe_both_z,[s2_bin_min,s2_bin_max]) & inrange(s2radius,[rcut_min,rcut_max]) ...
    & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);

%%%%variables for uber turbo boost!
s2_phe_both_z_2=s2_phe_both_z(time_energy_cut);
s2x_2=s2x(time_energy_cut);
s2y_2=s2y(time_energy_cut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for x_bin=S2_xbin_min:S2xybinsize:(S2_xbin_max-S2xybinsize);
      %fprintf('X = %d \n',x_bin)
    for y_bin=S2_ybin_min:S2xybinsize:(S2_ybin_max-S2xybinsize);          
      x_min = x_bin; x_max = x_bin+S2xybinsize;
      y_min = y_bin; y_max = y_bin+S2xybinsize;      
      
        x_count=int32(1+x_bin/S2xybinsize+S2_xbin_max/S2xybinsize);
        y_count=int32(1+y_bin/S2xybinsize+S2_ybin_max/S2xybinsize); 
      
      
        bin_cut = inrange(s2x_2,[x_min,x_max]) & inrange(s2y_2,[y_min,y_max]);            
      
       hist_bin = hist(s2_phe_both_z_2(bin_cut) ,x2)'/bin_s2; 
       Count_S2_both(y_count,x_count)=sum(hist_bin)*bin_s2; % X and Y flip ... because MATLAB
       
       %%%%%%%%Truncate variables after binning is complete
        s2_phe_both_z_2=s2_phe_both_z_2(~bin_cut);
        s2x_2=s2x_2(~bin_cut);
        s2y_2=s2y_2(~bin_cut);
        clear bin_cut;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
      if Count_S2_both(y_count,x_count)>350;                         
                
                   amp_start=max(hist_bin);
                   mean_start=sum(hist_bin.*x2)/sum(hist_bin);
                   sigma_start=std(hist_bin.*x2);    
                Fit_S2_both_xy= fit(x2(2:end-1),hist_bin(2:end-1),'gauss1'); %remove edges
                
                %Peak location (mean)
                mean_S2_both_xy(y_count,x_count)=Fit_S2_both_xy.b1; % X and Y flip ... because MATLAB
                
                %1-sigma of peak position
                Sig_b=Fit_S2_both_xy.c1/sqrt(2)/sqrt(Count_S2_both(y_count,x_count)); %sigma/sqrt(N)
   

                    %Fit error checking
                        if(isnan(Sig_b))
                         mean_S2_both_xy(y_count,x_count)=0;
                         Sig_b=0;
                        end

                    %uncertainty in mean
                    Sigma_S2_both(y_count,x_count)=Sig_b;% X and Y flip ... because MATLAB
       
      else
          
        mean_S2_both_xy(y_count,x_count)=0;  %don't so gaussin fit if less than 30 events
       
      end             
    end     
   end
   
   mean_S2_both_xy(isnan(mean_S2_both_xy)) = 0;
       
   
%%Plot the mean S2_XY both %%%%%%%%%%%%%%%%

s2xbins=(S2_xbin_min+S2xybinsize/2):S2xybinsize:(S2_xbin_max-S2xybinsize/2);
s2ybins=(S2_ybin_min+S2xybinsize/2):S2xybinsize:(S2_ybin_max-S2xybinsize/2);


%Calculate corrections matrix and 1 sigma corrections matrix
% if inpaint_on==1;  mean_S2_both_xy(mean_S2_both_xy==0)=nan; mean_S2_both_xy=inpaint_nans(mean_S2_both_xy,3); end;
center=interp2(s2xbins,s2ybins,mean_S2_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_S2_all=center./mean_S2_both_xy; 
norm_S2_all(isinf(norm_S2_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_all(isnan(norm_S2_all))=1;

% if inpaint_on==1;  Sigma_S2_both(Sigma_S2_both==0)=nan; Sigma_S2_both=inpaint_nans(Sigma_S2_both,3); end;
sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_both,x_center,y_center,'spline');%Normalize to the center (x=y=0)
sigma_norm_S2_all=sqrt((sigma_center./mean_S2_both_xy).^2+(Sigma_S2_both.*center./mean_S2_both_xy.^2).^2);
sigma_norm_S2_all(isinf(sigma_norm_S2_all))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S2_all(isnan(sigma_norm_S2_all))=1;%no nan. Set Sigma=1, which is 100% uncertainty


s2_both_center_xyz=center;
s2_both_center_xyz_error=sigma_center;


        
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S1 in XY (after correction for S1 Z-dependence). Both PMT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear hist_bin x

   
    %set up energy bins for the gaussian fit
    bin_s1=8;%5
    s1_bin_min=100;
    s1_bin_max=500;
    x=s1_bin_min:bin_s1:s1_bin_max; %150:bin_s1:500;
    x=x';

%set up matricies to be filled
mean_S1_both_xy = zeros(floor(s1_xbins),floor(s1_ybins));
Sigma_S1_both = zeros(floor(s1_xbins),floor(s1_ybins));
Count_S1_both = zeros(floor(s1_xbins),floor(s1_ybins));

% s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
s1_z_correction=polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);%normalize s1 to 160 ns.

s1_phe_both_z=s1_phe_both.*s1_z_correction;
% s2_phe_both_z=s2_phe_both.*exp(drift_distance_\mus/e_lifetime); %Correct S2 for lifeimte  %s2_phe_both_z[5000,18000]
% s2_phe_both_z=s2_phe_both.*exp(drift_time/e_lifetime_both); %Correct S2 for lifeimte  %s2_phe_both_z[5000,18000]



time_energy_cut = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[s1_bin_min,s1_bin_max]) & inrange(s2_phe_both_z,[1000,s2area_bound_max]) & inrange(s2radius,[rcut_min,rcut_max]) ...
    &inrange(s2_width,[100 1000]) & inrange(s1_width,[30 400]);


%%%%variables for uber turbo boost!
s1_phe_both_z_2=s1_phe_both_z(time_energy_cut);
s2x_2=s2x(time_energy_cut);
s2y_2=s2y(time_energy_cut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for x_bin=S2_xbin_min:S2xybinsize:(S2_xbin_max-S2xybinsize);
      %fprintf('X = %d \n',x_bin)
    for y_bin=S2_ybin_min:S2xybinsize:(S2_ybin_max-S2xybinsize);          
      x_min = x_bin; x_max = x_bin+S2xybinsize;
      y_min = y_bin; y_max = y_bin+S2xybinsize;      
      
        x_count=int32(1+x_bin/S2xybinsize+S2_xbin_max/S2xybinsize);
        y_count=int32(1+y_bin/S2xybinsize+S2_ybin_max/S2xybinsize); 
      
      
      bin_cut = inrange(s2x_2,[x_min,x_max]) & inrange(s2y_2,[y_min,y_max]);            
      
       hist_bin = hist(s1_phe_both_z_2(bin_cut) ,x)'/bin_s1; 
       Count_S1_both(y_count,x_count)=sum(hist_bin)*bin_s1; % X and Y flip ... because MATLAB
       
       %%%%%%%%Truncate variables after binning is complete
        s1_phe_both_z_2=s1_phe_both_z_2(~bin_cut);
        s2x_2=s2x_2(~bin_cut);
        s2y_2=s2y_2(~bin_cut);
        clear bin_cut;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
      if Count_S1_both(y_count,x_count)>350; 
                        
                    amp_start=max(hist_bin);
                    mean_start=sum(hist_bin.*x)/sum(hist_bin);
                    sigma_start=std(hist_bin.*x);                     
                Fit_S1_both_xy= fit(x(2:end-1),hist_bin(2:end-1),'gauss1'); %remove edges
                
                %Peak location (mean)
                mean_S1_both_xy(y_count,x_count)=Fit_S1_both_xy.b1; % X and Y flip ... because MATLAB
                
                %1-sigma of peak position
                Sig_b=Fit_S1_both_xy.c1/sqrt(2)/sqrt(Count_S1_both(y_count,x_count)); %sigma/sqrt(N)
   
                    %Fit error checking
                        if(isnan(Sig_b))
                          mean_S1_both_xy(y_count,x_count)=0;
                         Sig_b=0;
                        end
                                
                    %uncertainty in mean
                    Sigma_S1_both(y_count,x_count)=Sig_b;% X and Y flip ... because MATLAB   
     
      else
          
        mean_S1_both_xy(y_count,x_count)=0;  
       
      end            
    end    
   end
   
   mean_S1_both_xy(isnan(mean_S1_both_xy)) = 0;
   
   
%%Plot the correction S1_XY both %%%%%%%%%%%%

s1xbins=(S1_xbin_min+S1xybinsize/2):S1xybinsize:(S1_xbin_max-S1xybinsize/2);
s1ybins=(S1_ybin_min+S1xybinsize/2):S1xybinsize:(S1_ybin_max-S1xybinsize/2);

%Calculate corrections matrix and 1 sigma corrections matrix
% if inpaint_on==1;  mean_S1_both_xy(mean_S1_both_xy==0)=nan; mean_S1_both_xy=inpaint_nans(mean_S1_both_xy,3); end;
center=interp2(s1xbins,s1ybins,mean_S1_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_S1_all=center./mean_S1_both_xy; 
norm_S1_all(isinf(norm_S1_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_all(isnan(norm_S1_all))=1;

% if inpaint_on==1;  Sigma_S1_both(Sigma_S1_both==0)=nan; Sigma_S1_both=inpaint_nans(Sigma_S1_both,3); end;
sigma_center=interp2(s1xbins,s1ybins,Sigma_S1_both,x_center,y_center,'spline');%Normalize to the center (x=y=0)
sigma_norm_S1_all=sqrt((sigma_center./mean_S1_both_xy).^2+(Sigma_S1_both.*center./mean_S1_both_xy.^2).^2);
sigma_norm_S1_all(isinf(sigma_norm_S1_all))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S1_all(isnan(sigma_norm_S1_all))=1;%no nan. Set Sigma=1, which is 100% uncertainty

s1_both_center_xyz=center;
s1_both_center_xyz_error=sigma_center;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S1 in XYZ. 3D Map, caluculate and Submit IQ
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% if Kr_events>150000
%    
% %     LUXkrypCal_S1_3D(s1_phe_both,s1_phe_bottom,s2_phe_both,s2x,s2y,drift_distance_\mus,user_name,file_id_cp,dp_version,submit);
%     LUXkrypCal_S1_3D_Run04_local_2p14FixedZHisto(s1_phe_both,s1_phe_bottom,s2_phe_both,s2x,s2y,drift_time,x_center,y_center,z_center,det_edge,inpaint_on,user_name,file_id_cp,dp_version, algorithm_name, submit);
% 
% 
% end    
    

%% Apply corrections to CH3T and calculate the ER band shift

%First load the data

clear s1_phe_both s2_phe_both
load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p18\CH3T_Sep2015_NoCorrections.mat');


% Apply the corrections
%-----------------------------------
% Finding the electron lifetime
%-----------------------------------

%  electron_lifetime = e_lifetime_both;

%-----------------------------------
% Finding the S1 Z-dep values
%-----------------------------------

z_dep_both_values = P_s1_both;
z_dep_par_all = z_dep_both_values;      

%------------------------------------------
% Finding the S2 xy correction map values
%------------------------------------------

   s2_x_bins = s2xbins;  
   s2_y_bins = s2ybins; 
   s2_map_all = norm_S2_all; 
   
%------------------------------------------
% Finding the S1 xy correction map values
%------------------------------------------

      
   s1_x_bins = s1xbins; 
   s1_y_bins = s1ybins; 
   s1_map_all = norm_S1_all; 

   
% %------------------------------------------
% % Finding the S1 xyz correction map values
% %------------------------------------------
   
% 
% if nargin == 3      
%    s1_xyz_x_bins = s1xbins_3D;
%    s1_xyz_y_bins = s1ybins_3D;
%    s1_xyz_z_bins = s1zbins_3D;
%    
% %    xx = size(s1_xyz_x_bins); yy = size(s1_xyz_y_bins); zz = size(s1_xyz_z_bins);
%    
%    s1_xyz_map_all = norm_s1_both_xyz; %lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_both_xyz; s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
%    s1_xyz_map_bottom = norm_s1_bottom_xyz; %lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_bottom_xyz;   s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));
% end



%% Calculating corrections

%% Finding the drift time, x and y associated with S1 - this isn't how they are actually applied in DP, since DP estimates the S1 and here we know the exact S1.  
% It's the best we can do for ease here and it shouldn't matter much


s1_drift_ns = drift_time.*1000;
s1_ref_z_ns=z_center*1000;
s1_x_cm = S2_x;
s1_y_cm = S2_y;

%--------------------------------------------------------------------------
% Calculate Z corrections
%--------------------------------------------------------------------------

% Calculate S1 Z-correction
s1_z_correction = polyval(z_dep_par_all,s1_ref_z_ns./1000)./polyval(z_dep_par_all,s1_drift_ns./1000);
s1_z_correction(find(isnan(s1_z_correction))) = 1.0;

s2_z_correction = RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);

% % Calculate electron lifetime correction (S2 Z-correction)
% electron_lifetime_correction = exp(drift_time./electron_lifetime);
% electron_lifetime_correction(find(isnan(electron_lifetime_correction))) = 1.0;
 
%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------

% Calculate S1 XY-only corrections
s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'spline',1);
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

% Calculate S2 XY corrections
s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,S2_x,S2_y,'spline',1);
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

s1xyz_correction=s1xy_correction.*s1_z_correction;

%% add RQs of correction factors

% d.correction_electron_lifetime = single(electron_lifetime_correction);
d.correction_s1_z_dependence = single(s1_z_correction);
d.correction_s1_xyz_dependence = single(s1xyz_correction); %this is either 3D or 1D+2D
d.correction_s1_xy_dependence = single(s1xy_correction);
d.correction_s2_xy_dependence = single(s2xy_correction);

% d.admin.corrections.electron_lifetime = electron_lifetime;
d.admin.corrections.s1_z_dependence.all = z_dep_par_all;

%--------------------------------------------------------------------------
% Apply Z corrections
%--------------------------------------------------------------------------

s1_phe_both_z=s1_phe_both.*s1_z_correction;
s2_phe_both_z=s2_phe_both.*s2_z_correction;


%% Calculate Energy spectra
s1_phe_both_xyz=s1_phe_both.*s1xyz_correction;
s2_phe_both_xyz=s2_phe_both.*s2_z_correction.*s2xy_correction;

   h3_cut= inrange(drift_time,[40 300]) & inrange(s2radius_c,[0,18]) & inrange(s1_phe_both,[0 120]) ...
    & (log10(s1_phe_both)+0.5*log10(s2_phe_both_xyz)<4) & s2_phe_both>100 & s2_phe_both<16000 ...
    & drift_time<-2.*s2radius.^2+380; %used to be drift time 100 to 250
    
h3_s1_c_2015=s1_phe_both_xyz(h3_cut);
h3_s2_c_2015=s2_phe_both_xyz(h3_cut);


   kr_cut= inrange(drift_time,[40 300]) & inrange(s2radius_c,[0,18]) & inrange(s1_phe_both,[30 1000]) ...
    & (log10(s1_phe_both)+0.5*log10(s2_phe_both_xyz)>4) & s2_phe_both>3000 & s2_phe_both<60000 ...
    & drift_time<-2.*s2radius.^2+380; %used to be drift time 100 to 250

kr_s1_c_2015=s1_phe_both_xyz(kr_cut);
kr_s2_c_2015=s2_phe_both_xyz(kr_cut);
kr_2015_drift_time=drift_time(kr_cut);
kr_2015_s2x=S2_x(kr_cut);
kr_2015_s2y=S2_y(kr_cut);

 clearvars -except param_counter a1_list b1_list a2_list b2_list total_chi2_list h3_chi2_list kr_chi2_list g1_list g2_list total_chi2_map h3_chi2_map kr_chi2_map chi2_midtobot_list...
    best_chi2_midtotop besta1_midtotop besta2_midtotop bestb1_midtotop bestb2_midtotop best_chi2_midtobot...
    besta1_midtobot besta1_midtotop bestb1_midtobot bestb1_midtotop a1 b1 a2 b2 chi2_midtotop_map chi2_midtobot_map b_counter a_counter...
    besta1_combined bestb1_combined besta2_combined bestb2_combined best_chi2_combined chi2_combined_map chi2_combined_list b1_counter a1_counter a2_counter b2_counter...
    a1_counter b1_counter a2_counter b2_counter c1_list c2_list bestc1_midtobot bestc2_midtobot bestc1_midtotop bestc2_midtotop bestc1_combined bestc2_combined...
    besta1_combined besta2_combined bestb1_combined bestb2_combined besta1 bestb1 bestc1 besta2 bestb2 bestc2 bestg1 bestg2 best_total_chi2 NEST_CH3T_Spectrum inpaint_on...
    best_kr_chi2 best_h3_chi2 electron_lifetime_list electron_lifetime_map a2_start a1_start b1_start b2_start ...
    total_h3_chi2_list total_kr_chi2_list total_chi2_list h3_chi2_2014_list h3_chi2_2015_list kr_chi2_2014_list kr_chi2_2015_list ...
    electron_lifetime_2014_list electron_lifetime_2015_list total_h3_chi2_map total_kr_chi2_map total_chi2_map h3_chi2_2014_map h3_chi2_2015_map ...
    kr_chi2_2014_map kr_chi2_2015_map electron_lifetime_2014_map electron_lifetime_2015_map NEST_CH3T_Spectrum g1 g2 bestg1 bestg2 g1_counter g2_counter ...
    firsta1 firstb1 firsta2 firstb2 firstg1 firstg2 kr_means_2015_map kr_means_2014_map kr_sigma_2015_map kr_sigma_2014_map kr_staterr_2015_map kr_staterr_2014_map ...
    kr_means_2015_list kr_means_2014_list kr_2015_num_events kr_2014_num_events kr_sigma_2015_list kr_sigma_2014_list kr_staterr_2015_list kr_staterr_2014_list ...
    pseudo_lifetime_2014 pseudo_lifetime_2015 kr_s1_c_2015 kr_s2_c_2015 h3_s1_c_2015 h3_s2_c_2015 s1ab_s2_xyz_fit_map s1ab_s1_xyz_fit_map ...
    kr_2015_drift_time kr_2015_s2x kr_2015_s2y kr_chi2_2014_zslice_list kr_chi2_2015_zslice_list kr_chi2_zslice_zcenter_list ...
    full_kr_chi2_2014_list full_kr_chi2_2015_list best_total_kr_chi2 best_total_h3_chi2 alternate_total_chi2_map alternate_total_chi2_list ...
    best_total_alternate_chi2 best_total_alternate_kr_chi2 best_total_alternate_h3_chi2 besta2_alternate bestb2_alternate bestg1_alternate bestg2_alternate bestc2_alternate c2 ...
    full_kr_2014_resolution_list full_kr_2015_resolution_list full_kr_2014_resolution_err_list full_kr_2015_resolution_err_list ...
    full_kr_2014_resolution_map full_kr_2015_resolution_map full_kr_2014_resolution_err_map full_kr_2015_resolution_err_map besta1_alternate besta2_alternate



             


%%%%%%%%%%%%%%%%%%%%%
%% REPEAT FOR 2014 DATA
%%%%%%%%%%%%%%%%%%%%%

clear s1_phe_both s2_phe_both
load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p18\MappingParametersFit_Start_2014_CutDown.mat')

 %If we can use 3D map, do so.  If there aren't enough events do 1D then 2D
%If we can use 3D map, do so.  If there aren't enough events do 1D then 2D
 if length(s1ab_z)>100000; %require 100,000 events to do 3D binning
%Removing zeros from map with nearest neighbor method
s1ab_xyz_mean_old=s1ab_xyz_mean;
temp_map=s1ab_xyz_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s1ab_xyz_mean=temp_map; temp_mapQ=temp_map(q); s1ab_xyz_mean(~q)=temp_mapQ(r);s1ab_xyz_mean=reshape(s1ab_xyz_mean,size(s1ab_xyz_mean_old,1),size(s1ab_xyz_mean_old,2),size(s1ab_xyz_mean_old,3));

%interpolating with cubic so that we get NaN when extrapolating
 s1as1b_values=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1ab_xyz_mean,s2x,s2y,drift_time,'cubic');
 s1as1b_at_center=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1ab_xyz_mean,x_center,y_center,z_center,'cubic');
 s2_3Dfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_values+s1ab_s2_xyz_fit_map.p2)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_at_center+s1ab_s2_xyz_fit_map.p2); 
 s1_3Dfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_values+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_at_center+s1ab_s1_xyz_fit_map.p2); 
 
%Filling in NaN values with the s1a/s1b Z map polynomial 
  s1as1b_zvalues= polyval(s1ab_P,drift_time);
  s1as1b_z_at_center=polyval(s1ab_P,z_center);   
  s2_zfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zvalues+s1ab_s2_xyz_fit_map.p2)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_z_at_center+s1ab_s2_xyz_fit_map.p2);
  s1_zfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zvalues+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_z_at_center+s1ab_s1_xyz_fit_map.p2); 
  s2_3Dfieldremoval(isnan(s2_3Dfieldremoval))=s2_zfieldremoval(isnan(s2_3Dfieldremoval));
  s1_3Dfieldremoval(isnan(s1_3Dfieldremoval))=s1_zfieldremoval(isnan(s1_3Dfieldremoval));
 
s2_phe_bottom=s2_phe_bottom./s2_3Dfieldremoval;
s2_phe_both=s2_phe_both./s2_3Dfieldremoval;
s1_phe_bottom=s1_phe_bottom./s1_3Dfieldremoval;
s1_phe_both=s1_phe_both./s1_3Dfieldremoval;


 else 
%Removing zeros from map with nearest neighbor method
s1ab_xy_mean_old=s1ab_xy_mean;
temp_map=s1ab_xy_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s1ab_xy_mean=temp_map; temp_mapQ=temp_map(q); s1ab_xy_mean(~q)=temp_mapQ(r);s1ab_xy_mean=reshape(s1ab_xy_mean,size(s1ab_xy_mean_old,1),size(s1ab_xy_mean_old,2));

%interpolating with cubic so that we get NaN when extrapolating
 s1as1b_zvalues= polyval(s1ab_P,drift_time);
 s1as1b_zcorr_xyvalues=interp2(s1ab_xbins,s1ab_ybins,s1ab_xy_mean,s2x,s2y,'cubic');
 s1as1b_z_at_center=polyval(s1ab_P,z_center); 
 s1as1b_xy_at_center=interp2(s1ab_xbins,s1ab_ybins,s1ab_xy_mean,x_center,y_center,'cubic');
 
 s2_zfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zvalues+s1ab_s2_xyz_fit_map.p2)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_z_at_center+s1ab_s2_xyz_fit_map.p2);
 s2_xyfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zcorr_xyvalues+s1ab_s2_xyz_fit_map.p2)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_xy_at_center+s1ab_s2_xyz_fit_map.p2); 
 s1_zfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zvalues+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_z_at_center+s1ab_s1_xyz_fit_map.p2); 
 s1_xyfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zcorr_xyvalues+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_xy_at_center+s1ab_s1_xyz_fit_map.p2); 

 %Filling in NaN values with no XY correction.  (Z dependence should never be NaN)
 s2_xyfieldremoval(isnan(s2_xyfieldremoval))=1;
 s2_zfieldremoval(isnan(s2_zfieldremoval))=1;
 s1_xyfieldremoval(isnan(s1_xyfieldremoval))=1;
 s1_zfieldremoval(isnan(s1_zfieldremoval))=1;

 
s2_phe_bottom=s2_phe_bottom./(s2_zfieldremoval.*s2_xyfieldremoval);
s2_phe_both=s2_phe_both./(s2_zfieldremoval.*s2_xyfieldremoval);
s1_phe_bottom=s1_phe_bottom./(s1_zfieldremoval.*s1_xyfieldremoval);
s1_phe_both=s1_phe_both./(s1_zfieldremoval.*s1_xyfieldremoval);



 end 
    
    clear y_s1_both y_s1_bottom Fit_s1_both Fit_s1_bottom mean_s1_both mean_s1_both
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both x xfit S
    
    %set up energy bins for the gaussian fit
    bin_s1=8;%5
    s1_bin_min=100;
    s1_bin_max=500;
    x=s1_bin_min:bin_s1:s1_bin_max; %150:bin_s1:500;
    x=x';
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s1_phe_both,[s1_bin_min,s1_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);   
    
    A = sort(size(s1_phe_both(cut))); %%sort, incase 1xn or nx1 vector is defined.
    bins = ceil(A(2)./min_evts_per_bin); %get the number of points, adjust binning if stats are low.
    if bins>max_num_bins
       bins = max_num_bins; % usually we will have 65 bin. 0:5:325
    end
%     bin_size = (488/bins);%span the entire detector. ~49.0 cm. Bin centers every 4\mus
     bin_size = (floor(det_edge/bins*4)/4);%span the entire detector. 325us. typically = 5us. Round to nearest 0.25.

   
    hist_s1_both = zeros(length(x),bins);
    Fit_s1_both = cell(1,bins);
    means = zeros(bins,1);
    means_error = zeros(bins,1);
    mean_index = zeros(bins,1);
    
      
    for bin=1:bins

     binStart = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   

%      time_cut = inrange(drift_distance_\mus,[binStart,binEnd]) & cut==1;
       time_cut = inrange(drift_time,[binStart,binEnd]) & cut==1;
     
        hist_s1_both(:,bin)= hist(s1_phe_both(time_cut),x)'/bin_s1;
        temp_hist=hist_s1_both(:,bin);   
            amp_start=max(temp_hist);
            mean_start=sum(temp_hist.*x)/sum(temp_hist);
            sigma_start=std(temp_hist.*x);
        Fit_s1_both{bin}=fit(x(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges

     means(bin)=Fit_s1_both{bin}.b1;
     means_error(bin) = Fit_s1_both{bin}.c1/sqrt(2)/sqrt(sum(hist_s1_both(:,bin))*bin_s1); % == sigma/sqrt(N);
     mean_index(bin) = (binStart+binEnd)/2;     

    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
    dT_range= inrange( mean_index,[zcut_min, zcut_max] );
    dT_fit=0:1:det_edge;    
    [P_s1_both S] = polyfit(mean_index(dT_range),means(dT_range),2);
    sigma_P_s1_both = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';    

    s1_both_means=means;
    s1_both_means_sigma=means_error;
    s1_both_means_bin=mean_index;
    
    
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Calculating the S2 Z-dependence (i.e. electron lifetime)
%   Using Both PMT arrays
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
    s1_z_correction=polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);
    
    %s1_both_z_correction=interp1(s1_both_means_bin,s1_both_means,160)./interp1(s1_both_means_bin,s1_both_means,drift_time);
    s1_phe_both_z=s1_phe_both.*s1_z_correction;

    clear y_s1_both y_s1_both Fit_s1_both Fit_EL xEL yfitEL xfit
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both dT_range dT_fit S x2
    
    %set up energy bins for the gaussian fit. Using both S2 only
   bin_s2=500;
   s2_bin_min=1000;
   s2_bin_max=40000;
   x2=s2_bin_min:bin_s2:s2_bin_max;
   x2=x2';
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s2_phe_both,[s2_bin_min,s2_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]) ;
    A = sort(size(s2_phe_both(cut))); %%sort, incase 1xn or nx1 vector is defined.
    bins = ceil(A(2)./min_evts_per_bin); %get the number of points, adjust binning if stats are low.
    if bins>max_num_bins
       bins = max_num_bins; % usually we will have 65 bin. 0:5:325
    end
%     bin_size = (488/bins);%span the entire detector. 325us. typically = 5us
    bin_size = (floor(det_edge/bins*4)/4);%span the entire detector. 325us. typically = 5us. Round to nearest 0.25.
   
        hist_s2_both = zeros(length(x2),bins);
        Fit_s2_both = cell(1,bins);
        means=zeros(bins,1);
        means_error = zeros(bins,1);
        mean_index= zeros(bins,1);
        
    for bin=1:bins

     binStart = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   

%      time_cut = inrange(drift_distance_\mus,[binStart,binEnd]) & cut==1;
       time_cut = inrange(drift_time,[binStart,binEnd]) & cut==1; 
       
        hist_s2_both(:,bin)= hist(s2_phe_both(time_cut),x2)'/bin_s2;
        temp_hist=hist_s2_both(:,bin);
            amp_start=max(temp_hist);
            mean_start=sum(temp_hist.*x2)/sum(temp_hist);
            sigma_start=std(temp_hist.*x2);
        Fit_s2_both{bin}=fit(x2(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges
          

         means(bin)=Fit_s2_both{bin}.b1;
         means_error(bin) = Fit_s2_both{bin}.c1/sqrt(2)/sqrt(sum(hist_s2_both(:,bin))*bin_s2); % == sigma/sqrt(N);
         mean_index(bin) = (binStart+binEnd)/2;     
     
    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
    dT_range= inrange( mean_index,[5+bin_size, 500] );%from 50 to 300us for the fit, based on when s1as1b values get strange

s2_both_norm_z_index=mean_index(dT_range);
s2_both_norm_z_means=means(dT_range);
    s2_at_top=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4); %normalize to right below the gate
    s2_both_norm_z=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,4)./s2_both_norm_z_means;

    s2_both_means=means;
    s2_both_means_sigma=means_error;
    s2_both_means_bin=mean_index;      
    
    %%Tracker for pseudo-lifetime
    s2_at_bottom=RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z_means,320);
    pseudo_lifetime_2014=-320/log(s2_at_bottom/s2_at_top);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S2 in XY (after correcting for Z-dep). Both PMT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s2_phe_both_z=s2_phe_both.*RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);

    clear x2 hist_bin

    %set up energy bins for the gaussian fit. Using both S2 only
   bin_s2=500;
   s2_bin_min=1000;
   s2_bin_max=40000;
   x2=s2_bin_min:bin_s2:s2_bin_max;
   x2=x2';

%set up matricies to be filled
mean_S2_both_xy = zeros(floor(s2_xbins),floor(s2_ybins));
Sigma_S2_both = zeros(floor(s2_xbins),floor(s2_ybins));
Count_S2_both = zeros(floor(s2_xbins),floor(s2_ybins));

time_energy_cut = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[100,s1area_bound_max]) & inrange(s2_phe_both_z,[s2_bin_min,s2_bin_max]) & inrange(s2radius,[rcut_min,rcut_max]) ...
    & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);

%%%%variables for uber turbo boost!
s2_phe_both_z_2=s2_phe_both_z(time_energy_cut);
s2x_2=s2x(time_energy_cut);
s2y_2=s2y(time_energy_cut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for x_bin=S2_xbin_min:S2xybinsize:(S2_xbin_max-S2xybinsize);
      %fprintf('X = %d \n',x_bin)
    for y_bin=S2_ybin_min:S2xybinsize:(S2_ybin_max-S2xybinsize);          
      x_min = x_bin; x_max = x_bin+S2xybinsize;
      y_min = y_bin; y_max = y_bin+S2xybinsize;      
      
        x_count=int32(1+x_bin/S2xybinsize+S2_xbin_max/S2xybinsize);
        y_count=int32(1+y_bin/S2xybinsize+S2_ybin_max/S2xybinsize); 
      
      
        bin_cut = inrange(s2x_2,[x_min,x_max]) & inrange(s2y_2,[y_min,y_max]);            
      
       hist_bin = hist(s2_phe_both_z_2(bin_cut) ,x2)'/bin_s2; 
       Count_S2_both(y_count,x_count)=sum(hist_bin)*bin_s2; % X and Y flip ... because MATLAB
       
       %%%%%%%%Truncate variables after binning is complete
        s2_phe_both_z_2=s2_phe_both_z_2(~bin_cut);
        s2x_2=s2x_2(~bin_cut);
        s2y_2=s2y_2(~bin_cut);
        clear bin_cut;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
      if Count_S2_both(y_count,x_count)>350;                         
                
                   amp_start=max(hist_bin);
                   mean_start=sum(hist_bin.*x2)/sum(hist_bin);
                   sigma_start=std(hist_bin.*x2);    
                Fit_S2_both_xy= fit(x2(2:end-1),hist_bin(2:end-1),'gauss1'); %remove edges
                
                %Peak location (mean)
                mean_S2_both_xy(y_count,x_count)=Fit_S2_both_xy.b1; % X and Y flip ... because MATLAB
                
                %1-sigma of peak position
                Sig_b=Fit_S2_both_xy.c1/sqrt(2)/sqrt(Count_S2_both(y_count,x_count)); %sigma/sqrt(N)
   

                    %Fit error checking
                        if(isnan(Sig_b))
                         mean_S2_both_xy(y_count,x_count)=0;
                         Sig_b=0;
                        end

                    %uncertainty in mean
                    Sigma_S2_both(y_count,x_count)=Sig_b;% X and Y flip ... because MATLAB
       
      else
          
        mean_S2_both_xy(y_count,x_count)=0;  %don't so gaussin fit if less than 30 events
       
      end             
    end     
   end
   
   mean_S2_both_xy(isnan(mean_S2_both_xy)) = 0;
       
   
%%Plot the mean S2_XY both %%%%%%%%%%%%%%%%

s2xbins=(S2_xbin_min+S2xybinsize/2):S2xybinsize:(S2_xbin_max-S2xybinsize/2);
s2ybins=(S2_ybin_min+S2xybinsize/2):S2xybinsize:(S2_ybin_max-S2xybinsize/2);


%Calculate corrections matrix and 1 sigma corrections matrix
% if inpaint_on==1;  mean_S2_both_xy(mean_S2_both_xy==0)=nan; mean_S2_both_xy=inpaint_nans(mean_S2_both_xy,3); end;
center=interp2(s2xbins,s2ybins,mean_S2_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_S2_all=center./mean_S2_both_xy; 
norm_S2_all(isinf(norm_S2_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_all(isnan(norm_S2_all))=1;

% if inpaint_on==1;  Sigma_S2_both(Sigma_S2_both==0)=nan; Sigma_S2_both=inpaint_nans(Sigma_S2_both,3); end;
sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_both,x_center,y_center,'spline');%Normalize to the center (x=y=0)
sigma_norm_S2_all=sqrt((sigma_center./mean_S2_both_xy).^2+(Sigma_S2_both.*center./mean_S2_both_xy.^2).^2);
sigma_norm_S2_all(isinf(sigma_norm_S2_all))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S2_all(isnan(sigma_norm_S2_all))=1;%no nan. Set Sigma=1, which is 100% uncertainty


s2_both_center_xyz=center;
s2_both_center_xyz_error=sigma_center;


        
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S1 in XY (after correction for S1 Z-dependence). Both PMT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear hist_bin x

   
    %set up energy bins for the gaussian fit
    bin_s1=8;%5
    s1_bin_min=100;
    s1_bin_max=500;
    x=s1_bin_min:bin_s1:s1_bin_max; %150:bin_s1:500;
    x=x';

%set up matricies to be filled
mean_S1_both_xy = zeros(floor(s1_xbins),floor(s1_ybins));
Sigma_S1_both = zeros(floor(s1_xbins),floor(s1_ybins));
Count_S1_both = zeros(floor(s1_xbins),floor(s1_ybins));

% s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
s1_z_correction=polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);%normalize s1 to 160 ns.

s1_phe_both_z=s1_phe_both.*s1_z_correction;
% s2_phe_both_z=s2_phe_both.*exp(drift_distance_\mus/e_lifetime); %Correct S2 for lifeimte  %s2_phe_both_z[5000,18000]
% s2_phe_both_z=s2_phe_both.*exp(drift_time/e_lifetime_both); %Correct S2 for lifeimte  %s2_phe_both_z[5000,18000]



time_energy_cut = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[s1_bin_min,s1_bin_max]) & inrange(s2_phe_both_z,[1000,s2area_bound_max]) & inrange(s2radius,[rcut_min,rcut_max]) ...
    &inrange(s2_width,[100 1000]) & inrange(s1_width,[30 400]);


%%%%variables for uber turbo boost!
s1_phe_both_z_2=s1_phe_both_z(time_energy_cut);
s2x_2=s2x(time_energy_cut);
s2y_2=s2y(time_energy_cut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for x_bin=S2_xbin_min:S2xybinsize:(S2_xbin_max-S2xybinsize);
      %fprintf('X = %d \n',x_bin)
    for y_bin=S2_ybin_min:S2xybinsize:(S2_ybin_max-S2xybinsize);          
      x_min = x_bin; x_max = x_bin+S2xybinsize;
      y_min = y_bin; y_max = y_bin+S2xybinsize;      
      
        x_count=int32(1+x_bin/S2xybinsize+S2_xbin_max/S2xybinsize);
        y_count=int32(1+y_bin/S2xybinsize+S2_ybin_max/S2xybinsize); 
      
      
      bin_cut = inrange(s2x_2,[x_min,x_max]) & inrange(s2y_2,[y_min,y_max]);            
      
       hist_bin = hist(s1_phe_both_z_2(bin_cut) ,x)'/bin_s1; 
       Count_S1_both(y_count,x_count)=sum(hist_bin)*bin_s1; % X and Y flip ... because MATLAB
       
       %%%%%%%%Truncate variables after binning is complete
        s1_phe_both_z_2=s1_phe_both_z_2(~bin_cut);
        s2x_2=s2x_2(~bin_cut);
        s2y_2=s2y_2(~bin_cut);
        clear bin_cut;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
      if Count_S1_both(y_count,x_count)>350; 
                        
                    amp_start=max(hist_bin);
                    mean_start=sum(hist_bin.*x)/sum(hist_bin);
                    sigma_start=std(hist_bin.*x);                     
                Fit_S1_both_xy= fit(x(2:end-1),hist_bin(2:end-1),'gauss1'); %remove edges
                
                %Peak location (mean)
                mean_S1_both_xy(y_count,x_count)=Fit_S1_both_xy.b1; % X and Y flip ... because MATLAB
                
                %1-sigma of peak position
                Sig_b=Fit_S1_both_xy.c1/sqrt(2)/sqrt(Count_S1_both(y_count,x_count)); %sigma/sqrt(N)
   
                    %Fit error checking
                        if(isnan(Sig_b))
                          mean_S1_both_xy(y_count,x_count)=0;
                         Sig_b=0;
                        end
                                
                    %uncertainty in mean
                    Sigma_S1_both(y_count,x_count)=Sig_b;% X and Y flip ... because MATLAB   
     
      else
          
        mean_S1_both_xy(y_count,x_count)=0;  
       
      end            
    end    
   end
   
   mean_S1_both_xy(isnan(mean_S1_both_xy)) = 0;
   
   
%%Plot the correction S1_XY both %%%%%%%%%%%%

s1xbins=(S1_xbin_min+S1xybinsize/2):S1xybinsize:(S1_xbin_max-S1xybinsize/2);
s1ybins=(S1_ybin_min+S1xybinsize/2):S1xybinsize:(S1_ybin_max-S1xybinsize/2);

%Calculate corrections matrix and 1 sigma corrections matrix
% if inpaint_on==1;  mean_S1_both_xy(mean_S1_both_xy==0)=nan; mean_S1_both_xy=inpaint_nans(mean_S1_both_xy,3); end;
center=interp2(s1xbins,s1ybins,mean_S1_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_S1_all=center./mean_S1_both_xy; 
norm_S1_all(isinf(norm_S1_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_all(isnan(norm_S1_all))=1;

% if inpaint_on==1;  Sigma_S1_both(Sigma_S1_both==0)=nan; Sigma_S1_both=inpaint_nans(Sigma_S1_both,3); end;
sigma_center=interp2(s1xbins,s1ybins,Sigma_S1_both,x_center,y_center,'spline');%Normalize to the center (x=y=0)
sigma_norm_S1_all=sqrt((sigma_center./mean_S1_both_xy).^2+(Sigma_S1_both.*center./mean_S1_both_xy.^2).^2);
sigma_norm_S1_all(isinf(sigma_norm_S1_all))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S1_all(isnan(sigma_norm_S1_all))=1;%no nan. Set Sigma=1, which is 100% uncertainty

s1_both_center_xyz=center;
s1_both_center_xyz_error=sigma_center;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply corrections to 2014 Kr data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear s1_phe_both s2_phe_both
load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p18\MappingParametersFit_Start_2014_CutDown.mat')



% Apply the corrections
%-----------------------------------
% Finding the electron lifetime
%-----------------------------------

%  electron_lifetime = e_lifetime_both;

%-----------------------------------
% Finding the S1 Z-dep values
%-----------------------------------

z_dep_both_values = P_s1_both;
z_dep_par_all = z_dep_both_values;      

%------------------------------------------
% Finding the S2 xy correction map values
%------------------------------------------

   s2_x_bins = s2xbins;  
   s2_y_bins = s2ybins; 
   s2_map_all = norm_S2_all; 
   
%------------------------------------------
% Finding the S1 xy correction map values
%------------------------------------------

      
   s1_x_bins = s1xbins; 
   s1_y_bins = s1ybins; 
   s1_map_all = norm_S1_all; 

   
% %------------------------------------------
% % Finding the S1 xyz correction map values
% %------------------------------------------
   
% 
% if nargin == 3      
%    s1_xyz_x_bins = s1xbins_3D;
%    s1_xyz_y_bins = s1ybins_3D;
%    s1_xyz_z_bins = s1zbins_3D;
%    
% %    xx = size(s1_xyz_x_bins); yy = size(s1_xyz_y_bins); zz = size(s1_xyz_z_bins);
%    
%    s1_xyz_map_all = norm_s1_both_xyz; %lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_both_xyz; s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
%    s1_xyz_map_bottom = norm_s1_bottom_xyz; %lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_bottom_xyz;   s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));
% end



%% Calculating corrections

%% Finding the drift time, x and y associated with S1 - this isn't how they are actually applied in DP, since DP estimates the S1 and here we know the exact S1.  
% It's the best we can do for ease here and it shouldn't matter much


s1_drift_ns = drift_time.*1000;
s1_ref_z_ns=z_center*1000;
s1_x_cm = s2x;
s1_y_cm = s2y;

%--------------------------------------------------------------------------
% Calculate Z corrections
%--------------------------------------------------------------------------

% Calculate S1 Z-correction
s1_z_correction = polyval(z_dep_par_all,s1_ref_z_ns./1000)./polyval(z_dep_par_all,s1_drift_ns./1000);
s1_z_correction(find(isnan(s1_z_correction))) = 1.0;

s2_z_correction = RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);

% % Calculate electron lifetime correction (S2 Z-correction)
% electron_lifetime_correction = exp(drift_time./electron_lifetime);
% electron_lifetime_correction(find(isnan(electron_lifetime_correction))) = 1.0;
%  
%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------

% Calculate S1 XY-only corrections
s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'spline',1);
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

% Calculate S2 XY corrections
s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,s2x,s2y,'spline',1);
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

s1xyz_correction=s1xy_correction.*s1_z_correction;

%% add RQs of correction factors

% d.correction_electron_lifetime = single(electron_lifetime_correction);
d.correction_s1_z_dependence = single(s1_z_correction);
d.correction_s1_xyz_dependence = single(s1xyz_correction); %this is either 3D or 1D+2D
d.correction_s1_xy_dependence = single(s1xy_correction);
d.correction_s2_xy_dependence = single(s2xy_correction);

% d.admin.corrections.electron_lifetime = electron_lifetime;
d.admin.corrections.s1_z_dependence.all = z_dep_par_all;

%--------------------------------------------------------------------------
% Apply Z corrections
%--------------------------------------------------------------------------

s1_phe_both_z=s1_phe_both.*s1_z_correction;
s2_phe_both_z=s2_phe_both.*s2_z_correction;


%% Calculate ER Band shift
s1_phe_both_xyz=s1_phe_both.*s1xyz_correction;
s2_phe_both_xyz=s2_phe_both.*s2_z_correction.*s2xy_correction;

    liquid_cut= inrange(drift_time,[40 300]) & inrange(s2radius,[0,18]) & inrange(s1_phe_both_xyz,[30 1000]) ...
    & (log10(s1_phe_both_xyz)+0.5*log10(s2_phe_both_xyz)>4.2)  & s2_phe_both>3000 & s2_phe_both<60000 ... ;%used to be drift time 100 to 250
    & drift_time<-1.4.*s2radius.^2+400;

kr_s1_c_2014=s1_phe_both_xyz(liquid_cut);
kr_s2_c_2014=s2_phe_both_xyz(liquid_cut);
kr_2014_drift_time=drift_time(liquid_cut);
kr_2014_s2x=s2x(liquid_cut);
kr_2014_s2y=s2y(liquid_cut);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Apply Corrections to 2014 CH3T Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear s1_phe_both s2_phe_both
load('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p18\CH3T_Sep2014_NoCorrections.mat');


% Apply the corrections
%-----------------------------------
% Finding the electron lifetime
%-----------------------------------

%  electron_lifetime = e_lifetime_both;

%  electron_lifetime_2014=e_lifetime_both;
%-----------------------------------
% Finding the S1 Z-dep values
%-----------------------------------

z_dep_both_values = P_s1_both;
z_dep_par_all = z_dep_both_values;      

%------------------------------------------
% Finding the S2 xy correction map values
%------------------------------------------

   s2_x_bins = s2xbins;  
   s2_y_bins = s2ybins; 
   s2_map_all = norm_S2_all; 
   
%------------------------------------------
% Finding the S1 xy correction map values
%------------------------------------------

      
   s1_x_bins = s1xbins; 
   s1_y_bins = s1ybins; 
   s1_map_all = norm_S1_all; 

   
% %------------------------------------------
% % Finding the S1 xyz correction map values
% %------------------------------------------
   
% 
% if nargin == 3      
%    s1_xyz_x_bins = s1xbins_3D;
%    s1_xyz_y_bins = s1ybins_3D;
%    s1_xyz_z_bins = s1zbins_3D;
%    
% %    xx = size(s1_xyz_x_bins); yy = size(s1_xyz_y_bins); zz = size(s1_xyz_z_bins);
%    
%    s1_xyz_map_all = norm_s1_both_xyz; %lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_both_xyz; s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
%    s1_xyz_map_bottom = norm_s1_bottom_xyz; %lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_bottom_xyz;   s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));
% end



%% Calculating corrections

%% Finding the drift time, x and y associated with S1 - this isn't how they are actually applied in DP, since DP estimates the S1 and here we know the exact S1.  
% It's the best we can do for ease here and it shouldn't matter much


s1_drift_ns = drift_time.*1000;
s1_ref_z_ns=z_center*1000;
s1_x_cm = S2_x;
s1_y_cm = S2_y;

%--------------------------------------------------------------------------
% Calculate Z corrections
%--------------------------------------------------------------------------

% Calculate S1 Z-correction
s1_z_correction = polyval(z_dep_par_all,s1_ref_z_ns./1000)./polyval(z_dep_par_all,s1_drift_ns./1000);
s1_z_correction(find(isnan(s1_z_correction))) = 1.0;

s2_z_correction = RK_1DCubicInterp_LinearExtrap(s2_both_norm_z_index,s2_both_norm_z, drift_time);

% % Calculate electron lifetime correction (S2 Z-correction)
% electron_lifetime_correction = exp(drift_time./electron_lifetime);
% electron_lifetime_correction(find(isnan(electron_lifetime_correction))) = 1.0;
 
%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------

% Calculate S1 XY-only corrections
s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'spline',1);
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

% Calculate S2 XY corrections
s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,S2_x,S2_y,'spline',1);
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

s1xyz_correction=s1xy_correction.*s1_z_correction;

%% add RQs of correction factors

% d.correction_electron_lifetime = single(electron_lifetime_correction);
d.correction_s1_z_dependence = single(s1_z_correction);
d.correction_s1_xyz_dependence = single(s1xyz_correction); %this is either 3D or 1D+2D
d.correction_s1_xy_dependence = single(s1xy_correction);
d.correction_s2_xy_dependence = single(s2xy_correction);

% d.admin.corrections.electron_lifetime = electron_lifetime;
d.admin.corrections.s1_z_dependence.all = z_dep_par_all;

%--------------------------------------------------------------------------
% Apply Z corrections
%--------------------------------------------------------------------------

s1_phe_both_z=s1_phe_both.*s1_z_correction;
s2_phe_both_z=s2_phe_both.*s2_z_correction;


%% Calculate ER Band shift
s1_phe_both_xyz=s1_phe_both.*s1xyz_correction;
s2_phe_both_xyz=s2_phe_both.*s2_z_correction.*s2xy_correction;


s1_min=0;
s1_max=120;
  
  %% Redo with r<10 cm


    liquid_cut= clean_cut & inrange(drift_time,[40 300]) & inrange(s2radius_c,[0,18]) & inrange(s1_phe_both_xyz,[s1_min s1_max]) ...
    & (log10(s1_phe_both_xyz)+0.5*log10(s2_phe_both_xyz)<4.2) & s2_phe_both>100  & s2_phe_both<16000 ... ;%used to be drift time 100 to 250
    & drift_time<-1.4.*s2radius.^2+400;

h3_s1_c_2014=s1_phe_both_xyz(liquid_cut);
h3_s2_c_2014=s2_phe_both_xyz(liquid_cut);

 clearvars -except param_counter a1_list b1_list a2_list b2_list total_chi2_list h3_chi2_list kr_chi2_list g1_list g2_list total_chi2_map h3_chi2_map kr_chi2_map chi2_midtobot_list...
    best_chi2_midtotop besta1_midtotop besta2_midtotop bestb1_midtotop bestb2_midtotop best_chi2_midtobot...
    besta1_midtobot besta1_midtotop bestb1_midtobot bestb1_midtotop a1 b1 a2 b2 chi2_midtotop_map chi2_midtobot_map b_counter a_counter...
    besta1_combined bestb1_combined besta2_combined bestb2_combined best_chi2_combined chi2_combined_map chi2_combined_list b1_counter a1_counter a2_counter b2_counter...
    a1_counter b1_counter a2_counter b2_counter c1_list c2_list bestc1_midtobot bestc2_midtobot bestc1_midtotop bestc2_midtotop bestc1_combined bestc2_combined...
    besta1_combined besta2_combined bestb1_combined bestb2_combined besta1 bestb1 bestc1 besta2 bestb2 bestc2 bestg1 bestg2 best_total_chi2 NEST_CH3T_Spectrum inpaint_on...
    best_kr_chi2 best_h3_chi2 electron_lifetime_list electron_lifetime_map a2_start a1_start b1_start b2_start ...
    total_h3_chi2_list total_kr_chi2_list total_chi2_list h3_chi2_2014_list h3_chi2_2015_list kr_chi2_2014_list kr_chi2_2015_list ...
    electron_lifetime_2014_list electron_lifetime_2015_list total_h3_chi2_map total_kr_chi2_map total_chi2_map h3_chi2_2014_map h3_chi2_2015_map ...
    kr_chi2_2014_map kr_chi2_2015_map electron_lifetime_2014_map electron_lifetime_2015_map NEST_CH3T_Spectrum g1 g2 bestg1 bestg2 g1_counter g2_counter ...
    firsta1 firstb1 firsta2 firstb2 firstg1 firstg2 kr_means_2015_map kr_means_2014_map kr_sigma_2015_map kr_sigma_2014_map kr_staterr_2015_map kr_staterr_2014_map ...
    kr_means_2015_list kr_means_2014_list kr_2015_num_events kr_2014_num_events kr_sigma_2015_list kr_sigma_2014_list kr_staterr_2015_list kr_staterr_2014_list ...
    pseudo_lifetime_2014 pseudo_lifetime_2015 kr_s1_c_2015 kr_s2_c_2015 h3_s1_c_2015 h3_s2_c_2015 kr_s1_c_2014 kr_s2_c_2014 h3_s1_c_2014 h3_s2_c_2014 ...
    kr_2015_drift_time kr_2015_s2x kr_2015_s2y  kr_2014_drift_time kr_2014_s2x kr_2014_s2y ...
    kr_chi2_2014_zslice_list kr_chi2_2015_zslice_list kr_chi2_zslice_zcenter_list ...
    full_kr_chi2_2014_list full_kr_chi2_2015_list best_total_kr_chi2 best_total_h3_chi2 alternate_total_chi2_map alternate_total_chi2_list ...
    best_total_alternate_chi2 best_total_alternate_kr_chi2 best_total_alternate_h3_chi2 besta2_alternate bestb2_alternate bestg1_alternate bestg2_alternate bestc2_alternate c2 ...
    full_kr_2014_resolution_list full_kr_2015_resolution_list full_kr_2014_resolution_err_list full_kr_2015_resolution_err_list ...
    full_kr_2014_resolution_map full_kr_2015_resolution_map full_kr_2014_resolution_err_map full_kr_2015_resolution_err_map besta1_alternate besta2_alternate






%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO CHI2 TO ENEGY SPECTRA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g1_counter=1;
for g1=0.095:0.0025:0.125;
  g2_counter=1;
    for g2=16.5:0.1:21;
         h3_energy_2014=(1/73).*(h3_s1_c_2014./g1 + h3_s2_c_2014./g2);
         kr_energy_2014=(1/73).*(kr_s1_c_2014./g1 + kr_s2_c_2014./g2);
         h3_energy_2015=(1/73).*(h3_s1_c_2015./g1 + h3_s2_c_2015./g2);
         kr_energy_2015=(1/73).*(kr_s1_c_2015./g1 + kr_s2_c_2015./g2);
         
         
         
        %make CH3T data spectrum
        [h3_hist_2014,h3_hist_bin_2014]=hist(h3_energy_2014,NEST_CH3T_Spectrum(1,:));
        [h3_hist_2015,h3_hist_bin_2015]=hist(h3_energy_2015,NEST_CH3T_Spectrum(1,:));
        
        %scale the simulation to the data
        scale_factor_2014=sum(h3_hist_2014(17:80))./sum(NEST_CH3T_Spectrum(2,17:80));
        h3_scaled_sim_2014=scale_factor_2014.*NEST_CH3T_Spectrum(2,:);
        
        scale_factor_2015=sum(h3_hist_2015(17:80))./sum(NEST_CH3T_Spectrum(2,17:80));
        h3_scaled_sim_2015=scale_factor_2015.*NEST_CH3T_Spectrum(2,:);

        %Calculate chi2 for h3 spectrum
        h3_chi2_2014=(1/61).*sum( (h3_scaled_sim_2014(17:80) - h3_hist_2014(17:80)).^2 ./(h3_hist_2014(17:80)));
        h3_chi2_2015=(1/61).*sum( (h3_scaled_sim_2015(17:80) - h3_hist_2015(17:80)).^2 ./(h3_hist_2015(17:80)));
        
        dT_step=20;
        kr_mean_2014=zeros(length(40+dT_step/2:dT_step:300-dT_step/2),1);
        kr_mean_2015=zeros(length(40+dT_step/2:dT_step:300-dT_step/2),1);
        dT_centers=zeros(length(40+dT_step/2:dT_step:300-dT_step/2),1);
        
        full_kr_fit_res_2014=fit([0:1:1000].',hist(kr_energy_2014,[0:1:1000]).','gauss1');
        full_kr_fit_res_2015=fit([0:1:1000].',hist(kr_energy_2015,[0:1:1000]).','gauss1');
        full_kr_chi2_2014=((41.55-full_kr_fit_res_2014.b1)./0.2395).^2; %SHOULD BE 0.2395 BASED ON RUN03
        full_kr_chi2_2015=((41.55-full_kr_fit_res_2015.b1)./0.2395).^2;        
        
        full_kr_2014_resolution=(full_kr_fit_res_2014.c1/sqrt(2))/full_kr_fit_res_2014.b1;
        full_kr_2015_resolution=(full_kr_fit_res_2015.c1/sqrt(2))/full_kr_fit_res_2015.b1;
        
        temp_res_confint_2014=confint(full_kr_fit_res_2014,0.68);
        temp_res_confint_2015=confint(full_kr_fit_res_2015,0.68);
        width_err_2014=abs(full_kr_fit_res_2014.c1-temp_res_confint_2014(1,3))/sqrt(2);
        width_err_2015=abs(full_kr_fit_res_2014.c1-temp_res_confint_2015(1,3))/sqrt(2);
        mean_err_2014=abs(full_kr_fit_res_2014.b1-temp_res_confint_2014(1,2));
        mean_err_2015=abs(full_kr_fit_res_2015.b1-temp_res_confint_2015(1,2));

        full_kr_2014_resolution_err=sqrt( (width_err_2014./full_kr_fit_res_2014.b1).^2 + ( (mean_err_2014.*(full_kr_fit_res_2014.c1/sqrt(2)))/(full_kr_fit_res_2014.b1.^2)).^2 );
        full_kr_2015_resolution_err=sqrt( (width_err_2015./full_kr_fit_res_2015.b1).^2 + ( (mean_err_2015.*(full_kr_fit_res_2015.c1/sqrt(2)))/(full_kr_fit_res_2015.b1.^2)).^2 );
        
         clear kr_fit_res_2014 kr_fit_res_2015 kr_conf_2014 kr_conf_2015 dT_max
            kk=1;
            for dT_max=40+dT_step/2:dT_step:300-dT_step/2;
                  clear kr_fit_res_2014 kr_fit_res_2015
                  kr_fit_res_2014=fit([0:1:1000].',hist(kr_energy_2014(inrange(kr_2014_drift_time,[dT_max-dT_step dT_max])),[0:1:1000]).','gauss1');
                  kr_fit_res_2015=fit([0:1:1000].',hist(kr_energy_2015(inrange(kr_2015_drift_time,[dT_max-dT_step dT_max])),[0:1:1000]).','gauss1');     
                  kr_mean_2014(kk)=kr_fit_res_2014.b1;
                  kr_mean_2015(kk)=kr_fit_res_2015.b1;
                  dT_centers(kk)=dT_max-dT_step/2;
                  kk=kk+1;
            end

            clear cut
            cut=kr_mean_2015==0 | kr_mean_2014==0;
            kr_mean_2014(cut)=[];
            kr_mean_2015(cut)=[];
            dT_centers(cut)=[];
            
    
        kr_chi2_2014=((41.55-kr_mean_2014)./0.2395).^2; %SHOULD BE 0.2395 BASED ON RUN03
        kr_chi2_2015=((41.55-kr_mean_2015)./0.2395).^2;        

        total_h3_chi2=(h3_chi2_2014+h3_chi2_2015)/2;

        total_kr_chi2_2014=mean(kr_chi2_2014);
        total_kr_chi2_2015=mean(kr_chi2_2015);
     
        total_kr_chi2=(total_kr_chi2_2014+total_kr_chi2_2015)./2;
        total_chi2=(total_kr_chi2+total_h3_chi2)/2;        
        
        alternate_total_chi2=(sum( (h3_scaled_sim_2014(17:80) - h3_hist_2014(17:80)).^2 ./(h3_hist_2014(17:80)))+sum(kr_chi2_2014)+sum(kr_chi2_2015))/(61+length(kr_chi2_2014)+length(kr_chi2_2015)-4);
 
           if total_chi2 < best_total_chi2;
            best_total_chi2=total_chi2;
            best_total_kr_chi2=total_kr_chi2;
            best_total_h3_chi2=total_h3_chi2;
            besta1=a1;
            bestb1=b1;
%             bestc1=c1;
            besta2=a2;
            bestb2=b2;
%             bestc2=c2;    
            bestg1=g1;
            bestg2=g2;
           end


            if alternate_total_chi2 < best_total_alternate_chi2;
            best_total_alternate_chi2=alternate_total_chi2;
            best_total_alternate_kr_chi2=total_kr_chi2;
            best_total_alternate_h3_chi2=total_h3_chi2;
            besta1_alternate=a1;
            bestb1_alternate=b1;
%             bestc1=c1;
            besta2_alternate=a2;
            bestb2_alternate=b2;
%             bestc2_alternate=c2;    
            bestg1_alternate=g1;
            bestg2_alternate=g2;
           end                   
 

a2_list(param_counter)=a2;
b2_list(param_counter)=b2;
a1_list(param_counter)=a1;
b1_list(param_counter)=b1;
% c2_list(param_counter)=c2;
g1_list(param_counter)=g1;
g2_list(param_counter)=g2;
total_h3_chi2_list(param_counter)=total_h3_chi2;
total_kr_chi2_list(param_counter)=total_kr_chi2;
total_chi2_list(param_counter)=total_chi2;
alternate_total_chi2_list(param_counter)=alternate_total_chi2;
h3_chi2_2014_list(param_counter)=h3_chi2_2014;
h3_chi2_2015_list(param_counter)=h3_chi2_2015;
kr_chi2_2014_list(param_counter)=total_kr_chi2_2014;
kr_chi2_2015_list(param_counter)=total_kr_chi2_2015;
full_kr_chi2_2014_list(param_counter)=full_kr_chi2_2014;
full_kr_chi2_2015_list(param_counter)=full_kr_chi2_2015;
electron_lifetime_2014_list(param_counter)=pseudo_lifetime_2014;
electron_lifetime_2015_list(param_counter)=pseudo_lifetime_2015;

kr_chi2_2014_zslice_list(param_counter,:)=kr_chi2_2014;
kr_chi2_2015_zslice_list(param_counter,:)=kr_chi2_2015;
kr_chi2_zslice_zcenter_list(param_counter,:)=dT_centers;

total_h3_chi2_map(a1_counter,a2_counter,g1_counter,g2_counter)=total_h3_chi2;
total_kr_chi2_map(a1_counter,a2_counter,g1_counter,g2_counter)=total_kr_chi2;
total_chi2_map(a1_counter,a2_counter,g1_counter,g2_counter)=total_chi2;
alternate_total_chi2_map(a1_counter,a2_counter,g1_counter,g2_counter)=alternate_total_chi2;
h3_chi2_2014_map(a1_counter,a2_counter,g1_counter,g2_counter)=h3_chi2_2014;
kr_chi2_2014_map(a1_counter,a2_counter,g1_counter,g2_counter)=total_kr_chi2_2014;
h3_chi2_2015_map(a1_counter,a2_counter,g1_counter,g2_counter)=h3_chi2_2015;
kr_chi2_2015_map(a1_counter,a2_counter,g1_counter,g2_counter)=total_kr_chi2_2015;
electron_lifetime_2014_map(a1_counter,a2_counter,g1_counter,g2_counter)=pseudo_lifetime_2014;
electron_lifetime_2015_map(a1_counter,a2_counter,g1_counter,g2_counter)=pseudo_lifetime_2015;

full_kr_2014_resolution_list(param_counter)=full_kr_2014_resolution;
full_kr_2015_resolution_list(param_counter)=full_kr_2015_resolution;
full_kr_2014_resolution_err_list(param_counter)=full_kr_2014_resolution_err;
full_kr_2015_resolution_err_list(param_counter)=full_kr_2015_resolution_err;

full_kr_2014_resolution_map(a1_counter,a2_counter,g1_counter,g2_counter)=full_kr_2014_resolution;
full_kr_2015_resolution_map(a1_counter,a2_counter,g1_counter,g2_counter)=full_kr_2015_resolution;
full_kr_2014_resolution_err_map(a1_counter,a2_counter,g1_counter,g2_counter)=full_kr_2014_resolution_err;
full_kr_2015_resolution_err_map(a1_counter,a2_counter,g1_counter,g2_counter)=full_kr_2015_resolution_err;

param_counter=param_counter+1;
g2_counter=g2_counter+1;
    end
g1_counter=g1_counter+1;
end


catch 
 %Tried the b2 loop, but failed   
sprintf('Caught Error, so skipped. a2=%0.3f b2=%0.3f', a2, b2)

a2_list(param_counter)=a2;
b2_list(param_counter)=b2;
a1_list(param_counter)=a1;
b1_list(param_counter)=b1;
g1_list(param_counter)=-1;
g2_list(param_counter)=-1;
total_h3_chi2_list(param_counter)=-1;
total_kr_chi2_list(param_counter)=-1;
total_chi2_list(param_counter)=-1;
alternate_total_chi2_list(param_counter)=-1;
h3_chi2_2014_list(param_counter)=-1;
h3_chi2_2015_list(param_counter)=-1;
kr_chi2_2014_list(param_counter)=-1;
kr_chi2_2015_list(param_counter)=-1;
electron_lifetime_2014_list(param_counter)=-1;
electron_lifetime_2015_list(param_counter)=-1;
full_kr_chi2_2014_list(param_counter)=-1;
full_kr_chi2_2015_list(param_counter)=-1;

kr_chi2_2014_zslice_list(param_counter,:)=-1;
kr_chi2_2015_zslice_list(param_counter,:)=-1;
kr_chi2_zslice_zcenter_list(param_counter,:)=-1;


total_h3_chi2_map(a1_counter,a2_counter,:,:)=-1;
total_kr_chi2_map(a1_counter,a2_counter,:,:)=-1;
total_chi2_map(a1_counter,a2_counter,:,:)=-1;
alternate_total_chi2_map(a1_counter,a2_counter,:,:)=-1;
h3_chi2_2014_map(a1_counter,a2_counter,:,:)=-1;
kr_chi2_2014_map(a1_counter,a2_counter,:,:)=-1;
h3_chi2_2015_map(a1_counter,a2_counter,:,:)=-1;
kr_chi2_2015_map(a1_counter,a2_counter,:,:)=-1;
electron_lifetime_2014_map(a1_counter,a2_counter,:,:)=-1;
electron_lifetime_2015_map(a1_counter,a2_counter,:,:)=-1;

full_kr_2014_resolution_list(param_counter)=-1;
full_kr_2015_resolution_list(param_counter)=-1;
full_kr_2014_resolution_err_list(param_counter)=-1;
full_kr_2015_resolution_err_list(param_counter)=-1;

full_kr_2014_resolution_map(a1_counter,a2_counter,:,:)=-1;
full_kr_2015_resolution_map(a1_counter,a2_counter,:,:)=-1;
full_kr_2014_resolution_err_map(a1_counter,a2_counter,:,:)=-1;
full_kr_2015_resolution_err_map(a1_counter,a2_counter,:,:)=-1;

param_counter=param_counter+1;
end

g1_counter=1;
% g1=firstg1;

%Done with b2, so save and incriment counter
clearvars -except param_counter a1_list b1_list a2_list b2_list total_chi2_list h3_chi2_list kr_chi2_list g1_list g2_list total_chi2_map h3_chi2_map kr_chi2_map chi2_midtobot_list...
    best_chi2_midtotop besta1_midtotop besta2_midtotop bestb1_midtotop bestb2_midtotop best_chi2_midtobot...
    besta1_midtobot besta1_midtotop bestb1_midtobot bestb1_midtotop a1 b1 a2 b2 chi2_midtotop_map chi2_midtobot_map b_counter a_counter...
    besta1_combined bestb1_combined besta2_combined bestb2_combined best_chi2_combined chi2_combined_map chi2_combined_list b1_counter a1_counter a2_counter b2_counter...
    a1_counter b1_counter a2_counter b2_counter c1_list c2_list bestc1_midtobot bestc2_midtobot bestc1_midtotop bestc2_midtotop bestc1_combined bestc2_combined...
    besta1_combined besta2_combined bestb1_combined bestb2_combined besta1 bestb1 bestc1 besta2 bestb2 bestc2 bestg1 bestg2 best_total_chi2 NEST_CH3T_Spectrum inpaint_on...
    best_kr_chi2 best_h3_chi2 electron_lifetime_list electron_lifetime_map a2_start a1_start b1_start b2_start ...
    total_h3_chi2_list total_kr_chi2_list total_chi2_list h3_chi2_2014_list h3_chi2_2015_list kr_chi2_2014_list kr_chi2_2015_list ...
    electron_lifetime_2014_list electron_lifetime_2015_list total_h3_chi2_map total_kr_chi2_map total_chi2_map h3_chi2_2014_map h3_chi2_2015_map ...
    kr_chi2_2014_map kr_chi2_2015_map electron_lifetime_2014_map electron_lifetime_2015_map NEST_CH3T_Spectrum g1 g2 bestg1 bestg2 g1_counter g2_counter ...
    firsta1 firstb1 firsta2 firstb2 firstg1 firstg2 kr_means_2015_map kr_means_2014_map kr_sigma_2015_map kr_sigma_2014_map kr_staterr_2015_map kr_staterr_2014_map ...
    kr_means_2015_list kr_means_2014_list kr_2015_num_events kr_2014_num_events kr_sigma_2015_list kr_sigma_2014_list kr_staterr_2015_list kr_staterr_2014_list ...
    pseudo_lifetime_2014 pseudo_lifetime_2015 kr_2015_drift_time kr_2015_s2x kr_2015_s2y  kr_2014_drift_time kr_2014_s2x kr_2014_s2y ...
    kr_chi2_2014_zslice_list kr_chi2_2015_zslice_list kr_chi2_zslice_zcenter_list ...
    full_kr_chi2_2014_list full_kr_chi2_2015_list best_total_kr_chi2 best_total_h3_chi2 alternate_total_chi2_map alternate_total_chi2_list ...
    best_total_alternate_chi2 best_total_alternate_kr_chi2 best_total_alternate_h3_chi2 besta2_alternate bestb2_alternate bestg1_alternate bestg2_alternate bestc2_alternate c2 ...
    full_kr_2014_resolution_list full_kr_2015_resolution_list full_kr_2014_resolution_err_list full_kr_2015_resolution_err_list ...
    full_kr_2014_resolution_map full_kr_2015_resolution_map full_kr_2014_resolution_err_map full_kr_2015_resolution_err_map besta1_alternate besta2_alternate


save('C:\Program Files\MATLAB\R2012a\bin\Run04_ERBand\Kr2p18\DirectMeasurement\S1aS1bParamScan_EnergyChi2_AllFloat_Linear.mat') 

    
    a2_counter=a2_counter+1;    

    end
  
  a1_counter=a1_counter+1;    
  a2_counter=1;
 end
    
