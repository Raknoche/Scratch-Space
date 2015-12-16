%The desired outputs are the "field" corrections,s2z_normalization_fit, norm_S2_xy (and s2_xy_xbins s2_xy_ybins),
%mean_fit_s1z and norm_s1_xy (and s1_xy_xbins s1_xy_ybins)

%Apply corrections with
% s2_phe_bottom_xyz_z=s2_phe_bottom_xyz.*polyval(s2z_normalization_fit,0)./polyval(s2z_normalization_fit,drift_time)
% s2_phe_bottom_xyz_xyz=s2_phe_bottom_xyz_z.*interp2(s2_xy_xbins,s2_xy_ybins,norm_S2_xy,s2x,s2y,'spline');
% s1_phe_both_xyz_z=s1_phe_both_xyz.*polyval(mean_fit_s1z,(det_edge-4)/2)./polyval(mean_fit_s1z,drift_time);
% s1_phe_both_xyz_xyz=s1_phe_both_xyz_z.*interp2(s1_xy_xbins,s1_xy_ybins,norm_S1_both,s2x,s2y,'spline');


%%%%%%%%%%%%%%%%%
%First measure the "true" lifetime from CH3T, assuming no field dependence
%%%%%%%%%%%%%%%%%

load('CH3T_Sep2014_TAGCB_1_7_1_085_2_DP21');


g1=0.123;%1.01
g2_bot=6.47; %DP 2.0 (no vuv gain correction) -- 5.63;
g2_both=15.0; %DP 2.0 -- 13.75;

s2radius=(s2x.^2+s2y.^2).^(1/2);
s2_phe_bottom=s2_phe_bottom(s2radius<25);
s1_phe_both=s1_phe_both(s2radius<25);
s1_phe_bottom=s1_phe_bottom(s2radius<25);
s2_phe_both=s2_phe_both(s2radius<25);
s2x=s2x(s2radius<25);
s2y=s2y(s2radius<25);
drift_time=drift_time(s2radius<25);
s2radius=(s2x.^2+s2y.^2).^(1/2);

h3_s1_phe_both=s1_phe_both(s2radius<25 & (log10(s1_phe_both)+0.5*log10(s2_phe_bottom) < 3.8));
h3_s1_phe_bottom=s1_phe_bottom(s2radius<25 & (log10(s1_phe_both)+0.5*log10(s2_phe_bottom) < 3.8));
h3_s2_phe_bottom=s2_phe_bottom(s2radius<25 & (log10(s1_phe_both)+0.5*log10(s2_phe_bottom) < 3.8));
h3_s2_phe_both=s2_phe_both(s2radius<25 & (log10(s1_phe_both)+0.5*log10(s2_phe_bottom) < 3.8));
h3_drift_time=drift_time(s2radius<25 & (log10(s1_phe_both)+0.5*log10(s2_phe_bottom) < 3.8));
h3_radius=s2radius(s2radius<25 & (log10(s1_phe_both)+0.5*log10(s2_phe_bottom) < 3.8));
h3_energy=(1/73)*(h3_s1_phe_both./g1 + h3_s2_phe_bottom./g2_bot);
h3_energy_both=(1/73)*(h3_s1_phe_both./g1 + h3_s2_phe_both./g2_both);
h3_x=s2x(s2radius<25 & (log10(s1_phe_both)+0.5*log10(s2_phe_bottom) < 3.8));
h3_y=s2y(s2radius<25 & (log10(s1_phe_both)+0.5*log10(s2_phe_bottom) < 3.8));



s2_limit=3500; %For TAGCB=1,7,1,8.5,2 s2_limit=2500, FOR TAGCB=1,7,2,10,2 s2_limit=1200
s2_min=400; %For TAGCB=1,7,1,8.5,2 = s2_min=100, FOR TAGCB=1,7,2,10,2 s2_min=s2_min

dT_step=15;
dT_max=30+dT_step;
i=1;

for dT_max=30+dT_step:dT_step:310-dT_step;
h3_mean(i)=histfitlandau(h3_s2_phe_bottom(inrange(h3_drift_time,[dT_max-dT_step dT_max])), 30, 400, (3000+(285-dT_max)*2.0408));
dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

h3_mean_fit=fit(dT_loc.',h3_mean.','exp1');
fit_points=30:2:310;
h3_fit_line=h3_mean_fit.a.*exp(h3_mean_fit.b.*fit_points);

figure
scatter(dT_loc,h3_mean,'.k');
hold on;
plot(fit_points,h3_fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('Uncorrected CH3T S2Bot Peak (phe)'); myfigview(16);
title('Electron Lifetime Extracted from Landau Fits to CH3T');
e_lifetime_landau=-1/h3_mean_fit.b;
b=confint(h3_mean_fit, 0.683);%1 sigma matrix
sigma_e_lifetime=(1/b(1,2)-1/b(2,2))/2;
legend('S2\_bottom Mean', strcat( '\lambda= ', num2str(e_lifetime_landau,4), ' \pm ', num2str(sigma_e_lifetime,2), ' \mus' )...
        ,'location','northeast');

    
%Show some slices for Landau
dT_step=60;
figure
hold on;
colors={'r','b','g'};
colors2={[1 0.6 0.6],[0.6 0.6 1],[0.6 1 0.6]};
i=1;
    while i <3
for dT_max=30+dT_step:dT_step:310-dT_step;
val=h3_s2_phe_bottom(inrange(h3_drift_time,[dT_max-dT_step dT_max]));
passo=30;
inizio=s2_min;
fine=(3000+(285-dT_max)*2.0408);%s2_limit;
[y,x]=hist_(val,passo,min(val),fine,fine,1);
step(val,[inizio:passo:fine],colors2{i})
x=x+passo/2;
inz=find(x<inizio);
    if isempty(inz)
        inz=1;
    end;
y(1:inz(end))=0;
a0=max(y);
mpv0=min(x(find(y==a0)));
sigma0=std(val);
ftype=fittype('a*landau(x,mpv,sigma,0)');
fitres=fit(x,y,ftype,'StartPoint',[a0 mpv0 sigma0]);
bound=confint(fitres);
xfit=x(1):(passo/10):x(end);
yfit=fitres(xfit);
temp=line(xfit,yfit,'LineWidth',2,'LineStyle','-','Color',colors{i});
hAnnotation = get(temp,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off')
picco=max(yfit);
mpv=xfit(find(yfit==picco));
mpv=mpv(end);
mv=mean(val);
temp=line([mpv mpv],[0 picco],'LineWidth',2,'LineStyle','--','Color',colors{i});
hAnnotation = get(temp,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off')
i=i+1;
    end
    end
legend('30-90 uSec','90-150 uSec','150-210 uSec');
xlabel('Uncorrected S2Bot (phe)');ylabel('Counts'); title('Landau Fits to S2Bot Spectra');
myfigview(16);

electron_lifetime=e_lifetime_landau;
h3_s2_phe_bottom_z=h3_s2_phe_bottom.*exp(h3_drift_time./electron_lifetime);


%Showing that the correction worked
dT_step=15;
dT_max=30+dT_step;
i=1;

for dT_max=30+dT_step:dT_step:310-dT_step;
h3_mean_test(i)=histfitlandau(h3_s2_phe_bottom_z(inrange(h3_drift_time,[dT_max-dT_step dT_max])), 30, 400, (3000+(285-dT_max)*2.0408));
dT_loc_test(i)=dT_max-dT_step/2;
i=i+1;
end
figure
scatter(dT_loc_test,h3_mean_test,'.k');

% clearvars -except h3_s2_phe_bottom_z h3_radius h3_drift_time s2_limit s2_min h3_s2_phe_bottom e_lifetime_landau h3_energy h3_x h3_y

%% Get lifetime from S2_both

dT_step=15;
dT_max=30+dT_step;
i=1;

for dT_max=30+dT_step:dT_step:310-dT_step;
h3_mean(i)=histfitlandau(h3_s2_phe_both(inrange(h3_drift_time,[dT_max-dT_step dT_max])), 30, 400, (5000+(285-dT_max)*2.0408));
dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

h3_mean_fit=fit(dT_loc.',h3_mean.','exp1');
fit_points=30:2:310;
h3_fit_line=h3_mean_fit.a.*exp(h3_mean_fit.b.*fit_points);

figure
scatter(dT_loc,h3_mean,'.k');
hold on;
plot(fit_points,h3_fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('Uncorrected CH3T S2Both Peak (phe)'); myfigview(16);
title('Electron Lifetime Extracted from Landau Fits to CH3T');
e_lifetime_landau=-1/h3_mean_fit.b;
b=confint(h3_mean_fit, 0.683);%1 sigma matrix
sigma_e_lifetime=(1/b(1,2)-1/b(2,2))/2;
legend('S2\_both Mean', strcat( '\lambda= ', num2str(e_lifetime_landau,4), ' \pm ', num2str(sigma_e_lifetime,2), ' \mus' )...
        ,'location','northeast');

    
%Show some slices for Landau
dT_step=60;
figure
hold on;
colors={'r','b','g'};
colors2={[1 0.6 0.6],[0.6 0.6 1],[0.6 1 0.6]};
i=1;
    while i <3
for dT_max=30+dT_step:dT_step:310-dT_step;
val=h3_s2_phe_both(inrange(h3_drift_time,[dT_max-dT_step dT_max]));
passo=30;
inizio=s2_min;
fine=(5000+(285-dT_max)*2.0408);%s2_limit;
[y,x]=hist_(val,passo,min(val),fine,fine,1);
step(val,[inizio:passo:fine],colors2{i})
x=x+passo/2;
inz=find(x<inizio);
    if isempty(inz)
        inz=1;
    end;
y(1:inz(end))=0;
a0=max(y);
mpv0=min(x(find(y==a0)));
sigma0=std(val);
ftype=fittype('a*landau(x,mpv,sigma,0)');
fitres=fit(x,y,ftype,'StartPoint',[a0 mpv0 sigma0]);
bound=confint(fitres);
xfit=x(1):(passo/10):x(end);
yfit=fitres(xfit);
temp=line(xfit,yfit,'LineWidth',2,'LineStyle','-','Color',colors{i});
hAnnotation = get(temp,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off')
picco=max(yfit);
mpv=xfit(find(yfit==picco));
mpv=mpv(end);
mv=mean(val);
temp=line([mpv mpv],[0 picco],'LineWidth',2,'LineStyle','--','Color',colors{i});
hAnnotation = get(temp,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off')
i=i+1;
    end
    end
legend('30-90 uSec','90-150 uSec','150-210 uSec');
xlabel('Uncorrected S2Both (phe)');ylabel('Counts'); title('Landau Fits to S2Both Spectra');
myfigview(16);

electron_lifetime_both=e_lifetime_landau;
h3_s2_phe_both_z=h3_s2_phe_both.*exp(h3_drift_time./electron_lifetime_both);


%Showing that the correction worked
dT_step=15;
dT_max=30+dT_step;
i=1;

for dT_max=30+dT_step:dT_step:310-dT_step;
h3_mean_test(i)=histfitlandau(h3_s2_phe_both_z(inrange(h3_drift_time,[dT_max-dT_step dT_max])), 30, 400, (5000+(285-dT_max)*2.0408));
dT_loc_test(i)=dT_max-dT_step/2;
i=i+1;
end
figure
scatter(dT_loc_test,h3_mean_test,'.k');

%%
%%%%%%%%%%%%%%%%%%%
% Get S2 XY geomtric corrections based on CH3T
%%%%%%%%%%%%%%%%%%%
clear x2 hist_bin

xx_step=3;
yy_step=3;
%set up matricies to be filled
h3_s2_bottom_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s2_phe_bottom_z(inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 600;
            h3_s2_bottom_xymeans(j,i)=histfitlandau(h3_s2_phe_bottom_z(inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 30, 400, 3000);
             
           else 
          
        h3_s2_bottom_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s2_bottom_xymeans(isnan(h3_s2_bottom_xymeans)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;

old_s2xbins=s2x_bins;
old_s2ybins=s2y_bins;
    
    color_range_max=max(max(h3_s2_bottom_xymeans));
    color_range_min=min(min(h3_s2_bottom_xymeans(h3_s2_bottom_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_bottom_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s2_bottom_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);


    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,h3_s2_bottom_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S2_bot=center./h3_s2_bottom_xymeans; 
norm_S2_bot(isinf(norm_S2_bot))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_bot(isnan(norm_S2_bot))=1;
% 
% sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_bottom,0,0,'spline',1);%Normalize to the center (x=y=0)
% sigma_norm_S2_bot=sqrt((sigma_center./mean_S2_bottom_xy).^2+(Sigma_S2_bottom.*center./mean_S2_bottom_xy.^2).^2);
% sigma_norm_S2_bot(isinf(sigma_norm_S2_bot))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
% sigma_norm_S2_bot(isnan(sigma_norm_S2_bot))=1;%no nan. Set Sigma=1, which is 100% uncertainty
% 
h3_s2_phe_bottom_xyz=h3_s2_phe_bottom_z.*interp2(s2x_bins,s2y_bins,norm_S2_bot,h3_x,h3_y,'spline');

%Showing that the correction worked
xx_step=2;
yy_step=3;
h3_s2_bottom_xymeans_test = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s2_phe_bottom_z(inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 600;
            h3_s2_bottom_xymeans_test(j,i)=histfitlandau(h3_s2_phe_bottom_xyz(inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 30, 400, 3000);
             
           else 
          
        h3_s2_bottom_xymeans_test(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s2_bottom_xymeans_test(isnan(h3_s2_bottom_xymeans_test)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
    
    color_range_max=max(max(h3_s2_bottom_xymeans_test));
    color_range_min=min(min(h3_s2_bottom_xymeans_test(h3_s2_bottom_xymeans_test>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_bottom_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s2_bottom_xymeans_test,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('XYZ Corrected S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
%% Getting S2 XY for using S2 both

%%%%%%%%%%%%%%%%%%%
% Get S2 XY geomtric corrections based on CH3T
%%%%%%%%%%%%%%%%%%%
clear x2 hist_bin

xx_step=3;
yy_step=3;
%set up matricies to be filled
h3_s2_both_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s2_phe_both_z(inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 600;
            h3_s2_both_xymeans(j,i)=histfitlandau(h3_s2_phe_both_z(inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 40, 400, 5000);
             
           else 
          
        h3_s2_both_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s2_both_xymeans(isnan(h3_s2_both_xymeans)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;

old_s2xbins_both=s2x_bins;
old_s2ybins_both=s2y_bins;
    
    color_range_max=max(max(h3_s2_both_xymeans));
    color_range_min=min(min(h3_s2_both_xymeans(h3_s2_both_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_both_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s2_both_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S2 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);

    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,h3_s2_both_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S2_both=center./h3_s2_both_xymeans; 
norm_S2_both(isinf(norm_S2_both))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_both(isnan(norm_S2_both))=1;

h3_s2_phe_both_xyz=h3_s2_phe_both_z.*interp2(s2x_bins,s2y_bins,norm_S2_both,h3_x,h3_y,'spline');

%Showing that the correction worked
xx_step=3;
yy_step=3;
h3_s2_both_xymeans_test = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s2_phe_both_z(inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 600;
            h3_s2_both_xymeans_test(j,i)=histfitlandau(h3_s2_phe_both_xyz(inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 40, 400, 5000);
             
           else 
          
        h3_s2_both_xymeans_test(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s2_both_xymeans_test(isnan(h3_s2_both_xymeans_test)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
    
    color_range_max=max(max(h3_s2_both_xymeans_test));
    color_range_min=min(min(h3_s2_both_xymeans_test(h3_s2_both_xymeans_test>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_both_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s2_both_xymeans_test,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('XYZ Corrected S2 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);



%%
%%%%%%%
% GET S1 Z Dependence
%%%%%%%%%%%

dT_step=15;
dT_max=30+dT_step;
i=1;

for dT_max=30+dT_step:dT_step:310-dT_step;
% h3_s1_mean(i)=histfitlandau(h3_s1_phe_both(inrange(h3_energy,[3 5]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])), 0.5, 1, 100);
h3_s1_mean_fit=fit([0:0.5:50]',hist(h3_s1_phe_both(inrange(h3_energy,[3 4.5]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])),[0:0.5:50])','gauss1'); 
h3_s1_mean(i)=h3_s1_mean_fit.b1;
s1z_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

h3_mean_fit_s1z=polyfit(s1z_dT_loc,h3_s1_mean,2);
%     sigma_P_s1_both = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';    
fit_points=30:2:310;
h3_fit_line=polyval(h3_mean_fit_s1z,fit_points);

figure
scatter(s1z_dT_loc,h3_s1_mean,'.k');
hold on;
plot(fit_points,h3_fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('Uncorrected CH3T S1 Peak (phe)'); myfigview(16);
title('S1 Z Dependence Extracted from Landau Fits to CH3T');

%     
% %Show some slices for Landau
% dT_step=60;
% figure
% hold on;
% colors={'r','b','g'};
% colors2={[1 0.6 0.6],[0.6 0.6 1],[0.6 1 0.6]};
% i=1;
%     while i <3
% for dT_max=30+dT_step:dT_step:310-dT_step;
% val=h3_s1_phe_both(inrange(h3_energy,[2.5 4]) & inrange(h3_drift_time,[dT_max-dT_step dT_max]));
% passo=0.5;
% inizio=1;
% fine=100;%s2_limit;
% [y,x]=hist_(val,passo,min(val),fine,fine,1);
% step(val,[inizio:passo:fine],colors2{i})
% x=x+passo/2;
% inz=find(x<inizio);
%     if isempty(inz)
%         inz=1;
%     end;
% y(1:inz(end))=0;
% a0=max(y);
% mpv0=min(x(find(y==a0)));
% sigma0=std(val);
% ftype=fittype('a*landau(x,mpv,sigma,0)');
% fitres=fit(x,y,ftype,'StartPoint',[a0 mpv0 sigma0]);
% bound=confint(fitres);
% xfit=x(1):(passo/10):x(end);
% yfit=fitres(xfit);
% temp=line(xfit,yfit,'LineWidth',2,'LineStyle','-','Color',colors{i});
% hAnnotation = get(temp,'Annotation');
% hLegendEntry = get(hAnnotation','LegendInformation');
% set(hLegendEntry,'IconDisplayStyle','off')
% picco=max(yfit);
% mpv=xfit(find(yfit==picco));
% mpv=mpv(end);
% mv=mean(val);
% temp=line([mpv mpv],[0 picco],'LineWidth',2,'LineStyle','--','Color',colors{i});
% hAnnotation = get(temp,'Annotation');
% hLegendEntry = get(hAnnotation','LegendInformation');
% set(hLegendEntry,'IconDisplayStyle','off')
% i=i+1;
%     end
%     end
% legend('30-90 uSec','90-150 uSec','150-210 uSec');
% xlabel('Uncorrected S1 (phe)');ylabel('Counts'); title('Landau Fits to S1 Spectra');
% myfigview(16);

det_edge=320;
h3_s1_phe_both_z=h3_s1_phe_both.*polyval(h3_mean_fit_s1z,(det_edge-4)/2)./polyval(h3_mean_fit_s1z,h3_drift_time);

% Showing that the correction worked

dT_step=15;
dT_max=30+dT_step;
i=1;
for dT_max=30+dT_step:dT_step:310-dT_step;
% h3_s1_mean_test(i)=histfitlandau(h3_s1_phe_both_z(inrange(h3_energy,[3 5]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])), 0.5, 1, 100);
h3_s1_mean_fit=fit([0:0.5:50]',hist(h3_s1_phe_both_z(inrange(h3_energy,[3 4.5]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])),[0:0.5:50])','gauss1'); 
h3_s1_mean_test(i)=h3_s1_mean_fit.b1;
i=i+1;
end

figure
scatter(s1z_dT_loc,h3_s1_mean_test,'.k')

%% S1 Z dependence using S1 bot

%%
%%%%%%%
% GET S1 Z Dependence
%%%%%%%%%%%

dT_step=15;
dT_max=30+dT_step;
i=1;

for dT_max=30+dT_step:dT_step:det_edge-dT_step;
% h3_s1_mean(i)=histfitlandau(h3_s1_phe_both(inrange(h3_energy,[3 5]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])), 0.5, 1, 100);
h3_s1_mean_fit=fit([0:0.5:50]',hist(h3_s1_phe_bottom(inrange(h3_energy,[3 4.5]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])),[0:0.5:50])','gauss1'); 
h3_s1_mean(i)=h3_s1_mean_fit.b1;
s1z_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

h3_mean_fit_s1z_bottom=polyfit(s1z_dT_loc,h3_s1_mean,2);
%     sigma_P_s1_both = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';    
fit_points=30:2:310;
h3_fit_line=polyval(h3_mean_fit_s1z_bottom,fit_points);

figure
scatter(s1z_dT_loc,h3_s1_mean,'.k');
hold on;
plot(fit_points,h3_fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('Uncorrected CH3T S1 Peak (phe)'); myfigview(16);
title('S1 Z Dependence Extracted from Landau Fits to CH3T');

det_edge=320;
h3_s1_phe_bottom_z=h3_s1_phe_bottom.*polyval(h3_mean_fit_s1z_bottom,(det_edge-4)/2)./polyval(h3_mean_fit_s1z_bottom,h3_drift_time);

% Showing that the correction worked

dT_step=15;
dT_max=30+dT_step;
i=1;
for dT_max=30+dT_step:dT_step:det_edge-dT_step;
% h3_s1_mean_test(i)=histfitlandau(h3_s1_phe_both_z(inrange(h3_energy,[3 5]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])), 0.5, 1, 100);
h3_s1_mean_fit=fit([0:0.5:50]',hist(h3_s1_phe_bottom_z(inrange(h3_energy,[3 4.5]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])),[0:0.5:50])','gauss1'); 
h3_s1_mean_test(i)=h3_s1_mean_fit.b1;
i=i+1;
end

figure
scatter(s1z_dT_loc,h3_s1_mean_test,'.k')


%%
%%%%%%%%%%%%%%%%%%%
% Get S1 XY geomtric corrections based on CH3T
%%%%%%%%%%%%%%%%%%%
clear x2 hist_bin

xx_step=3;
yy_step=3;
%set up matricies to be filled
h3_s1_both_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s1_phe_both_z(inrange(h3_energy,[3 4.5]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 100;
            h3_s1_both_xymeans_fit=fit([0:0.5:50]',hist(h3_s1_phe_both_z(inrange(h3_energy,[3 4.5]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])),[0:0.5:50])','gauss1'); %histfitlandau(h3_s1_phe_both_z(inrange(h3_energy,[2 18]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 1, 0.5, 100);
            h3_s1_both_xymeans(j,i)=h3_s1_both_xymeans_fit.b1;
           else 
          
        h3_s1_both_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s1_both_xymeans(isnan(h3_s1_both_xymeans)) = 0;
     
   
%%Plot the correction S1_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
old_s1xbins=s2x_bins;
old_s1ybins=s2y_bins;
    
    color_range_max=max(max(h3_s1_both_xymeans));
    color_range_min=min(min(h3_s1_both_xymeans(h3_s1_both_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s1_xy_both_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s1_both_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S1 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
% 
% %%Plot the 1 sigma of the mean
% 
%     color_range_max=max(max(Sigma_S2_bottom));
%     color_range_min=min(min(Sigma_S2_bottom(Sigma_S2_bottom>0)));
%     vc_step=(color_range_max-color_range_min)/50;
% 
%     vc=color_range_min:vc_step:color_range_max;
% 
%     sigma_s2_xy_bottom_fig = figure;
%     contourf(s2xbins,s2ybins,Sigma_S2_bottom,vc,'LineColor','none');
%     xlabel('x (cm)','fontsize', 18);
%     ylabel('y (cm)','fontsize', 18);
%     title(strcat(file_id_cp, '. 1-Sigma S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
%     caxis([color_range_min color_range_max])
%     g=colorbar;
%     ylabel(g,'Phe','FontSize',16)
%     set(g,'FontSize',18)
%     hold on
%     LUXPlotTopPMTs
%     axis([-25 25 -25 25])
%     myfigview(16);

    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,h3_s1_both_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S1_both=center./h3_s1_both_xymeans; 
norm_S1_both(isinf(norm_S1_both))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_both(isnan(norm_S1_both))=1;
% 
% sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_bottom,0,0,'spline',1);%Normalize to the center (x=y=0)
% sigma_norm_S2_bot=sqrt((sigma_center./mean_S2_bottom_xy).^2+(Sigma_S2_bottom.*center./mean_S2_bottom_xy.^2).^2);
% sigma_norm_S2_bot(isinf(sigma_norm_S2_bot))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
% sigma_norm_S2_bot(isnan(sigma_norm_S2_bot))=1;%no nan. Set Sigma=1, which is 100% uncertainty
% 
h3_s1_phe_both_xyz=h3_s1_phe_both_z.*interp2(s2x_bins,s2y_bins,norm_S1_both,h3_x,h3_y,'spline');


%Prove the correction worked
xx_step=3;
yy_step=3;
h3_s1_both_xymeans_test = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s1_phe_both_z(inrange(h3_energy,[3 4.5]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 100;
            h3_s1_both_xymeans_fit=fit([0:0.5:50]',hist(h3_s1_phe_both_xyz(inrange(h3_energy,[3 4.5]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])),[0:0.5:50])','gauss1'); %histfitlandau(h3_s1_phe_both_z(inrange(h3_energy,[2 18]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 1, 0.5, 100);
            h3_s1_both_xymeans_test(j,i)=h3_s1_both_xymeans_fit.b1;
           else 
          
        h3_s1_both_xymeans_test(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s1_both_xymeans_test(isnan(h3_s1_both_xymeans_test)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
    
    color_range_max=max(max(h3_s1_both_xymeans_test));
    color_range_min=min(min(h3_s1_both_xymeans_test(h3_s1_both_xymeans_test>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s1_xy_both_test_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s1_both_xymeans_test,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('Corrected S1 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);

%%
%%%%%%%%%%%%%%%%%%%
% Get S1 XY geomtric corrections based on CH3T using bot only
%%%%%%%%%%%%%%%%%%%
clear x2 hist_bin

xx_step=3;
yy_step=3;
%set up matricies to be filled
h3_s1_bot_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s1_phe_bottom_z(inrange(h3_energy,[3 4.5]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 100;
            h3_s1_bottom_xymeans_fit=fit([0:0.5:50]',hist(h3_s1_phe_bottom_z(inrange(h3_energy,[3 4.5]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])),[0:0.5:50])','gauss1'); %histfitlandau(h3_s1_phe_bottom_z(inrange(h3_energy,[2 18]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 1, 0.5, 100);
            h3_s1_bottom_xymeans(j,i)=h3_s1_bottom_xymeans_fit.b1;
           else 
          
        h3_s1_bottom_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s1_bottom_xymeans(isnan(h3_s1_bottom_xymeans)) = 0;
     
   
%%Plot the correction S1_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
old_s1xbins_bottom=s2x_bins;
old_s1ybins_bottom=s2y_bins;
    
    color_range_max=max(max(h3_s1_bottom_xymeans));
    color_range_min=min(min(h3_s1_bottom_xymeans(h3_s1_bottom_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s1_xy_bottom_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s1_bottom_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S1 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
% 
% %%Plot the 1 sigma of the mean
% 
%     color_range_max=max(max(Sigma_S2_bottom));
%     color_range_min=min(min(Sigma_S2_bottom(Sigma_S2_bottom>0)));
%     vc_step=(color_range_max-color_range_min)/50;
% 
%     vc=color_range_min:vc_step:color_range_max;
% 
%     sigma_s2_xy_bottom_fig = figure;
%     contourf(s2xbins,s2ybins,Sigma_S2_bottom,vc,'LineColor','none');
%     xlabel('x (cm)','fontsize', 18);
%     ylabel('y (cm)','fontsize', 18);
%     title(strcat(file_id_cp, '. 1-Sigma S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
%     caxis([color_range_min color_range_max])
%     g=colorbar;
%     ylabel(g,'Phe','FontSize',16)
%     set(g,'FontSize',18)
%     hold on
%     LUXPlotTopPMTs
%     axis([-25 25 -25 25])
%     myfigview(16);

    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,h3_s1_bottom_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S1_bottom=center./h3_s1_bottom_xymeans; 
norm_S1_bottom(isinf(norm_S1_bottom))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_bottom(isnan(norm_S1_bottom))=1;
% 
% sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_bottom,0,0,'spline',1);%Normalize to the center (x=y=0)
% sigma_norm_S2_bot=sqrt((sigma_center./mean_S2_bottom_xy).^2+(Sigma_S2_bottom.*center./mean_S2_bottom_xy.^2).^2);
% sigma_norm_S2_bot(isinf(sigma_norm_S2_bot))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
% sigma_norm_S2_bot(isnan(sigma_norm_S2_bot))=1;%no nan. Set Sigma=1, which is 100% uncertainty
% 
h3_s1_phe_bottom_xyz=h3_s1_phe_bottom_z.*interp2(s2x_bins,s2y_bins,norm_S1_bottom,h3_x,h3_y,'spline');


%Prove the correction worked
xx_step=3;
yy_step=3;
h3_s1_bottom_xymeans_test = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(h3_s1_phe_bottom_z(inrange(h3_energy,[3 4.5]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step]))) > 100;
            h3_s1_bottom_xymeans_fit=fit([0:0.5:50]',hist(h3_s1_phe_bottom_xyz(inrange(h3_energy,[3 4.5]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])),[0:0.5:50])','gauss1'); %histfitlandau(h3_s1_phe_bottom_z(inrange(h3_energy,[2 18]) & inrange(h3_x,[x_bin-xx_step x_bin+xx_step]) & inrange(h3_y,[y_bin-yy_step y_bin+yy_step])), 1, 0.5, 100);
            h3_s1_bottom_xymeans_test(j,i)=h3_s1_bottom_xymeans_fit.b1;
           else 
          
        h3_s1_bottom_xymeans_test(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    h3_s1_bottom_xymeans_test(isnan(h3_s1_bottom_xymeans_test)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
    
    color_range_max=max(max(h3_s1_bottom_xymeans_test));
    color_range_min=min(min(h3_s1_bottom_xymeans_test(h3_s1_bottom_xymeans_test>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s1_xy_bottom_test_fig = figure;
    contourf(s2x_bins,s2y_bins,h3_s1_bottom_xymeans_test,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('Corrected S1 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
    %%

%%%%%%%%%%%%%%%%%%%%%%
%Correct the Kr data with the lifetime found above
%%%%%%%%%%%%%%%%%%%%%%
dir_path='C:\Users\Richard\Desktop\LUX Analysis\';
user_name='RichardKnoche';
dp_version='2.0';
algorithm_name='LUXkrypCal';
use_ccv_flag=1;

folders_to_load{1}='lux10_20140903T1918_cp12334';


rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing','golden'...
   ,'z_corrected_pulse_area_all_phe','xyz_corrected_pulse_area_all_phe'...
   ,'xyz_corrected_pulse_area_bot_phe'...
   ,'x_corrected','y_corrected','selected_s1_s2','top_bottom_ratio','x_cm','y_cm','z_corrected_pulse_area_bot_phe','xyz_corrected_pulse_area_bot_phe','full_evt_area_phe'};
    
path=strcat(dir_path,'/',folders_to_load{1});         
d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);                
                  
%Correcting x and y
 myname = 'Corrections_PositionCorrection';           
 position_correction_path = 'C:\Program Files\MATLAB\R2012a\bin\LuxCode\Trunk\DataProcessing\MatlabModules\Corrections_PositionCorrection\Corrections_PositionCorrection.m';
 IQs = dir([position_correction_path(1:(end-numel(myname)-2)) '*Mercury.mat']);
 table_corrections = load(IQs(4).name);

 d = Corrections_PositionCorrection_FunctionUseThis(d,table_corrections);
%%  Defining cuts and variables

submit=0;
delay_min=0;
grid_size=2; %cm XY plane

rcut_min = 0;
rcut_max = 25;%cm
s1area_bound_min = 20;%100
s1area_bound_max = 600;%600

s2area_bound_min = 100;%200 %change to cut out CH3T if needed
s2area_bound_max = 35000;%30000

min_evts_per_bin = 200;
max_num_bins = 65;

S1xybinsize = grid_size;%cm
S1_xbin_min = -25;
S1_xbin_max = 25;
S1_ybin_min = -25;
S1_ybin_max = 25;

S2xybinsize = grid_size;%cm
S2_xbin_min = -25;
S2_xbin_max = 25;
S2_ybin_min = -25;
S2_ybin_max = 25;


% Pulse Classification
 Pulse_Class=1;
  
    d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN   
    events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10
          
    s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
    s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);
    
    s1_class=(d.pulse_classification==1 )& s1_area_cut ; %area cut for Kr events
    s2_class=(d.pulse_classification==2) & s2_area_cut ;   
    
    s1_single_cut =logical( (s1_class & d.golden & d.selected_s1_s2).*repmat(sum(s2_class & d.golden & d.selected_s1_s2)==1,events(1),1) ); % was s1_before_s2_cut before using golden
    s2_single_cut =logical( (s2_class & d.golden & d.selected_s1_s2).*repmat(sum(s1_class & d.golden & d.selected_s1_s2)==1,events(1),1) );
    
    
    drift_time = d.z_drift_samples(s2_single_cut)/100;  % us
        
    d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area
    
    s1_phe_both = d.pulse_area_phe(s1_single_cut);
    s1_phe_bottom = d.phe_bottom(s1_single_cut);
    
    s2_phe_both = d.pulse_area_phe(s2_single_cut);
    s2_phe_bottom = d.phe_bottom(s2_single_cut);
    
    s2x = d.x_cm(s2_single_cut);
    s2y = d.y_cm(s2_single_cut);   
    
    d.livetime_sec=sum(d.livetime_end_samples-d.livetime_latch_samples)/1e8;
    evt_cut=logical(sum(s2_single_cut));%Cut for all the events passing the single S1 & S2 cut
    event_number=d.event_number(evt_cut);
    event_timestamp_samples=d.event_timestamp_samples(evt_cut);
    
    time_wait_cut=event_timestamp_samples/1e8/60 > delay_min; %allow x min for Kr mixing

    
    s2radius = (s2x.^2+s2y.^2).^(0.5);
    clean_cut=s2radius<25;

         s1_phe_bottom=s1_phe_bottom(clean_cut);
         s1_phe_both=s1_phe_both(clean_cut);
         s2_phe_bottom=s2_phe_bottom(clean_cut);
         s2_phe_both=s2_phe_both(clean_cut);
         drift_time=drift_time(clean_cut);
         s2x=s2x(clean_cut);
         s2y=s2y(clean_cut);
         s2radius = (s2x.^2+s2y.^2).^(0.5);
         event_timestamp_samples=event_timestamp_samples(clean_cut);
         
         
%%%%%%%%%%%%%%%%%
%Correctins S2 Z and finds remaining "field" dependence
%%%%%%%%%%%%%%%%%

         %%%%%%%% CORRECTING S2 SIGNAL BASED ON CH3T %%%%%%%%%%%%%
         clear hist_s2_bottom Fit_s2_bottom means means_error mean_index mean_ratios s2_phe_bottom_zslice
%          electron_lifetime=e_lifetime_landau;
         s2_phe_bottom_z=s2_phe_bottom.*exp(drift_time./electron_lifetime);
         s2_phe_bottom_xyz=s2_phe_bottom_z.*interp2(old_s2xbins,old_s2ybins,norm_S2_bot,s2x,s2y,'spline');
         s2_phe_both_z=s2_phe_both.*exp(drift_time./electron_lifetime_both);
         s2_phe_both_xyz=s2_phe_both_z.*interp2(old_s2xbins_both,old_s2ybins_both,norm_S2_both,s2x,s2y,'spline');
   
%Calculating remaining S2 Z dependence
clear dT_loc s2z_mean s2z_fit
dT_step=10;
i=1;
for dT_max=10+dT_step:dT_step:320-dT_step;
s2z_fit=fit([0:100:16000]',hist(s2_phe_bottom_xyz(inrange(drift_time,[dT_max-dT_step dT_max])),[0:100:16000])','gauss1');
s2z_mean(i)=s2z_fit.b1;
znorm_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

s2z_normalization_fit=polyfit(znorm_dT_loc,s2z_mean,6);
fit_points=0:2:330;
s2z_fit_line=polyval(s2z_normalization_fit,fit_points);

figure
scatter(znorm_dT_loc,s2z_mean,'.k');
hold on;
plot(fit_points,s2z_fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('CH3T Corrected S2 Peak (phe)'); myfigview(16);
title('Field Normalization from CH3T Corrected Kr S2 Peaks');



% z_top=interp1(znorm_dT_loc,s2z_mean,10,'spline');%Normalize to the center (x=y=0)
% norm_S2_z=z_top./s2z_mean; 
% norm_S2_z(isinf(norm_S2_z))=1;%no infinity! Do not correct events outside r=25cm
% norm_S2_z(isnan(norm_S2_z))=1;
s2_phe_bottom_xyz_z=s2_phe_bottom_xyz.*polyval(s2z_normalization_fit,0)./polyval(s2z_normalization_fit,drift_time);

%%Showing Z correction worked
i=1;
for dT_max=10+dT_step:dT_step:320-dT_step;
s2z_fit=fit([0:100:16000]',hist(s2_phe_bottom_xyz_z(inrange(drift_time,[dT_max-dT_step dT_max])),[0:100:16000])','gauss1');
s2z_mean_z(i)=s2z_fit.b1;
i=i+1;
end

figure
scatter(znorm_dT_loc,s2z_mean_z,'.k')
%%
   
%Calculating remaining S2 Z dependence using S2 both
clear dT_loc s2z_mean s2z_fit
dT_step=10;
i=1;
for dT_max=10+dT_step:dT_step:320-dT_step;
s2z_fit_both=fit([0:200:40000]',hist(s2_phe_both_xyz(inrange(drift_time,[dT_max-dT_step dT_max])),[0:200:40000])','gauss1');
s2z_mean_both(i)=s2z_fit_both.b1;
znorm_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

s2z_normalization_fit_both=polyfit(znorm_dT_loc,s2z_mean_both,6);
fit_points=0:2:330;
s2z_fit_line=polyval(s2z_normalization_fit_both,fit_points);

figure
scatter(znorm_dT_loc,s2z_mean_both,'.k');
hold on;
plot(fit_points,s2z_fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('CH3T Corrected S2Both Peak (phe)'); myfigview(16);
title('Field Normalization from CH3T Corrected Kr S2Both Peaks');



% z_top=interp1(znorm_dT_loc,s2z_mean,10,'spline');%Normalize to the center (x=y=0)
% norm_S2_z=z_top./s2z_mean; 
% norm_S2_z(isinf(norm_S2_z))=1;%no infinity! Do not correct events outside r=25cm
% norm_S2_z(isnan(norm_S2_z))=1;
s2_phe_both_xyz_z=s2_phe_both_xyz.*polyval(s2z_normalization_fit_both,0)./polyval(s2z_normalization_fit_both,drift_time);

%%Showing Z correction worked
i=1;
for dT_max=10+dT_step:dT_step:320-dT_step;
s2z_fit=fit([0:200:40000]',hist(s2_phe_both_xyz_z(inrange(drift_time,[dT_max-dT_step dT_max])),[0:200:40000])','gauss1');
s2z_mean_z(i)=s2z_fit.b1;
i=i+1;
end

figure
scatter(znorm_dT_loc,s2z_mean_z,'.k')


%%
%Corrects S2 XY and finds remaining "field" dependence

clear s2_phe_bottom_xyz_z_xymeans
%Calculating remaining S2 XY dependence
xx_step=2;
yy_step=2;
s2_phe_bottom_xyz_z_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s2_phe_bottom_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 300;
            s2_phe_bottom_xyz_z_xymeans_fit=fit([0:100:16000]',hist(s2_phe_bottom_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:100:16000])','gauss1');
            s2_phe_bottom_xyz_z_xymeans(j,i)=s2_phe_bottom_xyz_z_xymeans_fit.b1;
      else 
            s2_phe_bottom_xyz_z_xymeans(j,i)=0;      
      end
      j=j+1;
    end
    i=i+1;
   end
   
s2xbins=-25+xx_step/2:xx_step:25-xx_step/2;
s2ybins=-25+yy_step/2:yy_step:25-yy_step/2;
xy_center=interp2(s2xbins,s2ybins,s2_phe_bottom_xyz_z_xymeans,0,0,'spline');%Normalize to the center (x=y=0)
norm_S2_xy=xy_center./s2_phe_bottom_xyz_z_xymeans; 
norm_S2_xy(isinf(norm_S2_xy))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_xy(isnan(norm_S2_xy))=1;

s2_phe_bottom_xyz_xyz=s2_phe_bottom_xyz_z.*interp2(s2xbins,s2ybins,norm_S2_xy,s2x,s2y,'spline');
s2_xy_xbins=s2xbins;
s2_xy_ybins=s2ybins;
   
%Showing that XY correction worked
xx_step=2;
yy_step=2;
s2_phe_bottom_xyz_xyz_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s2_phe_bottom_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 300;
            s2_phe_bottom_xyz_xyz_xymeans_fit=fit([0:100:16000]',hist(s2_phe_bottom_xyz_xyz(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:100:16000])','gauss1');
            s2_phe_bottom_xyz_xyz_xymeans(j,i)=s2_phe_bottom_xyz_xyz_xymeans_fit.b1;
      else 
            s2_phe_bottom_xyz_xyz_xymeans(j,i)=0;      
      end
      j=j+1;
    end
    i=i+1;
   end
s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
    
    color_range_max=max(max(s2_phe_bottom_xyz_xyz_xymeans));
    color_range_min=min(min(s2_phe_bottom_xyz_xyz_xymeans(s2_phe_bottom_xyz_xyz_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_bottom_fig = figure;
    contourf(s2x_bins,s2y_bins,s2_phe_bottom_xyz_xyz_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
%%
%Corrects S2 Both XY and finds remaining "field" dependence

clear s2_phe_both_xyz_z_xymeans
%Calculating remaining S2 XY dependence
xx_step=2;
yy_step=2;
s2_phe_both_xyz_z_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s2_phe_both_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 300;
            s2_phe_both_xyz_z_xymeans_fit=fit([0:200:40000]',hist(s2_phe_both_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:200:40000])','gauss1');
            s2_phe_both_xyz_z_xymeans(j,i)=s2_phe_both_xyz_z_xymeans_fit.b1;
      else 
            s2_phe_both_xyz_z_xymeans(j,i)=0;      
      end
      j=j+1;
    end
    i=i+1;
   end
   
s2xbins=-25+xx_step/2:xx_step:25-xx_step/2;
s2ybins=-25+yy_step/2:yy_step:25-yy_step/2;
xy_center=interp2(s2xbins,s2ybins,s2_phe_both_xyz_z_xymeans,0,0,'spline');%Normalize to the center (x=y=0)
norm_S2_xy_both=xy_center./s2_phe_both_xyz_z_xymeans; 
norm_S2_xy_both(isinf(norm_S2_xy_both))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_xy_both(isnan(norm_S2_xy_both))=1;

s2_phe_both_xyz_xyz=s2_phe_both_xyz_z.*interp2(s2xbins,s2ybins,norm_S2_xy_both,s2x,s2y,'spline');
s2_xy_xbins_both=s2xbins;
s2_xy_ybins_both=s2ybins;
   
%Showing that XY correction worked
xx_step=2;
yy_step=2;
s2_phe_both_xyz_xyz_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s2_phe_both_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 300;
            s2_phe_both_xyz_xyz_xymeans_fit=fit([0:200:40000]',hist(s2_phe_both_xyz_xyz(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:200:40000])','gauss1');
            s2_phe_both_xyz_xyz_xymeans(j,i)=s2_phe_both_xyz_xyz_xymeans_fit.b1;
      else 
            s2_phe_both_xyz_xyz_xymeans(j,i)=0;      
      end
      j=j+1;
    end
    i=i+1;
   end
s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
    
    color_range_max=max(max(s2_phe_both_xyz_xyz_xymeans));
    color_range_min=min(min(s2_phe_both_xyz_xyz_xymeans(s2_phe_both_xyz_xyz_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_both_fig = figure;
    contourf(s2x_bins,s2y_bins,s2_phe_both_xyz_xyz_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S2 Mean (both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
    

    
    %% Correcting s1 for z and xy based on CH3T and calculating remaining field dependence
    
             %%%%%%%% CORRECTING S1 SIGNAL BASED ON CH3T %%%%%%%%%%%%%
         clear hist_s2_bottom Fit_s2_bottom means means_error mean_index mean_ratios s2_phe_bottom_zslice
         s1_phe_both_z=s1_phe_both.*polyval(h3_mean_fit_s1z,(det_edge-4)/2)./polyval(h3_mean_fit_s1z,drift_time);
         s1_phe_both_xyz=s1_phe_both_z.*interp2(old_s1xbins,old_s1ybins,norm_S1_both,s2x,s2y,'spline');
         s1_phe_bottom_z=s1_phe_bottom.*polyval(h3_mean_fit_s1z_bottom,(det_edge-4)/2)./polyval(h3_mean_fit_s1z_bottom,drift_time);
         s1_phe_bottom_xyz=s1_phe_bottom_z.*interp2(old_s1xbins_bottom,old_s1ybins_bottom,norm_S1_bottom,s2x,s2y,'spline');
         
         
%%
%%%%%%%
% GET S1 Z Dependence
%%%%%%%%%%%

dT_step=10;
% dT_max=30+dT_step;
i=1;

for dT_max=30+dT_step:dT_step:det_edge-dT_step;
s1xyz_mean_fit=fit([0:3:500]',hist(s1_phe_both_xyz(inrange(drift_time,[dT_max-dT_step dT_max])),[0:3:500])','gauss1'); 
s1xyz_mean(i)=s1xyz_mean_fit.b1;
s1xyz_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

mean_fit_s1z=polyfit(s1xyz_dT_loc,s1xyz_mean,6);
%     sigma_P_s1_both = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';    
fit_points=30:2:det_edge;
fit_line=polyval(mean_fit_s1z,fit_points);

figure
scatter(s1xyz_dT_loc,s1xyz_mean,'.k');
hold on;
plot(fit_points,fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('CH3T Corrected Kr S1 Peak (phe)'); myfigview(16);
title('Kr S1 Z Dependence After CH3T Corrections');

det_edge=320;
s1_phe_both_xyz_z=s1_phe_both_xyz.*polyval(mean_fit_s1z,(det_edge-4)/2)./polyval(mean_fit_s1z,drift_time);

% Showing that the correction worked

dT_step=10;
dT_max=30+dT_step;
i=1;
for dT_max=30+dT_step:dT_step:det_edge-dT_step;
% h3_s1_mean_test(i)=histfitlandau(h3_s1_phe_both_z(inrange(h3_energy,[3 5]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])), 0.5, 1, 100);
s1xyz_mean_fit_test=fit([0:3:500]',hist(s1_phe_both_xyz_z(inrange(drift_time,[dT_max-dT_step dT_max])),[0:3:500])','gauss1'); 
s1xyz_mean_test(i)=s1xyz_mean_fit_test.b1;
s1z_dT_loc_test(i)=dT_max-dT_step/2;
i=i+1;
end

figure
scatter(s1z_dT_loc_test,s1xyz_mean_test,'.k')


%%
%%%%%%%
% GET S1 Z Dependence with S1 bottom
%%%%%%%%%%%

dT_step=10;
dT_max=30+dT_step;
i=1;

for dT_max=30+dT_step:dT_step:det_edge-dT_step;
s1xyz_mean_fit_bottom=fit([0:3:500]',hist(s1_phe_bottom_xyz(inrange(drift_time,[dT_max-dT_step dT_max])),[0:3:500])','gauss1'); 
s1xyz_mean_bottom(i)=s1xyz_mean_fit_bottom.b1;
s1xyz_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

mean_fit_s1z_bottom=polyfit(s1xyz_dT_loc,s1xyz_mean_bottom,6);
%     sigma_P_s1_both = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';    
fit_points=30:2:det_edge;
fit_line=polyval(mean_fit_s1z_bottom,fit_points);

figure
scatter(s1xyz_dT_loc,s1xyz_mean_bottom,'.k');
hold on;
plot(fit_points,fit_line,'-r','LineWidth',2);
xlabel('Drift Time (uSec)');ylabel('CH3T Corrected KrBottom S1 Peak (phe)'); myfigview(16);
title('Kr S1Bottom Z Dependence After CH3T Corrections');

det_edge=320;
s1_phe_bottom_xyz_z=s1_phe_bottom_xyz.*polyval(mean_fit_s1z_bottom,(det_edge-4)/2)./polyval(mean_fit_s1z_bottom,drift_time);

% Showing that the correction worked

dT_step=10;
dT_max=30+dT_step;
i=1;
for dT_max=30+dT_step:dT_step:det_edge-dT_step;
% h3_s1_mean_test(i)=histfitlandau(h3_s1_phe_both_z(inrange(h3_energy,[3 5]) & inrange(h3_drift_time,[dT_max-dT_step dT_max])), 0.5, 1, 100);
s1xyz_mean_fit_test=fit([0:3:500]',hist(s1_phe_bottom_xyz_z(inrange(drift_time,[dT_max-dT_step dT_max])),[0:3:500])','gauss1'); 
s1xyz_mean_test(i)=s1xyz_mean_fit_test.b1;
i=i+1;
end

figure
scatter(s1z_dT_loc_test,s1xyz_mean_test,'.k')

%%
%%%%%%%%%%%%%%%%%%%
% Get remaining S1 XY field dependence
%%%%%%%%%%%%%%%%%%%
clear x2 hist_bin

xx_step=3;
yy_step=3;
s1_both_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s1_phe_both_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 100;
            s1_both_xymeans_fit=fit([0:3:500]',hist(s1_phe_both_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:3:500])','gauss1'); 
            s1_both_xymeans(j,i)=s1_both_xymeans_fit.b1;
           else 
          
        s1_both_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    s1_both_xymeans(isnan(s1_both_xymeans)) = 0;
     
   
%%Plot the correction S1_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
s1_xy_xbins=s2x_bins;
s1_xy_ybins=s2y_bins;
    
    color_range_max=max(max(s1_both_xymeans));
    color_range_min=min(min(s1_both_xymeans(s1_both_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s1_xy_fig = figure;
    contourf(s2x_bins,s2y_bins,s1_both_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S1 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
% 
% %%Plot the 1 sigma of the mean
% 
%     color_range_max=max(max(Sigma_S2_bottom));
%     color_range_min=min(min(Sigma_S2_bottom(Sigma_S2_bottom>0)));
%     vc_step=(color_range_max-color_range_min)/50;
% 
%     vc=color_range_min:vc_step:color_range_max;
% 
%     sigma_s2_xy_bottom_fig = figure;
%     contourf(s2xbins,s2ybins,Sigma_S2_bottom,vc,'LineColor','none');
%     xlabel('x (cm)','fontsize', 18);
%     ylabel('y (cm)','fontsize', 18);
%     title(strcat(file_id_cp, '. 1-Sigma S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
%     caxis([color_range_min color_range_max])
%     g=colorbar;
%     ylabel(g,'Phe','FontSize',16)
%     set(g,'FontSize',18)
%     hold on
%     LUXPlotTopPMTs
%     axis([-25 25 -25 25])
%     myfigview(16);

    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,s1_both_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S1_both=center./s1_both_xymeans; 
norm_S1_both(isinf(norm_S1_both))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_both(isnan(norm_S1_both))=1;
% 
% sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_bottom,0,0,'spline',1);%Normalize to the center (x=y=0)
% sigma_norm_S2_bot=sqrt((sigma_center./mean_S2_bottom_xy).^2+(Sigma_S2_bottom.*center./mean_S2_bottom_xy.^2).^2);
% sigma_norm_S2_bot(isinf(sigma_norm_S2_bot))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
% sigma_norm_S2_bot(isnan(sigma_norm_S2_bot))=1;%no nan. Set Sigma=1, which is 100% uncertainty
% 
s1_phe_both_xyz_xyz=s1_phe_both_xyz_z.*interp2(s1_xy_xbins,s1_xy_ybins,norm_S1_both,s2x,s2y,'spline');


%Prove the correction worked
xx_step=3;
yy_step=3;
s1_both_xymeans_test = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s1_phe_both_xyz_xyz(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 100;
            s1_both_xymeans_fit_test=fit([0:3:500]',hist(s1_phe_both_xyz_xyz(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:3:500])','gauss1'); 
            s1_both_xymeans_test(j,i)=s1_both_xymeans_fit_test.b1;
           else 
          
        s1_both_xymeans_test(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    s1_both_xymeans_test(isnan(s1_both_xymeans_test)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;


    color_range_max=max(max(s1_both_xymeans_test));
    color_range_min=min(min(s1_both_xymeans_test(s1_both_xymeans_test>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s1_xy_test_fig = figure;
    contourf(s2x_bins,s2y_bins,s1_both_xymeans_test,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('Corrected S1 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
    
%%
%%%%%%%%%%%%%%%%%%%
% Get remaining S1 XY field dependence using S1 bot
%%%%%%%%%%%%%%%%%%%
clear x2 hist_bin

xx_step=3;
yy_step=3;
s1_bottom_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s1_phe_bottom_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 100;
            s1_bottom_xymeans_fit=fit([0:3:500]',hist(s1_phe_bottom_xyz_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:3:500])','gauss1'); 
            s1_bottom_xymeans(j,i)=s1_bottom_xymeans_fit.b1;
           else 
          
        s1_bottom_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    s1_bottom_xymeans(isnan(s1_bottom_xymeans)) = 0;
     
   
%%Plot the correction S1_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;
s1_xy_xbins_bottom=s2x_bins;
s1_xy_ybins_bottom=s2y_bins;
    
    color_range_max=max(max(s1_bottom_xymeans));
    color_range_min=min(min(s1_bottom_xymeans(s1_bottom_xymeans>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s1_xy_fig = figure;
    contourf(s2x_bins,s2y_bins,s1_bottom_xymeans,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('S1 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
% 
% %%Plot the 1 sigma of the mean
% 
%     color_range_max=max(max(Sigma_S2_bottom));
%     color_range_min=min(min(Sigma_S2_bottom(Sigma_S2_bottom>0)));
%     vc_step=(color_range_max-color_range_min)/50;
% 
%     vc=color_range_min:vc_step:color_range_max;
% 
%     sigma_s2_xy_bottom_fig = figure;
%     contourf(s2xbins,s2ybins,Sigma_S2_bottom,vc,'LineColor','none');
%     xlabel('x (cm)','fontsize', 18);
%     ylabel('y (cm)','fontsize', 18);
%     title(strcat(file_id_cp, '. 1-Sigma S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
%     caxis([color_range_min color_range_max])
%     g=colorbar;
%     ylabel(g,'Phe','FontSize',16)
%     set(g,'FontSize',18)
%     hold on
%     LUXPlotTopPMTs
%     axis([-25 25 -25 25])
%     myfigview(16);

    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,s1_bottom_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S1_bottom=center./s1_bottom_xymeans; 
norm_S1_bottom(isinf(norm_S1_bottom))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_bottom(isnan(norm_S1_bottom))=1;
% 
% sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_bottom,0,0,'spline',1);%Normalize to the center (x=y=0)
% sigma_norm_S2_bot=sqrt((sigma_center./mean_S2_bottom_xy).^2+(Sigma_S2_bottom.*center./mean_S2_bottom_xy.^2).^2);
% sigma_norm_S2_bot(isinf(sigma_norm_S2_bot))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
% sigma_norm_S2_bot(isnan(sigma_norm_S2_bot))=1;%no nan. Set Sigma=1, which is 100% uncertainty
% 
s1_phe_bottom_xyz_xyz=s1_phe_bottom_xyz_z.*interp2(s1_xy_xbins,s1_xy_ybins,norm_S1_bottom,s2x,s2y,'spline');


%Prove the correction worked
xx_step=3;
yy_step=3;
s1_bottom_xymeans_test = zeros(floor(50/yy_step),floor(50/xx_step));
i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s1_phe_bottom_xyz_xyz(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 100;
            s1_bottom_xymeans_fit_test=fit([0:3:500]',hist(s1_phe_bottom_xyz_xyz(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:3:500])','gauss1'); 
            s1_bottom_xymeans_test(j,i)=s1_bottom_xymeans_fit_test.b1;
           else 
          
        s1_bottom_xymeans_test(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
    s1_bottom_xymeans_test(isnan(s1_bottom_xymeans_test)) = 0;
     
   
%%Plot the correction S2_XY bottom %%%%%%%%%%%%%%%

s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;


    color_range_max=max(max(s1_bottom_xymeans_test));
    color_range_min=min(min(s1_bottom_xymeans_test(s1_bottom_xymeans_test>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s1_xy_test_fig = figure;
    contourf(s2x_bins,s2y_bins,s1_bottom_xymeans_test,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat('Corrected S1 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
    
%%

% 
%     
%Comparing spectra
figure
step(s2_phe_bottom_xyz_z,[0:100:16000],'b')
hold on;
step(s2_phe_bottom_xyz_xyz,[0:100:16000],'g')
step(s2_phe_bottom_xyz,[0:100:16000],'r')
step(s2_phe_bottom_z,[0:100:16000],'m')
step(s2_phe_bottom,[0:100:16000],'k')

figure
step(s1_phe_both_xyz_z,[0:2:500],'b')
hold on;
step(s1_phe_both_xyz_xyz,[0:2:500],'g')
step(s1_phe_both_xyz,[0:2:500],'r')
step(s1_phe_both_z,[0:2:500],'m')
step(s1_phe_both,[0:2:500],'k')

save('CH3T_field_normalization_DP21','s2z_normalization_fit','s2z_normalization_fit_both','norm_S1_both','norm_S1_bottom','norm_S2_xy','norm_S2_xy_both','s2_xy_xbins','s2_xy_xbins_both','s2_xy_ybins',...
    's2_xy_ybins_both','mean_fit_s1z','mean_fit_s1z_bottom','s1_xy_xbins','s1_xy_xbins_bottom','s1_xy_ybins','s1_xy_ybins_bottom');
