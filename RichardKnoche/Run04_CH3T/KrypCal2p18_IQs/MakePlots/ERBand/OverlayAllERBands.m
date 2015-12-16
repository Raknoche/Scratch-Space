x=[0:0.1:50];

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\Sep2015_2p18_ERFit.mat')
y_Sep2015_lower=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
y_Sep2015_middle=ER_mean_power_fit.a.*x.^ER_mean_power_fit.b;
y_Sep2015_upper=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;


load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\Sep2014_2p18_ERFit.mat')
y_Sep2014_lower=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
y_Sep2014_middle=ER_mean_power_fit.a.*x.^ER_mean_power_fit.b;
y_Sep2014_upper=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;


load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\Nov2014_2p18_ERFit.mat')
y_Nov2014_lower=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
y_Nov2014_middle=ER_mean_power_fit.a.*x.^ER_mean_power_fit.b;
y_Nov2014_upper=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;


load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\Feb2015_2p18_ERFit.mat')
y_Feb2015_lower=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
y_Feb2015_middle=ER_mean_power_fit.a.*x.^ER_mean_power_fit.b;
y_Feb2015_upper=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;

figure
hold on;
plot(x,y_Sep2015_middle,'-k','LineWidth',2,'Color',[0 0 0]);
plot(x,y_Feb2015_middle,'-k','LineWidth',2,'Color',[1 0 0]);
plot(x,y_Nov2014_middle,'-k','LineWidth',2,'Color',[0 1 0]);
plot(x,y_Sep2014_middle,'-k','LineWidth',2,'Color',[0 0 1]);
xlabel('Corrected S1 (phe)');
ylabel('Corrected log10(S2/S1)');
legend('Sep 2015','Feb 2015','Nov 2014','Sep 2014');
myfigview(16);
plot(x,y_Sep2015_lower,'--k','LineWidth',2,'Color',[0.6 0.6 0.6]);
plot(x,y_Feb2015_lower,'--k','LineWidth',2,'Color',[1 0.6 0.6]);
plot(x,y_Nov2014_lower,'--k','LineWidth',2,'Color',[0.6 1 0.6]);
plot(x,y_Sep2014_lower,'--k','LineWidth',2,'Color',[0.6 0.6 1]);
plot(x,y_Sep2015_upper,'--k','LineWidth',2,'Color',[0.6 0.6 0.6]);
plot(x,y_Feb2015_upper,'--k','LineWidth',2,'Color',[1 0.6 0.6]);
plot(x,y_Nov2014_upper,'--k','LineWidth',2,'Color',[0.6 1 0.6]);
plot(x,y_Sep2014_upper,'--k','LineWidth',2,'Color',[0.6 0.6 1]);
ylim([1.8 3])