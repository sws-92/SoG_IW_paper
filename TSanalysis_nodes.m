%% Time-series analysis of IW currents
addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));
%% Load ONC nodes
clear
central=load('/ocean/kstankov/ADCP/central/ADCP_central.mat');
DDL=load('/ocean/kstankov/ADCP/ddl/ADCP_ddl.mat');
east=load('/ocean/kstankov/ADCP/east/ADCP_east.mat');
coastang=-30;
%% Take average and rotate
idx=central.chartdepth>50 & central.chartdepth<150;
central.mnV=nanmean(central.vtrue(idx,:),1)/100;
central.mnU=nanmean(central.utrue(idx,:),1)/100;
central.lat=49.0404366;
central.lon=-123.425795;
uv_i=central.mnU+1i*central.mnV;
central.across=imag(uv_i.*exp(-1i*coastang*pi/180));
central.along=real(uv_i.*exp(-1i*coastang*pi/180));

idx=DDL.chartdepth>50 & DDL.chartdepth<150;
DDL.mnV=nanmean(DDL.vtrue(idx,:),1)/100;
DDL.mnU=nanmean(DDL.utrue(idx,:),1)/100;
DDL.lat=49.085095;
DDL.lon=-123.330155;
uv_i=DDL.mnU+1i*DDL.mnV;
DDL.across=imag(uv_i.*exp(-1i*coastang*pi/180));
DDL.along=real(uv_i.*exp(-1i*coastang*pi/180));

idx=east.chartdepth>50 & east.chartdepth<150;
east.mnV=nanmean(east.vtrue(idx,:),1)/100;
east.mnU=nanmean(east.utrue(idx,:),1)/100;
east.lat=49.042835;
east.lon=-123.317265;
uv_i=east.mnU+1i*east.mnV;
east.across=imag(uv_i.*exp(-1i*coastang*pi/180));
east.along=real(uv_i.*exp(-1i*coastang*pi/180));

%% Create yearly climatology

% Central
tmp_day=day(datetime(datestr(central.mtime)),'dayofyear');
central.along_clim=NaN(365,1);
for i=1:365
    central.along_clim(i)=nanmean(central.along(tmp_day==i));
end

% East
tmp_day=day(datetime(datestr(east.mtime)),'dayofyear');
east.along_clim=NaN(365,1);
for i=1:365
    east.along_clim(i)=nanmean(east.along(tmp_day==i));
end

% DDL
tmp_day=day(datetime(datestr(DDL.mtime)),'dayofyear');
DDL.along_clim=NaN(365,1);
for i=1:365
    DDL.along_clim(i)=nanmean(DDL.along(tmp_day==i));
end

%% Look at AC functions and spectrograms
figure
subplot(2,3,1)
tmp=DDL.along_clim-nanmean(DDL.along_clim);
autocorr([tmp;tmp],'NumLags',366,'NumSTD',2)
ylim([-1 1]);
title('DDL');
N = length([tmp;tmp;tmp]);
PER = abs(fft([tmp;tmp;tmp])).^2/N;
Freq = [0:N-1]/N;
subplot(2,3,4)
plot(Freq(1:N/2),PER(1:N/2));
axis tight; grid on
x=[0.00818918108189182;0.0191080891910809;0.0263873612638736;0.0609639036096391;0.0791620837916208;0.0973602639736027;0.232026797320268];
y=[0.0244602418607606;0.0120006012346359;0.00924636488570306;0.00714789909603996;0.00557404975379262;0.00347558396412953;0.00150827228632036];
T=round(1./x);
for i=1:length(T)
    txt = ['\leftarrow ' num2str(T(i)) ' days'];
    text(x(i),y(i),txt)
end

subplot(2,3,2)
tmp=east.along_clim-nanmean(east.along_clim);
autocorr([tmp;tmp],'NumLags',366,'NumSTD',2)
ylim([-1 1]);
title('East');
N = length([tmp;tmp;tmp]);
PER = abs(fft([tmp;tmp;tmp])).^2/N;
Freq = [0:N-1]/N;
subplot(2,3,5)
plot(Freq(1:N/2),PER(1:N/2));
axis tight; grid on
x=[0.0228310502283104;0.0465753424657535;0.0739726027397260;0.106849315068493;0.126940639269406;0.150684931506849];
y=[0.00521305191426108;0.00313685807827831;0.00503251331982780;0.00670249531833568;0.00200849186307029;0.00119606818812051];
T=round(1./x);
for i=1:length(T)
    txt = ['\leftarrow ' num2str(T(i)) ' days'];
    text(x(i),y(i),txt)
end

subplot(2,3,3)
tmp=central.along_clim-nanmean(central.along_clim);
autocorr([tmp;tmp],'NumLags',366,'NumSTD',2)
title('Central');
ylim([-1 1]);
N = length([tmp;tmp;tmp]);
PER = abs(fft([tmp;tmp;tmp])).^2/N;
Freq = [0:N-1]/N;
subplot(2,3,6)
plot(Freq(1:N/2),PER(1:N/2));
axis tight; grid on
x=[0.0118721461187217;0.0356164383561646;0.0520547945205481;0.0684931506849316;0.0831050228310504;0.121461187214612;0.172602739726027];
y=[0.00757880340530981;0.00258962248018423;0.00363497467402006;0.00268465449780567;0.00344491063877718;0.00154427028634839;0.00135420625110551];
T=round(1./x);
for i=1:length(T)
    txt = ['\leftarrow ' num2str(T(i)) ' days'];
    text(x(i),y(i),txt)
end

%% Find longest uniterrupted sections of data and run FFT on those
npos=find(isnan(DDL.along));
[~,idx]=max(diff(npos));
cut=DDL.along(npos(idx)+1:npos(idx+1)-1);

figure
subplot(3,3,1)
[ac,nlags]=autocorr(cut,'NumLags',2*24*140,'NumSTD',2);
plot(nlags/48,ac);
ylim([-1 1]);
title('DDL');
ylabel('ACF');

N = length(cut);
PER = abs(fft(cut)).^2/N;
Freq = [0:N-1]/N;
subplot(3,3,4)
plot(Freq(1:N/2),PER(1:N/2));
axis tight; grid on; xlim([0 0.1]);
ylabel('Spectral coeff.');
x=[0.0206140350877193;0.0407894736842105;0.0609649122807018;0.0811403508771930];
y=[0.0563824101076666;0.323664589650271;0.0313247057755474;0.0146195695541346];
T=round(0.5./x);
for i=1:length(T)
    txt = ['\leftarrow ' num2str(T(i)) ' days'];
    text(x(i),y(i),txt)
end

subplot(3,3,7)
loglog(Freq(1:N/2),PER(1:N/2))
ylabel('log(spectral coeff.)');
xt=get(gca,'xtick');

npos=find(isnan(east.along));
[~,idx]=max(diff(npos));
cut=east.along(npos(idx)+1:npos(idx+1)-1);

subplot(3,3,2)
[ac,nlags]=autocorr(cut,'NumLags',2*24*140,'NumSTD',2);
plot(nlags/48,ac);
ylim([-1 1]);
xlabel('Lag (days)')
title('East');

N = length(cut);
PER = abs(fft(cut)).^2/N;
Freq = [0:N-1]/N;
subplot(3,3,5)
plot(Freq(1:N/2),PER(1:N/2));
axis tight; grid on;xlim([0 0.1]);
xlabel('f (cps)')
% x=[0.000438596491228072;0.0206140350877193;0.0407894736842105;0.0609649122807018;0.0811403508771930];
% y=[0.181670931768263;0.0563824101076666;0.323664589650271;0.0313247057755474;0.0146195695541346];
% T=round(1./x);
% for i=1:length(T)
%     txt = ['\leftarrow ' num2str(T(i)) ' days'];
%     text(x(i),y(i),txt)
% end

subplot(3,3,8)
loglog(Freq(1:N/2),PER(1:N/2))
xlabel('log(f)')
set(gca,'xtick',xt);
for i=1:2
    txt = ['\leftarrow ' num2str(T(i)) ' days'];
    text(x(i),y(i),txt)
end

npos=find(isnan(central.along));
[~,idx]=max(diff(npos));
cut=central.along(npos(idx)+1:npos(idx+1)-1);

subplot(3,3,3)
[ac,nlags]=autocorr(cut,'NumLags',2*24*140,'NumSTD',2);
plot(nlags/48,ac);
ylim([-1 1]);
title('Central');
N = length(cut);
PER = abs(fft(cut)).^2/N;
Freq = [0:N-1]/N;
subplot(3,3,6)
plot(Freq(1:N/2),PER(1:N/2));
axis tight; grid on;xlim([0 0.1]);
% x=[0.000438596491228072;0.0206140350877193;0.0407894736842105;0.0609649122807018;0.0811403508771930];
% y=[0.181670931768263;0.0563824101076666;0.323664589650271;0.0313247057755474;0.0146195695541346];
% T=round(1./x);
% for i=1:length(T)
%     txt = ['\leftarrow ' num2str(T(i)) ' days'];
%     text(x(i),y(i),txt)
% end
subplot(3,3,9)
loglog(Freq(1:N/2),PER(1:N/2))
set(gca,'xtick',xt);

%%
set(gcf,'color','w');
export_fig /ocean/sstevens/IW_project/figures/paper/nodeTS_analysis.pdf