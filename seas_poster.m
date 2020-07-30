clc
clear
addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));
surf=load('sinfit_temp_data_plotting.mat');
core=load('sinfit_temp_data_plotting_core.mat');
deep=load('sinfit_temp_data_plotting_deep.mat');
pth=load('pathXY.mat');
load pathXY_DEEP

%% BC Speed

surf.S_idx=strcmp(surf.rel_station_names,'SN-3');
surf.N_idx=strcmp(surf.rel_station_names,'IOS_44');
EBCtime=surf.cold_day(surf.S_idx)-surf.cold_day(surf.N_idx); %days
EBCdist=sum(gsw_distance(pth.Xlon,pth.Xlat)/1000);%KM
surf.EBCspeed=(EBCdist*100000)/(EBCtime*60*60*24); %cm s-1

fprintf('SURFACE: EBC distance= %2.0f km \n EBC time= %2.0f days \n EBC speed= %2.2f cm s-1 \n\n',...
    EBCdist,EBCtime,surf.EBCspeed);

core.S_idx=strcmp(core.rel_station_names,'SN-3');
core.N_idx=strcmp(core.rel_station_names,'IOS_44');
EBCtime=core.cold_day(core.S_idx)-core.cold_day(core.N_idx); %days
EBCdist=sum(gsw_distance(pth.Xlon,pth.Xlat)/1000);%KM
core.EBCspeed=(EBCdist*100000)/(EBCtime*60*60*24); %cm s-1

fprintf('CORE: EBC distance= %2.0f km \n EBC time= %2.0f days \n EBC speed= %2.2f cm s-1 \n\n',...
    EBCdist,EBCtime,core.EBCspeed);

deep.S_idx=strcmp(deep.rel_station_names,'SN-5');
deep.N_idx=strcmp(deep.rel_station_names,'IOS_45');
EBCtime=deep.cold_day(deep.S_idx)-deep.cold_day(deep.N_idx); %days
EBCdist=sum(gsw_distance(pth_DEEP.Xlon,pth_DEEP.Xlat)/1000);%KM
deep.EBCspeed=(EBCdist*100000)/(EBCtime*60*60*24); %cm s-1

fprintf('CORE: EBC distance= %2.0f km \n EBC time= %2.0f days \n EBC speed= %2.2f cm s-1 \n\n',...
    EBCdist,EBCtime,deep.EBCspeed);

%% Plotting - Phase shift


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface

% Set up XY projection
lat_limV=[min(surf.grid_lat(:)) max(surf.grid_lat(:))];
lon_limV=[min(surf.grid_lon(:)) max(surf.grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(surf.station_lon,surf.station_lat);
[Xgrid,Ygrid]=m_ll2xy(surf.grid_lon,surf.grid_lat);

dday = variogram([X' Y'],surf.cold_day','maxdist',49000,'plotit',false);

[a,c,~,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,X',Y',surf.cold_day',Xgrid,Ygrid);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([surf.station_lon';-125.25],...
    [surf.station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,surf.grid_lon,surf.grid_lat);
zi(~inside_idx)=NaN;
zi(surf.Zint>-10)=NaN; % 10m SHALLOW WATER MASK
zi_age=zi-40;

lat_lim=[50.42283746 48.252093478];
lon_lim=[-125.82436392 -122.252345235];

% figure('units','centimeters','outerposition',[0 0 19 23],'color','w'); % paper
figure('units','centimeters','outerposition',[0 0 17 22],'color','w');% poster
ax1=axes('Position',[0.035 0.05 0.3 0.8]);
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
hold on

[CS,CH]=m_contourf(surf.grid_lon,surf.grid_lat,zi_age,[0:10:20 30:15:130],...
        'linestyle','none');
colormap(ax1,m_colmap('jet'));
m_gshhs_f('patch',rgb('light grey'),'edgecolor','k');
s=m_scatter(surf.station_lon,surf.station_lat,5,'x','k');
m_scatter(surf.station_lon([find(surf.N_idx) find(surf.S_idx)]),...
    surf.station_lat([find(surf.N_idx) find(surf.S_idx)]),7,...
    'markerfacecolor','none','markeredgecolor','w');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',7);

text(0.9,0.025,'a)','units','normalized','fontsize',8,'Fontweight','bold');
m_plot(pth.Xlon,pth.Xlat,'w:','linewidth',2);
m_text(-123.0690,49.0852,sprintf('%2.1f cm s^{-1}',surf.EBCspeed),...
    'fontsize',6,'color','w','fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Core

% Set up XY projection
lat_limV=[min(core.grid_lat(:)) max(core.grid_lat(:))];
lon_limV=[min(core.grid_lon(:)) max(core.grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(core.station_lon,core.station_lat);
[Xgrid,Ygrid]=m_ll2xy(core.grid_lon,core.grid_lat);

dday = variogram([X' Y'],core.cold_day','maxdist',49000,'plotit',false);

[a,c,~,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,X',Y',core.cold_day',Xgrid,Ygrid);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([core.station_lon';-125.25],...
    [core.station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,core.grid_lon,core.grid_lat);
zi(~inside_idx)=NaN;
zi(core.Zint>-10)=NaN; % 10m SHALLOW WATER MASK
zi_age=zi-40;

ax2=axes('Position',[0.36 0.05 0.3 0.8]);
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
hold on
RGB=rgb('light grey');

[CS,CH]=m_contourf(core.grid_lon,core.grid_lat,zi_age,[0:10:20 30:15:130],...
        'linestyle','none');
colormap(ax2,m_colmap('jet'));
m_gshhs_f('patch',rgb('light grey'),'edgecolor','k');
s=m_scatter(core.station_lon,core.station_lat,5,'x','k');
m_scatter(core.station_lon([find(core.N_idx) find(core.S_idx)]),...
    core.station_lat([find(core.N_idx) find(core.S_idx)]),7,...
    'markerfacecolor','none','markeredgecolor','w');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xticklabels',[],'yticklabels',[]);
m_northarrow(-125.8184,50.0591,0.25,'type',2);
% [ax,h]=m_contfbar(ax2,[-.9 1.9],1,CS,CH,'endpiece','no','axfrac',.03,'fontsize',8);% paper
[ax,h]=m_contfbar(ax2,[-.9 1.9],1.05,CS,CH,'endpiece','no','axfrac',.03,'fontsize',8,...
    'xtick',[0:10:20 30:15:130]);
ylabel(ax,{'Age (days)'},'fontsize',8,'color','k','fontweight','bold','rotation',0);
text(0.9,0.025,'b)','units','normalized','fontsize',8,'Fontweight','bold');
m_plot(pth.Xlon,pth.Xlat,'w:','linewidth',2);
m_text(-123.0690,49.0852,sprintf('%2.1f cm s^{-1}',core.EBCspeed),...
    'fontsize',6,'color','w','fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deep

% Set up XY projection
lat_limV=[min(deep.grid_lat(:)) max(deep.grid_lat(:))];
lon_limV=[min(deep.grid_lon(:)) max(deep.grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(deep.station_lon,deep.station_lat);
[Xgrid,Ygrid]=m_ll2xy(deep.grid_lon,deep.grid_lat);

dday = variogram([X' Y'],deep.cold_day','maxdist',49000,'plotit',false);

[a,c,~,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,X',Y',deep.cold_day',Xgrid,Ygrid);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([deep.station_lon';-125.25],...
    [deep.station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,deep.grid_lon,deep.grid_lat);
zi(~inside_idx)=NaN;
zi(deep.Zint>-10)=NaN; % 10m SHALLOW WATER MASK
zi_age=zi-40;

ax3=axes('Position',[0.685 0.05 0.3 0.8]);
% ax3=subplot(2,2,1);
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
hold on
RGB=rgb('light grey');

[CS,CH]=m_contourf(deep.grid_lon,deep.grid_lat,zi_age,[0:10:20 30:15:130],...
        'linestyle','none');
colormap(ax3,m_colmap('jet'));
m_gshhs_f('patch',rgb('light grey'),'edgecolor','k');
s=m_scatter(deep.station_lon,deep.station_lat,5,'x','k');
m_scatter(deep.station_lon([find(deep.N_idx) find(deep.S_idx)]),...
    deep.station_lat([find(deep.N_idx) find(deep.S_idx)]),7,...
    'markerfacecolor','none','markeredgecolor','w');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xticklabels',[],'yticklabels',[]);
text(0.9,0.025,'c)','units','normalized','fontsize',8,'Fontweight','bold');
m_plot(pth_DEEP.Xlon,pth_DEEP.Xlat,'w:','linewidth',2);
m_text(-123.0690,49.0852,sprintf('%2.1f cm s^{-1}',deep.EBCspeed),...
    'fontsize',6,'color','w','fontweight','bold');

%%
set(gcf,'color','w');
export_fig /ocean/sstevens/IW_project/figures/paper/seas_CD_poster_V2.png -png -m3
print /ocean/sstevens/IW_project/figures/paper/vec/seas_CD_poster_V2.pdf -dpdf -painters
%% Plotting - Amplitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface

% Set up XY projection
lat_limV=[min(surf.grid_lat(:)) max(surf.grid_lat(:))];
lon_limV=[min(surf.grid_lon(:)) max(surf.grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(surf.station_lon,surf.station_lat);
[Xgrid,Ygrid]=m_ll2xy(surf.grid_lon,surf.grid_lat);

damp = variogram([X' Y'],surf.amp','maxdist',40000,'plotit',false);

[a,c,~,vstruct] = variogramfit(damp.distance,damp.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,X',Y',surf.amp',Xgrid,Ygrid);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([surf.station_lon';-125.25],...
    [surf.station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,surf.grid_lon,surf.grid_lat);
zi(~inside_idx)=NaN;
zi(surf.Zint>-10)=NaN; % 10m SHALLOW WATER MASK

lat_lim=[50.42283746 48.252093478];
lon_lim=[-125.82436392 -122.252345235];

% figure('units','centimeters','outerposition',[0 0 19 23],'color','w');% paper
figure('units','centimeters','outerposition',[0 0 17 22],'color','w');% poster
ax1=axes('Position',[0.035 0.05 0.3 0.8]);
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
hold on

[CS,CH]=m_contourf(surf.grid_lon,surf.grid_lat,zi,[0:0.2:0.8 1:.25:1.5],...
    'linestyle','none');
colormap(ax1,flipud(m_colmap('jet')));
m_gshhs_f('patch',rgb('light grey'),'edgecolor','k');
s=m_scatter(surf.station_lon,surf.station_lat,5,'x','k');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',7);
caxis([0 1.5]);
text(0.9,0.025,'a)','units','normalized','fontsize',8,'Fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Core

% Set up XY projection
lat_limV=[min(core.grid_lat(:)) max(core.grid_lat(:))];
lon_limV=[min(core.grid_lon(:)) max(core.grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(core.station_lon,core.station_lat);
[Xgrid,Ygrid]=m_ll2xy(core.grid_lon,core.grid_lat);

damp = variogram([X' Y'],core.amp','maxdist',40000,'plotit',false);

[a,c,~,vstruct] = variogramfit(damp.distance,damp.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,X',Y',core.amp',Xgrid,Ygrid);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([core.station_lon';-125.25],...
    [core.station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,core.grid_lon,core.grid_lat);
zi(~inside_idx)=NaN;
zi(core.Zint>-10)=NaN; % 10m SHALLOW WATER MASK

ax2=axes('Position',[0.36 0.05 0.3 0.8]);
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
hold on
RGB=rgb('light grey');

[CS,CH]=m_contourf(core.grid_lon,core.grid_lat,zi,[0:0.2:0.8 1:.25:1.5],...
        'linestyle','none');
colormap(ax2,flipud(m_colmap('jet')));
caxis([0 1.5]);
m_gshhs_f('patch',rgb('light grey'),'edgecolor','k');
s=m_scatter(core.station_lon,core.station_lat,5,'x','k');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xticklabels',[],'yticklabels',[]);
m_northarrow(-125.8184,50.0591,0.25,'type',2);
% [ax,h]=m_contfbar(ax2,[-.9 1.9],1,CS,CH,'endpiece','no','axfrac',.03,'fontsize',8);%paper
[ax,h]=m_contfbar(ax2,[-.9 1.9],1.05,CS,CH,'endpiece','no','axfrac',.03,'fontsize',8,...
    'xtick',[0:0.2:0.8 1:.25:1.5]); %poster
ylabel(ax,{'Amplitude';'(^\circC)'},'fontsize',8,'color','k','fontweight','bold','rotation',0);
text(0.9,0.025,'b)','units','normalized','fontsize',8,'Fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deep

% Set up XY projection
lat_limV=[min(deep.grid_lat(:)) max(deep.grid_lat(:))];
lon_limV=[min(deep.grid_lon(:)) max(deep.grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(deep.station_lon,deep.station_lat);
[Xgrid,Ygrid]=m_ll2xy(deep.grid_lon,deep.grid_lat);

damp = variogram([X' Y'],deep.amp','maxdist',40000,'plotit',false);

[a,c,~,vstruct] = variogramfit(damp.distance,damp.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,X',Y',deep.amp',Xgrid,Ygrid);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([deep.station_lon';-125.25],...
    [deep.station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,deep.grid_lon,deep.grid_lat);
zi(~inside_idx)=NaN;
zi(deep.Zint>-10)=NaN; % 10m SHALLOW WATER MASK

ax3=axes('Position',[0.685 0.05 0.3 0.8]);
% ax3=subplot(2,2,1);
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
hold on
RGB=rgb('light grey');

[CS,CH]=m_contourf(deep.grid_lon,deep.grid_lat,zi,[0:0.2:0.8 1:.25:1.5],...
        'linestyle','none');
colormap(ax3,flipud(m_colmap('jet')));
caxis([0 1.5]);
m_gshhs_f('patch',rgb('light grey'),'edgecolor','k');
s=m_scatter(deep.station_lon,deep.station_lat,5,'x','k');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xticklabels',[],'yticklabels',[]);
text(0.9,0.025,'c)','units','normalized','fontsize',8,'Fontweight','bold');
%%
export_fig /ocean/sstevens/IW_project/figures/paper/seas_AMP_poster_V2.png -png -m3
print /ocean/sstevens/IW_project/figures/paper/vec/seas_AMP_poster_V2.pdf -dpdf -painters
