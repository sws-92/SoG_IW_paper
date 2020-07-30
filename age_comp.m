%% AGE comparison between methods
% addpath(genpath('/ocean/sstevens/'));
% addpath(genpath('/ocean/rich/home/matlab/m_map/'));

%% Load
clear
clear global
load NE_path.mat
DP_path=load('DP_path.mat');
Eu=load('eulerianM.mat');
load PT_AS_inlets.mat

%% Eulerian
fid=fopen('/ocean/sstevens/IW_project/code/paper/age_table.txt','w');

% NE path distance
NEdist=gsw_distance([ones(length(NE_path.points(1:20:end)),1)*...
    Eu.mod_lon(NE_path.points(1)) Eu.mod_lon(NE_path.points(1:20:end))],...
    [ones(length(NE_path.points(1:20:end)),1)*...
    Eu.mod_lat(NE_path.points(1)) Eu.mod_lat(NE_path.points(1:20:end))]);

% NE path along strait speeds
tmp=Eu.mean_v; tmp(tmp<0)=NaN;
[l,ll]=ind2sub(size(Eu.mod_lon),NE_path.points);

[llu,lli]=unique(ll);
linspac=diff(llu);

mnV=[];
spac=[];
NElon_list=[]; NElat_list=[];
for i=1:160%length(llu)-1
    
    idx=lli(i):lli(i+1);
    
    if length(idx)==2
        idx=idx(1);
    else
        idx(end)=[];
    end
    
    if linspac(i)==1
        mnV(i)=nanmean(tmp(l(idx),llu(i)));
        spac(i)=500;
        NElon_list=[NElon_list;Eu.mod_lon(l(idx),llu(i))];
        NElat_list=[NElat_list;Eu.mod_lat(l(idx),llu(i))];
    else
        mnV(i)=nanmean(tmp(l(idx),llu(i:i+1)),'all');
        spac(i)=1000;
        lats=Eu.mod_lat(l(idx),llu(i:i+1));
        lons=Eu.mod_lon(l(idx),llu(i:i+1));
        NElon_list=[NElon_list;lons(:)];
        NElat_list=[NElat_list;lats(:)];
    end
    
end
NE_path.lon_list=NElon_list;
NE_path.lat_list=NElat_list;
save /ocean/sstevens/IW_project/data/NE_path.mat NE_path

% Isolate boundary current
% mnV(212:end)=[];spac(212:end)=[];
NEtt=nansum(spac./mnV)/(60*60*24); % days
NE_path.dist=sum(spac)/1000;

fprintf(fid,'Eulerian Southern inflow\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',sum(spac)/1000);
fprintf(fid,'Speed=%4.2f cm/s\r\n',nanmean(mnV)*100);
fprintf(fid,'Transit time=%4.0f days\r\n\r\n',NEtt);

% DP path distance
DPdist=gsw_distance([ones(length(DP_path.points(1:20:end)),1)*...
    Eu.mod_lon(DP_path.points(1)) Eu.mod_lon(DP_path.points(1:20:end))],...
    [ones(length(DP_path.points(1:20:end)),1)*...
    Eu.mod_lat(DP_path.points(1)) Eu.mod_lat(DP_path.points(1:20:end))]);

% DP path along strait speeds
tmp=Eu.mean_v;tmp(tmp>0)=NaN;
[l,ll]=ind2sub(size(Eu.mod_lon),DP_path.points);

[llu,lli]=unique(ll);
linspac=diff(llu);
mnV=[];
spac=[];
for i=1:length(llu)-1
    
    idx=lli(i:i+1);
    
    if length(idx)==2
        idx=idx(1);
    else
        idx(end)=[];
    end
    
    if linspac(i)==1
        
        mnV(i)=nanmean(tmp(l(idx),llu(i)));
        spac(i)=500;
        
    else
        mnV(i)=nanmean(tmp(l(idx),llu(i:i+1)));
        spac(i)=1000;
    end
    
end

DPtt=nansum(spac./mnV*-1)/(60*60*24); % days
DP_path.dist=sum(spac(~isnan(spac)))/1000;

fprintf(fid,'Eulerian Northern inflow\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',sum(spac(~isnan(spac)))/1000);
fprintf(fid,'Speed=%4.2f cm/s\r\n',nanmean(mnV)*100);
fprintf(fid,'Transit time=%4.0f days\r\n\r\n',DPtt);

DPlon_list=Eu.mod_lon(DP_path.points);
DPlat_list=Eu.mod_lat(DP_path.points);
clearvars -except NE_path DP_path fid Eu NElon_list NElat_list...
    DPlon_list DPlat_list

%% Eulerian whole strait
load model_IW_outputUV.mat
load PT_AS_inlets.mat
%%
station_bound=alphaShape(AS.lon,AS.lat);
station_bound.Alpha=station_bound.Alpha*2;
qlon=double(model_uv.lon);qlon(isnan(qlon))=0;
qlat=double(model_uv.lat);qlat(isnan(qlat))=0;
inside_idx=inShape(station_bound,qlon,qlat);

straitV=Eu.mean_v;straitV(~inside_idx)=NaN;

generalV=nanmean(straitV(straitV>0));% ms-1
generalTS=250e3/(generalV*60*60*24); % days

maxV=NaN(size(straitV,2),1);
mnV=NaN(size(straitV,2),1);
for i=1:size(straitV,2)
    tmpV=straitV(:,i);
    maxV(i)=max(tmpV);
    mnV(i)=nanmean(tmpV(tmpV>0)); 
end
maxV(maxV<=0)=NaN;
maxTS=nansum(ones(length(maxV),1)*500./maxV)/(60*60*24); % days
mnTS=nansum(ones(length(mnV),1)*500./mnV)/(60*60*24); % days

fprintf(fid,'Eulerian whole strait\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',sum(~isnan(mnV))/2);
fprintf(fid,'Speed=%4.2f cm/s\r\n',nanmean(mnV)*100);
fprintf(fid,'Transit time=%4.0f days\r\n\r\n',mnTS);

%% Lagrangian
clearvars -except NE_path DP_path fid Eu NElon_list NElat_list...
    DPlon_list DPlat_list
load lagrangianM
%%
% Find first and last point on path
idx=dsearchn([BP.lon_grid(:)/100 BP.lat_grid(:)/100],...
    [[NElon_list(1);NElon_list(end)] [NElat_list(1);NElat_list(end)]]);

% find mean around start and end point
[row,col]=ind2sub(size(BP.mean_map),idx);
first=mean(BP.mean_map(row(1)-2:row(1)+2,col(1)-2:col(1)+2),'all');
last=mean(BP.mean_map(row(2)-2:row(2)+2,col(2)-2:col(2)+2),'all');
NEtt=last-first; %days

fprintf(fid,'Lagrangian Southern inflow\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',NE_path.dist);
fprintf(fid,'Speed=%4.2f cm/s\r\n',NE_path.dist/NEtt/86400*1e5); %cm/s
fprintf(fid,'Transit time=%4.0f days\r\n\r\n',NEtt);

% Find first and last point on path
idx=dsearchn([BP.lon_grid(:)/100 BP.lat_grid(:)/100],...
    [[DPlon_list(1);DPlon_list(end)] [DPlat_list(1);DPlat_list(end)]]);

% find mean around start and end point
[row,col]=ind2sub(size(DP.mean_map),idx);
first=nanmean(DP.mean_map(row(1)-2:row(1)+2,col(1)-2:col(1)+2),'all');
last=nanmean(DP.mean_map(row(2)-2:row(2)+2,col(2)-2:col(2)+2),'all');
DPtt=last-first; %days

fprintf(fid,'Lagrangian Northern inflow\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',DP_path.dist);
fprintf(fid,'Speed=%4.2f cm/s\r\n',DP_path.dist/DPtt/86400*1e5); %cm/s
fprintf(fid,'Transit time=%4.0f days\r\n\r\n',DPtt);

% Find all ages in strait
load PT_AS_inlets.mat
station_bound=alphaShape(AS.lon,AS.lat);
station_bound.Alpha=station_bound.Alpha*2;
inside_idx=inShape(station_bound,BP.lon_grid/100,BP.lat_grid/100);

straitA=BP.mean_map(inside_idx);
straitLa=BP.lat_grid(inside_idx)/100;
straitLo=BP.lat_grid(inside_idx)/100;

% create lines
center_line=[linspace(-125.15,-122.75,1000);linspace(50.15,48.4,1000)];
load centre_line_dist.mat

distance_idx=NaN(1,length(straitA));
distance_N=NaN(1,length(straitA));
for i=1:length(straitA)
    [~,idxb]=min(abs(center_line(2,:)-straitLa(i)));
    distance_N(i)=center_line_dist(idxb);
end

% dist=gsw_distance([ones(length(straitA),1)*-122.4364 straitLo],...
%     [ones(length(straitA),1)*48.1429 straitLa]);
% dist=(dist-min(dist))/1000;

distance_N(isnan(straitA))=[];
straitA(isnan(straitA))=[];

[coeffs(1), coeffs(2)] = TheilSen([distance_N'/1000,straitA]);
fittedy=polyval(coeffs,distance_N/1000);

figure;
subplot(1,3,1);
scatter(straitA,distance_N/1000,30,'r','filled',...
    'markeredgecolor','none','markerfacealpha',0.1)
hold on
plot(fittedy,distance_N/1000,'k');
xlabel('Age (days)');
ylabel('Distance (km)');
title('Lagrangian');

s=TheilSen([straitA,distance_N'/1000]); % km/day
[r,p]=corrcoef(distance_N'/1000,straitA);

fprintf(fid,'Lagrangian Whole inflow\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',max(distance_N/1000));
fprintf(fid,'Speed=%4.2f cm/s\r\n',s/86400*1e5); %cm/s
fprintf(fid,'Transit time=%4.0f days\r\n',max(fittedy));
fprintf(fid,'r=%4.2f, p=%4.2f\r\n\r\n',r(2),p(2));

%% Tracer- observations
clearvars -except NE_path DP_path fid Eu NElon_list NElat_list...
    DPlon_list DPlat_list
core=load('sinfit_temp_data_plotting_core.mat');
%% kriging

% Set up XY projection
lat_limV=[min(core.grid_lat(:)) max(core.grid_lat(:))];
lon_limV=[min(core.grid_lon(:)) max(core.grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(core.station_lon,core.station_lat);
[Xgrid,Ygrid]=m_ll2xy(core.grid_lon,core.grid_lat);

dday = variogram([X' Y'],core.cold_day','maxdist',49000,'plotit',false);

[a,c,~,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,X',Y',core.cold_day',Xgrid,Ygrid);

% Find first and last point on path
idx=dsearchn([core.grid_lon(:) core.grid_lat(:)],...
    [[NElon_list(1);NElon_list(end)] [NElat_list(1);NElat_list(end)]]);

% find mean around start and end point
[row,col]=ind2sub(size(zi),idx);
first=mean(zi(row(1)-2:row(1)+2,col(1)-2:col(1)+2),'all');
last=mean(zi(row(2)-2:row(2)+2,col(2)-2:col(2)+2),'all');
NEtt=last-first; %days

fprintf(fid,'Tracer observations Southern inflow\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',NE_path.dist);
fprintf(fid,'Speed=%4.2f cm/s\r\n',NE_path.dist/NEtt/86400*1e5); %cm/s
fprintf(fid,'Transit time=%4.0f days\r\n\r\n',NEtt);

% Find first and last point on path
idx=dsearchn([core.grid_lon(:) core.grid_lat(:)],...
    [[DPlon_list(1);DPlon_list(end)] [DPlat_list(1);DPlat_list(end)]]);

% find mean around start and end point
[row,col]=ind2sub(size(zi),idx);
first=nanmean(zi(row(1)-2:row(1)+2,col(1)-2:col(1)+2),'all');
last=nanmean(zi(row(2)-2:row(2)+2,col(2)-2:col(2)+2),'all');
DPtt=last-first; %days

fprintf(fid,'Tracer observations Northern inflow\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',DP_path.dist);
fprintf(fid,'Speed=%4.2f cm/s\r\n',DP_path.dist/DPtt/86400*1e5); %cm/s
fprintf(fid,'Transit time=%4.0f days\r\n',DPtt);
fprintf(fid,'!PROBABLY NOT VALID!\r\n\r\n');

% Find youngest water in Haro Strait
HS_AS=load('HS_AS.mat');
AS=alphaShape(HS_AS.l,HS_AS.ll);
msk=inShape(AS,core.Dsort_lon,core.Dsort_lat);
age=core.Dsort_coldday-min(core.Dsort_coldday(msk));

[coeffs(1), coeffs(2)] = TheilSen([core.distance_N'/1000 age']);
fittedy=polyval(coeffs,core.distance_N/1000);
[r,p]=corrcoef(core.distance_N'/1000,age);

subplot(1,3,2);
scatter(age,core.distance_N/1000,30,'r','filled',...
    'markeredgecolor','none','markerfacealpha',0.1)
hold on
plot(fittedy,core.distance_N/1000,'k');
title('Tracer (observations)');

s=TheilSen([age',core.distance_N'/1000]); % km/day

fprintf(fid,'Tracer observations Whole inflow\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',max(core.distance_N/1000));
fprintf(fid,'Speed=%4.2f cm/s\r\n',s/86400*1e5); %cm/s
fprintf(fid,'Transit time=%4.0f days\r\n',max(fittedy));
fprintf(fid,'r=%4.2f, p=%4.2f\r\n\r\n',r(2),p(2));

%% Tracer- model
clearvars -except NE_path DP_path fid Eu NElon_list NElat_list...
    DPlon_list DPlat_list
m_core=load('sinfit_SSC_modelres_core.mat');

%% kriging
HS_AS=load('HS_AS.mat');
AS=alphaShape(HS_AS.l,HS_AS.ll);

qlon=m_core.grid_lon;qlon(isnan(qlon))=0;
qlat=m_core.grid_lat;qlat(isnan(qlat))=0;
msk=inShape(AS,qlon,qlat);
m_core.cold_day=m_core.cold_day-min(m_core.cold_day(msk));

% Find first and last point on path
idx=dsearchn([m_core.grid_lon(:) m_core.grid_lat(:)],...
    [[NElon_list(1);NElon_list(end)] [NElat_list(1);NElat_list(end)]]);

% find mean around start and end point
[row,col]=ind2sub(size(m_core.cold_day),idx);
first=mean(m_core.cold_day(row(1)-2:row(1)+2,col(1)-2:col(1)+2),'all');
last=mean(m_core.cold_day(row(2)-2:row(2)+2,col(2)-2:col(2)+2),'all');
NEtt=last-first; %days

fprintf(fid,'Tracer model Southern inflow\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',NE_path.dist);
fprintf(fid,'Speed=%4.2f cm/s\r\n',NE_path.dist/NEtt/86400*1e5); %cm/s
fprintf(fid,'Transit time=%4.0f days\r\n\r\n',NEtt);

% Find first and last point on path
idx=dsearchn([m_core.grid_lon(:) m_core.grid_lat(:)],...
    [[DPlon_list(1);DPlon_list(end)] [DPlat_list(1);DPlat_list(end)]]);

% find mean around start and end point
[row,col]=ind2sub(size(m_core.cold_day),idx);
first=nanmean(m_core.cold_day(row(1)-2:row(1)+2,col(1)-2:col(1)+2),'all');
last=nanmean(m_core.cold_day(row(2)-2:row(2)+2,col(2)-2:col(2)+2),'all');
DPtt=last-first; %days

fprintf(fid,'Tracer model Northern inflow\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',DP_path.dist);
fprintf(fid,'Speed=%4.2f cm/s\r\n',DP_path.dist/DPtt/86400*1e5); %cm/s
fprintf(fid,'Transit time=%4.0f days\r\n',DPtt);
fprintf(fid,'!PROBABLY NOT VALID!\r\n\r\n');

% Find all ages in strait
load PT_AS_inlets.mat
station_bound=alphaShape(AS.lon,AS.lat);
station_bound.Alpha=station_bound.Alpha*2;
inside_idx=inShape(station_bound,qlon,qlat);

straitA=m_core.cold_day(inside_idx);
straitLa=m_core.grid_lat(inside_idx);
msk=isnan(straitA); 
straitA(msk)=[];
straitLa(msk)=[];

% create lines
center_line=[linspace(-125.15,-122.75,1000);linspace(50.15,48.4,1000)];
load centre_line_dist.mat

distance_idx=NaN(1,length(straitA));
distance_N=NaN(1,length(straitA));
for i=1:length(straitA)
    [~,idxb]=min(abs(center_line(2,:)-straitLa(i)));
    distance_N(i)=center_line_dist(idxb);
end

% [coeffs(1), coeffs(2)] = tsreg(distance_N'/1000,straitA',...
%     length(distance_N));
% fittedy=polyval(coeffs,distance_N/1000);

[coeffs(1), coeffs(2)] = TheilSen([distance_N'/1000,straitA]);
fittedy=polyval(coeffs,distance_N/1000);

subplot(1,3,3);
scatter(straitA,distance_N/1000,30,'r','filled',...
    'markeredgecolor','none','markerfacealpha',0.1)
hold on
plot(fittedy,distance_N/1000,'k');
[r,p]=corrcoef(distance_N'/1000,straitA);
title('Tracer (model)');

s=TheilSen([straitA,distance_N'/1000]);

fprintf(fid,'Tracer model whole inflow\r\n');
fprintf(fid,'Distance=%4.1f km\r\n',max(distance_N/1000));
fprintf(fid,'Speed=%4.2f cm/s\r\n',s/86400*1e5); %cm/s
fprintf(fid,'Transit time=%4.0f days\r\n',max(fittedy));
fprintf(fid,'r=%4.2f, p=%4.2f\r\n\r\n',r(2),p(2));

export_fig /ocean/sstevens/IW_project/figures/all_TT_trend.png