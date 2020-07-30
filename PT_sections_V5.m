clear
addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));
%% Load bathy data
fname='/ocean/rich/more/mmapbase/noaa_bc3/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');
lon=ncread(fname,'lon');
lat_lim=[48 50.3];
lon_lim=[-125.99298287 -122.22345235];
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
Z_lon=lon(ilon);
Z_lat=lat(ilat);
% [Z_lon, Z_lat]=meshgrid(Z_lon,Z_lat);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);

res=1/24; % 1-hr traj positions

HS_box=[[-123.2727;-122.8011;-123.1087;-123.5899] ...
    [48.7141;48.9287;49.0628;48.7810]];
HS_as=alphaShape(HS_box(:,1),HS_box(:,2));
SB_box=[[HS_box(3:4,1);-123.9969;-124.3668] ...
    [HS_box(3:4,2);49.5783;49.3051]];
SB_as=alphaShape(SB_box(:,1),SB_box(:,2));
NB_box=[[SB_box(3:4,1);-125.3255;-124.6681] ...
    [SB_box(3:4,2);50.0033;50.1915]];
NB_as=alphaShape(NB_box(:,1),NB_box(:,2));

allB_as=alphaShape([[NB_box(:,1);HS_box(:,1);SB_box(:,1)] ...
    [NB_box(:,2);HS_box(:,2);SB_box(:,2)]]);

m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection

dir_list=dir('/ocean/sstevens/ariane_old/results_BP/2016*365d');
tic
SB_pass_IN=NaN(30000,4);SB_IN_count=1;
SB_pass_OUT=NaN(30000,4);SB_OUT_count=1;
NB_pass_IN=NaN(30000,4);NB_IN_count=1;
NB_pass_OUT=NaN(30000,4);NB_OUT_count=1;

count2=0;
for i=1:length(dir_list)
    
    % Load data
    path = [dir_list(i).folder '/'  dir_list(i).name '/'];
    filestr=[path 'TrajTBL'];
    load(filestr);
    disp(i)
    toc
    
    for j=1:length(TrajTBL)
        count=0;
        
        tmt=datenum(2016,i,ceil(TrajTBL(j).initHour/24));
        DOY=day(datetime(datestr(tmt)),'dayofyear');
        
        mn_traj=[];
        mn_traj(:,1)=movmean(TrajTBL(j).lon,3*24,'omitnan');
        mn_traj(:,2)=movmean(TrajTBL(j).lat,3*24,'omitnan');
        mn_traj(:,3)=movmean(TrajTBL(j).z_m,3*24,'omitnan');
        mn_traj(:,4)=ones(size(mn_traj,1),1)*DOY;
                
        msk=isnan(mn_traj(:,1));
        if sum(msk)
            mn_traj(msk,:)=0;
        end
        
         % Does it go into SoG?
         SB_msk=inShape(SB_as,mn_traj(:,1),mn_traj(:,2));
         
         if sum(SB_msk) % if it enter SoG, where does it enter?
%              figure;
%              m_gshhs_f('patch',rgb('light grey'),'edgecolor','k');
%              hold on;
%              m_plot(mn_traj(:,1),mn_traj(:,2),'color','r');
             
             HS_msk=inShape(HS_as,mn_traj(:,1),mn_traj(:,2));
             NB_msk=inShape(NB_as,mn_traj(:,1),mn_traj(:,2));
             all_msk=~inShape(allB_as,mn_traj(:,1),mn_traj(:,2));
             
             transit=zeros(length(HS_msk),1);
             transit(all_msk)=NaN;
             transit(HS_msk)=2;
             transit(SB_msk)=3;
             transit(NB_msk)=5;

             transit=diff(transit);
             
             % THIS ONLY FINDS THE FIRST CROSSING
             % SHOULD DO TEST WITH ALL CROSSINGS TOO
             % 1= up-strait crossing at SB
             idx=find(transit==1,1,'first');
             if ~isempty(idx)
                 SB_pass_IN(SB_IN_count,:)=mn_traj(idx,:);
                 SB_IN_count=SB_IN_count+1;
                 
             end
             
             %-1= down-strait crossing at SB
             idx=find(transit==-1,1,'first');
             if ~isempty(idx)
                 SB_pass_OUT(SB_OUT_count,:)=mn_traj(idx,:);
                 SB_OUT_count=SB_OUT_count+1;
             end

             % 2= up-strait crossing at NB
             idx=find(transit==2,1,'first');
             if ~isempty(idx)
                 NB_pass_IN(NB_IN_count,:)=mn_traj(idx,:);
                 NB_IN_count=NB_IN_count+1;
             end
             
             %-2= down-strait crossing at NB
             idx=find(transit==-2,1,'first');
             if ~isempty(idx)
                 NB_pass_OUT(NB_OUT_count,:)=mn_traj(idx,:);
                 NB_OUT_count=NB_OUT_count+1;
             end
         end
    end
end
keyboard
save('/ocean/sstevens/ariane_old/results_BP/PT_sections.mat');

m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
figure
hold on
m_gshhs_f('patch',rgb('light grey'),'edgecolor','k');
[lo,la]=m_ll2xy(HS_box(:,1),HS_box(:,2));
tmp_as=alphaShape(lo,la);
plot(tmp_as,'FaceColor','red','FaceAlpha',0.1,'edgecolor','none');
[lo,la]=m_ll2xy(SB_box(:,1),SB_box(:,2));
tmp_as=alphaShape(lo,la);
plot(tmp_as,'FaceColor','b','FaceAlpha',0.1,'edgecolor','none');
[lo,la]=m_ll2xy(NB_box(:,1),NB_box(:,2));
tmp_as=alphaShape(lo,la);
plot(tmp_as,'FaceColor','g','FaceAlpha',0.1,'edgecolor','none');
[lo,la]=m_ll2xy([NB_box(:,1);HS_box(:,1);SB_box(:,1)],...
    [NB_box(:,2);HS_box(:,2);SB_box(:,2)]);
tmp_as=alphaShape(lo,la);
plot(tmp_as,'FaceColor','g','FaceAlpha',0.0,'edgecolor','k');
m_scatter(SB_pass_IN(:,1),SB_pass_IN(:,2),1,'k');
m_scatter(SB_pass_OUT(:,1),SB_pass_OUT(:,2),1,'r');
m_scatter(NB_pass_OUT(:,1),NB_pass_OUT(:,2),1,'r');
m_scatter(NB_pass_IN(:,1),NB_pass_IN(:,2),1,'k');

% Find cross-strait distance from western side of gate
idx=find(isnan(SB_pass_IN(:,1)),1)-1;
SB_pass_IN(1:idx,end+1)=gsw_distance([ones(idx,1)*HS_box(4,1)...
    SB_pass_IN(1:idx,1)],[ones(idx,1)*HS_box(4,2) SB_pass_IN(1:idx,2)])/1000;

idx=find(isnan(SB_pass_OUT(:,1)),1)-1;
SB_pass_OUT(1:idx,end+1)=gsw_distance([ones(idx,1)*HS_box(4,1)...
    SB_pass_OUT(1:idx,1)],[ones(idx,1)*HS_box(4,2) SB_pass_OUT(1:idx,2)])/1000;

idx=find(isnan(NB_pass_IN(:,1)),1)-1;
NB_pass_IN(1:idx,end+1)=gsw_distance([ones(idx,1)*SB_box(4,1)...
    NB_pass_IN(1:idx,1)],[ones(idx,1)*SB_box(4,2) NB_pass_IN(1:idx,2)])/1000;

idx=find(isnan(NB_pass_OUT(:,1)),1)-1;
NB_pass_OUT(1:idx,end+1)=gsw_distance([ones(idx,1)*SB_box(4,1)...
    NB_pass_OUT(1:idx,1)],[ones(idx,1)*SB_box(4,2) NB_pass_OUT(1:idx,2)])/1000;

SB_pass_IN(end,:)=NaN;
SB_pass_OUT(end,:)=NaN;
NB_pass_IN(end,:)=NaN;
NB_pass_IN(end,:)=NaN;

figure;
hist(NB_pass_OUT(:,5))
hold on
hist(NB_pass_IN(:,5))

figure;
hist(SB_pass_OUT(:,5))
hold on
hist(SB_pass_IN(:,5))

save('/ocean/sstevens/ariane_old/results_BP/PT_sections_V3.mat');

%% Plot stats
clear
load /ocean/sstevens/ariane_old/results_BP/PT_sections_V3.mat
%%
SB_pass_IN(SB_pass_IN(:,5)==0,5)=NaN;
NB_pass_IN(NB_pass_IN(:,5)==0,5)=NaN;
SB_pass_OUT(SB_pass_OUT(:,5)==0,5)=NaN;
NB_pass_OUT(NB_pass_OUT(:,5)==0,5)=NaN;

% Find section bathymetry
[X,Y]=meshgrid(Z_lat,Z_lon);
SB_sec=[linspace(HS_box(4,1),HS_box(3,1),1000);...
    linspace(HS_box(4,2),HS_box(3,2),1000)];
NB_sec=[linspace(SB_box(4,1),SB_box(3,1),1000);...
    linspace(SB_box(4,2),SB_box(3,2),1000)];

SB_bathy=interp2(X,Y,Z,SB_sec(2,:),SB_sec(1,:))*-1;
NB_bathy=interp2(X,Y,Z,NB_sec(2,:),NB_sec(1,:))*-1;
SB_bathy(isnan(SB_bathy))=0;SB_bathy=smooth(SB_bathy,20);
NB_bathy(isnan(NB_bathy))=0;NB_bathy=smooth(NB_bathy,20);

idx=length(SB_sec);
SB_secd=gsw_distance([ones(idx,1)*HS_box(4,1)...
    SB_sec(1,:)'],[ones(idx,1)*HS_box(4,2) SB_sec(2,:)'])/1000;
NB_secd=gsw_distance([ones(idx,1)*SB_box(4,1)...
    NB_sec(1,:)'],[ones(idx,1)*SB_box(4,2) NB_sec(2,:)'])/1000;

SB_linsec(1,:)=interp1(SB_secd,SB_sec(1,:),0:0.5:max(SB_secd));
SB_linsec(2,:)=interp1(SB_secd,SB_sec(2,:),0:0.5:max(SB_secd));

NB_linsec(1,:)=interp1(NB_secd,NB_sec(1,:),0:0.5:max(NB_secd));
NB_linsec(2,:)=interp1(NB_secd,NB_sec(2,:),0:0.5:max(NB_secd));


SB_IN_section=histcounts2(SB_pass_IN(:,3)*-1,SB_pass_IN(:,5),...
    0:5:450,0:0.5:max(SB_secd));
    
SB_OUT_section=histcounts2(SB_pass_OUT(:,3)*-1,SB_pass_OUT(:,5)...
    ,0:5:450,0:0.5:max(SB_secd));

NB_IN_section=histcounts2(NB_pass_IN(:,3)*-1,NB_pass_IN(:,5),...
    0:5:450,0:0.5:max(NB_secd));

NB_OUT_section=histcounts2(NB_pass_OUT(:,3)*-1,NB_pass_OUT(:,5),...
    0:5:450,0:0.5:max(NB_secd));

IW_count=sum(sum(SB_IN_section(50/5:150/5,:)));
IW_pcnt=IW_count/sum(SB_IN_section(:));

IW_count=sum(sum(SB_OUT_section(1:50/5,:)));
IW_pcnt=IW_count/sum(SB_OUT_section(:));

B_pcnt=sum(NB_IN_section(:))/sum(SB_IN_section(:));

IW_count=sum(sum(SB_OUT_section(1:50/5,:)));
IW_pcnt=IW_count/sum(SB_OUT_section(:));

Ross=[sqrt((9.81*200)/gsw_f(49)) sqrt((9.81*400)/gsw_f(49))];

%% 
IN_cm=cmocean('amp');IN_cm(1,:)=[1 1 1];
OUT_cm=cmocean('-ice');OUT_cm(1,:)=[1 1 1];

f1=figure('units','centimeters','outerposition',[0 0 19 16],'color','w');
map_ax=axes('Position',[0.065 0.025 0.35 0.35]);
lon_lim=[-125.5 -122.522345235];
lat_lim=[50.1 48.5];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.25)
hold on
m_gshhs_f('patch',rgb_x('light grey'),'edgecolor','none');
m_northarrow(-125.32,49.75,0.15,'type',2,'linewidth',.5);
m_grid('linestyle','none','tickdir','out','xaxisloc','bottom','yaxisloc','left','fontsize',8,...
    'xtick',-126:-122,'ytick',48:51);
m_scatter(HS_box(3:4,1),HS_box(3:4,2),10,'k','filled');
m_scatter(SB_box(3:4,1),SB_box(3:4,2),10,'k','filled');
m_patch([-124.2723;-123.6640;-124.1356;-124.7425],...
    [49.0594;49.5403;49.7926;49.3092],'m','facealpha',0.1,'edgecolor','k');
m_patch([-123.4397;-122.8289;-123.2230;-123.8327],...
    [48.6058;49.0823;49.3003;48.8217],'m','facealpha',0.1,'edgecolor','k');
text(0.7,0.1,'g)','units','normalized','fontweight','bold','fontsize',8);

slice_ax1=axes('Position',[0.385 0.2 0.45 0.2]);
lon_lim=[-123.6640 -124.7425];
lat_lim=[49.7926 49.0594];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','horiz','aspect',0.3)
hold on
m_gshhs_f('patch',rgb_x('light grey'),'edgecolor','none');
m_grid('linestyle','none','tickdir','out','xaxisloc','bottom','yaxisloc','right','fontsize',8,...
    'xtick',[],'ytick',[]);
m_scatter(HS_box(3:4,1),HS_box(3:4,2),10,'k','filled');
% m_plot(HS_box(3:4,1),HS_box(3:4,2),'k');
% m_plot(SB_box(3:4,1),SB_box(3:4,2),'k');
text(0.9,0.1,'h)','units','normalized','fontweight','bold','fontsize',8);

Bhist3(:,1)=histcounts(NB_pass_IN(:,5),0:0.5:max(NB_secd));
Bhist3(:,2)=histcounts(NB_pass_OUT(:,5),0:0.5:max(NB_secd));
scle=max(Bhist3(:));

for i=1:length(NB_linsec)-1
    if Bhist3(i,1)>0
        [reclon,reclat]=m_ll2xy(NB_linsec(1,i),NB_linsec(2,i));
        [reclon2,~]=m_ll2xy(NB_linsec(1,i+1),NB_linsec(2,i+1));
        rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist3(i,1)/scle)/400],...
            'facecolor',[IN_cm(end/2,:) 0.4],'edgecolor','none');
    end
    
     if Bhist3(i,2)>0
        [reclon,reclat]=m_ll2xy(NB_linsec(1,i),NB_linsec(2,i));
        [reclon2,~]=m_ll2xy(NB_linsec(1,i+1),NB_linsec(2,i+1));
        a=rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist3(i,2)/scle)/400],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        rectangle('Position',[reclon,reclat-(a.Position(4)),reclon2-reclon,a.Position(4)],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        delete(a)
    end
end

m_scatter(nanmean(NB_pass_IN(:,1)),nanmean(NB_pass_IN(:,2)),5,'filled',...
    'markerfacecolor',IN_cm(end/2,:),'markeredgecolor','none');
m_scatter(nanmean(NB_pass_OUT(:,1)),nanmean(NB_pass_OUT(:,2)),5,'filled',...
    'markerfacecolor',OUT_cm(end/2,:),'markeredgecolor','none');

[reclon,reclat]=m_ll2xy(NB_linsec(1,end-1),NB_linsec(2,end-1));
[reclon2,~]=m_ll2xy(NB_linsec(1,end),NB_linsec(2,end));
a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*2,(350/scle)/400],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*2,a.Position(4)],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*2,(700/scle)/400],...
    'facecolor',[0 0 0 0.4],'edgecolor','none');
rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*2,a.Position(4)],...
    'facecolor',[0 0 0 0.4],'edgecolor','none');

m_text([-124.1263;-124.0548;-123.9797;-123.9120;-123.8408],...
    [49.6771;49.6326;49.5859;49.5435;49.4989],...
    {700;350;0;350;700});

slice_ax2=axes('Position',[0.385 0. 0.45 0.2]);
lon_lim=[-122.8289 -123.8327];
lat_lim=[49.3003 48.6058];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','horiz','aspect',0.3)
hold on
m_gshhs_f('patch',rgb_x('light grey'),'edgecolor','none');
m_grid('linestyle','none','tickdir','out','xaxisloc','bottom','yaxisloc','right','fontsize',8,...
    'xtick',[],'ytick',[]);
% m_scatter(HS_box(3:4,1),HS_box(3:4,2),10,'k','filled');
text(0.9,0.1,'i)','units','normalized','fontweight','bold','fontsize',8);

Bhist4(:,1)=histcounts(SB_pass_IN(:,5),0:0.5:max(SB_secd));
Bhist4(:,2)=histcounts(SB_pass_OUT(:,5),0:0.5:max(SB_secd));
scle=max(Bhist3(:));

for i=1:length(SB_linsec)-1
    if Bhist4(i,1)>0

        [reclon,reclat]=m_ll2xy(SB_linsec(1,i),SB_linsec(2,i));
        [reclon2,~]=m_ll2xy(SB_linsec(1,i+1),SB_linsec(2,i+1));
        rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist4(i,1)/scle)/650],...
            'facecolor',[IN_cm(end/2,:) 0.4],'edgecolor','none');
    end
    
     if Bhist4(i,2)>0
        [reclon,reclat]=m_ll2xy(SB_linsec(1,i),SB_linsec(2,i));
        [reclon2,~]=m_ll2xy(SB_linsec(1,i+1),SB_linsec(2,i+1));
        a=rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist4(i,2)/scle)/650],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        rectangle('Position',[reclon,reclat-(a.Position(4)),reclon2-reclon,a.Position(4)],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        delete(a)
    end
end

m_scatter(nanmean(SB_pass_IN(:,1)),nanmean(SB_pass_IN(:,2)),5,'filled',...
    'markerfacecolor',IN_cm(end/2,:),'markeredgecolor','none');
m_scatter(nanmean(SB_pass_OUT(:,1)),nanmean(SB_pass_OUT(:,2)),5,'filled',...
    'markerfacecolor',OUT_cm(end/2,:),'markeredgecolor','none');

[reclon,reclat]=m_ll2xy(SB_linsec(1,end-1),SB_linsec(2,end-1));
[reclon2,~]=m_ll2xy(SB_linsec(1,end),SB_linsec(2,end));
a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*1.5,(450/scle)/650],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*1.5,a.Position(4)],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*1.5,(900/scle)/650],...
    'facecolor',[0 0 0 0.4],'edgecolor','none');
rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*1.5,a.Position(4)],...
    'facecolor',[0 0 0 0.4],'edgecolor','none');
a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*1.5,(1350/scle)/650],...
    'facecolor',[0 0 0 0.1],'edgecolor','none');
rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*1.5,a.Position(4)],...
    'facecolor',[0 0 0 0.1],'edgecolor','none');

m_text([-123.2579;-123.2076;-123.1539;-123.0970;-123.0402;-122.9868],...
    [49.1762;49.1448;49.1114;49.0758;49.0403;49.0067],...
    {1350;900;450;0;450;900});

ax(1)=axes('Position',[0.225 0.7 0.25 0.15]);
hold on
[X,Y]=meshgrid([0:0.5:max(SB_secd)-0.5]',0:5:450-5);
[CS,h]=contourf(X,Y,SB_IN_section,0:5:max([NB_IN_section(:);...
    SB_IN_section(:)]),'linestyle','none');
caxis([0 max([NB_IN_section(:);SB_IN_section(:)])]);
set(gca,'ydir','reverse','xticklabel',[]);
colormap(ax(1),IN_cm);
area(SB_secd,SB_bathy,400,'FaceColor',rgb_x('light grey'),'edgecolor','none');
ylim([0 200]);
ylabel('Depth (m)');
title('Into Basin');
text(0.9,0.1,'a)','units','normalized','fontweight','bold','fontsize',8);
text(-0.2,0.5,'SB','units','normalized','fontweight','bold','fontsize',8);

ax(2)=axes('Position',[0.5 0.7 0.25 0.15]);
hold on
contourf(X,Y,SB_OUT_section,0:5:max([NB_IN_section(:);...
    SB_OUT_section(:)]),'linestyle','none');
set(gca,'ydir','reverse','xticklabel',[],'yticklabel',[]);
caxis([0 max([NB_OUT_section(:);SB_OUT_section(:)])]);
colormap(ax(2),OUT_cm);
area(SB_secd,SB_bathy,400,'FaceColor',rgb_x('light grey'),'edgecolor','none');
ylim([0 200]);
title('Out of Basin');
text(0.9,0.1,'b)','units','normalized','fontweight','bold','fontsize',8);

ax(3)=axes('Position',[0.225 0.5 0.25 0.15]);
hold on
[X,Y]=meshgrid([0:0.5:max(NB_secd)-0.5]',0:5:450-5);
[CS,h]=contourf(X,Y,NB_IN_section,0:5:max([NB_IN_section(:);...
    SB_IN_section(:)]),'linestyle','none');
set(gca,'ydir','reverse');
caxis([0 max([NB_IN_section(:);SB_IN_section(:)])]);
colormap(ax(3),IN_cm);
cm_ax1=m_contfbar(ax(3),[0 1],-0.425,CS,h,'endpiece','no','axfrac',.1,'fontsize',8);
area(NB_secd,NB_bathy,430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
area([max(NB_secd) max(SB_secd)],[0 0],430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
linkaxes([ax(1) ax(3)],'x');
ylim([0 430]);
xlabel('Across-Strait Distance (km)');
text(0.9,0.1,'d)','units','normalized','fontweight','bold','fontsize',8);
text(-0.2,0.5,'NB','units','normalized','fontweight','bold','fontsize',8);

ax(4)=axes('Position',[0.5 0.5 0.25 0.15]);
hold on
[CS,h]=contourf(X,Y,NB_OUT_section,0:5:max([NB_OUT_section(:);...
    SB_OUT_section(:)]),'linestyle','none');
set(gca,'ydir','reverse','yticklabel',[]);
caxis([0 max([SB_OUT_section(:);NB_OUT_section(:)])]);
colormap(ax(4),OUT_cm);
area(NB_secd,NB_bathy,430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
% xlim([0 max(SB_secd)]);
area([max(NB_secd) max(SB_secd)],[0 0],430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
cm_ax2=m_contfbar(ax(4),[0 1],-0.425,CS,h,'endpiece','no','axfrac',.1,'fontsize',8);
linkaxes([ax(2) ax(4)],'x');
ylim([0 430]);
text(0.9,0.1,'e)','units','normalized','fontweight','bold','fontsize',8);

% Add mini histograms
axh(1)=axes('Position',[0.76 0.7 0.075 0.15]);
hold on
Bhist1(:,1)=histcounts(SB_pass_IN(:,3)*-1,0:5:200);
Bhist1(:,2)=histcounts(SB_pass_OUT(:,3)*-1,0:5:200);
b1=barh(0:5:200-5,Bhist1(:,1),'histc');
set(b1,'facecolor',IN_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
b2=barh(0:5:200-5,Bhist1(:,2),'histc');
set(b2,'facecolor',OUT_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
linkaxes([ax(1) ax(2) axh(1)],'y');
set(gca,'ydir','reverse','yticklabel',[],'xtick',0:500:1000,...
    'xticklabel',[]);
% plot([0 max(xlim)],[nanmean(SB_pass_IN(:,3)*-1) ...
%     nanmean(SB_pass_IN(:,3))*-1],':','color',IN_cm(end/2,:));
% plot([0 max(xlim)],[nanmean(SB_pass_OUT(:,3)*-1) ...
%     nanmean(SB_pass_OUT(:,3))*-1],':','color',OUT_cm(end/2,:));
scatter(0,nanmean(SB_pass_IN(:,3)*-1),15,'filled','markerfacecolor',...
    IN_cm(end/2,:),'markeredgecolor','none');
scatter(0,geomean(SB_pass_IN(~isnan(SB_pass_IN(:,3)),3)*-1,'all'),15,...
    'markeredgecolor',IN_cm(end/2,:));
scatter(0,nanmean(SB_pass_OUT(:,3)*-1),15,'filled','markerfacecolor',...
    OUT_cm(end/2,:),'markeredgecolor','none');
scatter(0,geomean(SB_pass_OUT(~isnan(SB_pass_OUT(:,3)),3)*-1,'all'),15,...
    'markeredgecolor',OUT_cm(end/2,:));
grid on
text(0.9,0.1,'c)','units','normalized','fontweight','bold','fontsize',8);

axh(2)=axes('Position',[0.76 0.5 0.075 0.15]);
hold on
Bhist2(:,1)=histcounts(NB_pass_IN(:,3)*-1,0:5:430);
Bhist2(:,2)=histcounts(NB_pass_OUT(:,3)*-1,0:5:430);
b1=barh(0:5:430-5,Bhist2(:,1),'histc');
set(b1,'facecolor',IN_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
b2=barh(0:5:430-5,Bhist2(:,2),'histc');
set(b2,'facecolor',OUT_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
set(gca,'ydir','reverse','yticklabel',[],'xtick',0:500:1000,...
    'xticklabelrotation',45);
linkaxes([ax(3) ax(4) axh(2)],'y');
linkaxes([axh(1) axh(2)],'x');
scatter(0,nanmean(NB_pass_IN(:,3)*-1),15,'filled','markerfacecolor',...
    IN_cm(end/2,:),'markeredgecolor','none');
scatter(0,nanmean(NB_pass_OUT(:,3)*-1),15,'filled','markerfacecolor',...
    OUT_cm(end/2,:),'markeredgecolor','none');
grid on
xlabel('Particle Count');
text(0.9,0.1,'f)','units','normalized','fontweight','bold','fontsize',8);

set(findall(gcf,'-property','tickdir'),'tickdir','out');
set(findall(gcf,'-property','ticklength'),'ticklength',[0.025 0.02]);
set(findall(gcf,'-property','fontsize'),'fontsize',6);

%%
% export_fig '/ocean/sstevens/IW_project/figures/paper/PT_sections.pdf' -dpdf
export_fig '/ocean/sstevens/IW_project/figures/paper/PT_sections.png' -png -m3

%% Plot stats- SUMMER
clear
load /ocean/sstevens/ariane_old/results_BP/PT_sections_V3.mat
%%
SB_pass_IN(SB_pass_IN(:,5)==0,5)=NaN;
NB_pass_IN(NB_pass_IN(:,5)==0,5)=NaN;
SB_pass_OUT(SB_pass_OUT(:,5)==0,5)=NaN;
NB_pass_OUT(NB_pass_OUT(:,5)==0,5)=NaN;

MAY=day(datetime(datestr(datenum(2016,5,1))),'dayofyear');
SEPT=day(datetime(datestr(datenum(2016,9,1))),'dayofyear');

idx=SB_pass_IN(:,4)>=MAY & SB_pass_IN(:,4)<SEPT;SB_pass_IN(~idx,:)=NaN;
idx=NB_pass_IN(:,4)>=MAY & NB_pass_IN(:,4)<SEPT;NB_pass_IN(~idx,:)=NaN;
idx=SB_pass_OUT(:,4)>=MAY & SB_pass_OUT(:,4)<SEPT;SB_pass_OUT(~idx,:)=NaN;
idx=NB_pass_OUT(:,4)>=MAY & NB_pass_OUT(:,4)<SEPT;NB_pass_OUT(~idx,:)=NaN;

% Find section bathymetry
[X,Y]=meshgrid(Z_lat,Z_lon);
SB_sec=[linspace(HS_box(4,1),HS_box(3,1),1000);...
    linspace(HS_box(4,2),HS_box(3,2),1000)];
NB_sec=[linspace(SB_box(4,1),SB_box(3,1),1000);...
    linspace(SB_box(4,2),SB_box(3,2),1000)];

SB_bathy=interp2(X,Y,Z,SB_sec(2,:),SB_sec(1,:))*-1;
NB_bathy=interp2(X,Y,Z,NB_sec(2,:),NB_sec(1,:))*-1;
SB_bathy(isnan(SB_bathy))=0;SB_bathy=smooth(SB_bathy,20);
NB_bathy(isnan(NB_bathy))=0;NB_bathy=smooth(NB_bathy,20);

idx=length(SB_sec);
SB_secd=gsw_distance([ones(idx,1)*HS_box(4,1)...
    SB_sec(1,:)'],[ones(idx,1)*HS_box(4,2) SB_sec(2,:)'])/1000;
NB_secd=gsw_distance([ones(idx,1)*SB_box(4,1)...
    NB_sec(1,:)'],[ones(idx,1)*SB_box(4,2) NB_sec(2,:)'])/1000;

SB_linsec(1,:)=interp1(SB_secd,SB_sec(1,:),0:0.5:max(SB_secd));
SB_linsec(2,:)=interp1(SB_secd,SB_sec(2,:),0:0.5:max(SB_secd));

NB_linsec(1,:)=interp1(NB_secd,NB_sec(1,:),0:0.5:max(NB_secd));
NB_linsec(2,:)=interp1(NB_secd,NB_sec(2,:),0:0.5:max(NB_secd));


SB_IN_section=histcounts2(SB_pass_IN(:,3)*-1,SB_pass_IN(:,5),...
    0:5:450,0:0.5:max(SB_secd));
    
SB_OUT_section=histcounts2(SB_pass_OUT(:,3)*-1,SB_pass_OUT(:,5)...
    ,0:5:450,0:0.5:max(SB_secd));

NB_IN_section=histcounts2(NB_pass_IN(:,3)*-1,NB_pass_IN(:,5),...
    0:5:450,0:0.5:max(NB_secd));

NB_OUT_section=histcounts2(NB_pass_OUT(:,3)*-1,NB_pass_OUT(:,5),...
    0:5:450,0:0.5:max(NB_secd));

IW_count=sum(sum(SB_IN_section(50/5:200/5,:)));
IW_pcnt=IW_count/sum(SB_IN_section(:));

%% 
IN_cm=cmocean('amp');IN_cm(1,:)=[1 1 1];
OUT_cm=cmocean('-ice');OUT_cm(1,:)=[1 1 1];

f1=figure('units','centimeters','outerposition',[0 0 19 16],'color','w');
map_ax=axes('Position',[0.065 0.025 0.35 0.35]);
lon_lim=[-125.5 -122.522345235];
lat_lim=[50.1 48.5];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.25)
hold on
m_gshhs_f('patch',rgb_x('light grey'),'edgecolor','none');
m_northarrow(-125.32,49.75,0.15,'type',2,'linewidth',.5);
m_grid('linestyle','none','tickdir','out','xaxisloc','bottom','yaxisloc','left','fontsize',8,...
    'xtick',-126:-122,'ytick',48:51);
m_scatter(HS_box(3:4,1),HS_box(3:4,2),10,'k','filled');
m_scatter(SB_box(3:4,1),SB_box(3:4,2),10,'k','filled');
m_patch([-124.2723;-123.6640;-124.1356;-124.7425],...
    [49.0594;49.5403;49.7926;49.3092],'m','facealpha',0.1,'edgecolor','k');
m_patch([-123.4397;-122.8289;-123.2230;-123.8327],...
    [48.6058;49.0823;49.3003;48.8217],'m','facealpha',0.1,'edgecolor','k');
text(0.7,0.1,'g)','units','normalized','fontweight','bold','fontsize',8);

slice_ax1=axes('Position',[0.385 0.2 0.45 0.2]);
lon_lim=[-123.6640 -124.7425];
lat_lim=[49.7926 49.0594];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','horiz','aspect',0.3)
hold on
m_gshhs_f('patch',rgb_x('light grey'),'edgecolor','none');
m_grid('linestyle','none','tickdir','out','xaxisloc','bottom','yaxisloc','right','fontsize',8,...
    'xtick',[],'ytick',[]);
m_scatter(HS_box(3:4,1),HS_box(3:4,2),10,'k','filled');
% m_plot(HS_box(3:4,1),HS_box(3:4,2),'k');
% m_plot(SB_box(3:4,1),SB_box(3:4,2),'k');
text(0.9,0.1,'h)','units','normalized','fontweight','bold','fontsize',8);

Bhist3(:,1)=histcounts(NB_pass_IN(:,5),0:0.5:max(NB_secd));
Bhist3(:,2)=histcounts(NB_pass_OUT(:,5),0:0.5:max(NB_secd));
scle=max(Bhist3(:));

for i=1:length(NB_linsec)-1
    if Bhist3(i,1)>0
        [reclon,reclat]=m_ll2xy(NB_linsec(1,i),NB_linsec(2,i));
        [reclon2,~]=m_ll2xy(NB_linsec(1,i+1),NB_linsec(2,i+1));
        rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist3(i,1)/scle)/400],...
            'facecolor',[IN_cm(end/2,:) 0.4],'edgecolor','none');
    end
    
     if Bhist3(i,2)>0
        [reclon,reclat]=m_ll2xy(NB_linsec(1,i),NB_linsec(2,i));
        [reclon2,~]=m_ll2xy(NB_linsec(1,i+1),NB_linsec(2,i+1));
        a=rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist3(i,2)/scle)/400],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        rectangle('Position',[reclon,reclat-(a.Position(4)),reclon2-reclon,a.Position(4)],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        delete(a)
    end
end

m_scatter(nanmean(NB_pass_IN(:,1)),nanmean(NB_pass_IN(:,2)),5,'filled',...
    'markerfacecolor',IN_cm(end/2,:),'markeredgecolor','none');
m_scatter(nanmean(NB_pass_OUT(:,1)),nanmean(NB_pass_OUT(:,2)),5,'filled',...
    'markerfacecolor',OUT_cm(end/2,:),'markeredgecolor','none');

[reclon,reclat]=m_ll2xy(NB_linsec(1,end-1),NB_linsec(2,end-1));
[reclon2,~]=m_ll2xy(NB_linsec(1,end),NB_linsec(2,end));
a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*2,(250/scle)/400],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*2,a.Position(4)],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
% a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*2,(700/scle)/400],...
%     'facecolor',[0 0 0 0.4],'edgecolor','none');
% rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*2,a.Position(4)],...
%     'facecolor',[0 0 0 0.4],'edgecolor','none');

m_text([-124.1263;-124.0548;-123.9797;-123.9120;-123.8408],...
    [49.6771;49.6326;49.5859;49.5435;49.4989],...
    {250;' ';0;' ';250});

slice_ax2=axes('Position',[0.385 0. 0.45 0.2]);
lon_lim=[-122.8289 -123.8327];
lat_lim=[49.3003 48.6058];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','horiz','aspect',0.3)
hold on
m_gshhs_f('patch',rgb_x('light grey'),'edgecolor','none');
m_grid('linestyle','none','tickdir','out','xaxisloc','bottom','yaxisloc','right','fontsize',8,...
    'xtick',[],'ytick',[]);
% m_scatter(HS_box(3:4,1),HS_box(3:4,2),10,'k','filled');
text(0.9,0.1,'i)','units','normalized','fontweight','bold','fontsize',8);

Bhist4(:,1)=histcounts(SB_pass_IN(:,5),0:0.5:max(SB_secd));
Bhist4(:,2)=histcounts(SB_pass_OUT(:,5),0:0.5:max(SB_secd));
scle=max(Bhist3(:));

for i=1:length(SB_linsec)-1
    if Bhist4(i,1)>0

        [reclon,reclat]=m_ll2xy(SB_linsec(1,i),SB_linsec(2,i));
        [reclon2,~]=m_ll2xy(SB_linsec(1,i+1),SB_linsec(2,i+1));
        rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist4(i,1)/scle)/650],...
            'facecolor',[IN_cm(end/2,:) 0.4],'edgecolor','none');
    end
    
     if Bhist4(i,2)>0
        [reclon,reclat]=m_ll2xy(SB_linsec(1,i),SB_linsec(2,i));
        [reclon2,~]=m_ll2xy(SB_linsec(1,i+1),SB_linsec(2,i+1));
        a=rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist4(i,2)/scle)/650],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        rectangle('Position',[reclon,reclat-(a.Position(4)),reclon2-reclon,a.Position(4)],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        delete(a)
    end
end

m_scatter(nanmean(SB_pass_IN(:,1)),nanmean(SB_pass_IN(:,2)),5,'filled',...
    'markerfacecolor',IN_cm(end/2,:),'markeredgecolor','none');
m_scatter(nanmean(SB_pass_OUT(:,1)),nanmean(SB_pass_OUT(:,2)),5,'filled',...
    'markerfacecolor',OUT_cm(end/2,:),'markeredgecolor','none');

[reclon,reclat]=m_ll2xy(SB_linsec(1,end-1),SB_linsec(2,end-1));
[reclon2,~]=m_ll2xy(SB_linsec(1,end),SB_linsec(2,end));
a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*1.5,(250/scle)/650],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*1.5,a.Position(4)],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
% a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*1.5,(900/scle)/650],...
%     'facecolor',[0 0 0 0.4],'edgecolor','none');
% rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*1.5,a.Position(4)],...
%     'facecolor',[0 0 0 0.4],'edgecolor','none');
% a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*1.5,(1350/scle)/650],...
%     'facecolor',[0 0 0 0.1],'edgecolor','none');
% rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*1.5,a.Position(4)],...
%     'facecolor',[0 0 0 0.1],'edgecolor','none');

m_text([-123.2579;-123.2076;-123.1539;-123.0970;-123.0402;-122.9868],...
    [49.1762;49.1448;49.1114;49.0758;49.0403;49.0067],...
    {' ';' ';250;0;250;' '});

ax(1)=axes('Position',[0.225 0.7 0.25 0.15]);
hold on
[X,Y]=meshgrid([0:0.5:max(SB_secd)-0.5]',0:5:450-5);
[CS,h]=contourf(X,Y,SB_IN_section,0:2:max([NB_IN_section(:);...
    SB_IN_section(:)]),'linestyle','none');
caxis([0 max([NB_IN_section(:);SB_IN_section(:)])]);
set(gca,'ydir','reverse','xticklabel',[]);
colormap(ax(1),IN_cm);
area(SB_secd,SB_bathy,400,'FaceColor',rgb_x('light grey'),'edgecolor','none');
ylim([0 200]);
ylabel('Depth (m)');
title('Into Basin');
text(0.9,0.1,'a)','units','normalized','fontweight','bold','fontsize',8);
text(-0.2,0.5,'SB','units','normalized','fontweight','bold','fontsize',8);

ax(2)=axes('Position',[0.5 0.7 0.25 0.15]);
hold on
contourf(X,Y,SB_OUT_section,0:2:max([NB_IN_section(:);...
    SB_OUT_section(:)]),'linestyle','none');
set(gca,'ydir','reverse','xticklabel',[],'yticklabel',[]);
caxis([0 max([NB_OUT_section(:);SB_OUT_section(:)])]);
colormap(ax(2),OUT_cm);
area(SB_secd,SB_bathy,400,'FaceColor',rgb_x('light grey'),'edgecolor','none');
ylim([0 200]);
title('Out of Basin');
text(0.9,0.1,'b)','units','normalized','fontweight','bold','fontsize',8);

ax(3)=axes('Position',[0.225 0.5 0.25 0.15]);
hold on
[X,Y]=meshgrid([0:0.5:max(NB_secd)-0.5]',0:5:450-5);
[CS,h]=contourf(X,Y,NB_IN_section,0:2:max([NB_IN_section(:);...
    SB_IN_section(:)]),'linestyle','none');
set(gca,'ydir','reverse');
caxis([0 max([NB_IN_section(:);SB_IN_section(:)])]);
colormap(ax(3),IN_cm);
cm_ax1=m_contfbar(ax(3),[0 1],-0.425,CS,h,'endpiece','no','axfrac',.1,'fontsize',8);
area(NB_secd,NB_bathy,430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
area([max(NB_secd) max(SB_secd)],[0 0],430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
linkaxes([ax(1) ax(3)],'x');
ylim([0 430]);
xlabel('Across-Strait Distance (km)');
text(0.9,0.1,'d)','units','normalized','fontweight','bold','fontsize',8);
text(-0.2,0.5,'NB','units','normalized','fontweight','bold','fontsize',8);

ax(4)=axes('Position',[0.5 0.5 0.25 0.15]);
hold on
[CS,h]=contourf(X,Y,NB_OUT_section,0:2:max([NB_OUT_section(:);...
    SB_OUT_section(:)]),'linestyle','none');
set(gca,'ydir','reverse','yticklabel',[]);
caxis([0 max([SB_OUT_section(:);NB_OUT_section(:)])]);
colormap(ax(4),OUT_cm);
area(NB_secd,NB_bathy,430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
% xlim([0 max(SB_secd)]);
area([max(NB_secd) max(SB_secd)],[0 0],430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
cm_ax2=m_contfbar(ax(4),[0 1],-0.425,CS,h,'endpiece','no','axfrac',.1,'fontsize',8);
linkaxes([ax(2) ax(4)],'x');
ylim([0 430]);
text(0.9,0.1,'e)','units','normalized','fontweight','bold','fontsize',8);

% Add mini histograms
axh(1)=axes('Position',[0.76 0.7 0.075 0.15]);
hold on
Bhist1(:,1)=histcounts(SB_pass_IN(:,3)*-1,0:5:200);
Bhist1(:,2)=histcounts(SB_pass_OUT(:,3)*-1,0:5:200);
b1=barh(0:5:200-5,Bhist1(:,1),'histc');
set(b1,'facecolor',IN_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
b2=barh(0:5:200-5,Bhist1(:,2),'histc');
set(b2,'facecolor',OUT_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
linkaxes([ax(1) ax(2) axh(1)],'y');
% plot([0 max(xlim)],[nanmean(SB_pass_IN(:,3)*-1) ...
%     nanmean(SB_pass_IN(:,3))*-1],':','color',IN_cm(end/2,:));
% plot([0 max(xlim)],[nanmean(SB_pass_OUT(:,3)*-1) ...
%     nanmean(SB_pass_OUT(:,3))*-1],':','color',OUT_cm(end/2,:));
scatter(0,nanmean(SB_pass_IN(:,3)*-1),15,'filled','markerfacecolor',...
    IN_cm(end/2,:),'markeredgecolor','none');
scatter(0,nanmean(SB_pass_OUT(:,3)*-1),15,'filled','markerfacecolor',...
    OUT_cm(end/2,:),'markeredgecolor','none');
set(gca,'ydir','reverse','yticklabel',[],'xtick',0:250:500,...
    'xticklabel',[]);
xlim([0 500]);
grid on
text(0.9,0.1,'c)','units','normalized','fontweight','bold','fontsize',8);

axh(2)=axes('Position',[0.76 0.5 0.075 0.15]);
hold on
Bhist2(:,1)=histcounts(NB_pass_IN(:,3)*-1,0:5:430);
Bhist2(:,2)=histcounts(NB_pass_OUT(:,3)*-1,0:5:430);
b1=barh(0:5:430-5,Bhist2(:,1),'histc');
set(b1,'facecolor',IN_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
b2=barh(0:5:430-5,Bhist2(:,2),'histc');
set(b2,'facecolor',OUT_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
linkaxes([ax(3) ax(4) axh(2)],'y');
linkaxes([axh(1) axh(2)],'x');
scatter(0,nanmean(NB_pass_IN(:,3)*-1),15,'filled','markerfacecolor',...
    IN_cm(end/2,:),'markeredgecolor','none');
scatter(0,nanmean(NB_pass_OUT(:,3)*-1),15,'filled','markerfacecolor',...
    OUT_cm(end/2,:),'markeredgecolor','none');
grid on
xlabel('Particle Count');
xlim([0 500]);
set(gca,'ydir','reverse','yticklabel',[],'xtick',0:250:500,...
    'xticklabelrotation',45);
text(0.9,0.1,'f)','units','normalized','fontweight','bold','fontsize',8);

set(findall(gcf,'-property','tickdir'),'tickdir','out');
set(findall(gcf,'-property','ticklength'),'ticklength',[0.025 0.02]);
set(findall(gcf,'-property','fontsize'),'fontsize',6);

%%
% export_fig '/ocean/sstevens/IW_project/figures/paper/PT_sections.pdf' -dpdf
export_fig '/ocean/sstevens/IW_project/figures/paper/PT_sections_SUMMER.png' -png -m3

%% Plot stats- WINTER
clear
load /ocean/sstevens/ariane_old/results_BP/PT_sections_V3.mat
%%
SB_pass_IN(SB_pass_IN(:,5)==0,5)=NaN;
NB_pass_IN(NB_pass_IN(:,5)==0,5)=NaN;
SB_pass_OUT(SB_pass_OUT(:,5)==0,5)=NaN;
NB_pass_OUT(NB_pass_OUT(:,5)==0,5)=NaN;

DEC=day(datetime(datestr(datenum(2016,12,1))),'dayofyear');
APR=day(datetime(datestr(datenum(2016,4,1))),'dayofyear');

idx=SB_pass_IN(:,4)>=DEC | SB_pass_IN(:,4)<APR;SB_pass_IN(~idx,:)=NaN;
idx=NB_pass_IN(:,4)>=DEC | NB_pass_IN(:,4)<APR;NB_pass_IN(~idx,:)=NaN;
idx=SB_pass_OUT(:,4)>=DEC | SB_pass_OUT(:,4)<APR;SB_pass_OUT(~idx,:)=NaN;
idx=NB_pass_OUT(:,4)>=DEC | NB_pass_OUT(:,4)<APR;NB_pass_OUT(~idx,:)=NaN;

% Find section bathymetry
[X,Y]=meshgrid(Z_lat,Z_lon);
SB_sec=[linspace(HS_box(4,1),HS_box(3,1),1000);...
    linspace(HS_box(4,2),HS_box(3,2),1000)];
NB_sec=[linspace(SB_box(4,1),SB_box(3,1),1000);...
    linspace(SB_box(4,2),SB_box(3,2),1000)];

SB_bathy=interp2(X,Y,Z,SB_sec(2,:),SB_sec(1,:))*-1;
NB_bathy=interp2(X,Y,Z,NB_sec(2,:),NB_sec(1,:))*-1;
SB_bathy(isnan(SB_bathy))=0;SB_bathy=smooth(SB_bathy,20);
NB_bathy(isnan(NB_bathy))=0;NB_bathy=smooth(NB_bathy,20);

idx=length(SB_sec);
SB_secd=gsw_distance([ones(idx,1)*HS_box(4,1)...
    SB_sec(1,:)'],[ones(idx,1)*HS_box(4,2) SB_sec(2,:)'])/1000;
NB_secd=gsw_distance([ones(idx,1)*SB_box(4,1)...
    NB_sec(1,:)'],[ones(idx,1)*SB_box(4,2) NB_sec(2,:)'])/1000;

SB_linsec(1,:)=interp1(SB_secd,SB_sec(1,:),0:0.5:max(SB_secd));
SB_linsec(2,:)=interp1(SB_secd,SB_sec(2,:),0:0.5:max(SB_secd));

NB_linsec(1,:)=interp1(NB_secd,NB_sec(1,:),0:0.5:max(NB_secd));
NB_linsec(2,:)=interp1(NB_secd,NB_sec(2,:),0:0.5:max(NB_secd));


SB_IN_section=histcounts2(SB_pass_IN(:,3)*-1,SB_pass_IN(:,5),...
    0:5:450,0:0.5:max(SB_secd));
    
SB_OUT_section=histcounts2(SB_pass_OUT(:,3)*-1,SB_pass_OUT(:,5)...
    ,0:5:450,0:0.5:max(SB_secd));

NB_IN_section=histcounts2(NB_pass_IN(:,3)*-1,NB_pass_IN(:,5),...
    0:5:450,0:0.5:max(NB_secd));

NB_OUT_section=histcounts2(NB_pass_OUT(:,3)*-1,NB_pass_OUT(:,5),...
    0:5:450,0:0.5:max(NB_secd));

%% 
IN_cm=cmocean('amp');IN_cm(1,:)=[1 1 1];
OUT_cm=cmocean('-ice');OUT_cm(1,:)=[1 1 1];

f1=figure('units','centimeters','outerposition',[0 0 19 16],'color','w');
map_ax=axes('Position',[0.065 0.025 0.35 0.35]);
lon_lim=[-125.5 -122.522345235];
lat_lim=[50.1 48.5];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.25)
hold on
m_gshhs_f('patch',rgb_x('light grey'),'edgecolor','none');
m_northarrow(-125.32,49.75,0.15,'type',2,'linewidth',.5);
m_grid('linestyle','none','tickdir','out','xaxisloc','bottom','yaxisloc','left','fontsize',8,...
    'xtick',-126:-122,'ytick',48:51);
m_scatter(HS_box(3:4,1),HS_box(3:4,2),10,'k','filled');
m_scatter(SB_box(3:4,1),SB_box(3:4,2),10,'k','filled');
m_patch([-124.2723;-123.6640;-124.1356;-124.7425],...
    [49.0594;49.5403;49.7926;49.3092],'m','facealpha',0.1,'edgecolor','k');
m_patch([-123.4397;-122.8289;-123.2230;-123.8327],...
    [48.6058;49.0823;49.3003;48.8217],'m','facealpha',0.1,'edgecolor','k');
text(0.7,0.1,'g)','units','normalized','fontweight','bold','fontsize',8);

slice_ax1=axes('Position',[0.385 0.2 0.45 0.2]);
lon_lim=[-123.6640 -124.7425];
lat_lim=[49.7926 49.0594];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','horiz','aspect',0.3)
hold on
m_gshhs_f('patch',rgb_x('light grey'),'edgecolor','none');
m_grid('linestyle','none','tickdir','out','xaxisloc','bottom','yaxisloc','right','fontsize',8,...
    'xtick',[],'ytick',[]);
m_scatter(HS_box(3:4,1),HS_box(3:4,2),10,'k','filled');
% m_plot(HS_box(3:4,1),HS_box(3:4,2),'k');
% m_plot(SB_box(3:4,1),SB_box(3:4,2),'k');
text(0.9,0.1,'h)','units','normalized','fontweight','bold','fontsize',8);

Bhist3(:,1)=histcounts(NB_pass_IN(:,5),0:0.5:max(NB_secd));
Bhist3(:,2)=histcounts(NB_pass_OUT(:,5),0:0.5:max(NB_secd));
scle=max(Bhist3(:));

for i=1:length(NB_linsec)-1
    if Bhist3(i,1)>0
        [reclon,reclat]=m_ll2xy(NB_linsec(1,i),NB_linsec(2,i));
        [reclon2,~]=m_ll2xy(NB_linsec(1,i+1),NB_linsec(2,i+1));
        rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist3(i,1)/scle)/400],...
            'facecolor',[IN_cm(end/2,:) 0.4],'edgecolor','none');
    end
    
     if Bhist3(i,2)>0
        [reclon,reclat]=m_ll2xy(NB_linsec(1,i),NB_linsec(2,i));
        [reclon2,~]=m_ll2xy(NB_linsec(1,i+1),NB_linsec(2,i+1));
        a=rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist3(i,2)/scle)/400],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        rectangle('Position',[reclon,reclat-(a.Position(4)),reclon2-reclon,a.Position(4)],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        delete(a)
    end
end

m_scatter(nanmean(NB_pass_IN(:,1)),nanmean(NB_pass_IN(:,2)),5,'filled',...
    'markerfacecolor',IN_cm(end/2,:),'markeredgecolor','none');
m_scatter(nanmean(NB_pass_OUT(:,1)),nanmean(NB_pass_OUT(:,2)),5,'filled',...
    'markerfacecolor',OUT_cm(end/2,:),'markeredgecolor','none');

[reclon,reclat]=m_ll2xy(NB_linsec(1,end-1),NB_linsec(2,end-1));
[reclon2,~]=m_ll2xy(NB_linsec(1,end),NB_linsec(2,end));
a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*2,(250/scle)/400],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*2,a.Position(4)],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
% a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*2,(700/scle)/400],...
%     'facecolor',[0 0 0 0.4],'edgecolor','none');
% rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*2,a.Position(4)],...
%     'facecolor',[0 0 0 0.4],'edgecolor','none');

m_text([-124.1263;-124.0548;-123.9797;-123.9120;-123.8408],...
    [49.6771;49.6326;49.5859;49.5435;49.4989],...
    {' ';250;0;250;' '});

slice_ax2=axes('Position',[0.385 0. 0.45 0.2]);
lon_lim=[-122.8289 -123.8327];
lat_lim=[49.3003 48.6058];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','horiz','aspect',0.3)
hold on
m_gshhs_f('patch',rgb_x('light grey'),'edgecolor','none');
m_grid('linestyle','none','tickdir','out','xaxisloc','bottom','yaxisloc','right','fontsize',8,...
    'xtick',[],'ytick',[]);
% m_scatter(HS_box(3:4,1),HS_box(3:4,2),10,'k','filled');
text(0.9,0.1,'i)','units','normalized','fontweight','bold','fontsize',8);

Bhist4(:,1)=histcounts(SB_pass_IN(:,5),0:0.5:max(SB_secd));
Bhist4(:,2)=histcounts(SB_pass_OUT(:,5),0:0.5:max(SB_secd));
scle=max(Bhist3(:));

for i=1:length(SB_linsec)-1
    if Bhist4(i,1)>0

        [reclon,reclat]=m_ll2xy(SB_linsec(1,i),SB_linsec(2,i));
        [reclon2,~]=m_ll2xy(SB_linsec(1,i+1),SB_linsec(2,i+1));
        rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist4(i,1)/scle)/650],...
            'facecolor',[IN_cm(end/2,:) 0.4],'edgecolor','none');
    end
    
     if Bhist4(i,2)>0
        [reclon,reclat]=m_ll2xy(SB_linsec(1,i),SB_linsec(2,i));
        [reclon2,~]=m_ll2xy(SB_linsec(1,i+1),SB_linsec(2,i+1));
        a=rectangle('Position',[reclon,reclat,reclon2-reclon,(Bhist4(i,2)/scle)/650],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        rectangle('Position',[reclon,reclat-(a.Position(4)),reclon2-reclon,a.Position(4)],...
            'facecolor',[OUT_cm(end/2,:) 0.4],'edgecolor','none');
        delete(a)
    end
end

m_scatter(nanmean(SB_pass_IN(:,1)),nanmean(SB_pass_IN(:,2)),5,'filled',...
    'markerfacecolor',IN_cm(end/2,:),'markeredgecolor','none');
m_scatter(nanmean(SB_pass_OUT(:,1)),nanmean(SB_pass_OUT(:,2)),5,'filled',...
    'markerfacecolor',OUT_cm(end/2,:),'markeredgecolor','none');

[reclon,reclat]=m_ll2xy(SB_linsec(1,end-1),SB_linsec(2,end-1));
[reclon2,~]=m_ll2xy(SB_linsec(1,end),SB_linsec(2,end));
a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*1.5,(250/scle)/650],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*1.5,a.Position(4)],...
    'facecolor',[0 0 0 0.8],'edgecolor','none');
% a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*1.5,(900/scle)/650],...
%     'facecolor',[0 0 0 0.4],'edgecolor','none');
% rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*1.5,a.Position(4)],...
%     'facecolor',[0 0 0 0.4],'edgecolor','none');
% a=rectangle('Position',[reclon,reclat,(reclon2-reclon)*1.5,(1350/scle)/650],...
%     'facecolor',[0 0 0 0.1],'edgecolor','none');
% rectangle('Position',[reclon,reclat-(a.Position(4)),(reclon2-reclon)*1.5,a.Position(4)],...
%     'facecolor',[0 0 0 0.1],'edgecolor','none');

m_text([-123.2579;-123.2076;-123.1539;-123.0970;-123.0402;-122.9868],...
    [49.1762;49.1448;49.1114;49.0758;49.0403;49.0067],...
    {' ';' ';250;0;250;' '});

ax(1)=axes('Position',[0.225 0.7 0.25 0.15]);
hold on
[X,Y]=meshgrid([0:0.5:max(SB_secd)-0.5]',0:5:450-5);
[CS,h]=contourf(X,Y,SB_IN_section,0:2:max([NB_IN_section(:);...
    SB_IN_section(:)]),'linestyle','none');
caxis([0 max([NB_IN_section(:);SB_IN_section(:)])]);
set(gca,'ydir','reverse','xticklabel',[]);
colormap(ax(1),IN_cm);
area(SB_secd,SB_bathy,400,'FaceColor',rgb_x('light grey'),'edgecolor','none');
ylim([0 200]);
ylabel('Depth (m)');
title('Into Basin');
text(0.9,0.1,'a)','units','normalized','fontweight','bold','fontsize',8);
text(-0.2,0.5,'SB','units','normalized','fontweight','bold','fontsize',8);

ax(2)=axes('Position',[0.5 0.7 0.25 0.15]);
hold on
contourf(X,Y,SB_OUT_section,0:2:max([NB_IN_section(:);...
    SB_OUT_section(:)]),'linestyle','none');
set(gca,'ydir','reverse','xticklabel',[],'yticklabel',[]);
caxis([0 max([NB_OUT_section(:);SB_OUT_section(:)])]);
colormap(ax(2),OUT_cm);
area(SB_secd,SB_bathy,400,'FaceColor',rgb_x('light grey'),'edgecolor','none');
ylim([0 200]);
title('Out of Basin');
text(0.9,0.1,'b)','units','normalized','fontweight','bold','fontsize',8);

ax(3)=axes('Position',[0.225 0.5 0.25 0.15]);
hold on
[X,Y]=meshgrid([0:0.5:max(NB_secd)-0.5]',0:5:450-5);
[CS,h]=contourf(X,Y,NB_IN_section,0:2:max([NB_IN_section(:);...
    SB_IN_section(:)]),'linestyle','none');
set(gca,'ydir','reverse');
caxis([0 max([NB_IN_section(:);SB_IN_section(:)])]);
colormap(ax(3),IN_cm);
cm_ax1=m_contfbar(ax(3),[0 1],-0.425,CS,h,'endpiece','no','axfrac',.1,'fontsize',8);
area(NB_secd,NB_bathy,430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
area([max(NB_secd) max(SB_secd)],[0 0],430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
linkaxes([ax(1) ax(3)],'x');
ylim([0 430]);
xlabel('Across-Strait Distance (km)');
text(0.9,0.1,'d)','units','normalized','fontweight','bold','fontsize',8);
text(-0.2,0.5,'NB','units','normalized','fontweight','bold','fontsize',8);

ax(4)=axes('Position',[0.5 0.5 0.25 0.15]);
hold on
[CS,h]=contourf(X,Y,NB_OUT_section,0:2:max([NB_OUT_section(:);...
    SB_OUT_section(:)]),'linestyle','none');
set(gca,'ydir','reverse','yticklabel',[]);
caxis([0 max([SB_OUT_section(:);NB_OUT_section(:)])]);
colormap(ax(4),OUT_cm);
area(NB_secd,NB_bathy,430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
% xlim([0 max(SB_secd)]);
area([max(NB_secd) max(SB_secd)],[0 0],430,'FaceColor',rgb_x('light grey'),'edgecolor','none');
cm_ax2=m_contfbar(ax(4),[0 1],-0.425,CS,h,'endpiece','no','axfrac',.1,'fontsize',8);
linkaxes([ax(2) ax(4)],'x');
ylim([0 430]);
text(0.9,0.1,'e)','units','normalized','fontweight','bold','fontsize',8);

% Add mini histograms
axh(1)=axes('Position',[0.76 0.7 0.075 0.15]);
hold on
Bhist1(:,1)=histcounts(SB_pass_IN(:,3)*-1,0:5:200);
Bhist1(:,2)=histcounts(SB_pass_OUT(:,3)*-1,0:5:200);
b1=barh(0:5:200-5,Bhist1(:,1),'histc');
set(b1,'facecolor',IN_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
b2=barh(0:5:200-5,Bhist1(:,2),'histc');
set(b2,'facecolor',OUT_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
linkaxes([ax(1) ax(2) axh(1)],'y');
set(gca,'ydir','reverse','yticklabel',[],'xtick',0:250:500,...
    'xticklabel',[]);
% plot([0 max(xlim)],[nanmean(SB_pass_IN(:,3)*-1) ...
%     nanmean(SB_pass_IN(:,3))*-1],':','color',IN_cm(end/2,:));
% plot([0 max(xlim)],[nanmean(SB_pass_OUT(:,3)*-1) ...
%     nanmean(SB_pass_OUT(:,3))*-1],':','color',OUT_cm(end/2,:));
scatter(0,nanmean(SB_pass_IN(:,3)*-1),15,'filled','markerfacecolor',...
    IN_cm(end/2,:),'markeredgecolor','none');
scatter(0,nanmean(SB_pass_OUT(:,3)*-1),15,'filled','markerfacecolor',...
    OUT_cm(end/2,:),'markeredgecolor','none');
grid on
xlim([0 500]);
text(0.9,0.1,'c)','units','normalized','fontweight','bold','fontsize',8);

axh(2)=axes('Position',[0.76 0.5 0.075 0.15]);
hold on
Bhist2(:,1)=histcounts(NB_pass_IN(:,3)*-1,0:5:430);
Bhist2(:,2)=histcounts(NB_pass_OUT(:,3)*-1,0:5:430);
b1=barh(0:5:430-5,Bhist2(:,1),'histc');
set(b1,'facecolor',IN_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
b2=barh(0:5:430-5,Bhist2(:,2),'histc');
set(b2,'facecolor',OUT_cm(end/2,:),'facealpha',0.4,'edgecolor','none');
set(gca,'ydir','reverse','yticklabel',[],'xtick',0:250:500,...
    'xticklabelrotation',45);
linkaxes([ax(3) ax(4) axh(2)],'y');
linkaxes([axh(1) axh(2)],'x');
scatter(0,nanmean(NB_pass_IN(:,3)*-1),15,'filled','markerfacecolor',...
    IN_cm(end/2,:),'markeredgecolor','none');
scatter(0,nanmean(NB_pass_OUT(:,3)*-1),15,'filled','markerfacecolor',...
    OUT_cm(end/2,:),'markeredgecolor','none');
grid on
xlim([0 500]);
xlabel('Particle Count');
text(0.9,0.1,'f)','units','normalized','fontweight','bold','fontsize',8);

set(findall(gcf,'-property','tickdir'),'tickdir','out');
set(findall(gcf,'-property','ticklength'),'ticklength',[0.025 0.02]);
set(findall(gcf,'-property','fontsize'),'fontsize',6);

%%
% export_fig '/ocean/sstevens/IW_project/figures/paper/PT_sections.pdf' -dpdf
export_fig '/ocean/sstevens/IW_project/figures/paper/PT_sections_WINTER.png' -png -m3





