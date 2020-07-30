% OCEAN /results/SalishSea/ and /results2/
% WEATHER /results/forcing/atmospheric/GEM2.5/operational/ops*
% Indivudual UV and temp files are used to reduce file size
clear
addpath(genpath('/ocean/sstevens/'));

% ncdisp('/results/SalishSea/hindcast.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc');
obs=load('sinfit_temp_data.mat');
% addpath(genpath('/results/SalishSea/hindcast.201812/'));
% addpath(genpath('/results2/SalishSea/nowcast-green.201806/'));

% Load in metadata
model_temp.lat=ncread('/results2/SalishSea/nowcast-green.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc',...
    'nav_lat');
model_temp.lon=ncread('/results2/SalishSea/nowcast-green.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc',...
    'nav_lon');
model_temp.depth=ncread('/results2/SalishSea/nowcast-green.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc',...
    'deptht');

% Replicate metadata in all structs
model_uv.lat=model_temp.lat; model_uv.lon=model_temp.lon;
model_uv.depth=model_temp.depth;
model_ssh.lat=model_temp.lat; model_ssh.lon=model_temp.lon;
model_ssh.depth=model_temp.depth;
model_uv.tdis=gsw_distance([ones(size(model_temp.lon,1),1).*model_temp.lon(1,400) ...
    model_temp.lon(:,400)],...
    [ones(size(model_temp.lon,1),1).*model_temp.lat(1,400) ...
    model_temp.lat(:,400)]);
model_uv.tdis(model_uv.tdis>1e6)=NaN;
model_uv.tcoords=[model_temp.lon(:,400) model_temp.lat(:,400)];

% Find specific depth level to start at
% [~,startloc]=min(abs(model_temp.depth-50));
% [~,endloc]=min(abs(model_temp.depth-100));

loc=find(model_temp.depth>=60 & model_temp.depth<=80);

if length(loc)==1
    step=1;
else
    step=loc(end)-loc(1);
end


% Make list of file dirs
moddir=[];
% moddir=dir('/results/SalishSea/hindcast.201812/*1*');
moddir=[moddir;dir('/results2/SalishSea/nowcast-green.201812/*1*')];
    
% Loop through all files and pull out temperature from desired depth range
model_temp.t=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_ssh.ssh=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_uv.u=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_uv.v=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_uv.v_trans=NaN(size(model_temp.lon,1),size(model_temp.depth,1),length(moddir));
model_temp.mtime=NaN(length(moddir));

cd('/results2/SalishSea/nowcast-green.201812/')
h = waitbar(0,'Parsing SSC model data...');
tic
errdat=[];


for i=1:length(moddir)
    try
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_T.nc']);
        model_temp.t(:,:,i)=mean(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'votemper',...
            [1 1 loc(1) 1],...
            [Inf Inf step Inf]),3);
        keyboard
        model_temp.mtime(i)=datenum('1900-01-01 00:00:00')+...
            ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'time_counter')/86400;
        model_ssh.ssh(:,:,i)=ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'sossheig');
        
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_U.nc']);
        model_uv.u(:,:,i)=mean(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'vozocrtx',...
            [1 1 loc(1) 1],...
            [Inf Inf step Inf]),3);
        
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_V.nc']);
        model_uv.v(:,:,i)=mean(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'vomecrty',...
            [1 1 loc(1) 1],...
            [Inf Inf step Inf]),3);
        
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_V.nc']);
        model_uv.v_trans(:,:,i)=squeeze(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'vomecrty',...
            [1 400 1 1],...
            [Inf 1 Inf Inf]));
        
    catch
        model_uv.v(:,:,i)=NaN;
        model_uv.v(:,:,i)=NaN;
        errdat=[errdat;model_temp.mtime(i)];
    end
    waitbar(i/length(moddir),h);
end

close
toc

[model_temp.mtime,idx]=sort(model_temp.mtime(:,1));
model_temp.t=model_temp.t(:,:,idx);
model_uv.u=model_uv.u(:,:,idx);
model_uv.v=model_uv.v(:,:,idx);

model_uv.mtime=model_temp.mtime;
model_ssh.mtime=model_temp.mtime;

clearvars -except model_uv model_temp model_ssh
save('/ocean/sstevens/IW_project/data/model_IW_outputT_surf.mat','model_temp','-v7.3');
save('/ocean/sstevens/IW_project/data/model_IW_outputUV_surf.mat','model_uv','-v7.3');
save('/ocean/sstevens/IW_project/data/model_IW_outputssh_surf.mat','model_ssh','-v7.3');

%%
clear
addpath(genpath('/ocean/sstevens/'));

% ncdisp('/results/SalishSea/hindcast.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc');
obs=load('sinfit_temp_data.mat');
% addpath(genpath('/results/SalishSea/hindcast.201812/'));
% addpath(genpath('/results2/SalishSea/nowcast-green.201806/'));

% Load in metadata
model_temp.lat=ncread('/results2/SalishSea/nowcast-green.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc',...
    'nav_lat');
model_temp.lon=ncread('/results2/SalishSea/nowcast-green.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc',...
    'nav_lon');
model_temp.depth=ncread('/results2/SalishSea/nowcast-green.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc',...
    'deptht');

% Replicate metadata in all structs
model_uv.lat=model_temp.lat; model_uv.lon=model_temp.lon;
model_uv.depth=model_temp.depth;
model_ssh.lat=model_temp.lat; model_ssh.lon=model_temp.lon;
model_ssh.depth=model_temp.depth;
model_uv.tdis=gsw_distance([ones(size(model_temp.lon,1),1).*model_temp.lon(1,400) ...
    model_temp.lon(:,400)],...
    [ones(size(model_temp.lon,1),1).*model_temp.lat(1,400) ...
    model_temp.lat(:,400)]);
model_uv.tdis(model_uv.tdis>1e6)=NaN;
model_uv.tcoords=[model_temp.lon(:,400) model_temp.lat(:,400)];

% Find specific depth level to start at
% [~,startloc]=min(abs(model_temp.depth-50));
% [~,endloc]=min(abs(model_temp.depth-100));

loc=find(model_temp.depth>=90 & model_temp.depth<=110);

if length(loc)==1
    step=1;
else
    step=loc(end)-loc(1);
end


% Make list of file dirs
moddir=[];
% moddir=dir('/results/SalishSea/hindcast.201812/*1*');
moddir=[moddir;dir('/results2/SalishSea/nowcast-green.201812/*1*')];
    
% Loop through all files and pull out temperature from desired depth range
model_temp.t=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_ssh.ssh=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_uv.u=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_uv.v=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_uv.v_trans=NaN(size(model_temp.lon,1),size(model_temp.depth,1),length(moddir));
model_temp.mtime=NaN(length(moddir));

cd('/results2/SalishSea/nowcast-green.201812/')
h = waitbar(0,'Parsing SSC model data...');
tic
errdat=[];


for i=1:length(moddir)
    try
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_T.nc']);
        model_temp.t(:,:,i)=mean(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'votemper',...
            [1 1 loc(1) 1],...
            [Inf Inf step Inf]),3);
        
        model_temp.mtime(i)=datenum('1900-01-01 00:00:00')+...
            ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'time_counter')/86400;
        model_ssh.ssh(:,:,i)=ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'sossheig');
        
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_U.nc']);
        model_uv.u(:,:,i)=mean(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'vozocrtx',...
            [1 1 loc(1) 1],...
            [Inf Inf step Inf]),3);
        
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_V.nc']);
        model_uv.v(:,:,i)=mean(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'vomecrty',...
            [1 1 loc(1) 1],...
            [Inf Inf step Inf]),3);
        
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_V.nc']);
        model_uv.v_trans(:,:,i)=squeeze(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'vomecrty',...
            [1 400 1 1],...
            [Inf 1 Inf Inf]));
        
    catch
        model_uv.v(:,:,i)=NaN;
        model_uv.v(:,:,i)=NaN;
        errdat=[errdat;model_temp.mtime(i)];
    end
    waitbar(i/length(moddir),h);
end

close
toc

[model_temp.mtime,idx]=sort(model_temp.mtime(:,1));
model_temp.t=model_temp.t(:,:,idx);
model_uv.u=model_uv.u(:,:,idx);
model_uv.v=model_uv.v(:,:,idx);

model_uv.mtime=model_temp.mtime;
model_ssh.mtime=model_temp.mtime;

clearvars -except model_uv model_temp model_ssh
save('/ocean/sstevens/IW_project/data/model_IW_outputT_core.mat','model_temp','-v7.3');
save('/ocean/sstevens/IW_project/data/model_IW_outputUV_core.mat','model_uv','-v7.3');
save('/ocean/sstevens/IW_project/data/model_IW_outputssh_core.mat','model_ssh','-v7.3');

%% 
clear
addpath(genpath('/ocean/sstevens/'));

% ncdisp('/results/SalishSea/hindcast.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc');
obs=load('sinfit_temp_data.mat');
% addpath(genpath('/results/SalishSea/hindcast.201812/'));
% addpath(genpath('/results2/SalishSea/nowcast-green.201806/'));

% Load in metadata
model_temp.lat=ncread('/results2/SalishSea/nowcast-green.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc',...
    'nav_lat');
model_temp.lon=ncread('/results2/SalishSea/nowcast-green.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc',...
    'nav_lon');
model_temp.depth=ncread('/results2/SalishSea/nowcast-green.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc',...
    'deptht');

% Replicate metadata in all structs
model_uv.lat=model_temp.lat; model_uv.lon=model_temp.lon;
model_uv.depth=model_temp.depth;
model_ssh.lat=model_temp.lat; model_ssh.lon=model_temp.lon;
model_ssh.depth=model_temp.depth;
model_uv.tdis=gsw_distance([ones(size(model_temp.lon,1),1).*model_temp.lon(1,400) ...
    model_temp.lon(:,400)],...
    [ones(size(model_temp.lon,1),1).*model_temp.lat(1,400) ...
    model_temp.lat(:,400)]);
model_uv.tdis(model_uv.tdis>1e6)=NaN;
model_uv.tcoords=[model_temp.lon(:,400) model_temp.lat(:,400)];

% Find specific depth level to start at
% [~,startloc]=min(abs(model_temp.depth-50));
% [~,endloc]=min(abs(model_temp.depth-100));

loc=find(model_temp.depth>=120 & model_temp.depth<=140);

if length(loc)==1
    step=1;
else
    step=loc(end)-loc(1);
end


% Make list of file dirs
moddir=[];
% moddir=dir('/results/SalishSea/hindcast.201812/*1*');
moddir=[moddir;dir('/results2/SalishSea/nowcast-green.201812/*1*')];
    
% Loop through all files and pull out temperature from desired depth range
model_temp.t=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_ssh.ssh=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_uv.u=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_uv.v=NaN(size(model_temp.lon,1),size(model_temp.lon,2),length(moddir));
model_uv.v_trans=NaN(size(model_temp.lon,1),size(model_temp.depth,1),length(moddir));
model_temp.mtime=NaN(length(moddir));

cd('/results2/SalishSea/nowcast-green.201812/')
h = waitbar(0,'Parsing SSC model data...');
tic
errdat=[];


for i=1:length(moddir)
    try
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_T.nc']);
        model_temp.t(:,:,i)=mean(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'votemper',...
            [1 1 loc(1) 1],...
            [Inf Inf step Inf]),3);
        
        model_temp.mtime(i)=datenum('1900-01-01 00:00:00')+...
            ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'time_counter')/86400;
        model_ssh.ssh(:,:,i)=ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'sossheig');
        
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_U.nc']);
        model_uv.u(:,:,i)=mean(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'vozocrtx',...
            [1 1 loc(1) 1],...
            [Inf Inf step Inf]),3);
        
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_V.nc']);
        model_uv.v(:,:,i)=mean(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'vomecrty',...
            [1 1 loc(1) 1],...
            [Inf Inf step Inf]),3);
        
        file=dir([moddir(i).folder,'/',moddir(i).name,'/SalishSea_1d*grid_V.nc']);
        model_uv.v_trans(:,:,i)=squeeze(ncread([moddir(i).folder,'/',moddir(i).name,'/',file.name],'vomecrty',...
            [1 400 1 1],...
            [Inf 1 Inf Inf]));
        
    catch
        model_uv.v(:,:,i)=NaN;
        model_uv.v(:,:,i)=NaN;
        errdat=[errdat;model_temp.mtime(i)];
    end
    waitbar(i/length(moddir),h);
end

close
toc

[model_temp.mtime,idx]=sort(model_temp.mtime(:,1));
model_temp.t=model_temp.t(:,:,idx);
model_uv.u=model_uv.u(:,:,idx);
model_uv.v=model_uv.v(:,:,idx);

model_uv.mtime=model_temp.mtime;
model_ssh.mtime=model_temp.mtime;

clearvars -except model_uv model_temp model_ssh
save('/ocean/sstevens/IW_project/data/model_IW_outputT_deep.mat','model_temp','-v7.3');
save('/ocean/sstevens/IW_project/data/model_IW_outputUV_deep.mat','model_uv','-v7.3');
save('/ocean/sstevens/IW_project/data/model_IW_outputssh_deep.mat','model_ssh','-v7.3');