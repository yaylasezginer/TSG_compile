% Process TSG underway data 
%
% Step 1: Read in raw data from csv files
% Step 2: Complete data calibrations
% Step 3: Save compiled and cleaned data as .mat file
% Step 4: Plotting
% Steps 1-3 not required after first run
% Requires additional function O2sol.m and m_map package from Rich

clc
clearvars
close all

% %% Step 1: Read in raw data from csv files
% 
% 
% % Change to local directory
% filedir = '/Users/yaylasez/Desktop/Chile/rti_ulagos/UW';
% 
% fn = dir(filedir);
% lvm = contains({fn.name}, '.lvm');
% csv = contains({fn.name}, '.csv');
% datafile = any([lvm', csv'],2);
% fn = fn(datafile);
% 
% data = []; 
% % Open each file in loop
% for i = 1:numel(fn)
%     
%     file = [filedir '/' fn(i).name];
%     opts = detectImportOptions(file, 'FileType','text','NumHeaderLines',22);
%     T = readtable(file,opts);
% 
%     fields = T.Properties.VariableNames;
%     col = true(numel(fields),1);
%     out = [];
%     for j = 1:numel(fields)
%         type = class(T.(fields{j}));
%         
%         % Prevent empty cell arrays of inconsistent size from being added
%         if contains(type,'cell')
%             col(j) = false;   
%         else
%         % Remove outliers
% %         spread = mad(T.(fields{j}));
% %         medio = median(T.(fields{j}));
% %         out(:,j) = abs(T.(fields{j}) - medio) > 3*spread;
%           out(:,j) = isoutlier(T.(fields{j}), 'movmedian',500);
%           % Remove points when TSG drops to 0
%           out(:,j) = any([out(:,j), T.(fields{j}) == 0],2);
%         end
%     end
%     
%     outlier = any(out,2);
%     Tkeep = T(~outlier,col);
%     Tkeep.Properties.VariableNames = [{'Julian_Day'},{'Hour'},{'Minute'},{'Second'}...
%     {'Lat'},{'Long'},{'SST'},{'salinity'},{'Conductivity'},{'O2_sat'},{'O2_umol_L'},{'Optode_T'}...
%     {'Transmisometry'},{'ChlF_mgm3'},{'ChlF_cal'}];
% 
%     data = [data; Tkeep];
%     
% end
% 
% out2 = isoutlier(data.Lat);
% data = data(~out2,:);
% 
% %% Step 2: Complete data calibrations
% % Calibration values from Huinay tests
% 
% data.SST = data.SST * 0.8943505898124280 + 2.1076824587830300;
% data.salinity = data.salinity * 1.0207059862045200 + 0.9528476239631850;
% data.Transmisometry = data.Transmisometry * -0.0014375910015595100 + 4.587121985071780;
% data.ChlF_mgm3 = data.ChlF_mgm3 * 0.0007876896646601650 + 0.10200416997295000;
% o2sol = O2sol(data.salinity,data.SST);
% data.O2_sat = (data.O2_umol_L - o2sol)./o2sol * 100;
% 
% tsg = struct();
% var = Tkeep.Properties.VariableNames;
% for i = 1:numel(var)
%     tsg.(var{i}) = table2array(data(:,i));
% end
% 
% % convert date stamp to UTC - manage offset with uw computer
% 
% [n,v] = size(data);
% uv = ones(n,1);
% tsg.mdate = datenum(2023*uv, 0*uv,data.Julian_Day);
% offset = datenum(0,0,0,0,3,0);
% tsg.mdate = tsg.mdate - offset;
% 
% % Load PAR data and align data with UW data
% load '/Users/yaylasez/Desktop/Chile/rti_ulagos/PAR/PAR.mat'
% tsg.par = PAR.par_clean;
% 
% %% Save data as .mat file
% % Change to local dir!
% save([filedir '/tsg.mat'], 'tsg')

%% Plotting

load '/Users/yaylasez/Desktop/Chile/rti_ulagos/UW/tsg.mat'
load '/Users/yaylasez/Desktop/Chile/rti_ulagos/PAR/PAR.mat'
 
% Set map bounds - 
% m_map path Change to local directory!
addpath('/Users/yaylasez/Documents/MATLAB/m_map')
% Bathymetry datafile location - change to local directory
addpath('/Users/yaylasez/Desktop/Chile')

latrange = [min(tsg.Lat) - 0.05, max(tsg.Lat) + 0.05];
lonrange = [min(tsg.Long) - 0.02, max(tsg.Long) + 0.05];

m_proj('lambert','lat',latrange,'lon',lonrange)
[maplong,maplat] = m_ll2xy(tsg.Long, tsg.Lat);

time_plot = [{'SST'},{'salinity'},{'O2_sat'},...
    {'Transmisometry'},{'ChlF_mgm3'}];


figure
sp = tight_subplot(6,1,[0.02 0.1],[.1 0.01],[0.1 0.01]);

for j = 1:5 %numel(time_plot)
    figure(1)
    sp(j) = subplot(sp(j));
    plot(tsg.mdate, tsg.(time_plot{j}),'k.')
    ylabel(time_plot{j},'interpreter','none')
    set(gca, 'XTickLabel',[])
    
    figure(j+1)
    %color land in dark grey. use previously saved coastline dat for speed
    m_usercoast('PatCoast','patch',[.3 .3 .3]); hold on 
    m_grid('color','k','fontsize',16)  %add grid around plot
    scatter(maplong,maplat, 30, tsg.(time_plot{j}),'filled')
    colormap(jet)
    rgb = colorbar;
    title(time_plot{j},'interpreter','none')
    set(gca, 'FontSize',14)
end

figure(1)
sp(6) = subplot(sp(6));
plot(PAR.mdate, tsg.par,'k.')
ylabel('PAR'); 
linkaxes(sp,'x')
datetick('x')


% Sampling times and locations
load '/Users/yaylasez/Desktop/Chile/rti_ulagos/StationData.mat'
[station_lon, station_lat] = m_ll2xy(meta.Long, meta.Lat);
figure
m_usercoast('PatCoast','patch',[.3 .3 .3]); hold on
m_grid('color','k','fontsize',16)  %add grid around plot
scatter(maplong,maplat, 30, tsg.Julian_Day - tsg.Julian_Day(1),'filled');
colormap(colorscale([18]));
rgb = colorbar;
plot(station_lon, station_lat, 'ko', 'markerfacecol','k')
title('Day','interpreter','none')
set(gca, 'FontSize',14)
hold off


% World context figure
figure
SAlat = [-55.9, 13.8]; SAlong = [-81.12, -33.94];
m_proj('lambert','lat',SAlat,'lon',SAlong)
[maplong,maplat] = m_ll2xy(tsg.Long, tsg.Lat);
m_coast('patch','k'); hold on
m_grid('linest','-', 'xticklabels', [], 'yticklabels', []);
plot([min(maplong), max(maplong), max(maplong), min(maplong), min(maplong)], [min(maplat), min(maplat), max(maplat), max(maplat), min(maplat)], 'r-','linewidth',4);


