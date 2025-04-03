

function [resultsG]  = Analysis_grid(self,cel,binSize)

% Parameters

Smooth_kernel = 2;
ratemap_binside = binSize; % define bin size for rate maps (cm)
xcorr_binsize = 5; % ms
velocity_thresh = [-1 -1]; % cm/sec
% plot_vec = zeros(25, 1); % which plots are active
% plot_vec(14) = 1; % header information is active
theta_skip = 0; % run theta skipping analysis with spike time autocorr



% set self.cel
self.cel = cel;
% get tetrode number and cell number
tet = cel(1,1);
cell = cel(1,2);

% get session and mouse ID
[~,session_name,~] = fileparts(self.name); % get session name
[~,mouse_ID_name,~] = fileparts(pwd); % get mouse ID
session = strrep(session_name,'_',' ');
mouse_ID = strrep(mouse_ID_name,'_',' ');
% laser ON or OFF session


% plot trajectory
fig_trajectory = figure;
self.plot_trajectory(cel);

% plot rate map
fig_rate_map = figure;
try
    rate_map = self.RateMap(cel, 'continuize_epochs', 1, 'supress_plot', 0, 'binside', ratemap_binside,'std_smooth_kernel', 3); % get rate map
    %rate_map = self.RateMap(cel, 'continuize_epochs', 1, 'supress_plot',
    %0, 'binside', ratemap_binside); % get rate map Original

catch
end

% plot rate map autocorrelogram
fig_rate_map_ac = figure;
try
    self.plot_rate_map_ac(cel, rate_map); % plot rate map autocorrelation
    axis xy
catch
end


AcorrBNT = spatial_autocorrelation(rate_map);

try % in case gridness score cannot be calculated
    [gridness_score, props] = self.Gridness(cel, 'binside', ratemap_binside, 'continuize_epochs', 1, 'grid3', 1,'std_smooth_kernel',2);

catch
end


%% 

[gscoreAcorBNT, varargoutAcorBNT] = analyses.gridnessScore(AcorrBNT)

if ~isempty(varargoutAcorBNT.spacing) 
      Grid_spacing = sum(varargoutAcorBNT.spacing)/size(varargoutAcorBNT.spacing,1);
else
      Grid_spacing = 0;
end


 
%% plot gridness score 3
fig_grid3 = figure;
if exist('gridness_score','var') && ~isempty(gridness_score) % if analysis valid
    line(props.periodicity(:,1), props.periodicity(:,2), 'Color', 'k', 'LineWidth', 1.5), hold on;
    text(85, .92, {['Gridness3: ' num2str(gridness_score)];['BNTAc: ' num2str(gscoreAcorBNT)];['GridSpaci: ' num2str(Grid_spacing)]}, 'FontSize', 10);
    title('Periodicity of Correlation of Rotated AutoCorr');
    xlabel('Rotation Angle'); ylabel('Correlation');
    xlim([0 180]);
    ylim([-.5, 1.1]);
else
    axis off
    text(.1, .4, 'Invalid Gridness: No Peaks in Autocorrelogram');
end
axis square
ax = gca; % get axis object
ax.FontSize = 18; % change FontSize of axes and title
ax.LineWidth = 1.5; % change line width of axes outline, tick marks, and grid lines
ax.TickLength = [0.03 0.1]; % change width and length of tick marks
ax.XMinorTick = 'on'; % turn on minor tick marks for x axis
ax.YMinorTick = 'on'; % turn on minor tick marks for y axis
ax.TitleFontWeight = 'normal'; % change font of title from 'bold' to 'normal'
ax.TitleFontSizeMultiplier = 0.8;

% plot waveform
fig_waveform = figure;
try % try-catch, because waveform currently doesn't work for merged sessions
    ind = find(sum(repmat(self.cel,size(self.cells,1),1) == self.cells,2) == 2);
    m = {self.user_def.waveform(ind,:).mean};
    s = {self.user_def.waveform(ind,:).std};
    m=cellfun(@(x) x(:)', m,'UniformOutput',0);
    s=cellfun(@(x) x(:)', s,'UniformOutput',0);
    hold on
    for i = 1:length(m)
        t = ((i-1)*length(m{i})+11:(i)*length(m{i})+10) + (i-1)*5;
        t = [t fliplr(t)];
        patch(t,[m{i}+s{i} fliplr(m{i}-s{i})],[.8 .8 .8],'EdgeColor',[.8 .8 .8]);
        t = ((i-1)*length(m{i})+11:(i)*length(m{i})+10) + (i-1)*5;
        plot(t,m{i},'r','LineWidth',3)
    end
    title('Waveforms (1ms each)')
    ylim([min(cellfun(@(x) min(x),m)) max(cellfun(@(x) max(x),m))])
    xlim([0 t(end)+5])
    ylabel('Voltage (µV)')
    ax = gca; % get handle for graphics axes
    ax.XTickLabel = []; % turn of x labels
    ax.XTick = []; % turn of X ticks
    ax.FontSize = 18; % change FontSize of axes and title
    ax.LineWidth = 1.5; % change line width of axes outline, tick marks, and grid lines
    ax.TickLength = [0.03 0.1]; % change width and length of tick marks
    ax.XMinorTick = 'on'; % turn on minor tick marks for x axis
    ax.YMinorTick = 'on'; % turn on minor tick marks for y axis
    ax.TitleFontWeight = 'normal'; % change font of title from 'bold' to 'normal'
    ax.TitleFontSizeMultiplier = 0.8;
catch
end

%%



% plot spiking rate across session
try % try loop to deal with root.epoch being either a cell (meaning the session has been epoched) or a matrix
    if strcmp(laser, 'none')%self.epoch(1) == self.ts(1) && self.epoch(2) == self.ts(end)
        fig_firingRate = figure;
        load(session_name,'Inh','Inh_total','lightOFF','laserON','manipulation')
        if exist('Inh','var') % for non-merged sessions
            if isa(Inh,'struct')
                Inh_total = Inh.S1;
            end
            [stimts] = Inh_total(:,1)-30;
            stop_stimts = stimts + 145;
        end
        if exist('lightOFF','var') % for lightVSdarkness sessions
            stimts = lightOFF(:,1);
            stop_stimts = lightOFF(:,2);
        end
        if exist('laserON','var') % for rhythmic laser stimulation epochs
            stimts = laserON(:,1);
            stop_stimts = laserON(:,2);
        end
        if exist('manipulation','var') % for baselineVSmanipulation sessions
            stimts = manipulation(:,1);
            stop_stimts = manipulation(:,2);
        end
        if stop_stimts(end) > self.ts(end)
            stop_stimts(end) = [];
        end
        All_ts = vertcat(stimts,stimts,stop_stimts,stop_stimts);
        All_ts = sort(All_ts);
        s = zeros(length(All_ts),1);
        s(1:4:end) = 0;
        s(2:4:end) = 1.3;
        s(3:4:end) = 1.3;
        s(4:4:end) = 0;
        All_ts = [All_ts;All_ts(end)];
        s = [s;0];
        % smoothed histogram figure
        smooth_factor = 10; % smoothing in s
        hold on
        edges = 0:1:self.cel_ts{1}(end);
        [N,edges] = histcounts(self.cel_ts{1},edges);
        edges_corrected = edges(1:end-1)+diff(edges);
        N_smooth = smooth(N,smooth_factor);
        fill(All_ts,s./1.3.*max(N_smooth),'g')
        alpha(0.5)
        plot(edges_corrected,N_smooth,'k')
        ylabel('Frequency (Hz)')
        xlabel('Time (s)')
        title(strcat(session,sprintf(', T%dC%d',tet,cell)))
        
        % average for non-inhibition and inhibition epochs
        boundaries = vertcat(stimts,stop_stimts);
        boundaries = sort(boundaries);
        boundaries = round(boundaries);
        avg = nan(length(boundaries),1);
        for b = 1:length(boundaries)
            if b == 1
                avg(b) = mean(N(1:boundaries(b)));
            else
                try
                    avg(b) = mean(N(boundaries(b-1):boundaries(b)));
                catch
                    avg(b) = mean(N(boundaries(b-1):end));
                end
            end
        end
        x = boundaries(1:end-1) + round(diff(boundaries)/2);
        x = vertcat(round(boundaries(1)/2),x);
        plot(x,avg,'linewidth',2,'Color','b')
        hold off
    end

end

resultsG.mouse_ID = mouse_ID;
resultsG.session = session;
resultsG.BinSize_smooth = ['binSize' num2str(ratemap_binside)  ' Smooth' num2str(Smooth_kernel)];
resultsG.tetrode = ['t' num2str(tet)  ' cel_' num2str(cell)];
resultsG.cell = cel;
if exist('gridness_score','var') && ~isempty(gridness_score)
    resultsG.gridness3 = gridness_score;
else
    resultsG.gridness3 = NaN;
end

resultsG.rate_map=props.rate_map;
resultsG.gscoreAcorBNT= gscoreAcorBNT;
resultsG.Grid_spacing = Grid_spacing;

end

