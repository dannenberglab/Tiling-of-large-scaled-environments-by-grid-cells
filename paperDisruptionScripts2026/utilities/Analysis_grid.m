

function [resultsG]  = Analysis_grid(self,cel,binSize)

% Parameters
Smooth_kernel = 2;
ratemap_binside = binSize; % define bin size for rate maps (cm)
xcorr_binsize = 5; % ms
velocity_thresh = [-1 -1]; % cm/sec
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
   catch
end

% plot rate map autocorrelogram
fig_rate_map_ac = figure;
try
    self.plot_rate_map_ac(cel, rate_map); % plot rate map autocorrelation
    axis xy
catch
end

try % in case gridness score cannot be calculated
    [gridness_score, props] = self.Gridness(cel, 'binside', ratemap_binside, 'continuize_epochs', 1, 'grid3', 1,'std_smooth_kernel',2);

catch
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

%%  grid score number 2

[HDgridScore]=get_HDGridScore_v2_PowePoint(props)

%% 
AcorrBNT = spatial_autocorrelation(rate_map);
[gscoreAcorBNT, varargoutAcorBNT] = analyses.gridnessScore(AcorrBNT)

if ~isempty(varargoutAcorBNT.spacing) 
      Grid_spacing = sum(varargoutAcorBNT.spacing)/size(varargoutAcorBNT.spacing,1);
else
      Grid_spacing = 0;
end

%%  SpatialInformation

        try
        SI_c = root.SpatialInformation2([tet cell],...
        'binside', ratemap_binside, 'occupation_thresh',...
        0.25,'continuize_epochs', 1,'std_smooth_kernel', 2);
        catch
            SI_c = nan;
        end

        Spatial_Information = SI_c;

%% stability field by two half cros-corr

epoch1 = self.epoch;
epoch2 = epoch1;
        root = self;
        r2Cor = root;
        r2CorB = root;
        % cel = [tet cell]; 

        r2Cor.epoch = [1 round(root.b_ts(end)/2)];
        r2Cor.cel = cel;

        r2CorB.epoch = [round(r2Cor.b_ts(end)/2) r2Cor.b_ts(end)];
        r2CorB.cel = cel;

        S_Stability = r2Cor.CorrelateRateMaps(cel,r2CorB,cel);


%%  results

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
resultsG.S_Stability = S_Stability;
resultsG.Spatial_Information = Spatial_Information;

end

