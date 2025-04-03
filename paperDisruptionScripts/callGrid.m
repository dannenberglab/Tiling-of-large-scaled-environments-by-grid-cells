clear all
close all


%% load data file
load('20230827_20min_45x45cmHighWalls')

binSize = 1.5;  % define bin size for rate maps (cm)
self.epoch = [-Inf Inf];
self = root;
cells = self.cells;
cel = [4 2]; % tetrode number and cell number
self.cel = cel;
Grid = Analysis_grid(self,cel,binSize);


%% load data file
load('20230827_35min_65x65x45cmHighWalls')
binSize = 2.05;  % define bin size for rate maps (cm)
self.epoch = [-Inf Inf];
self = root;
cells = self.cells;
cel = [4 2]; % tetrode number and cell number
self.cel = cel;
GridScaled = Analysis_grid(self,cel,binSize);

%% calculate spatial crosscorrelation

mat1 = Grid.rate_map;
mat2 = GridScaled.rate_map;
mat1b = mat1(1:end-1,1:end-1);
mat2b = mat2(1:end-1,1:end-1);
out_mat = Cross_Correlation(mat1b,mat2b);
out_mat_grid2=out_mat;
gridCrossCorr2cell = out_mat_grid2(3:end,3:end); %% removing nan of borders 

figure('color','white')
imagesc(gridCrossCorr2cell)
axis xy
xline(size(gridCrossCorr2cell,2)/2, 'Color', 'k', 'LineWidth', 2); % Draw line for Y axis.
yline(size(gridCrossCorr2cell,1)/2, 'Color', 'k', 'LineWidth', 2); % Draw line for X axis.
set(gca,'XTick',[], 'YTick', [])
colormap('parula')

%% calculate distance and angle of the nearest peak to the origin
SPcroscorr = out_mat_grid2;
resulCC_Grid = CC_DistAng2(SPcroscorr);





