function spatial_information = SpatialInformation2(self, cel, varargin)
% spatial_information = root.SpatialInformation(cel);
%
% Computes the spatial information of a cell
%
% Will return the spatial information score in bits/spike for
% tetrode cel(1), cell cel(2). Can return continuized epochs, or vectorizes 
% multiple epochs. Note that if xdim and ydim are not passed, the common xdim 
% and y dim between epochs if selected, which may affect information score. 
% For example, without xdim and ydim, running this function with epochs A=>SIa1 and B=SIb1, 
% and then with just A=>SIa2, SIa1 is not necessarily the same as SIa2. 
%
% See parameters below.
%
% See Cacucci et al 2007 Methods
%
% OPTIONAL PARAMETERS
%   
%   xdim                vector of bin edges along x dimension
%   ydim                vector of bin edges along y dimension
%   continuize_epochs   0 or 1 (0). If 0, ratemap is calculated for each
%                       epoch (adds 3rd dim to rate_map output, if 1, ratemap
%                       is calculated across all epochs
%   std_smooth_kernel   STD of the gaussian kernel to smooth the rate map  
%   binside             The length in cm of the side of a bin when
%                       calculating the rate map
%   occupation_thresh   Bins that had less than this number of seconds of
%                       occupation are not included in the score. (0)
%   ratemap             1x2 Cell array with a pre-computed ratemap in the
%                       first cell and the occupancy in the second cell.
% andrew 14 mat 2010
% enewman 20131004 added pre-computed ratemap parameter
% wchapman 20140514 recalculated 'F' by overall firing rate

p = inputParser;

p.addRequired('self')
p.addRequired('cel', @isnumeric)
p.addParamValue('xdim', [], @isnumeric);
p.addParamValue('ydim', [], @isnumeric);
p.addParamValue('continuize_epochs', 0, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('std_smooth_kernel', 0, @isnumeric);
p.addParamValue('binside', 3, @isnumeric);
p.addParamValue('occupation_thresh', 0, @isnumeric);
p.addParamValue('ratemap', [], @iscell)

p.parse(self, cel, varargin{:});

self = p.Results.self;
cel = p.Results.cel;
xdim = p.Results.xdim;
ydim = p.Results.ydim;
continuize_epochs = p.Results.continuize_epochs;
std_smooth_kernel = p.Results.std_smooth_kernel;
binside = p.Results.binside;
occupation_thresh = p.Results.occupation_thresh;
ratemap = p.Results.ratemap;

self.cel = cel;
n_thresh = 50; % spikes... how many before information score is meaningless

% if isempty(ratemap)
%   [rate_map, ~, ~, occupancy] = self.RateMap('cel',cel, 'continuize_epochs', continuize_epochs, 'std_smooth_kernel', std_smooth_kernel,'supress_plot',0, 'binside', binside, 'xdim', xdim, 'ydim', ydim);
% else
if isempty(ratemap)
  [rate_map, ~, ~, occupancy] = self.RateMap('cel',cel, 'continuize_epochs', continuize_epochs, 'std_smooth_kernel', std_smooth_kernel, 'binside', binside, 'xdim', xdim, 'ydim', ydim);
else

  rate_map = ratemap{1};
  occupancy = ratemap{2};
end

n_spikes = self.cel_ts;

if iscell(n_spikes)
    n_spikes = cellfun(@length, n_spikes);
    
    if continuize_epochs, n_spikes = sum(n_spikes); end
    
else
    n_spikes = length(n_spikes);
end

% 1) threshold in seconds (raw occupancy)

% ####  update 20251227 with occupation_thresh
% Occupancy threshold (seconds) 
bad = isfinite(occupancy) & (occupancy <= occupation_thresh);

rate_map(bad) = NaN;
occupancy(bad) = NaN;

% Normalize occupancy to probability (per epoch slice)
for i = 1:size(occupancy, 3)
    occSum = nansum(occupancy(:,:,i), 'all');
    if occSum > 0
        p_occupancy(:,:,i) = occupancy(:,:,i) ./ occSum;
    else
        p_occupancy(:,:,i) = NaN; % nothing valid
    end
end

% Robustness: avoid log2(0) 
rate_map(rate_map <= 0) = NaN;

%  ##########

% for i = 1:size(occupancy, 3)
%     occupancy(:, :, i) = occupancy(:, :, i) / sum(sum(occupancy(:, :, i))); % normalize to probability
% end
   
%rate_map(occupancy <= occupation_thresh) = NaN;

%%F = nanmean(rate_map(:));

F = (cellfun(@(x) numel(x),self.cel_i)) ./ diff(self.epoch,1); %firing rate

% ####  update 20251227
F = F(:)';          % ensure 1 x K
% #####

spatial_information = nansum(nansum( p_occupancy .* rate_map .* log2(rate_map ./ repmat(F, size(rate_map, 1), size(rate_map, 2)) ), 1),2) ./ F;

spatial_information = spatial_information(:);

spatial_information(n_spikes<n_thresh) = NaN; % where spiking was too low, get rid of spatial information score

end




