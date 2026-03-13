function [correlation_matrix,correlation_matrixRaw] = spatial_autocorrelation(map, varargin)
%SPATIAL_AUTOCORRELATION Calculate 2D autocorrelation.
%   Calculate 2D autocorrelation and remove edge artifacts from results. We
%   perform a correlation of the map with itself, and define an outer
%   boundry which we do not return. To return the entire autocorrelation
%   map set `output_region` to `1`.
%
%  Args:
%                   map: A 2D ratemap.
%         output_region: Fraction of map to keep.
%               inpaint: If true, inpaints NaN values in ratemap.
%
%  Returns:
%    correlation_matrix: A 2D autocorrelation matrix of the map.
%
%  Dependencies:
%    `normxcorr2_masked()` [1].
%
%  Example of use:
%    ```
%    ratemap = randn(151);
%    acorr = analyse.spatial_autocorrelation(ratemap);
%    figure;
%    imagesc(acorr);
%    ```
%
%  Ref:
%    [1] "Masked object registration in the Fourier domain"
%         Padfield D, 2012, IEEE Transactions on Image Processing
%         https://doi.org/10.1109/TIP.2011.2181402
%
%
% 2018-06-10 frolov, BNT.
% 2020-06-25 waaga, new correlation function.

%%   taken from DOI 10.5281/zenodo.6181066.


ip = inputParser();
ip.addRequired('map', @(x) validateattributes(x, ...
    {'numeric'}, {'real', 'nonsparse', '2d'}));
ip.addParameter('output_region', .9);
ip.addParameter('inpaint', false);
ip.addParameter('mask_edge', true);
ip.parse(map, varargin{:});
p = ip.Results;

% Fraction of output region: Edge-areas produce artefacts, so we only
% return the center of the map. Should be a value in range 0..1.
output_region = p.output_region;

% handle nan values; we mask the nan values before correlating.
nan_idx = isnan(map);
if p.inpaint
    map = external.inpaintn(map);
    
end
mask = ~nan_idx;
map(nan_idx) = 0;

new_size_v = round((size(map, 1) + size(map, 1)) * output_region);
new_size_h = round((size(map, 2) + size(map, 2)) * output_region);

if mod(new_size_v, 2) == 0 && new_size_v > 0
    new_size_v = new_size_v - 1;
end

if mod(new_size_h, 2) == 0 && new_size_h > 0
    new_size_h = new_size_h - 1;
end

% Calculate correlation matrix [1].
correlation_matrix = normxcorr2_masked(map, map, mask, mask);
[correlation_matrixRaw] = correlation_matrix;
[n_row, n_col] = size(correlation_matrix);

offset_v = round((n_row - new_size_v) / 2 + 1);
offset_h = round((n_col - new_size_h) / 2 + 1);

if p.mask_edge
    cm_size = size(correlation_matrix);
    [hi_col, hi_row] = meshgrid(1:cm_size(1), 1:cm_size(2));
    radius = (cm_size(1)-offset_v-offset_h+1) / 2;
    
%     circle_mask = (hi_row - offset_v - radius).^2 + ...
%         (hi_col - offset_h - radius).^2 <= radius.^2;

    circle_mask = (hi_row - offset_v - radius).^2 + ...
        (hi_col - offset_h - radius).^2 <= radius.^2;  

    circle_mask2= circle_mask'; %% added jjhp 2/19/2024

%     correlation_matrix3 = correlation_matrix;
%     correlation_matrix(~circle_mask) = 0; %% original
    correlation_matrix(~circle_mask2) = 0;
%     correlation_matrix2(~circle_mask2) = 0;

end

correlation_matrix = correlation_matrix( ...
    offset_v:end-offset_v+1, ...
    offset_h:end-offset_h+1);

% correlation_matrix = correlation_matrix3( ...
%     offset_v:end-offset_v+1, ...
%     offset_h:end-offset_h+1);

end

