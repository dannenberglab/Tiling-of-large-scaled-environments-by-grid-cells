function [ correlation_matrix, correlation_matrix_Square, n_overlap_masked_pixel ] = ...
    normxcorr2_masked( fixed_image, moving_image, fixed_mask, moving_mask )
%NORMXCORR2_MASKED Masked normalized 2D cross-correlation
%   This function calculates the Masked NCC using FFTs instead of 
%   spatial correlation.  It is therefore much faster for larger
%   structuring elements.
%
%   The masked normalized cross-correlation of fixed_image and
%   moving_image are calculated using masks fixed_mask and moving_mask.
%
%   fixed_mask and moving_mask should consist of only 1s and 0s, where 1
%   indicates locations of useful information in the corresponding image,
%   and 0 indicates locations that should be masked (ignored).
%
%   fixed_image and moving_image need not be the same size, but fixed_mask
%   must be the same size as fixed_image, and moving_mask must be the same
%   size as moving_image.
%
%   The resulting matrix map contains correlation coefficients and its
%   values range from -1.0 to 1.0.
%
%   If you use this code in your work, please reference the following:
%   D. Padfield. "Masked Object Registration in the Fourier Domain".
%   Transactions on Image Processing.
%   D. Padfield. "Masked FFT registration". In Proc. Computer Vision and
%   Pattern Recognition, 2010.
%
%   Author: Dirk Padfield, GE Global Research, padfield@research.ge.com

ip = inputParser();
ip.addRequired('fixed_image', @(x) validateattributes(x, ...
    {'numeric'}, {'real','nonempty', 'nonsparse', '2d', 'finite'}));
ip.addRequired('moving_image', @(x) validateattributes(x, ...
    {'numeric'}, {'real','nonempty', 'nonsparse', '2d', 'finite'}));
ip.addRequired('fixed_mask', @(x) validateattributes(x, ...
    {'logical','numeric'}, {'real', 'nonsparse', '2d', 'finite'}));
ip.addRequired('moving_mask', @(x) validateattributes(x, ...
    {'logical','numeric'}, {'real', 'nonsparse', '2d', 'finite'}));
ip.parse(fixed_image, moving_image, fixed_mask, moving_mask);
p = ip.Results;

% If either fixed_image or moving_image has a minimum value which is 
% negative, we need to shift the array so all values are positive to 
% ensure numerically robust results for the normalized cross-correlation.
fixed_image  = shift_data(p.fixed_image);
moving_image = shift_data(p.moving_image);
fixed_mask   = double(p.fixed_mask);
moving_mask  = double(p.moving_mask);

% Ensure that the masks consist of only 0s and 1s.  Anything less than or
% equal to 0 gets set to 0, and everything else gets set to 1.
fixed_mask(  fixed_mask <= 0) = 0;
fixed_mask(  fixed_mask >  0) = 1;
moving_mask(moving_mask <= 0) = 0;
moving_mask(moving_mask >  0) = 1;

% The fixed and moving images need to be masked for the equations below to
% work correctly.
fixed_image(  fixed_mask == 0) = 0;
moving_image(moving_mask == 0) = 0;

% Flip the moving image and mask in both dimensions so that its correlation
% can be more easily handled.
rotated_moving_image = rot90(moving_image, 2);
rotated_moving_mask  = rot90(moving_mask,  2);

% Calculate all of the FFTs that will be needed.
fixed_image_size  = size(fixed_image);
moving_image_size = size(rotated_moving_image);
combined_size     = fixed_image_size + moving_image_size - 1;

% Find the next largest size that is a multiple of a combination of 2, 3,
% and/or 5.  This makes the FFT calculation much faster.
optimal_size(1) = find_closest_valid_dimension(combined_size(1));
optimal_size(2) = find_closest_valid_dimension(combined_size(2));

% Only 6 FFTs are needed.
fixed_fft = fft2( ...
    fixed_image, ...
    optimal_size(1), ...
    optimal_size(2));

rotated_moving_fft = fft2( ...
    rotated_moving_image, ...
    optimal_size(1), ...
    optimal_size(2));

fixed_mask_fft = fft2( ...
    fixed_mask, ...
    optimal_size(1), ...
    optimal_size(2));

rotated_moving_mask_fft = fft2( ...
    rotated_moving_mask, ...
    optimal_size(1), ...
    optimal_size(2));

% Only 6 IFFTs are needed.
% Compute and save these results rather than computing them multiple times.
n_overlap_masked_pixel = real(ifft2( ... 
    rotated_moving_mask_fft .* fixed_mask_fft));

n_overlap_masked_pixel = round(n_overlap_masked_pixel);
n_overlap_masked_pixel = max(n_overlap_masked_pixel, eps);

mask_correlated_fixed_fft = real(ifft2( ...
    rotated_moving_mask_fft .* fixed_fft));

mask_correlated_rotated_moving_fft = real(ifft2( ...
    fixed_mask_fft .* rotated_moving_fft));

numerator = real(ifft2(rotated_moving_fft .* fixed_fft)) - ...
    mask_correlated_fixed_fft .* ...
    mask_correlated_rotated_moving_fft ./ ...
    n_overlap_masked_pixel;

fixed_squared_fft = fft2( ...
    fixed_image .* fixed_image, ...
    optimal_size(1), ...
    optimal_size(2));

fixed_denom = real(ifft2(rotated_moving_mask_fft.*fixed_squared_fft)) - ...
    mask_correlated_fixed_fft.^2 ./ ...
    n_overlap_masked_pixel;

fixed_denom = max(fixed_denom, 0);

rotated_moving_squared_fft = fft2( ...
    rotated_moving_image .* rotated_moving_image, ...
    optimal_size(1), ...
    optimal_size(2));

moving_denom = real(ifft2( ...
    fixed_mask_fft .* rotated_moving_squared_fft)) - ...
    mask_correlated_rotated_moving_fft.^2 ./ ...
    n_overlap_masked_pixel;
moving_denom = max(moving_denom,0);

denom = sqrt(fixed_denom .* moving_denom);

% denom is the sqrt of the product of positive numbers so it must be
% positive or zero.  Therefore, the only danger in dividing the numerator
% by the denominator is when dividing by zero.
% Since the correlation value must be between -1 and 1, we therefore
% saturate at these values.
correlation_matrix = zeros(size(numerator));
tol = 1000 * eps(max(abs(denom(:))));
i_nonzero = find(denom > tol);
correlation_matrix(i_nonzero) = numerator(i_nonzero) ./ denom(i_nonzero);
correlation_matrix(correlation_matrix < -1)  = -1;
correlation_matrix(correlation_matrix > 1)   = 1;

% Crop out the correct size.
correlation_matrix = correlation_matrix( ...
    1 : combined_size(1), 1 : combined_size(2));
n_overlap_masked_pixel = ...
    n_overlap_masked_pixel(1 : combined_size(1), 1 : combined_size(2));


correlation_matrix_Square = MakeSquare(correlation_matrix);


end

function matout = MakeSquare(matin)
% no longer pads, but instead cuts

a = diff(size(matin));

    if a<0 % if there are more rows than cols

        matout = matin(1+floor(-a/2):end-ceil(-a/2), :);
        
    elseif a>0 % if there are more cols than rows

        matout = matin(:, 1+floor(a/2):end-ceil(a/2));

    else
        matout = matin;
        
    end
    
end

function b = shift_data(x)
b = double(x);
is_unsigned = isa(x,'uint8') || isa(x,'uint16') || isa(x,'uint32');
if ~is_unsigned 
    min_b = min(b(:));
    if min_b < 0
        b = b - min_b;
    end
end
end


function [ new_n ] = find_closest_valid_dimension( n )
% Find the closest valid dimension above the desired dimension.  This
% will be a combination of 2s, 3s, and 5s.

% Incrementally add 1 to the size until
% we reach a size that can be properly factored.
new_n = n;
result = 0;
new_n = new_n - 1;
while(result ~= 1)
    new_n  = new_n + 1;
    result = factorize_number(new_n);
end
end


function [ n ] = factorize_number(n)

for ifac = [2 3 5]
    while(rem(n,ifac) == 0)
        n = n / ifac;
    end
end
end

