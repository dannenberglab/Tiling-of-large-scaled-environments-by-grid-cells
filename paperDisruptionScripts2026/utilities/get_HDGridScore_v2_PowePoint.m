function [HDgridScore]=get_HDGridScore_v2_PowePoint(props);
%this function computes an alternative grid score based on the
%powerspectrum of the ellipticity-corrected ('grid3',1)
%ring-autocorrelogram  computed by the CMBHOME gridness function.

% % % ratemap_binsize=1;
% % % 
% % % %set epoch
% % % root.epoch=epoch;
% % % [gridness3Score, props] = root.Gridness(cel, 'binside', ratemap_binsize, 'continuize_epochs', 1, 'grid3', 1);
map=props.auto_corr;
map(~props.auto_corr_mask)=nan; 

if size(map,1)~=size(map,2)
    disp('Autocorrelogram is not quadratic.')
    return
end

CenterX=ceil(length(map)/2);
CenterY=ceil(length(map)/2);

angle=nan(1,size(map,1)*2+size(map,2)*2-3); %initialize
profileLines=cell(size(angle)); %initialize
profileLines_woNan = cell(size(angle));
profile=nan(size(angle)); %initialize
%first quadrant
c=0;
% first half of first quadrant
for i=CenterY:size(map,2) %
    c=c+1; %counter
    x=[CenterX size(map,1)]; y=[CenterY i];
    angle(c)=atan((y(1,2)-y(1,1))/(x(1,2)-x(1,1)));
    profileLines{c}=improfile(map,x,y);
    %%%
    profileLines_woNan{c} = profileLines{c}(~isnan(profileLines{c}));
    % get midpoint
    midpoint = length(profileLines_woNan{c})/2;
    if isinteger(midpoint)
        profile(c) = mean(profileLines_woNan{c}(midpoint:midpoint + 1));
    else
%         cellfun(@(profileLines_woNan) nansum(x,ndims(x)+1), raw_data, 'UniformOutput', false)
%         profileLines_woNan(isnan(profileLines_woNan))=0;
        profile(c) = profileLines_woNan{c}(round(midpoint));
 
    end
    %%%
end
%second half of first quadrant
for i=1:CenterX-1 %
    c=c+1; %counter
    x=[CenterX size(map,1)-i]; y=[CenterY size(map,2)];
    angle(c)=atan((y(1,2)-y(1,1))/(x(1,2)-x(1,1)));
    profileLines{c}=improfile(map,x,y);
    %%%
    profileLines_woNan{c} = profileLines{c}(~isnan(profileLines{c}));
    % get midpoint
    midpoint = length(profileLines_woNan{c})/2;
    if isinteger(midpoint)
        profile(c) = mean(profileLines_woNan{c}(midpoint:midpoint + 1));
    else
        profile(c) = profileLines_woNan{c}(round(midpoint));
    end
    %%%
end
% first half of second quarter
for i=CenterX-1:size(map,1)-1 %
    if i==CenterX-1 % for correct computation of angle
        c=c+1; %counter
        x=[CenterX size(map,1)-i]; y=[CenterY size(map,2)];
        angle(c)=atan((y(1,2)-y(1,1))/(x(1,2)-x(1,1)));
        profileLines{c}=improfile(map,x,y);
    %%%
    profileLines_woNan{c} = profileLines{c}(~isnan(profileLines{c}));
    % get midpoint
    midpoint = length(profileLines_woNan{c})/2;
    if isinteger(midpoint)
        profile(c) = mean(profileLines_woNan{c}(midpoint:midpoint + 1));
    else
        profile(c) = profileLines_woNan{c}(round(midpoint));
    end
    %%%
    else
        c=c+1; %counter
        x=[CenterX size(map,1)-i]; y=[CenterY size(map,2)];
        angle(c)=atan((y(1,2)-y(1,1))/(x(1,2)-x(1,1)))+pi;
        profileLines{c}=improfile(map,x,y);
    %%%
    profileLines_woNan{c} = profileLines{c}(~isnan(profileLines{c}));
    % get midpoint
    midpoint = length(profileLines_woNan{c})/2;
    if isinteger(midpoint)
        profile(c) = mean(profileLines_woNan{c}(midpoint:midpoint + 1));
    else
        profile(c) = profileLines_woNan{c}(round(midpoint));
    end
    %%%
    end
end
% second half of second quarter
for i=1:floor(size(map,2)/2) %
    c=c+1; %counter
    x=[CenterX 1]; y=[CenterY size(map,2)-i];
    angle(c)=atan((y(1,2)-y(1,1))/(x(1,2)-x(1,1)))+pi;
    profileLines{c}=improfile(map,x,y);
    %%%
    profileLines_woNan{c} = profileLines{c}(~isnan(profileLines{c}));
    % get midpoint
    midpoint = length(profileLines_woNan{c})/2;
    if isinteger(midpoint)
        profile(c) = mean(profileLines_woNan{c}(midpoint:midpoint + 1));
    else
        profile(c) = profileLines_woNan{c}(round(midpoint));
    end
    %%%
end
%first half of third quarter
for i=CenterY:size(map,2)-1 %
    c=c+1; %counter
    x=[CenterX 1]; y=[CenterY size(map,2)-i];
    angle(c)=atan((y(1,2)-y(1,1))/(x(1,2)-x(1,1)))+pi;
    profileLines{c}=improfile(map,x,y);
    %%%
    profileLines_woNan{c} = profileLines{c}(~isnan(profileLines{c}));
    % get midpoint
    midpoint = length(profileLines_woNan{c})/2;
    if isinteger(midpoint)
        profile(c) = mean(profileLines_woNan{c}(midpoint:midpoint + 1));
    else
        profile(c) = profileLines_woNan{c}(round(midpoint));
    end
    %%%
end
%second half of third quarter
for i=2:ceil(size(map,2)/2)-1 %
    c=c+1; %counter
    x=[CenterX i]; y=[CenterY 1];
    angle(c)=atan((y(1,2)-y(1,1))/(x(1,2)-x(1,1)))+pi;
    profileLines{c}=improfile(map,x,y);
    %%%
    profileLines_woNan{c} = profileLines{c}(~isnan(profileLines{c}));
    % get midpoint
    midpoint = length(profileLines_woNan{c})/2;
    if isinteger(midpoint)
        profile(c) = mean(profileLines_woNan{c}(midpoint:midpoint + 1));
    else
        profile(c) = profileLines_woNan{c}(round(midpoint));
    end
    %%%
end
%first half of fourth quarter
for i=CenterX:size(map,2) %
    c=c+1; %counter
    x=[CenterX i]; y=[CenterY 1];
    angle(c)=atan((y(1,2)-y(1,1))/(x(1,2)-x(1,1)))+2*pi;
    profileLines{c}=improfile(map,x,y);
    %%%
    profileLines_woNan{c} = profileLines{c}(~isnan(profileLines{c}));
    % get midpoint
    midpoint = length(profileLines_woNan{c})/2;
    if isinteger(midpoint)
        profile(c) = mean(profileLines_woNan{c}(midpoint:midpoint + 1));
    else
        profile(c) = profileLines_woNan{c}(round(midpoint));
    end
    %%%
end
%second half of fourth quarter
for i=2:floor(size(map,2)/2) %
    c=c+1; %counter
    x=[CenterX size(map,2)]; y=[CenterY i];
    angle(c)=atan((y(1,2)-y(1,1))/(x(1,2)-x(1,1)))+2*pi;
        profileLines{c}=improfile(map,x,y);
    %%%
    profileLines_woNan{c} = profileLines{c}(~isnan(profileLines{c}));
    % get midpoint
    midpoint = length(profileLines_woNan{c})/2;
    if isinteger(midpoint)
        profile(c) = mean(profileLines_woNan{c}(midpoint:midpoint + 1));
    else
        profile(c) = profileLines_woNan{c}(round(midpoint));
    end
    %%%
end

%% interpolate profile
% angleDeg=linspace(0,360,length(profile)+1);
x=[rad2deg(angle) 360];
[~, indexXUnique]=unique(x);
v=[profile profile(1)];
temp=interp1(x(indexXUnique),v(indexXUnique),0:360);
profileDeg=temp(1:end-1); %delete last value because it is a duplicate (0 and 360 degree is the same).
% %plot figure
% figure
% plot(angleDeg(mm,:),profileDeg(mm,:))

%% compute FFT of profile
% Fs=360; %Sampling frequency
% T=1/Fs; %Sampling period
L=360; %Length of signal
% t=(0:L-1)*T; %Time vector
Y=fft(profileDeg-mean(profileDeg)); % perform FFT
P2=abs(Y/L);
P1=P2(1:L/2+1);
P1(2:end-1)=2*P1(2:end-1);
% f=Fs*(0:(L/2))/L;
% figure
% plot(f,P1)

% pull out result
HDgridScore=P1(7);

