clear all;
clc;

load('pointTargetData.mat');

data = veraStrct.data;
fs = 20e6;
[rows_data,col_data, z_data] = size(data);
time = [0:1:rows_data-1]*(1/fs);
time_upsample = [0:1/4:rows_data-1]*(1/fs);

for kk = 1:128;
        data_upsample(:,:,kk) = interp1(time,data(:,:,kk),time_upsample,'linear');
end

fs_upsample = 4*fs;
speed = 1540; %m/s in body
pixel_size_through_depth = speed/fs_upsample; 
time_zero = 80;

for dd = 0.04:0.01:0.08
    
depth = dd; %m
for ii = 1:64
    xe(ii) = 0.1953e-3*ii;
    d(ii) = (xe(ii)^2+depth^2)^0.5;
    time_to_point(ii) = d(ii)/speed;
end

ttp_neg_xe = fliplr(time_to_point);
ttp_all_xe = [ttp_neg_xe time_to_point];
time_from_zero = ttp_all_xe(64);
time_from_zero_v = ones(1,length(ttp_all_xe))*time_from_zero;
time_delay = ttp_all_xe - time_from_zero_v;

time_array = [0:1/fs_upsample:max(size(data_upsample))];

for bb = 1:128
for aa = 1:128
    delay = time_delay(aa);
    time_array_delayed(:,aa,bb) = time_array(:,aa,bb)*delay;
    delayed_channel(:,aa,bb) = interp1(time_array_delayed(:,aa,bb),data_upsample(:,aa,bb),time_array,'linear');
end
end

% for zz = 1:128
%     delayed_channel(:,:,zz) = interp1(time_array_delayed,data_upsample(:,:,zz),time_array,'linear');
% end
    
figure;
min_data = min(min(min(delayed_channel)));
max_data = max(max(max(delayed_channel)));
imagesc(delayed_channel(:,:),[min_data max_data])
colormap('gray');
if dd == 0.04
    title('Channel data with delays (pointTargetData.mat), focus at 0.04 m');
elseif dd == 0.05
    title('Channel data with delays (pointTargetData.mat), focus at 0.05 m');
elseif dd == 0.06
    title('Channel data with delays (pointTargetData.mat), focus at 0.06 m');
elseif dd == 0.07
    title('Channel data with delays (pointTargetData.mat), focus at 0.07 m');
elseif dd == 0.08
    title('Channel data with delays (pointTargetData.mat), focus at 0.08 m');
end

summed_channels = sum(delayed_channel,2);
figure;
imagesc(20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
if dd == 0.04
    title('Compressed B-mode image (pointTargetData.mat), focus at 0.04 m');
elseif dd == 0.05
    title('Compressed B-mode image (pointTargetData.mat), focus at 0.05 m');
elseif dd == 0.06
    title('Compressed B-mode image (pointTargetData.mat), focus at 0.06 m');
elseif dd == 0.07
    title('Compressed B-mode image (pointTargetData.mat), focus at 0.07 m');
elseif dd == 0.08
    title('Compressed B-mode image (pointTargetData.mat), focus at 0.08 m');
end

end


load('anecoicCystData.mat');

data = veraStrct.data;
fs = 20e6;
[rows_data,col_data, z_data] = size(data);
time = [0:1:rows_data-1]*(1/fs);
time_upsample = [0:1/4:rows_data-1]*(1/fs);

for kk = 1:128;
        data_upsample(:,:,kk) = interp1(time,data(:,:,kk),time_upsample,'linear');
end

fs_upsample = 4*fs;
speed = 1540; %m/s in body
pixel_size_through_depth = speed/fs_upsample; 
time_zero = 80;

for dd = 0.04:0.01:0.08
    
depth = dd; %m
for ii = 1:64
    xe(ii) = 0.1953e-3*ii;
    d(ii) = (xe(ii)^2+depth^2)^0.5;
    time_to_point(ii) = d(ii)/speed;
end

ttp_neg_xe = fliplr(time_to_point);
ttp_all_xe = [ttp_neg_xe time_to_point];
time_from_zero = ttp_all_xe(64);
time_from_zero_v = ones(1,length(ttp_all_xe))*time_from_zero;
time_delay = ttp_all_xe - time_from_zero_v;

time_array = [0:1/fs_upsample:max(size(data_upsample))];

for bb = 1:128
for aa = 1:128
    delay = time_delay(aa);
    time_array_delayed(:,aa,bb) = time_array(:,aa,bb)*delay;
    delayed_channel(:,aa,bb) = interp1(time_array_delayed(:,aa,bb),data_upsample(:,aa,bb),time_array,'linear');
end
end

% for zz = 1:128
%     delayed_channel(:,:,zz) = interp1(time_array_delayed,data_upsample(:,:,zz),time_array,'linear');
% end
    
figure;
min_data = min(min(min(delayed_channel)));
max_data = max(max(max(delayed_channel)));
imagesc(delayed_channel(:,:),[min_data max_data])
colormap('gray');
if dd == 0.04
    title('Channel data with delays (anecoicCystData.mat), focus at 0.04 m');
elseif dd == 0.05
    title('Channel data with delays (anecoicCystData.mat), focus at 0.05 m');
elseif dd == 0.06
    title('Channel data with delays (anecoicCystData.mat), focus at 0.06 m');
elseif dd == 0.07
    title('Channel data with delays (anecoicCystData.mat), focus at 0.07 m');
elseif dd == 0.08
    title('Channel data with delays (anecoicCystData.mat), focus at 0.08 m');
end

summed_channels = sum(delayed_channel,2);
figure;
imagesc(20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
if dd == 0.04
    title('Compressed B-mode image (anecoicCystData.mat), focus at 0.04 m');
elseif dd == 0.05
    title('Compressed B-mode image (anecoicCystData.mat), focus at 0.05 m');
elseif dd == 0.06
    title('Compressed B-mode image (anecoicCystData.mat), focus at 0.06 m');
elseif dd == 0.07
    title('Compressed B-mode image (anecoicCystData.mat), focus at 0.07 m');
elseif dd == 0.08
    title('Compressed B-mode image (anecoicCystData.mat), focus at 0.08 m');
end

end





