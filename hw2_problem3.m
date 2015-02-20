clear all;
clc;

load('pointTargetData.mat');

data = veraStrct.data;
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = speed/fs; 

for ii = 1:max(size(data))
    time_array_all(ii) = ii/fs;
end

channel = [[-63.5:1:63.5]];

zone1 = data([1:485],:,:);
zone2 = data([1:485*2],:,:);
zone3 = data([1:485*3],:,:);
zone4 = data([1:485*4],:,:);
zone5 = data(:,:,:);
zone = struct('zone_1',zone1, 'zone_2',zone2 ,'zone_3',zone3, 'zone_4',zone4, 'zone_5',zone5);

for z = 1:5

    if z == 1
        data_matrix = zone.zone_1;
    elseif z == 2
        data_matrix = zone.zone_2;
    elseif z == 3
        data_matrix = zone.zone_3;
    elseif z == 4
        data_matrix = zone.zone_4;
    elseif z == 5
        data_matrix = zone.zone_5;
    end
last_row = max(size(data_matrix));
depth_row = last_row - 150;
depth = depth_row*pixel_size_through_depth; %m

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    d(ii) = (xe(ii)^2+depth^2)^0.5;
    time_to_point(ii) = d(ii)/speed;
end
time_from_zero = time_to_point(64);
time_from_zero_v = ones(1,length(time_to_point))*time_from_zero;
time_delay = time_to_point - time_from_zero_v;

time_array = time_array_all(1:max(size(data_matrix)));


for aa = 1:128
    delay = time_delay(aa);
    time_array_delayed = time_array+delay;
    delayed_channel([1:max(size(data_matrix))],aa,:) = interp1(time_array,data_matrix([1:max(size(data_matrix))],aa,:),time_array_delayed,'linear');
end
tic
for cc = 1:last_row
    if z ==1
        delayed_data([1:485],:,:) = delayed_channel([1:485],:,:);
    elseif z==2
        delayed_data([486:485*2],:,:) = delayed_channel([486:485*2],:,:);
    elseif z==3
        delayed_data([485*2+1:485*3],:,:) = delayed_channel([485*2+1:485*3],:,:);
    elseif z==4
        delayed_data([485*3+1:485*4],:,:) = delayed_channel([485*3+1:485*4],:,:);
    elseif z==5
        delayed_data([485*4+1:max(size(delayed_channel))],:,:) = delayed_channel([485*4+1:max(size(delayed_channel))],:,:);
    end
end
toc
end

for ll = 1:numel(delayed_data)
    if isnan(delayed_data(ll))==1
        delayed_data(ll) = 0;
    end
end

figure;
min_data = min(min(min(delayed_data)));
max_data = max(max(max(delayed_data)));
imagesc(delayed_data(:,:),[min_data, max_data])
colormap('gray');
title('Channel data with delays (pointTargetData.mat)');

summed_channels = sum(delayed_data,2);
figure;
imagesc(20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
title('Compressed B-mode image (pointTargetData.mat)');


load('anecoicCystData.mat');

data = veraStrct.data;
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = speed/fs; 

for ii = 1:max(size(data))
    time_array_all(ii) = ii/fs;
end

channel = [[-64:0] [0:64]];

zone1 = data([1:485],:,:);
zone2 = data([1:485*2],:,:);
zone3 = data([1:485*3],:,:);
zone4 = data([1:485*4],:,:);
zone5 = data(:,:,:);
zone = struct('zone_1',zone1, 'zone_2',zone2 ,'zone_3',zone3, 'zone_4',zone4, 'zone_5',zone5);

for z = 1:5

    if z == 1
        data_matrix = zone.zone_1;
    elseif z == 2
        data_matrix = zone.zone_2;
    elseif z == 3
        data_matrix = zone.zone_3;
    elseif z == 4
        data_matrix = zone.zone_4;
    elseif z == 5
        data_matrix = zone.zone_5;
    end
last_row = max(size(data_matrix));
depth_row = last_row - 150;
depth = depth_row*pixel_size_through_depth; %m
    
for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    d(ii) = (xe(ii)^2+depth^2)^0.5;
    time_to_point(ii) = d(ii)/speed;
end
time_from_zero = time_to_point(64);
time_from_zero_v = ones(1,length(time_to_point))*time_from_zero;
time_delay = time_to_point - time_from_zero_v;

time_array = time_array_all(1:max(size(data_matrix)));

for aa = 1:128
    delay = time_delay(aa);
    time_array_delayed = time_array+delay;
    delayed_channel([1:max(size(data_matrix))],aa,:) = interp1(time_array_delayed,data_matrix([1:max(size(data_matrix))],aa,:),time_array,'linear');
end

for cc = 1:last_row
    if z ==1
        delayed_data([1:485],:,:) = delayed_channel([1:485],:,:);
    elseif z==2
        delayed_data([486:485*2],:,:) = delayed_channel([486:485*2],:,:);
    elseif z==3
        delayed_data([485*2+1:485*3],:,:) = delayed_channel([485*2+1:485*3],:,:);
    elseif z==4
        delayed_data([485*3+1:485*4],:,:) = delayed_channel([485*3+1:485*4],:,:);
    elseif z==5
        delayed_data([485*4+1:max(size(delayed_channel))],:,:) = delayed_channel([485*4+1:max(size(delayed_channel))],:,:);
    end
end
end

for ll = 1:numel(delayed_data)
    if isnan(delayed_data(ll))==1
        delayed_data(ll) = 0;
    end
end

figure;
min_data = min(min(min(delayed_data)));
max_data = max(max(max(delayed_data)));
imagesc(delayed_data(:,:),[min_data, max_data])
colormap('gray');
title('Channel data with delays (anecoicCystData.mat)');

summed_channels = sum(delayed_data,2);
figure;
imagesc(20*log10(abs(hilbert(summed_channels(:,:)))));
colormap('gray');
title('Compressed B-mode image (anecoicCystData.mat)');






