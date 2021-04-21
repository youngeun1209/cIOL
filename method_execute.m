%% Offline version for continuous data example
% cap_CNT = cap_cnt{9,4};
% IMU_CNT = IMU_cnt{9,4};
load('example_data')
X = cap_CNT.x;
Y = IMU_CNT.x;


% execute
[filt_cntX,ref_cICA] = cIOL(X, Y, 100, 'idxRef', 3, 'learningRate',0.001, ...
    'window_time', 2000,'moving_time', 500, 'flag_PCA', true);

% plot the results
clab = cap_CNT.clab;
figure(1); plot_each_channel(X, 100, [100 110], 'nTrial',200, 'fs', 100, 'chanName', clab, 'channels',{'Cz','Pz','POz','Oz'}, 'scale',100)
figure(2); plot_each_channel(filt_cntX, 100, [100 110], 'nTrial',200, 'fs', 100, 'chanName', clab, 'channels',{'Cz','Pz','POz','Oz'}, 'scale',100)

%% Pseudo online version for continuous data example
X = cap_CNT.x;
Y = IMU_CNT.x;
fs = cap_CNT.fs;

filt_cntX = [];
online_time = 2000; % 2 seconds
online_size = online_time/1000*fs;
% tic
for i = 1:size(X,1)/online_size
%     toc
%     tic
    seg = (1:online_size) + online_size*(i-1);
    
    % execute
    [filt_cntX_temp,ref_cICA_temp] = cIOL( X(seg,:), Y(seg,:), fs, 'idxRef', 3, 'learningRate',0.001, ...
        'window_time', 1000,'moving_time', 200, 'flag_PCA', true);
    filt_cntX = [filt_cntX; filt_cntX_temp];
    
end
filt_cntX = [filt_cntX; X(length(filt_cntX)+1:end,:)];

% plot the results
clab = cap_CNT.clab;
figure(1); plot_each_channel(X, 100, [100 110], 'nTrial',200, 'fs', 100, 'chanName', clab, 'channels',{'Cz','Pz','POz','Oz'}, 'scale',100)
figure(2); plot_each_channel(filt_cntX, 100, [100 110], 'nTrial',200, 'fs', 100, 'chanName', clab, 'channels',{'Cz','Pz','POz','Oz'}, 'scale',100)

%% Offline version for epoched data example
cap_EPO = cap_epo{9,4};
IMU_EPO = IMU_epo{9,4};
X = cap_EPO.x;
Y = IMU_EPO.x;

% execute
[filt_cntX,ref_cICA] = cIOL(X, Y, 100, 'idxRef', 1:3, 'learningRate',0.001, ...
    'window_time', 200,'moving_time', 100);

% plot the results
clab = cap_EPO.clab;
figure(1); plot_each_channel(X, 100, [0 1], 'nTrial',200, 'fs', 100, 'chanName', clab, 'channels',{'Cz','Pz','POz','Oz'}, 'scale',100)
figure(2); plot_each_channel(filt_cntX, 100, [0 1], 'nTrial',200, 'fs', 100, 'chanName', clab, 'channels',{'Cz','Pz','POz','Oz'}, 'scale',100)

