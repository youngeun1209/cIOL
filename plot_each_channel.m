function plot_each_channel(data, fs, segtime, varargin)
% plot each channel in segment time
% example: 
%           plot_each_channel(EPOX, 100, [0 10], 'chanName', clab, ...
%               'nTrial',1, 'channels', {'Cz','Pz','Oz'}, 'scale', 200, ...
%               'str_title', 'epoched X')
%
% input:    EPOX :      continuous or epoched data [time X channel] or [time X channel X trial]
%           fs   :      sampling frequency of the data [int]
%           segtime :   selected time segment in second [time_start, time_end]
%           varargin:   
%               chanName -  the total Name of channel [cell]
%               nTrial   -  selected trial index [int]
%               channels -  selected channel name [cell]
%               scale    -  scale of the plotting for each channel [double]
%               str_title - the title name of plot [str]

if isempty(data)
    error('check your SMT data')
end

if isempty(fs)
    fs = 500;
end

if isempty(segtime)
    t_end = size(data,1)/fs;
    segtime=[0 t_end];
end

varargin = reshape(varargin,2,[])';
opt = cell2struct(varargin(:,2), varargin(:,1));

if isfield(opt,'chanName')
    clab = opt.chanName;
else
    clab = 1:size(data,2);
end

if isfield(opt,'nTrial')
    nTrial = opt.nTrial;
else
    nTrial =1;
end

if isfield(opt, 'channels')
    selected_ch = opt.channels;
    
    if ~prod(ismember(selected_ch, clab))
        error('check selected channels')
    end
    
    ch_idx = ismember(clab,selected_ch);
    data = data(:,ch_idx,:);
    clab = clab(ch_idx);
end

str_title = 'time domain plot';
if isfield(opt,'title')
    str_title = opt.title;
end

if isfield(opt,'en_text')
    en_text = opt.en_text;
else
    en_text = true;
end


% time segment
t=segtime(1)+1/fs:1/fs:segtime(2);
size_data = size(data);
if length(size_data) == 3
    if length(t) > size_data(1)
        error('Selected segment should be shorter than input data segment')
    end
    chVec = squeeze(data(int64(t*fs),:,nTrial));
elseif length(size_data) == 2
    if length(clab) == size_data(2)
        chVec = squeeze(data(int64(t*fs),:));
    else
        chVec = squeeze(data(int64(t*fs),nTrial));
    end
end
    
% scale for plotting
windowSize = fs*0.01; b = (1/windowSize)*ones(1,windowSize); a = 1; 
t_dat_filt = filter(b,a,chVec);
plot_max = max(reshape(abs(t_dat_filt),1,[]));
if plot_max == 0
    plot_scale = 1;
    warning('plot scale is 0')
elseif plot_max < 10 && plot_max >= 4
    plot_scale = 10;
elseif  plot_max < 5 && plot_max >= 1
    plot_scale = 5;
elseif  plot_max < 2 && plot_max >=0.5
    plot_scale = 1;
elseif plot_max < 0.5
    plot_scale = plot_max*1.5;
else
    plot_scale = round((plot_max*1.5)/10)*10;
end

if isfield(opt,'scale')
    plot_scale = opt.scale;
end

% bias
bias = plot_scale * 2;
bias = cumsum(repmat(bias, 1, size(chVec,2)));

chVec = chVec + flip(bias);

% figure;
plot(t, chVec,'k')
% ylim([0-bias(1) bias(end)+bias(1)*2])
ylim([0 bias(end)+bias(1)])

xlim(segtime)
yticks(bias)
yticklabels(flip(clab))
title(str_title);
ylabel('channels')
xlabel('time [s]')
str = sprintf('scale: %d', plot_scale);  % 
if en_text == true
annotation('textbox',[.01 .68 .3 .3],'String',str,'LineStyle','none');
end
grid on
grid minor
