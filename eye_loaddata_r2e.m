function [eegT, eegNT, varargout] = eye_loaddata_r2e(path, epoch_size, left_border, expTitle, latency, btn, varargin)
%
% make target and non-target EEG trials from data in files
%

if nargin >= 7
    eventsT = varargin{1};
else
    eventsT = {'msgbuttonPressed'};
end

% ---------- load data -------------------
json = fileread([path '/meta.json']);
[str ~] = parse_json(json);

sRate = 500;
epoch_size = epoch_size - left_border;
eegT = [];
eegNT = [];
labelsT = [];
for i = 1:length(str{1}.valid_files)
    edf_fname = str{1}.valid_files{i}.name_edf;
    ascii_fname = ['eventsLatencies' edf_fname '.ascii'];
    eeg_fname = str{1}.valid_files{i}.name_eeg;
    
    if ~isempty(latency)
        if (str2num(str{1}.valid_files{i}.latency) ~= latency)
            continue;
        end
    end
    
    if ~isempty(btn)
        if ~strcmp(str{1}.valid_files{i}.btn, btn)
            continue;
        end
    end

    if strcmp(expTitle, '06--01') || strcmp(expTitle, '07--01')
        load([path '/' eeg_fname '.raw.mat']);
        signal = raw;
        tmp = find(diff(raw(:,end) == 65282) == 1);    
        signal = signal(tmp(2)+1:end, :); % second label sync
    elseif strcmp(expTitle, 'test10')        
        load([path '/' eeg_fname '.mat']);
        signal = EEG;
    elseif strcmp(expTitle, 'e401') || strcmp(expTitle, 'e402') || strcmp(expTitle, 'e403') || strcmp(expTitle, 'e404') || strcmp(expTitle, 'e405') || strcmp(expTitle, 'e406') || ...
           strcmp(expTitle, 'e407') || strcmp(expTitle, 'e408') || strcmp(expTitle, 'e409') || strcmp(expTitle, 'e410') || strcmp(expTitle, 'e411')
        load([path '/' eeg_fname '.mat']);
        signal = EEG;
    end
    
    fid = fopen([path '/' ascii_fname]);
    tmp = fgetl(fid); % read header
    B = textscan(fid, '%f %s'); % B{1} - times, B{2} - labels
    fclose(fid);

    [eegTi, eegNTi, labelsTi] = loaddata1(signal, B, epoch_size, sRate, left_border, eventsT);    
    eegT = cat(3, eegT, eegTi);
    eegNT = cat(3, eegNT, eegNTi);
    labelsT = [labelsT; labelsTi];
end

varargout{1} = labelsT;
end

function [eegT, eegNT, labelsT] = loaddata1(signal, B, epoch_size, sRate, left_border, eventsT)
% ----------------------------------------

msgbuttonPressed_idx = find(strcmp(B{2}, 'msgbuttonPressed'));
msgbuttonPressed_times = B{1}(msgbuttonPressed_idx);
%msgbuttonPressed_times = B{1}(msgbuttonPressed_idx) / 2; % /2 - only for test10!
msgbuttonPressed_t = round(msgbuttonPressed_times * sRate);

msgballChosen_idx = find(strcmp(B{2}, 'msgballChosen'));
msgballChosen_times = B{1}(msgballChosen_idx);
msgballChosen_t = round(msgballChosen_times * sRate);

msgBallMoved_idx = find(strcmp(B{2}, 'msgBallMoved'));
msgBallMoved_times = B{1}(msgBallMoved_idx);
msgBallMoved_t = round(msgBallMoved_times * sRate);

msgClickedInBlockMode_idx = find(strcmp(B{2}, 'msgClickedInBlockMode'));
msgClickedInBlockMode_times = B{1}(msgClickedInBlockMode_idx);
msgClickedInBlockMode_t = round(msgClickedInBlockMode_times * sRate);

msgBallClickedInBlockedMode_idx = find(strcmp(B{2}, 'msgBallClickedInBlockedMode'));
msgBallClickedInBlockedMode_times = B{1}(msgBallClickedInBlockedMode_idx);
%msgBallClickedInBlockedMode_times = B{1}(msgBallClickedInBlockedMode_idx) / 2; % /2 - only for test10!
msgBallClickedInBlockedMode_t = round(msgBallClickedInBlockedMode_times * sRate);

msgBoardClickedInBlockedMode_idx = find(strcmp(B{2}, 'msgBoardClickedInBlockedMode'));
msgBoardClickedInBlockedMode_times = B{1}(msgBoardClickedInBlockedMode_idx);
msgBoardClickedInBlockedMode_t = round(msgBoardClickedInBlockedMode_times * sRate);

fixationDuration_t = round((epoch_size/1000) * sRate);

%eventsT_t = [msgbuttonPressed_t; msgballChosen_t; msgBallMoved_t];
%eventsNT_t = [msgClickedInBlockMode_t; msgBallClickedInBlockedMode_t; msgBoardClickedInBlockedMode_t];
%eventsT_t = [msgbuttonPressed_t];

eventsT_t = [];
labelsT = [];
for i = 1:length(eventsT)
    if strcmpi(eventsT{i}, 'msgbuttonPressed')
        eventsT_t = [eventsT_t; msgbuttonPressed_t];
        labelsT = [labelsT; 1*ones(length(msgbuttonPressed_t), 1)]; %'msgbuttonPressed'
    elseif strcmpi(eventsT{i}, 'msgballChosen')
        eventsT_t = [eventsT_t; msgballChosen_t];
        labelsT = [labelsT; 2*ones(length(msgballChosen_t), 1)]; %'msgballChosen'
    elseif strcmpi(eventsT{i}, 'msgBallMoved')
        eventsT_t = [eventsT_t; msgBallMoved_t];
        labelsT = [labelsT; 3*ones(length(msgBallMoved_t), 1)]; %'msgBallMoved'
    end
end

eventsNT_t = [msgBallClickedInBlockedMode_t];

% exclude trials with insufficient pre-onset
eventsT_t = eventsT_t(eventsT_t+left_border > 0);
eventsNT_t = eventsNT_t(eventsNT_t+left_border > 0);
labelsT = labelsT(eventsT_t+left_border > 0);

% ---------- make EEG trials ------------------------------
nChannels = size(signal, 2);
%signal = signal - repmat(mean(signal, 2), 1, nChannels);
left_border = left_border/1000*sRate;
eegT = zeros(fixationDuration_t, nChannels, length(eventsT_t));
for i = 1:length(eventsT_t)
    eegT(:, :, i) = signal(eventsT_t(i)+left_border+1:eventsT_t(i)+fixationDuration_t+left_border, :);    
end

eegNT = zeros(fixationDuration_t, nChannels, length(eventsNT_t));
for i = 1:length(eventsNT_t)
    eegNT(:, :, i) = signal(eventsNT_t(i)+left_border+1:eventsNT_t(i)+fixationDuration_t+left_border, :);    
end
i;
end