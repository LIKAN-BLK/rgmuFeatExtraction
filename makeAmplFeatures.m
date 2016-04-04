function [X1, X0] = makeAmplFeatures(eegT, eegNT, f_channels, fs, beg_time)
%
[eegT, eegNT] = preprocess(eegT, eegNT, fs, beg_time);
N0 = size(eegNT, 3);
N1 = size(eegT, 3);
nChannels = size(eegT, 2);
T = size(eegT, 1);

% without filtering
eegTfilt = eegT;
eegNTfilt = eegNT;

% features
times_beg = [0.2:0.02:0.45];
times_end = times_beg + 0.05;
ts_beg = round((times_beg - beg_time) .* fs);
ts_end = round((times_end - beg_time) .* fs);

X0 = [];
X1 = [];    
for i = 1:N1
    x = [];
    for t = 1:length(ts_beg)    
        x = [x; mean(eegTfilt(ts_beg(t):ts_end(t), f_channels, i), 1)];
    end
    X1(i, :) = x(:);
end

for i = 1:N0
    x = [];
    for t = 1:length(ts_beg)    
        x = [x; mean(eegNTfilt(ts_beg(t):ts_end(t), f_channels, i), 1)];
    end
    X0(i, :) = x(:);
end 
end


function [eegT, eegNT] = preprocess(eegT, eegNT, fs, beg_time)
%


% baseline correction
base_beg_time = 0.2;
base_end_time = 0.3;
base_beg = (base_beg_time - beg_time) * fs;
base_end = (base_end_time - beg_time) * fs;
eegT = eegT - repmat(mean(eegT(base_beg:base_end,:,:),1), [size(eegT,1), 1, 1]);
eegNT = eegNT - repmat(mean(eegNT(base_beg:base_end,:,:),1), [size(eegNT,1), 1, 1]);


end

