function eye_saveWaveAndAmple1(paths, expTitles, clf_type)
%
%
%

latency = 500; % fixation duration for training
btn = []; % all button positions 'R' - only right, 'L' - only left

fs = 500;
beg_time = -0.5; % seconds
end_time = 1; % seconds

if strcmp(clf_type, 'clfA')
    eventsT = {'msgbuttonPressed'}; %clfA
elseif strcmp(clf_type, 'clfB')
    eventsT = {'msgbuttonPressed' 'msgballChosen' 'msgBallMoved'}; %clfB
end

chNames = {'FZ','F3','F4','Cz','C3','C4','PZ','P3','P4','P1','P2','PO7','PO8','PO3','PO4','Oz','O1','O2','POz','FP2'}; % all channel set
chNamesAll = [chNames, 'VEOG', 'HEOGL', 'HEOGR'];

chans = {'FZ','F3','F4','Cz','C3','C4','PZ','P3','P4','P1','P2','PO7','PO8','PO3','PO4','Oz','O1','O2','POz'}; %work channel set
f_channels = zeros(length(chans), 1);
for i = 1:length(chans)
        f_channels(i) =  find(strcmpi(chNames, chans{i}));
end

exp = [1:8];%1:8; % number of experiment

for n=exp
    path = paths{n};
    expTitle = expTitles{n};
    expTitle = [expTitle btn];
    % load trials
    [eegT, eegNT, labelsT] = eye_loaddata_r2e(path, end_time*1000, beg_time*1000, expTitle, latency, btn, eventsT);
    
    [X1, X0] = makeAmplFeatures(eegT, eegNT, f_channels, fs, beg_time);
    
    mat_dir = '../mat/wLets/';
    mkdir([mat_dir expTitle(2:end) '/'])
    save([mat_dir expTitle(2:end) '/' 'aEEG_NT' expTitle],'X0');
    save([mat_dir expTitle(2:end) '/' 'aEEG_T' expTitle],'X1');
    eegT = eegT(275:500,:,:); %50-500ms
    eegNT = eegNT(275:500,:,:);


    %-----------------------------------------------------------------------------------------------
    

    freqs_lo = 5; % Hz
    freqs_hi = 30; % Hz
    freqs = freqs_lo:0.5:freqs_hi; % Hz
    scales = fliplr(fs ./ freqs);
    freqs = scal2frq(scales, 'cmor2-1', 1/fs); % double check

    x = eegT(:, 1, 1);
    s = abs(cwt(x-mean(x), scales, 'cmor2-1'));
    s1 = abs(cwt(x+1000000, scales, 'cmor2-1'));
    mask = abs(s-s1) < 10000;

    N1 = size(eegT, 3);
    N0 = size(eegNT, 3);

    w= zeros(size(scales,2),size(eegT, 1),size(f_channels,1),size(eegT, 3));
    for i = 1:N1
        disp(i);
        for ch = 1:size(f_channels, 1)
            x(:,ch) = eegT(:, f_channels(ch), i);
            s = abs(cwt(x(:,ch)-mean(x(:,ch)), scales, 'cmor2-1'));        
            w(:,:,ch,i) = log(s);        
        end
    %     save([mat_dir 'wEEG_' expTitle '#' num2str(i,'%03u')], 'w', 'x');
    end
    
    save([mat_dir expTitle(2:end) '/' 'wEEG_T' expTitle], 'w');



    % -------------------NONTarget---------------------------------
    nw = zeros(size(scales,2),size(eegNT, 1),size(f_channels,1),size(eegNT, 3));
    for i = 1:N0
        disp(i);
        for ch = 1:size(f_channels,1)
            x(:,ch) = eegNT(:, f_channels(ch), i);
            s = abs(cwt(x(:,ch)-mean(x(:,ch)), scales, 'cmor2-1'));        
            nw(:,:,ch,i) = log(s);        
        end
    %     save([mat_dir 'wEEG_' expTitle '#' num2str(N1+i,'%03u')], 'w', 'x');
    end
    save([mat_dir expTitle(2:end) '/' 'wEEG_NT' expTitle], 'nw');

    
%   mkdir([mat_dir expTitle(2:end) '/' '/mean' expTitle])
%   mkdir([mat_dir expTitle(2:end) '/' '/var' expTitle])
%   mkdir([mat_dir expTitle(2:end) '/' '/R2' expTitle])


%     % % -----------------MEAN-------------
% 
%     trial_mean = mean(w,4);
%     for i=1:23
%         im = imagesc(trial_mean(:,:,i));
%         saveas(im,[mat_dir expTitle(2:end) '/mean' expTitle '/' 'target' num2str(i,'%03u') '.png']);
%     end;
% 
%     %----------------------NONTARGET MEAN---------------------------
%     trial_mean = mean(nw,4);
%     for i=1:23
%         im = imagesc(trial_mean(:,:,i));
%         saveas(im,[mat_dir expTitle(2:end) '/mean' expTitle '/' 'NONtarget' num2str(i,'%03u') '.png']);
%     end;
%     %  ---------------------R2-----------------------
%     for ch = 1:size(nw,3)
%         for s = 1:size(nw,1)
%             for t = 1:size(nw,2)
%                 r2 = corrcoef(cat(4,nw(s,t,ch,:), w(s,t,ch,:)),[zeros(1,size(nw,4)) ones(1,size(w,4))]).^2;
%                 r2_mat(s,t,ch) = r2(1,2);
%             end
%         end
%         im = imagesc(r2_mat(:,:,ch));
%         colorbar;
%         saveas(im,[mat_dir expTitle(2:end) '/R2' expTitle '/' 'ch' num2str(ch,'%03u') '.png'])
%     end
%     save([mat_dir expTitle(2:end) '/' 'R2' expTitle], 'r2_mat');


    % % -----------------VAR-------------

    % for sc = 1:size(w,1)
    %     for ch = 1:size(w,3)
    %         for trial = 1:size(w,4)
    %             dw(sc,:,ch,trial) = decimate(w(sc,:,ch,trial),10);
    %         end
    %     end
    % end
    % w = dw;
    % trial_var = var(w,0,4);
    % for i=1:23
    %     im = imagesc(trial_var(:,:,i));
    %     saveas(im,[mat_dir expTitle(2:end) '/var' expTitle '/' 'target' num2str(i,'%03u') '.png']);
    % end;
    % %----------------------NONTARGET VAR---------------------------
    % 
    % for sc = 1:size(nw,1)
    %     for ch = 1:size(nw,3)
    %         for trial = 1:size(nw,4)
    %             dnw(sc,:,ch,trial) = decimate(nw(sc,:,ch,trial),10);
    %         end
    %     end
    % end
    % nw = dnw;
    % trial_var = var(nw,0,4);
    % for i=1:23
    %     im = imagesc(trial_var(:,:,i));
    %     saveas(im,[mat_dir expTitle(2:end) '/var' expTitle '/' 'NONtarget' num2str(i,'%03u') '.png']);
    % end;



    times = beg_time:1/fs:end_time;
    times = times(2:end);
    labels = [labelsT; zeros(N0, 1)];
    save([mat_dir expTitle(2:end) '/wEEG_info'], 'mask', 'times', 'freqs', 'labels', 'chNamesAll');
end
