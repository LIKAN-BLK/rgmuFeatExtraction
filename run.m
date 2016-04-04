paths = {'../../data/01/','../../data/02/','../../data/06/','../../data/07/','../../data/08/','../../data/09/','../../data/10/','../../data/11/'};
expTitles={'e401','e402','e406','e407','e408','e409','e410','e411'};
eye_saveWaveAndAmple1(paths, expTitles, 'clfA')

% % -----------------DECIMATE-------------

for sc = 1:size(tmp,1)
    for t = 1:size(tmp,2)
        for ch = 1:size(tmp,3)
            for trial = 1:size(tmp,4)
                tmp = decimate(tmp(sc,:,ch,trial),10);
            end
        end
    end
end

% % -----------------MEAN-------------
tmp=load('wEEG_NTe406','-mat');
trial_mean = mean(tmp,4);
for i=1:23
    im = imagesc(trial_mean(:,1:10:750,i));
    saveas(im,['NONtarget' num2str(i,'%03u') '.png']);
end;

tmp = load('wEEG_Te406','-mat');
trial_mean = mean(tmp,4);
for i=1:23
    im = imagesc(trial_mean(:,1:10:750,i));
    saveas(im,['target' num2str(i,'%03u') '.png'])
end;

% % -----------------VARIANCE-------------
tmp=load('wEEG_NTe406','-mat');
% tmp = tmp.w;
for sc = 1:size(tmp,1)
    for ch = 1:size(tmp,3)
        for trial = 1:size(tmp,4)
            dtmp(sc,:,ch,trial) = decimate(tmp(sc,:,ch,trial),10);
        end
    end
end
tmp = dtmp;
trial_var = var(tmp,0,4);
for i=1:23
    im = imagesc(trial_var(:,:,i));
    saveas(im,['NONtarget' num2str(i,'%03u') '.png']);
end;

tmp = load('wEEG_Te406','-mat');
% tmp = tmp.w;
for sc = 1:size(tmp,1)
    for ch = 1:size(tmp,3)
        for trial = 1:size(tmp,4)
            dtmp(sc,:,ch,trial) = decimate(tmp(sc,:,ch,trial),10);
        end
    end    
end
tmp=dtmp;
trial_var = var(tmp,0,4);
for i=1:23
    im = imagesc(trial_var(:,:,i));
    colorbar;
    saveas(im,['target' num2str(i,'%03u') '.png'])
end;

%---------CALC R2 CUBE------------------
NT = load('wEEG_NTe406','-mat');
T = load('wEEG_Te406','-mat');
for ch = 1:size(NT,3)
    for s = 1:size(NT,1)
        for t = 1:size(NT,2)
            r2 = corrcoef(cat(4,NT(s,t,ch,:), T(s,t,ch,:)),[zeros(1,size(NT,4)) ones(1,size(T,4))]).^2;
            r2_mat(s,t,ch) = r2(1,2);
        end
    end
    im = imagesc(r2_mat(:,:,ch));
    colorbar;
    saveas(im,['401r2' num2str(ch,'%03u') '.png'])
end
save(['401R2'], 'r2_mat');
