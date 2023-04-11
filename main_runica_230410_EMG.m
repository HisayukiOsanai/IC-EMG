% This is an example code that shows how to extract EMG signals from
% skull-screw-referenced LFP data using independent component analysis (ICA).
%
% Input:
%   ICEMG_threshold: threshold of weight distribution uniformity for
%                    detecting IC-EMG.
%                    default: ICEMG_threshold = 0.1;
%   filename: path to LFP file (causion: ICA computaion may take long time for heavy input files)
%   chn:    number of LFP channels
%   Fs:     sampling rate [Hz]
%   badch:  abnormal ch# you want to avoid in computing ICA. Including
%           abnormal channels may cause unwanted results.
%           example:    badch = []; if all channels are healthy
%                       badch = [15, 27]; if you want to avoid ch#15 and #27
%
%   (optional) EMGfilename: path to EMG file
%
%
% Output:
%   EEG.icaact = ICs obtained by ICA from LFP
%   EEG.icaweights, EEG.icasphere: both are relating to IC weight distribution.
%       LFP and IC have relationship that  LFP = mixing*IC,
%         where  mixing = pinv(EEG.icaweights * EEG.icasphere).
%       Weight distribution of each IC is shown in each column of "mixing" matrix.
%
%   ICEMG:  IC-EMG, identified with uniformity of IC weight distribution.
%   ICEMGi: Index of IC-EMG in ICs.
% 
% Hisayuki Osanai, 04/10/2023, Hisayuki.Osanai@utsouthwestern.edu


close all
clear


%% settings

% threshold of weight distribution uniformity for detecting IC-EMG
% default: ICEMG_threshold = 0.1;
% increase (e.g. ICEMG_threshold = 0.2) if LFP data length is short
ICEMG_threshold = 0.1;

chn = 32; % number of LFP channels

Fs=2000; % sampling rate [Hz]

svpath = strcat(pwd, '\results\');
filepath = strcat(pwd, '\example\');

%% impport LFP files

fprintf('Reading files...\n');

filename{1} = strcat(filepath, 'LFP.mat');

badch{1} = [];
% badch{1} = [15 27];

%% ICA

for k = 1:length(filename)

    % load LFP data
    [path, file, ext] = fileparts(filename{k});
    file = strcat(file , ext);
    path = strcat(path , '\');

    S = load(filename{k});
    LFP = S.EEG.data; % example data contains 32-ch LFP data in "EEG.data" variable


    % delete bad channels
    badchannel = badch{k};
    if min(badchannel) > 0
        LFP(badchannel,:) = [];
    end


    %% run ICA

    Samples = LFP;
    LFP_beforecrop = Samples;

    % Use 4-min data for ICA if total LFP recording is longer than 5min to
    % reduce computational cost.
    % You can reduce the length if you think computation time is too long.
    % We recommend to use data of at least 10-sec data (~20000-sampling points).
    if length(LFP)/Fs/60 > 5
        startI = 1 + Fs*0*60;
        stopI = startI + Fs*4*60; % 4-min
        Samples = Samples(:, startI:stopI);
    end


    %EEG.data = LFP; % LFP after deleting abnormal channels

    tmpdata = Samples;
    tmpdata = double(tmpdata);
    tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); % zero mean

    % run ICA
    tic
    [weights,sphere,meanvar,bias,signs,lrates,data,y] = runica(tmpdata, 'lrate', 0.001, 'extended', 1);
    % in case if you may want to use PCA before ICA
    %[weights,sphere,meanvar,bias,signs,lrates,data,y] = runica(tmpdata, 'lrate', 0.001, 'extended', 1, 'pca', 20);
    toc

    tic
    unmixing_matrix = weights * sphere;
    IC_normalized = unmixing_matrix * LFP_beforecrop;

    fprintf('calculated ICA, file #%d \n', k);

    %% plot row LFP
    fprintf('plot LPFs\n');

    Wave = LFP;

    % Replace Nan to bad channel data
    a = Wave;
    b = NaN(1,length(a));
    c = a;
    for i = 1:1:length(badchannel)
        c = insertRows(c, b, badchannel(i));
    end
    Wave = c;

    fig_LFP = func_plotWaves(Wave, Fs);

    nCh = size(Wave, 1);
    chind = cell(1, nCh);
    for i = 1:1:nCh
        chind{i} = strcat('ch', num2str(i));
    end
    GCA = gca;
    GCA.YTickLabel =  fliplr( chind);
    title('LFP')


    %% plot ICs and their spatial weight distributions

    fprintf('plot ICs\n');
    fig_IC = func_plotWaves(IC_normalized, Fs);
    title('Normalized ICs')


    %before scaling. ICs values are origninally normalized.
    %EEG.etc.icaweights_beforerms = weights;
    %EEG.etc.icasphere_beforerms = sphere;

    %% scale ICA weights
    fprintf('Scaling ICs\n');

    IC_normalized = unmixing_matrix * LFP;
    EEG.data = LFP;

    [icaweights_new, icawinv, IC, fig_ICscaled] = func_plotScaledIC(EEG, Fs, weights, sphere);
    set(gcf,'Renderer', 'Painters');

    EEG.icaweights = icaweights_new;
    EEG.icasphere = sphere;
    EEG.icawinv = icawinv;
    EEG.icaact = IC;
    %EEG.icaact_normalized = IC_normalized; %normalized IC
    EEG.srate = Fs;


    spweight = pinv(icaweights_new*sphere)'; % spatial weight
    [fig_Weights, nweights] = func_plotWeight(chn, spweight, badchannel);

    %set(gcf,'Renderer', 'Painters');

    %% Identify IC-EMG (IC with most uniform spatial weight distribution)

    W = nweights;

    MW = mean(W, 2);
    SW = std(W, [], 2);

    fig_scatter = figure;
    hold on
    plot(MW, SW, '.', 'MarkerSize', 5 )

    xlim([0 1])
    ylim([0 1])
    xlabel('mean')
    ylabel('standard deviation')
    title('Scatter plot of mean weight and std')

    % get(gcf,'Position')
    set(gcf, 'Position', [1005         736         279         228]);
    ax = gca;
    ax.FontSize = 10;


    % identify IC-EMG as most-uniform IC with lower std than threshold
    [Min_SW, ICEMGi] = min(SW);
    % change ylabel color if IC-EMG was detected.
    if Min_SW < ICEMG_threshold
        set(0, 'currentfigure', fig_IC);
        GCA = gca;
        GCA.YTickLabel{size(GCA.YTickLabel,1)+1 - ICEMGi} = ['\color{red}' GCA.YTickLabel{size(GCA.YTickLabel,1)+1 - ICEMGi}];

        set(0, 'currentfigure', fig_ICscaled);
        GCA = gca;
        GCA.YTickLabel{size(GCA.YTickLabel,1)+1 - ICEMGi} = ['\color{red}' GCA.YTickLabel{size(GCA.YTickLabel,1)+1 - ICEMGi}];

        set(0, 'currentfigure', fig_Weights);
        AX = fig_Weights.Children;
        AX_ICEMG = AX(size(AX,1)+1 - ICEMGi);
        AX_ICEMG.YLabel.Color = 'r';
    end

    ICEMG = IC(ICEMGi,:);

    %% IC back-projections to LFP channels


    mixing = spweight';
    recLFP = mixing(:,ICEMGi)*IC(ICEMGi,:);
    LFP_sub = LFP - recLFP;

    Wave = recLFP;
    a = Wave; b = NaN(1,length(a)); c = a;
    for i = 1:1:length(badchannel)
        c = insertRows(c, b, badchannel(i));
    end
    recLFP = c;

    Wave = LFP_sub;
    a = Wave; b = NaN(1,length(a)); c = a;
    for i = 1:1:length(badchannel)
        c = insertRows(c, b, badchannel(i));
    end
    LFP_sub = c;

    fig_recLFP = func_plotWaves(recLFP, Fs);
    title('Projection of uniform-IC (IC-EMG)')
    %set(gcf,'Renderer', 'Painters');

    fig_subLFP = func_plotWaves(LFP_sub, Fs);
    title('uniform-IC-deleted LFP')
    %set(gcf,'Renderer', 'Painters');


    %% plot EMGs (if available)

    %         EMGfilename = strcat(filepath, 'EMG.mat');
    %         E = load(EMGfilename);
    %
    %
    %         EMG = E.EEG.data;
    %
    %         Wave = EMG;
    %         Wave_EMG = [Wave; EMG(1,:) - EMG(2,:)];
    %         Wave_ICEMG = [ICEMG;  Wave_EMG];
    %
    %         %z-scoring
    %         Wave_z = Wave_ICEMG;
    %         for i = 1:4
    %             Wave_z(i,:) = (Wave_ICEMG(i,:) - mean(Wave_ICEMG(i,:)) ) ./ std( Wave_ICEMG(i,:) );
    %         end
    %         fig_EMG = func_plotWaves(Wave_z, Fs);
    %
    %         pos = [ 1004          362         484         276];
    %         set(gcf,'Position',pos);
    %
    %         YLabels = {'IC-EMG', 'R-EMG', 'L-EMG', '(R-L)-EMG'};
    %         GCA = gca;
    %         GCA.YTickLabel =  fliplr( YLabels);
    %         title('EMG')
    %
    %         %     set(gcf,'Renderer', 'Painters');
    %         s = strcat(svpath, 'ICEMG.png');
    %         saveas(fig_EMG, s);


    %% save figures and .mat file

    s = strcat(svpath, 'RawLFP.png');
    saveas(fig_LFP, s);

    s = strcat(svpath, 'ICs_normalized.png');
    saveas(fig_IC, s);

    s = strcat(svpath, 'ICs.png');
    saveas(fig_ICscaled, s);

    s = strcat(svpath, 'ICweights.png');
    saveas(fig_Weights, s);

    s = strcat(svpath, 'uniformIC projection.png');
    saveas(fig_recLFP, s);

    s = strcat(svpath, 'uniformIC-deleted LFP.png');
    saveas(fig_subLFP, s);

    s = strcat(svpath, 'scatter of weight distribution.png');
    saveas(fig_scatter, s);

    s = strcat(svpath, 'IC-EMG results.mat');
    save(s ,'EEG', 'ICEMG', 'ICEMGi', '-v7.3'); %https://jp.mathworks.com/matlabcentral/answers/220925-save

    %close all

    toc

end



%% functions

%% plot LFP/ICA waveforms
function fig = func_plotWaves(Wave, Fs)

% Wave = LFP;
nCh = size(Wave, 1);
fig = figure;
%get(gcf,'Position')
pos = [ 1104          62         484         876];
set(gcf,'Position',pos);
hold on

% display 0-5 sec waveforms
ts = 0; %[s]
te = ts + 5;

ts_i = ts*Fs + 1; %time[s] -> index
te_i = te*Fs + 1;
Wave = Wave(:, ts_i:te_i); % crop LFP into display epoch time

Time = 0:1/Fs: (length(Wave) - 1)/Fs;
Time = Time(:, ts_i:te_i);

ymax = max(max(abs( Wave))) * 1.2; % setting y axis scale
for k = 1:1:nCh %length(ch_i)
    %     plot(t, Wave(k,:) - (k-1)*ymax ,'k', 'LineWidth', 0.5)
    plot(Time, Wave(k,:) - (k-1)*ymax ,'k', 'LineWidth', 0.2)
end


ylim([-ymax*nCh,  ymax*1])
% yLimits = get(gca,'YLim');

GCA = gca;
GCA.YTick = -ymax*(nCh-1) : ymax : 0;

chind = cell(1,nCh);
for k = 1:1:nCh
    chind{k} = strcat('c', num2str(k));
end

GCA.YTickLabel =  fliplr( chind);
GCA.XLabel.String = 'time [s]';


end


%%  plot ICA spatial weights
function [fig, nweights] = func_plotWeight(chn, spweight, badchannel)

chi = 1:1:chn;

fig = figure;
pos = [940    75   488   866];
set(gcf,'Position',pos);


size(spweight,1);
Nh = ceil(size(spweight,1)/2); Nw = 2; gap = [.009 .08]; marg_h =[.03 .03]; marg_w = [.1 .1]; %marg_h =[.03 .01]; marg_w = [.02 .02];
[ha, ~] = tight_subplot(Nh, Nw, gap, marg_h, marg_w); %https://www.mathworks.com/matlabcentral/fileexchange/27991-tight-subplot-nh--nw--gap--marg-h--marg-w-
X = [1, chn];

% normalize weight distribution
nweights = spweight ./max( abs(spweight  ), [], 2);

%inserse if weight distribution tends to be negative
N = size(nweights, 1);
for i = 1:1:N
    if mean(nweights) <0
        nweights(i,:) = - nweights(i,:);
    end
end

% plot
if min(badchannel) > 0
    chi(badchannel) = [];
end
i = 0;
for k = 1:1:N
    if k <= length(chi)
        i = i+1;
        plot(ha(k), chi, nweights(i,:) ,  'k-', 'Linewidth', 1);
    end

    ylim(ha(k), [-1, 1])
    xlim(ha(k), X)
    ha(k).XTick=[]; ha(k).XTickLabel=[];
    ha(k).YTick=[-1 0 1]; %ha(k).YTickLabel=[strcat('c', num2str(ici(k)))];
    hold(ha(k),'on')
    plot(ha(k), [0 chn], [0 0], '--', 'Color',[0 0 0]+0.5);
    ha(k).Box = 'off';
    ha(k).YLabel.String = strcat('c', num2str(k));
end

Nh = chn - length(badchannel);
ha(Nh).XAxis.Color = 'k';
ha(Nh).XLim = [1 chn]; ha(Nh).YLim = [-1 1];
ha(Nh).XTick= [1, round(chn/2), chn]; ha(Nh).XTickLabel= [1, round(chn/2), chn];
ha(Nh).YTick=[-1, 0, 1]; ha(Nh).YTickLabel=[-1, 0, 1];

ha(Nh-1).XAxis.Color = 'k';
ha(Nh-1).XLim = [1 chn]; ha(Nh-1).YLim = [-1 1];
ha(Nh-1).XTick= [1, round(chn/2), chn]; ha(Nh-1).XTickLabel= [1, round(chn/2), chn];
ha(Nh-1).YTick=[-1, 0, 1]; ha(Nh-1).YTickLabel=[-1, 0, 1];

ha(1).Title.String = 'ICA spatial weights';

ha(Nh-1).XLabel.String = 'channel';
ha(Nh).XLabel.String = 'channel';

end


%% calculate scaled ICs
function [icaweights_new, icawinv, ICAcomponents_scaled, fig] = func_plotScaledIC(EEG, Fs, weights, sphere)


icawinv_temp   = pinv(weights * sphere); % a priori same result as inv
scaling = repmat(sqrt(mean(icawinv_temp.^2))', [1 size(weights,2)]);

icaweights_new = weights .* scaling;
unmixing_matrix = icaweights_new * sphere;
ICAcomponents = unmixing_matrix * EEG.data;

icawinv = icawinv_temp;

fig = func_plotWaves(ICAcomponents, Fs);
ICAcomponents_scaled = ICAcomponents;
title('Scaled ICs')

end


%% Insert NaN to bad channels

function [ result ] = insertRows( inMat , doMat , numRow )
[ rowsDo , colsDo ] = size( doMat );
[ rowsIn , colsIn ] = size( inMat );

if colsDo > colsIn;    error('too many rows'); end

insMat = zeros( rowsDo , colsIn );
insMat( 1:rowsDo , 1:colsDo ) = doMat;
result = cat( 1 , inMat( 1:numRow -1 , : ) , insMat , inMat( numRow : rowsIn , : ) );

% ref: https://yurufuwa-engineer.blogspot.com/2013/12/matlab.html
end