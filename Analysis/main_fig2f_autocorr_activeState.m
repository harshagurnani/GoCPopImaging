%% Autocorrelation in modes from dff
% HG. Updated May 2020

% path to all folders
if ~exist( 'datapath', 'var'),   datapath = 'D:\Work\OneDrive - University College London\pubs and work\Golgi in vivo imaging\Paper\Datasets\'; end
% will save summary data into ..\FigureData if doSave==true
if ~exist( 'doSave', 'var' ), doSave = false; end
if ~exist( 'doPlot', 'var' ), doPlot = true; end

crus  = dir( [datapath,'Crus\*.mat'] );
lob45 = dir( [datapath,'Lob4_5\*.mat'] );
nCrus = length(crus);
nLob = length(lob45);

%% parameters and get autocorr
bin = 50; maxlag = 5e3;%ms
nPlot = maxlag/bin;
nTBins = 1+ 2*nPlot;
nPC = 5;
nShuffle = 1;

autocorr_ts = nan( nTBins, nPC, nCrus+nLob ); % lagTimepoints x PC# x Session
autocorr_ts_shuff = nan( nTBins, nPC, nCrus+nLob ); % for shuffled data

locoSD = nan( nCrus+nLob, 2); % duration variation, activity variation
locoMean = nan( nCrus+nLob, 2); % duration mean, activity mean

% all animals
for jj=1:nCrus+nLob
    if jj>nCrus, roi = jj-nCrus;  f = load( [datapath,'Lob4_5\',lob45(roi).name], 'allAnalysed', 'allData', 'allEvents');
    else, roi = jj; f = load( [datapath,'Crus\',crus(roi).name], 'allAnalysed', 'allData', 'allEvents'); end
    
    numModes = min( nPC, size(f.allAnalysed.PCA.dff.all.proj,1) ); %dff
    x = f.allAnalysed.PCA.dff.all.proj(1:numModes,:); % dynamics along PMs
    t = nanmean( f.allData.neurons.time,2);
 
    % crop during loco periods (else no activity dominates)
    loco = f.allEvents.locomotion;
    dur = loco(:,2)-loco(:,1);
    locoSD(jj,:) = [nanstd(dur), nanstd(loco(:,3))];
    locoMean(jj,:) = [nanmean(dur), nanmean(loco(:,3))];    
    if ~isempty(loco),   [x,t] = crop_data( x,t,loco, maxlag);    end
    tq = (1+floor(min(t)/bin))*bin:bin:max(t);
    for kk=1:numModes
        % resample only for active data
        vq = interp1( t, x(kk,:), tq );
        
        % lag correlations
        cc = xcorr( vq, maxlag/bin,'coeff' );
        autocorr_ts(:,kk,jj) = cc/cc(nPlot+1);
        
        shuffcc = [];
        for reps = 1:nShuffle
            id = randperm( length(vq) );
            vqs = vq(id);
            cc = xcorr( vqs-nanmean(vqs), maxlag/bin , 'coeff' );
            shuffcc = cat(1,shuffcc, cc/cc(nPlot+1));
        end
        % shuffled correlations
        autocorr_ts_shuff(:,kk,jj) = nanmean(shuffcc,1);
    end

end
if doSave
   fname = '..\FigureData\data_fig2f_autocorr_dff.mat';
   save( fname, 'autocorr_ts', 'autocorr_ts_shuff', '-v7.3' );
end


%% Plotting
if doPlot
   fig2e = plot_fig2e_autocorr_activeState( autocorr_ts, autocorr_ts_shuff, bin, maxlag, 'dff'); 
end


%% Stats
alpha = 0.05;
lags = 0:bin:maxlag;
TestLags = [500, maxlag]; nLags = length(TestLags);
for jj=1:nLags
    timeID = find( lags>=TestLags(jj),1);
    PM1 = squeeze( autocorr_ts( timeID+nPlot, 1, :) );
    nonPM1 = squeeze( nanmean(autocorr_ts( timeID+nPlot, 2:nPC, :),2) );
    [p,~] = signrank( PM1, nonPM1);
    fprintf(' ------------- Difference between autocorr of PM1 and other modes ---------- \n')
    if p<alpha
       fprintf( 'At lag of %d ms, significant difference \ndifference = %2.2f, p-value = %f \n', lags(timeID), nanmean(PM1-nonPM1), p);
    else
       fprintf( 'At lag of %d ms, no significant difference \ndifference =  %2.2f, p-value = %f \n', lags(timeID), nanmean(PM1-nonPM1), p);
    end
end




%% ---------------- helper
function [x,t] = crop_data( x,t,loco, extratime)

    loco = [loco(:,1)-.7*extratime, loco(:,2)+.7*extratime ];
    newx = [];
    newt = [];
    for jj=1:size(loco,1)
        t1 = find(t>loco(jj,1),1); if isempty(t1), t1 = 1; end
        t2 = find(t>loco(jj,2),1); if isempty(t2), t2 = length(t); end
        newx = cat(2, newx, x(:,t1:t2) ); newt = cat(1, newt, t(t1:t2));
    end
    bint = nanmean(diff(t));
    t = [1:size(newx,2)]*bint;
    x=newx;
end


