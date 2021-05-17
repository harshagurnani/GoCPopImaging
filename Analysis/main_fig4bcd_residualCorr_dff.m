%% Residual correlations 
% - after projecting out PM1
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

%% load data for both lobules

pw_dist.all = []; pw_dist.pos = []; pw_dist.neg = [];
TotalCorr.all = [];
[ResCorr.all, ResCorr.pos, ResCorr.neg]  = deal([]); % corr, region id

TotalCorr.shuffCorr = []; TotalCorr.shuffCorr = [];
ResCorr.shuffCorr = []; ResCorr.shuffCorr = [];
blockLen = 0.5;

sessMeanCorr = [];

for roi=1:nCrus+nLob
    
    if roi>nCrus
        jj=roi-nCrus; f = load( [datapath,'Lob4_5\',lob45(jj).name], 'allAnalysed', 'allData'); reg=2;
    else
        jj=roi; f = load( [datapath,'Crus\',crus(jj).name], 'allAnalysed', 'allData'); reg=1;
    end
    dist = f.allAnalysed.pw_distance;
    x = f.allAnalysed.correlations.noPC1.pw.corr;
    xTotal = f.allAnalysed.correlations.all.pw.corr;
    sig = f.allAnalysed.correlations.noPC1.pw.sig;
    
    nNeu = size(dist,1);
    ids = find(triu(ones(nNeu),1)==1); % all pairs - counted once
    ResCorr.all = cat(1, ResCorr.all, [x(ids(:)), reg*ones(length(ids),1)] ); % all pairs - rescorr
    TotalCorr.all = cat(1, TotalCorr.all, [xTotal(ids(:)), reg*ones(length(ids),1)] ); % all pairs - totalcorr
    pw_dist.all = cat(1, pw_dist.all, dist(ids(:)) );
    sessMeanCorr = cat( 1, sessMeanCorr, [nanmean(x(ids(:))), reg] ); % session mean
    
    ids = find(sig>0 & triu(ones(nNeu),1)==1);  % significantly positive
    ResCorr.pos = cat(1, ResCorr.pos, [x(ids(:)), reg*ones(length(ids),1)] );
    pw_dist.pos = cat(1, pw_dist.pos, dist(ids(:)) );
    
    ids = find(sig<0 & triu(ones(nNeu),1)==1);  % significantly negative
    ResCorr.neg = cat(1, ResCorr.neg, [x(ids(:)), reg*ones(length(ids),1)] );
    pw_dist.neg = cat(1, pw_dist.neg, dist(ids(:)) );
    
    
    % shuffle distribution
    x = f.allData.neurons.f;    shuffX = x;
    y = subspace_svd(x, -1);     shuffY = y;
    tm = f.allData.neurons.time;
    [T, nNeu] = size(x);
    maxT = max(tm(:));
    try
        rate = f.allData.params.acquisition_rate;
    catch
        dt = nanmean(diff(tm(:,1)));
        rate = 1000/dt;
    end
    
    for neu=1:nNeu
        id = block_shuffle_time(T, rate, blockLen);
        shuffX(:, neu) = x(id, neu);
        shuffY(:, neu) = y(id, neu);        
    end
    
    [~, corr, ~] = all_period_lagcorr_pairwise_dff( [0, maxT], shuffX, tm, 20, 0 );
    corr = squeeze(corr);
    
    [~, corr2, ~] = all_period_lagcorr_pairwise_dff( [0, maxT], shuffY, tm, 20, 0 );
    corr2 = squeeze(corr2);
    
    ids = find(triu(ones(nNeu),1)==1);  % upper block     
    
    TotalCorr.shuffCorr = cat(1, TotalCorr.shuffCorr, corr(ids(:))); 
    ResCorr.shuffCorr = cat(1, ResCorr.shuffCorr, corr2(ids(:))); 


    
end

if doSave
   fname = '..\FigureData\data_fig4bcd_residualcorr_dff.mat';
   save( fname, 'ResCorr', 'sessMeanCorr', 'pw_dist', 'TotalCorr', '-v7.3' );
end


%% Plotting
if doPlot
    fig4bcd = plot_fig4bcd_residualCorr( ResCorr, TotalCorr, pw_dist, sessMeanCorr, 'dff' );
end