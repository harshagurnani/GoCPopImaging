%% Pairwise correlations  - dF/F
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


%% load/save summary data

% all pairs                  %pw distance & sessID       % significance             % mean session corr 
totalCorr.crus.pwCorr = [];    totalCorr.crus.pwDist = [];   totalCorr.crus.corrSig = [];   totalCorr.crus.sessCorr = [];     
totalCorr.lob45.pwCorr = [];   totalCorr.lob45.pwDist = [];  totalCorr.lob45.corrSig = [];  totalCorr.lob45.sessCorr = [];    

for exp=1:nCrus + nLob
    if exp> nCrus
        jj=exp-nCrus; reg = 'lob45'; 
        f = load( [datapath,'Lob4_5\',lob45(jj).name], 'allAnalysed'); 
    else
        jj=exp; reg = 'crus'; 
        f = load( [datapath,'Crus\',crus(jj).name], 'allAnalysed'); 
    end
    x = f.allAnalysed.correlations.all.pw.corr;     % all corr
    sig = f.allAnalysed.correlations.all.pw.sig;    % significance (outside 95% CI of shuffle control)
    d = f.allAnalysed.pw_distance;                  % distance between GoCs
    
    nNeu = size(x,1);
    ids = find(triu(ones(nNeu),1)==1);  % upper block
    totalCorr.(reg).pwCorr  = cat(1, totalCorr.(reg).pwCorr, x(ids(:)));
    totalCorr.(reg).pwDist  = cat(1, totalCorr.(reg).pwDist, [d(ids(:)), jj*ones(length(ids),1)]);
    totalCorr.(reg).corrSig   = cat(1, totalCorr.(reg).corrSig, sig(ids(:)));
    totalCorr.(reg).sessCorr = cat(1, totalCorr.(reg).sessCorr, nanmean(x(ids(:))) );

end

if doSave
    fname = '..\FigureData\data_fig1fg_totalCorr_dff.mat';
    save( fname, 'totalCorr', '-v7.3' );
end

%%      Plotting
%------------------------------------------------------------------------%
if doPlot
    fig1fg = plot_fig1fg_s2cd_pairwise_corr( totalCorr, 'dff');
end