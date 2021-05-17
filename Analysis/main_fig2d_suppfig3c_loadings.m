% Factor loadings 
% Uses scatter_kde.m from File Exchange by Nils Haentjens
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


%% load data

loadings.dff = [];
loadings.fr = [];
reg = [];

for jj=1:nCrus + nLob    
    if jj>nCrus
        f = load( [datapath,'Lob4_5\',lob45(jj-nCrus).name], 'allAnalysed'); rn = 2;
    else
        f = load( [datapath,'Crus\',crus(jj).name], 'allAnalysed'); rn=1;
    end
    tmp = f.allAnalysed.PCA.dff.all.eigvec;
    % fix for session 11 and 6 (pc1 has negative dynamics)
    if jj==11 || jj==6, tmp(:,1) = -tmp(:,1); end
    loadings.dff =  cat(1, loadings.dff, tmp(:,1:3));   % pcs for top 3 modes  -dff
    reg = cat(1,reg,rn*ones(size(tmp,1),1));
    tmp = f.allAnalysed.PCA.fr.all.eigvec;
    loadings.fr =  cat(1, loadings.fr, tmp(:,1:3));     % pcs for top 3 modes  -events
end
if doSave
   fname = '..\FigureData\data_fig2d_s3g_loadings.mat';
   save( fname, 'loadings', '-v7.3' );
end

%% Plotting
if doPlot
    plot_fig2d_suppfig3c_loadings( loadings )
end