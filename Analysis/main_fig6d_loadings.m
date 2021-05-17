% Factor loadings 
% Uses scatter_kde.m from File Exchange by Nils Haentjens
% HG. Updated May 2020

% path to all folders
if ~exist( 'datapath', 'var'),   datapath = 'D:\Work\OneDrive - University College London\pubs and work\Golgi in vivo imaging\Paper\Datasets\'; end
% will save summary data into ..\FigureData if doSave==true

crus  = dir( [datapath,'Crus\*.mat'] );
lob45 = dir( [datapath,'Lob4_5\*.mat'] );
nCrus = length(crus);
nLob = length(lob45);


%% load data

loadings.dff = [];
loadings.fr = [];
loadings.var = [];
loadings.mean = [];
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
    loadings.dff =  cat(1, loadings.dff, tmp(:,1:5));   % pcs for top 3 modes  -dff
    loadings.var =  cat(1, loadings.var, nanstd(tmp(:,1:5),[],1));   % pcs for top 3 modes  -dff
    loadings.mean =  cat(1, loadings.mean, nanmean(tmp(:,1:5),1));   % pcs for top 3 modes  -dff
    
    tmp = f.allAnalysed.PCA.fr.all.eigvec;
    loadings.fr =  cat(1, loadings.fr, tmp(:,1:5));     % pcs for top 3 modes  -events

    reg = cat(1,reg,jj*ones(size(tmp,1),1));
end

if doSave
   fname = '..\FigureData\data_fig6d_loadings.mat';
   save( fname, 'loadings', '-v7.3' );
end

%% Plotting
if doPlot
    


fig6d = figure; hold on;
xt = [0,2:5];
clr =  [ 12, 185 192; 204, 43, 83]/255;

bb=bar( xt, nanmean(loadings.var,1), 'k');
bb.FaceAlpha=0;
for jj=1:nCrus+nLob
    scatter( xt, loadings.var(jj,:), 20, 'markerfacecolor', clr(reg(jj),:), 'markeredgecolor', clr(reg(jj),:)); 
end

errorbar( xt, nanmean(loadings.var,1), nanstd(loadings.var,[], 1), 'ko')

xticks(xt);
xticklabels({'PM1','PM2', 'PM3',  'PM4', 'PM5'});
ylabel('Variability of loadings');

title('Fig 6D - Variability of loadings (dFF)')

end