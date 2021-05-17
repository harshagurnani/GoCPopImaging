%% Eigenvalue Spectrum, and cross-validated Explained variance by PM1
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

nPC = 40;

% explained variance by PM1
[ExpVar_PM1.crus.dff, ExpVar_PM1.crus.fr] = deal(nan(nCrus,1)); %dff//fr
[ExpVar_PM1.lob.dff, ExpVar_PM1.lob.fr]   = deal(nan(nLob,1));
% normalised eigenspectrum
[Spectrum.dff, Spectrum.fr] = deal(nan(nCrus+nLob,nPC)); %dff//fr [Sess x nPC ]

for exp=1:nCrus+nLob
    if exp>nCrus
        jj = exp-nCrus; reg = 'lob'; 
        f = load( [datapath,'Lob4_5\',lob45(jj).name], 'allAnalysed');
    else
        jj = exp; reg = 'crus';
        f = load( [datapath,'Crus\',crus(jj).name], 'allAnalysed');
    end
    
    % Explained variance by PM1 (cross-validated)
    tmp = f.allAnalysed.PCA.dff.crossval.res; %dff
    if size(tmp,1)>size(tmp,2), tmp=permute(tmp,[2,1,3]); end
    allcv = nanmean(nanmean(tmp,3),2);
    ExpVar_PM1.(reg).dff(jj) = allcv(1);
    
    tmp=f.allAnalysed.PCA.fr.crossval.res; %events
    if size(tmp,1)>size(tmp,2), tmp=permute(tmp,[2,1,3]); end
    allcv = nanmean(nanmean(tmp,3),2);
    ExpVar_PM1.(reg).fr(jj) = allcv(1);
    
    % Normalised eigenspectrum
    allpc = f.allAnalysed.PCA.dff.all.eig_val; allpc = allpc/allpc(1);  %dff
    nc = min(nPC, length(allpc));
    Spectrum.dff(exp,1:nc) = allpc(1:nc);
    
    allpc = f.allAnalysed.PCA.fr.all.eig_val; allpc = allpc/allpc(1); %events
    nc = min(nPC, length(allpc));
    Spectrum.fr(exp,1:nc) = allpc(1:nc);
end
if doSave
    fname = '..\FigureData\data_fig2b_d_ModesSpectrum.mat';
    save( fname, 'ExpVar_PM1', 'Spectrum', '-v7.3' ); 
end

%%      Plotting - Run from here if you have already loaded data
%------------------------------------------------------------------------%
if doPlot
    fig2bc = plot_fig2bc_suppfig3ab_eigenSpectrum( ExpVar_PM1, Spectrum);
end