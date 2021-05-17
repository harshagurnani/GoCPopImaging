%% Cross Val SVD
% Cross-validated explained variance versus number of modes
% peak of curve as dimensionality estimate
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

%% cross-validated expvar and dim

nModes = 25;
%cvSVD curves
[cvExpVar.crus.dff, cvExpVar.crus.fr] = deal(nan(nModes,nCrus));
[cvExpVar.lob.dff, cvExpVar.lob.fr]   = deal(nan(nModes,nLob));
% dimensionality
[dim.peak.dff, dim.spectral.dff] = deal(nan( nCrus+nLob,1));
[dim.peak.fr, dim.spectral.fr]   = deal(nan( nCrus+nLob,1));

for roi=1:nCrus+nLob    
    if roi>nCrus
        jj=roi-nCrus; f = load( [datapath,'Lob4_5\',lob45(jj).name], 'allAnalysed', 'allData'); reg = 'lob';
    else
        jj=roi; f = load( [datapath,'Crus\',crus(jj).name], 'allAnalysed', 'allData');  reg = 'crus';
    end
    % for dff
    dtype = 'dff';
    x = nanmean( nanmean(f.allAnalysed.PCA.(dtype).crossval.res,3), 1);
    npc = min( nModes, length(x) );
    cvExpVar.(reg).(dtype)(1:npc,jj) = x(1:npc);
    dim.peak.(dtype)(roi) = f.allAnalysed.Dimensionality.(dtype).cv(1); 
    dim.spectral.(dtype)(roi) = f.allAnalysed.Dimensionality.(dtype).Spectral(1); 
    
    % for events
    dtype = 'fr';
    x = nanmean( nanmean(f.allAnalysed.PCA.(dtype).crossval.res,3), 1);
    npc = min( nModes, length(x) );
    cvExpVar.(reg).(dtype)(1:npc,jj) = x(1:npc);
    dim.peak.(dtype)(roi) = f.allAnalysed.Dimensionality.(dtype).cv(1); 
    dim.spectral.(dtype)(roi) = f.allAnalysed.Dimensionality.(dtype).Spectral(1); 
    
    nNeu_dim.dff(roi) = length(f.allData.ROI.keep_f);
    nNeu_dim.fr(roi) = length(f.allData.ROI.keep_spikes);
end

if doSave
   fname = '..\FigureData\data_fig6abc_dimensionality.mat';
   save( fname, 'cvExpVar', 'dim', 'nNeu_dim', '-v7.3' );
end

%% Plotting
if doPlot
    [fig4ab, lm.dff, lm.fr] = plot_fig6abc_Dimensionality( cvExpVar, dim, nNeu_dim );
end
