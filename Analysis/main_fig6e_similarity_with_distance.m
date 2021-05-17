%% Distribution of modes
% dot product of loadings versus distance
%
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


%% load data and compute dot product between modes
SignalSim = []; pw_dist = [];

for roi=1:nCrus+nLob
    if roi>nCrus
        jj=roi-nCrus; f = load( [datapath,'Lob4_5\',lob45(jj).name], 'allAnalysed');
    else
        jj=roi; f = load( [datapath,'Crus\',crus(jj).name], 'allAnalysed'); 
    end
    use_pc = f.allAnalysed.Dimensionality.dff.cv(1); use_pc = 1:use_pc; % only use signal modes
    dist = f.allAnalysed.pw_distance;
    x = f.allAnalysed.PCA.dff.all.eigvec(:,use_pc);
    nNeu = size(dist,1);
    ids = find(triu(ones(nNeu),1)==1);    % upper tri block
    cc = x*x';
    % all pairs
    SignalSim = cat(1, SignalSim, cc(ids(:)) );
    pw_dist = cat(1, pw_dist, dist(ids(:)) );
end

if doSave
   fname = '..\FigureData\data_fig6e_similarity.mat';
   save( fname, 'SignalSim', 'pw_dist', '-v7.3' );
end

%% Plotting
if doPlot
   [fig4d, linMdl ] = plot_fig6e_similarity_with_distance ( SignalSim, pw_dist);
end

% Linear fir
fprintf(' ------ Linear fit to similarity versus distance --------- \n');
fprintf('Adj. Rsquared = %2.2f,  Slope = %2.3f/300um \n', linMdl.Rsquared.Adjusted, linMdl.Coefficients.Estimate(2)*300);