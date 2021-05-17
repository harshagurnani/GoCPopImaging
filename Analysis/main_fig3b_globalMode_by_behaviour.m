%% pc1 - explained variance by behaviour during active state 
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
flds = { 'lincomb_state', 'state', 'wmi', 'loco', 'wamp','wangle','wsp'};
labels = {'LinComb', 'State', 'WMI', 'Loco', 'WAmp', 'Angle', 'WSP'};
nBeh = length(flds);

[PC1EVBeh.crus.data, PC1EVBeh.crus.shuffle] = deal(nan(nCrus,nBeh)); % all sessions and behaviours, for real data and shuffled control
[PC1EVBeh.lob.data, PC1EVBeh.lob.shuffle] = deal(nan(nLob,nBeh));

for exp=1:nCrus+nLob
    if exp> nCrus
        jj=exp-nCrus; reg = 'lob'; 
        f = load( [datapath,'Lob4_5\',lob45(jj).name], 'allAnalysed', 'allData', 'allEvents'); 
    else
        jj=exp; reg = 'crus'; 
        f = load( [datapath,'Crus\',crus(jj).name], 'allAnalysed', 'allData', 'allEvents' ); 
    end
    mdl = get_lin_model_active( f.allData, f.allAnalysed, f.allEvents);
    for ff=1:length(flds)
        beh = flds{ff};
        if isfield( mdl.pcs, beh )   % maximal cv exp var for PM1
            PC1EVBeh.(reg).data(jj,ff) = max(mdl.pcs.(beh).test(1,:)); %cross-validated
            PC1EVBeh.(reg).shuffle(jj,ff) = max(mdl.pcs.(beh).shuff(1,:)); 
        end
    end
end
    
if doSave
   fname = '..\FigureData\data_fig3b_globalMode_dff_withbehaviour.mat';
   save( fname, 'PC1EVBeh', 'flds', '-v7.3' );
end


%% Plotting
allExp.data = cat(1, PC1EVBeh.crus.data, PC1EVBeh.lob.data); %combined
allExp.shuffle = cat(1, PC1EVBeh.crus.shuffle, PC1EVBeh.lob.shuffle);
mc = abs(allExp.data);
for sessn=1:size(allExp.data,1)
   for beh = 1:size(allExp.data,2)
      if  allExp.data(sessn,beh)< abs(allExp.shuffle(sessn,beh)), mc(sessn,beh)=NaN; end % worse than shuffle -> Nan (terrible fit)
   end
end


flds = { 'lincomb_state', 'state', 'wmi', 'loco', 'wamp','wangle','wsp'};
labels = {'LinComb', 'State', 'WMI', 'Loco', 'WAmp', 'Angle', 'WSP'};
nBeh = length(flds);

if doPlot
    fig2g = plot_fig2g_globalMode_by_behaviour( mc, nCrus, labels, 'dff' );
end

%% Stats
alphaF = 0.05;
fprintf( '----- Improvement to Linear Model --------- \n' )
alpha = alphaF/(size(allExp.data,2)-1);  %Bonferroni correction by #beh
for beh = 2:size(allExp.data,2)
    id = find( ~isnan( squeeze(allExp.data(:,1)) ) & ~isnan( squeeze(allExp.data(:,beh)) ) );
    [p,~] = signrank( squeeze(allExp.data(id,1)), squeeze(allExp.data(id,beh)), 'tail','right' );
    
    if p<alpha
       fprintf( 'Significant increase from %s to Linear model, \n p-value = %f (One-sided Signed Rank test) \n', labels{beh}, p );
    else
       fprintf( 'No significant increase from %s to Linear model, \n p-value = %f (One-sided Signed Rank test) \n', labels{beh}, p );
    end
    
end

alpha = alphaF/2;  %Bonferroni correction
fprintf( '\n----- Comparision across Lobule --------- \n' )
for beh = [3,4]
    id1 = find( ~isnan(squeeze( PC1EVBeh.crus.data(:,beh))) ) ;
    id2 = find( ~isnan(squeeze(PC1EVBeh.lob.data(:,beh))) );
    [p,~] = ranksum( squeeze( PC1EVBeh.crus.data(id1,beh)), squeeze(PC1EVBeh.lob.data(id2,beh)) );
    if p<alpha
       fprintf( 'Significant difference for %s, \n p-value = %f (Signed Rank test) \n', labels{beh}, p );
    else
       fprintf( 'No significant difference for %s, \n p-value = %f (Signed Rank test) \n', labels{beh}, p );
    end
end

alpha = alphaF/2;  %Bonferroni correction
fprintf( '\n----- Comparision with state --------- \n' )
for beh = [3,5]
    id = find( ~isnan( squeeze(allExp.data(:,2)) ) & ~isnan( squeeze(allExp.data(:,beh)) ) );
    [p,~] = signrank( squeeze(allExp.data(id,2)), squeeze(allExp.data(id,beh)), 'tail','left' );
    
    if p<alpha
       fprintf( 'Significant increase from state to %s, \n p-value = %f (One-sided Signed Rank test) \n', labels{beh}, p );
    else
       fprintf( 'No significant increase from state to %s, \n p-value = %f (One-sided Signed Rank test) \n', labels{beh}, p );
    end
    
end