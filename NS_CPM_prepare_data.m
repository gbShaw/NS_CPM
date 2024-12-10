function Data = NS_CPM_prepare_data(CPMcfg)
% input:  CPMcfg > configuration struct
% output: Data


%% read FC file list
f = dir([CPMcfg.path.FC filesep '*.txt']);
fileNames_FC = cat(1,{f(:).name})';
[fileNames_FC,ia] = sortrows(fileNames_FC);
f = f(ia);




%% read behav
try
    Data.scores = load(CPMcfg.path.scores);
    fileNames_score = fileNames_FC;
catch    
    behav = textread(CPMcfg.path.scores,'%s','delimiter', '\t' ) ;
    behav = reshape(behav,2,[])';
    behav = sortrows(behav,1);
    fileNames_score = behav(:,1);    
    Data.scores   = cellfun(@str2double,behav(:,2));
end

%% read covariable
if logical(exist(CPMcfg.path.covar,'file'))
    try
        Data.covar = load(CPMcfg.path.covar);
        fileNames_covar=fileNames_FC;    
    catch    
        behav = textread(CPMcfg.path.covar,'%[^\n]', 1 )  ;
        nc = length(strsplit(behav{1})) ;
        behav = textread(CPMcfg.path.covar,'%s','delimiter', '\t' ) ;    
        behav = reshape(behav,nc,[])';
        behav = sortrows(behav,1);
        fileNames_covar = behav(:,1);    
        Data.covar   = cellfun(@str2double,behav(:,2:end));         
    end
else
    Data.covar=[];
    fileNames_covar = fileNames_FC;
end

%
if ~isequal(fileNames_covar,fileNames_score)
    fileNames_covar
    fileNames_score
    error('covariable and score not consistent')
end


%% read FC

[~, ia,ib ] = intersect(fileNames_score,fileNames_FC );
if ~isequal(sortrows(fileNames_score), fileNames_FC(ib))
    error('Input not consistent!')
end
f = f(ib);
Data.scores = Data.scores(ia);
if ~isempty(Data.covar)
    Data.covar  = Data.covar (ia,:);
end
fileNames_FC = cat(1,{f(:).name})';
fileNames_score = fileNames_score(ia);


for nSub =1 :length(f)
    fprintf('Reading data : %s\n', f(nSub).name);
    Data.FCmat(:,:,nSub) = load( [CPMcfg.path.FC filesep f(nSub).name]);
end

%% reshape data
FC_mat       = squeeze( Data.FCmat(:,:,1));
Data.emptyMask = FC_mat*0;
Data.feature_idx  = find(triu(ones(size(FC_mat)),1));
feature_delete  = ~(triu(ones(size(FC_mat)),1));
Data.FCvcts  = reshape(Data.FCmat, [], nSub)';
Data.FCvcts(:,feature_delete) =[];
Data.nSub  = nSub;

% m = Data.emptyMask;
% size(m(Data.feature_idx))
% size(Data.FCvcts(1,:))
% m(Data.feature_idx) = Data.FCvcts(1,:);
% whos m
% subplot(1,3,1);
% imagesc(FC_mat)
% subplot(1,3,2);
% imagesc(m);
% subplot(1,3,3);
% imagesc(FC_mat-m);


%% confirm
Data.confirm = [fileNames_FC num2cell(Data.scores) num2cell(Data.covar) ];







