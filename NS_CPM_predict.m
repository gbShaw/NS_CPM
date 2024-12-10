function [predicted, Performance ] = NS_CPM_predict(Data_pre, cfg, DataOrder)


if nargin>2
    Data_pre.scores=Data_pre.scores(DataOrder);
end

cfg.nSub = Data_pre.nSub;
FCvcts = Data_pre.FCvcts;
behav  = Data_pre.scores;
covar  = Data_pre.covar;

if cfg.paralle==0
    [predicted, Performance ]  = predict_singlecore(FCvcts, behav, covar, cfg ) ;
else
    [predicted, Performance ]  = predict_parallel(FCvcts, behav, covar, cfg ) ;    
end
Performance.Permutation.p_negative=nan;
Performance.Permutation.p_positive=nan;



function [predicted, Performance ] =predict_singlecore(FCvcts, behav, covar ,cfg)
predicted_neg = zeros(length(behav),cfg.cvrp);   predicted_pos = predicted_neg;
% cfg.cvrp
% fprintf('Current repetition: ')
for rp=1:cfg.cvrp
%     fprintf('%d', rp);
     if cfg.modelTypeIdx==1
         [predicted_pos(:,rp), predicted_neg(:,rp) ] =  core_predict_cpm(FCvcts, behav , covar, cfg  ) ;
     else
         [predicted_pos(:,rp), predicted_neg(:,rp) ] =  core_predict_lasso (FCvcts, behav , covar, cfg  ) ;
     end
end
%  fprintf('\n');
predicted.positive = predicted_pos;
predicted.negative = predicted_neg;
 
[Performance.positive.r,Performance.positive.p ] = corr(mean(predicted_pos,2), behav) ;
[Performance.negative.r,Performance.negative.p ] = corr(mean(predicted_neg,2), behav) ;



function [predicted, Performance ] =predict_parallel (FCvcts, behav, covar ,cfg)


predicted_neg = zeros(length(behav),cfg.cvrp);   predicted_pos = predicted_neg;
modelTypeIdx = cfg.modelTypeIdx;
 parfor rp=1:cfg.cvrp
     if modelTypeIdx==1
         [predicted_pos(:,rp), predicted_neg(:,rp) ] =  core_predict_cpm(FCvcts, behav , covar, cfg  ) ;
     else
         [predicted_pos(:,rp), predicted_neg(:,rp) ] =  core_predict_lasso (FCvcts, behav , covar, cfg  ) ;
     end
 end
predicted.positive = predicted_pos;
predicted.negative = predicted_neg;
 
[Performance.positive.r,Performance.positive.p ] = corr(mean(predicted_pos,2), behav) ;
[Performance.negative.r,Performance.negative.p ] = corr(mean(predicted_neg,2), behav) ;



function [pos_edges, neg_edges]= core_featureSelection(vcts, behav, covariables ,cfg )
if ~isempty(covariables )
     [r_mat,p_mat] = partialcorr (vcts,behav ,covariables ) ;
%              [r_mat,p_mat] = corr (vcts,behav ,'type','Spearman'  );
else
     [r_mat,p_mat] = corr (vcts,behav ) ;%,'type','Spearman'  );
end


if  cfg.modelTypeIdx==1
    pos_edges = (r_mat > 0 & p_mat < cfg.thresh);
    neg_edges = (r_mat < 0 & p_mat <  cfg.thresh);
else
    pos_edges =  p_mat < cfg.thresh ;
    neg_edges = pos_edges;
end
    
    
function [predicted_pos, predicted_neg ] =  core_predict_cpm(all_vcts, all_behav , covariables  , cfg)
warning off
predicted_pos = zeros(cfg.nSub,1); predicted_neg = predicted_pos ;
indices = crossvalind('Kfold',1:cfg.nSub, cfg.kfold) ;
for ki = 1: max(indices)
    testInd = (indices == ki); 
    trainInd = ~testInd;

    train_vcts  =all_vcts(trainInd,:);
    train_behav =all_behav(trainInd,:);
    if ~isempty(covariables)
        train_cov = covariables(trainInd,:);    
    else
        train_cov = [];    
    end
    test_vct  =all_vcts(testInd,:);


     [pos_edges, neg_edges] = core_featureSelection(train_vcts, train_behav, train_cov ,cfg);

    if sum(pos_edges+neg_edges)==0,error('no correlation feature found');end

    train_sum = sum( train_vcts(: , pos_edges) ,2)   ;
    test_sum  = sum( test_vct( :, pos_edges ) , 2 )  ;

    Beta = polyfit(train_sum, train_behav,1) ;
    predicted_pos(testInd,1)= Beta(1)*test_sum + Beta(2) ;

    train_sum = sum( train_vcts(: , neg_edges) ,2)   ;
    test_sum  = sum( test_vct( :, neg_edges ) , 2 )  ;

    Beta = polyfit(train_sum, train_behav,1) ;
    predicted_neg(testInd,1) = Beta(1)*test_sum + Beta(2) ;                        
%         FitInfo.Beta = Beta ;
end        
warning on
    
 function [predicted_pos, predicted_neg ] =  core_predict_lasso(all_vcts, all_behav , covariables  , cfg)
indices = crossvalind('Kfold',1:cfg.nSub, cfg.kfold) ;
predicted_pos = zeros(cfg.nSub,1); predicted_neg = predicted_pos ;
    for ki = 1: max(indices)
        testInd = (indices == ki); 
        trainInd = ~testInd;
        
        train_vcts  =all_vcts(trainInd,:);
        train_behav =all_behav(trainInd,:);
        train_cov = covariables(trainInd,:);        
        test_vct  =all_vcts(testInd,:);

        
         [pos_edges, neg_edges] = core_featureSelection(train_vcts, train_behav, train_cov ,cfg.lassoAlpha);
        
        if sum(pos_edges+neg_edges)==0,error('no correlation feature found');end
        
        [B,FitInfo] = lasso(train_vcts(: , pos_edges) ,train_behav,'CV', 10 ,'Alpha',1  );        
        b  = B(:,FitInfo.IndexMinMSE);
        b0 = FitInfo.Intercept(FitInfo.IndexMinMSE);
        predicted_pos(testInd,1)  = sum(test_vct( :, pos_edges)' .* b , 1) + b0 ;     

       
        predicted_neg(testInd,1)  = predicted_pos(testInd,1)  ;

    end    
 
    
    
function [tInd,eInd] = crossvalind(method,N,varargin)
%CROSSVALIND generates cross-validation indices
%
%   INDICES = CROSSVALIND('Kfold',N,K) returns randomly generated indices
%   for a K-fold cross-validation of N observations. INDICES contains equal
%   (or approximately equal) proportions of the integers 1 through K that
%   define a partition of the N observations into K disjoint subsets.
%   Repeated calls return different randomly generated partitions. K
%   defaults to 5 when omitted. In K-fold cross-validation, K-1 folds are
%   used for training and the last fold is used for evaluation. This
%   process is repeated K times, leaving one different fold for evaluation
%   each time.
%
%   [TRAIN,TEST] = CROSSVALIND('HoldOut',N,P) returns logical index vectors
%   for cross-validation of N observations by randomly selecting P*N
%   (approximately) observations to hold out for the evaluation set. P must
%   be a scalar between 0 and 1. P defaults to 0.5 when omitted,
%   corresponding to holding 50% out. Using holdout cross-validation within
%   a loop is similar to K-fold cross-validation one time outside the loop,
%   except that non-disjointed subsets are assigned to each evaluation.
%
%   [TRAIN,TEST] = CROSSVALIND('LeaveMOut',N,M), where M is an integer,
%   returns logical index vectors for cross-validation of N observations by
%   randomly selecting M of the observations to hold out for the evaluation
%   set. M defaults to 1 when omitted. Using LeaveMOut cross-validation
%   within a loop does not guarantee disjointed evaluation sets. Use K-fold
%   instead.
%
%   [TRAIN,TEST] = CROSSVALIND('Resubstitution',N,[P,Q]) returns logical
%   index vectors of indices for cross-validation of N observations by
%   randomly selecting P*N observations for the evaluation set and Q*N
%   observations for training. Sets are selected in order to minimize the
%   number of observations that are used in both sets. P and Q are scalars
%   between 0 and 1. Q=1-P corresponds to holding out (100*P)%, while P=Q=1
%   corresponds to full resubstitution. [P,Q] defaults to [1,1] when omitted.
%
%   [...] = CROSSVALIND(METHOD,GROUP,...) takes the group structure of the
%   data into account. GROUP is a grouping vector that defines the class for
%   each observation. GROUP can be a numeric vector, a string array, or a
%   cell array of strings. The partition of the groups depends on the type
%   of cross-validation: For K-fold, each group is divided into K subsets,
%   approximately equal in size. For all others, approximately equal
%   numbers of observations from each group are selected for the evaluation
%   set. In both cases the training set will contain at least one
%   observation from each group.
%
%   [...] = CROSSVALIND(METHOD,GROUP,...,'CLASSES',C) restricts the
%   observations to only those values specified in C.  C can be a numeric
%   vector, a string vector, or a cell array of character vectors, but it is
%   of the same form as GROUP. If one output argument is specified, it will
%   contain the value 0 for observations belonging to excluded classes. If
%   two output arguments are specified, both will contain the logical value
%   false for observations belonging to excluded classes.
%
%   [...] = CROSSVALIND(METHOD,GROUP,...,'MIN',MIN) sets the minimum number
%   of observations that each group has in the training set. MIN defaults
%   to 1. Setting a large value for MIN can help to balance the training
%   groups, but adds partial resubstitution when there are not enough
%   observations. You cannot set MIN when using K-fold cross-validation.
%
%   Examples:
%
%      % Create a 10-fold cross-validation to compute classification error.
%      load fisheriris
%      indices = crossvalind('Kfold',species,10);
%      cp = classperf(species);
%      for i = 1:10
%          test = (indices == i); train = ~test;
%          class = classify(meas(test,:),meas(train,:),species(train,:));
%          classperf(cp,class,test)
%      end
%      cp.ErrorRate
%
%      % Approximate a leave-one-out prediction error estimate.
%      load carbig
%      x = Displacement; y = Acceleration;
%      N = length(x);
%      sse = 0;
%      for i = 1:100
%          [train,test] = crossvalind('LeaveMOut',N,1);
%          yhat = polyval(polyfit(x(train),y(train),2),x(test));
%          sse = sse + sum((yhat - y(test)).^2);
%      end
%      CVerr = sse / 100
%
%      % Divide cancer data 60/40 without using the 'Benign' observations.
%      % Assume groups are the true labels of the observations.
%      labels = {'Cancer','Benign','Control'};
%      groups = labels(ceil(rand(100,1)*3));
%      [train,test] = crossvalind('holdout',groups,0.6,'classes',...
%          {'Control','Cancer'});
%      sum(test) % Total groups allocated for testing
%      sum(train) % Total groups allocated for training
%
%   See also CLASSPERF, CLASSIFY, GRP2IDX, FITCKNN, CLASSIFICATIONSVM.

%   References:
%   [1] Hastie, T. Tibshirani, R, and Friedman, J. (2001) The Elements of
%       Statistical Learning, Springer, pp. 214-216.
%   [2] Theodoridis, S. and Koutroumbas, K.  (1999) Pattern Recognition,
%       Academic Press, pp. 341-342.

% Copyright 2003-2016 The MathWorks, Inc.
% edited by NS_CPM, 20211214 remove convertStringsToChars

% set defaults
% if nargin > 0
% %     method = convertStringsToChars(method);
% end
% 
% if nargin > 1
% %     N = convertStringsToChars(N);
% end
% 
% if nargin > 2
% %     [varargin{:}] = convertStringsToChars(varargin{:});
% end

classesProvided = false;
MG = 1;   % default for minimum number of observations for every training group
P = 0.5;  % default value for holdout
K = 5;    % default value for Kfold
M = 1;    % default value for leave-M-out
Q = [1 1];% default value for resubstitution

% get and validate the method (first input)
if ischar(method) && size(method,1)==1
    validMethods = {'holdout','kfold','resubstitution','leavemout'};
    method = strmatch(lower(method),validMethods); 
    if isempty(method)
        error(message('bioinfo:crossvalind:NotValidMethod'))
    end
    method = validMethods{method};
else
    error(message('bioinfo:crossvalind:NotValidTypeForMethod'))
end

if nargout>1 && isequal(method,'kfold')
    error(message('bioinfo:crossvalind:TooManyOutputArgumentsForKfold'))
end

% take P,K,Q, or M if provided by the third input (first varargin) and
% validate it
if numel(varargin) && isnumeric(varargin{1})
    S = varargin{1};
    varargin(1)=[];
    switch method
        case 'holdout'
            if numel(S)==1 && S>0 && S<1
                P = S;
            else
                error(message('bioinfo:crossvalind:InvalidThirdInputP'));
            end
        case 'kfold'
            if  numel(S)==1 && S>=1
                K = round(S);
            else
                error(message('bioinfo:crossvalind:InvalidThirdInputK'));
            end
        case 'leavemout'
            if  numel(S)==1 && S>=1
                M = round(S);
            else
                error(message('bioinfo:crossvalind:InvalidThirdInputM'));
            end
        case 'resubstitution'
            if numel(S)==2 && all(S>0) && all(S<=1)
                Q = S(:);
            else
                error(message('bioinfo:crossvalind:InvalidThirdInputQ'));
            end
    end %switch
end

% read optional paired input arguments in
if numel(varargin)
    if rem(numel(varargin),2)
        error(message('bioinfo:crossvalind:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'classes','min'};
    for j=1:2:numel(varargin)
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:crossvalind:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:crossvalind:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % classes
                    classesProvided = true;
                    classes = pval;
                case 2 % min
                    MG = round(pval(1));
                    if MG<0
                        error(message('bioinfo:crossvalind:NotValidMIN'))
                    end
            end
        end
    end
end

if isscalar(N) && isnumeric(N)
    if N<1 || N~=floor(N)
        error(message('bioinfo:crossvalind:NNotPositiveInteger'));
    end
    group = ones(N,1);
else
    [group, groupNames] = grp2idx(N); % at this point group is numeric only
    N = numel(group);
end

if classesProvided
    orgN = N;
    % change classes to same type as groups
    [~,classes]=grp2idx(classes);
    validGroups = intersect(classes,groupNames);
    if isempty(validGroups)
        error(message('bioinfo:crossvalind:EmptyValidGroups'))
    end
    selectedGroups = ismember(groupNames(group),validGroups);
    group = grp2idx(group(selectedGroups)); % group idxs are reduced to only the sel groups
    N = numel(group);     % the new size of the reduced vector
end

nS = accumarray(group(:),1);
if min(nS)<MG
    error(message('bioinfo:crossvalind:MissingObservations'))
end

switch method
    case {'leavemout','holdout','resubstitution'}
        switch method
            case 'leavemout'
                % number of samples for holdout in every group
                nSE = repmat(M,numel(nS),1);
                % at least there is MG sample(s) for training in every group
                nST = max(nS-nSE,MG);
            case 'holdout'
                % computes the number of samples for holdout in every group
                nSE = floor(nS*P);
                % at least there is MG sample(s) for training in every group
                nST = max(nS-nSE,MG);
            case 'resubstitution'
                % computes the number of samples for training and evaluation
                nSE = floor(nS*Q(1));
                nST = floor(nS*Q(2));
                % at least there is MG sample(s) for training in every group
                nST = max(nST,MG);
        end
        % Initializing the outputs
        tInd = false(N,1);
        eInd = false(N,1);
        % for every group select randomly the samples for both sets
        for g = 1:numel(nS)
            h = find(group==g);
            randInd = randperm(nS(g));
            tInd(h(randInd(1:nST(g))))=true;
            eInd(h(randInd(end-nSE(g)+1:end)))=true;
        end
    case 'kfold'
        tInd = zeros(N,1);
        for g = 1:numel(nS)
            h = find(group==g);
            % compute fold id's for every observation in the  group
            q = ceil(K*(1:nS(g))/nS(g));
            % and permute them to try to balance among all groups
            pq = randperm(K);
            % randomly assign the id's to the observations of this group
            randInd = randperm(nS(g));
            tInd(h(randInd))=pq(q);
        end
end

if classesProvided
    if isequal(method,'kfold')
        temp = zeros(orgN,1);
        temp(selectedGroups) = tInd;
        tInd = temp;
    else
        temp = false(orgN,1);
        temp(selectedGroups) = tInd;
        tInd = temp;
        temp = false(orgN,1);
        temp(selectedGroups) = eInd;
        eInd = temp;
    end
end






    
    
    
    
    
    
