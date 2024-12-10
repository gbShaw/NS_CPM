function [Prediction_r, Performance ] = NS_CPM_permutationtest(Data_pre, cfg, Performance)

cfg.cvrp = 1;

if cfg.paralle>0
     [Prediction_r ]= permutation_it_prarallel(Data_pre, cfg);
else
      [Prediction_r ]= permutation_it_single(Data_pre, cfg);
end

% Performance
Prediction_r(1,:) = [Performance.positive.r Performance.negative.r];

Prediction_r(isnan(Prediction_r)) = 0; 

 sorted_prediction_r_pos = sort(Prediction_r(:,1),'descend') ;
 position_pos                   = find(sorted_prediction_r_pos==Performance.positive.r);
 if Performance.positive.r==0
     Performance.Permutation.p_positive = nan;
 else              
     Performance.Permutation.p_positive = position_pos/cfg. PermutationN;
 end
 
 sorted_prediction_r_neg = sort(Prediction_r(:,2),'descend') ;
 position_neg                   = find(sorted_prediction_r_neg==Performance.negative.r);
 if Performance.negative.r==0
     Performance.Permutation.p_negative = nan;
 else              
     Performance.Permutation.p_negative = position_neg/cfg. PermutationN;
 end



        
        
function [Prediction_r ]= permutation_it_prarallel(Data_pre, cfg)

try 
    parpool(cfg.paralle)
catch
    gp = gcp;
    if gp.NumWorkers~=cfg.paralle
        error('parallel open fialed')
    end
end

no_sub = Data_pre.nSub;

cfg.paralle=0;

Prediction_r = zeros(cfg.PermutationN,2);
fprintf('Permutation for %d times\n',cfg.PermutationN)
TBX_parfor_progress(cfg.PermutationN);
parfor it =2:cfg.PermutationN
%      fprintf( '\n Iteration %d out of %d.. ', it , cfg.PermutationN  ) 
%         disp permutation.
  
        DataOrder  = randperm(no_sub);
        [~, Performance ] = NS_CPM_predict(Data_pre, cfg, DataOrder);        
        r_pos = Performance.positive.r;
        r_neg = Performance.negative.r;
        Prediction_r(it,:) = [r_pos r_neg];
        TBX_parfor_progress;
end
TBX_parfor_progress(0);
        
        
function [Prediction_r ]= permutation_it_single(Data_pre, cfg)
no_sub = Data_pre.nSub;

Prediction_r = zeros(cfg.PermutationN,2);
for it =2:cfg.PermutationN
        fprintf( '\n Iteration %d out of %d.. ', it , cfg.PermutationN  ) 
        
        DataOrder  = randperm(no_sub);
        [~, Performance ] = NS_CPM_predict(Data_pre, cfg, DataOrder);
        r_pos = Performance.positive.r;
        r_neg = Performance.negative.r;
        Prediction_r(it,:) = [r_pos r_neg];
end
  
        
        
        