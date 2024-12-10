function Features = NS_CPM_export_Features(Data, CPMcfg)

CPMcfg.nSub = Data.nSub;
FCvcts = Data.FCvcts;
behav  = Data.scores;
covar  = Data.covar;
for i = 1:length(CPMcfg.thresh)
    cfg = CPMcfg;
    cfg.thresh= CPMcfg.thresh(i);
    fprintf('\nExporting consensus features for threshold %f\n',cfg.thresh)
    if CPMcfg.paralle==0
        [  pos_edges{i}, neg_edges{i} ]  = predict_singlecore(FCvcts, behav, covar, cfg ) ;
    else
        [  pos_edges{i}, neg_edges{i} ]  = predict_parallel(FCvcts, behav, covar, cfg ) ;    
    end
    
    Features.consensus95.pos_mask{i} =  revert2mask( pos_edges{i}, Data, 0.95);
    Features.consensus95.neg_mask{i} =  revert2mask( neg_edges{i}, Data, 0.95);
    dlmwrite(  [CPMcfg.path.output filesep sprintf( 'Consensus_Feature_95Overlap_mask_positive_%f.txt' , cfg.thresh) ] , Features.consensus95.pos_mask{i},'delimiter','\t')
    dlmwrite(  [CPMcfg.path.output filesep sprintf( 'Consensus_Feature_95Overlap_mask_negative_%f.txt' , cfg.thresh)],  Features.consensus95.neg_mask{i},  'delimiter','\t')
    
    Features.consensus100.pos_mask{i} =  revert2mask( pos_edges{i}, Data, 1);
    Features.consensus100.neg_mask{i} =  revert2mask( neg_edges{i}, Data, 1);
    dlmwrite(  [CPMcfg.path.output filesep sprintf( 'Consensus_Feature_100Overlap_mask_positive_%f.txt' , cfg.thresh) ], Features.consensus100.pos_mask{i},'delimiter','\t')
    dlmwrite(  [CPMcfg.path.output filesep sprintf( 'Consensus_Feature_100Overlap_mask_negative_%f.txt' , cfg.thresh)],  Features.consensus100.neg_mask{i},  'delimiter','\t')
end
Features.pos_edges = pos_edges;
Features.neg_edges = neg_edges;

save ([CPMcfg.path.output filesep 'NS_CPM_Features.mat'],'Features')


function mask = revert2mask(featureIdx, Data, overlapThresh)
featureIdx = reshape(featureIdx, size(featureIdx,1),[]);
featureIdx = sum(featureIdx,2)>=overlapThresh *size(featureIdx,2);
mask = Data.emptyMask;
mask(Data.feature_idx)=featureIdx;

function writeNiifile(mask,CPMcfg,filename )
v.dim     = [size(mask,1) size(mask,2) 1];
v.fname = [CPMcfg.path.output filesep filename ];
v.mat = [-3,0,0,93;0,3,0,-129;0,0,3,-75;0,0,0,1];
v.dt = [4 0];
spm_write_vol(v,mask)






function [pos_edges, neg_edges ]   =predict_singlecore(FCvcts, behav, covar ,CPMcfg)
fprintf('Current repetition:');
for rp=1:CPMcfg.cvrp
          fprintf('%d', rp);
          [pos_edges(:,:,rp), neg_edges(:,:,rp) ]=  core_featureSelectionCV(FCvcts, behav , covar, CPMcfg  ) ;
end
 fprintf('\n');

 
function  [pos_edges, neg_edges ]  =predict_parallel (FCvcts, behav, covar ,CPMcfg)

TBX_parfor_progress(CPMcfg.cvrp);
 parfor rp=1:CPMcfg.cvrp

        [pos_edges(:,:,rp), neg_edges(:,:,rp)  ] =  core_featureSelectionCV(FCvcts, behav , covar, CPMcfg  ) ;
        TBX_parfor_progress;
 end
TBX_parfor_progress(0);


 
    
function [pos_edges, neg_edges ] =  core_featureSelectionCV(all_vcts, all_behav , covariables  , CPMcfg)

indices = crossvalind('Kfold',1:CPMcfg.nSub, CPMcfg.kfold) ;
Features.pos_edges = ones(size(all_vcts,2),max(indices));
Features.neg_edges =Features.pos_edges ;

for ki = 1: max(indices)
    testInd = (indices == ki); 
    trainInd = ~testInd;

    train_vcts  =all_vcts(trainInd,:);
    train_behav =all_behav(trainInd,:);
    if ~isempty(covariables)
        train_cov = covariables(trainInd,:);        
    else
        train_cov =[];
    end
    test_vct  =all_vcts(testInd,:);


     [pos_edges(:,ki), neg_edges(:,ki)] = core_featureSelection(train_vcts, train_behav, train_cov ,CPMcfg);
          
   
end      





function [pos_edges, neg_edges]= core_featureSelection(vcts, behav, covariables ,CPMcfg )
if ~isempty(covariables )
     [r_mat,p_mat] = partialcorr (vcts,behav ,covariables ) ;
%              [r_mat,p_mat] = corr (vcts,behav ,'type','Spearman'  );
else
     [r_mat,p_mat] = corr (vcts,behav ) ;%,'type','Spearman'  );
end

if  CPMcfg.modelTypeIdx==1
    pos_edges = (r_mat > 0 & p_mat < CPMcfg.thresh);
    neg_edges = (r_mat < 0 & p_mat <  CPMcfg.thresh);
else
    pos_edges =  p_mat < CPMcfg.thresh ;
    neg_edges = pos_edges;
end    
    
    
    
    
    
