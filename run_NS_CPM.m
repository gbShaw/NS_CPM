function [ALLPrediction, ALLPerformance ]= run_NS_CPM(Data,CPMcfg)
% Data
% CPMcfg
disp('running...')


%% open parpool
if CPMcfg.paralle>0
    try 
        parpool(CPMcfg.paralle)
    catch
        gp = gcp;
        if gp.NumWorkers~=CPMcfg.paralle
            error('parallel open fialed')
        end
    end
end


%% compute
Datapre  = Data;
Datapre  = rmfield(Datapre,{'FCmat' 'confirm'});
cfgpre   = CPMcfg;
for i = 1:length(CPMcfg.thresh)
    cfgpre.thresh = CPMcfg.thresh(i); 
    
    %==  perdiction
    fprintf('Perform prediction: thresh = %f\n',cfgpre.thresh); 
    [ALLPrediction{i}, ALLPerformance{i} ] = NS_CPM_predict(Datapre, cfgpre);   % predict    
    fprintf(['Performance: \n',...
                        '           positive :  r =%.3f; p = %.3f \n' ,...
                        '           negetive :  r =%.3f; p = %.3f \n' ],...
                        ALLPerformance{i} .positive.r, ALLPerformance{i} .positive.p,...
                        ALLPerformance{i} .negative.r,ALLPerformance{i} .negative.p)                            
    
    %== permutation
    warning off
    if CPMcfg.isPermutation==1
        [ALLPrediction{i}.Permutation_r, ALLPerformance{i} ] = NS_CPM_permutationtest(Datapre, cfgpre, ALLPerformance{i} );
    end
    warning on  %
    
    %== Export consensus features
    if CPMcfg.isExportConsesus==1
        NS_CPM_export_Features(Data, cfgpre);
    end

    %== save txt format results file
    txtReport(cfgpre, ALLPerformance{i})
end


save([CPMcfg.path.output filesep 'NS_CPM_results.mat'], 'ALLPerformance','ALLPrediction','CPMcfg')
disp('Congratulation, work finished!')

function txtReport(CPMcfg,Performance)
fid = fopen([CPMcfg.path.output filesep 'NS_CPM_reports_' num2str(CPMcfg.thresh) '.txt'],'w');
fprintf(fid,'Reports for feather selection thresh : %f \n',CPMcfg.thresh);
fprintf(fid,'\n========================');
fprintf(fid,'Prediction:\n\n  For positive features:\n    r=%0.5f\n    p=%0.5f\n', ...
                Performance.positive.r,Performance.positive.p);
fprintf(fid,'\n  For negative features:\n    r=%0.5f\n    p=%0.5f\n', ...
                Performance.negative.r,Performance.negative.p);            

fprintf(fid,'\n========================');
fprintf(fid,'Permutation test:\n\n  For positive features:\n    p=%0.5f\n', ...
                Performance.Permutation.p_positive );
fprintf(fid,'\n  For negative features:\n    p=%0.5f\n', ...
                Performance.Permutation.p_negative);                  
fprintf(fid,'\nTips:  Only the results of the CPM model distinguish between positive and negative \n'       );
    
fclose(fid);


 