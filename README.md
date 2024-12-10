# Brain Network Machine Learning Toolkit NS_CPM

NS_CPM is a GUI toolkit developed for machine learning analysis of connectivity data. It uses univariate linear models, lasso regression, or relevance vector regression (RVR) to predict continuous behavioral variables and has established a widely used k-fold cross-validation framework and permutation testing. NS_CPM can output overlapping features for feature generalization. It provides concise visual operations as well as highly flexible scripting options. All its functions support parallel computation, allowing for rapid data analysis on high-performance servers, and it has been tested successfully in MATLAB 2010a and above (WIN/LINUX).

## Feature Input

Connectivity data serves as common features for machine learning in the field of cognitive neuroscience, which was the initial design focus of this tool. Connectivity networks have a solid research foundation and potential, holding significant research value for an extended period in the future. Both functional connectivity and structural connectivity data can be used as feature inputs. The author may also incorporate other modalities into NS_CPM in the future.

## Model Construction

CPM (Shen et al., 2017). As a method based on univariate linear models, it results in fewer false positives and stronger interpretability while being widely applied.
LASSO Regression. LASSO regularization shares a key similarity with connectivity networks: sparsity. It is commonly used for feature extraction in machine learning but can also be utilized for model construction.
Elastic-net. The script mode also supports elastic net regression by adjusting the ratio of L1 to L2 regularization to modify features within the model. For simplicity's sake, this model has not yet been included in the GUI.
RVR. Relevance Vector Regression algorithm.

## Model Evaluation

NS_CPM provides a cross-validation framework. Leave-one-out methods may lead to overestimation of models and are not recommended for cognitive neuroscience model evaluation (Poldrack et al., 2020). Currently, only k-fold cross-validation is incorporated into NS_CPM.

Permutation testing can further evaluate models by randomly permuting labels followed by repeated model construction to obtain distributions of performance metrics; p-value is defined as the percentage of permutation scores that exceed those from original classifications.

## Result Output

NS_CPM uses correlation coefficients between predicted values and true values as evaluation metrics; p-values are calculated through permutation tests to enhance result credibility. Additionally, NS_CPM can output overlapping features for feature generalization.





**This toolkit is a graphical version of the code in this article:**

Xiao, Z., Chen, Z., Chen, W., Gao, W., He, L., Wang, Q., Lei, X., Qiu, J., Feng, T., Chen, H., Turel, O., Bechara, A., & He, Q. (2022). Maladaptive changes in delay discounting in males during the COVID-19 pandemic: the predictive role of functional connectome. Cerebral Cortex, 32(20), 4605–4618.

https://academic.oup.com/cercor/article/32/20/4605/6512152 





 

## References

### Published papers:

- Xiao, Z., Chen, Z., Chen, W., Gao, W., He, L., Wang, Q., Lei, X., Qiu, J., Feng, T., Chen, H., Turel, O., Bechara, A., & He, Q. (2022). Maladaptive changes in delay discounting in males during the COVID-19 pandemic: the predictive role of functional connectome. Cerebral Cortex, 32(20), 4605–4618.

- Chen, Z., Tang, Y., Liu, X., Li, W., Hu, Y., Hu, B., Xu, T., Zhang, R., Xia, L., Zhang, J.-X., Xiao, Z., Chen, J., Feng, Z., Zhou, Y., He, Q., Qiu, J., Lei, X., Chen, H., Qin, S., & Feng, T. (2024). Edge-centric connectome-genetic markers of bridging factor to comorbidity between depression and anxiety. Nature Communications, 15(1), 10560.



### Model references:

- Shen, X., Finn, E. S., Scheinost, D., Rosenberg, M. D., Chun, M. M., Papademetris, X., & Constable, R. T. (2017). Using connectome-based predictive modeling to predict individual behavior from brain connectivity. nature protocols, 12(3), 506-518.

- Poldrack, R. A., Huckins, G., & Varoquaux, G. (2020). Establishment of best practices for evidence for prediction: a review. JAMA psychiatry, 77(5), 534-540.


## Release note
1.0.31_beta 20211223
1.fix bug: missing TBX_parfor_progress

1.0.3_beta 20211214
1.add online update
2.add crossvalind as subfunction to avoid conflicts


1.0.2_beta 20211130
1.Add elastic-net methods for batch mode
2.fix bug: missing TBX_parfor_progress
3.fix bug: missing crossvalind


1.0.1_beta
add support for matlab2010a