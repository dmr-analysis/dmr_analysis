import numpy as np
import pandas as pd
import os
import matplotlib as mlt
import matplotlib.pyplot as plt
#from statistics import NormalDist
#import scipy.stats as st
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix
from scipy import stats

#in-house function
#import plot_dmr
#import dmr_data_analysis
from .plot_dmr import get_parameters_for_analysis
from .dmr_data_analysis import accuracy 

#backend plot is for writing to file,
#mlt.use('Agg')
#exec(open('rank_dmr.py').read())

#GUI backends plot
#mlt.use('TkAgg')
mlt.use('Agg')

EPS=np.finfo(np.float32).eps
def logRegress_feature_score(in_dmr_df):
   ''' use logsitic regression to fit feature vectors from DMR dataframe, then compute a score logReg_score for ranking all DMRs
   '''
   tmp_dmr_df=in_dmr_df.copy()
   tmp_dmr_df['log10_gcb_vs_grpsDist_pval']=-1*np.log10(tmp_dmr_df['gcb_vs_grpsDist_pval']+EPS)
   tmp_dmr_df['log10_tumor_vs_grpsDist_pval']=-1*np.log10(tmp_dmr_df['tumor_vs_grpsDist_pval']+EPS)
  
   #min-max normalization of log10_values
   #max_tval=np.amax([tmp_dmr_df['log10_gcb_vs_grpsDist_pval'].max(),tmp_dmr_df['log10_tumor_vs_grpsDist_pval'].max()])
   #min_tval=np.amin([tmp_dmr_df['log10_gcb_vs_grpsDist_pval'].min(),tmp_dmr_df['log10_tumor_vs_grpsDist_pval'].min()])
   
   max_tval=[tmp_dmr_df['log10_gcb_vs_grpsDist_pval'].max(),tmp_dmr_df['log10_tumor_vs_grpsDist_pval'].max()]
   min_tval=[tmp_dmr_df['log10_gcb_vs_grpsDist_pval'].min(),tmp_dmr_df['log10_tumor_vs_grpsDist_pval'].min()]


   #add minMax normalized log10_pval to the data frame
   tmp_dmr_df['log10_gcb_vs_grpsDist_pval_minMaxNorm']=(tmp_dmr_df['log10_gcb_vs_grpsDist_pval'].to_numpy()-min_tval[0])/(max_tval[0]-min_tval[0])
   tmp_dmr_df['log10_tumor_vs_grpsDist_pval_minMaxNorm']=(tmp_dmr_df['log10_tumor_vs_grpsDist_pval'].to_numpy()-min_tval[1])/(max_tval[1]-min_tval[1])
   #calculate weighted score of 5 feature vectors for a potential dmr ranking
   tmp_dmr_df['dmr_weight_score']=(
                  0.8*(tmp_dmr_df['log10_gcb_vs_grpsDist_pval_minMaxNorm'].to_numpy()+tmp_dmr_df['log10_tumor_vs_grpsDist_pval_minMaxNorm'].to_numpy()) + \
                  0.1*tmp_dmr_df['percent_data_passed_ttest'].to_numpy() + 0.6*tmp_dmr_df['cluster_accuracy'].to_numpy() + \
                  1.0*(tmp_dmr_df['high_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['high_positive_tumor_vs_gcb_percent'].to_numpy()) \
                  )/5

   #jbd check and replace nan in dataframe
   tmp_dmr_df=tmp_dmr_df.fillna(0)
   #print(tmp_dmr_df.isnull().sum())

   len_of_df=tmp_dmr_df.shape[0]
   percentage_cutoff=0
   p_cutoff=0.05
   #use logistic regression to esimate beta parameters for 7 featuer vectors then compute a logReg_score for the final dmr ranking  
   if True:
       #condition 1, 30 for top 400, 100% for top 800, this condition is stable 
       tmp_dmr_df['percent_data_passed_ttest_gt_pval']=tmp_dmr_df['percent_data_passed_ttest'].copy()
       tmp_dmr_df.loc[tmp_dmr_df['percent_data_passed_ttest_gt_pval']>percentage_cutoff,'percent_data_passed_ttest_gt_pval']=1
       tmp_dmr_df.loc[tmp_dmr_df['percent_data_passed_ttest_gt_pval']<=percentage_cutoff,'percent_data_passed_ttest_gt_pval']=0
       x=np.concatenate((tmp_dmr_df['log10_gcb_vs_grpsDist_pval_minMaxNorm'].to_numpy().reshape(-1,1), \
                  tmp_dmr_df['log10_tumor_vs_grpsDist_pval_minMaxNorm'].to_numpy().reshape(-1,1), \
                  tmp_dmr_df['percent_data_passed_ttest_gt_pval'].to_numpy().reshape(-1,1) , tmp_dmr_df['cluster_accuracy'].to_numpy().reshape(-1,1), \
                  (tmp_dmr_df['high_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['high_positive_tumor_vs_gcb_percent'].to_numpy() ).reshape(-1,1), \
                  (tmp_dmr_df['median_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['median_positive_tumor_vs_gcb_percent'].to_numpy() ).reshape(-1,1) ), axis=1)
#                  (tmp_dmr_df['low_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['low_positive_tumor_vs_gcb_percent'].to_numpy() ).reshape(-1,1) ), axis=1)

   if False:
       #condition 2, 31 for top 400, 100% for top 800
       tmp_dmr_df['percent_data_passed_ttest_gt_pval']=tmp_dmr_df['percent_data_passed_ttest'].copy()
       tmp_dmr_df.loc[tmp_dmr_df['percent_data_passed_ttest_gt_pval']>percentage_cutoff,'percent_data_passed_ttest_gt_pval']=1
       tmp_dmr_df.loc[tmp_dmr_df['percent_data_passed_ttest_gt_pval']<=percentage_cutoff,'percent_data_passed_ttest_gt_pval']=0
       x=np.concatenate((tmp_dmr_df[['log10_gcb_vs_grpsDist_pval_minMaxNorm','log10_tumor_vs_grpsDist_pval_minMaxNorm']].max(axis=1).to_numpy().reshape(-1,1), \
                  tmp_dmr_df['percent_data_passed_ttest_gt_pval'].to_numpy().reshape(-1,1) , tmp_dmr_df['cluster_accuracy'].to_numpy().reshape(-1,1), \
                  (tmp_dmr_df['high_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['high_positive_tumor_vs_gcb_percent'].to_numpy() ).reshape(-1,1), \
                  (tmp_dmr_df['median_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['median_positive_tumor_vs_gcb_percent'].to_numpy() ).reshape(-1,1)) , axis=1) #\
#                  (tmp_dmr_df['low_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['low_positive_tumor_vs_gcb_percent'].to_numpy() ).reshape(-1,1) ), axis=1)


#   if False:
#       #condition 3, 27 for top 400, 100% for top 900, this condition is not stable when change the p_cutoff??
#       tmp_dmr_df['percent_data_passed_ttest_gt_pval']=tmp_dmr_df['percent_data_passed_ttest'].copy()
#       tmp_dmr_df.loc[tmp_dmr_df['percent_data_passed_ttest_gt_pval']>percentage_cutoff,'percent_data_passed_ttest_gt_pval']=1
#       tmp_dmr_df.loc[tmp_dmr_df['percent_data_passed_ttest_gt_pval']<=percentage_cutoff,'percent_data_passed_ttest_gt_pval']=0
#  
#       tmp_dmr_df['tumor_or_gcb_grpsDist_ls_pval']=np.zeros(tmp_dmr_df[ 'gcb_vs_grpsDist_pval'].shape)
#       tmp_dmr_df.loc[ ~( (tmp_dmr_df['gcb_vs_grpsDist_pval'] >=p_cutoff) & (tmp_dmr_df[ 'tumor_vs_grpsDist_pval'] >=p_cutoff)),'tumor_or_gcb_grpsDist_ls_pval' ]=1
#
#       x=np.concatenate((tmp_dmr_df['tumor_or_gcb_grpsDist_ls_pval'].to_numpy().reshape(-1,1), \
#                  tmp_dmr_df['percent_data_passed_ttest_gt_pval'].to_numpy().reshape(-1,1) , tmp_dmr_df['cluster_accuracy'].to_numpy().reshape(-1,1), \
#                  (tmp_dmr_df['high_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['high_positive_tumor_vs_gcb_percent'].to_numpy() ).reshape(-1,1), \
#                  (tmp_dmr_df['median_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['median_positive_tumor_vs_gcb_percent'].to_numpy() ).reshape(-1,1), \
#                  (tmp_dmr_df['low_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['low_positive_tumor_vs_gcb_percent'].to_numpy() ).reshape(-1,1) ), axis=1)


   #assume the predefined is_DMR as a binary classifiy for DMR
   y=tmp_dmr_df.is_DMR.copy()
   y.iloc[y=='D']=1
   y.iloc[y=='U']=0
   y=np.array(y,dtype=int)

   model = LogisticRegression(solver='liblinear', C=10,random_state=0).fit(x, y)
   #rand_model=LogisticRegression(solver='liblinear',C=10, random_state=0).fit(rand_x,y)
   # the larger value of C means weaker regularization, or weaker penalization related to high values of b.0 and b.i
   
   #compute final score for all DMR based on estimated beta values
   #log-odds score = log(P/1-P) for class 1
   score=np.matmul(x, model.coef_.T)+ model.intercept_

   #predicted probability for class 1
   predicted_x=model.predict_proba(x)[:,1] 
   print('DMR class ')
   print(model.classes_)
   print('beta coefficient for DMR features :\n gcb_vs_grpsDist, tumor_vs_grpsDist, %_passed_ttest, cluster_accuracy, high_(t-g)%, median_(t-g)%')
   print(model.coef_)
   print('interceptic of Logistic Regression model ')
   print(model.intercept_)
   #print('Prediction accuracy in logReg')
   #print(model.score(x,y))

#   if False:
#      rand_score=np.zeros(score.shape)
#      for i in range(10):
#         rand_x=np.concatenate(( np.random.permutation( tmp_dmr_df['log10_gcb_vs_grpsDist_pval_minMaxNorm'].to_numpy()).reshape(-1,1), \
#                  np.random.permutation(tmp_dmr_df['log10_tumor_vs_grpsDist_pval_minMaxNorm'].to_numpy()).reshape(-1,1), \
#                  np.random.permutation(tmp_dmr_df['percent_data_passed_ttest_gt_pval'].to_numpy()).reshape(-1,1) ,np.random.permutation(tmp_dmr_df['cluster_accuracy'].to_numpy()).reshape(-1,1), \
#                  np.random.permutation((tmp_dmr_df['high_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['high_positive_tumor_vs_gcb_percent'].to_numpy() )).reshape(-1,1), \
#                  np.random.permutation((tmp_dmr_df['median_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['median_positive_tumor_vs_gcb_percent'].to_numpy() )).reshape(-1,1) , \
#                  np.random.permutation((tmp_dmr_df['low_negative_tumor_vs_gcb_percent'].to_numpy() + tmp_dmr_df['low_positive_tumor_vs_gcb_percent'].to_numpy() )).reshape(-1,1) ), axis=1)
#        
#        # rand_x=np.concatenate(( np.random.random((1,len_of_df)).reshape(-1,1), \
#        #                       np.random.random((1,len_of_df)).reshape(-1,1), \
#        #                       np.random.random((1,len_of_df)).reshape(-1,1), np.random.random((1,len_of_df)).reshape(-1,1), \
#        #                       np.random.random((1,len_of_df)).reshape(-1,1), np.random.random((1,len_of_df)).reshape(-1,1),np.random.random((1,len_of_df)).reshape(-1,1) ),axis=1)
#
#         rand_model=LogisticRegression(solver='liblinear',C=10, random_state=0).fit(rand_x,y)
#         rand_score=np.matmul(rand_x,rand_model.coef_.T)+rand_score
#      rand_score=rand_score/10
#      tmp_dmr_df['logReg_randScore']=list(rand_score)
#      tmp_dmr_df['logReg_randScore']=tmp_dmr_df.logReg_randScore.astype(float)


   #add new scorse to data frame
   tmp_dmr_df['logReg_score']=list(score)
   tmp_dmr_df['logReg_score']=tmp_dmr_df.logReg_score.astype(float)

   tmp_dmr_df['logReg_predicted_dmr']=list(predicted_x)
   tmp_dmr_df['logReg_predicted_dmr']=tmp_dmr_df.logReg_predicted_dmr.astype(float)

   #calculate prediction accuracy
   conf_matrix=confusion_matrix(y, model.predict(x))
   #prediction_accuracy=dmr_data_analysis.accuracy(conf_matrix)

   return tmp_dmr_df, conf_matrix

if __name__ == '__main__':

  in_folder='out/DMR/chrY/plots/'
  in_dmr_file='chrY_maxDist_250_minSize_5_DMR_clusterAccuracy_gt_0.5_miniMethyChange_gt_0.07_0.15_0.2_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0_1756_all'

  #in_folder='out/DMR/chr1/plots/'
  #in_dmr_file='chr1_maxDist_250_minSize_5_DMR_clusterAccuracy_gt_0.5_miniMethyChange_gt_0.07_0.15_0.2_high_miniPercentChange_gt_0.0001_Pcutoff_0.05_isSmooth_2_isModTest_0_93816_all'



  in_dmr_file=os.path.join(in_folder, in_dmr_file)
  print(in_dmr_file)
  in_dmr_df=pd.read_csv(in_dmr_file,sep='\t')
  print(in_dmr_df.shape)

  #get results parameters
  data_start_col=3
  in_data_df=pd.DataFrame()
  gcb_col_idx, tumor_col_idx, percent_cutoff, low_median_high_cutoff, in_chrm, \
      max_read_length, mini_block_size, P_cutoff, isSmooth, is_modT, mini_percentage_cutoff = get_parameters_for_analysis(data_start_col, in_data_df, in_dmr_file)

  tmp_dmr_df, conf_matrix =logRegress_feature_score(in_dmr_df)
  prediction_accuracy=accuracy(conf_matrix)
  print('Prediction accuracy and confusion matrix')
  print(prediction_accuracy)
  print(conf_matrix) 
  
  tmp_dmr_df=tmp_dmr_df.sort_values(by='logReg_score')
  tmp_dmr_df=tmp_dmr_df.reset_index(drop=True)
  selected_df=tmp_dmr_df[(tmp_dmr_df.is_DMR=='D') & (tmp_dmr_df.logReg_predicted_dmr==0)]
 
  #conver to zscore
  #tmp_dmr_df['logReg_score2zscore']=stats.zscore(tmp_dmr_df['logReg_score'])
  #tmp_dmr_df['pval_zscore']=stats.norm.sf(abs(tmp_dmr_df.logReg_score2zscore))
 
  #random permutation test
  #num_of_draw=int(np.ceil(tmp_dmr_df.shape[0]*0.9))
  #tmp_score=tmp_dmr_df.logReg_score.copy()
  #tmp_index=tmp_score.index.tolist()
  #tmp_expt_pval=[]
  #for idx,rr in  tmp_dmr_df.iterrows():
  #    tmp_index2=set(tmp_index)-set([idx])
  #    tmp_idx=np.random.permutation(list(tmp_index2))
  #    tmp_score2=tmp_score[tmp_idx[0:num_of_draw]]
  #    tmp_expt_pval.append(np.where(rr['logReg_score']>tmp_score2)[0].shape[0]/num_of_draw)
  
  #mini_logReg_score= tmp_dmr_df.loc[tmp_dmr_df[tmp_dmr_df.is_DMR=='D' ].index.min(),].logReg_score

  #plot score
  #nn, bins, patches=plt.hist(x=tmp_dmr_df.logReg_score,bins=[np.floor(tmp_dmr_df.logReg_score.min()-1),0, 1, 2,4,8,16,32,np.ceil(tmp_dmr_df.logReg_score.max()+5)])   
  #plt.show()

  #export data
  out_file=in_dmr_file + "_test"
  print(out_file)
  tmp_dmr_df.to_csv(out_file, index=False, sep='\t')








