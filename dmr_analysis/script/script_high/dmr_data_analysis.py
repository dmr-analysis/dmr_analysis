#Data analysis functions in DMR pipeline
import numpy as np
import pandas as pd
import logging
import os
import warnings
warnings.filterwarnings("ignore")

from .others.ranksum import ranksum 
from scipy.stats import ks_2samp
from scipy.stats import ttest_ind
from scipy.stats import zscore
from scipy.stats import t as t_pdf
from scipy.stats import norm 

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import confusion_matrix
from statsmodels.stats.multitest import fdrcorrection

from sklearn.manifold import TSNE
#import json
from sklearn.cluster import KMeans

#in-house function
from .dmr_plots import plot_raw_tumor_data, plot_raw_gcb_data, plot_3D_PCA, plot_2D_tSNE, plot_smoothed_data, plot_hist_of_difference_in_groups, plot_hist_for_block_length_and_data_size
 
#gloabl variables for mininum values
EPS=np.finfo(np.float32).eps
#np.seterr(all='raise')
np.seterr(all='ignore')

def too_much_zeros(rd_smooth1, logger_n):
  #check column samples with zeros
  tmp_df1=pd.DataFrame(data=np.array(rd_smooth1).T)
  #added may jbw
  tmp_sample_size=tmp_df1.shape[1]
  non_zero_tmp_df1=tmp_df1.loc[:,tmp_df1.all()]
  non_zero_cols=non_zero_tmp_df1.shape[1]
  if non_zero_cols/tmp_sample_size < 0.3:
     logger_n.info('Too much zeros in datagrame, number of nonzero columns is %g', non_zero_cols)
     return True
  else:
     return False 

def check_zeros_in_df(rd_smooth1, rd_smooth2,logger_n):
  #check zerors in dataframe
  tmp_df1=pd.DataFrame(data=np.array(rd_smooth1).T)
  tmp_df2=pd.DataFrame(data=np.array(rd_smooth2).T)
  non_zero_tmp_df1=tmp_df1.loc[:,tmp_df1.all()]
  non_zero_tmp_df2=tmp_df2.loc[:,tmp_df2.all()]
  if non_zero_tmp_df1.shape[1]+ non_zero_tmp_df2.shape[1]<3 and non_zero_tmp_df1.shape[1]>=1:
     pop_col_name=tmp_df1.loc[:,tmp_df1.all()].columns.to_list()[0]
     #put random small value to avoid underflow
     rd_small_values=np.linspace(10**(-3),10**(-1),tmp_df1.shape[0]*10)
     rd_small_values2=np.random.RandomState(seed=2).permutation(rd_small_values)[0:tmp_df1.shape[0]]
     pop_col_data= np.abs(tmp_df1.pop(pop_col_name) -  rd_small_values2)
     insert_col_idx=int(np.ceil(tmp_df1.shape[1]*0.6))
     tmp_df1.insert(insert_col_idx,pop_col_name,pop_col_data)
     al_smooth=np.concatenate((tmp_df1,tmp_df2),axis=1).T
     logger_n.info('Too much zeros in datagrame, reorder group1 sample columns with data to the middle %s, %s', pop_col_name, str(insert_col_idx) )
  elif non_zero_tmp_df1.shape[1]+ non_zero_tmp_df2.shape[1]<3 and non_zero_tmp_df2.shape[1]>=1:
     pop_col_name=tmp_df2.loc[:,tmp_df2.all()].columns.to_list()[0]
     #put random small value to avoid underflow
     rd_small_values=np.linspace(10**(-3),10**(-1),tmp_df2.shape[0]*10)
     rd_small_values2=np.random.RandomState(seed=2).permutation(rd_small_values)[0:tmp_df2.shape[0]]
     pop_col_data= np.abs(tmp_df2.pop(pop_col_name) - rd_small_values2)
     insert_col_idx=int(np.ceil(tmp_df2.shape[1]*0.6))
     tmp_df2.insert(insert_col_idx,pop_col_name,pop_col_data)
     al_smooth=np.concatenate((tmp_df1,tmp_df2),axis=1).T
     logger_n.info('Too much zeros in datagrame, reorder group2 sample columns with data to the middle %s, %s', pop_col_name,str(insert_col_idx)  )
  else:
     al_smooth=np.concatenate([np.array(rd_smooth1),np.array(rd_smooth2)],axis=0)
  return al_smooth

def precision(label, confusion_matrix):
    col = confusion_matrix[:, label]
    return confusion_matrix[label, label] / col.sum()

def recall(label, confusion_matrix):
    row = confusion_matrix[label, :]
    return confusion_matrix[label, label] / row.sum()

def accuracy(confusion_matrix):
    ''' calculate accurancy of clusteirng based on confusion matrix
    '''
    confusion_matrix=confusion_matrix.astype(float)
    diagonal_sum = confusion_matrix.trace()
    sum_of_all_elements = confusion_matrix.sum()
    return diagonal_sum / sum_of_all_elements

def replace_nan_by_median(selected_data,tmp_id):
   ''' replace nan values by median of data in the dataframe
   '''
   if np.product(selected_data[np.isnan(selected_data)].shape) >0:
     logging.info( 'Replace nan by median in %s ' , tmp_id)
     median_of_data=np.nanmedian(selected_data, axis=1)
     for j in range(selected_data.shape[1]):
          where_is_nan=np.isnan(selected_data[:,j])
          selected_data[where_is_nan,j]=median_of_data[where_is_nan]
     return selected_data
   else:
     return selected_data

def sort_selected_data(selected_pos, selected_data,selected_cols):
   ''' sort data frame by selected_cols
   '''
   new_data=np.concatenate((selected_pos.reshape(selected_pos.shape[0],1), selected_data), axis=1)
   new_columns=['position']+list(selected_cols)
   tmp_df=pd.DataFrame(data=new_data,columns=new_columns)
   sorted_df=tmp_df.sort_values(by='position')
   return sorted_df

def calculate_mean_and_std(record_smooth1):
   ''' compute mean and standard deviation of matrix
   '''
   tumor_smooth=np.mean(np.array(record_smooth1).T,axis=1)
   tumor_smooth_median=np.median(np.array(record_smooth1).T,axis=1)
   #here compute 95% confidence interval 1.96, 98% confidence interval is 2.33, 99% confidence interval is 2.58
   tumor_ci=1.96*np.std(np.array(record_smooth1).T,axis=1)/np.sqrt(len(record_smooth1) )

   return tumor_smooth, tumor_ci, tumor_smooth_median

def do_test_between_two_groups(gcb_smooth, tumor_smooth, smooth_or_interpolated):
   ''' do pair wise comparsion between the two groups mean by using Ranksum, KS-test, T-test
   '''
   p_ranksum=ranksum(gcb_smooth, tumor_smooth, 'Pranksum')
   try :
      ks_test= ks_2samp(gcb_smooth,tumor_smooth)
   except Exception as inst:
      print(inst)
      print(gcb_smooth)
      print(tumor_smooth)
      print(smooth_or_interpolated)
      print('Error in ks_test : do_test_between_two_groups')
      ks_test[0]=0
      ks_test[1]=1
      pass
 
   t_test= ttest_ind(gcb_smooth,tumor_smooth)
   logging.info( 'Data analysis on mean %s ' , smooth_or_interpolated )
   logging.info( 'Ranksum %g ' , p_ranksum)
   logging.info( 'KS-test %g ', ks_test[1])
   logging.info( 'T-test %s ', np.array2string(t_test[1]))
   ks_test_smooth=ks_test[1]
   return p_ranksum, ks_test, t_test, ks_test_smooth

def matrix4eucdist(X, Y):
  '''
    input two matrix and compute euclian-distance between the two matrix elements
    input matrix X, Y row is the sample name, column is the data points or positions
    return squred euclian-distance sum(X-Y)^2
  '''
  #X=gcb_data.T
  #Y=gcb_data.T #tumor_data.T

  U=np.zeros(Y.shape)
  U[~np.isnan(Y)]=1
  Y[np.isnan(Y)]=0

  V=np.zeros(X.shape)
  V[~np.isnan(X)]=1
  X[np.isnan(X)]=0

  d=np.matmul(np.power(X,2),U.T) + np.matmul(V, np.power(Y.T,2)) -2*np.matmul(X,Y.T)
  return d

def distance_matrix2vector(d):
  '''
     input matrix must be symmetric matrix
     get off-diagnoal elements from a matrix, remove duplicated elements and take a squrare root of the remaining elements
     reutrn a vector from an input matrix
  '''
  if d.shape[0]==d.shape[1] :
     #symmetric matrix, within cluster distance
     ##d2=d[np.where(~np.eye(d.shape[0],dtype=bool))]
     #d2=np.sqrt(np.unique(d2))
     ##d2=np.unique(d2)

     #extract upper undiagonal elements from the matrix
     d2= d[np.triu_indices(d.shape[0],k=1)]
  else:
     #unsymmetric matrix, between clusters distance
     #d2=np.sqrt(d)
     d2=d
     d2=d2.reshape(1,np.prod(d2.shape))
  return d2

class jbw_rratios:
   statistic=0.0
   pvalue=0.0
   def __init__(self, dist1,dist2):
       self.dist1=dist1
       self.dist2=dist2
       jbw_rratios.statistic=dist1/dist2
       jbw_rratios.pvalue=norm.sf(abs(dist1/dist2))*2.0

def difference_of_eucdist(X,Y):
   '''
     Input two matrices X and Y and  compute the euclain distance between them and within each group.
     Then test the significant of difference of euclian distances between the within group and the across groups by T-test
     Here assumes X is gcb_data, Y is the tumor_data 
     and matrice X and Y shall have the same column numbers, but the row number can be different
     return T-valie and P-value of the signifance between across group distance and within group distance.
   '''
   #X=gcb_data.T
   #Y=tumor_data.T

   dist2_gcb=matrix4eucdist(X, X)
   dist2_tumor=matrix4eucdist(Y,Y)
   dist2_gcb_tumor=matrix4eucdist(X,Y)
   gcb_dist=distance_matrix2vector(dist2_gcb)
   tumor_dist=distance_matrix2vector(dist2_tumor)
   gcb2tumor_dist=distance_matrix2vector(dist2_gcb_tumor)

   gcb_eu_distance=gcb_dist
   tumor_eu_distance=tumor_dist
   #within_distance=np.concatenate([gcb_eu_distance, tumor_eu_distance], axis=0)
   #within_distance=within_distance.reshape(1,within_distance.shape[0])
   #within_distance=within_distance.reshape(np.prod(within_distance.shape),)

   gcb2tumor_dist=gcb2tumor_dist.reshape(np.prod(gcb2tumor_dist.shape),)

   gcb2tumor_d=np.log2(gcb2tumor_dist+EPS)
   gcb_d=np.log2(gcb_eu_distance+ EPS)
   tumor_d=np.log2(tumor_eu_distance +EPS)
   try:
      if len(gcb2tumor_d)>1 and len(gcb_d)>1 :
         t_test4gcb= ttest_ind(gcb2tumor_d, gcb_d,equal_var=False)
      else:
         t_test4gcb=jbw_rratios(gcb2tumor_d[0],gcb_d[0])
   except:
      print(gcb2tumor_d)
      print(gcb_d)
   if len(gcb2tumor_d)>1 and len(tumor_d)>1:
      t_test4tumor= ttest_ind(gcb2tumor_d, tumor_d,equal_var=False)
   else:
      t_test4tumor=jbw_rratios(gcb2tumor_d[0],tumor_d[0])
   return t_test4gcb, t_test4tumor


def moderated_ttest(group1, group2,is_FDR=False):
   '''
     input two numpy data matrices row is data points, column is samples
     perform moderated t-test between the two groups/matrices
     https://compgenomr.github.io/book/how-to-test-for-differences-between-samples.html#moderated-t-tests-using-information-from-multiple-comparisons
     Smyth GK (2004) "Linear Models and Empirical Bayes Methods for Assessing Differential Expression in Microarray Experiments,
     " Statistical Applications in Genetics and Molecular Biology: Vol. 3: Iss. 1, Article 3. (Pertains to the Moderated t-test.)
   '''
   #data2=np.concatenate((X,Y),axis=0)
   #group1=X.T
   #group2=Y.T
   #count number of samples in each group
   n1=group1.shape[1]
   n2=group2.shape[1]

   dx=np.mean(group1,axis=1)-np.mean(group2,axis=1)
   stderr= np.sqrt( (np.var(group1,axis=1,ddof=1)*(n1-1) + np.var(group2,axis=1,ddof=1)*(n2-1))/(n1+n2-2)*(1/n1+1/n2) )
   mod_stderr= (stderr+ np.median(stderr))/2
   t_mod = dx / mod_stderr
   p_mod = 2*t_pdf.sf( np.abs(t_mod), n1+n2-2 )

   if is_FDR:
      is_rejected, corrected_p= fdrcorrection(p_mod)
   else:
      corrected_p= p_mod

   #startdard T-test
   #t=dx/stderr
   #p = 2*t_pdf.sf( np.abs(t), n1+n2-2 )
   return t_mod, corrected_p

def differential_analaysis_in_two_groups(group1, group2, P_cutoff,is_modT):
   '''input numpy matrix, row is sample name, column is the positions
      return the percentage of data points/positions in the matrics are differentially expressed between the two groups
      is_modT=0,1,2 for T-test, moderate T-test, KS-test respectively
   '''
   #print('P<',P_cutoff)
   #print(group1.shape)
   #print(group2.shape)
   #ttval,tpval=ttest_ind(group1,group2, equal_var=False)
   #print(tpval)
   #is_rejected, corrected_tpval=fdrcorrection(tpval)
   #print(corrected_tpval)

   #mod_tval, mod_tpval=moderated_ttest(group1.T, group2.T)
   #print(mod_tpval)

   if is_modT==1:
     logging.info('Moderated T-test !')
     mod_tval, mod_tpval=moderated_ttest(group1.T, group2.T)
     corrected_tpval=mod_tpval
   elif is_modT==0:
     ttval,tpval=ttest_ind(group1,group2, equal_var=False)
     logging.info('Standard T-test !')
     corrected_tpval=tpval
     #print(group1.shape)
     #print(group2.shape)
     #print(corrected_tpval)
     #print(ttval)
   elif is_modT==2:
       rows,cols=group1.shape
       corrected_tpval=np.zeros((cols,))
       for i in range(cols):
          #print(group1[:,i])
          #print(group2[:,i])
          ks_test= ks_2samp(group1[:,i],group2[:,i])
          #print(ks_test)
          corrected_tpval[i]=ks_test[1]
       #print(corrected_tpval)

   passed_data=np.where(corrected_tpval<P_cutoff)
   percentage_of_passed_data=corrected_tpval[passed_data].shape[0]/corrected_tpval.shape[0]
   #print('Percentage of data points passed FDR of p-values filtering is: ', percentage_of_passed_data)
   return percentage_of_passed_data

def do_PCA_analysis(all_smooth,gcb_col_idx,tumor_col_idx, tmp_column_names) :
   ''' Perform PCA analysis on all_smooth, gcb_col_idx and tumor_col_idx are column index of gcb and tumor group
       Output: PCA dataframe and column names
       ALL input data are numpy array format
   '''
   #do PCA analysis on data 
   #the first 4 is gcb, and the rest is tumor 
   X= StandardScaler().fit_transform(all_smooth)
   col_name= np.concatenate([tmp_column_names[gcb_col_idx[0]],tmp_column_names[tumor_col_idx[0]]],axis=0)
   #feat_cols= pd.DataFrame(X )

   #PCA analysis
   #print(all_smooth,col_name)
   pca_X=PCA(n_components=3)
   principalComponents_X=pca_X.fit_transform(X)
   pca_df=pd.DataFrame(data=principalComponents_X, columns=['PCA 1','PCA 2','PCA 3'])
   logging.info('Explained variation per principal component: {}'.format(pca_X.explained_variance_ratio_))
   return pca_df, col_name, pca_X

def do_tSNE_analysis(all_smooth, gcb_col_idx, tumor_col_idx ):
    ''' Perform t-distributed Stochastic Neighbor Embedding analysis on all_smooth 
    '''
    #tSNE analysis is changed in new version of sklearn >1.2.0 which requires preplexity smaller than samples size
    if  min(len(gcb_col_idx),len(tumor_col_idx)) >4:
        num2perplexity= min(len(gcb_col_idx),len(tumor_col_idx))
    else:
        num2perplexity= len(gcb_col_idx)+len(tumor_col_idx)-1

    tsne = TSNE(n_components=2,learning_rate=1000,perplexity = num2perplexity,early_exaggeration = 12,init = 'random',  random_state=2019)
    X_tsne = tsne.fit_transform(all_smooth)
    tsne_df=pd.DataFrame(data=X_tsne, columns=['PCA_1','PCA_2'])
    #cluster of tSNE then compare to known cluster label
    return tsne_df

def cacl_neg_pos_percent_changes(hist_df,percent_cutoff):
   '''input a dataframe contains percentage column from histrgram, and a percentage_cutoff 
      to calculate percentage of negative or positiave changes greater than the percentage_cutoff
   '''
   #changed after use manual bins in historgram plot
   #negative_percent=hist_df['percentage'][hist_df.start_bin<=-percent_cutoff].sum()
   #positive_percent=hist_df['percentage'][hist_df.end_bin>=percent_cutoff].sum()
   negative_percent=hist_df['percentage'][hist_df.end_bin<=-percent_cutoff].sum()
   positive_percent=hist_df['percentage'][hist_df.start_bin>=percent_cutoff].sum()
   total_percent=negative_percent+positive_percent
   return negative_percent, positive_percent, total_percent


def do_kmeans_on_tSNE(all_tSNE, num_of_clusters,num_of_init,record_smooth2,record_smooth1):
   '''
     input a dataframe results of tSNE, and number of cluster, number of inital values
     gcb matrix smooth2, tumor matrix smooth1
     return a clustering accuracy
   '''
   kmeans = KMeans(init='k-means++', n_clusters=num_of_clusters, n_init=num_of_init)
   #all_tSNE=tSNE_df.to_numpy()
   kmeans.fit(all_tSNE)
   P= kmeans.predict(all_tSNE)

   #evaluate accurancy
   P_true1=np.concatenate([np.zeros(len(record_smooth2)), np.ones(len(record_smooth1))])
   P_true2=np.concatenate([np.ones(len(record_smooth2)), np.zeros(len(record_smooth1))])

   cm1=confusion_matrix(P_true1, list(P))
   cm2=confusion_matrix(P_true2, list(P))
   cluster_accuracy=max(accuracy(cm1), accuracy(cm2))
   return cluster_accuracy

def export_selected_dmr_data(tmp_chr,tmp_id,block_methylation_columns2, gcb_col_idx2,
                tumor_col_idx2, position2data, all_smooth2,dotOrUnderscore, out_folder, out_file_string='.dat' ):
    out_col_name= np.concatenate([ block_methylation_columns2[gcb_col_idx2[0]], block_methylation_columns2[tumor_col_idx2[0]]],axis=0)
    new_col_name=[]
    for i in out_col_name:
        if dotOrUnderscore==1:
          #changed jbw march 2022
          tmp_i=i.split('_')
          new_col_name.append('_'.join([tmp_i[i] for i in (1,3,-1) ]))
        else:
          #changed jbw april 2022, here is the bug which needs to find a better way for exporting new column names??
          tmp_i=i.split('.')[0]
          if tmp_chr in tmp_i:
            new_col_name.append(tmp_i)
          else:
            new_col_name.append('_'.join([tmp_i,tmp_chr]))

    #print(all_smooth*100, out_col_name,new_col_name)
    out_data=all_smooth2*100
    out_dmr_df=pd.DataFrame(data=out_data.T, columns=new_col_name)
    out_dmr_df['position']=position2data
    out_dmr_df=out_dmr_df.set_index('position')
    #out_dmr_file='_'.join([tmp_i[-1],tmp_id+out_file_string])
    out_dmr_file='_'.join([tmp_chr,tmp_id+out_file_string])

    #print(tmp_id,out_dmr_df)
    #output MRdata path already set up in do_block_finding at dmr_analysis_block, 
    #here is a double check  
    out_folder_data=os.path.join(out_folder, tmp_chr,'data')
    if not os.path.exists(out_folder_data):
       os.mkdir(out_folder_data)
       print('Create: ' , out_folder_data)
    out_dmr_file=os.path.join(out_folder_data,out_dmr_file)
    out_dmr_df.to_csv(out_dmr_file,index=True,sep='\t',float_format='%g',compression='gzip')


def do_DMR_analysis_in_a_block(logger_n, percent_cutoff, selected_data2tumor,
                selected_pos, selected_data2gcb,tmp_len, tmp_id, tmp_num_of_data,
                gcb_col_idx, tumor_col_idx, block_methylation_columns,fig,isSmooth,isExport,dotOrUnderscore,P_cutoff, 
                is_modT, low_median_high_cutoff, wildType_fileString,out_folder,tmp_chr):
    '''
        isSmooth=2 then use smoothed data for DMR analysis and return P-values for smoothed data
        otherwise, if isSmooth==1 use interpolated data, 0 for using raw data to calculate difference between groups
    '''
    if isSmooth==1:
       dataStr='interpolate'
    elif isSmooth==2:
       dataStr='smooth'
    elif isSmooth==0:
       dataStr='raw'
    num_of_gcb=len(gcb_col_idx[0])
    num_of_tumor=len(tumor_col_idx[0])
    #plot raw tumo data
    sub_plot_num=231
    xnew, record_smooth1, record_interpolated1, num_of_data_size= plot_raw_tumor_data(fig,sub_plot_num,selected_pos,selected_data2tumor,tmp_len,tmp_id,tmp_num_of_data,dataStr)

    #plot raw gcb data
    sub_plot_num=232
    xnew2, record_smooth2, record_interpolated2, num_of_data_size= plot_raw_gcb_data(fig,sub_plot_num,selected_pos,selected_data2gcb,tmp_len,tmp_id,tmp_num_of_data,dataStr)

    #jbw debuge
    #if tmp_id=='mr7623':
    #    print('Position :', selected_pos)
    #    print('ID : ', tmp_id)
    #    print('data : ', tmp_num_of_data)

    #compute mean and standard deviation of smoothed data
    tumor_smooth, tumor_ci, tumor_smooth_median=[],[],[]
    gcb_smooth, gcb_ci, gcb_smooth_median=[],[],[]

    #return mean, index, median values for smoothed data to plot confidence interval
    tumor_smooth, tumor_ci,tumor_smooth_median = calculate_mean_and_std(record_smooth1)
    gcb_smooth, gcb_ci, gcb_smooth_median = calculate_mean_and_std(record_smooth2)

    #logging.info('Percentage of data points passed T-test filtering on '+ dataStr+' data')
    if isSmooth==1:
       #logging.info('Percentage of data points passed T-test filtering on interpolated data')
       gcb_data=np.array(record_interpolated2)
       tumor_data=np.array(record_interpolated1)
    elif isSmooth==2:
       #logging.info('Percentage of data points passed T-test filtering on smoothed data')
       gcb_data=np.array(record_smooth2)
       tumor_data=np.array(record_smooth1)
    elif isSmooth==0:
       #logging.info('Percentage of data points passed T-test filtering on raw data')
       gcb_data=np.array(selected_data2gcb.T.tolist())
       tumor_data=np.array(selected_data2tumor.T.tolist())

    #differential data analysis between samples on smoothed data
    if is_modT==0:
       testStr='T-test'
    elif is_modT==1:
       testStr='modT-test'
    elif is_modT==2:
       testStr='KS-test'

    #set a mininum percentage in data before log2 transforming percentage for using T-test to see how many positions are passed P_cutoff
    gcb_dd=np.log2((gcb_data+EPS)*100)
    tumor_dd=np.log2((tumor_data+EPS)*100)
    #print('gcb\n', gcb_dd.shape)
    #print(gcb_dd)
    #print('tumor\n', tumor_dd.shape)
    #print(tumor_dd)
    percentage_of_data_passed_filtering=differential_analaysis_in_two_groups(gcb_dd, tumor_dd, P_cutoff,is_modT)
    logger_n.info('Percentage of '+dataStr+' data points passed '+ testStr + ' P<'+ str(P_cutoff) + ' is: %g', percentage_of_data_passed_filtering)

    #calculate euclidean distance between two groups
    gcb_euDist_ttest, tumor_euDist_ttest=difference_of_eucdist(gcb_data,tumor_data)

    tumor_raw, tumor_rw_ci,tumor_raw_median=[],[],[]
    gcb_raw,gcb_rw_ci,gcb_raw_median=[],[],[]
    tumor_raw,tumor_rw_ci,tumor_raw_median=calculate_mean_and_std(selected_data2tumor.T.tolist())
    gcb_raw,gcb_rw_ci,gcb_raw_median=calculate_mean_and_std(selected_data2gcb.T.tolist())

    tumor_interpolated, tumor_int_ci, tumor_interpolated_median=[],[],[]
    tumor_interpolated, tumor_int_ci, tumor_interpolated_median =  calculate_mean_and_std(record_interpolated1)
    gcb_interpolated, gcb_int_ci, gcb_interpolated_median=[],[],[]
    gcb_interpolated, gcb_int_ci, gcb_interpolated_median=  calculate_mean_and_std(record_interpolated2)

    #do ranksum test between mean on smoothed data
    p_ranksum, ks_test, t_test, ks_test_smooth= do_test_between_two_groups(gcb_smooth, tumor_smooth,'smoothed_data')

    #do ranksum test between mean on interpolated data
    #jbw debug
    #if tmp_id=='mr2009':
    #   print(gcb_raw)
    #   print(tumor_raw)
    #   print(gcb_interpolated)
    #   print(tumor_interpolated)
    p_ranksum_int, ks_test_int, t_test_int, ks_test_interpolated= do_test_between_two_groups(gcb_interpolated, tumor_interpolated,'interpolated_data')

    #do PCA analysis on smoothed or interpolated data and plot the first 3 PCA
    ##first 4 is gcb, and the rest is tumor 
    #here force all PCA and clustering on smoothed Data!
    logger_n.info( 'Perform PCA , t-SNE analysis on '+dataStr + ' data !')
    if isSmooth==1:
        all_smooth=np.concatenate([np.array(record_interpolated2),np.array(record_interpolated1)], axis=0)
        #all_smooth=check_zeros_in_df(record_interpolated2, record_interpolated1,logger_n)
        #logging.info( 'Perform PCA , t-SNE analysis on interpolated data !')
    elif isSmooth==2 :
        #jbw 2021 for avoiding underflow in pca
        all_smooth=np.concatenate([np.array(record_smooth2),np.array(record_smooth1)],axis=0)
        #all_smooth=check_zeros_in_df(record_smooth2, record_smooth1,logger_n)
        #logging.info('Perform PCA, tSNE, and clustering analysis on smoothed data !')
    elif isSmooth==0:
        all_smooth=np.concatenate([np.array(selected_data2gcb.T.tolist()),np.array(selected_data2tumor.T.tolist())],axis=0)
        #all_smooth=check_zeros_in_df(selected_data2gcb.T.tolist(), selected_data2tumor.T.tolist(), logger_n)
        #logging.info('Perform PCA, tSNE analysis on raw data !')
    all_raw_data=np.concatenate([np.array(selected_data2gcb.T.tolist()),np.array(selected_data2tumor.T.tolist())],axis=0)

    #do PCA analysis Debug
    #print(record_smooth2,record_smooth1)
    #logger_n.info(all_smooth)
    #logger_n.info(xnew)
    #if sample order is reordered in order to avoid underflowing in PCA, then sample labels in the same group may be wrong!!
    #this can be a bug for plot figures. Thus we skip the PCA calculation if there is too much columns with zeros.
    if not too_much_zeros(all_smooth, logger_n):
      pca_df, all_smooth_columns, pca_X= do_PCA_analysis(all_smooth,gcb_col_idx,tumor_col_idx, block_methylation_columns )

      #plot 3D PCA
      sub_plot_num=233
      plot_3D_PCA(fig,sub_plot_num,pca_df, wildType_fileString,num_of_gcb,num_of_tumor)
    else:
      logger_n.info('Skip PCA calculation because of too much columns with zeros in dataframe!')

    #do tSNE analysis
    tSNE_df= do_tSNE_analysis(all_smooth,gcb_col_idx[0], tumor_col_idx[0] )

    #do clustering on tSNE data
    num_of_clusters=2
    num_of_init=10
    all_tSNE=tSNE_df.to_numpy()
    cluster_accuracy=do_kmeans_on_tSNE(all_tSNE, num_of_clusters,num_of_init,record_smooth2,record_smooth1)
    logger_n.info( dataStr +' clusteirng accuracy of tSNE %g', cluster_accuracy)

    #plot 2D TSNE
    sub_plot_num=235
    plot_2D_tSNE(fig,sub_plot_num,tSNE_df, cluster_accuracy,dataStr,num_of_gcb,num_of_tumor)

    #plot smoothed data for the two groups
    sub_plot_num=234
    plot_smoothed_data(fig,sub_plot_num,selected_pos,xnew,tumor_smooth,tumor_ci,xnew2,gcb_smooth,gcb_ci,percentage_of_data_passed_filtering,dataStr,testStr,P_cutoff)

    #export data for DMR if it is needed
    if isExport:
      export_selected_dmr_data(tmp_chr,tmp_id,block_methylation_columns, gcb_col_idx,tumor_col_idx,list(xnew), all_smooth,dotOrUnderscore, 
                               out_folder, out_file_string='.dat.gz' )
      export_selected_dmr_data(tmp_chr,tmp_id,block_methylation_columns, gcb_col_idx,tumor_col_idx,selected_pos.tolist(), all_raw_data, dotOrUnderscore,
                            out_folder, out_file_string='_raw.dat.gz' )

    #compute the mean/median difference between the two group samples
    isMean=True
    if isMean:
       delt_of_raw=tumor_raw - gcb_raw
       delt_of_smooth = tumor_smooth - gcb_smooth
       delt_of_interpolated = tumor_interpolated -gcb_interpolated
    else:
       delt_of_raw=tumor_raw_median - gcb_raw_median
       delt_of_smooth = tumor_smooth_median - gcb_smooth_median
       delt_of_interpolated = tumor_interpolated_median -gcb_interpolated_median


    #at least 5% methylation level has changes in a data point betwween tumor and normal groups, Now need do U<=0.05, 0.05<L<=0.1, 0.1<M<=0.2, H>0.2 DMR  ?? 
    #print(percent_cutoff)
    logger_n.info('Minimum precentage of the mean methylation changes between the two groups (e.g., mean(tumor) -mean(gcb)) in each position : Low=%g, Median=%g, High=%g', percent_cutoff[0],percent_cutoff[1],percent_cutoff[2])

    #plot historm of delt_of_smooth
    #hist_df2, results2 are interpolated data
    #only plot the High percent_cutoff in percent_cutoff[L,M,H]
    #there is a darray warning bug in this function
    #jbw debug
    #if tmp_id=='mr7623':
    #   print(delt_of_raw)
    #   print(delt_of_smooth)
      
    hist_df, hist_df2, hist_df3,  results, results2 , results3 = plot_hist_of_difference_in_groups(fig, delt_of_raw,delt_of_smooth, delt_of_interpolated,percent_cutoff, low_median_high_cutoff, gcb_euDist_ttest.pvalue, tumor_euDist_ttest.pvalue)

    negative_percent,positive_percent,total_percent=np.zeros(3),np.zeros(3),np.zeros(3)
    negative_percent2,positive_percent2,total_percent2=np.zeros(3),np.zeros(3),np.zeros(3)
    negative_percent3,positive_percent3,total_percent3=np.zeros(3),np.zeros(3),np.zeros(3)
    if isSmooth==2:
       #delt of median smooth
       for loop in range(3):
              negative_percent[loop],positive_percent[loop],total_percent[loop]=  cacl_neg_pos_percent_changes(hist_df,percent_cutoff[loop])

       logger_n.info( 'Difference of the mean smooth data , Negative percentage= %g,  %g,  %g ' , negative_percent[0],negative_percent[1],negative_percent[2])
       logger_n.info( ' Low, Median, High                   Positive percentae= %g,  %g,  %g ', positive_percent[0], positive_percent[1], positive_percent[2])
       logger_n.info( '                                     Total perecntage= %g,  %g,  %g', total_percent[0],total_percent[1], total_percent[2] )
    elif isSmooth==1:
       #delt of median interpolated
       for loop in range(3):
           negative_percent2[loop],positive_percent2[loop],total_percent2[loop]= cacl_neg_pos_percent_changes(hist_df2,percent_cutoff[loop])

       logger_n.info( 'Difference of the mean interpolate data, Negative2 percentage= %g,  %g,  %g ' , negative_percent2[0],negative_percent2[1],negative_percent2[2])
       logger_n.info( 'Low, Median, High                        Positive2 percentage= %g,  %g,  %g ', positive_percent2[0],  positive_percent2[1], positive_percent2[2])
       logger_n.info( '                                         Total percengate= %g,  %g,  %g', total_percent2[0],total_percent2[1], total_percent2[2])
    elif isSmooth==0:
       #delt of median raw
       for loop in range(3):
          negative_percent3[loop], positive_percent3[loop],total_percent3[loop] = cacl_neg_pos_percent_changes(hist_df3,percent_cutoff[loop])

       logger_n.info( 'Difference of the mean raw data,     Negative3 percentage= %g,  %g,  %g ' , negative_percent3[0],negative_percent3[1],negative_percent3[2])
       logger_n.info( 'Low, Mediam, High                    Positive3 percentage= %g,  %g,  %g ', positive_percent3[0], positive_percent3[1], positive_percent3[2])
       logger_n.info( '                                     Total percengate= %g,  %g,  %g', total_percent3[0],total_percent3[1], total_percent3[2] )

    #record data
    #use interpolated to record
    #export both smoothed and interpolated P-value
    statistic_test_p=[ks_test_smooth,ks_test_interpolated, percentage_of_data_passed_filtering, \
                        gcb_euDist_ttest.pvalue, tumor_euDist_ttest.pvalue, gcb_euDist_ttest.statistic, tumor_euDist_ttest.statistic]

    #export percentage based on either smoothed or interpolated data points
    if isSmooth==2:
        #return smooth
        return statistic_test_p, cluster_accuracy, negative_percent, positive_percent,  results, total_percent
    elif isSmooth==1:
         #return interpolated
        return statistic_test_p, cluster_accuracy, negative_percent2, positive_percent2, results2, total_percent2
    else:
        #return raw
        return statistic_test_p, cluster_accuracy, negative_percent3, positive_percent3, results3, total_percent3
 
