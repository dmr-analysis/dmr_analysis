#plot function in DMR pipeline
#plot historgram of block length, and block data points size
import matplotlib as mlt
import pandas as pd
import numpy as np
mlt.use('Agg')
import matplotlib.pyplot as plt
#load library for plot 3D
from mpl_toolkits.mplot3d import Axes3D
from .dmr_utility import make_new_bin_range 
import logging
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d

#np.seterr(all='raise')
np.seterr(all='ignore')

def interpolate_and_smooth_selected_data(selected_pos2, selected_data2):
   ''' data interpolatation of selected pos and selected data
      and data smooth for interpolated data points
      output: new_x, smoothed data points, interpolated data points, num_of_data_size in interpolated data points
   '''
   x=selected_pos2
   y=selected_data2

   #num_of_data_size=100
   delt_length=max(selected_pos2)-min(selected_pos2)

   #based on lenght of MR to estimate the step size of data interpolation
   if delt_length> 10000:
        step_size=250
   elif delt_length>5000 and delt_length<=10000:
        step_size=50
   elif delt_length>100 and delt_length <= 5000:
        step_size=10
   elif delt_length>20  and delt_length<=100:
        step_size=2
   else:
        step_size=1
   num_of_data_size=int(delt_length/step_size)
   #print num_of_data_size
   xnew=np.linspace(min(selected_pos2),max(selected_pos2),num_of_data_size)
   xnew=xnew.astype(int)

   #plot interpoloated data points
   #record tumor data
   record_smooth1=[]
   record_interpolated1=[]
   #debuge
   #print(x)
   #print(y)
   for i in range(y.shape[1]):
       str2kind='linear'
       #str2kind='nearest'
       f1 = interp1d(x, y[:,i],kind = str2kind, assume_sorted=False)
       interpolated_y=f1(xnew)
       #plt.plot(xnew,interpolated_y,'g-')
       #try to smooth the curve
       smooth=gaussian_filter1d(interpolated_y, sigma=3)
       record_smooth1.append(list(smooth))
       record_interpolated1.append(list(interpolated_y))

   return xnew, record_smooth1, record_interpolated1, num_of_data_size

def plot_hist_for_block_length_and_data_size(block_length,block_size,max_block_length,out_fig_name,logger_num):
   ''' Plot historgram of block_lenght and block_size (number of data points in a block) distributions
   '''
   mr_len=[]
   mr_size=[]
   for kk in block_length.keys():
     mr_len.append(block_length[kk])
     mr_size.append(block_size[kk])

   logger_num.info('minimum MR length %g ', min(mr_len))
   logger_num.info('maximum MR length %g ', max(mr_len))
   plt.subplot(1,2,1)
   if max_block_length<= 2000 and  max_block_length > 1000  :
     #2000 
     bin_range=[0, 500, 1000, 5000, 10000, 50000 ,100000, 500000,max(mr_len) +100 ]
   elif max_block_length <= 1000 and max_block_length > 500 :
     #1000 length
     if max(mr_len)> 100000:
        bin_range=[min(mr_len), 100,500,1000,5000,10000,50000, 100000,  max(mr_len) + 100]
     else:
        bin_range=[min(mr_len), 100,500,1000,5000,10000,50000,  max(mr_len) + 100]
   elif max_block_length<= 500 and max_block_length > 200 :
     #500 length
     bin_range=[min(mr_len), 100,500,1000,5000,10000,50000,  max(mr_len) +100]
   elif max_block_length<= 200 :
     #200 length
     bin_range=[min(mr_len), 100,500,1000,5000,10000,  max(mr_len) + 100]
   else:
     bin_range=[0, 500, 1000, 5000, 10000, 50000 ,100000, 500000, max(mr_len) +100 ]

   #python3 reshape the bin_range if maximum value is smaller than the other range in the list
   new_bin_range=make_new_bin_range(bin_range)
   print(new_bin_range)

   n, bins, patches= plt.hist(x=mr_len,bins=new_bin_range)
   logger_num.info( 'Maximum length of adjancey CpG sites in a block %s ', str(max_block_length))
   logger_num.info( 'Hist plot n %s ' , np.array2string(n.astype(int)))
   logger_num.info( '         bins %s ', np.array2string(bins.astype(int)))
   plt.title('MR Length')

   plt.subplot(1,2,2)
   #2000
   #ln,lbins,lpatches=plt.hist(x=mr_size,bins=[0,50,100,500,1000,5000,10000,50000])

   #1000 length
   #bin_range=[min(mr_size), 100,500,1000,5000,max(mr_size)]

   #500/200 length
   bin_range=[min(mr_size), 100,500,1000,max(mr_size)]
   new_bin_range=make_new_bin_range(bin_range)
   ln,lbins,lpatches=plt.hist(x=mr_size,bins= new_bin_range ,range=(min(mr_size),max(mr_size)))
   logger_num.info('mininum MR data size %g', min(mr_size))
   logger_num.info('maximum MR data size %g', max(mr_size))
   #print ln.astype(int) , lbins.astype(int)
   plt.title('MR data size')
   #plt.show()
   plt.savefig(out_fig_name,dpi=100)


def plot_3D_PCA(fig,sub_plot_num, pca_df,wildType_fileString,num_of_gcb,num_of_tumor ):
   ''' 3-D plot of the first 3 PCAs
   '''
   #first is GCB or wildtype group , then comes Tumor group
   #plot 3D PCA
   ax=fig.add_subplot(sub_plot_num,projection='3d')
   my_color=[1]*num_of_gcb + [2]*num_of_tumor
   ax1= ax.scatter(pca_df['PCA 1'][0:num_of_gcb], pca_df['PCA 2'][0:num_of_gcb], pca_df['PCA 3'][0:num_of_gcb], c='green', s=80,label=wildType_fileString.replace('_','') )
   ax2= ax.scatter(pca_df['PCA 1'][num_of_gcb:num_of_tumor+num_of_gcb], pca_df['PCA 2'][num_of_gcb:num_of_tumor+num_of_gcb], pca_df['PCA 3'][num_of_gcb:num_of_tumor+num_of_gcb], c='red', s=80,label='tumor/KO')

   # make simple, bare axis lines through space:
   xAxisLine = ((min(pca_df['PCA 1']), max(pca_df['PCA 1'])), (0, 0), (0,0))
   ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'orange')
   yAxisLine = ((0, 0), (min(pca_df['PCA 2']), max(pca_df['PCA 2'])), (0,0))
   ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'orange')
   zAxisLine = ((0, 0), (0,0), (min(pca_df['PCA 3']), max(pca_df['PCA 3'])))
   ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'orange')

   #add label for plotted data points
   if False:
     for idx, row in pca_df.iterrows():
          ax.text(row['PCA 1'], row['PCA 2'], row['PCA 3'], idx)

   # label the axes
   ax.set_xlabel("PC1")
   ax.set_ylabel("PC2")
   ax.set_zlabel("PC3")
   ax.set_title("PCA ")
   ax.legend()


   #ax=fig.add_subplot(235,projection='3d')
   #ax1=ax.scatter(pca_df['PCA 1'][0:4], pca_df['PCA 2'][0:4],  pca_df['PCA 3'][0:4], c='green', s=80, label='gcb')
   #ax2=ax.scatter(pca_df['PCA 1'][4:12], pca_df['PCA 2'][4:12],  pca_df['PCA 3'][4:12], c='red', s=80,  label='tumor')
   #xAxisLine = ((min(pca_df['PCA 1']), max(pca_df['PCA 1'])), (0, 0), (0,0))
   #ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r-')
   #yAxisLine = ((0, 0), (min(pca_df['PCA 2']), max(pca_df['PCA 2'])), (0,0))
   #ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r-')
   #yAxisLine = ((0, 0), (min(pca_df['PCA 3']), max(pca_df['PCA 3'])), (0,0))
   #ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r-')
   return ax


def plot_2D_PCA(fig,sub_plot_num, pca_df,wildType_fileString, num_of_gcb,num_of_tumor,tumorType_fileString='tumor/KO',color2wildType='green' , color2tumorType='red'):
   ''' 2-D plot of the first 3 PCAs
   	first is GCB or wildtype group , then comes Tumor group
   	plot 2D PCA
   	color2wildType default is green
   	color2tumorTyp default is red
   '''
   ax=fig.add_subplot(sub_plot_num)
   my_color=[1]*num_of_gcb + [2]*num_of_tumor
   ax1= ax.scatter(pca_df['PCA 1'][0:num_of_gcb], pca_df['PCA 2'][0:num_of_gcb], c=color2wildType, s=80,label=wildType_fileString )
   ax2= ax.scatter(pca_df['PCA 1'][num_of_gcb:num_of_tumor+num_of_gcb], pca_df['PCA 2'][num_of_gcb:num_of_tumor+num_of_gcb], c=color2tumorType, s=80,label=tumorType_fileString)

   # make simple, bare axis lines through space:
   #xAxisLine = ((min(pca_df['PCA 1']), max(pca_df['PCA 1'])), (0, 0), (0,0))
   #ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'blue')
   #yAxisLine = ((0, 0), (min(pca_df['PCA 2']), max(pca_df['PCA 2'])), (0,0))
   #ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'blue')

   #add label for plotted data points
   if False:
     for idx, row in pca_df.iterrows():
          ax.text(row['PCA 1'], row['PCA 2'], idx)
   # label the axes
   ax.set_xlabel("PC1",fontsize=12, fontweight='bold')
   ax.set_ylabel("PC2",fontsize=12, fontweight='bold')
   #ax.set_title("PCA ", fontsize=14, fontweight='bold')
   ax.legend()
   return ax


def plot_2D_tSNE(fig,sub_plot_num,pca_df,cluster_accuracy,dataStr,num_of_gcb,num_of_tumor) :
    '''a 2-D plot tSNE
    '''
    #plot 2D
    ax=fig.add_subplot(sub_plot_num)
    ax.scatter(pca_df['PCA_1'][0:num_of_gcb], pca_df['PCA_2'][0:num_of_gcb], c='green', s=80, label='gcb')
    ax.scatter(pca_df['PCA_1'][num_of_gcb:num_of_gcb+num_of_tumor], pca_df['PCA_2'][num_of_gcb:num_of_gcb+num_of_tumor],  c='red', s=80,  label='tumor')
    xAxisLine = ((min(pca_df['PCA_1']), max(pca_df['PCA_1'])), (0, 0), (0,0))
    ax.plot(xAxisLine[0], xAxisLine[1],  'orange')
    yAxisLine = ((0, 0), (min(pca_df['PCA_2']), max(pca_df['PCA_2'])), (0,0))
    ax.plot(yAxisLine[0], yAxisLine[1],  'orange')

    # label the axes
    ax.set_xlabel("tsne1")
    ax.set_ylabel("tsne2")

    #if isSmooth==1:
    #   dataStr='interpolate'
    #elif isSmooth==2:
    #   dataStr='smooth'
    #elif isSmooth==0:
    #   dataStr='raw'
    ax.set_title("t-SNE (cluster " + dataStr + " "+ str('{:{width}.{prec}f}'.format(cluster_accuracy,width=5,prec=3)) + ')' )
    #ax.legend()


def plot_hist_of_difference_in_groups(fig, delt_of_raw, delt_of_smooth, delt_of_interpolated, all_percent_cutoff, low_median_high_level, gcb_euDist_ttest_pvalue, tumor_euDist_ttest_pvalue ):
   ''' Plot a historgram of the differential methylation between the two samples
       hist_df for smooth, hist_df2 for interpolated, hist_df3 for raw data
   '''
   percent_cutoff=all_percent_cutoff[low_median_high_level]
   ax=fig.add_subplot(236)

   #here to made bin range based on low,median, high level of delt(tumor-gcb) changes
   xbins1=np.sort(np.array(all_percent_cutoff))
   xbins2=xbins1.copy()
   xbins=np.concatenate([xbins1,xbins2*(-1)])

   xx=np.array([np.array(delt_of_smooth,dtype="object"), np.array(delt_of_interpolated,dtype="object"), np.array(delt_of_raw,dtype="object")],dtype="object").T
   abs_xx=np.abs(xx)
   tmp_max=[np.amax(abs_xx[0]), np.amax(abs_xx[1]) ,np.amax(abs_xx[2]) ]
   max_bin=np.amax(tmp_max) +0.2
   if max_bin >1 :
      max_bin=1
   step_size= np.median(np.abs(xbins1[0:-1]-xbins1[1:]))
   #extend bin to positive side
   max_bin_range1=np.arange(xbins1[-1], max_bin, step_size)
   #jbw
   if len(max_bin_range1)>0:
     max_bin_range1=np.delete(max_bin_range1,0)

   #extend bin to negative side
   max_bin_range2=max_bin_range1.copy()
   max_bin_range2=max_bin_range2*(-1)
   #merge all bins
   xbins=xbins.tolist()+ max_bin_range1.tolist() + max_bin_range2.tolist() +[0]
   xbins.sort()
   #make minimum and maximum as -1 and 1 for the range
   #jbw december 2021
   if xbins[0]>-1:
      xbins[0]=-1
   if xbins[-1]<1:
      xbins[-1]=1

   #this line has a warning bug for narray but use dtype=object shall pay attention to the shape of array 
   #the label is the column names of the numpy array in x
   #therefore has to conver the np array to row represents data points, columns represent three types of data
   dn,dbins,dpatches=ax.hist(x=xx,bins=xbins,label=['smooth','interpolate','raw'])

   #smooth results
   #jbw debug 
   #if np.isnan(np.sum(dn[0])):
   #  print(xx)
   #  print(dn)

   results= np.array([dn[0], dn[0]/np.sum(dn[0]),dbins[0:-1],dbins[1:]]).T
   #interpolated results
   results2= np.array([dn[1], dn[1]/np.sum(dn[1]),dbins[0:-1],dbins[1:]]).T
   #raw results
   results3=np.array([dn[2], dn[2]/np.sum(dn[2]),dbins[0:-1],dbins[1:]]).T

   #only plot the High cut off of changes 
   ax.axvline(x=percent_cutoff,color='y')
   ax.axvline(x=-percent_cutoff,color='y')

   ax.set_title('(tum-gcb)' + str(' [gcb vs gps {:{width}.{prec}e}, tumor vs gps {:{width}.{prec}e}]'.format(gcb_euDist_ttest_pvalue, tumor_euDist_ttest_pvalue,width=3,prec=2)) )
   ax.legend()

   hist_df=pd.DataFrame(data=results, columns=['counts','percentage','start_bin','end_bin'])
   hist_df2=pd.DataFrame(data=results2, columns=['counts','percentage','start_bin','end_bin'])
   hist_df3=pd.DataFrame(data=results3, columns=['counts','percentage','start_bin','end_bin'])
   #smooth, interpolated
   #result are -.... om smoothed data
   #results2 are histgram bin and counts of differental methylation in interpolated data (count, percentage, bin_pos_start, bin_pos_end) 

   #compute total number of percentage > percent_cutoff
   #results2 (count, percentage changes, start bin, end bin)
   #sum_of_negative_count=results[results[:,2]<=-percent_cutoff,0].sum()
   #sum_of_positive_count=results[results[:,3]>=percent_cutoff,0].sum()
   #total_counts=results[:,0].sum()
   #total_percentage_gt_cutoff=(sum_of_negative_count+sum_of_positive_count)/total_counts

   #sum_of_negative_count=results2[results2[:,2]<=-percent_cutoff,0].sum()
   #sum_of_positive_count=results2[results2[:,3]>=percent_cutoff,0].sum()
   #total_counts=results2[:,0].sum()
   #total_percentage_gt_cutoff2=(sum_of_negative_count+sum_of_positive_count)/total_counts

   # *2 for interpolated data results
   # *3 for raw data results

   return hist_df, hist_df2, hist_df3, results, results2,  results3

def plot_raw_tumor_data(fig,sub_plot_num,selected_pos,selected_data2tumor,tmp_len,tmp_id,tmp_num_of_data,dataStr):
    fig.add_subplot(sub_plot_num)
    #plot raw tumo data
    plt.plot(selected_pos,selected_data2tumor ,'o',ms=1.5,color='orange')

    #plot tumor interpoloated1 and smoothed data1
    xnew, record_smooth1, record_interpolated1 = [],[],[]
    xnew, record_smooth1,record_interpolated1,num_of_data_size=interpolate_and_smooth_selected_data(selected_pos, selected_data2tumor)
    for i in range(len(record_smooth1)):
         plt.plot(xnew,record_smooth1[i],'-',color='red')

    plt.title( ' Length='+ str(tmp_len) +' '+ tmp_id + ' #=' + str(tmp_num_of_data) + ' #' + dataStr+'= '+ str(num_of_data_size))
    plt.axis([min(selected_pos)-50, max(selected_pos)+50, 0 ,1 ])
    return xnew, record_smooth1, record_interpolated1, num_of_data_size
 

def plot_raw_gcb_data(fig,sub_plot_num,selected_pos,selected_data2gcb,tmp_len,tmp_id,tmp_num_of_data,dataStr):
    fig.add_subplot(sub_plot_num)
    #plot raw gcb data
    #jbw
    #print(selected_pos)
    #print(selected_data2gcb)
    if selected_data2gcb.size==0:
       print('No data of selected MR is found, please check input option for "wild type file string", I STOP!')
       exit(1)
    else:
       plt.plot(selected_pos,selected_data2gcb ,'o',ms=1.5,color='y')

    #plot gcb interpolated2 and smoothed data2
    record_smooth2, record_interpolated2, xnew2=[],[],[]
    xnew2, record_smooth2,record_interpolated2,num_of_data_size=interpolate_and_smooth_selected_data(selected_pos, selected_data2gcb)

    for i in range(len(record_smooth2)):
         plt.plot(xnew2,record_smooth2[i],'-', color='green')
    plt.title( ' Length='+ str(tmp_len) +' '+ tmp_id + ' #=' + str(tmp_num_of_data) + ' #'+dataStr+'= '+ str(num_of_data_size))
    plt.axis([min(selected_pos)-50, max(selected_pos)+50, 0 ,1 ])
    return xnew2, record_smooth2, record_interpolated2, num_of_data_size
    

def plot_smoothed_data(fig,sub_plot_num,selected_pos,xnew,tumor_smooth,tumor_ci,xnew2,gcb_smooth,gcb_ci,percentage_of_data_passed_filtering,dataStr,testStr,P_cutoff):
    #plot smoothed data and their confidence interval
    fig.add_subplot(sub_plot_num)
    plt.plot(xnew2,gcb_smooth,'g-',label='gcb')
    plt.plot(xnew,tumor_smooth,'r-',label='tumor')
    plt.fill_between(xnew2, (gcb_smooth - gcb_ci), (gcb_smooth + gcb_ci), color='g', alpha=0.2)
    plt.fill_between(xnew, (tumor_smooth -tumor_ci), (tumor_smooth + tumor_ci), color='r', alpha=0.2)
    plt.axis([min(selected_pos)-50, max(selected_pos)+50, 0 ,1 ])
    #plt.title('KS-test P value ' + str(  '{:{width}.{prec}e}'.format(ks_test_interpolated, width=5, prec=3) ) )
    plt.title(str('{:{width}.{prec}g}'.format(percentage_of_data_passed_filtering*100,width=5,prec=3)) + '% '+ dataStr + ' '+ testStr +' P<' + str('{:{width}.{prec}}'.format(P_cutoff,width=5,prec=3)))







