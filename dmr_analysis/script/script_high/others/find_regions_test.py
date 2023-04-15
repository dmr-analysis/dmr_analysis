#python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d

def find_methylation_blocks(tmp_position,tmp_methylation, tmp_column_names, max_read_length,mini_block_size):
  '''
    Input: tmp_position is a numpy array of chromosome positions,
           tmp_methylation is a numpy data matrix 
           tmp_column_names is a list of column names/sample name in the data matrix
           max_read_length is the maximum length between two adjacent data points that can be included in the same block/region
           mini_block_size is the minimum number of data points that are needed for a valid block 
    Output: blcoks is predicted block
            block_size is the number of data points in the block 
            block_length is the length of block
            block_methylation is the  data matrix within the block
            All results in dictionary format
    Find all methylation blocks passed filterring conditions 
  '''
  #distance between two consecutive points
  delt= abs(tmp_position[0:-1]- tmp_position[1:])
  #find index of points > max_length
  is_break_index= np.where(delt>=max_read_length)[0]
  #print(delt)
  #print(np.where(delt>=max_read_length))
  #print(np.max(delt))
  #print(is_break_index)

  #find idex for blocks
  #plus 2 is for python does not read the last position in slice, and one position move back in delt position calculation.
  if is_break_index.size>0:
    block_position=np.zeros((len(is_break_index)+2,2))
    loop=0
    for i in is_break_index :
      if loop==0:
        block_position[loop,0]=0
        block_position[loop,1]=i+1
      else:
        block_position[loop,0]=block_position[loop-1,1]
        block_position[loop,1]=i+1
      loop += 1

    #add the last point
    if is_break_index[-1]<len(tmp_position):
       block_position[loop,0]=block_position[loop-1,1]
       block_position[loop,1]=len(tmp_position)
       loop +=1
  else:
    #here assue no block is found and all input data belong to the same block
    block_position=np.zeros((1,2))
    block_position[0,0] = 0
    block_position[0,1] = len(tmp_position)
   
  #remove rows with all zeros
  block_position=block_position[~np.all(block_position==0, axis=1)]
  block_position=block_position.astype(int)
 
  #extract position from each block
  blocks={}
  blocks_size={}
  blocks_methylation={}
  blocks_length={}
  minimum_size=mini_block_size
  loop = 0
  for ii in range(block_position.shape[0]):
    tmp_len= block_position[ii,1]-block_position[ii,0]
    if tmp_len >= minimum_size:
       tmp_id='mr'+str(loop)
      #blocks.append(tmp_position[ block_position[ii,0]: block_position[ii,1] ])
      #blocks_size.append( tmp_len)
      #blocks_methylation.append(tmp_methylation[ block_position[ii,0]: block_position[ii,1] ,:])
       blocks[tmp_id]=tmp_position[ block_position[ii,0]: block_position[ii,1] ]
       blocks_length[tmp_id]=max(blocks[tmp_id])-min(blocks[tmp_id])
       blocks_size[tmp_id]=tmp_len
       if len(tmp_methylation.shape)>1:
         blocks_methylation[tmp_id] = tmp_methylation[ block_position[ii,0]: block_position[ii,1] ,:]
       else:
         blocks_methylation[tmp_id] = tmp_methylation[ block_position[ii,0]: block_position[ii,1] ,]
       loop += 1

  return blocks, blocks_size, blocks_methylation, blocks_length



if __name__=='__main__':
 ##########################################
 #tmp_position=np.load('tmp_position.npy')
 #tmp_methylation=np.load('tmp_methylation.npy')
 #tmp_column_names=np.load('tmp_column_names.npy',allow_pickle=True)

 test_index=None
 #tmp_position=np.loadtxt('tmp_position.csv',delimiter=',')[:test_index]
 #tmp_methylation=np.loadtxt('tmp_methylation.csv',delimiter=',')[:test_index]
 #tmp_column_names=np.loadtxt('tmp_column_names.csv',delimiter=',',dtype=object)

 tmp_pd=pd.read_csv('tmp_out_pd.csv',sep='\t')
 tmp_pd2=tmp_pd[0:test_index]
 tmp_position=tmp_pd2.Starts.values
 tmp_methylation=tmp_pd2.iloc[:,3:15].values
 tmp_column_names=tmp_pd2.columns[3:].values
 max_read_length=500
 mini_block_size= 5
 print('Select ', test_index, ' rows from input')
 print('Blocks with distance greater than ', max_read_length, ' and minimum ', mini_block_size, ' data points in a block')
 block, block_size, block_methylation,block_length = find_methylation_blocks(tmp_position,tmp_methylation, tmp_column_names, max_read_length,mini_block_size)
 block_methylation_columns=tmp_column_names

 #find gcb in columns
 tmp_cols=block_methylation_columns.astype(str)
 gcb_col_idx=np.where(np.char.find(tmp_cols,'gcb')>=0)
 tumor_col_idx=np.where(np.char.find(tmp_cols,'gcb')<0)

 #for each MR perfrom data interpretation
 #for kk in block.keys():
 if True:
   kk=block.keys()[0]
   tmp_id=kk
   tmp_data=block_methylation[kk]
   tmp_pos=block[kk]
   tmp_len=block_size[kk]
   plt.plot(tmp_pos,tmp_data[:,gcb_col_idx[0]],'ro')
   x=tmp_pos
   y=tmp_data[:,gcb_col_idx[0]]
   flank_space=50
   xnew=np.linspace(min(tmp_pos),max(tmp_pos),50)
   xnew=xnew.astype(int)
   for i in range(y.shape[1]):
     f1 = interp1d(x, y[:,i],kind = 'linear')
     plt.plot(xnew,f1(xnew),'g-')
#   plt.plot(tmp_pos,tmp_data[:,tumor_col_idx[0]],'g+-')
   plt.title(str(tmp_len) +' '+ tmp_id)
   plt.show()
   raw_input(" ")






