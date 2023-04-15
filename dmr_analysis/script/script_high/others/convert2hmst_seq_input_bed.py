#this script is used to convert common bed methylation file to HMST-seq-abalyze supported BED format 
#especially to export files to each chrom folder separately.
#for python3
import gzip
import os
import pandas as pd
from multiprocessing import Process
import glob

def runInParallel(fn,chrs,output_folder, output_name,out_df):
  proc=[]
  for tchr in chrs:
     p=Process(target=fn,args=(tchr,output_folder,output_name,out_df))
     p.start()
     proc.append(p)
  for p in proc:
     p.join()

def read_bed_gz_files(f1):
  ''' Here the default is 6 columns in a file such as chr,start, end, methy, count, strand ,
    here read a bed files and return a dataframe of the input file
  '''
  record_data=[]
  record_chrs=[]
  with gzip.open(f1,'rt') as inf1:
    try:
      for line in inf1:
          line=line.rstrip('\n')
          lines=line.split('\t')
          if len(lines)>6:
             #compute total counts
             tmp_total=int(float(lines[4])+float(lines[6]))
          else:
             #assume total counts exist
             tmp_total=int(float(lines[4]))
          if lines[0]=='23':
             lines[0]='X'
          elif lines[0]=='24':
             lines[0]='Y'
          if 'chr' not in lines[0]:
              tmp_data=['chr'+lines[0], int(lines[1]),int(lines[2]),lines[3],str(tmp_total),lines[5] ]
          else:
              tmp_data=[lines[0], int(lines[1]),int(lines[2]),lines[3],str(tmp_total),lines[5] ]
          record_data.append(tmp_data)
          if lines[0] not in record_chrs:
             record_chrs.append(lines[0])
    except gzip.BadGzipFile:
      print(f1)
      print('input_file is not a valid gzip file by BadGzipFile')
  df=pd.DataFrame(record_data,columns=['Chrs','Starts','Ends','Methylation','Total_counts','Strand'])
  #map object to float or int or string
  #python2
  #df['Chrs']=map(str,df['Chrs'])
  #df['Starts']=map(int,df['Starts'])
  #df['Ends']=map(int,df['Ends'])
  #df['Methylation']=map(float,df['Methylation'])
  #df['Total_counts']=map(int,df['Total_counts'])
  #df['Strand']=map(str,df['Strand'])

  #python3
  df['Chrs']=df['Chrs'].astype(str)
  df['Starts']=df['Starts'].astype(int)
  df['Ends']=df['Ends'].astype(int)
  df['Methylation']=df['Methylation'].astype(float)
  df['Total_counts']=df['Total_counts'].astype(int)
  df['Strand']=df['Strand'].astype(str)
 
  print(len(record_data))
  print(record_chrs)
  return df

def call_export_files(i,output_folder,output_name,out_df):
   if i<23:
      tmp_chr=str(i)
   elif i==23:
      tmp_chr='X'
   elif i==24:
      tmp_chr='Y'
   out_folder=os.path.join(output_folder,'chr'+tmp_chr)

   if not os.path.isdir(out_folder):
       os.system('mkdir '+ out_folder)
       print("Create output file folder " + out_folder)
   else:
       print("Output file folder " + out_folder + " exist")
   out_file=output_name+'_chr'+tmp_chr+'.bed'
   out_file2=os.path.join(out_folder,out_file)
   #print "Export data at " , out_file2
   #print tmp_chr
   tmp_out_data=out_df.copy()
   tmp_out_data=tmp_out_data[tmp_out_data['Chrs'] == 'chr'+tmp_chr]
   #if tmp_chr=='1':
   #   print len(tmp_out_data)
   if len(tmp_out_data)>0:
      #sort data
      print("Export to ", out_file2)
      sorted_tmp_out_data=tmp_out_data.sort_values(['Starts'],ascending=True)
      sorted_tmp_out_data.to_csv(out_file2,sep='\t',index=False,header=False)
      os.system('gzip -f ' + out_file2)
   else:
      print("Not exist " , out_file2)

if __name__ == "__main__":
  #pool=mp.Pool(mp.cpu_count())

  input_folder='extracted_data/'
  #input_name='gcb_meth1_4160735_4160735'
  #input_name='gcb_meth1_4118819_4118819'
#  input_name='gcb_meth1_4122131_4122131'
  input_names=['gcb_meth1_4174884_4174884','meth_4159170_4159170','meth_4158726_4158726','meth_4121361_4121361',\
                'meth_4189200_4189200','meth_4134005_4134005','meth_4188900_4188900','meth_4175837_4175837',\
                 'meth_4105105_4105105']
  for input_name in input_names:   
    input_file_name=input_name+'_WGBS_GRCh37.bed.gz'
    input_file=os.path.join(input_folder,input_file_name)
    print("Read ", input_file)

    output_folder='in_data/WGBS-data/'
    output_name=input_name

    out_df=read_bed_gz_files(input_file)
    all_chrs=[]

    #maximum threads 12
    for i in range(1,7):
       all_chrs.append(i)  
    runInParallel(call_export_files,all_chrs,output_folder, output_name,out_df)
 
    all_chrs=[]
    for i in range(7,13):
       all_chrs.append(i)
    runInParallel(call_export_files,all_chrs,output_folder, output_name,out_df)
 
    all_chrs=[]
    for i in range(13,19):
       all_chrs.append(i)
    runInParallel(call_export_files,all_chrs,output_folder, output_name,out_df)

    all_chrs=[]
    for i in range(19,25):
       all_chrs.append(i)
    runInParallel(call_export_files,all_chrs,output_folder, output_name,out_df)




