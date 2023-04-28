import argparse
import os
import pandas as pd
import sys
# from .scripts.preprocess_reference import main as pr_main
# from .scripts.create_bed import main as cb_main

def make_names(chr, start, end, rna_name, x, y, reg, gene_name, strand,
               orig_start, orig_end):
    """Creates the string names and returns it. The format of names is:
    chr:start_pos:end_pos:rna_name||region_name:X:Y||gene_name:strand:gene_start:gene_end

    Keyword arguments:
    row -- row of dataframe passed to make_TSS_TES_gene_5dist()
    start -- start position
    end -- end position
    x -- X
    y -- Y
    name -- region name"""

    names = ('%s:%s:%s:%s||%s:%d:%d||%s:%c:%s:%s' % (
             chr,                     #chr
             start,                   #start_pos depending on type(TSS, TES, gene, 5dist)
             end,                     #end_pos of type(TSS, TES, gene, 5dist)
             rna_name,                #rna_name
             reg,                     #TSS, TES, gene or 5dist
             x,                       #X
             y,                       #Y
             gene_name,               #gene_name
             strand,                  #strand
             orig_start,              #start_pos
             orig_end))               #end_pos
    return names



def combine_non_unique(df):
    """Combines similar genomic regions starting and ending at the same site, but
    having different rna_names"""
    pd.options.display.max_colwidth = 100
    #print(df)
    groups = df.groupby(['chr', 'start_pos', 'end_pos', 'x', 'y', 'gene_name', 'reg'], sort=False)
    all_rows = []
    for g in groups:
        group_df = g[1].reset_index()

        chr = group_df.chr[0]
        start_pos = group_df.start_pos[0]
        end_pos = group_df.end_pos[0]
        X = group_df.x[0]
        Y = group_df.y[0]
        gene_name = group_df.gene_name[0]
        reg = group_df.reg[0]
        strand = group_df.strand[0]

        if len(group_df) > 1:
            rna_name = '&'.join(group_df.rna_name.tolist())
            orig_start = group_df.orig_start.drop_duplicates().astype(str).tolist()
            orig_end = group_df.orig_end.drop_duplicates().astype(str).tolist()
            if len(orig_start) > 1:
                orig_start = '&'.join(orig_start)
            else:
                orig_start = orig_start[0]
            if len(orig_end) > 1:
                orig_end = '&'.join(orig_end)
            else:
                orig_end = orig_end[0]
        else:
            rna_name = group_df.rna_name[0]
            orig_start = group_df.orig_start[0].astype(str)
            orig_end = group_df.orig_end[0].astype(str)

        names = make_names(chr, start_pos, end_pos, rna_name, X, Y, reg, gene_name,
                           strand, orig_start, orig_end)
        new_row = [chr, start_pos, end_pos, names]
        all_rows.append(new_row)

    new_df = pd.DataFrame(all_rows, columns=['chr', 'start_pos', 'end_pos', 'names'])
    return new_df

def make_TSS_TES_gene_5dist(df, X, Y, M, N, rem):
    """Creates four BED files: TSS, TES, gene and 5dist from the reference

    df -- sorted dataframe, columns=['chr', 'start_pos', 'end_pos', 'gene_name:rna_name',
                                     'dot', 'strand']
    X -- number of upstream bp, TSS, TES, gene
    Y -- number of downstream bp, TSS, TES, gene
    M -- number of bp from gene start site, 5dist
    N -- number of bp from gene start site, 5dist
    rem -- True if regions should be removed, False if TSS, TES and 5dist should be kept
    """
    removed = []
    # added jbw
    TSS = [];
    TES = [];
    gene = [];
    dist5 = [];
    dist5D = []
    size = len(df.values)

    for i, row in enumerate(df.values):
        gene_start_pos = int(row[1]) + Y
        gene_end_pos = int(row[2]) - Y
        rna_gene = row[3].split(':')
        rna_name = rna_gene[1]
        gene_name = rna_gene[0]
        if gene_start_pos < gene_end_pos:
            gene.append({'chr': row[0],
                         'start_pos': gene_start_pos,
                         'end_pos': gene_end_pos,
                         'x': X,
                         'y': Y,
                         'gene_name': gene_name,
                         'rna_name': rna_name,
                         'orig_start': row[1],
                         'orig_end': row[2],
                         'strand': row[5],
                         'reg': 'gene'})
            # 'names': make_names(row, gene_start_pos, gene_end_pos,
            #                    X, Y, 'gene')})

            TSS_start_pos = row[1] - X if row[5] == '+' else row[2] - Y
            TSS_end_pos = row[1] + Y if row[5] == '+' else row[2] + X
            TSS.append({'chr': row[0],
                        'start_pos': TSS_start_pos,
                        'end_pos': TSS_end_pos,
                        'x': X,
                        'y': Y,
                        'gene_name': gene_name,
                        'rna_name': rna_name,
                        'orig_start': row[1],
                        'orig_end': row[2],
                        'strand': row[5],
                        'reg': 'TSS'})
            # 'names': make_names(row, TSS_start_pos, TSS_end_pos, X, Y,
            #                    'TSS')})

            # added wang
            if TSS_start_pos < 1:
                TSS_start_pos = 1

            TES_start_pos = row[1] - X if row[5] == '-' else row[2] - Y
            TES_end_pos = row[1] + Y if row[5] == '-' else row[2] + X
            TES.append({'chr': row[0],
                        'start_pos': TES_start_pos,
                        'end_pos': TES_end_pos,
                        'x': X,
                        'y': Y,
                        'gene_name': gene_name,
                        'rna_name': rna_name,
                        'orig_start': row[1],
                        'orig_end': row[2],
                        'strand': row[5],
                        'reg': 'TES'})
            # 'names': make_names(row, TES_start_pos, TES_end_pos, X, Y,
            #                    'TES')})
            # added wang
            if TES_start_pos < 1:
                TES_start_pos = 1

            # 5'distance upstream
            dist5_start_pos = row[1] - N if row[5] == '+' else row[2] + M
            dist5_end_pos = row[1] - M if row[5] == '+' else row[2] + N
            if dist5_start_pos < 1:
                dist5_start_pos = 1

            # added wang
            if dist5_end_pos < 1:
                dist5_end_pos = 1

            if not (dist5_end_pos == 1 & dist5_start_pos == 1):
                dist5.append({'chr': row[0],
                              'start_pos': dist5_start_pos,
                              'end_pos': dist5_end_pos,
                              'x': M,
                              'y': N,
                              'gene_name': gene_name,
                              'rna_name': rna_name,
                              'orig_start': row[1],
                              'orig_end': row[2],
                              'strand': row[5],
                              'reg': '5dist'})
                # 'names': make_names(row, dist5_start_pos, dist5_end_pos, M,
                #                    N, '5dist')})

            # added jbw for 5distD Downstream
            dist5D_start_pos = row[1] + M if row[5] == '+' else row[2] - N
            dist5D_end_pos = row[1] + N if row[5] == '+' else row[2] - M
            if dist5D_start_pos < 1:
                dist5D_start_pos = 1

            if dist5D_end_pos < 1:
                dist5D_end_pos = 1

            if not (dist5D_end_pos == 1 & dist5D_start_pos == 1):
                dist5D.append({'chr': row[0],
                               'start_pos': dist5D_start_pos,
                               'end_pos': dist5D_end_pos,
                               'x': M,
                               'y': N,
                               'gene_name': gene_name,
                               'rna_name': rna_name,
                               'orig_start': row[1],
                               'orig_end': row[2],
                               'strand': row[5],
                               'reg': '5distD'})
        else:
            if rem:  # remove TSS, TES, 5dist regions even though geneBody too small
                removed.append({'chr': row[0],
                                'start_pos': row[1],
                                'end_pos': row[2],
                                'x': 0,
                                'y': 0,
                                'gene_name': gene_name,
                                'rna_name': rna_name,
                                'orig_start': row[1],
                                'orig_end': row[2],
                                'strand': row[5],
                                'reg': 'removed_from_geneBody_TSS_TES_5dist'})
                # 'names': make_names(row, row[1], row[2], 0, 0,
                #                    'removed_from_geneBody_TSS_TES_5dist')})
            else:
                removed.append({'chr': row[0],
                                'start_pos': row[1],
                                'end_pos': row[2],
                                'x': 0,
                                'y': 0,
                                'gene_name': gene_name,
                                'rna_name': rna_name,
                                'orig_start': row[1],
                                'orig_end': row[2],
                                'strand': row[5],
                                'reg': 'removed_only_from_geneBody'})
                # 'names': make_names(row, row[1], row[2], 0, 0,
                #                    'removed_only_from_geneBody')})

                TSS_start_pos = row[1] - X if row[5] == '+' else row[2] - Y
                TSS_end_pos = row[1] + Y if row[5] == '+' else row[2] + X
                TSS.append({'chr': row[0],
                            'start_pos': TSS_start_pos,
                            'end_pos': TSS_end_pos,
                            'x': X,
                            'y': Y,
                            'gene_name': gene_name,
                            'rna_name': rna_name,
                            'orig_start': row[1],
                            'orig_end': row[2],
                            'strand': row[5],
                            'reg': 'TSS'})
                # 'names': make_names(row, TSS_start_pos, TSS_end_pos, X, Y,
                #                    'TSS')})

                # added wang
                if TSS_start_pos < 1:
                    TSS_start_pos = 1

                TES_start_pos = row[1] - X if row[5] == '-' else row[2] - Y
                TES_end_pos = row[1] + Y if row[5] == '-' else row[2] + X
                TES.append({'chr': row[0],
                            'start_pos': TES_start_pos,
                            'end_pos': TES_end_pos,
                            'x': X,
                            'y': Y,
                            'gene_name': gene_name,
                            'rna_name': rna_name,
                            'orig_start': row[1],
                            'orig_end': row[2],
                            'strand': row[5],
                            'reg': 'TES'})
                # 'names': make_names(row, TES_start_pos, TES_end_pos, X, Y,
                #                    'TES')})

                # added wang
                if TES_start_pos < 1:
                    TES_start_pos = 1

                # 5'distance upstream
                dist5_start_pos = row[1] - N if row[5] == '+' else row[2] + M
                dist5_end_pos = row[1] - M if row[5] == '+' else row[2] + N
                if dist5_start_pos < 1:
                    dist5_start_pos = 1
                # added wang
                if dist5_end_pos < 1:
                    dist5_end_pos = 1

                if not (dist5_start_pos == 1 & dist5_end_pos == 1):
                    dist5.append({'chr': row[0],
                                  'start_pos': dist5_start_pos,
                                  'end_pos': dist5_end_pos,
                                  'x': M,
                                  'y': N,
                                  'row': row,
                                  'gene_name': gene_name,
                                  'rna_name': rna_name,
                                  'orig_start': row[1],
                                  'orig_end': row[2],
                                  'strand': row[5],
                                  'reg': '5dist'})
                    # 'names': make_names(row, dist5_start_pos, dist5_end_pos, M,
                    #                    N, '5dist')})

                # added jbw 5distance downstream
                dist5D_start_pos = row[1] + M if row[5] == '+' else row[2] - N
                dist5D_end_pos = row[1] + N if row[5] == '+' else row[2] - M
                if dist5D_start_pos < 1:
                    dist5D_start_pos = 1

                if dist5D_end_pos < 1:
                    dist5D_end_pos = 1

                if not (dist5D_start_pos == 1 & dist5D_end_pos == 1):
                    dist5D.append({'chr': row[0],
                                   'start_pos': dist5D_start_pos,
                                   'end_pos': dist5D_end_pos,
                                   'x': M,
                                   'y': N,
                                   'row': row,
                                   'gene_name': gene_name,
                                   'rna_name': rna_name,
                                   'orig_start': row[1],
                                   'orig_end': row[2],
                                   'strand': row[5],
                                   'reg': '5distD'})

    TSS_df = pd.DataFrame(TSS)
    TSS_df = combine_non_unique(TSS_df)

    TES_df = pd.DataFrame(TES)
    TES_df = combine_non_unique(TES_df)

    gene_df = pd.DataFrame(gene)
    gene_df = combine_non_unique(gene_df)

    dist5_df = pd.DataFrame(dist5)
    dist5_df = combine_non_unique(dist5_df)

    removed_df = pd.DataFrame(removed)
    if not removed_df.empty:
        removed_df = combine_non_unique(removed_df)

    # added jbw
    dist5D_df = pd.DataFrame(dist5D)
    dist5D_df = combine_non_unique(dist5D_df)

    return TSS_df, TES_df, gene_df, dist5_df, removed_df, dist5D_df


def make_intergenic_region(genome_file, reference, bed_files, min_intergenic_reg_len, max_intergenic_reg_len, X, Y,
                           out_folder, intergenic_between_genes):
    """Finds intergenic regions with regions bigger than min_intergenic_reg_len,
    and returns dataframe.

    genome_file -- genome file name
    bed_files -- BED filenames: [TSS_filename, TES_filename, gene_filename]
    min_intergenic_reg_len -- minimum size of intergenic region
    out_folder -- folder name where output should go
    intergenic_between_genes -- complement on gene only or TSS, TES and gene? find intergenic
                  regions between gene regions, or include TSS and TES as well?"""

    if intergenic_between_genes:
        # Runs the command line which creates intergenic regions BED file
        os.system('bedtools complement -i '+reference+' -g '
                  +genome_file+' > '+out_folder+'/'+'intergenic_regions.bed')
    else:
        os.system('cat '+bed_files[0]+' '+ bed_files[1]+' '+bed_files[2]+\
                  ' > '+out_folder+'/'+'combined.bed')

        #added jbw to check whether negative values in bed file and replace them as 1
        tmp_file=out_folder+'/'+'combined.bed'
        tmp_df=pd.read_csv(tmp_file,sep='\t',header=None)
        tmp_df.loc[tmp_df[1].apply(lambda x: x<0),1]=1
        tmp_df.to_csv(tmp_file,sep='\t',index=False, header=None)
        #end change jbw

        os.system('bedtools sort -g '+genome_file+' -i '+out_folder+'/'+\
                  'combined.bed > '+out_folder+'/'+'combined_sorted.bed')

        os.system('bedtools merge -i '+out_folder+'/'+'combined_sorted.bed > \
                  '+out_folder+'/'+'combined_sorted_merged.bed')

        os.system('bedtools complement -i '+out_folder+'/'+'combined_sorted_merged.bed -g ' +
                  genome_file + ' > '+out_folder+'/'+'intergenic_regions.bed')
        #jbw test
        os.system('rm '+out_folder+'/'+'combined*.bed')

    inter_gen = pd.read_csv(out_folder+'/'+'intergenic_regions.bed', sep='\t',
                names=['chr', 'start_pos', 'end_pos'])
    #jbw test
    os.system('rm '+out_folder+'/'+'intergenic_regions.bed')

    # Removes intergenetic regions shorter than parameter
    inter_gen = inter_gen[((inter_gen['end_pos'] - inter_gen['start_pos']) >
                min_intergenic_reg_len) & ((inter_gen['end_pos'] - inter_gen['start_pos']) <
                max_intergenic_reg_len)]

    inter_gen['names'] = inter_gen.chr.astype(str) + ':' +\
                         inter_gen.start_pos.astype(str) + ':' +\
                         inter_gen.end_pos.astype(str) + '||' +\
                         'intergenic:' +\
                         ('None' if intergenic_between_genes else str(X)) + ':' +\
                         ('None' if intergenic_between_genes else str(Y))
    return inter_gen


def cb_main(reference, genomeFile, X, Y, M, N, intergenic_between_genes, min_intergenic_len, max_intergenic_len,
         remove_short, out_folder):
    if reference[-4:].lower() != '.bed':
        print(
            '***** ERROR in create_bed.py: The uniquely sorted reference must be a BED formatted file. Exiting. *****')
        exit(1)

    if intergenic_between_genes in ('yes', 'y'):
        intergenic_between_genes = True
    elif intergenic_between_genes in ('no', 'n'):
        intergenic_between_genes = False
    else:
        print('***** ERROR in create_bed.py: intergenic_between_genes must be either [yes/y] or [no/n]. Exiting. *****')
        exit(1)

    if remove_short in ('yes', 'y'):
        remove_short = True
    elif remove_short in ('no', 'n'):
        remove_short = False
    else:
        print('***** ERROR in create_bed.py: removeShort must be either [yes/y] or [no/n]. Exiting. *****')
        exit(1)

    ref_df = pd.read_csv(reference, index_col=False, header=0, sep='\t')
    reference_name = reference.split('/')[-1][:-4]

    dfs = make_TSS_TES_gene_5dist(ref_df, X, Y, M, N, remove_short)
    # added jbw
    TSS_df, TES_df, gene_df, dist5_df, removed_df, dist5D_df = dfs

    if remove_short:
        TSS_filename = out_folder + '/TSS_Up' + str(X) + '_Down' + str(Y) + '_removedShort' + '.bed'
        TES_filename = out_folder + '/TES_Up' + str(X) + '_Down' + str(Y) + 'removedShort' + '.bed'
        gene_filename = out_folder + '/gene_Up' + str(X) + '_Down' + str(Y) + 'removedShort' + '.bed'
        dist5_filename = out_folder + '/5dist_Up' + str(N) + '_Up' + str(M) + 'removedShort' + '.bed'
        removed_filename = out_folder + '/removed_regions_all_TSS_TES_5dist_geneBodyLessThan0.bed'
        # added jbw
        dist5D_filename = out_folder + '/5dist_Down' + str(N) + '_Down' + str(M) + 'removedShort' + '.bed'
    else:
        TSS_filename = out_folder + '/TSS_Up' + str(X) + '_Down' + str(Y) + '.bed'
        TES_filename = out_folder + '/TES_Up' + str(X) + '_Down' + str(Y) + '.bed'
        gene_filename = out_folder + '/gene_Up' + str(X) + '_Down' + str(Y) + 'removedShort' + '.bed'
        dist5_filename = out_folder + '/5dist_Up' + str(N) + '_Up' + str(M) + '.bed'
        removed_filename = out_folder + '/removed_regions_geneBodyLessThan0.bed'
        # added jbw
        dist5D_filename = out_folder + '/5dist_Down' + str(N) + '_Down' + str(M) + '.bed'

    order = ['chr', 'start_pos', 'end_pos', 'names']
    TSS_df = TSS_df[order]
    TES_df = TES_df[order]
    gene_df = gene_df[order]
    dist5_df = dist5_df[order]
    if not removed_df.empty:
        removed_df = removed_df[order]
    # added jbw
    dist5D_df = dist5D_df[order]

    # added wang check negative postion in exported files
    TSS_df.loc[TSS_df['start_pos'].apply(lambda x: x < 0), 'start_pos'] = 1
    TES_df.loc[TES_df['start_pos'].apply(lambda x: x < 0), 'start_pos'] = 1

    TSS_df.to_csv(TSS_filename, sep='\t', index=False, header=None)
    TES_df.to_csv(TES_filename, sep='\t', index=False, header=None)
    gene_df.to_csv(gene_filename, sep='\t', index=False, header=None)
    dist5_df.to_csv(dist5_filename, sep='\t', index=False, header=None)
    if not removed_df.empty:
        removed_df.to_csv(removed_filename, sep='\t', index=False, header=None)
    # added jbw
    dist5D_df.to_csv(dist5D_filename, sep='\t', index=False, header=None)

    # Creates filtered intergenic region BED file
    bed_files = [TSS_filename, TES_filename, gene_filename]
    inter_gen = make_intergenic_region(genomeFile, reference, bed_files, min_intergenic_len, max_intergenic_len,
                                       X, Y, out_folder, intergenic_between_genes)
    inter_gen = inter_gen[order]

    if intergenic_between_genes:
        intergen_filename = out_folder + '/' + 'intergenic_uniqueSorted_betweenGenes_minLen' + str(
            min_intergenic_len) + '.bed'
    else:
        intergen_filename = out_folder + '/' + 'intergenic_uniqueSorted_betweenTSS_TES_genes_minLen' + str(
            min_intergenic_len) + '.bed'
    inter_gen.to_csv(intergen_filename, sep='\t', index=False, header=None)

    # added jbw
    return [TSS_filename, TES_filename, gene_filename, dist5_filename, intergen_filename, removed_filename,
            dist5D_filename]

def make_BED_of_ref(df):
    """Creates BED formatted dataframe from input file and returns the dataframe

    Keyword arguments:
    df -- dataframe of reference"""

    df['gene_rna_name'] = df['gene_name'] + ':' + df['rna_name']
    df['dot'] = '.'
    df_BED = pd.DataFrame(df[['chr', 'start_pos', 'end_pos', 'gene_rna_name',
                              'dot', 'strand']])
    return df_BED
def get_num_chr(human, string):
    """Returns the numbers of the respective chromosome names

    Keyword arguments:
    human -- is the data from human sample or mouse or rat ? Yes/mouse/rat
    string -- X, Y or M"""
    if human =='human':
        if string == 'X':
            return '23'
        elif string == 'Y':
            return '24'
        elif string == 'M':
            return '25'
    elif human =='mouse ':   # mouse
        if string == 'X':
            return '20'
        elif string == 'Y':
            return '21'
        elif string == 'M':
            return '22'
    elif human =='rat': #rat
        if string == 'X':
            return '21'
        elif string == 'Y':
            return '22'
        elif string == 'M':
            return '23'

def find_uniq_genes(df, num_chr_name, human, remove_mir):
    """Filters and orders data and writes to file

    Keyword arguments:
    df -- pandas dataframe
    num_chr_name -- numerical chromosome name?
    human -- data from human? if False -> mouse or rat
    remove_mir -- remove genes with names starting with mir? (default True)
    """
    df_copy = df.copy()

    if remove_mir:
        #removes rows starting with case insensitive 'mir':
        df_copy = df_copy[df_copy.gene_name.str.lower().str.startswith('mir') == False]

    #removes rows where chr contains '_':
    df_copy = df_copy[df_copy.chr.str.contains('_') == False]

    #removes duplicates
    df_copy = df_copy.drop_duplicates(subset=['chr', 'strand', 'TSS', 'TES'])


    chrs = df_copy.chr.str[3:]

    X_mask = chrs == 'X'
    chrs.loc[X_mask] = get_num_chr(human, 'X')

    Y_mask = chrs == 'Y'
    chrs.loc[Y_mask] = get_num_chr(human, 'Y')

    M_mask = chrs == 'M'
    chrs.loc[M_mask] = get_num_chr(human, 'M')

    df_copy['order_by_chr'] = chrs

    if num_chr_name:
        df_copy.chr = pd.to_numeric(df_copy.order_by_chr)

        #sorts first by chr, then TSS, then rna_name
        df_copy = df_copy.sort_values(['chr', 'TSS', 'rna_name'], ascending=True)

    else:
        df_copy.order_by_chr = pd.to_numeric(df_copy.order_by_chr)
        df_copy = df_copy.sort_values(['order_by_chr', 'TSS', 'rna_name'], ascending=True)

    return df_copy


def pr_main(reference, human, numerical_chr_name, remove_mir, folder_out):
    """Cleans reference file (refFlat)"""

    if human.lower() in ('yes', 'y'):
        human = 'human'
    elif human.lower() in ('mouse', 'm'):
        human = 'mouse'
    elif human.lower() in ('rat', 'r'):
        human = 'rat'
    else:
        print('***** ERROR in preprocess_reference.py: No info whether it is human or mouse. If human sample, enter [yes/y]. If mouse sample, enter [no/n]. Exiting. *****')
        exit(1)

    if numerical_chr_name.lower() in ('yes', 'y'):
        numerical_chr_name = True
    elif numerical_chr_name.lower() in ('no', 'n'):
        numerical_chr_name = False
    else:
        print('***** ERROR in preprocess_reference.py: numerical chromosome name must be whether [yes/y] or [no/n]. Exiting. *****')
        exit(1)

    if remove_mir.lower() in ('yes', 'y'):
        remove_mir = True
    elif remove_mir.lower() in ('no', 'n'):
        remove_mir = False
    else:
        print(remove_mir)
        print('***** ERROR: removeMir must be either [yes/y] or [no/n]. Exiting. *****')
        exit(1)
    outfilename = reference.split('/')[-1][:-4]


    df = pd.read_csv(reference, sep='\t', usecols=[0, 1, 2, 3, 4, 5],
                names=['gene_name', 'rna_name', 'chr', 'strand', 'TSS', 'TES'])

    new_ref = find_uniq_genes(df, numerical_chr_name, human, remove_mir)

    new_ref.rename(columns={'TSS': 'start_pos', 'TES': 'end_pos'}, inplace=True)
    ref_BED_df = make_BED_of_ref(new_ref)

    # Makes sure the columns are in this exact order
    new_ref = new_ref[['rna_name', 'chr', 'strand', 'start_pos', 'end_pos', 'gene_name']]
    ref_BED_df = ref_BED_df[['chr', 'start_pos', 'end_pos', 'gene_rna_name', 'dot', 'strand']]

    # Keeps header
    new_ref.to_csv(folder_out+'/'+outfilename+'_clean_sorted.txt', sep='\t', index=False)

    out_bed = folder_out+'/'+outfilename+'_clean_sorted.bed'
    # No header
    ref_BED_df.to_csv(out_bed, sep='\t', index=False, header=None)

    return out_bed



def main(reference_file, genome_file, human, numerical_chr_name, remove_mir, X, Y,
         M, N, intergenic_between_genes, min_intergenic_len, max_intergenic_len, remove_short, folder_out):
    """Calls preprocess_reference and create_bed, which cleans the reference and
    extracts genomic regions from the genes of the reference"""

    out_list_file = folder_out+'/list_region_files.txt'
    folder_out += '/data'

    ref_bed = pr_main(reference_file, human, numerical_chr_name, remove_mir, folder_out)

    #added wang
    region_files = cb_main(ref_bed, genome_file, X, Y, M, N, intergenic_between_genes,
                           min_intergenic_len, max_intergenic_len, remove_short, folder_out)


    outfiles = open(out_list_file, 'w')
    #added jbw
    TSS, TES, gene, dist5, intergenic, removed , dist5D = region_files
    outfiles.write(TSS+'\n')
    outfiles.write(gene+'\n')
    outfiles.write(TES+'\n')
    outfiles.write(dist5+'\n')
    outfiles.write(intergenic+'\n')
    #added jbw
    outfiles.write(dist5D +'\n')
    outfiles.close()

    return out_list_file

def my_parser(parser):
    required = parser.add_argument_group("Required")
    required.add_argument("-r", "--referenceFile",
                        metavar='file',
                        required=True,
                        help="reference file")
    required.add_argument("-g", "--genomeFile",
                        metavar='file',
                        required=True,
                        help="the genome file, one column with chromosome names and one column with the length of each chromosome")
    required.add_argument("-hu", "--human",
                        metavar='y/m/r',
                        required=True,
                        help="is sample from human? [yes/y] or [mouse/m] or [rat/r]. if no, it must be from either mouse/m or rat/r genome")
    required.add_argument("-n", "--numericalChr",
                        metavar='y/n',
                        required=True,
                        help="will you operate with numerical chromosome names or not, [yes/y] or [no/n]. the same as chosen here must be as in the genomeFile, default=no")

    optional = parser.add_argument_group("Optional, has default values")
    optional.add_argument("-rm", "--removeMir",
                        metavar='',
                        default='yes',
                        help="remove genes with names starting with 'Mir'? [yes/y] or [no/n]. default=yes")
    optional.add_argument("-X",
                        metavar='',
                        type=int,
                        default=1000,
                        help="the number of upstream bp, TSS, TES, gene. default=1000")
    optional.add_argument("-Y",
                        metavar='',
                        type=int,
                        default=1000,
                        help="the number of downstream bp, TSS, TES, gene. default=1000.")
    optional.add_argument("-M",
                        metavar='',
                        type=int,
                        default=10000,
                        help="the number of bp from gene start site, 5dist. default=10000.")
    optional.add_argument("-N",
                        metavar='',
                        type=int,
                        default=1000000,
                        help="the number of bp from gene start site, 5dist. default=1000000.")
    optional.add_argument("-i", "--intergenicBTGenes",
                        metavar='',
                        default='yes',
                        help="intergenic regions is between gene body regions [yes/y], or between TSS and TES [no/n]?. default=yes")
    optional.add_argument("-l", "--minIntergenicLen",
                        metavar='',
                        default=2000,
                        type=int,
                        help="minimum intergenic region distance. default=2000")
    #added wang
    optional.add_argument("-xL","--maxIntergenicLen",
                          metavar='',
                          default=10000000,
                          type=int,
                          help="maximum intergenic region distance, default=10000000")
    optional.add_argument("-rem", "--removeShort",
                        metavar='',
                        default='yes',
                        help="should regions be removed from TSS, TES and 5distance regions if gene body is removed for its size being less than 0? [yes/y] or [no/n]. default=yes")
    optional.add_argument("-F", "--folderOut",
                        metavar='',
                        default='out',
                        help="what should be the name of the out-folder. default=out")
    return parser

def run(args):
    """Calls the main function with the unpacked arguments"""
    if args.folderOut[-1] == '/':
        args.folderOut = args.folderOut[:-1]

    if not os.path.exists(args.folderOut):
        os.system('mkdir '+args.folderOut)

    if not os.path.exists(args.folderOut+'/data'):
        os.system('mkdir '+args.folderOut+'/data')

    #added wang
    print("5' distance regions ", args.M, args.N)
    print("minimum , maximum length of intergenic region ", args.minIntergenicLen, args.maxIntergenicLen)

    main(args.referenceFile, args.genomeFile, args.human, args.numericalChr, args.removeMir,
         args.X, args.Y, args.M, args.N, args.intergenicBTGenes, args.minIntergenicLen,args.maxIntergenicLen,
         args.removeShort, args.folderOut)

if __name__ == "__main__":
    args = my_parser(argparse.ArgumentParser('python gene_annotation.py')).parse_args()

    run(args)
