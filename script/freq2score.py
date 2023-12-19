import glob
import pandas as pd

def get_expression_score(df, rep):
    exp_weight = df[rep+'_t0_gate1_freq']*0.25 + df[rep+'_t0_gate2_freq']*0.5 + \
                 df[rep+'_t0_gate3_freq']*0.75  + df[rep+'_t0_gate4_freq']*1
    exp_norm_factor = df[rep+'_t0_gate1_freq'] + df[rep+'_t0_gate2_freq'] + \
                 df[rep+'_t0_gate3_freq'] + df[rep+'_t0_gate4_freq']
    df[rep+"_exp_score"] = exp_weight/exp_norm_factor
    return df 

def get_expression_pos_neg(df, rep):
    df[rep+'_exp_pos/neg'] = df[rep+'_exp_pos_freq']/df[rep+'_exp_neg_freq']
    return df

def get_binding_score(df, rep):
    binding_weight = df[rep+'_t1_gate1_freq']*0.25 + df[rep+'_t1_gate2_freq']*0.5 + \
                 df[rep+'_t1_gate3_freq']*0.75  + df[rep+'_t1_gate4_freq']*1
    binding_norm_factor = df[rep+'_t1_gate1_freq'] + df[rep+'_t1_gate2_freq'] + \
                 df[rep+'_t1_gate3_freq'] + df[rep+'_t1_gate4_freq']
    df[rep+"_binding_score"] = binding_weight/binding_norm_factor
    return df 

def get_norm_binding_score(df, rep):
    neg_ctrl_pep = 'DSAKEALDKYFKNH'
    sars2_pep    = 'DSFKEELDKYFKNH'
    s2p6_wt = 'QIVHLL-MMRN'
    neg_ctrl_df = df[(df['SH_pep'] == neg_ctrl_pep) & (df['totalfreq']>0.00006)]
    avg_neg_ctrl_binding_score = neg_ctrl_df[rep + '_binding_score'].mean()
    print(avg_neg_ctrl_binding_score)
    s2p6_wt_to_sars2_binding_df = df[(df['SH_pep'] == sars2_pep) & (df['mut_ID'] == s2p6_wt)]
    s2p6_wt_to_sars2_binding_df = s2p6_wt_to_sars2_binding_df.reset_index()
    print(s2p6_wt_to_sars2_binding_df)
    s2p6_wt_to_sars2_binding_score = s2p6_wt_to_sars2_binding_df.loc[0, rep+'_binding_score']
    print(s2p6_wt_to_sars2_binding_score)
    df[rep + '_norm_binding_score'] = (df[rep + '_binding_score'] - avg_neg_ctrl_binding_score) / (s2p6_wt_to_sars2_binding_score - avg_neg_ctrl_binding_score)
    return df

def get_norm_binding_avg(df):
    df['avg_norm_binding_score'] = (df['Rep1_norm_binding_score'] + df['Rep2_norm_binding_score'])/2
    return df

def main():
  inputfile = 'result/mut_freq.tsv'
  outfile = 'result/mut_scores.tsv'
  freq_df = pd.read_csv(inputfile, sep = '\t')

  
  freq_df = get_expression_score(freq_df, 'Rep1')
  freq_df = get_expression_score(freq_df, 'Rep2')

  freq_df =  get_expression_pos_neg(freq_df, 'Rep1')
  freq_df =  get_expression_pos_neg(freq_df, 'Rep2')

  freq_df = get_binding_score(freq_df, 'Rep1')
  freq_df = get_binding_score(freq_df, 'Rep2')  

  freq_df = get_norm_binding_score(freq_df, 'Rep1')
  freq_df = get_norm_binding_score(freq_df, 'Rep2')
  freq_df = get_norm_binding_avg(freq_df)

  print(len(freq_df))

  cols = ['SH_pep', 'mut_ID', 'totalfreq', 'Rep1_t1_freq', 'Rep2_t1_freq',
          'Rep1_t1_gate1_count','Rep1_t1_gate2_count','Rep1_t1_gate3_count','Rep1_t1_gate4_count',
          'Rep2_t1_gate1_count','Rep2_t1_gate2_count','Rep2_t1_gate3_count','Rep2_t1_gate4_count',
	  'Rep1_exp_score', 'Rep1_exp_pos/neg', 'Rep2_exp_score', 'Rep2_exp_pos/neg', 'Rep1_binding_score', 'Rep2_binding_score',
	  'Rep1_norm_binding_score', 'Rep2_norm_binding_score', 'avg_norm_binding_score']
  
  freq_df = freq_df[cols]
  freq_df = freq_df.fillna(0)
  freq_df.to_csv(outfile, sep="\t", index = False)



if __name__ == "__main__":
  main()
