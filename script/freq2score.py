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
    df[rep+"_binding_score"] = (binding_weight/binding_norm_factor)/df[rep+"_exp_score"]
    return df 


def main():
  inputfile = 'result/mut_freq_2*10-5.tsv'
  outfile = 'result/mut_scores_2*10-5.tsv'
  freq_cutoff = 0.00001
  freq_df = pd.read_csv(inputfile, sep = '\t')

  
  freq_df = get_expression_score(freq_df, 'Rep1')
  freq_df = get_expression_score(freq_df, 'Rep2')

  freq_df =  get_expression_pos_neg(freq_df, 'Rep1')
  freq_df =  get_expression_pos_neg(freq_df, 'Rep2')

  freq_df = get_binding_score(freq_df, 'Rep1')
  freq_df = get_binding_score(freq_df, 'Rep2')  

  cols = ['SH_pep', 'mut_ID', 
          'Rep1_exp_score', 'Rep1_exp_pos/neg', 'Rep2_exp_score', 'Rep2_exp_pos/neg', 'Rep1_binding_score', 'Rep2_binding_score']
  
  freq_df = freq_df[cols]
  freq_df = freq_df.fillna(0)
  freq_df.to_csv(outfile, sep="\t", index = False)



if __name__ == "__main__":
  main()
