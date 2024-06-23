### usage: compute_consensus.py /path/to/out/dir

import sys
import result_converter as rc

out_dir=sys.argv[1]

# collect results from multiple tools in dataframe
results_df = rc.combine_results(out_dir)

# if exist drop columns
if 'DQA1_1' in results_df:
    results_df = results_df.drop(columns=['DQA1_1'])
if 'DQA1_2' in results_df:
    results_df = results_df.drop(columns=['DQA1_2'])

# create text files with hla-allele input for netMHC
rc.input_netMHC(results_df, out_dir)

# create dataframe with additional consensus row
df_with_cons = rc.results_with_consensus(results_df, out_dir)