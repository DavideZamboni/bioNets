import numpy as np
import pandas as pd
# Rec
'''true_vs_computed = {
    'HBonds':[
        '/home/davide/storage/lab_dati/risultati_paragone/rec/hb-graph_1micros.dat', 
        '/home/davide/storage/lab_dati/psns_out/cli_outs/rec_grk1/Adjacency/named_interactors/Hbonds_adjacency.dat'
              ],
    'HC':[
        '/home/davide/storage/lab_dati/risultati_paragone/rec/hc-graph_1micros.dat', 
        '/home/davide/storage/lab_dati/psns_out/cli_outs/rec_grk1/Adjacency/named_interactors/HC_adjacency.dat'
              ],
    'SB':[
        '/home/davide/storage/lab_dati/risultati_paragone/rec/sb-graph_1micros.dat', 
        '/home/davide/storage/lab_dati/psns_out/cli_outs/rec_grk1/Adjacency/named_interactors/SB_adjacency.dat'
              ]
}
'''
# Pyinteraph 1
true_vs_computed = {
    'HBonds':[
        '/home/davide/storage/lab_dati/risultati_paragone/pyinter_test/hb-graph.dat', 
        '/home/davide/storage/lab_dati/psns_out/cli_outs/pyinteraph_testcase/Adjacency/named_interactors/Hbonds_adjacency.dat'
              ],
    'HC':[
        '/home/davide/storage/lab_dati/risultati_paragone/pyinter_test/hc-graph.dat', 
        '/home/davide/storage/lab_dati/psns_out/cli_outs/pyinteraph_testcase/Adjacency/named_interactors/HC_adjacency.dat'
              ],
    'SB':[
        '/home/davide/storage/lab_dati/risultati_paragone/pyinter_test/sb-graph.dat', 
        '/home/davide/storage/lab_dati/psns_out/cli_outs/pyinteraph_testcase/Adjacency/named_interactors/SB_adjacency.dat'
              ]
}


sele_rows   = []   # one of (truth, computed) is zero, other not
params_rows = []   # both non-zero, but values differ
op = '/home/davide/storage/lab_dati/'
for inttype, (truth_p, comp_p) in true_vs_computed.items():
    true_tot = 0
    computed_tot = 0
    diff = 0

    true     = np.loadtxt(truth_p)
    computed = np.loadtxt(comp_p)

    if true.shape != computed.shape:
        print(f'{inttype}: SHAPE MISMATCH {true.shape} vs {computed.shape}')
        continue

    if np.allclose(true, computed, atol=1e-2):
        print(f'{inttype}: OK')
        continue

    diff_mask  = ~np.isclose(true, computed, atol=1e-2)
    same_mask  = np.isclose(true, computed, atol=1e-2)
    same_rows, same_cols = np.where(same_mask)
    rows, cols = np.where(diff_mask)

    for r, c in zip(rows, cols):
        if r <= c: 
            continue

        t, k = true[r, c], computed[r, c]
        mydiff = abs(k - t)
        diff += mydiff
        true_tot += t
        computed_tot += k


        t_zero = np.equal(t, 0.0)
        k_zero = np.equal(k, 0.0)

        if t_zero ^ k_zero:                       # XOR → exactly one is zero
            if t_zero:
                present = 'computed'
            else:
                present = 'true'
            sele_rows.append({
                'row': r,
                'col': c,
                'intertype': inttype,
                'present_in': present,
                'non-0 val': mydiff
            })
        elif not t_zero and not k_zero:           # both non-zero, values differ
            params_rows.append({
                'row': r,
                'col': c,
                'intertype': inttype,
                'difference': t - k,
            })
    
    for r, c in zip(same_rows, same_cols):
        if r <= c: 
            continue
        t, k = true[r, c], computed[r, c]
        mydiff = abs(k - t)
        diff += mydiff
        true_tot += t
        computed_tot += k

    print(f'  true:     min={true.min():.4g}  max={true.max():.4g}  mean(>0)={true[true>0].mean():.4g}')
    print(f'  computed: min={computed.min():.4g}  max={computed.max():.4g}  mean(>0)={computed[computed>0].mean():.4g}')

    print(f'{inttype}: {len(rows)} differing entries - {len(same_rows)} similar entries')
    print(f'True total: {true_tot}, computed total: {computed_tot} observed difference abs(computed - true): {diff}\n% difference: {(diff / true_tot) * 100}')


pd.DataFrame(sele_rows,
             columns=['row', 'col', 'intertype', 'present_in', 'non-0 val']
            ).to_csv( op + 'possible_sele_mismatch.csv', index=False)

pd.DataFrame(params_rows,
             columns=['row', 'col', 'intertype', 'difference']
            ).to_csv( op + 'possible_params_mismatch.csv', index=False)

