import numpy as np
import pandas as pd
true_vs_computed = {
    'HBonds':[
        '/home/davide/storage/lab_dati/risultati_rec/hb-graph_1micros.dat', 
        '/home/davide/storage/lab_dati/psns_out/cli_outs/rec_grk1/Adjacency/named/Hbonds_adjacency.dat'
              ],
    'HC':[
        '/home/davide/storage/lab_dati/risultati_rec/hc-graph_1micros.dat', 
        '/home/davide/storage/lab_dati/psns_out/cli_outs/rec_grk1/Adjacency/named/HC_adjacency.dat'
              ],
    'SB':[
        '/home/davide/storage/lab_dati/risultati_rec/sb-graph_1micros.dat', 
        '/home/davide/storage/lab_dati/psns_out/cli_outs/rec_grk1/Adjacency/named/SB_adjacency.dat'
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

    print(f'{inttype}: {len(rows)} differing entries')
    print(f'True total: {true_tot}, coomputed total: {computed_tot} observed difference abs(computed - true): {diff}\n% difference: {(diff / true_tot) * 100}')


pd.DataFrame(sele_rows,
             columns=['row', 'col', 'intertype', 'present_in', 'non-0 val']
            ).to_csv( op + 'possible_sele_mismatch.csv', index=False)

pd.DataFrame(params_rows,
             columns=['row', 'col', 'intertype', 'difference']
            ).to_csv( op + 'possible_params_mismatch.csv', index=False)

