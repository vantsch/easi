
#!/usr/bin/env python3
"""Summarise MAJIQ VOILA Modulizer TSVs for EASI.

Outputs:
  MAJIQ_summary.xlsx
  MAJIQ_gene_summary.tsv

Cassette sheet: one row per distinct alternative-exon interval within a gene.
Other classified sheets: one row per native MAJIQ event ID.
Junctions sheet: visible in Excel, excluded from gene-level event counts.

IMPORTANT: displayed PSI/dPSI are for the selected representative junction,
not a mathematically inferred whole-event PSI. The native event_changing flag
is reported separately from junction probability.
"""

from pathlib import Path
import pandas as pd

INPUT_DIR = Path('.')
OUTPUT_EXCEL = Path('MAJIQ_summary.xlsx')
OUTPUT_GENES = Path('MAJIQ_gene_summary.tsv')
PROBABILITY_THRESHOLD = 0.95
DIRECTION_DPSI_THRESHOLD = 0.10

EVENT_FILES = {
    'alt3prime', 'alt3and5prime', 'alt5prime',
    'alternate_first_exon', 'alternate_last_exon',
    'alternative_intron', 'cassette', 'junctions',
    'multi_exon_spanning', 'mutually_exclusive',
    'orphan_junction', 'tandem_cassette',
}
SKIP_FILES = {
    'heatmap', 'other', 'p_alt3prime', 'p_alt5prime',
    'p_alternate_first_exon', 'p_alternate_last_exon',
}
EXCLUDE_FROM_GENE_COUNT = {'junctions'}
SHEET_NAMES = {
    'alt3prime': "Alternative 3' SS",
    'alt3and5prime': "Alternative 3' and 5' SS",
    'alt5prime': "Alternative 5' SS",
    'alternate_first_exon': 'Alternative First Exon',
    'alternate_last_exon': 'Alternative Last Exon',
    'alternative_intron': 'Alternative Intron',
    'cassette': 'Exon Skipping',
    'junctions': 'Junctions',
    'multi_exon_spanning': 'Multi-Exon Spanning',
    'mutually_exclusive': 'Mutually Exclusive Exons',
    'orphan_junction': 'Orphan Junctions',
    'tandem_cassette': 'Tandem Exon Skipping',
}
SUMMARY_COLUMNS = [
    'gene_name', 'gene_id', 'seqid', 'strand', 'module_id',
    'module_event_combination', 'num_events', 'denovo_juncs',
]
INCLUSION_JUNCTIONS = {'C1_A', 'C2_A'}
SKIPPING_JUNCTIONS = {'C1_C2', 'C2_C1'}


def clean(value):
    if pd.isna(value):
        return ''
    s = str(value).strip()
    return '' if s.lower() in {'nan', 'none', 'na', '<na>'} else s


def chromosome(value):
    s = clean(value)
    return '' if not s else (s if s.startswith('chr') else 'chr' + s)


def coord(chrom, value):
    s = clean(value)
    if not s:
        return ''
    if ':' in s:
        c, pos = s.split(':', 1)
        return chromosome(c) + ':' + pos
    c = chromosome(chrom)
    return (c + ':' if c else '') + s


def lsv_coord(row):
    s = clean(row.get('lsv_id', ''))
    return coord(row.get('seqid', ''), s.rsplit(':', 1)[-1]) if s else ''


def junction_coord(row):
    return coord(row.get('seqid', ''), row.get('junction_coord', ''))


def bool_series(df, name):
    if name not in df.columns:
        return pd.Series(False, index=df.index)
    return df[name].astype(str).str.lower().str.strip().isin({'true', '1', 'yes'})


def joined(values):
    return '|'.join(sorted({clean(x) for x in values if clean(x)}))


def analysis_columns(df):
    psi = [x for x in df if x.endswith('_median_psi')]
    dpsi = [x for x in df if x.endswith('_median_dpsi')]
    prob = [x for x in df if x.endswith('_probability_changing')]
    if len(psi) < 2 or len(dpsi) != 1 or len(prob) != 1:
        return None
    return psi[:2], dpsi[0], prob[0]


def add_event_keys(df, event_type):
    work = df.copy()
    if event_type != 'cassette':
        work['alternative_exon_coord'] = ''
        work['canonical_event_key'] = (
            event_type + '|' + work['gene_id'].astype(str)
            + '|' + work['event_id'].astype(str)
        )
        return work

    required = {'spliced_with', 'spliced_with_coord', 'junction_name'}
    missing = required.difference(work.columns)
    if missing:
        raise ValueError(f'Cassette file missing columns: {sorted(missing)}')

    # The A rows describe the same alternative exon from either side.
    a = work.loc[work['junction_name'].isin(INCLUSION_JUNCTIONS)]
    exon_coords = a.groupby('event_id')['spliced_with_coord'].agg(
        lambda values: sorted({clean(v) for v in values if clean(v)})
    )
    valid = exon_coords[exon_coords.map(len).eq(1)].map(lambda x: x[0])
    work['alternative_exon_coord'] = work['event_id'].map(valid).fillna('')
    work['canonical_event_key'] = work.apply(
        lambda r: (
            f"cassette|{clean(r['gene_id'])}|{chromosome(r['seqid'])}|"
            f"{clean(r['strand'])}|{clean(r['alternative_exon_coord'])}"
            if clean(r['alternative_exon_coord'])
            else f"cassette-unresolved|{clean(r['gene_id'])}|{clean(r['event_id'])}"
        ), axis=1
    )
    unresolved = work.loc[work['alternative_exon_coord'].eq(''), 'event_id'].nunique()
    if unresolved:
        print(f'  Warning: {unresolved} cassette IDs lack unambiguous exon coordinates')
    return work


def choose_representative(group, dpsi_col, prob_col, cassette=False):
    candidates = group
    if cassette:
        inclusion = group.loc[group['junction_name'].isin(INCLUSION_JUNCTIONS)]
        if not inclusion.empty:
            candidates = inclusion
    valid = candidates.copy()
    valid['_dpsi'] = pd.to_numeric(valid[dpsi_col], errors='coerce')
    valid['_prob'] = pd.to_numeric(valid[prob_col], errors='coerce')
    valid = valid.dropna(subset=['_dpsi', '_prob'])
    if valid.empty:
        return None
    valid['_abs'] = valid['_dpsi'].abs()
    valid = valid.sort_values(['_prob', '_abs'], ascending=[False, False], kind='stable')
    return valid.iloc[0]


def cassette_direction(group, dpsi_col, prob_col):
    """Translate all supported cassette junction changes into exon direction.

    C1_A and C2_A: positive dPSI means more inclusion.
    C1_C2 and C2_C1: positive dPSI means more skipping.

    dPSI comparison must be LMNA minus WT for the descriptions below.
    Only junctions meeting BOTH thresholds vote on the direction.
    Multiple event IDs can reuse the same junction; count it once.
    """
    evidence = group.loc[
        group['junction_name'].isin(INCLUSION_JUNCTIONS | SKIPPING_JUNCTIONS)
    ].copy()
    evidence['_dpsi'] = pd.to_numeric(evidence[dpsi_col], errors='coerce')
    evidence['_prob'] = pd.to_numeric(evidence[prob_col], errors='coerce')
    evidence = evidence.loc[
        evidence['_dpsi'].abs().ge(DIRECTION_DPSI_THRESHOLD)
        & evidence['_prob'].ge(PROBABILITY_THRESHOLD)
    ].copy()
    if evidence.empty:
        return 'Insufficient directional evidence'

    # Keep alternative splice-graph roles separate: an identical genomic
    # junction could be labelled differently in different event structures.
    evidence = evidence.drop_duplicates(
        subset=['gene_id', 'seqid', 'strand', 'junction_name', 'junction_coord', '_dpsi']
    )
    roles = evidence['junction_name'].isin(INCLUSION_JUNCTIONS)
    biological_sign = evidence['_dpsi'].where(roles, -evidence['_dpsi'])
    if biological_sign.gt(0).all():
        return 'More exon inclusion'
    if biological_sign.lt(0).all():
        return 'More exon skipping'
    return 'Discordant junction evidence'


def describe_change(row, event_type, dpsi_col):
    d = float(row[dpsi_col])
    if d == 0:
        return 'No directional change'
    prefix = 'Increased' if d > 0 else 'Decreased'
    labels = {
        'alt3prime': "alternative 3' splice-site junction usage",
        'alt5prime': "alternative 5' splice-site junction usage",
        'alt3and5prime': 'alternative splice-site junction usage',
        'alternative_intron': 'intron-associated junction usage',
        'alternate_first_exon': 'first-exon-associated junction usage',
        'alternate_last_exon': 'last-exon-associated junction usage',
        'tandem_cassette': 'tandem-cassette junction usage',
        'mutually_exclusive': 'selected exon-associated junction usage',
        'junctions': 'junction usage',
        'orphan_junction': 'junction usage',
        'multi_exon_spanning': 'multi-exon-spanning junction usage',
    }
    return f"{prefix} {labels.get(event_type, 'event-associated junction usage')}"


def summarise(df, event_type):
    cols = analysis_columns(df)
    if cols is None or not {'event_id', 'gene_id', 'seqid', 'strand', 'junction_coord'}.issubset(df):
        return pd.DataFrame()
    psi_cols, dpsi_col, prob_col = cols
    work = add_event_keys(df, event_type)
    rows = []
    for key, group in work.groupby('canonical_event_key', dropna=False):
        rep = choose_representative(group, dpsi_col, prob_col, event_type == 'cassette')
        if rep is None:
            continue
        dpsi = float(rep[dpsi_col])
        prob = float(rep[prob_col])
        junctions = sorted({coord(rep['seqid'], x) for x in group['junction_coord'] if clean(x)})
        if event_type == 'cassette':
            exon = clean(rep.get('alternative_exon_coord', ''))
            coordinates = coord(rep['seqid'], exon) if exon else lsv_coord(rep)
            change = cassette_direction(group, dpsi_col, prob_col)
        else:
            coordinates = lsv_coord(rep)
            change = describe_change(rep, event_type, dpsi_col)

        # A group-level 'Significant' means at least one junction exceeds
        # the posterior probability criterion, NOT native event significance.
        probs = pd.to_numeric(group[prob_col], errors='coerce')
        is_supported = bool(probs.ge(PROBABILITY_THRESHOLD).any())
        rows.append({
            'Gene': clean(rep.get('gene_name', '')),
            'Gene ID': clean(rep.get('gene_id', '')),
            'Chr': chromosome(rep.get('seqid', '')),
            'Strand': clean(rep.get('strand', '')),
            'Coordinates': coordinates,
            'Junction coordinates': junction_coord(rep),
            'Canonical Event ID': key,
            'Event ID': clean(rep.get('event_id', '')),
            'All MAJIQ Event IDs': joined(group['event_id']),
            'Module IDs': joined(group['module_id']) if 'module_id' in group else '',
            'Change': change,
            psi_cols[0]: round(float(rep[psi_cols[0]]), 3) if pd.notna(rep[psi_cols[0]]) else None,
            psi_cols[1]: round(float(rep[psi_cols[1]]), 3) if pd.notna(rep[psi_cols[1]]) else None,
            'dPSI': round(dpsi, 3),
            '|ΔPSI|': round(abs(dpsi), 3),
            'Probability changing': round(prob, 4),
            'Junction changing': 'Yes' if bool_series(group, 'junction_changing').any() else 'No',
            'Event changing': 'Yes' if bool_series(group, 'event_changing').any() else 'No',
            'Significant': 'Yes' if is_supported else 'No',
            'Complex': 'Yes' if bool_series(group, 'complex').any() else 'No',
            'Number of MAJIQ IDs': group['event_id'].nunique(),
            'Number of unique junctions': len(junctions),
            'All junctions': '|'.join(junctions),
        })
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return result.sort_values(
        ['Significant', '|ΔPSI|', 'Gene'], ascending=[False, False, True]
    ).reset_index(drop=True)


def prepare_summary(df):
    columns = [c for c in SUMMARY_COLUMNS if c in df]
    if not columns:
        return pd.DataFrame()
    out = df[columns].copy()
    if 'module_event_combination' in out:
        out['module_event_combination'] = (
            out['module_event_combination'].astype(str)
            .str.replace('tandem_cassette', 'tandem exon skipping', regex=False)
            .str.replace('cassette', 'exon skipping', regex=False)
        )
    return out


def determine_effect_direction(group):
    significant = group.loc[group['significant']]
    selected = significant if not significant.empty else group
    effects = pd.to_numeric(selected['dPSI'], errors='coerce').dropna()
    effects = effects.loc[effects.ne(0)]
    if effects.empty:
        return 'none'
    if effects.gt(0).any() and effects.lt(0).any():
        return 'mixed'
    return 'increased' if effects.gt(0).any() else 'decreased'


def prepare_gene_summary(event_tables):
    columns = [
        'gene_id', 'gene_name', 'method', 'detected', 'significant',
        'n_events', 'n_significant_events', 'n_majiq_changing_events',
        'max_abs_effect', 'effect_direction', 'event_types',
    ]
    if not event_tables:
        return pd.DataFrame(columns=columns)
    events = pd.concat(event_tables, ignore_index=True)
    if events.empty:
        return pd.DataFrame(columns=columns)
    events['is_significant'] = events['Significant'].eq('Yes')
    events['is_event_changing'] = events['Event changing'].eq('Yes')
    rows = []
    for (gene_id, gene_name, canonical), group in events.groupby(
        ['Gene ID', 'Gene', 'Canonical Event ID'], dropna=False
    ):
        effects = pd.to_numeric(group['dPSI'], errors='coerce')
        if effects.notna().any():
            idx = effects.abs().idxmax()
            effect = float(effects.loc[idx])
            max_effect = float(effects.abs().max())
        else:
            effect, max_effect = 0.0, 0.0
        rows.append({
            'Gene ID': gene_id, 'Gene': gene_name,
            'significant': bool(group['is_significant'].any()),
            'event_changing': bool(group['is_event_changing'].any()),
            'dPSI': effect, 'max_abs_effect': max_effect,
            'event_types': joined(group['event_type']),
        })
    collapsed = pd.DataFrame(rows)
    output = []
    for (gene_id, gene_name), group in collapsed.groupby(['Gene ID', 'Gene'], dropna=False):
        effect_max = pd.to_numeric(group['max_abs_effect'], errors='coerce').max()
        types = sorted({x for entry in group['event_types'] for x in entry.split('|') if x})
        output.append({
            'gene_id': gene_id, 'gene_name': gene_name,
            'method': 'MAJIQ', 'detected': True,
            'significant': bool(group['significant'].any()),
            'n_events': len(group),
            'n_significant_events': int(group['significant'].sum()),
            'n_majiq_changing_events': int(group['event_changing'].sum()),
            'max_abs_effect': round(float(effect_max), 3) if pd.notna(effect_max) else None,
            'effect_direction': determine_effect_direction(group),
            'event_types': '|'.join(types),
        })
    result = pd.DataFrame(output, columns=columns)
    return result.sort_values(
        ['significant', 'n_significant_events', 'max_abs_effect', 'gene_name'],
        ascending=[False, False, False, True]
    ).reset_index(drop=True)


def main():
    sheets = {}
    event_tables = []
    for file in sorted(INPUT_DIR.glob('*.tsv')):
        name = file.stem.lower()
        if name in SKIP_FILES or (name != 'summary' and name not in EVENT_FILES):
            continue
        print(f'Processing {file.name}')
        df = pd.read_csv(file, sep='\t', comment='#')
        if df.empty:
            continue
        if name == 'summary':
            table = prepare_summary(df)
            if not table.empty:
                sheets['Summary'] = table
            continue
        table = summarise(df, name)
        if table.empty:
            print('  Skipped: no reportable events')
            continue
        print(f"  Native MAJIQ IDs: {df['event_id'].nunique()} -> reported rows: {len(table)}")
        table['event_type'] = name
        sheets[SHEET_NAMES[name][:31]] = table.drop(columns=['event_type'])
        if name not in EXCLUDE_FROM_GENE_COUNT:
            event_tables.append(table)
    if not sheets:
        raise RuntimeError('No supported MAJIQ TSV files found.')
    with pd.ExcelWriter(OUTPUT_EXCEL, engine='openpyxl') as writer:
        for name, sheet in sheets.items():
            sheet.to_excel(writer, sheet_name=name, index=False)
            ws = writer.sheets[name]
            ws.freeze_panes = 'A2'
            ws.auto_filter.ref = ws.dimensions
            for cells in ws.columns:
                letter = cells[0].column_letter
                longest = max(len(str(c.value)) if c.value is not None else 0 for c in cells[:200])
                ws.column_dimensions[letter].width = min(max(longest + 2, 12), 45)
    genes = prepare_gene_summary(event_tables)
    genes.to_csv(OUTPUT_GENES, sep='\t', index=False)
    print(f'Finished: {OUTPUT_EXCEL}; {OUTPUT_GENES}; {len(genes)} genes')


if __name__ == '__main__':
    main()
