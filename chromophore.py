import pandas as pd
from rdkit import Chem

def shorten_string(s, max_len=10):
    return s if len(s) <= max_len else s[:max_len] + '...'

df = pd.read_csv('/prep1.csv')
df['Molecule'] = df['Chromophore'].apply(lambda x: Chem.MolFromSmiles(str(x)))
df.dropna(subset=['Molecule'], inplace=True)

total_unique_solvents = df['Solvent'].nunique()

df_stats = (
    df.groupby('Chromophore')['Solvent']
      .agg(total_count='size', unique_solvents='nunique')
      .reset_index()
)
df_stats.sort_values('total_count', ascending=False, inplace=True)

top10 = df_stats.head(10).reset_index(drop=True)

df_unique = top10[['Chromophore', 'unique_solvents']].copy()
df_unique.rename(columns={'unique_solvents': 'UniqueSolvents'}, inplace=True)
df_unique['UniquePct'] = (df_unique['UniqueSolvents'] / total_unique_solvents * 100).round(1)

accumulated = set()
coverage = []
for chromo in top10['Chromophore']:
    solvents = set(df.loc[df['Chromophore'] == chromo, 'Solvent'])
    accumulated.update(solvents)
    coverage.append({
        'Chromophore': chromo,
        'CoverageSolvents': len(accumulated),
        'CoveragePct': round(len(accumulated) / total_unique_solvents * 100, 1)
    })
df_coverage = pd.DataFrame(coverage)

df_merged = pd.merge(
    df_unique,
    df_coverage,
    on='Chromophore',
    how='inner'
)[['Chromophore', 'UniqueSolvents', 'UniquePct', 'CoverageSolvents', 'CoveragePct']]

df_merged['Chromophore'] = df_merged['Chromophore'].apply(lambda x: shorten_string(x, 10))

print("=== Merged Top-10 Chromophores Table ===")
print(df_merged.to_string(index=False))
