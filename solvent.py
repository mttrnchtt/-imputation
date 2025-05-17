import pandas as pd
from rdkit import Chem
import plotly.graph_objects as go
from plotly.subplots import make_subplots

df = pd.read_csv('prep1.csv')
df['Molecule'] = df['Chromophore'].apply(lambda x: Chem.MolFromSmiles(str(x)))
df.dropna(subset=['Molecule'], inplace=True)

total_rows = len(df)
total_unique_chromos = df['Chromophore'].nunique()

df_stats = (
    df.groupby('Solvent')['Chromophore']
      .agg(total_count='size', unique_chromophores='nunique')
      .reset_index()
)
df_stats.sort_values('total_count', ascending=False, inplace=True)
top10 = df_stats.head(10).reset_index(drop=True)

df_total = top10[['Solvent', 'total_count']].copy()
df_total.rename(columns={'total_count': 'Count'}, inplace=True)
df_total['Percentage'] = (df_total['Count'] / total_rows * 100).round(1)

print("=== Top-10 Solvents by Total Occurrences ===")
print(df_total[['Solvent', 'Count', 'Percentage']])

df_unique = top10[['Solvent', 'unique_chromophores']].copy()
df_unique.rename(columns={'unique_chromophores': 'Count'}, inplace=True)
df_unique['Percentage'] = (df_unique['Count'] / total_unique_chromos * 100).round(1)

print("\n=== Top-10 Solvents by Unique Chromophores ===")
print(df_unique[['Solvent', 'Count', 'Percentage']])

accumulated = set()
coverage = []
for solvent in df_unique['Solvent']:
    chromos = set(df.loc[df['Solvent'] == solvent, 'Chromophore'])
    accumulated.update(chromos)
    cnt = len(accumulated)
    pct = round(cnt / total_unique_chromos * 100, 1)
    coverage.append({'Solvent': solvent, 'Count': cnt, 'Percentage': pct})

df_coverage = pd.DataFrame(coverage)

print("\n=== Cumulative Coverage by Top-10 Solvents ===")
print(df_coverage[['Solvent', 'Count', 'Percentage']])

df_unique_renamed = df_unique.rename(
    columns={'Count': 'UniqueCount', 'Percentage': 'UniquePct'}
)
df_coverage_renamed = df_coverage.rename(
    columns={'Count': 'CoverageCount', 'Percentage': 'CoveragePct'}
)

df_merged = pd.merge(
    df_unique_renamed,
    df_coverage_renamed,
    on='Solvent',
    how='inner'
)

print("\n=== Merged Top-10 Solvents Table ===")
print(df_merged[['Solvent', 'UniqueCount', 'UniquePct', 'CoverageCount', 'CoveragePct']])

fig = make_subplots(specs=[[{"secondary_y": True}]])

fig.add_trace(
    go.Bar(
        x=df_unique['Solvent'],
        y=df_unique['Percentage'],
        name='% Unique Chromophores',
        text=df_unique['Count'],
        textposition='outside'
    ),
    secondary_y=False
)

fig.add_trace(
    go.Scatter(
        x=df_coverage['Solvent'],
        y=df_coverage['Percentage'],
        name='Cumulative Coverage %',
        mode='lines+markers'
    ),
    secondary_y=True
)

fig.update_layout(
    title="Top-10 Solvents: Unique Chromophores & Cumulative Coverage",
    xaxis_title="Solvent",
    legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1)
)
fig.update_yaxes(title_text="% Unique Chromophores", secondary_y=False)
fig.update_yaxes(title_text="% Cumulative Coverage", secondary_y=True)

fig.show()

df_stats['PercentageUnique'] = (df_stats['unique_chromophores'] / total_unique_chromos * 100).round(1)

accumulated = set()
coverage_all = []
for solvent in df_stats['Solvent']:
    chromos = set(df.loc[df['Solvent'] == solvent, 'Chromophore'])
    accumulated.update(chromos)
    cnt = len(accumulated)
    pct = round(cnt / total_unique_chromos * 100, 1)
    coverage_all.append({'Solvent': solvent, 'CoverageCount': cnt, 'CoveragePct': pct})

df_coverage_all = pd.DataFrame(coverage_all)
df_full = pd.merge(
    df_stats[['Solvent', 'PercentageUnique']],
    df_coverage_all[['Solvent', 'CoverageCount', 'CoveragePct']],
    on='Solvent',
    how='outer'
)

df_full.to_csv('/Users/diana/Desktop/RDKit/full_data_all_solvents.csv', index=False)

print("\n=== Full Data (All Solvents) Saved to CSV ===")
