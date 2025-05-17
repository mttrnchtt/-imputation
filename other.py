import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('prep1.csv')

df = df.rename(columns={
    'Chromophore': 'Chromophore_SMILES',
    'Absorption max (nm)': 'Absorption_nm',
    'Emission max (nm)': 'Emission_nm'
})

df['Molecule'] = df['Chromophore_SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(str(smi)) if pd.notnull(smi) else None
)
invalid = df[df['Molecule'].isnull()]
print(f"Found {len(invalid)} invalid SMILES examples:")
print(invalid['Chromophore_SMILES'].head().to_list())

df = df.dropna(subset=['Molecule']).reset_index(drop=True)

def analyze_exclusivity(df):
    solvent_counts = (
        df.groupby('Chromophore_SMILES')['Solvent']
          .nunique()
          .reset_index(name='num_solvents')
    )
    exclusive_chromos = solvent_counts[solvent_counts['num_solvents'] == 1]
    sol_stats = (
        df[df['Chromophore_SMILES'].isin(exclusive_chromos['Chromophore_SMILES'])]
        .groupby('Solvent')['Chromophore_SMILES']
        .nunique()
        .reset_index(name='exclusive_chromophores')
        .sort_values('exclusive_chromophores', ascending=False)
    )
    return sol_stats

exclusive_solvents = analyze_exclusivity(df)
print("\nTop-5 solvents by number of exclusive chromophores:")
print(exclusive_solvents.head())

def add_molecular_features(df):
    df['MolWeight'] = df['Molecule'].apply(Descriptors.MolWt)
    df['LogP']      = df['Molecule'].apply(Descriptors.MolLogP)
    df['NumAtoms']  = df['Molecule'].apply(lambda m: m.GetNumAtoms())
    return df

df = add_molecular_features(df)

plt.figure(figsize=(12, 6))
sns.scatterplot(
    data=df,
    x='Absorption_nm',
    y='Emission_nm',
    hue='Solvent',
    palette='tab10',
    s=80,
    edgecolor='w'
)
plt.title("Chromophore Absorption vs. Emission Maxima")
plt.xlabel("Absorption max (nm)")
plt.ylabel("Emission max (nm)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

def detect_anomalies(df):
    spec_anoms = df[(df['Absorption_nm'] < 300) | (df['Absorption_nm'] > 500)]
    mw_stats = df['MolWeight'].describe()
    mw_anoms = df[(df['MolWeight'] < mw_stats['25%']) | (df['MolWeight'] > mw_stats['75%'])]
    return spec_anoms, mw_anoms

spec_anoms, mw_anoms = detect_anomalies(df)
print("\nSpectral anomalies (first 5):")
print(spec_anoms[['Chromophore_SMILES','Absorption_nm','Emission_nm']].head())
print("\nMolecular weight anomalies (first 5):")
print(mw_anoms[['Chromophore_SMILES','MolWeight']].head())
