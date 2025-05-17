import pandas as pd
from rdkit import Chem

df = pd.read_csv('prep1.csv')

df['Molecule'] = df['Chromophore'].apply(lambda x: Chem.MolFromSmiles(str(x)))
df.dropna(subset=['Molecule'], inplace=True)

df_solvent_count = (
    df.groupby('Chromophore')['Solvent']
      .nunique()
      .reset_index(name='NumSolvents')
)

df_solvent_count['NumSolvents'] = df_solvent_count['NumSolvents'].astype(int)

output_path = '/Users/diana/Desktop/RDKit/chromophore_solvent_count.csv'
df_solvent_count.to_csv(output_path, index=False)

print(f"Saved CSV with {len(df_solvent_count)} unique chromophores to:")
print(output_path)
