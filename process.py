import pandas as pd

def serial_dilution(start_value, dilution_fract, num_dilutions):
    concentrations = []
    for i in range(1,num_dilutions):
        if i == 1:
            new_value = start_value
        else:
            new_value = concentrations[i-2] * dilution_fract
        concentrations.append(new_value)
    conc_strs = [str(val) for val in concentrations]
    return conc_strs

def normalize_lowest_and_reorder(df):    
    cols = [float(i) for i in list(df.columns)[:-2]]
    min_col = str(min(cols))
    min_col_mean = df[min_col].mean()
    vel_mean = df['velcade'].mean()
    df_norm = (df-vel_mean)/(min_col_mean-vel_mean)*100

    # Reorder low to high
    columns = df_norm.columns.tolist()
    columns = columns[::-1]
    df_norm = df_norm[columns]

    return df_norm

def extract_f114_and_dmso_data(input_df):
    suffixes = ('GFP+', 'PI-A- , GFP-A+', 'GFP-A+ , PI-A-')
    f114 = []
    dmso = []
    for n, sample in enumerate(input_df['Name']):
        if 'F114' in sample:
            if sample.endswith(suffixes):
                f114.append(input_df['Statistic'][n])
        if 'DMSO' in sample:
            if sample.endswith(suffixes):
                dmso.append(input_df['Statistic'][n])

    df = pd.DataFrame()
    df['dmso'] = dmso
    df['f114'] = f114
    return df

def normalize_to_dmso(df):
    dmso_mean = df['dmso'].mean()
    df_norm = df/dmso_mean
    return df_norm