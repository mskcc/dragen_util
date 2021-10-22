import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.pyplot import figure 
import glob

figure(figsize=(16,8), dpi=80)
TARGET_X_AXIS = 'control_mapq'
TARGET_Y_AXIS = 'ct_diff'

X_LABEL = "BWA-MEM MAPQ"
Y_LABEL = "Median( BWA MEM MAPQ - DRAGEN MAPQ )"

files = glob.glob('*_bam_differences.csv')

for csv_f in files:
    basename = csv_f.split('.')[0]

    df = pd.read_csv(csv_f)

    grouped = df.groupby(TARGET_X_AXIS).median()

    # Write intermediate grouping results to file
    grouped_filename = "{}___grouped.csv".format(basename)
    grouped_data = grouped.to_csv(index=True)
    grouped_file = open(grouped_filename, "a")
    grouped_file.write(grouped_data)
    grouped_file.close()

    # Graph intermediate results as a bar graph
    grouped_df = pd.read_csv(grouped_filename)
    plt.bar(grouped_df[TARGET_X_AXIS], grouped_df[TARGET_Y_AXIS], color='#1f77b4') 
    plt.title("{} vs. {}".format(X_LABEL, Y_LABEL))
    plt.xlabel(X_LABEL)                    
    plt.ylabel(Y_LABEL)
    plt.xlim([-2, 62])                          # Helpful for plotting MAPQ, which has range [1-60]

    plt.savefig("{}.pdf".format(basename))
    plt.close()
