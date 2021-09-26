import numpy as np
import pandas as pd

df = pd.DataFrame(np.arange(6).reshape(3,2), columns=list('AB'))
df[df['B'] == 3].A
# 1    2
# Name: A, dtype: int64

print(df.loc[df['B'] == 3, 'A'])
print(df)
