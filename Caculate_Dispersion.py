import numpy as np
import pandas as pd
# calculate the dispersion
# wavelength: 波长；D：色散参量；L：纤芯长度(km)
def Dispersion(wavelength, D, L):
    c = 300000
    beta = -D*(wavelength**2/(2*np.pi*c))
    beta_final = beta*L
    return beta_final, beta


beta2_EDF_end, beta2_EDF = Dispersion(1550, -22, 0.0005)
beat_SMF_end, beta2_SMF = Dispersion(1550, 15, 0.015)
beta2_DCF_end, beta2_DCF = Dispersion(1550, -58, 0.0003)
beta2_DSF_end, beta2_DSF = Dispersion(1550, -4, 0.050)

dispersion_data = np.ones((4, 2))
dispersion_data[0, 0] = beta2_EDF_end
dispersion_data[0, 1] = beta2_EDF
dispersion_data[1, 0] = beat_SMF_end
dispersion_data[1, 1] = beta2_SMF
dispersion_data[2, 0] = beta2_DCF_end
dispersion_data[2, 1] = beta2_DCF
dispersion_data[3, 0] = beta2_DSF_end
dispersion_data[3, 1] = beta2_DSF

dispersion_data = pd.DataFrame(dispersion_data, columns=['二阶色散(ps²)', '二阶色散系数(ps²/km)'])
name_fiber = ['EDF', 'SMF', 'DCF', 'DSF']
name_pd = pd.DataFrame(name_fiber, columns=['纤芯名称'])
dispersion_data = pd.concat([name_pd, dispersion_data], axis=1)
print('二阶色散数据如下：\n', dispersion_data)
final_dispersion = np.sum(dispersion_data['二阶色散(ps²)'])
print('最终二阶色散为：', final_dispersion, 'ps²')
dispersion_data['总色散(ps²)'] = final_dispersion
pd.DataFrame.to_csv(dispersion_data, 'd:/Onedrive/OneDrive - CQU/Code/Python/Active_Regenerative_Mode_Lock/data/dispersion_data.csv')