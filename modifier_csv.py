import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

la, lb, lc = 0.27, 0.315, 0.08
def jacobien(theta1, theta4):
    D = (la + lb + lc / 2) / 3
    r1, r2 = la / D, lb / D
    s1, s4 = np.sin(theta1), np.sin(theta4)
    c1, c4 = np.cos(theta1), np.cos(theta4)
    A = r1 * s1 - r1 * s4
    B = 2 * (lc / 2) / D + r1 * c1 + r1 * c4
    C = np.pi / 2 + np.arctan(A / B)
    D_val = 8 * r2 * np.sqrt(1 - (B * 2 + A * 2) / (4 * r2 ** 2))
    E = r2 * np.sin(C) / (1 + A * 2 / B * 2) * np.sqrt(1 - (B * 2 + A * 2) / (4 * r2 ** 2))
    Ep = r2 * np.cos(C) / (1 + A * 2 / B * 2) * np.sqrt(1 - (B * 2 + A * 2) / (4 * r2 ** 2))
    J11 = -r1 * s1 / 2 - 2 * np.cos(C) * (A * r1 * c1 - B * r1 * s1) / D_val - E / B ** 2 * (B * r1 * c1 + A * r1 * s1)
    J12 = -r1 * s4 / 2 + 2 * np.cos(C) * (A * r1 * c4 + B * r1 * s4) / D_val - E / B ** 2 * (-B * r1 * c4 + A * r1 * s4)
    J21 = -r1 * c1 / 2 - 2 * np.sin(C) * (A * r1 * c1 - B * r1 * s1) / D_val + Ep / B ** 2 * (B * r1 * c1 + A * r1 * s1)
    J22 = -r1 * c4 / 2 + 2 * np.sin(C) * (A * r1 * c4 + B * r1 * s4) / D_val + Ep * B ** 2 * (-B * r1 * c4 + A * r1 * s4)
    return J11, J12, J21, J22

def dynamique(theta1, theta4, tau1, tau2):
    J_T = np.transpose(jacobien(theta1, theta4))
    F = np.linalg.solve(J_T, np.array([tau1, tau2]))
    return F[0], F[1]
def calcul_dynamique(row):
    t1, t4 = row['Theta1_Encodeur'], row['Theta4_Encodeur']
    tau1, tau4 = row['Tau1'], row['Tau4']
    
    a11, a12, a21, a22= jacobien(t1, t4)
    
    A = np.array([[a11, a12], [a21, a22]])
    
    # det(A) est critique pour la stabilité (Type II)
    detA = np.linalg.det(A)
    
    if np.abs(detA) < 1e-2: # Seuil d'instabilité
        return 0.0, 0.0, detA
    
    # J = A^-1 * B  => J^T = B^T * (A^-1)^T
    J_T = np.transpose(A)
    
    try:
        F = np.linalg.solve(J_T, np.array([tau1, tau4]))
        return F[0], F[1], detA
    except np.linalg.LinAlgError:
        return 0.0, 0.0, detA
CHEMIN_CSV = "1rectangle_force_externe.csv"

# 1. Chargement des données
df = pd.read_csv(CHEMIN_CSV)

ids = df['Essai_ID'].values

# Découpage de tous les vecteurs
tours_XR = df['X_Encodeur'].values
tours_YR = df['Y_Encodeur'].values 
tours_FX = df['Force1'].values
tours_FY = df['Force4'].values 
tours_detJ = df['Jacobien'].values 
tours_torque1 = df['Tau1'].values
tours_torque4 = df['Tau4'].values
tours_theta1 = df['Theta1_Encodeur'].values
tours_theta4 = df['Theta4_Encodeur'].values 

#df['Tau1'] = df['Tau1']/np.abs(df['Tau1']) * np.abs((df['Courant1']-0.62)/0.92)
results = df.apply(calcul_dynamique, axis=1)
df['Force1'], df['Force4'], df['Jacobien'] = zip(*results)
new_file_name = CHEMIN_CSV.split('.')[0] + '_MODIFIÉ.csv'
# for i in range(len(ids)):
#     if np.abs(tours_detJ[i]) < 0.05:
#         df.at[i, 'Force1'] = 0
#         df.at[i, 'Force4'] = 0
#         df.at[i, 'Tau1'] = 0
#         df.at[i, 'Tau4'] = 0
df.to_csv(new_file_name, index= False)
print("✅ fichier csv modifié")