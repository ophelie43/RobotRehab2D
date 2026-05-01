import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

la, lb, lc = 0.27, 0.315, 0.08
def filter_force(force: np.ndarray, f_cutoff, dt):

    alpha = 1-np.exp(-dt * 2 *np.pi * f_cutoff)
    Fm_n = force[-2] * (1-alpha) + force[-1] * alpha
    return Fm_n
def jacobien(theta1, theta4):
    D = (la + lb + lc / 2) / 3
    r1, r2 = la / D, lb / D
    s1, s4 = np.sin(theta1), np.sin(theta4)
    c1, c4 = np.cos(theta1), np.cos(theta4)
    A = r1 * s1 - r1 * s4
    B = 2 * (lc / 2) / D + r1 * c1 + r1 * c4
    C = np.pi / 2 + np.arctan(A / B)
    D_val = 8 * r2 * np.sqrt(1 - (B ** 2 + A ** 2) / (4 * r2 ** 2))
    E = r2 * np.sin(C) / (1 + A ** 2 / B ** 2) * np.sqrt(1 - (B ** 2 + A ** 2) / (4 * r2 ** 2))
    Ep = r2 * np.cos(C) / (1 + (A ** 2 / B ** 2)) * np.sqrt(1 - (B ** 2 + A ** 2) / (4 * r2 ** 2))
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
    
    if np.abs(detA) < 5e-2: # Seuil d'instabilité
        return 0.0, 0.0, 0.1
    
    # J = A^-1 * B  => J^T = B^T * (A^-1)^T
    J_T = np.transpose(A)
    
    try:
        F = np.linalg.solve(J_T, np.array([tau1, tau4]))
        return F[0], F[1], detA
    except np.linalg.LinAlgError:
        return 0.0, 0.0, 0.1


def appliquer_filtre_df(df, f_cutoff, dt):
    # Calcul du coefficient alpha
    alpha = 1 - np.exp(-dt * 2 * np.pi * f_cutoff)
    
    # Initialisation des colonnes filtrées
    df['Force1_Filtree'] = df['Force1'].copy()
    df['Force4_Filtree'] = df['Force4'].copy()
    
    # Application récursive (Filtre IIR)
    for i in range(1, len(df)):
        df.at[i, 'Force1_Filtree'] = df.at[i-1, 'Force1_Filtree'] * (1 - alpha) + df.at[i, 'Force1'] * alpha
        df.at[i, 'Force4_Filtree'] = df.at[i-1, 'Force4_Filtree'] * (1 - alpha) + df.at[i, 'Force4'] * alpha
    
    return df

# --- Dans ton code principal ---

CHEMIN_CSV = "PID/CSV/10cercles_barre.csv"

# 1. Chargement des données
df = pd.read_csv(CHEMIN_CSV)
results = df.apply(calcul_dynamique, axis=1)
df['Force1'], df['Force4'], df['Jacobien'] = zip(*result