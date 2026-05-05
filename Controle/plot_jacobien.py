import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# --- CONFIGURATION ---
CHEMIN_CSV = "rectangle.csv"
def generer_rapport_visuel(chemin_fichier):
    if not os.path.exists(chemin_fichier):
        print(f"Erreur : Le fichier {chemin_fichier} est introuvable.")
        return

    # 1. Chargement des données
    df = pd.read_csv(chemin_fichier)
    
    # Palette de couleurs
    couleurs_palette = ['#FFD700', '#FF0000', '#00FF00', '#0000FF', '#FF00FF', 
                        '#00FFFF', '#FFA500', '#800080', '#A52A2A', '#008080']

    # --- PRÉPARATION DES DONNÉES (Une seule fois) ---
    # On utilise l'ID d'essai pour séparer les tours (ou une remise à zéro d'un compteur)
    ids = df['Essai_ID'].values
    points_de_coupure = np.where(np.diff(ids) < 0)[0] + 1

    # Découpage de tous les vecteurs
    tours_XR = np.split(df['X_Encodeur'].values, points_de_coupure)
    tours_YR = np.split(df['Y_Encodeur'].values, points_de_coupure)
    tours_FX = np.split(df['Force1'].values, points_de_coupure)
    tours_FY = np.split(df['Force4'].values, points_de_coupure)
    # tours_FX_Filtree = np.split(df['Force1_Filtree'].values, points_de_coupure)
    # tours_FY_Filtree = np.split(df['Force4_Filtree'].values, points_de_coupure)
    tours_detJ = np.split(df['Jacobien'].values, points_de_coupure)
    tours_torque1 = np.split(df['Tau1'].values, points_de_coupure)
    tours_torque4 = np.split(df['Tau4'].values, points_de_coupure)
    tours_theta1 = np.split(df['Theta1_Encodeur'].values, points_de_coupure)
    tours_theta4 = np.split(df['Theta4_Encodeur'].values, points_de_coupure)
    # 2. Création de la figure
    fig = plt.figure(figsize=(16, 12))
    fig.canvas.manager.set_window_title(f"Analyse Capstone - {os.path.basename(chemin_fichier)}")

    ax1 = plt.subplot2grid((4, 2), (0, 0), rowspan=2) 
    ax_fx = plt.subplot2grid((4, 2), (0, 1))           
    ax_fy = plt.subplot2grid((4, 2), (1, 1))           
    ax_detJ = plt.subplot2grid((4, 2), (2, 1))
    ax_torque = plt.subplot2grid((4,2), (2,0))
    ax_theta = plt.subplot2grid((4,2), (3,0))

    # Config axes
    for ax, titre in [(ax_fx, "Répétabilité Force X (Force1)"), (ax_fy, "Répétabilité Force Y (Force4)")]:
        ax.set_title(titre, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='black', linewidth=0.8, linestyle='--')

    ax_detJ.set_title("Déterminant de la Jacobienne (Zones critiques en rouge)", fontweight='bold')
    ax_detJ.grid(True, alpha=0.3)
    ax_detJ.axhline(0, color='black', linewidth=0.8, linestyle='--')

    ax_torque.set_title("Force X Filtree", fontweight='bold')
    ax_torque.grid(True, alpha=0.3)
    ax_torque.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax_torque.set_ylabel("Force [N]")

    ax_theta.set_title("Force Y Filtre", fontweight = 'bold')
    ax_theta.grid(True, alpha = 0.3)
    ax_theta.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax_theta.set_ylabel("Force [N]")

    # 3. Tracé de la consigne (Background)
    ax1.plot(df['X_Voulu'], df['Y_Voulu'], 'b--', label='Cible (Voulue)', linewidth=1, alpha=0.3)

    # 4. Boucle de tracé par TOUR
    for i in range(len(tours_XR)):
        couleur = couleurs_palette[i % len(couleurs_palette)]
        
        xr = tours_XR[i]
        yr = tours_YR[i]
        fx = tours_FX[i]
        fy = tours_FY[i]
        # fx_f = tours_FX_Filtree[i]
        # fy_f = tours_FY_Filtree[i]
        dj = tours_detJ[i]
        t1 = tours_torque1[i]
        t4 = tours_torque4[i]
        theta1 = tours_theta1[i]
        theta4 = tours_theta4[i]
        indices_local = np.arange(len(xr))
        
        # A. Trajectoire et Vecteurs Quiver
        pas = max(1, len(xr) // 15) 
        ax1.plot(xr, yr, color=couleur, alpha=0.4, linewidth=1)
        ax1.quiver(xr[::pas], yr[::pas], fx[::pas], fy[::pas], 
                   color=couleur, alpha=0.7, width=0.004, scale=50)

        # B. Forces
        ax_fx.plot(indices_local, fx, color=couleur, alpha=0.8, label=f'Tour {i+1}')
        ax_fy.plot(indices_local, fy, color=couleur, alpha=0.8)

        # C. Jacobienne et Percentile 2% (Zones de singularité)
        dj_clip = np.clip(dj, -0.2, 0.2)
        ax_detJ.plot(indices_local, dj_clip, color="black", alpha=0.2)

        ax_torque.plot(indices_local,t1, color = couleur, alpha = 0.8)
        # Filtrage des 2% plus proches de zéro (en valeur absolue)
        seuil = np.percentile(np.abs(dj), 2)
        masque = np.abs(dj) <= seuil
        
        # On utilise NaN pour ne pas tracer les points hors seuil
        vecteur_rouge = np.full_like(dj, np.nan)
        vecteur_rouge[masque] = dj[masque]
        
        ax_detJ.scatter(indices_local, dj, color=couleur, s=10, zorder=5)

        ax_theta.plot(indices_local, t4, color = couleur, alpha = 0.8)
    # 5. Finalisation
    ax1.set_title("Trajectoire et Vecteurs de Force", fontweight='bold')
    ax1.axis('equal')
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='lower left', fontsize='small')

    ax_fx.legend(loc='upper right', fontsize='7', ncol=2)
    
    # Ajout d'une entrée de légende manuelle pour le rouge
    from matplotlib.lines import Line2D
    securite_legende = Line2D([0], [0], marker='o', color='w', label='Proche Singularité (2%)',
                              markerfacecolor='red', markersize=5)
    handles, labels = ax_detJ.get_legend_handles_labels()
    handles.append(securite_legende)
    ax_detJ.legend(handles=handles, loc='upper right', fontsize='7')

    plt.tight_layout()
    print(f"Jarvis : Analyse terminée. {len(tours_XR)} tours détectés.")
    plt.show()

if __name__ == "__main__":
    generer_rapport_visuel(CHEMIN_CSV)