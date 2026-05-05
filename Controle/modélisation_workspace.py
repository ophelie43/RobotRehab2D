import numpy as np 
from math import atan, sqrt, degrees, radians
import matplotlib.pyplot as plt
import matplotlib.patches as patches

l1, l2, l5 = 0.265, 0.3, 0.08

def cinematique(theta1, theta4): 
    E = 2 * l2 * (l5 + l1 * (np.cos(theta4) - np.cos(theta1)))
    F = 2 * l1 * l2 * (np.sin(theta4) - np.sin(theta1))
    G = l5 ** 2 + 2 * l1 ** 2 + 2 * l5 * l1 * np.cos(theta4) - 2 * l5 * l1 * np.cos(theta1) - 2 * l1 ** 2 * np.cos(theta4 - theta1)
    xc = l5 + l1 * np.cos(theta4) + l2 * np.cos(2*np.atan((-F-(E**2+F**2-G**2)**0.5)/(G-E)))
    yc = l1 * np.sin(theta4) + l2 * np.sin(2*np.atan((-F-(E**2+F**2-G**2)**0.5)/(G-E)))
    return xc, yc

def cinematique_inverse(xc,yc): 
    E1 = -2 * l1 * xc
    F1 = - 2 * l1 * yc
    G1 = l1 ** 2 - l2 ** 2 + xc ** 2 + yc ** 2 
    E4 = 2 * l1 * (-xc + l5)
    F4 = -2 * l1 * yc
    G4 = l5 ** 2 + l1 ** 2 - l2 ** 2 + xc ** 2 + yc ** 2 - 2 * l5 * xc
    theta1 = 2 * np.atan((-F1+(E1**2+F1**2-G1**2)**0.5) / (G1-E1))
    theta4 = 2 * np.atan((-F4-(E4**2+F4**2-G4**2)**0.5) / (G4-E4))
    return degrees(theta1), degrees(theta4)

def jacobien(theta1, theta4):
    D = (l1 + l2 + l5 / 2) / 3
    r1, r2 = l1 / D, l2 / D
    s1, s4 = np.sin(theta1), np.sin(theta4)
    c1, c4 = np.cos(theta1), np.cos(theta4)
    A = r1 * s1 - r1 * s4
    B = 2 * (l5 / 2) / D + r1 * c1 + r1 * c4
    C = np.pi / 2 + np.arctan(A / B)
    D_val = 8 * r2 * np.sqrt(1 - (B * 2 + A * 2) / (4 * r2 ** 2))
    E = r2 * np.sin(C) / (1 + A * 2 / B * 2) * np.sqrt(1 - (B * 2 + A * 2) / (4 * r2 ** 2))
    Ep = r2 * np.cos(C) / (1 + A * 2 / B * 2) * np.sqrt(1 - (B * 2 + A * 2) / (4 * r2 ** 2))
    J11 = -r1 * s1 / 2 - 2 * np.cos(C) * (A * r1 * c1 - B * r1 * s1) / D_val - E / B ** 2 * (B * r1 * c1 + A * r1 * s1)
    J12 = -r1 * s4 / 2 + 2 * np.cos(C) * (A * r1 * c4 + B * r1 * s4) / D_val - E / B ** 2 * (-B * r1 * c4 + A * r1 * s4)
    J21 = -r1 * c1 / 2 - 2 * np.sin(C) * (A * r1 * c1 - B * r1 * s1) / D_val + Ep / B ** 2 * (B * r1 * c1 + A * r1 * s1)
    J22 = -r1 * c4 / 2 + 2 * np.sin(C) * (A * r1 * c4 + B * r1 * s4) / D_val + Ep * B ** 2 * (-B * r1 * c4 + A * r1 * s4)
    return np.array([[J11, J12], [J21, J22]])

def dynamique(theta1, theta4, tau1, tau2):
    J_T = np.transpose(jacobien(theta1, theta4))
    F = np.linalg.solve(J_T, np.array([tau1, tau2]))
    return F[0], F[1]

def compute_theta1(xc, yc):
    E1, F1, G1 = -2 * l1 * xc, -2 * l1 * yc, l1**2 - l2**2 + xc**2 + yc**2
    root = sqrt(max(0.0, E1**2 + F1**2 - G1**2))
    return 2 * atan((-F1 + root) / (G1 - E1))

def compute_theta4(xc, yc):
    E4, F4, G4 = 2 * l1 * (-xc + l5), -2 * l1 * yc, l5**2 + l1**2 - l2**2 + xc**2 + yc**2 - 2 * l5 * xc
    root = sqrt(max(0.0, E4**2 + F4**2 - G4**2))
    return 2 * atan((-F4 - root) / (G4 - E4))

print(cinematique_inverse(0.1, 0.3))

print(cinematique(radians(132.9666468467352), radians(22.494068109637446)))

# 1. Génération des angles de 0 à 180 degrés (0 à pi radians)
# On augmente la résolution (150 points) pour un graphique plus dense
angles_t1 = np.linspace(0, np.pi, 150)
angles_t4 = np.linspace(0, np.pi, 150)

x_workspace = []
y_workspace = []

def cinematique_securisee(theta1, theta4): 
    E = 2 * l2 * (l5 + l1 * (np.cos(theta4) - np.cos(theta1)))
    F = 2 * l1 * l2 * (np.sin(theta4) - np.sin(theta1))
    G = l5 ** 2 + 2 * l1 ** 2 + 2 * l5 * l1 * np.cos(theta4) - 2 * l5 * l1 * np.cos(theta1) - 2 * l1 ** 2 * np.cos(theta4 - theta1)
    
    # Vérification de l'assemblage physique (la racine doit être positive)
    val_sous_racine = E**2 + F**2 - G**2
    if val_sous_racine < 0:
        return None, None  # Position physiquement inatteignable
        
    # Calcul avec votre configuration de signes
    angle_intermediaire = 2 * np.atan((-F - (val_sous_racine)**0.5) / (G - E))
    xc = l5 + l1 * np.cos(theta4) + l2 * np.cos(angle_intermediaire)
    yc = l1 * np.sin(theta4) + l2 * np.sin(angle_intermediaire)
    
    return xc, yc

# 2. Balayage de toutes les combinaisons possibles
print("Calcul de l'espace de travail en cours...")
for t1 in angles_t1:
    for t4 in angles_t4:
        x, y = cinematique_securisee(t1, t4)
        if x is not None and y is not None:
            # On ne garde que les points avec un Y positif (le robot travaille vers l'avant)
            if y >= 0:
                x_workspace.append(x)
                y_workspace.append(y)

# 3. Visualisation des données
plt.figure(figsize=(10, 8))
plt.scatter(x_workspace, y_workspace, s=2, c='#007acc', alpha=0.6, label='Zone atteignable')

# Ajout des repères des moteurs
plt.plot(0, 0, 'ko', markersize=8, label='Moteur 1 (Origine)')
plt.plot(l5, 0, 'ro', markersize=8, label='Moteur 2 (lc)')

plt.title('Espace de travail du mécanisme à 5 barres', fontsize=14)
plt.xlabel('Axe X (m)', fontsize=12)
plt.ylabel('Axe Y (m)', fontsize=12)

# L'axe 'equal' est crucial en mécanique pour ne pas déformer le visuel
plt.axis('equal')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(loc='lower right')

plt.show()

# --- Définition de la Zone Utile (Sécuritaire) ---
x_min, x_max = -0.20, 0.20
y_min, y_max = 0.23, 0.40

# --- Balayage ---
angles_t1 = np.linspace(0, np.pi, 150)
angles_t4 = np.linspace(0, np.pi, 150)

x_full = []
y_full = []
x_utile = []
y_utile = []

print("Calcul des espaces de travail...")
for t1 in angles_t1:
    for t4 in angles_t4:
        x, y = cinematique_securisee(t1, t4)
        if x is not None and y is not None and y >= 0:
            x_full.append(x)
            y_full.append(y)
            # Tri des points qui tombent dans la zone utile
            if (x_min <= x <= x_max) and (y_min <= y <= y_max):
                x_utile.append(x)
                y_utile.append(y)

# --- Affichage ---
fig, ax = plt.subplots(figsize=(10, 8))

# 1. Espace atteignable global (Gris clair)
ax.scatter(x_full, y_full, s=2, c='#cccccc', alpha=0.5, label='Espace atteignable complet (instable aux bords)')

# 2. Espace utile (Vert)
ax.scatter(x_utile, y_utile, s=2, c='#2ca02c', alpha=0.8, label='Zone de travail Utile (Fluide et stable)')

# 3. Dessin de la boîte de sécurité (Bounding box)
rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, 
                         linewidth=2, edgecolor='red', facecolor='none', label='Limites Logicielles (np.clip)')
ax.add_patch(rect)

# Repères des moteurs
ax.plot(0, 0, 'ko', markersize=8, label='Moteur 1 (S1)')
ax.plot(l5, 0, 'ro', markersize=8, label='Moteur 2 (S2)')

ax.set_title('Analyse de l\'Espace de Travail du Robot', fontsize=14)
ax.set_xlabel('Axe X (m)', fontsize=12)
ax.set_ylabel('Axe Y (m)', fontsize=12)
ax.axis('equal')
ax.grid(True, linestyle='--', alpha=0.7)
ax.legend(loc='lower right')

plt.show()