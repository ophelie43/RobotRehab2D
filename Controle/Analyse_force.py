import serial
import time
import threading
import numpy as np
import csv
import os
from math import atan, sqrt, degrees, radians, cos, sin, atan2
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- 1. Constantes et Modèle Cinématique ---
la, lb, lc = 0.27, 0.315, 0.08
NOM_FICHIER_CSV = "rectangle.csv"

def cinematique_directe(theta1, theta4):
    t1, t4 = radians(theta1), radians(theta4)
    try:
        E = 2 * lb * (lc + la * (cos(t4) - cos(t1)))
        F = 2 * la * lb * (sin(t4) - sin(t1))
        G = lc**2 + 2 * la**2 + 2 * lc * la * cos(t4) - 2 * lc * la * cos(t1) - 2 * la**2 * cos(t4 - t1)
        arg_racine = max(0, E**2 + F**2 - G**2)
        angle_interne = 2 * atan2((-F - sqrt(arg_racine)), (G - E))
        xc = lc + la * cos(t4) + lb * cos(angle_interne)
        yc = la * sin(t4) + lb * sin(angle_interne)
        return xc, yc
    except:
        return None, None

def cinematique_inverse(xc, yc):
    E1, F1 = -2 * la * xc, -2 * la * yc
    G1 = la**2 - lb**2 + xc**2 + yc**2 
    E4, F4 = 2 * la * (-xc + lc), -2 * la * yc
    G4 = lc**2 + la**2 - lb**2 + xc**2 + yc**2 - 2 * lc * xc 
    theta1 = 2 * np.atan((-F1 + np.sqrt(max(0, E1**2 + F1**2 - G1**2))) / (G1 - E1)) 
    theta4 = 2 * np.atan((-F4 - np.sqrt(max(0, E4**2 + F4**2 - G4**2))) / (G4 - E4)) 
    return degrees(theta1), degrees(theta4)

def jacobien(theta1, theta4):
    # Sécurisation mathématique totale des divisions par zéro
    t1, t4 = np.radians(theta1), np.radians(theta4)
    D = (la + lb + lc / 2) / 3
    r1, r2 = la / D, lb / D
    s1, s4 = np.sin(t1), np.sin(t4)
    c1, c4 = np.cos(t1), np.cos(t4)
    A = r1 * s1 - r1 * s4
    B = 2 * (lc / 2) / D + r1 * c1 + r1 * c4
    
    if abs(B) < 1e-6: B = 1e-6 # Anti-crash division
    C = np.pi / 2 + np.arctan(A / B)
    
    arg_racine = abs(1 - (B ** 2 + A ** 2) / (4 * r2 ** 2))
    D_val = 8 * r2 * np.sqrt(arg_racine)
    if abs(D_val) < 1e-6: D_val = 1e-6 # Anti-crash division
    
    E = r2 * np.sin(C) / (1 + A ** 2 / B ** 2) * np.sqrt(arg_racine)
    Ep = r2 * np.cos(C) / (1 + A ** 2 / B ** 2) * np.sqrt(arg_racine)
    
    J11 = -r1 * s1 / 2 - 2 * np.cos(C) * (A * r1 * c1 - B * r1 * s1) / D_val - E / B ** 2 * (B * r1 * c1 + A * r1 * s1)
    J12 = -r1 * s4 / 2 + 2 * np.cos(C) * (A * r1 * c4 + B * r1 * s4) / D_val - E / B ** 2 * (-B * r1 * c4 + A * r1 * s4)
    J21 = -r1 * c1 / 2 - 2 * np.sin(C) * (A * r1 * c1 - B * r1 * s1) / D_val + Ep / B ** 2 * (B * r1 * c1 + A * r1 * s1)
    J22 = -r1 * c4 / 2 + 2 * np.sin(C) * (A * r1 * c4 + B * r1 * s4) / D_val + Ep * B ** 2 * (-B * r1 * c4 + A * r1 * s4)
    return np.array([[J11, J12], [J21, J22]])

def dynamique(theta1, theta4, tau1, tau2):
    J_T = np.transpose(jacobien(theta1, theta4))
    # Utilisation de pinv (pseudo-inverse) au lieu de solve : NE PLANTE JAMAIS !
    F = np.linalg.pinv(J_T).dot(np.array([tau1, tau2]))
    return F[0], F[1]

def current_to_torque(current, theta): 
    direction = 1.0 
    for i in range(1, len(theta)):
        index = -1 - i
        delta_theta = theta[-1] - theta[index] # Correction du signe mathématique
        if abs(delta_theta) > 0.01:
            direction = delta_theta / abs(delta_theta)
            break
    return abs((current - 0.218) / 2.011) * direction

# --- 2. CLASSE PID AVANCÉ ---
class PIDController:
    def __init__(self, Kp, Ki, Kd, alpha=0.05, limite_sortie=50.0):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.alpha = alpha  
        self.limite_sortie = limite_sortie 
        self.en_1 = 0.0
        self.In_1 = 0.0
        self.Dn_1 = 0.0
        self.dernier_temps = time.time()

    def calculer(self, consigne, valeur_mesuree):
        maintenant = time.time()
        dt = maintenant - self.dernier_temps
        if dt <= 0.0: dt = 0.001
        en = consigne - valeur_mesuree
        Pn = self.Kp * en
        C1 = (2.0 * self.alpha - dt) / (2.0 * self.alpha + dt)
        C2 = (2.0 * self.Kd) / (2.0 * self.alpha + dt)
        Dn = (C1 * self.Dn_1) + (C2 * (en - self.en_1))
        Un_temp = Pn + self.In_1 + Dn
        In = self.In_1
        if abs(Un_temp) < self.limite_sortie:
            In = self.In_1 + (dt / 2.0) * self.Ki * (en + self.en_1)
        Un = Pn + In + Dn
        sortie = max(min(Un, self.limite_sortie), -self.limite_sortie)
        self.en_1 = en
        self.In_1 = In
        self.Dn_1 = Dn
        self.dernier_temps = maintenant
        return sortie

    def reset(self):
        self.en_1 = 0.0
        self.In_1 = 0.0
        self.Dn_1 = 0.0
        self.dernier_temps = time.time()

pid_moteur1 = PIDController(Kp=0.8, Ki=2.0, Kd=0.05, alpha=0.03, limite_sortie=50.0)
pid_moteur2 = PIDController(Kp=0.8, Ki=2.0, Kd=0.05, alpha=0.03, limite_sortie=50.0)

# --- 3. Fonctions CSV ---
def initialiser_csv():
    if not os.path.exists(NOM_FICHIER_CSV):
        with open(NOM_FICHIER_CSV, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Horodatage", "Essai_ID", "X_Voulu", "Y_Voulu", "Theta1_Cible", "Theta4_Cible", "Correction_M1", "Correction_M2", "X_Encodeur", "Y_Encodeur", "Theta1_Encodeur", "Theta4_Encodeur", "Courant1", "Courant4", "Tau1", "Tau4", "Force1", "Force4"])

def sauvegarder_donnees(data_list):
    with open(NOM_FICHIER_CSV, mode='a', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data_list)

# --- 4. Affichage et Analyse (Matplotlib) ---
def tracer_analyse(donnees):
    if not donnees:
        print("Jarvis : Aucune donnée à tracer.")
        return

    # 1. Extraction des IDs d'essais uniques
    tous_les_ids = [ligne[1] for ligne in donnees]
    ids_uniques = []
    for id_essai in tous_les_ids:
        if id_essai not in ids_uniques:
            ids_uniques.append(id_essai)
            
    # Palette de couleurs étendue pour 10 essais
    couleurs_palette = ['#FFD700', '#FF0000', '#00FF00', '#0000FF', '#FF00FF', 
                        '#00FFFF', '#FFA500', '#800080', '#A52A2A', '#008080']

    # Extraction des données globales
    x_voulu = [ligne[2] for ligne in donnees]
    y_voulu = [ligne[3] for ligne in donnees]
    x_reel = [ligne[8] for ligne in donnees]
    y_reel = [ligne[9] for ligne in donnees]
    force_x_globale = [ligne[16] for ligne in donnees]
    force_y_globale = [ligne[17] for ligne in donnees]
    
    # Création d'une figure complexe : 2 lignes, 2 colonnes
    fig = plt.figure(figsize=(16, 10))
    fig.canvas.manager.set_window_title("Analyse Comparative Détaillée")

    # Ax1 : La trajectoire (occupe les deux lignes de la colonne 0)
    ax1 = plt.subplot2grid((2, 2), (0, 0), rowspan=2)
    
    # Ax2 : Force X (Ligne 0, Colonne 1)
    ax_fx = plt.subplot2grid((2, 2), (0, 1))
    
    # Ax3 : Force Y (Ligne 1, Colonne 1)
    ax_fy = plt.subplot2grid((2, 2), (1, 1))

    # Configuration des axes de force
    for ax, titre in [(ax_fx, "Répétabilité Force X"), (ax_fy, "Répétabilité Force Y")]:
        ax.set_title(titre)
        ax.set_xlabel("Itérations au sein du tour")
        ax.set_ylabel("Force (N)")
        ax.axhline(0, color='black', linewidth=1, linestyle='--')
        ax.grid(True)

    # --- Tracé global ---
    ax1.plot(x_voulu, y_voulu, 'b--', label='Cible (Voulue)', linewidth=2, alpha=0.4)
    ax1.plot(x_reel, y_reel, 'k-', label='Réel (Encodeurs)', linewidth=2, alpha=0.1)

    # --- Boucle par tour ---
    for i, id_actuel in enumerate(ids_uniques):
        indices = [idx for idx, val in enumerate(tous_les_ids) if val == id_actuel]
        xr_sub = [x_reel[idx] for idx in indices]
        yr_sub = [y_reel[idx] for idx in indices]
        fx_sub = [force_x_globale[idx] for idx in indices]
        fy_sub = [force_y_globale[idx] for idx in indices]
        
        couleur = couleurs_palette[i % len(couleurs_palette)]
        temps_local = range(len(indices))
        
        # 1. Vecteurs sur la trajectoire
        pas = max(1, len(xr_sub) // 12) 
        ax1.quiver(xr_sub[::pas], yr_sub[::pas], fx_sub[::pas], fy_sub[::pas], 
                   color=couleur, alpha=0.7, width=0.004)

        # 2. Courbe Force X
        ax_fx.plot(temps_local, fx_sub, color=couleur, label=f'Tour {i+1}', alpha=0.8)
        
        # 3. Courbe Force Y
        ax_fy.plot(temps_local, fy_sub, color=couleur, label=f'Tour {i+1}', alpha=0.8)

    # Finalisation Ax1
    ax1.set_title("Trajectoire et Vecteurs de Force")
    ax1.set_xlabel("X (m)")
    ax1.set_ylabel("Y (m)")
    ax1.axis('equal')
    ax1.grid(True)

    # Légendes
    ax_fx.legend(loc='upper right', fontsize='xx-small', ncol=2)
    ax_fy.legend(loc='upper right', fontsize='xx-small', ncol=2)
    
    plt.tight_layout()
    plt.show()

# --- 5. Configuration des Ports Séries ---
try:
    carte_servo1 = serial.Serial('COM6', 115200, timeout=0.05)
    carte_servo2 = serial.Serial('COM5', 115200, timeout=0.05)
    print("Jarvis : Connexion établie.")
except serial.SerialException as e:
    print(f"Jarvis : Erreur COM - {e}"); exit()

time.sleep(2) 

# --- 6. Variables Globales & Suivi ---
offsets_encodeurs = {"COM6 - S1": None, "COM5 - S2": None}
angles_consigne_init = {"COM6 - S1": 135.0 - 7.0, "COM5 - S2": 45.0 + 9.0}
angles_reels_calib = {"COM6 - S1": 135.0, "COM5 - S2": 45.0}
angles_actuels = {"COM6 - S1": [135.0], "COM5 - S2": [45.0]}
courants_actuels = {"COM6 - S1": 0.0, "COM5 - S2": 0.0}
commande_en_attente = None
mode_silencieux = True

def lire_et_afficher(port, nom_carte):
    global offsets_encodeurs, angles_actuels, courants_actuels, mode_silencieux
    
    # LE SECRET EST ICI : "while" au lieu de "if".
    # On vide tout le bouchon de données pour avoir la position en TEMPS RÉEL
    while port.in_waiting > 0:
        ligne = port.readline().decode('utf-8', errors='ignore').strip()
        if ligne and "ENCODEUR:" in ligne and "| COURANT:" in ligne:
            try:
                parties = ligne.split("|")
                angle = float(parties[0].split(":")[1])
                courant = float(parties[1].split(":")[1])
                
                if not mode_silencieux:
                    print(f"{nom_carte} >> ENCODEUR: {angle} | COURANT: {courant}")
                    
                if offsets_encodeurs[nom_carte] is None:
                    offsets_encodeurs[nom_carte] = angle

                # L'angle ajouté à la fin de la boucle sera le plus récent physiquement !
                angles_actuels[nom_carte].append((angle - offsets_encodeurs[nom_carte]) + angles_reels_calib[nom_carte])
                
                # On ne garde que les 100 dernières valeurs en mémoire pour ne pas faire exploser la RAM
                if len(angles_actuels[nom_carte]) > 100:
                    angles_actuels[nom_carte].pop(0)
                    
                courants_actuels[nom_carte] = courant
            except: pass
# --- 7. Fonctions de Mouvement ---
def executer_trajectoire():
    essai_id_base = time.strftime("CARRE_%H%M%S")
    # Sommets décalés de 4cm vers X+ : (x + 0.04, y)
    sommets = [(0.04, 0.4), (0.14, 0.4), (0.14, 0.3), (-0.06, 0.3), (-0.06, 0.4), (0.04, 0.4)]
    nb_pas = 60 
    donnees_a_sauver = []
    
    print(f"\nJarvis : Démarrage de la trajectoire CARRE décalée (Essai {essai_id_base})...")
    
    # Homing sur le premier point décalé
    point_depart = sommets[0]
    th1_init, th4_init = cinematique_inverse(point_depart[0], point_depart[1])
    carte_servo1.write(f"S:{th1_init - 7.0}\n".encode())
    carte_servo2.write(f"S:{th4_init + 9.0}\n".encode())
    
    print("Jarvis : Mise en position initiale. Attente de stabilisation...")
    time.sleep(2) 
    
    pid_moteur1.reset()
    pid_moteur2.reset()
    
    nb_essais = 10 
    for n in range(nb_essais):
        essai_tour_id = f"{essai_id_base}_T{n+1}"
        print(f"Jarvis : Carré {n+1}/{nb_essais} en cours...")

        for i in range(len(sommets) - 1):
            p_start = sommets[i]
            p_end = sommets[i+1]
            
            for pas in range(nb_pas):
                lire_et_afficher(carte_servo1, "COM6 - S1")
                lire_et_afficher(carte_servo2, "COM5 - S2")
                
                ratio = pas / nb_pas
                px = p_start[0] + (p_end[0] - p_start[0]) * ratio
                py = p_start[1] + (p_end[1] - p_start[1]) * ratio
                
                try:
                    th1_cible, th4_cible = cinematique_inverse(px, py)
                    ang1_real = angles_actuels["COM6 - S1"][-1]
                    ang4_real = angles_actuels["COM5 - S2"][-1]
                    current1 = courants_actuels["COM6 - S1"]
                    current4 = courants_actuels["COM5 - S2"]
                    tau1 = current_to_torque(current1, angles_actuels["COM6 - S1"])
                    tau4 = current_to_torque(current4, angles_actuels["COM5 - S2"])
                    force1, force4 = dynamique(ang1_real, ang4_real, tau1, tau4)

                    corr1 = pid_moteur1.calculer(th1_cible, ang1_real)
                    corr2 = pid_moteur2.calculer(th4_cible, ang4_real)
                    
                    carte_servo1.write(f"S:{th1_cible + corr1 - 7.0}\n".encode())
                    carte_servo2.write(f"S:{th4_cible + corr2 + 9.0}\n".encode())
                    
                    time.sleep(0.08)
                    
                    rx, ry = cinematique_directe(ang1_real, ang4_real)
                    if rx is not None:
                        donnees_a_sauver.append([
                            time.strftime("%H:%M:%S"), essai_tour_id, px, py, 
                            round(th1_cible, 2), round(th4_cible, 2),
                            round(corr1, 2), round(corr2, 2),
                            round(rx, 4), round(ry, 4), 
                            round(ang1_real, 2), round(ang4_real, 2), 
                            round(current1, 2), round(current4, 2),
                            round(tau1, 3), round(tau4, 3), 
                            round(force1, 3), round(force4, 3)
                        ])
                except Exception as e:
                    continue 

    sauvegarder_donnees(donnees_a_sauver)
    print(f"Jarvis : {nb_essais} carrés terminés.")
    tracer_analyse(donnees_a_sauver)


def executer_cercle():
    essai_id_base = time.strftime("CERCLE_%H%M%S")
    # Centre décalé de 4cm : 0.0 + 0.04
    centre_x, centre_y, rayon = 0.04, 0.35, 0.1 
    nb_points = 200 
    donnees_a_sauver = []
    
    print(f"\nJarvis : Génération CERCLE décalé (Essai {essai_id_base})...")
    
    # Positionnement initial décalé
    th1_init, th4_init = cinematique_inverse(centre_x + rayon, centre_y)
    carte_servo1.write(f"S:{th1_init - 7.0}\n".encode())
    carte_servo2.write(f"S:{th4_init + 9.0}\n".encode())
    time.sleep(2)
    
    pid_moteur1.reset()
    pid_moteur2.reset()

    nb_tours = 1
    for n in range(nb_tours):
        essai_tour_id = f"{essai_id_base}_T{n+1}"
        print(f"Jarvis : Tour {n+1}/{nb_tours} en cours...")
        
        for i in range(nb_points + 1):
            lire_et_afficher(carte_servo1, "COM6 - S1")
            lire_et_afficher(carte_servo2, "COM5 - S2")
            
            angle_cercle = (2 * np.pi * i) / nb_points
            # px utilise maintenant le centre_x décalé à 0.04
            px = centre_x + rayon * np.cos(angle_cercle)
            py = centre_y + rayon * np.sin(angle_cercle)
            
            try:
                th1_cible, th4_cible = cinematique_inverse(px, py)
                ang1_real = angles_actuels["COM6 - S1"][-1]
                ang4_real = angles_actuels["COM5 - S2"][-1]
                current1 = courants_actuels["COM6 - S1"]
                current4 = courants_actuels["COM5 - S2"]
                tau1 = current_to_torque(current1, angles_actuels["COM6 - S1"])
                tau4 = current_to_torque(current4, angles_actuels["COM5 - S2"])
                force1, force4 = dynamique(ang1_real, ang4_real, tau1, tau4)

                corr1 = pid_moteur1.calculer(th1_cible, ang1_real)
                corr2 = pid_moteur2.calculer(th4_cible, ang4_real)
                
                carte_servo1.write(f"S:{th1_cible + corr1 - 7.0}\n".encode())
                carte_servo2.write(f"S:{th4_cible + corr2 + 9.0}\n".encode())       
                
                time.sleep(0.08)
                
                rx, ry = cinematique_directe(ang1_real, ang4_real)
                if rx is not None:
                    donnees_a_sauver.append([
                        time.strftime("%H:%M:%S"), essai_tour_id, px, py, round(th1_cible, 2), round(th4_cible, 2),
                        round(corr1, 2), round(corr2, 2), round(rx, 4), round(ry, 4), round(ang1_real, 2), round(ang4_real, 2), 
                        round(current1, 2), round(current4, 2), round(tau1, 3), round(tau4, 3), round(force1, 3), round(force4, 3)
                    ])
            except Exception as e:
                continue

    sauvegarder_donnees(donnees_a_sauver)
    print(f"Jarvis : {nb_tours} cercles terminés.")
    tracer_analyse(donnees_a_sauver)


# --- 8. Interface & Boucle Principale ---
commande_en_attente = None

def interface_saisie():
    global commande_en_attente
    time.sleep(3) 
    while True:
        entree = input("\nJarvis : 'GO', 'CERCLE', 'FORCE X Y' ou 'X Y' : ").strip().upper()
        commande_en_attente = entree
        while commande_en_attente is not None:
            time.sleep(0.1)

threading.Thread(target=interface_saisie, daemon=True).start()
    
try:
    initialiser_csv()
    print("\nJarvis : Système prêt.")
    carte_servo1.write(f"S:{angles_consigne_init['COM6 - S1']}\n".encode())
    carte_servo2.write(f"S:{angles_consigne_init['COM5 - S2']}\n".encode())
    time.sleep(2)
    carte_servo1.reset_input_buffer(); carte_servo2.reset_input_buffer()
    
    while True:
        if commande_en_attente:
            entree = commande_en_attente
            if entree == "GO": executer_trajectoire()
            elif entree == "CERCLE": executer_cercle()
            elif entree.startswith("FORCE"):
                try:
                    c = entree.split()
                    observer_forces_continu(float(c[1]), float(c[2]))
                except: pass
            else:
                try:
                    c = entree.split()
                    th1, th4 = cinematique_inverse(float(c[0]), float(c[1]))
                    carte_servo1.write(f"S:{th1 - 7.0}\n".encode())
                    carte_servo2.write(f"S:{th4 + 9.0}\n".encode())
                except: pass
            commande_en_attente = None
        else:
            lire_et_afficher(carte_servo1, "COM6 - S1")
            lire_et_afficher(carte_servo2, "COM5 - S2")
            time.sleep(0.005)

except KeyboardInterrupt:
    print("\nFermeture.")
    carte_servo1.close(); carte_servo2.close()
