import serial
import time
import threading
import numpy as np
import csv
import os
from math import atan, sqrt, degrees, radians, cos, sin, atan2
import matplotlib.pyplot as plt

# --- 1. Constantes et Modèle Cinématique ---
la, lb, lc = 0.27, 0.315, 0.08

dt = 0.08
NOM_FICHIER_CSV = "/Users/opheliesenechal/Desktop/Bureau - MacBook Pro de Ophélie/UPIR Hiver 2026/Code/rectangle.csv"

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
    # Les angles mathématiques purs (pas d'offsets physiques ici !)
    return degrees(theta1), degrees(theta4)
def jacobien(theta1, theta4):
    D = (la + lb + lc / 2) / 3
    r1, r2 = la / D, lb / D
    s1, s4 = np.sin(theta1), np.sin(theta4)
    c1, c4 = np.cos(theta1), np.cos(theta4)
    A = r1 * s1 - r1 * s4
    B = 2 * (lc / 2) / D + r1 * c1 + r1 * c4
    C = np.pi / 2 + np.arctan(A / B)
    D_val = 8 * r2 * np.sqrt(1 - (B ** 2 + A ** 2) / (4 * r2 ** 2))
    E = r2 * np.sin(C) / (1 + (A ** 2 / B ** 2)) * np.sqrt(1 - (B ** 2 + A ** 2) / (4 * r2 ** 2))
    Ep = r2 * np.cos(C) / (1 + (A ** 2 / B ** 2)) * np.sqrt(1 - (B ** 2 + A ** 2) / (4 * r2 ** 2))

    J11 = -((r1 * s1 )/ 2) - 2 * np.cos(C) * (A * r1 * c1 - B * r1 * s1) / D_val - E / B ** 2 * (B * r1 * c1 + A * r1 * s1)
    J12 = -r1 * s4 / 2 + 2 * np.cos(C) * (A * r1 * c4 + B * r1 * s4) / D_val - E / B ** 2 * (-B * r1 * c4 + A * r1 * s4)
    J21 = -r1 * c1 / 2 - 2 * np.sin(C) * (A * r1 * c1 - B * r1 * s1) / D_val + Ep / B ** 2 * (B * r1 * c1 + A * r1 * s1)
    J22 = -r1 * c4 / 2 + 2 * np.sin(C) * (A * r1 * c4 + B * r1 * s4) / D_val + Ep * B ** 2 * (-B * r1 * c4 + A * r1 * s4)
    return np.array([[J11, J12], [J21, J22]])

def dynamique(theta1, theta4, tau1, tau2):
    J = jacobien(theta1, theta4)
    detJ = np.linalg.det(J)
    if np.abs(detJ) < 1e-2: # Seuil d'instabilité
        return 0.0, 0.0
    else:
        J_T = np.transpose(J)
        F = np.linalg.solve(J_T, np.array([tau1, tau2]))
        return F[0], F[1]
def current_to_torque(current, theta: np.array, nom_carte): 
    for i in range(1,len(theta)):
            index = -1 - i
            if (theta[index] - theta[-1]) != 0:
                delta_theta = theta[index] - theta[-1]
                direction = delta_theta/np.abs(delta_theta)
                break
    if nom_carte == "COM9 - S1":
        return np.abs((current-0.62)/0.92) * direction
    if nom_carte == "COM5- S2":
        return np.abs((current-0.218)/2.011) * direction

# --- 2. CLASSE PID AVANCÉ (Logique Industrielle / Tustin) ---
class PIDController:
    def __init__(self, Kp, Ki, Kd, alpha=0.05, limite_sortie=20.0):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.alpha = alpha  # Constante de temps du filtre passe-bas
        self.limite_sortie = limite_sortie 
        
        # Variables d'état (mémoires de l'instant n-1)
        self.en_1 = 0.0
        self.In_1 = 0.0
        self.Dn_1 = 0.0
        
        self.dernier_temps = time.time()

    def calculer(self, consigne, valeur_mesuree):
        maintenant = time.time()
        dt = maintenant - self.dernier_temps
        if dt <= 0.0:
            dt = 0.001

        # Calcul de l'erreur actuelle
        en = consigne - valeur_mesuree
        
        # 1. Action Proportionnelle
        Pn = self.Kp * en
        
        # 2. Action Dérivée avec filtre de Tustin
        C1 = (2.0 * self.alpha - dt) / (2.0 * self.alpha + dt)
        C2 = (2.0 * self.Kd) / (2.0 * self.alpha + dt)
        Dn = (C1 * self.Dn_1) + (C2 * (en - self.en_1))
        
        # 3. Anti-Windup Conditionnel
        Un_temp = Pn + self.In_1 + Dn
        In = self.In_1
        
        if abs(Un_temp) < self.limite_sortie:
            In = self.In_1 + (dt / 2.0) * self.Ki * (en + self.en_1)
            
        # 4. Commande finale non-saturée
        Un = Pn + In + Dn
        
        # 5. Saturation finale de sécurité
        sortie = max(min(Un, self.limite_sortie), -self.limite_sortie)
        
        # 6. Mise à jour des mémoires
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
class AdmittanceController:
    def __init__(self, De, Ki, limite_sortie = 5.0):
        self.De = De
        self.Ki = Ki
        self.limite_sortie = limite_sortie

        # Variables d'état
        self.dernier_temps = time.time()
    def calculer(self,consigne, valeur_mesurée):
        maintenant = time.time()
        dt = maintenant - self.dernier_temps
        if dt <- 0.0:
            dt = 0.001
        
        #Calcul de l'erreur de force
        Fe = consigne - valeur_mesurée

        # 1. Action proportionnelle
        Pn = (1/self.De) * Fe

        # 2. Action intégrale
        In = self.Ki * Fe * dt

        # 3. Sortie
        vr = Pn + In
        return vr
class LowPassFilter:
    def __init__(self, f_cutoff, dt):
        self.alpha = 1- np.exp(-dt * 2 * np.pi * f_cutoff)
        self.last_val = 0.0

    def update(self, val):
        filtered_val = (1 - self.alpha) * self.last_val + self.alpha * val
        self.last_val = filtered_val
        return filtered_val

# Initialisation des deux PID (Gains ajustés pour un mouvement calme et précis)
pid_moteur1 = PIDController(Kp=0.8, Ki=2.0, Kd=0.05, alpha=0.03, limite_sortie=15.0)
pid_moteur2 = PIDController(Kp=0.8, Ki=2.0, Kd=0.05, alpha=0.03, limite_sortie=15.0)

admitance_moteur1 = AdmittanceController(De = 1.0, Ki = 1.0)
admitance_moteur2 = AdmittanceController(De = 1.0, Ki = 1.0)

lp_f1 = LowPassFilter(f_cutoff = 5.0, dt = dt)
lp_f4 = LowPassFilter(f_cutoff = 5.0, dt = dt)
# --- 3. Fonctions CSV ---
def initialiser_csv():
    if not os.path.exists(NOM_FICHIER_CSV):
        with open(NOM_FICHIER_CSV, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                "Horodatage", "Essai_ID", "X_Voulu", "Y_Voulu", "Theta1_Cible", "Theta4_Cible",
                "Correction_M1", "Correction_M2", 
                "X_Encodeur", "Y_Encodeur", "Theta1_Encodeur", "Theta4_Encodeur",
                "Courant1", "Courant4", "Tau1", "Tau4", "Force1", "Force4", "Jacobien"
            ])

def sauvegarder_donnees(data_list):
    with open(NOM_FICHIER_CSV, mode='a', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data_list)

# --- 4. Affichage et Analyse (Matplotlib) ---

def tracer_analyse(donnees):
    if not donnees:
        print("Jarvis : Aucune donnée à tracer.")
        return

    # Désactiver l'affichage interactif pour éviter le crash de thread
    plt.switch_backend('Agg') 
    
    x_voulu = [ligne[2] for ligne in donnees]
    y_voulu = [ligne[3] for ligne in donnees]
    x_reel = [ligne[8] for ligne in donnees]
    y_reel = [ligne[9] for ligne in donnees]
    
    plt.figure(figsize=(10, 6))
    plt.plot(x_voulu, y_voulu, 'r--', label='Consigne')
    plt.plot(x_reel, y_reel, 'b-', label='Réel (Encodeurs)')
    plt.legend()
    plt.title("Analyse de la trajectoire")
    
    # On sauvegarde au lieu de plt.show()
    chemin_image = NOM_FICHIER_CSV.replace(".csv", ".png")
    plt.savefig(chemin_image)
    plt.close()
    print(f"Jarvis : Graphique sauvegardé sous {chemin_image}")
# --- 5. Configuration des Ports Séries ---
try:
    carte_servo1 = serial.Serial('/dev/tty.usbmodem101', 115200, timeout=0.05)# COM9 et + 9°
    carte_servo2 = serial.Serial('/dev/tty.usbmodem1101', 115200, timeout=0.05) # COM5 et -7°
    print("Jarvis : Connexion établie.")
except serial.SerialException as e:
    print(f"Jarvis : Erreur COM - {e}"); exit()

time.sleep(2) 

# --- 6. Variables Globales & Suivi ---
offsets_encodeurs = {"COM9 - S1": None, "COM5 - S2": None}
angles_consigne_init = {"COM9 - S1": 135.0 - 7.0, "COM5 - S2": 45.0 + 9.0}
angles_reels_calib = {"COM9 - S1": 135.0, "COM5 - S2": 45.0}    
angles_actuels = {"COM9 - S1": [135.0], "COM5 - S2": [45.0]}
courants_actuels = {"COM9 - S1": 0.0, "COM5 - S2": 0.0}
torques_actuels = {"COM9 - S1": 0.0, "COM5 - S2": 0.0}
forces_actuels = {"COM9 - S1": 0.0, "COM5 - S2": 0.0}

def lire_et_afficher(port, nom_carte):
    global offsets_encodeurs, angles_actuels
    if port.in_waiting > 0:
        ligne = port.readline().decode('utf-8', errors='ignore').strip()
        if ligne and "ENCODEUR:" in ligne and "| COURANT:" in ligne:
            
            try:
                parties = ligne.split("|")
                angle = float(parties[0].split(":")[1])
                courant = float(parties[1].split(":")[1])
                print(f"{nom_carte} >> ENCODEUR: {angle} | COURANT: {courant}")
                if offsets_encodeurs[nom_carte] is None:
                    offsets_encodeurs[nom_carte] = angle


                # L'angle calculé ici est la mathématique pure
                angles_actuels[nom_carte].append((angle - offsets_encodeurs[nom_carte]) + angles_reels_calib[nom_carte])
                courants_actuels[nom_carte] = courant
                
            except: pass

# --- 7. Fonctions de Mouvement ---
def executer_trajectoire():
    essai_id = 0 #time.strftime("%H%M%S")
    sommets = [(0.0, 0.4), (0.1, 0.4), (0.1, 0.3), (-0.1, 0.3), (-0.1, 0.4), (0.0, 0.4)]
    nb_pas = 60 # Plus de pas pour une trajectoire fluide à haute vitesse
    donnees_a_sauver = []
    
    print(f"\nJarvis : Démarrage de la trajectoire fluide avec PID (Essai {essai_id})...")
    
    # Approche en douceur du premier point (Homing)
    point_depart = sommets[0]
    th1_init, th4_init = cinematique_inverse(point_depart[0], point_depart[1])
    carte_servo1.write(f"S:{th1_init - 7.0}\n".encode())
    lire_et_afficher(carte_servo1, "COM9 - S1")

    carte_servo2.write(f"S:{th4_init + 9.0}\n".encode())
    lire_et_afficher(carte_servo2, "COM5 - S2")
    
    print("Jarvis : Mise en position initiale. Attente de stabilisation...")
    time.sleep(2) 
    
    pid_moteur1.reset()
    pid_moteur2.reset()
    
    for n in range(1):
        essai_id = 0
        for i in range(len(sommets) - 1):
            p_depart = sommets[i]
            p_arrivee = sommets[i+1]
            
            for pas in range(nb_pas):
                ratio = pas / nb_pas
                px = p_depart[0] + (p_arrivee[0] - p_depart[0]) * ratio
                py = p_depart[1] + (p_arrivee[1] - p_depart[1]) * ratio
                essai_id += 1
                
                try:
                    th1_cible, th4_cible = cinematique_inverse(px, py)
                    ang1_real = angles_actuels["COM9 - S1"][-1]
                    ang4_real = angles_actuels["COM5 - S2"][-1]

                    #courants_actuels[carte_servo1] = courant
                    current1 = courants_actuels["COM9 - S1"]
                    current4 = courants_actuels["COM5 - S2"]

                    tau1 = current_to_torque(current1, angles_actuels["COM9 - S1"])
                    tau4 = current_to_torque(current4, angles_actuels["COM5 - S2"])

                    force1, force4 = dynamique(ang1_real, ang4_real,tau1,tau4)
                    forces_actuels["COM9 - S1"] = force1
                    forces_actuels["COM5 - S2"] = force4

                    corr1 = pid_moteur1.calculer(th1_cible, ang1_real)
                    corr2 = pid_moteur2.calculer(th4_cible, ang4_real)
                    
                    J = jacobien(ang1_real, ang4_real)
                    detJ = np.linalg.det(J)
                
                    cmd_th1 = th1_cible + corr1 - 7.0
                    cmd_th4 = th4_cible + corr2 + 9.0
                    
                    carte_servo1.write(f"S:{cmd_th1}\n".encode())
                    carte_servo2.write(f"S:{cmd_th4}\n".encode())
                    


                    time.sleep(0.08) # Boucle rapide de 80ms
                    
                    rx, ry = cinematique_directe(ang1_real, ang4_real)
                    if rx is not None:
                        donnees_a_sauver.append([
                            time.strftime("%Y-%m-%d %H:%M:%S"), essai_id,
                            px, py, round(th1_cible, 2), round(th4_cible, 2),
                            round(corr1, 2), round(corr2, 2),
                            round(rx, 4), round(ry, 4), round(ang1_real, 2), round(ang4_real, 2), round(current1, 2), round(current4, 2),
                            round(tau1, 3), round(tau4, 3), round(force1, 3), round(force4, 3), round(detJ, 4)
                        ])
                except:
                    continue 

    sauvegarder_donnees(donnees_a_sauver)
    print(f"Jarvis : Trajectoire terminée. Fichier {NOM_FICHIER_CSV} mis à jour.")
    # tracer_analyse(donnees_a_sauver)


def executer_cercle():
    essai_id = 0
    centre_x, centre_y = 0.0, 0.35 
    rayon = 0.1 
    nb_points = 200 # Plus de points pour garder une vitesse constante
    donnees_a_sauver = []
    
    print(f"\nJarvis : Génération cercle avec PID (Essai {essai_id})...")
    
    th1_init, th4_init = cinematique_inverse(centre_x + rayon, centre_y)
    carte_servo1.write(f"S:{th1_init - 7.0}\n".encode())
    lire_et_afficher(carte_servo1, "COM9 - S1")
    carte_servo2.write(f"S:{th4_init + 9.0}\n".encode())
    lire_et_afficher(carte_servo2, "COM5 - S2")

    print("Jarvis : Mise en position pour le cercle. Attente de stabilisation...")
    time.sleep(2)
    carte_servo1.reset_input_buffer()
    carte_servo2.reset_input_buffer()

    pid_moteur1.reset()
    pid_moteur2.reset()
    
    for n in range(1):
        essai_id = 0
        for i in range(nb_points + 1):
            timeout = time.time() + 0.5
            while (carte_servo1.in_waiting == 0 or carte_servo2.in_waiting == 0) and time.time() < timeout:
                time.sleep(0.001)

            # On lit toutes les données disponibles pour avoir la plus récente
            while carte_servo1.in_waiting > 0:
                lire_et_afficher(carte_servo1, "COM9 - S1")
            while carte_servo2.in_waiting > 0:
                lire_et_afficher(carte_servo2, "COM5 - S2")

            angle_cercle = (2 * np.pi * i) / nb_points
            px = centre_x + rayon * np.cos(angle_cercle)
            py = centre_y + rayon * np.sin(angle_cercle)
            essai_id += 1
            
            try:
                th1_cible, th4_cible = cinematique_inverse(px, py)

                # Sécurité : on vérifie que les listes ne sont pas vides
                if not angles_actuels["COM9 - S1"] or not angles_actuels["COM5 - S2"]:
                    continue

                ang1_real = angles_actuels["COM9 - S1"][-1]
                ang4_real = angles_actuels["COM5 - S2"][-1]

                #courants_actuels[carte_servo1] = courant
                current1 = courants_actuels["COM9 - S1"]
                current4 = courants_actuels["COM5 - S2"]

                tau1 = current_to_torque(current1, angles_actuels["COM9 - S1"], "COM0 - S1")
                tau4 = current_to_torque(current4, angles_actuels["COM5 - S2"], "COM5 - S2")

                force1, force4 = dynamique(ang1_real, ang4_real,tau1,tau4)

                forces_actuels["COM9 - S1"] = force1
                forces_actuels["COM5 - S2"] = force4

                corr1 = pid_moteur1.calculer(th1_cible, ang1_real)
                corr2 = pid_moteur2.calculer(th4_cible, ang4_real)
                
                cmd_th1 = th1_cible + corr1 - 7.0
                cmd_th4 = th4_cible + corr2 + 9.0

                carte_servo1.write(f"S:{cmd_th1}\n".encode())
                carte_servo2.write(f"S:{cmd_th4}\n".encode())  

                J = jacobien(ang1_real, ang4_real)
                detJ = np.linalg.det(J)
                
                time.sleep(0.08) # Boucle rapide de 20ms
                
                rx, ry = cinematique_directe(ang1_real, ang4_real)
                if rx is not None:
                    donnees_a_sauver.append([
                        time.strftime("%Y-%m-%d %H:%M:%S"), essai_id,
                        px, py, round(th1_cible, 2), round(th4_cible, 2),
                        round(corr1, 2), round(corr2, 2),
                        round(rx, 4), round(ry, 4), round(ang1_real, 2), round(ang4_real, 2), round(current1, 2), round(current4, 2),
                        round(tau1, 3), round(tau4, 3), round(force1, 3), round(force4, 3), round(detJ, 4)
                    ])
            except:
                continue

    sauvegarder_donnees(donnees_a_sauver)
    print(f"Jarvis : Cercle terminé. {len(donnees_a_sauver)} points enregistrés.")
    # tracer_analyse(donnees_a_sauver)
    
def admittance_control():
    essai_id = 0
    centre_x, centre_y = 0.0, 0.35 
    rayon = 0.1 
    nb_points = 200 # Plus de points pour garder une vitesse constante
    donnees_a_sauver = []
    
    print(f"\nJarvis : Position initiale (Essai {essai_id})...")
    
    th1_init, th4_init = cinematique_inverse(centre_x + rayon, centre_y)
    carte_servo1.write(f"S:{th1_init - 7.0}\n".encode())
    carte_servo2.write(f"S:{th4_init + 9.0}\n".encode())

    print("Jarvis : Attente de stabilisation...")
    time.sleep(2)

    carte_servo1.reset_input_buffer()
    carte_servo2.reset_input_buffer()

    th1_cible, th4_cible = th1_init, th4_init
    print("Jarvis : Mise en position pour le cercle. Attente de stabilisation...")
    time.sleep(2)
    
    pid_moteur1.reset()
    pid_moteur2.reset()
    
    for _ in range(1):
        essai_id = 0
        for i in range(nb_points + 1):

            timeout = time.time() + 0.5
            while (carte_servo1.in_waiting == 0 or carte_servo2.in_waiting == 0) and time.time() < timeout:
                time.sleep(0.001)

            # On lit toutes les données disponibles pour avoir la plus récente
            while carte_servo1.in_waiting > 0:
                lire_et_afficher(carte_servo1, "COM9 - S1")

            while carte_servo2.in_waiting > 0:
                lire_et_afficher(carte_servo2, "COM5 - S2")

            ang1_real = angles_actuels["COM9 - S1"][-1]
            ang4_real = angles_actuels["COM5 - S2"][-1]

            #courants_actuels[carte_servo1] = courant
            current1 = courants_actuels["COM9 - S1"]
            current4 = courants_actuels["COM5 - S2"]

            tau1 = current_to_torque(current1, angles_actuels["COM9 - S1"])
            tau4 = current_to_torque(current4, angles_actuels["COM5 - S2"])

            force1_brut, force4_brut = dynamique(ang1_real, ang4_real,tau1,tau4)
            force1_clean = lp_f1.update(force1_brut)
            force4_clean = lp_f4.update(force4_brut)

            forces_actuels["COM9 - S1"] = force1_clean
            forces_actuels["COM5 - S2"] = force4_clean

            v1 = admitance_moteur1.calculer(0, force1_clean)
            v4 = admitance_moteur2.calculer(0, force4_clean)
            th1_cible = th1_cible + (v1 * dt) 
            th4_cible = th4_cible + (v4 * dt)
            essai_id += 1
            
            try:



                cmd_th1 = th1_cible - 7.0
                cmd_th4 = th4_cible + 9.0
                # Sécurité : on vérifie que les listes ne sont pas vides
                if not angles_actuels["COM9 - S1"] or not angles_actuels["COM5 - S2"]:
                    continue
                px, py = cinematique_directe(cmd_th1, cmd_th4)
                
                J = jacobien(ang1_real, ang4_real)
                detJ = np.linalg.det(J)
                
                carte_servo1.write(f"S:{cmd_th1}\n".encode())
                carte_servo2.write(f"S:{cmd_th4}\n".encode())       
                
                time.sleep(dt)
                
                rx, ry = cinematique_directe(ang1_real, ang4_real)
                if rx is not None:
                    donnees_a_sauver.append([
                        time.strftime("%Y-%m-%d %H:%M:%S"), essai_id,
                        px, py, round(th1_cible, 2), round(th4_cible, 2),
                        round(th1_cible, 2), round(th4_cible, 2),
                        round(rx, 4), round(ry, 4), round(ang1_real, 2), round(ang4_real, 2), round(current1, 2), round(current4, 2),
                        round(tau1, 3), round(tau4, 3), round(force1_clean, 3), round(force4_clean, 3), round(detJ, 4)
                    ])
            except:
                continue

    sauvegarder_donnees(donnees_a_sauver)
    print(f"Jarvis : Cercle terminé.")
    tracer_analyse(donnees_a_sauver)
# --- 8. Interface & Boucle Principale ---
def interface_saisie():
    time.sleep(3) 
    while True:
        entree = input("\nJarvis : 'GO' (carré), 'CERCLE', 'ADMITTANCE' ou 'X Y' : ").strip().upper()
        if entree == "GO":
            executer_trajectoire()
        elif entree == "CERCLE":
            executer_cercle()
        elif entree == 'ADMITTANCE':
            admittance_control()
        elif entree:
            try:
                c = entree.split()
                if len(c) == 2:
                    th1, th4 = cinematique_inverse(float(c[0]), float(c[1]))
                    carte_servo1.write(f"S:{th1 - 7.0}\n".encode())
                    carte_servo2.write(f"S:{th4 + 9.0}\n".encode())
            except: pass

threading.Thread(target=interface_saisie, daemon=True).start()
    
try:
    initialiser_csv()
    print("\nJarvis : Système prêt.")
    carte_servo1.write(f"S:{angles_consigne_init['COM9 - S1']}\n".encode())
    carte_servo2.write(f"S:{angles_consigne_init['COM5 - S2']}\n".encode())
    time.sleep(2)
    
    carte_servo1.reset_input_buffer()
    carte_servo2.reset_input_buffer()
    
    while True:
        try: 
            lire_et_afficher(carte_servo1, "COM9 - S1")
            lire_et_afficher(carte_servo2, "COM5 - S2")
        except:
                pass
        time.sleep(0.005)

except KeyboardInterrupt:
    print("\nFermeture."); carte_servo1.close(); carte_servo2.close()