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
NOM_FICHIER_CSV = "test_admittance.csv"

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
    t1, t4 = radians(theta1), radians(theta4)
    D = (la + lb + lc / 2) / 3
    r1, r2 = la / D, lb / D
    s1, s4 = np.sin(t1), np.sin(t4)
    c1, c4 = np.cos(t1), np.cos(t4)
    A = r1 * s1 - r1 * s4
    B = 2 * (lc / 2) / D + r1 * c1 + r1 * c4
    C = np.pi / 2 + np.arctan(A / B)
    D_val = 8 * r2 * np.sqrt(max(0, 1 - (B**2 + A**2) / (4 * r2**2)))
    sqrt_term = np.sqrt(max(0, 1 - (B**2 + A**2) / (4 * r2**2)))
    E  = r2 * np.sin(C) / (1 + (A**2 / B**2)) * sqrt_term
    Ep = r2 * np.cos(C) / (1 + (A**2 / B**2)) * sqrt_term

    J11 = -((r1 * s1) / 2) - 2 * np.cos(C) * (A * r1 * c1 - B * r1 * s1) / D_val - E / B**2 * (B * r1 * c1 + A * r1 * s1)
    J12 = -r1 * s4 / 2 + 2 * np.cos(C) * (A * r1 * c4 + B * r1 * s4) / D_val - E / B**2 * (-B * r1 * c4 + A * r1 * s4)
    J21 = -r1 * c1 / 2 - 2 * np.sin(C) * (A * r1 * c1 - B * r1 * s1) / D_val + Ep / B**2 * (B * r1 * c1 + A * r1 * s1)
    J22 = -r1 * c4 / 2 + 2 * np.sin(C) * (A * r1 * c4 + B * r1 * s4) / D_val + Ep * B**2 * (-B * r1 * c4 + A * r1 * s4)
    return np.array([[J11, J12], [J21, J22]])

def dynamique(theta1, theta4, tau1, tau2):
    if tau1 is None or tau2 is None:
        return 0.0, 0.0
    J = jacobien(theta1, theta4)
    detJ = np.linalg.det(J)
    if np.abs(detJ) < 1e-3:
        return 0.0, 0.0
    try:
        J_T = np.transpose(J)
        F = np.linalg.solve(J_T, np.array([float(tau1), float(tau2)]))
        return F[0], F[1]
    except:
        return 0.0, 0.0

def current_to_torque(current, theta):
    direction = 1.0
    if len(theta) > 2:
        for i in range(1, len(theta)):
            index = -1 - i
            if (theta[index] - theta[-1]) != 0:
                delta_theta = theta[index] - theta[-1]
                direction = delta_theta / np.abs(delta_theta)
                break
    try:
        return np.abs((float(current) - 0.62) / 0.92) * direction
    except:
        return 0.0


# --- 2. CLASSE PID AVANCÉ (Logique Industrielle / Tustin) ---
class PIDController:
    def __init__(self, Kp, Ki, Kd, alpha=0.05, limite_sortie=20.0):
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
        if dt <= 0.0:
            dt = 0.001
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


class AdmittanceController:
    def __init__(self, De, Ki, Bd, limite_sortie=5.0):
        self.De = De
        self.Ki = Ki
        self.Bd = Bd
        self.limite_sortie = limite_sortie
        self.dernier_temps = time.time()

    def calculer(self, consigne, valeur_mesuree):
        maintenant = time.time()
        dt = maintenant - self.dernier_temps
        if dt <= 0.0:
            dt = 0.001
        Fe = valeur_mesuree - consigne
        Pn = (1 / self.De) * Fe
        In = self.Ki * Fe * dt
        vr = Pn + In
        self.dernier_temps = maintenant
        return vr


class LowPassFilter:
    def __init__(self, f_cutoff, dt):
        self.alpha = 1 - np.exp(-dt * 2 * np.pi * f_cutoff)
        self.last_val = 0.0

    def update(self, val):
        filtered_val = (1 - self.alpha) * self.last_val + self.alpha * val
        self.last_val = filtered_val
        return filtered_val


# Initialisation des PID
pid_moteur1 = PIDController(Kp=0.8, Ki=2.0, Kd=0.05, alpha=0.03, limite_sortie=15.0)
pid_moteur2 = PIDController(Kp=0.8, Ki=2.0, Kd=0.05, alpha=0.03, limite_sortie=15.0)

lp_f1 = LowPassFilter(f_cutoff=2.0, dt=0.08)
lp_f4 = LowPassFilter(f_cutoff=2.0, dt=0.08)


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


# --- 4. Affichage et Analyse ---
def tracer_analyse(donnees):
    if not donnees:
        print("Jarvis : Aucune donnee a tracer.")
        return
    plt.switch_backend('Agg')
    x_voulu = [ligne[2] for ligne in donnees]
    y_voulu = [ligne[3] for ligne in donnees]
    x_reel  = [ligne[8] for ligne in donnees]
    y_reel  = [ligne[9] for ligne in donnees]
    plt.figure(figsize=(10, 6))
    plt.plot(x_voulu, y_voulu, 'r--', label='Consigne')
    plt.plot(x_reel, y_reel, 'b-', label='Reel (Encodeurs)')
    plt.legend()
    plt.title("Analyse de la trajectoire")
    chemin_image = NOM_FICHIER_CSV.replace(".csv", ".png")
    plt.savefig(chemin_image)
    plt.close()
    print(f"Jarvis : Graphique sauvegarde sous {chemin_image}")


# --- 5. Configuration des Ports Series ---
try:
    carte_servo1 = serial.Serial('COM6', 115200, timeout=0.05)  # COM6 et +9 deg
    carte_servo2 = serial.Serial('COM5', 115200, timeout=0.05)  # COM5 et -7 deg
    print("Jarvis : Connexion etablie.")
except serial.SerialException as e:
    print(f"Jarvis : Erreur COM - {e}"); exit()

time.sleep(2)

# --- 6. Variables Globales & Suivi ---
offsets_encodeurs    = {"COM9 - S1": None,  "COM5 - S2": None}
angles_consigne_init = {"COM9 - S1": 135.0 - 7.0, "COM5 - S2": 45.0 + 9.0}
angles_reels_calib   = {"COM9 - S1": 135.0, "COM5 - S2": 45.0}
angles_actuels       = {"COM9 - S1": [135.0], "COM5 - S2": [45.0]}
courants_actuels     = {"COM9 - S1": 0.0,   "COM5 - S2": 0.0}
torques_actuels      = {"COM9 - S1": 0.0,   "COM5 - S2": 0.0}
forces_actuels       = {"COM9 - S1": 0.0,   "COM5 - S2": 0.0}

def lire_et_afficher(port, nom_carte):
    global offsets_encodeurs, angles_actuels
    if port.in_waiting > 0:
        ligne = port.readline().decode('utf-8', errors='ignore').strip()
        if ligne and "ENCODEUR:" in ligne and "| COURANT:" in ligne:
            try:
                parties = ligne.split("|")
                angle   = float(parties[0].split(":")[1])
                courant = float(parties[1].split(":")[1])
                print(f"{nom_carte} >> ENCODEUR: {angle} | COURANT: {courant}")
                if offsets_encodeurs[nom_carte] is None:
                    offsets_encodeurs[nom_carte] = angle
                angles_actuels[nom_carte].append(
                    (angle - offsets_encodeurs[nom_carte]) + angles_reels_calib[nom_carte]
                )
                courants_actuels[nom_carte] = courant
            except:
                pass


# --- 7. Fonctions de Mouvement ---
def executer_trajectoire():
    essai_id = 0
    sommets  = [(0.0, 0.4), (0.1, 0.4), (0.1, 0.3), (-0.1, 0.3), (-0.1, 0.4), (0.0, 0.4)]
    nb_pas   = 60
    donnees_a_sauver = []

    print(f"\nJarvis : Demarrage de la trajectoire fluide avec PID (Essai {essai_id})...")

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
            p_depart  = sommets[i]
            p_arrivee = sommets[i + 1]
            for pas in range(nb_pas):
                ratio = pas / nb_pas
                px = p_depart[0] + (p_arrivee[0] - p_depart[0]) * ratio
                py = p_depart[1] + (p_arrivee[1] - p_depart[1]) * ratio
                essai_id += 1
                try:
                    th1_cible, th4_cible = cinematique_inverse(px, py)
                    ang1_real = angles_actuels["COM9 - S1"][-1]
                    ang4_real = angles_actuels["COM5 - S2"][-1]
                    current1  = courants_actuels["COM9 - S1"]
                    current4  = courants_actuels["COM5 - S2"]
                    tau1 = current_to_torque(current1, angles_actuels["COM9 - S1"])
                    tau4 = current_to_torque(current4, angles_actuels["COM5 - S2"])
                    force1, force4 = dynamique(ang1_real, ang4_real, tau1, tau4)
                    forces_actuels["COM9 - S1"] = force1
                    forces_actuels["COM5 - S2"] = force4
                    corr1 = pid_moteur1.calculer(th1_cible, ang1_real)
                    corr2 = pid_moteur2.calculer(th4_cible, ang4_real)
                    J    = jacobien(ang1_real, ang4_real)
                    detJ = np.linalg.det(J)
                    cmd_th1 = th1_cible + corr1 - 7.0
                    cmd_th4 = th4_cible + corr2 + 9.0
                    carte_servo1.write(f"S:{cmd_th1}\n".encode())
                    carte_servo2.write(f"S:{cmd_th4}\n".encode())
                    print(f"DEBUG: Angles={ang1_real:.2f}/{ang4_real:.2f} | DetJ={detJ:.6f} | F={force1:.2f}/{force4:.2f}")
                    time.sleep(0.08)
                    rx, ry = cinematique_directe(ang1_real, ang4_real)
                    if rx is not None:
                        donnees_a_sauver.append([
                            time.strftime("%Y-%m-%d %H:%M:%S"), essai_id,
                            px, py, round(th1_cible, 2), round(th4_cible, 2),
                            round(corr1, 2), round(corr2, 2),
                            round(rx, 4), round(ry, 4), round(ang1_real, 2), round(ang4_real, 2),
                            round(current1, 2), round(current4, 2),
                            round(tau1, 3), round(tau4, 3), round(force1, 3), round(force4, 3), round(detJ, 4)
                        ])
                except:
                    continue

    sauvegarder_donnees(donnees_a_sauver)
    print(f"Jarvis : Trajectoire terminee. Fichier {NOM_FICHIER_CSV} mis a jour.")


def executer_cercle():
    essai_id = 0
    centre_x, centre_y = 0.0, 0.35
    rayon    = 0.1
    nb_points = 200
    donnees_a_sauver = []

    print(f"\nJarvis : Generation cercle (Essai {essai_id})...")

    th1_init, th4_init = cinematique_inverse(centre_x + rayon, centre_y)
    carte_servo1.write(f"S:{th1_init - 7.0}\n".encode())
    carte_servo2.write(f"S:{th4_init + 9.0}\n".encode())

    print("Jarvis : Attente de stabilisation...")
    time.sleep(2)

    carte_servo1.reset_input_buffer()
    carte_servo2.reset_input_buffer()

    pid_moteur1.reset()
    pid_moteur2.reset()

    for _ in range(1):
        essai_id = 0
        for i in range(nb_points + 1):
            timeout = time.time() + 0.5
            while (carte_servo1.in_waiting == 0 or carte_servo2.in_waiting == 0) and time.time() < timeout:
                time.sleep(0.001)
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
                if not angles_actuels["COM9 - S1"] or not angles_actuels["COM5 - S2"]:
                    continue
                ang1_real = angles_actuels["COM9 - S1"][-1]
                ang4_real = angles_actuels["COM5 - S2"][-1]
                current1  = courants_actuels["COM9 - S1"]
                current4  = courants_actuels["COM5 - S2"]
                tau1 = current_to_torque(current1, np.array(angles_actuels["COM9 - S1"]))
                tau4 = current_to_torque(current4, np.array(angles_actuels["COM5 - S2"]))
                force1, force4 = dynamique(ang1_real, ang4_real, tau1, tau4)
                corr1 = pid_moteur1.calculer(th1_cible, ang1_real)
                corr2 = pid_moteur2.calculer(th4_cible, ang4_real)
                J    = jacobien(ang1_real, ang4_real)
                detJ = np.linalg.det(J)
                cmd_th1 = th1_cible + corr1 - 7.0
                cmd_th4 = th4_cible + corr2 + 9.0
                carte_servo1.write(f"S:{cmd_th1}\n".encode())
                carte_servo2.write(f"S:{cmd_th4}\n".encode())
                print(f"DEBUG: Angles={ang1_real:.2f}/{ang4_real:.2f} | DetJ={detJ:.6f} | F={force1:.2f}/{force4:.2f}")
                rx, ry = cinematique_directe(ang1_real, ang4_real)
                if rx is not None:
                    donnees_a_sauver.append([
                        time.strftime("%Y-%m-%d %H:%M:%S"), essai_id,
                        px, py, round(th1_cible, 2), round(th4_cible, 2),
                        round(corr1, 2), round(corr2, 2),
                        round(rx, 4), round(ry, 4), round(ang1_real, 2), round(ang4_real, 2),
                        round(current1, 4), round(current4, 4),
                        round(tau1, 3), round(tau4, 3), round(force1, 3), round(force4, 3), round(detJ, 4)
                    ])
                time.sleep(0.08)
            except Exception as e:
                print(f"Erreur au point {i}: {e}")
                continue

    sauvegarder_donnees(donnees_a_sauver)
    print(f"Jarvis : Cercle termine. {len(donnees_a_sauver)} points enregistres.")


def admittance_control():
    dt = 0.02
    force_threshold = 0.1
    max_cartesian_speed = 0.08

    admittance_x = AdmittanceController(De=0.3, Ki=0.0, Bd=0.1)
    admittance_y = AdmittanceController(De=0.3, Ki=0.0, Bd=0.1)

    fc = 0.1
    fv = 0.05

    print("\nJarvis : Mode Admittance actif")

    while True:
        try:
            # lecture capteurs
            while carte_servo1.in_waiting > 0:
                lire_et_afficher(carte_servo1, "COM9 - S1")
            while carte_servo2.in_waiting > 0:
                lire_et_afficher(carte_servo2, "COM5 - S2")

            ang1 = angles_actuels["COM9 - S1"][-1]
            ang4 = angles_actuels["COM5 - S2"][-1]

            q = np.array([ang1, ang4])

            # vitesse articulaire
            q_dot = np.zeros(2)

            # courant → torque
            tau1 = current_to_torque(courants_actuels["COM9 - S1"], angles_actuels["COM9 - S1"])
            tau4 = current_to_torque(courants_actuels["COM5 - S2"], angles_actuels["COM5 - S2"])

            tau_ext = np.array([tau1, tau4])

            # force
            f1, f2 = dynamique(ang1, ang4, tau_ext[0], tau_ext[1])
            f = np.array([f1, f2])

            # filtre
            f[0] = lp_f1.update(f[0])
            f[1] = lp_f4.update(f[1])

            # deadzone
            f[np.abs(f) < force_threshold] = 0.0

            # ADMITTANCE → vitesse
            vx = admittance_x.calculer(0, f[0])
            vy = admittance_y.calculer(0, f[1])
            v = np.array([vx, vy])

            # saturation vitesse
            norm = np.linalg.norm(v)
            if norm > max_cartesian_speed:
                v = v / norm * max_cartesian_speed

            # position réelle
            rx, ry = cinematique_directe(ang1, ang4)
            if rx is None:
                continue

            x_meas = np.array([rx, ry])

            # ⭐⭐ POINT CRITIQUE ⭐⭐
            x_cmd = x_meas + v * dt

            # sécurité workspace
            x_cmd[0] = np.clip(x_cmd[0], -0.2, 0.2)
            x_cmd[1] = np.clip(x_cmd[1], 0.2, 0.45)

            th1, th4 = cinematique_inverse(x_cmd[0], x_cmd[1])

            carte_servo1.write(f"S:{th1 - 7.0}\n".encode())
            carte_servo2.write(f"S:{th4 + 9.0}\n".encode())

            print(f"F=({f[0]:.2f},{f[1]:.2f}) V=({vx:.3f},{vy:.3f})")

            time.sleep(dt)

        except KeyboardInterrupt:
            break


# --- 8. Interface & Boucle Principale ---
def interface_saisie():
    time.sleep(3)
    while True:
        entree = input("\nJarvis : 'GO' (carre), 'CERCLE', 'ADMITTANCE' ou 'X Y' : ").strip().upper()
        if entree == "GO":
            executer_trajectoire()
        elif entree == "CERCLE":
            executer_cercle()
        elif entree == "ADMITTANCE":
            admittance_control()
        elif entree:
            try:
                c = entree.split()
                if len(c) == 2:
                    th1, th4 = cinematique_inverse(float(c[0]), float(c[1]))
                    carte_servo1.write(f"S:{th1 - 7.0}\n".encode())
                    carte_servo2.write(f"S:{th4 + 9.0}\n".encode())
            except:
                pass

threading.Thread(target=interface_saisie, daemon=True).start()

try:
    initialiser_csv()
    print("\nJarvis : Systeme pret.")
    carte_servo1.write(f"S:{angles_consigne_init['COM9 - S1']}\n".encode())
    carte_servo2.write(f"S:{angles_consigne_init['COM5 - S2']}\n".encode())
    time.sleep(2)

    carte_servo1.reset_input_buffer()
    carte_servo2.reset_input_buffer()

    while True:
        time.sleep(0.005)

except KeyboardInterrupt:
    print("\nFermeture.")
    carte_servo1.close()
    carte_servo2.close()
