import serial
import time
import threading
import numpy as np
import csv
from math import degrees, radians, cos, sin, atan2, sqrt

# --- 1. Constantes et Modèle Cinématique ---
la, lb, lc = 0.27, 0.315, 0.08
filename = "lecture_encodeur.csv"
data = []

def create_csv(filename):
    with open(filename, 'w', newline='') as csvfile:
        fieldnames = [
            'commande_x_envoyee', 'commande_y_envoyee', 
            'theta_1_envoye', 'theta_4_envoye', 
            'pos_réel_x', 'pos_réel_y', 
            'theta_1_reel', 'theta_4_reel'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)
        print(f"Jarvis : Fichier csv {filename} écrit avec succès.")

def cinematique_inverse(xc, yc):
    """Calcule les angles (consignes moteurs avec offsets) pour atteindre (x, y)."""
    E1, F1 = -2 * la * xc, -2 * la * yc
    G1 = la**2 - lb**2 + xc**2 + yc**2 
    
    E4, F4 = 2 * la * (-xc + lc), -2 * la * yc
    G4 = lc**2 + la**2 - lb**2 + xc**2 + yc**2 - 2 * lc * xc 
    
    # Utilisation de atan2 pour la robustesse
    theta1 = 2 * atan2((-F1 + sqrt(max(0, E1**2 + F1**2 - G1**2))), (G1 - E1)) 
    theta4 = 2 * atan2((-F4 - sqrt(max(0, E4**2 + F4**2 - G4**2))), (G4 - E4)) 
    
    # Retourne les consignes avec offsets (-7 et +9)
    return degrees(theta1) - 7, degrees(theta4) + 9

def cinematique(theta1, theta4): 
    # Conversion en radians
    t1, t4 = radians(theta1), radians(theta4)
    
    try:
        # Correction des puissances (**2)
        E = 2 * lb * (lc + la * (cos(t4) - cos(t1)))
        F = 2 * la * lb * (sin(t4) - sin(t1))
        G = lc**2 + 2 * la**2 + 2 * lc * la * cos(t4) - 2 * lc * la * cos(t1) - 2 * la**2 * cos(t4 - t1)
        
        # Calcul de l'argument de la racine avec sécurité
        arg = max(0, E**2 + F**2 - G**2)
        
        # Utilisation de atan2 pour éviter les divisions par zéro
        angle_interne = 2 * atan2((-F - sqrt(arg)), (G - E))
        
        xc = lc + la * cos(t4) + lb * cos(angle_interne)
        yc = la * sin(t4) + lb * sin(angle_interne)
        return xc, yc
    except Exception as e:
        print(f"Erreur cinématique directe : {e}")
        return None, None

# --- 2. Configuration des Ports Séries ---
try:
    carte_servo1 = serial.Serial('COM9', 115200, timeout=0.05)
    carte_servo2 = serial.Serial('COM5', 115200, timeout=0.05)
    print("Jarvis : Connexion aux deux cartes Servo 2040 établie.")
except serial.SerialException as e:
    print(f"Jarvis : Erreur matérielle sur les ports COM - {e}")
    exit()

time.sleep(2) 

# --- 3. Variables Globales & Calibration ---
offsets_encodeurs = {"COM9 - S1": None, "COM5 - S2": None}
angles_consigne = {"COM9 - S1": 135.0 - 7, "COM5 - S2": 45.0 + 9}
angles_reels = {"COM9 - S1": 135.0, "COM5 - S2": 45.0}
afficher_en_continu = False 

def lire_et_afficher(port, nom_carte):
    global offsets_encodeurs, angles_reels
    angle_relatif = None
    if port.in_waiting > 0:
        ligne = port.readline().decode('utf-8').strip()
        if ligne.startswith("ENCODEUR:"):
            try:
                angle_brut = float(ligne.split(":")[1])
                if offsets_encodeurs[nom_carte] is None:
                    offsets_encodeurs[nom_carte] = angle_brut
                    print(f"Jarvis : [CALIBRATION] {nom_carte} calé sur {angles_reels[nom_carte]}°.")
                
                angle_relatif = (angle_brut - offsets_encodeurs[nom_carte]) + angles_reels[nom_carte]
                if afficher_en_continu:
                    print(f"[{nom_carte}] Position : {angle_relatif:.1f}°")
            except ValueError: pass
    return angle_relatif

def envoyer_angle(port, angle):
    port.write(f"S:{angle}\n".encode('utf-8'))

def envoyer_xy(xc, yc):
    th1, th4 = cinematique_inverse(xc, yc)
    print(f"\nJarvis : Vise ({xc}, {yc}) -> S1:{th1:.1f}°, S2:{th4:.1f}°")
    envoyer_angle(carte_servo1, th1)
    envoyer_angle(carte_servo2, th4)
    return th1, th4

# --- 4. Interface ---
def interface_saisie():
    global afficher_en_continu
    time.sleep(3) 
    while True:
        entree = input("\nJarvis : 'X Y' ou 'STATUT' : ").strip().upper()
        if entree == "STATUT":
            afficher_en_continu = not afficher_en_continu
        elif entree:
            try:
                v = entree.split()
                if len(v) == 2: envoyer_xy(float(v[0]), float(v[1]))
            except ValueError: pass

threading.Thread(target=interface_saisie, daemon=True).start()

# --- 5. Boucle Principale ---
try:
    print("\nJarvis : Initialisation...")
    envoyer_angle(carte_servo1, angles_consigne["COM9 - S1"])
    envoyer_angle(carte_servo2, angles_consigne["COM5 - S2"])
    time.sleep(2)
    
    carte_servo1.reset_input_buffer()
    carte_servo2.reset_input_buffer()

    points_trajectoire = [(0.0, 0.4), (0.1, 0.4), (0.1, 0.3), (0.0, 0.3), (0.0, 0.4)]

    for x_cmd, y_cmd in points_trajectoire:
        # On récupère les consignes envoyées pour les logger
        th1_env, th4_env = envoyer_xy(x_cmd, y_cmd)
        time.sleep(1.5) 
        
        t1_reel, t4_reel = None, None
        while t1_reel is None: t1_reel = lire_et_afficher(carte_servo1, "COM9 - S1")
        while t4_reel is None: t4_reel = lire_et_afficher(carte_servo2, "COM5 - S2")
            
        xc, yc = cinematique(t1_reel, t4_reel)

        data.append({
            'commande_x_envoyee': x_cmd, 
            'commande_y_envoyee': y_cmd, 
            'theta_1_envoye': th1_env, 
            'theta_4_envoye': th4_env, 
            'pos_réel_x': xc, 
            'pos_réel_y': yc,
            'theta_1_reel': t1_reel,
            'theta_4_reel': t4_reel,
        })

    create_csv(filename)
    print("Jarvis : Prêt pour commandes manuelles.")
    
    while True:
        lire_et_afficher(carte_servo1, "COM9 - S1")
        lire_et_afficher(carte_servo2, "COM5 - S2")
        time.sleep(0.01)

except KeyboardInterrupt:
    print("\nJarvis : Fermeture des ports.")
    carte_servo1.close()
    carte_servo2.close()