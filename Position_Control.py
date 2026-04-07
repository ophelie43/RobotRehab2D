import serial
import time
import threading
import numpy as np
from math import atan, sqrt, degrees

# --- 1. Constantes et Modèle Cinématique ---
la, lb, lc = 0.27, 0.31, 0.08

def cinematique_inverse(xc, yc):
    """
    Calcule les angles moteurs (theta1, theta4) (rad) necessaires pour
    atteindre un point (xc, yc) (m) dans le plan. [cite: 699]
    """
    # Equations geometriques pour les deux cotes du mecanisme 
    E1 = -2 * la * xc 
    F1 = -2 * la * yc 
    G1 = la**2 - lb**2 + xc**2 + yc**2 
    
    E4 = 2 * la * (-xc + lc) 
    F4 = -2 * la * yc 
    G4 = lc**2 + la**2 - lb**2 + xc**2 + yc**2 - 2 * lc * xc 
    
    # Resolution par la methode de l arc tangente double 
    # Configuration "+ -" du rapport 
    theta1 = 2 * np.atan((-F1 + np.sqrt(E1**2 + F1**2 - G1**2)) / (G1 - E1)) 
    theta4 = 2 * np.atan((-F4 - np.sqrt(E4**2 + F4**2 - G4**2)) / (G4 - E4)) 
    
    return degrees(theta1), degrees(theta4)

# --- 2. Configuration des Ports Séries ---
try:
    carte_servo1 = serial.Serial('COM9', 115200, timeout=0.05)
    carte_servo2 = serial.Serial('COM5', 115200, timeout=0.05)
    print("Jarvis : Connexion aux deux cartes Servo 2040 établie avec succès.")
except serial.SerialException as e:
    print(f"Jarvis : Erreur matérielle sur les ports COM - {e}")
    exit()

time.sleep(2) 

# --- 3. Variables Globales & Calibration ---
offsets_encodeurs = {"COM9 - S1": None, "COM5 - S2": None}
# NOUVEAU : Définition des angles de départ de votre mécanisme
angles_depart = {"COM9 - S1": 135.0, "COM5 - S2": 45.0}

afficher_en_continu = False 

def lire_et_afficher(port, nom_carte):
    global offsets_encodeurs, afficher_en_continu, angles_depart
    
    if port.in_waiting > 0:
        ligne = port.readline().decode('utf-8').strip()
        
        if ligne.startswith("ENCODEUR:"):
            try:
                angle_brut = float(ligne.split(":")[1])
                
                # SÉQUENCE DE CALIBRATION
                if offsets_encodeurs[nom_carte] is None:
                    offsets_encodeurs[nom_carte] = angle_brut
                    print(f"Jarvis : [CALIBRATION] Encodeur {nom_carte} calé sur {angles_depart[nom_carte]}°.")
                
                # NOUVEAU CALCUL : On applique le décalage et on ajoute la vraie position physique
                angle_relatif = (angle_brut - offsets_encodeurs[nom_carte]) + angles_depart[nom_carte]
                
                if afficher_en_continu:
                    print(f"[{nom_carte}] Position Encodeur : {angle_relatif:.1f}°")
                    
            except ValueError:
                pass
                
        elif ligne:
            print(f"[{nom_carte}] {ligne}")

def envoyer_angle(port, angle):
    commande = f"S:{angle}\n"
    port.write(commande.encode('utf-8'))

def envoyer_xy(xc, yc):
    try:
        theta1, theta4 = cinematique_inverse(xc, yc)
        print(f"\nJarvis : Coordonnées cartésiennes ({xc}, {yc}) valides.")
        print(f"Jarvis : Consignes transmises -> Servo 1 (COM9) : {theta1:.1f}°, Servo 2 (COM5) : {theta4:.1f}°")
        
        envoyer_angle(carte_servo1, theta1)
        envoyer_angle(carte_servo2, theta4)
    except Exception as e:
        print(f"\nJarvis : ERREUR - La position demandée est mécaniquement inatteignable. ({e})")

# --- 4. Interface de Commande ---
def interface_saisie():
    global afficher_en_continu
    time.sleep(3) 
    
    while True:
        entree = input("\nJarvis : Entrez 'X Y' pour viser un point, ou 'STATUT' pour lire les encodeurs : ").strip().upper()
        
        if entree == "STATUT":
            afficher_en_continu = not afficher_en_continu
            etat = "ACTIVÉE" if afficher_en_continu else "DÉSACTIVÉE"
            print(f"Jarvis : Lecture continue des encodeurs {etat}.")
        elif entree:
            try:
                valeurs = entree.split()
                if len(valeurs) == 2:
                    x_cible = float(valeurs[0])
                    y_cible = float(valeurs[1])
                    envoyer_xy(x_cible, y_cible)
                else:
                    print("Jarvis : Format incorrect. Utilisez un espace entre les deux coordonnées.")
            except ValueError:
                print("Jarvis : Caractères invalides. Veuillez entrer des nombres (ex: 0 0.4).")

thread_saisie = threading.Thread(target=interface_saisie, daemon=True)
thread_saisie.start()

# --- 5. Boucle Principale de Communication ---
try:
    print("\nJarvis : Démarrage du système.")
    
    # On force les moteurs à rejoindre leurs angles de départ
    print("Jarvis : Envoi des positions initiales aux servomoteurs...")
    envoyer_angle(carte_servo1, angles_depart["COM9 - S1"])
    envoyer_angle(carte_servo2, angles_depart["COM5 - S2"])
    
    # --- AJOUTE CETTE LIGNE ICI ---
    # On attend 2 secondes que les bras finissent de bouger
    time.sleep(2) 
    
    # Maintenant que les bras sont à 135 et 45, on vide les vieux messages
    carte_servo1.reset_input_buffer()
    carte_servo2.reset_input_buffer()
    
    print("Jarvis : Tampons vidés. Calibration des encodeurs en cours...")
    
    while True:
        lire_et_afficher(carte_servo1, "COM9 - S1")
        lire_et_afficher(carte_servo2, "COM5 - S2")
        time.sleep(0.01)

except KeyboardInterrupt:
    print("\nJarvis : Extinction des systèmes demandée. Fermeture des ports de communication.")
    carte_servo1.close()
    carte_servo2.close()