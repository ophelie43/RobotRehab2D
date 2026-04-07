import serial
import time
import threading
import numpy as np
import matplotlib.pyplot as plt # Nouvelle bibliothèque pour le tracé
from math import degrees, radians

# --- 1. Constantes et Modèles ---
l1, l2, l5 = 0.27, 0.315, 0.08 

def cinematique_directe(theta1_deg, theta4_deg):
    t1, t4 = radians(theta1_deg), radians(theta4_deg)
    try:
        E = 2 * l2 * (l5 + l1 * (np.cos(t4) - np.cos(t1)))
        F = 2 * l1 * l2 * (np.sin(t4) - np.sin(t1))
        G = l5**2 + 2*l1**2 + 2*l5*l1*np.cos(t4) - 2*l5*l1*np.cos(t1) - 2*l1**2*np.cos(t4-t1)
        angle_interne = 2 * np.atan((-F - np.sqrt(E**2 + F**2 - G**2)) / (G - E))
        xc = l5 + l1 * np.cos(t4) + l2 * np.cos(angle_interne)
        yc = l1 * np.sin(t4) + l2 * np.sin(angle_interne)
        return xc, yc
    except:
        return None, None

def cinematique_inverse(xc, yc):
    la, lb, lc = l1, l2, l5
    E1, F1 = -2 * la * xc, -2 * la * yc
    G1 = la**2 - lb**2 + xc**2 + yc**2 
    E4, F4 = 2 * la * (-xc + lc), -2 * la * yc
    G4 = lc**2 + la**2 - lb**2 + xc**2 + yc**2 - 2 * lc * xc 
    th1 = 2 * np.atan((-F1 + np.sqrt(E1**2 + F1**2 - G1**2)) / (G1 - E1)) 
    th4 = 2 * np.atan((-F4 - np.sqrt(E4**2 + F4**2 - G4**2)) / (G4 - E4)) 
    return degrees(th1) - 7, degrees(th4) + 9

# --- 2. Ports Séries ---
try:
    s1_port = serial.Serial('/dev/cu.usbmodem1101', 115200, timeout=0.05)
    s2_port = serial.Serial('/dev/cu.usbmodem101' , 115200, timeout=0.05)
    print("Jarvis : Connexion établie.")
except Exception as e:
    print(f"Erreur COM : {e}"); exit()

# --- 3. Suivi ---
offsets_encodeurs = {"S1": None, "S2": None}
pos_actuelle_deg = {"S1": 135.0, "S2": 45.0}

def lire_encodeurs():
    global pos_actuelle_deg
    while True:
        for p, name, start_angle in [(s1_port, "S1", 135.0), (s2_port, "S2", 45.0)]:
            if p.in_waiting > 0:
                ligne = p.readline().decode('utf-8').strip()
                if ligne.startswith("ENCODEUR:"):
                    try:
                        val = float(ligne.split(":")[1])
                        if offsets_encodeurs[name] is None:
                            offsets_encodeurs[name] = val
                        pos_actuelle_deg[name] = (val - offsets_encodeurs[name]) + start_angle
                    except: pass
        time.sleep(0.01)

threading.Thread(target=lire_encodeurs, daemon=True).start()

# --- 4. Trajectoire et Tracé ---
def executer_trajectoire_custom():
    points = [(-0.1, 0.3), (0.2, 0.3), (0.2, 0.5), (-0.1, 0.5), (-0.1, 0.3)]
    
    # Listes pour stocker les données du graphique
    x_voulu, y_voulu = [], []
    x_reel, y_reel = [], []

    print(f"\n{'Point':<5} | {'X voulu':<8} | {'Y voulu':<8} | {'X réel':<8} | {'Y réel':<8} | {'Erreur (mm)':<10}")
    print("-" * 75)

    for px, py in points:
        try:
            th1_c, th4_c = cinematique_inverse(px, py)
            s1_port.write(f"S:{th1_c}\n".encode())
            s2_port.write(f"S:{th4_c}\n".encode())
            
            time.sleep(1.2)
            
            rx, ry = cinematique_directe(pos_actuelle_deg["S1"], pos_actuelle_deg["S2"])
            
            if rx is not None:
                err = np.sqrt((px-rx)**2 + (py-ry)**2) * 1000
                print(f"P | {px:<8.2f} | {py:<8.2f} | {rx:<8.2f} | {ry:<8.2f} | {err:<10.1f} mm")
                
                # Sauvegarde pour le graphique
                x_voulu.append(px); y_voulu.append(py)
                x_reel.append(rx); y_reel.append(ry)
        except:
            print(f"Point {px},{py} hors de portée.")

    # --- Génération du Graphique ---
    plt.figure(figsize=(8, 6))
    plt.plot(x_voulu, y_voulu, 'r--', label='Trajectoire Voulue', marker='o')
    plt.plot(x_reel, y_reel, 'b-', label='Trajectoire Réelle (Encodeurs)', marker='x')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.title('Comparaison Trajectoire Théorique vs Réelle')
    plt.legend()
    plt.grid(True)
    plt.axis('equal') # Pour garder les proportions du carré
    print("\nJarvis : Affichage du graphique... Fermez la fenêtre pour continuer.")
    plt.show()

# --- 5. Interface ---
try:
    print("Jarvis : Initialisation à 135°/45°...")
    s1_port.write(f"S:{135-7}\n".encode())
    s2_port.write(f"S:{45+9}\n".encode())
    time.sleep(2)
    
    while True:
        cmd = input("\nCommande ('CARRE' ou 'X Y') : ").strip().upper()
        if cmd == "CARRE":
            executer_trajectoire_custom()
        elif cmd:
            parts = cmd.split()
            if len(parts) == 2:
                try:
                    th1, th4 = cinematique_inverse(float(parts[0]), float(parts[1]))
                    s1_port.write(f"S:{th1}\n".encode()); s2_port.write(f"S:{th4}\n".encode())
                except: print("Point inatteignable.")

except KeyboardInterrupt:
    print("\nFermeture."); s1_port.close(); s2_port.close()