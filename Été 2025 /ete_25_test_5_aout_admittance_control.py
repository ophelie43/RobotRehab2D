import serial
import time
import numpy as np
from math import atan, sqrt, degrees
from datetime import datetime
import csv
import math
import os
import asyncio

# -------------------- Global setup --------------------

SERVO1_PORT = 'COM7'  # usba theta1
SERVO2_PORT = 'COM6'  # usbc theta4


l1, l2, l5 = 0.28, 0.31, 0.08

B_virtual = np.array([10.0, 10.0])
K_virtual = np.array([10.0, 10.0])
De = np.array([200.0, 200.0])
Ki = np.array([0.1, 0.1])

F_threshold = 5

time_vec = [0,0.01]
x_vec = [0,0]
y_vec = [0,0]
x_target_vec = [0,0]
y_target_vec = [0,0]
pos1_real_vec = [0,0]
pos4_real_vec = [0,0]

# -------------------- Fonctions outils --------------------

def current_to_torque(current): return ((current + 0.1182)/0.0137) * 10**(-2) * 9.81

def compute_derivative(vec,time_vec):
    derivative = (vec[-1] - 2*vec[-2] + vec[-3])/(time_vec[-1] - time_vec[-2])
    return derivative

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
    return theta1, theta4

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

def unitary_vector(Fx, Fy):
    norm = np.sqrt(Fx**2 + Fy**2)
    return (Fx / norm, Fy / norm) if norm != 0 else (0, 0)

def get_x_target(x, y, Fx_unitary, Fy_unitary):
    return x + 0.001 * Fx_unitary, y + 0.001 * Fy_unitary

# Pour stocker les positions

commanded_angles = []
measured_positions = []
def send_angle(servo, angle, max_attempts=10, tolerance_deg=5.0):
    # Envoie l’angle au servo
    # Active le mode STREAM une seule fois
    
    servo.write(f"SET {angle}\n".encode())
    print(f"Envoi de l'angle {angle}° au servo.")
    time.sleep(0.05)

    attempts = 0
    current = 0.0
    position = 0.0

    while attempts < max_attempts:
        line = servo.readline().decode(errors="replace").strip()
        if not line:
            attempts += 1
            continue
        try:
            current_str, position_str = line.split(",")
            current = float(current_str.strip())
            position = float(position_str.strip())
            return current, position
           # print(f"[{attempts+1}/{max_attempts}] Courant: {current:.3f} A, Position: {position:.2f}°")
            if abs(position - angle) <= tolerance_deg:
                # Sauvegarde pour le graphe
                commanded_angles.append(angle)
                measured_positions.append(position)
                

        except ValueError:
            print(f"⚠️ Ligne ignorée : '{line}'")
        attempts += 1

    print("❌ Position non atteinte dans les limites.")
    commanded_angles.append(angle)
    measured_positions.append(position)
    return current, position

# -------------------- Main control loop --------------------

x, y = 0.0, 0.3
print(f"\n--- Boucle contrôle ---")
print(f"Position actuelle: x={x:.3f}, y={y:.3f}")
    
servo1 = serial.Serial(SERVO1_PORT, 115200, timeout=1)
servo2 = serial.Serial(SERVO2_PORT, 115200, timeout=1)
exercise_duration = 10

servo1.write(b"STREAM\n")
servo2.write(b"STREAM\n")
time.sleep(0.1)

file_path = "C:/Users/maely/OneDrive - polymtl.ca/Stage été 2025/Admittance/admittance_control.csv"

if os.path.exists(file_path):
    os.remove(file_path)
    print(f"✅ Fichier précédent '{file_path}' supprimé.")

write_header = not os.path.exists(file_path)

with open(file_path, mode="a", newline="") as file:
    writer = csv.writer(file)
    if write_header:
        writer.writerow(["timestamp", "x", "y", "tau1", "tau2", "Fx_real", "Fy_real", "Fx_theo", "Fy_theo",
                             "Fx_error", "Fy_error", "x_cmd", "y_cmd", "current1", "current2", 
                             "pos1_real", "pos4_real", "pos1_sent", "pos4_sent", "current_x", "current_y"])

    t_start = time.time()
    while time.time() - t_start < exercise_duration:
        dt_loop = 0.005
        theta1 = compute_theta1(x, y)
        theta4 = compute_theta4(x, y)
        theta1_deg = degrees(theta1)
        theta4_deg = degrees(theta4)

        #print(f"Commande envoyée: theta1={degrees(theta1):.2f}, theta4={degrees(theta4):.2f}")

        current1, pos1_real = send_angle(servo1, theta1_deg)
        current2, pos4_real = send_angle(servo2, theta4_deg)
        pos1_real_vec.append(pos1_real)
        pos4_real_vec.append(pos4_real)
        current_x,current_y = cinematique(pos1_real * np.pi / 180, pos4_real * np.pi / 180)
        print(f"Position réelle: x={current_x:.3f}, y={current_y:.3f}")
        
        x_vec.append(current_x)
        y_vec.append(current_y)
    
        time_vec.append(time.time())

        angular_velocity1 = compute_derivative(pos1_real_vec, time_vec)
        angular_velocity4 = compute_derivative(pos4_real_vec, time_vec)
        if angular_velocity1 < 0:
            current2 = -current2
        if angular_velocity4 < 0:
            current1 = -current1

        tau1 = current_to_torque(current1)
        tau2 = current_to_torque(current2)

        Fx_real, Fy_real = dynamique(pos1_real*np.pi/180, pos4_real*np.pi/180, tau1, tau2)
        if np.sqrt(Fx_real**2 + Fy_real**2) < F_threshold:
            Fx_unitary, Fy_unitary = 0.0, 0.0
        else:
            Fx_unitary, Fy_unitary = unitary_vector(Fx_real, Fy_real)

        x_target, y_target = get_x_target(x, y, Fx_unitary, Fy_unitary)
        x_target_vec.append(x_target)
        y_target_vec.append(y_target)
        vx_target = compute_derivative(x_target_vec,time_vec)
        vy_target = compute_derivative(y_target_vec, time_vec)

        now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        vx_real = compute_derivative(x_vec, time_vec)
        vy_real = compute_derivative(y_vec, time_vec)
        
        Fx_theo = B_virtual[0] * (vx_target - vx_real) + K_virtual[0] * (x_target - x)
        Fy_theo = B_virtual[1] * (vy_target - vy_real) + K_virtual[1] * (y_target - y)

        Fx_error = Fx_theo - Fx_real
        Fy_error = Fy_theo - Fy_real

        x_cmd = np.clip(current_x+ (1 / De[0]) * Fx_error * dt_loop + Ki[0] * dt_loop * Fx_error * dt_loop, -0.2, 0.2)
        y_cmd = np.clip(current_y+ (1 / De[1]) * Fy_error * dt_loop + Ki[1] * dt_loop * Fy_error * dt_loop, 0.23, 0.6)

        writer.writerow([now_str, x, y, tau1, tau2, Fx_real, Fy_real, Fx_theo, Fy_theo,
                             Fx_error, Fy_error, x_cmd, y_cmd, current1, current2,
                             pos1_real, pos4_real, theta1_deg, theta4_deg, current_x, current_y])
        file.flush()
        
        print(current_x, current_y)

        x, y = x_cmd, y_cmd
        del x_vec[0]
        del y_vec[0]
        del x_target_vec[0]
        del y_target_vec[0]
        del time_vec[0]
        time.sleep(dt_loop)

servo1.close()
servo2.close()
print("Contrôle terminé. CSV généré :", os.path.abspath(file_path))
