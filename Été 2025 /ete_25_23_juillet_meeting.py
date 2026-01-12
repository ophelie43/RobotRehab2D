import serial
import time
import numpy as np
from math import atan, sqrt, degrees, fmod
from datetime import datetime
import csv
import math
import os
import asyncio
from bleak import BleakClient

# -------------------- Global setup --------------------

SERVO1_PORT = 'COM6'  # usbc theta1
SERVO2_PORT = 'COM7'  # usba theta4

l1, l2, l5 = 0.28, 0.31, 0.08

B_virtual = np.array([10.0, 10.0])
K_virtual = np.array([10.0, 10.0])
De = np.array([30.0, 30.0])
Ki = np.array([0.1, 0.1])

F_threshold = 3

# Global IMU values
gz1 = 0.0
gz2 = 0.0

# -------------------- Classes --------------------

class IMU_BLE:
    WIT_READ_CHAR_UUID = "0000ffe4-0000-1000-8000-00805f9a34fb"

    def __init__(self, address, label):
        self.address = address
        self.label = label.upper()
        self.client = BleakClient(self.address)
        self.gyro_z_history = []
        self.offset = 0

    def to_signed(self, val):
        return val - 65536 if val > 32767 else val

    def decode_packet(self, data):
        global gz1, gz2
        if len(data) != 20 or data[0] != 0x55 or data[1] != 0x61:
            return
        def get_val(lo, hi):
            return self.to_signed((hi << 8) | lo)
        gyro_z = get_val(data[12], data[13]) * (2000.0 / 32768.0)
        self.gyro_z_history.append(gyro_z)
        if self.label == "IMU1":
            gz1 = gyro_z
        elif self.label == "IMU2":
            gz2 = gyro_z

    def handle_data(self, sender, data):
        for i in range(0, len(data), 20):
            packet = data[i:i+20]
            self.decode_packet(packet)

    async def connect_and_start(self):
        await self.client.connect()
        await self.client.start_notify(self.WIT_READ_CHAR_UUID, self.handle_data)

    async def stop(self):
        await self.client.stop_notify(self.WIT_READ_CHAR_UUID)
        await self.client.disconnect()

# -------------------- Helper functions --------------------

def current_to_torque1(current): return 50 * current #86.511 * current #10
def current_to_torque2(current): return 35 * current #67.089 * current #5

def jacobien(theta1, theta4):
    D = (l1 + l2 + l5 / 2) / 3
    r1, r2 = l1 / D, l2 / D
    s1, s4 = np.sin(theta1), np.sin(theta4)
    c1, c4 = np.cos(theta1), np.cos(theta4)
    A = r1 * s1 - r1 * s4
    B = 2 * (l5 / 2) / D + r1 * c1 + r1 * c4
    C = np.pi / 2 + np.arctan(A / B)
    D_val = 8 * r2 * np.sqrt(1 - (B ** 2 + A ** 2) / (4 * r2 ** 2))
    E = r2 * np.sin(C) / (1 + A ** 2 / B ** 2) * np.sqrt(1 - (B ** 2 + A ** 2) / (4 * r2 ** 2))
    Ep = r2 * np.cos(C) / (1 + A ** 2 / B ** 2) * np.sqrt(1 - (B ** 2 + A ** 2) / (4 * r2 ** 2))
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

#SERVO1_PORT = 'COM6'  # usbc
#SERVO2_PORT = 'COM7'  # usba
#servo1 = serial.Serial(SERVO1_PORT, 115200, timeout=1)
#servo2 = serial.Serial(SERVO2_PORT, 115200, timeout=1)

def send_angle(servo, angle, max_attempts=10):
    servo.write(f"{angle}\n".encode())
    attempts = 0
    while attempts < max_attempts:
        line = servo.readline().decode(errors="replace").strip()
        try:
            # Attendu: ligne au format "courant, position"
            current_str, position_str = line.split(",")
            current = float(current_str.strip())
            position = float(position_str.strip())
            return current, position
        except:
            attempts += 1
    print(f"Warning: no valid response from servo after {max_attempts} tries.")
    return 0.0, 0.0  # valeur par défaut pour continuer le programme


# -------------------- Main control loop --------------------

async def run_control(imu1, imu2):
    x, y = 0.0, 0.3
    print(f"\n--- Boucle contrôle ---")
    print(f"Position actuelle: x={x:.3f}, y={y:.3f}")
    print(f"Gyro IMU1={gz1:.2f}°/s, Gyro IMU2={gz2:.2f}°/s")
    
    servo1 = serial.Serial(SERVO1_PORT, 115200, timeout=1)
    servo2 = serial.Serial(SERVO2_PORT, 115200, timeout=1)
    exercise_duration = 40

    file_path = "C:/Users/maely/OneDrive - polymtl.ca/Stage été 2025/Admittance/admittance_control.csv"

    if os.path.exists(file_path):
        os.remove(file_path)
        print(f"✅ Fichier précédent '{file_path}' supprimé.")
        
    write_header = not os.path.exists(file_path)

    with open(file_path, mode="a", newline="") as file:
        writer = csv.writer(file)
        if write_header:
            writer.writerow(["timestamp", "x", "y", "vx", "vy", "Fx_real", "Fy_real", "Fx_theo", "Fy_theo", "Fx_error", "Fy_error", "x_cmd", "y_cmd", "current1", "current2", "gz1", "gz2", "position1_sent", "position2_sent"])

        t_start = time.time()
        while time.time() - t_start < exercise_duration:
            dt_loop = 0.005
            theta1 = compute_theta1(x, y)
            theta4 = compute_theta4(x, y)
            theta1_sent = (degrees(theta1) - 13) / 2.03
            theta4_sent = (degrees(theta4) - 11) / 2.03

            print(f"Commande envoyée: theta1={degrees(theta1):.2f}, theta4={degrees(theta4):.2f}")

            # Dummy current read (replace if needed)
            current1, position1 = send_angle(servo1, theta1_sent)
            current2, position2 = send_angle(servo2, theta4_sent)
            print(f"Dernière commande: position1 = {(position1*2.03)+13:.4f}, position 2 = {(position2*2.03)+11:.4f}")
            if gz1 < 0: current2 = -current2
            if gz2 < 0: current1 = -current1
            tau1 = current_to_torque1(current1)
            tau2 = current_to_torque2(current2)

            # Calculate force and target
            Fx_real, Fy_real = dynamique(theta1, theta4, tau1, tau2)
            if np.sqrt(Fx_real**2 + Fy_real**2) < F_threshold:
                Fx_unitary, Fy_unitary = 0.0, 0.0
            else:
                Fx_unitary, Fy_unitary = unitary_vector(Fx_real, Fy_real)
            x_target, y_target = get_x_target(x, y, Fx_unitary, Fy_unitary)
            print(f"Fx_real={Fx_real:.2f}, Fy_real={Fy_real:.2f}")
            print(f"Nouvelle cible: x_target={x_target:.3f}, y_target={y_target:.3f}")

            # Logging
            now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            vx_real = (x - x_target) / dt_loop
            vy_real = (y - y_target) / dt_loop
            #Fx_theo = 0
            #Fy_theo = 0
            Fx_theo = B_virtual[0] * vx_real + K_virtual[0] * (x - x_target)
            Fy_theo = B_virtual[1] * vy_real + K_virtual[1] * (y - y_target)
            
            Fx_error = Fx_real - Fx_theo
            Fy_error = Fy_real - Fy_theo
            x_cmd = np.clip(x + (1 / De[0]) * Fx_error * dt_loop + Ki[0] * dt_loop * Fx_error, -0.2, 0.2)
            y_cmd = np.clip(y + (1 / De[1]) * Fy_error * dt_loop + Ki[1] * dt_loop * Fy_error, 0.23, 0.4)

            writer.writerow([now_str, x, y, vx_real, vy_real, Fx_real, Fy_real, Fx_theo, Fy_theo, Fx_error, Fy_error, x_cmd*1000, y_cmd, current1, current2, gz1, gz2,(position1*2.03)+13, (position2*2.03)+11])
            file.flush()

            x, y = x_cmd, y_cmd
            await asyncio.sleep(dt_loop)

    servo1.close()
    servo2.close()
    print("Contrôle terminé. CSV généré :", os.path.abspath(file_path))

# -------------------- Main async runner --------------------

async def main():
    imu1 = IMU_BLE("C3:98:10:C5:26:51", "IMU1") #Gauche
    imu2 = IMU_BLE("F7:E5:3A:2D:90:F9", "IMU2") #Droit
    await imu1.connect_and_start()
    await imu2.connect_and_start()
    print("IMU connectés. Démarrage du contrôle...")

    control_task = asyncio.create_task(run_control(imu1, imu2))
    await control_task

    await imu1.stop()
    await imu2.stop()

asyncio.run(main())
