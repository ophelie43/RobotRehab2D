import serial
import time
import matplotlib.pyplot as plt

SERVO_PORT = 'COM4'
servo = serial.Serial(SERVO_PORT, 115200, timeout=1)
time.sleep(2)  # Laisse le temps à la carte de s'initialiser

# Active le mode STREAM une seule fois
def read_stream(servo):
    servo.write(b"STREAM\n")
    while True:
        # Lit une seule ligne de données depuis le servo
        # Sans boucle ni tentative multiple
        try:
            line = servo.readline().decode(errors="replace").strip()
            if not line:
                print("⚠️ Aucune donnée reçue.")
                continue

            current_str, position_str = line.split(",")
            current = float(current_str.strip())
            position = float(position_str.strip())
            print(f"Courant: {current:.3f} A, Position: {position:.2f}°")

        except Exception as e:
            print(f"⚠️ Erreur lecture servo : {e}")
        
        time.sleep(0.001)  # Pause pour éviter de surcharger le port série

# Pour stocker les positions
now = time.time()
servo.write(f"SET {0}\n".encode())
time.sleep(0.2)
while time.time() - now < 10:  # Lit pendant 10 secondes
    
    read_stream(servo)
servo.write(b"STOP\n")  # Arrête le mode STREAM
