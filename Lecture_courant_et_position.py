import serial
import time

SERVO_PORT = '/dev/tty.usbmodem101' 

try:
    # Utilisation de pyserial
    servo = serial.Serial(SERVO_PORT, 115200, timeout=1)
    print(f"Connecté avec succès à {SERVO_PORT}")
except Exception as e:
    print(f"Erreur de connexion : {e}")
    exit()

time.sleep(2) 

# On active le flux de données
servo.write(b"STREAM\n")
print("Lecture commencée...")

start_time = time.time()

try:
    while time.time() - start_time < 10:  # Boucle de 10 secondes
        line = servo.readline().decode(errors="replace").strip()
        
        if line:
            try:
                # On sépare les données (format attendu : "courant, position")
                parts = line.split(",")
                if len(parts) == 2:
                    current = float(parts[0].strip())
                    position = float(parts[1].strip())
                    print(f"Courant: {current:.3f} A, Position: {position:.2f}°")
            except ValueError:
                # Ignore les lignes mal formées au démarrage
                continue
        
        time.sleep(0.01)

finally:
    # S'exécute même s'il y a une erreur pour arrêter la carte proprement
    servo.write(b"STOP\n")
    servo.close()
    print("Communication terminée.")