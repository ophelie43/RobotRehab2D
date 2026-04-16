from machine import Pin, PWM
import time
import sys
import uselect
import plasma
from servo import Servo, servo2040
from pimorini import Analog, AnalogMux

# --- Configuration ---
# Si le servo de COM5 est branché sur le port 2, changez la valeur ci-dessous par 2.
BROCHE_SERVO = 2
servo = PWM(Pin(BROCHE_SERVO))
servo.freq(50)

# ==========================================
# CONFIGURATION DES LEDS (MÉTHODE INFAILLIBLE)
# ==========================================
BROCHE_LEDS = 18  
NUM_LEDS = 6

leds = plasma.WS2812(NUM_LEDS, 0, 0, BROCHE_LEDS)
leds.start()

# --- Configuration Encodeur ---
pin_a = Pin(26, Pin.IN)
pin_b = Pin(27, Pin.IN)
pin_z = Pin(28, Pin.IN)

compteur = 0
derniere_val_a = pin_a.value()

def handle_encoder(pin):
    global compteur, derniere_val_a
    val_a = pin_a.value()
    val_b = pin_b.value()
    if val_a != derniere_val_a:
        if val_b != val_a:
            compteur -= 1  
        else:
            compteur += 1  
    derniere_val_a = val_a

def handle_index(pin):
    global compteur
    compteur = 0 
    print("MSG: Passage au point Zéro (Index)")

pin_a.irq(trigger=Pin.IRQ_RISING | Pin.IRQ_FALLING, handler=handle_encoder)
pin_z.irq(trigger=Pin.IRQ_RISING, handler=handle_index)

# --- 2. Configurer Courant -----
current_adc = Analog(servo2040.SHARED_ADC, servo2040. CURRENT_GAIN,
                      servo2040.SHUNT_RESISTOR, servo2040.CURRENT_OFFSET)
mux = AnalogMux(servo2040.ADC_ADDR_0, servo2040.ADC_ADDR_1, servo2040.ADC_ADR_2, 
                muxed_pin = Pin(servo2040.SHARED_ADC))
mux.select(servo2040.CURRENT_SENSE_ADDR)



def aller_a_angle(moteur, angle):
    angle = max(0, min(180, angle))
    angle_calibre = angle / 2.03 
    min_duty = 1638
    max_duty = 8192
    duty = int(min_duty + (angle_calibre / 180) * (max_duty - min_duty))
    moteur.duty_u16(duty)
    print(f"CONFIRM:Moteur a {angle}")

# Séquence d'initialisation
print("MSG: Démarrage de la séquence d'initialisation...")
# ALLUMAGE DES LEDS EN VERT !
for i in range(NUM_LEDS):
    leds.set_rgb(i, 0, 50, 0) 

leds.update() # On pousse l'information vers les LEDs
aller_a_angle(servo, 135) # Vous pouvez ajuster cet angle de départ (ex: 135 pour l'autre)
time.sleep(1)
print("MSG: Initialisation terminée !")

# Configuration de l'écoute USB
poll_obj = uselect.poll()
poll_obj.register(sys.stdin, uselect.POLLIN)

while True:


    # Si l'ordinateur envoie un message
    if poll_obj.poll(0):
        commande = sys.stdin.readline().strip()
        # On attend une commande sous la forme "S:90"
        if commande.startswith("S:"):
            try:
                angle_cible = float(commande.split(":")[1])
                aller_a_angle(servo, angle_cible)
            except ValueError:
                print("ERREUR: Commande invalide")

    # Envoi continu de la position de l'encodeur à l'ordinateur
        current = sum([current_adc.read_current() for _ in range(10)]) / 10
    angle_encodeur = (compteur / (300 * 2)) * 360 
    print(f"ENCODEUR:{angle_encodeur:.1f} | COURANT: {current:.4f}")
    time.sleep(0.05)
    
    

leds.start()
