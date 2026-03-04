from machine import Pin
import time

# Configuration des pins SENSORS
pin_a = Pin(26, Pin.IN) # Phase A
pin_b = Pin(27, Pin.IN) # Phase B
pin_z = Pin(28, Pin.IN) # Index Z

compteur = 0
derniere_val_a = pin_a.value()
tours_complets = 0

# Interruption pour le comptage (Phase A)
def handle_encoder(pin):
    global compteur, derniere_val_a
    val_a = pin_a.value()
    val_b = pin_b.value()
    if val_a != derniere_val_a:
        if val_b != val_a:
            compteur += 1
        else:
            compteur -= 1
    derniere_val_a = val_a

# Interruption pour l'Index (Z)
def handle_index(pin):
    global compteur, tours_complets
    # Option 1 : Réinitialiser le compteur à zéro pour une référence absolue
    compteur = 0
    tours_complets += 1
    print("--- Passage au point Zéro (Index) ---")

# Activation des interruptions
pin_a.irq(trigger=Pin.IRQ_RISING | Pin.IRQ_FALLING, handler=handle_encoder)
pin_z.irq(trigger=Pin.IRQ_RISING, handler=handle_index)

print("Système avec Index prêt. Tournez l'encodeur...")

while True:
    angle = (compteur / (300 * 2)) * 360 # Ajuste selon tes tests de précision
    print(f"Angle : {angle:.1f}° | Tours : {tours_complets}")
    time.sleep(0.2)