from machine import Pin
import time

# On force la pin GP28 (A2) à 0V au repos avec PULL_DOWN
pin_z = Pin(28, Pin.IN, Pin.PULL_DOWN)

index_trouve = False

def alerte_index(pin):
    global index_trouve
    index_trouve = True # On mémorise le passage

# On surveille la montée du signal (RISING)
pin_z.irq(trigger=Pin.IRQ_RISING, handler=alerte_index)

print("Recherche de l'index... Tourne l'axe TRES lentement sur 360°")

while True:
    if index_trouve:
        print("!!! INDEX CAPTURÉ PAR LE MATÉRIEL !!!")
        index_trouve = False # Reset pour le tour suivant
    
    time.sleep(0.01) # On tourne vite pour ne pas saturer la console