import time
from machine import ADC, Pin

# A0=GP26, A1=GP27, A2=GP28
phase_a = ADC(Pin(26))
phase_b = ADC(Pin(27))
index_z = ADC(Pin(28))

print("Mode diagnostic : Tournez l'encodeur très lentement...")

while True:
    val_a = (phase_a.read_u16() / 65535) * 3.3
    val_b = (phase_b.read_u16() / 65535) * 3.3
    val_z = (index_z.read_u16() / 65535) * 3.3
    
    # Format optimisé pour le Traceur de Thonny
    
    print(f"A:{val_a:.2f} B:{val_b:.2f} Index:{val_z:.2f}")
    
    # On réduit le délai au minimum pour ne pas rater le pic Z
    time.sleep(0.01)