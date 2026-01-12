# === servo_hybride.py ===
from servo import Servo, servo2040
from pimoroni import Analog, AnalogMux
from machine import Pin, ADC
import time
import sys
import select

SAMPLES = 20
TIME_BETWEEN = 0.001

# Initialisation du servo
servo = Servo(servo2040.SERVO_1)
servo.enable()
angle = 0
servo.value(angle / 2.03 - 8)

# Position via potentiomètre
position_adc = ADC(Pin(26))

# Courant via ADC partagé
current_adc = Analog(servo2040.SHARED_ADC,
                     servo2040.CURRENT_GAIN,
                     servo2040.SHUNT_RESISTOR,
                     servo2040.CURRENT_OFFSET)
mux = AnalogMux(servo2040.ADC_ADDR_0, servo2040.ADC_ADDR_1, servo2040.ADC_ADDR_2,
                muxed_pin=Pin(servo2040.SHARED_ADC))
mux.select(servo2040.CURRENT_SENSE_ADDR)
time.sleep(0.002)

def read_data():
    current = 0
    for _ in range(SAMPLES):
        current += current_adc.read_current()
        time.sleep(TIME_BETWEEN)
    current /= SAMPLES

    raw_adc = position_adc.read_u16()
    position = (raw_adc - 18219) / 102.76
    print(f"{current - 0.016:.4f}, {position:.4f}")

# Mode de streaming
streaming = False
last_stream_time = time.ticks_ms()

while True:
    # Vérifie si une commande a été reçue depuis le PC
    if sys.stdin in select.select([sys.stdin], [], [], 0)[0]:
        line = sys.stdin.readline()
        if not line:
            continue
        parts = line.strip().split()
        if not parts:
            continue

        cmd = parts[0].upper()

        if cmd == "SET" and len(parts) == 2:
            try:
                angle = float(parts[1])
                servo.value(angle / 2.03 - 8)
            except:
                print("err")
                sys.stdout.flush()

        elif cmd == "GET":
            read_data()

        elif cmd == "STREAM":

            streaming = True

        elif cmd == "STOP":
            streaming = False

    # Streaming en continu toutes les 50 ms
    if streaming and time.ticks_diff(time.ticks_ms(), last_stream_time) > 50:
        read_data()
        last_stream_time = time.ticks_ms()
