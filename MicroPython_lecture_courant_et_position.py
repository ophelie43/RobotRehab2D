# === servo_hybride.py ===
from servo import Servo, servo2040
from pimoroni import Analog, AnalogMux
from machine import Pin, ADC
import time
import math
import sys
import select

SAMPLES = 20
TIME_BETWEEN = 0.001

# Initialisation du servo
servo = Servo(servo2040.SERVO_1)
servo.enable()
# Lecture du courant

current_adc = Analog(servo2040.SHARED_ADC,
                     servo2040.CURRENT_GAIN,
                     servo2040.SHUNT_RESISTOR,
                     servo2040.CURRENT_OFFSET)
mux = AnalogMux(servo2040.ADC_ADDR_0, servo2040.ADC_ADDR_1, servo2040.ADC_ADDR_2,
                muxed_pin=Pin(servo2040.SHARED_ADC))
mux.select(servo2040.CURRENT_SENSE_ADDR)
time.sleep(0.002)
position_adc = ADC(Pin(26))

def read_data():
    # Moyenner le courant
    current = current_adc.read_current()
    #time.sleep(TIME_BETWEEN)

    raw_adc = position_adc.read_u16()
    position = (raw_adc - 18219) / 102.76
    print(f"{current - 0.016:.4f}, {position:.4f}")

SWEEPS = 3              # How many sweeps of the servo to perform
STEPS = 10              # The number of discrete sweep steps
STEPS_INTERVAL = 0.5    # The time in seconds between each step of the sequence
SWEEP_EXTENT = 90.0     # How far from zero to move the servo when sweeping

# Do a sine sweep
for _j in range(SWEEPS):
    for i in range(360):
        servo.value(math.sin(math.radians(i)) * SWEEP_EXTENT)
        read_data()
        time.sleep(0.02)


# Disable the servo
s.disable()



# Mode de streaming
streaming = False
last_stream_time = time.ticks_ms()