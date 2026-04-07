int counter = 0;
String dir = "";
unsigned long last_run = 0;
void setup() {
  Serial.begin(9600);
  attachInterrupt(digitalPinToInterrupt(3), shaft_moved, FALLING);
  pinMode(4, INPUT);
}

void loop() {
  Serial.print("counterP : ");
  Serial.print(counter);
  Serial.print( "direction : ");
  Serial.println(dir);
}
void shaft_moved(){
  if (digitalRead(4) == 1) {
    counter++;
    dir = "CW";
  }
  if (digitalRead(4) == 0) {
    counter--;
    dir = "CCW";
  }
  last_run = millis();
}
