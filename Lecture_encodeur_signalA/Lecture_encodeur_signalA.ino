#include <SPI.h>

// Pin Slave Select
const int slaveSelectPin = D7;

// Stop / Run
bool stopProgram = false;

void setup() {

  Serial.begin(9600);

  pinMode(slaveSelectPin, OUTPUT);
  digitalWrite(slaveSelectPin, HIGH);

  SPI.begin();
  Serial.println("SPI Ready sur XIAO nRF52840");
  Serial.println("s = stop | r = run");
}

void loop() {
  // Gestion stop/run (inchangée) [cite: 3, 4, 5, 6]
  if (Serial.available()) {
    char cmd = Serial.read();
    if (cmd == 's') stopProgram = true;
    if (cmd == 'r') stopProgram = false;
  }
  if (stopProgram) return;

  // Transaction SPI corrigée
  SPI.beginTransaction(SPISettings(100000, MSBFIRST, SPI_MODE0)); // Vitesse stable [cite: 7]
  digitalWrite(slaveSelectPin, LOW);

  // On lit 2 octets sans envoyer de commande 0x80
  byte byte1 = SPI.transfer(0x00); 
  byte byte2 = SPI.transfer(0x00); 

  digitalWrite(slaveSelectPin, HIGH);
  SPI.endTransaction();

  // Extraction du statut (2 premiers bits)
  byte status = (byte1 & 0xC0) >> 6;
  
  // Extraction de la pression brute (14 bits)
  int rawData = ((byte1 & 0x3F) << 8) | byte2;

  // Affichage [cite: 9, 10]
  Serial.print("Statut: "); Serial.print(status);
  Serial.print(" | Pression Brute: "); Serial.println(rawData);

  delay(500); 
}