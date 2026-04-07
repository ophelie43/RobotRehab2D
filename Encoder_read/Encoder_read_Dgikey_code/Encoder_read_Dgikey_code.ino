#include <SPI.h>

// On définit la pin SS (Slave Select) sur la pin 10
const int slaveSelectPin = 10;

void setup() {
  // Initialisation du port série pour voir les résultats sur l'ordinateur
  Serial.begin(9600);

  // Configuration de la pin SS en sortie
  pinMode(slaveSelectPin, OUTPUT);
  digitalWrite(slaveSelectPin, HIGH);

  // Initialisation du bus SPI
  SPI.begin();
  
  // Configuration des paramètres SPI (vitesse, ordre des bits, mode)
  // Ces paramètres dépendent de la fiche technique de ton capteur
  SPI.beginTransaction(SPISettings(1000000, MSBFIRST, SPI_MODE0));
  
  Serial.println("Initialisation SPI terminee.");
}

void loop() {
  byte reponse;

  // 1. On active le capteur (SS bas)
  digitalWrite(slaveSelectPin, LOW);

  // 2. On envoie une commande (ex: 0x80 pour lire un registre)
  // et on récupère la réponse simultanément
  reponse = SPI.transfer(0x80); 

  // 3. On désactive le capteur (SS haut)
  digitalWrite(slaveSelectPin, HIGH);

  // Affichage du résultat
  Serial.print("Valeur lue du capteur : ");
  Serial.println(reponse, HEX);

  delay(1000); // Pause d'une seconde entre les lectures
}