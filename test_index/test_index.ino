const int indexPin = A0;
void setup() {
  // put your setup code here, to run once:
  Serial.begin(115200); // Plus rapide pour ne pas ralentir les interruptions
  pinMode(indexPin, INPUT);
}

void loop() {
  // Lecture de la tension brute (0 à 1023)
  int valAnalogique = analogRead(indexPin);
  
  // Lecture de l'interprétation numérique de l'Arduino (0 ou 1)
  // On multiplie par 1000 pour que les deux courbes soient sur la même échelle
  int valNumerique = digitalRead(indexPin) * 1000;

  // Envoi au traceur série
  Serial.print("Brut_0_a_1023:");
  Serial.print(valAnalogique);
  Serial.print(" "); // Séparateur pour le traceur
  Serial.print("Interpretation_Numerique:");
  Serial.println(valNumerique);

  // Pas de delay() pour ne rater aucun pic de bruit !
}
