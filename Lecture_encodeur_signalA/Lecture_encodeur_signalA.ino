void setup() {
  // Une vitesse plus haute (115200) est préférable pour le plotter
  Serial.begin(115200); 
}

void loop() {
  // Lecture de la valeur analogique (0 à 1023) sur A0
  int rawValue = analogRead(A0);
  
  // Conversion en tension réelle (V)
  float voltage = rawValue * (5.0 / 1023.0);

  // --- Affichage pour le Serial Plotter ---
  // On affiche des "bornes" fixes pour stabiliser le graphique
  Serial.print(0.0);    // Limite basse
  Serial.print(" ");
  Serial.print(3.0);    // Limite haute
  Serial.print(" ");
  
  // La valeur mesurée
  Serial.println(voltage);

  // Un délai très court pour voir les transitions en "temps réel"
  delay(5); 
}