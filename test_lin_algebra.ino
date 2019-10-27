#include <LinearAlgebra.h>

// handle diagnostic informations given by assertion and abort program execution:
void __assert(const char *__func, const char *__file, int __lineno, const char *__sexp) {
    // transmit diagnostic informations through serial link. 
    Serial.println(__func);
    Serial.println(__file);
    Serial.println(__lineno, DEC);
    Serial.println(__sexp);
    Serial.flush();
    // abort program execution.
    abort();
}



void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);

  vec v = vec::ones(6);


  for(uint8_t i=0; i<6;++i) {
      v(i) = random(0,100)/10.0d;
  }
  String str;
 
  mat m = mat::zeros(6,6);
  for(uint8_t i=0; i<6;++i) {
    for (uint8_t j=0;j<6;++j) {
      m(i,j) = random(0,1000)/100.0d;
    }
  }

  mat n = mat::zeros(6,6);
  for(uint8_t i=0; i<6;++i) {
    for (uint8_t j=0;j<6;++j) {
      n(i,j) = random(0,1000)/100.0d;
    }
  }

  //m.print(str);
  long mil = millis();
  mat xp = (m * m) + (n * n);
    Serial.println(millis() - mil);
  xp.print(str);
  Serial.println(str);

  mat xp2 = m * (m + n) * n;
  xp2.print(str);
  Serial.println(str);  
  return;
  //long mil = millis();
  m.inv().print(str);
  Serial.println(str);

  m.inv2().print(str);
  Serial.println(str);
  mat diff = m.inv() - m.inv2();
  diff.print(str);
  Serial.println(str);
  /*
  long mil = millis();
  m.inv();
  Serial.println(millis() - mil);
  mil = millis();
  m.inv2();
  Serial.println(millis() - mil);
  
  //long diff = millis() - mil;
  //Serial.println(diff);
  /*
  //mat n = m - (m.inv().inv());
  double det = m.det();
  Serial.println(det);
  
  m.print(str);
  Serial.println(str);
  n.print(str);
  Serial.println(str);
  */
}

void loop() {
  // put your main code here, to run repeatedly:

}
