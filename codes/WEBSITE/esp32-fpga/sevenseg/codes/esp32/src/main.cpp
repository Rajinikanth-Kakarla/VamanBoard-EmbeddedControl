#include <Arduino.h>
#include <WiFi.h>
#include <AsyncTCP.h>
#include <ESPAsyncWebServer.h>
#include "esp32_eoss3_spi.h"
#include "credentials.h"

#define GPIO_OUTPUT_VAL_REG 0x40021004
#define GPIO_OUTPUT_DIR_REG 0x40021008
#define PIN_a 1
#define PIN_b 2
#define PIN_c 3
#define PIN_d 4
#define PIN_e 5
#define PIN_f 6
#define PIN_g 7

AsyncWebServer server(80);

const char* PARAM_DIGIT = "digit";
const char index_html[] PROGMEM = R"rawliteral(
<!DOCTYPE HTML><html><head>
  <title>Vaman LED Form</title>
  <meta name="viewport" content="width=device-width, initial-scale=1">
  </head>
  <body>
  <h1>Seven-Segment Display</h1>
  <form action="/sevenseg" method="POST">
    Enter digit: <input type="number" min="0" max="9" step="1" name="state"><br>
    <input type="submit" value="Submit">
  </form>
</body></html>)rawliteral";

void notFound(AsyncWebServerRequest *request) {
  request->send(404, "text/plain", "Not found");
}

void display(int digit) { 
    int D = digit & (1<<3);
    int C = digit & (1<<2);
    int B = digit & (1<<1);
    int A = digit & (1<<0);
	int a = (C&&!B&&!A)||(!D&&!C&&!B&&A);
	int b = (C&&!B&&A)||(C&&B&&!A);
	int c = !C&&B&&!A;
	int d = (!D&&!C&&!B&&A)||(C&&!B&&!A)||(C&&B&&A);
	int e = (!B||A)&&(C||B||A);
	int f = (!D&&!C&&A)||(!C&&B)||(B&&A);
	int g = (!D&&!C&&!B)||(C&&B&&A);
    uint32_t gpioval = 0;
    gpioval |= (a<<PIN_a);
    gpioval |= (b<<PIN_b);
    gpioval |= (c<<PIN_c);
    gpioval |= (d<<PIN_d);
    gpioval |= (e<<PIN_e);
    gpioval |= (f<<PIN_f);
    gpioval |= (g<<PIN_g);
    esp32_eoss3_spi_ahb_write(GPIO_OUTPUT_VAL_REG, (uint8_t *)&gpioval, 4);
}

void setup() {
  esp32_eoss3_spi_init();
  uint32_t dirval = 0;
  dirval |= (1<<PIN_a);
  dirval |= (1<<PIN_b);
  dirval |= (1<<PIN_c);
  dirval |= (1<<PIN_d);
  dirval |= (1<<PIN_e);
  dirval |= (1<<PIN_f);
  dirval |= (1<<PIN_g);
  esp32_eoss3_spi_ahb_write(GPIO_OUTPUT_DIR_REG, (uint8_t *)&dirval, 4);
  Serial.begin(115200);
  WiFi.mode(WIFI_STA);
  WiFi.begin(ssid, password);
  while (WiFi.status() != WL_CONNECTED) {
      Serial.print(".");
      delay(2000);
  }
  Serial.println();
  Serial.print("IP Address: ");
  Serial.println(WiFi.localIP());

  server.on("/sevenseg", HTTP_GET, [](AsyncWebServerRequest *request){
    request->send_P(200, "text/html", index_html);
  });

  server.on("/sevenseg", HTTP_POST, [](AsyncWebServerRequest *request) {
    if (request->hasParam(PARAM_DIGIT)) {
        display(request->getParam(PARAM_DIGIT)->value().toInt());
    }
    request->send(200, "text/html", index_html);
  });

  server.onNotFound(notFound);
  server.begin();
}

void loop() {}
