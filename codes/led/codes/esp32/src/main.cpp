#include <Arduino.h>
#include <WiFi.h>
#include <AsyncTCP.h>
#include <ESPAsyncWebServer.h>
#include "esp32_eoss3_spi.h"
#include "credentials.h"

#define GPIO_OUTPUT_VAL_REG 0x40021004
#define GPIO_OUTPUT_DIR_REG 0x40021008
#define PIN_BLUE 18
#define PIN_GREEN 21
#define PIN_RED 22
#define PIN_ALL (1<<PIN_GREEN) | (1<<PIN_BLUE) | (1<<PIN_RED)

AsyncWebServer server(80);

const char* PARAM_STATE = "state";
const char index_html[] PROGMEM = R"rawliteral(
<!DOCTYPE HTML><html><head>
  <title>Vaman LED Form</title>
  <meta name="viewport" content="width=device-width, initial-scale=1">
  </head>
  <body>
  <h1>LED State</h1>
  <form action="/led" method="POST">
    <input type="radio" value="ON" name="state">
    <label for="ON">ON</label>
    <input type="radio" value="OFF" name="state">
    <label for="OFF">OFF</label>
    <input type="submit" value="Submit">
  </form>
</body></html>)rawliteral";

void notFound(AsyncWebServerRequest *request) {
  request->send(404, "text/plain", "Not found");
}

void setup() {
  esp32_eoss3_spi_init();
  uint32_t dirval = (1<<PIN_GREEN) | (1<<PIN_BLUE) | (1<<PIN_RED);
  uint32_t gpioval = 0;
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

  server.on("/led", HTTP_GET, [](AsyncWebServerRequest *request){
    request->send_P(200, "text/html", index_html);
  });

  server.on("/led", HTTP_POST, [](AsyncWebServerRequest *request) {
    int params = request->params();
    for(int i=0;i<params;i++){
      AsyncWebParameter* p = request->getParam(i);
      if(p->isPost()){
        if (p->name() == PARAM_STATE) {
          if (p->value() == "ON") gpioval |= (PIN_ALL);
          else if (p->value() == "OFF") gpioval &= ~(PIN_ALL);
          Serial.println(gpioval, HEX);
          esp32_eoss3_spi_ahb_write(GPIO_OUTPUT_VAL_REG, (uint8_t *)&gpioval, 4);
        }
      }
    }
    request->send(200, "text/html", index_html);
  });

  server.onNotFound(notFound);
  server.begin();
}

void loop() {}
