#include <Arduino.h>
#ifdef ESP32
  #include <WiFi.h>
  #include <AsyncTCP.h>
#else
  #include <ESP8266WiFi.h>
  #include <ESPAsyncTCP.h>
#endif
#include <ESPAsyncWebServer.h>
#include <math.h>
#include "matfun.h"
#include "esp32_eoss3_spi.h"

#define GPIO_OUTPUT_VAL_REG 0x3FF44004
#define GPIO_OUTPUT_DIR_REG 0x3FF44008
#define PIN_BLUE 18
#define PIN_GREEN 21
#define PIN_RED 22
#define PIN_ALL ((1<<PIN_GREEN) | (1<<PIN_BLUE) | (1<<PIN_RED))

AsyncWebServer server(80);

const char* ssid = "Rajinikanth";
const char* password = "Beast@03";

const char* PARAM_AX = "ax";
const char* PARAM_AY = "ay";
const char* PARAM_BX = "bx";
const char* PARAM_BY = "by";
const char* PARAM_CX = "cx";
const char* PARAM_CY = "cy";

const char index_html[] PROGMEM = R"rawliteral(
<!DOCTYPE HTML><html><head>
  <title>Angle Bisector</title>
  <meta name="viewport" content="width=device-width, initial-scale=1">
  </head>
  <body>
  <h1>Angle Bisecctor</h1>
  <form action="/anglebisector">
    <label for="ax">Vertex A - x:</label>
    <input type="number" id="ax" name="ax" required>
    <label for="ay">Vertex A - y:</label>
    <input type="number" id="ay" name="ay" required>
    <label for="bx">Vertex B - x:</label>
    <input type="number" id="bx" name="bx" required>
    <label for="by">Vertex B - y:</label>
    <input type="number" id="by" name="by" required>
    <label for="cx">Vertex C - x:</label>
    <input type="number" id="cx" name="cx" required>
    <label for="cy">Vertex C - y:</label>
    <input type="number" id="cy" name="cy" required>
    <input type="submit" value="Calculate">
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
  
  while (WiFi.waitForConnectResult() != WL_CONNECTED) {
    Serial.println("Connecting to WiFi...");
    delay(1000);
  }

  Serial.println("Connected to WiFi");
  Serial.println("IP Address: ");
  Serial.println(WiFi.localIP());

  server.on("/", HTTP_GET, [](AsyncWebServerRequest *request) {
    request->send_P(200, "text/html", index_html);
  });

  server.on("/anglebisector", HTTP_GET, [](AsyncWebServerRequest *request) {
    int params = request->params();
    float ax = 0, ay = 0, bx = 0, by = 0, cx = 0, cy = 0;

    for (int i = 0; i < params; i++) {
      AsyncWebParameter* p = request->getParam(i);
      if (p->isPost()) {
        if (p->name() == PARAM_AX) {
          ax = p->value().toFloat();
        } else if (p->name() == PARAM_AY) {
          ay = p->value().toFloat();
        } else if (p->name() == PARAM_BX) {
          bx = p->value().toFloat();
        } else if (p->name() == PARAM_BY) {
          by = p->value().toFloat();
        } else if (p->name() == PARAM_CX) {
          cx = p->value().toFloat();
        } else if (p->name() == PARAM_CY) {
          cy = p->value().toFloat();
        }
      }
    }

    double **A, **B, **C, **eigen, a, b, c, **D, **E, **F;
    int m1 = 2, n1 = 1;
    A = createMat(m1, n1);
    B = createMat(m1, n1);
    C = createMat(m1, n1);
    eigen = createMat(m1, n1);

    if (A && B && C && eigen) {
      A[0][0] = ax;
      A[1][0] = ay;
      B[0][0] = bx;
      B[1][0] = by;
      C[0][0] = cx;
      C[1][0] = cy;

      double **diff_AB = Matsub(B, A, 2, 1);
      double distance_AB = Matnorm(diff_AB, 2);

      double **diff_AC = Matsub(C, A, 2, 1);
      double distance_AC = Matnorm(diff_AC, 2);

      double **diff_BC = Matsub(C, B, 2, 1);
      double distance_BC = Matnorm(diff_BC, 2);

      a = distance_BC;
      b = distance_AC;
      c = distance_AB;

      double l1 = (a + c - b) / 2;
      double l2 = (a + b - c) / 2;
      double l3 = (c + b - a) / 2;

      D = Matscale(Matadd(Matscale(C, 2, 1, l1), Matscale(B, 2, 1, l2), 2, 1), 2, 1, 1.0 / (l1 + l2));
      E = Matscale(Matadd(Matscale(A, 2, 1, l2), Matscale(C, 2, 1, l3), 2, 1), 2, 1, 1.0 / (l2 + l3));
      F = Matscale(Matadd(Matscale(B, 2, 1, l3), Matscale(A, 2, 1, l1), 2, 1), 2, 1, 1.0 / (l3 + l1));

      double **temp1 = Matscale(A, 2, 1, a);
      double **temp2 = Matscale(B, 2, 1, b);
      double **temp3 = Matscale(C, 2, 1, c);

      double **eigenvalues = Matadd(Matadd(temp1, temp2, 2, 1), temp3, 2, 1);
      double eigendenominator = a + b + c;

      eigen[0][0] = eigenvalues[0][0] / eigendenominator;
      eigen[1][0] = eigenvalues[1][0] / eigendenominator;

      String result =
          "D: " + String(D[0][0], 2) + ", " + String(D[1][0], 2) + "<br>"
          "E: " + String(E[0][0], 2) + ", " + String(E[1][0], 2) + "<br>"
          "F: " + String(F[0][0], 2) + ", " + String(F[1][0], 2) + "<br>"
          "distance of AB: " + String(c) + "<br>"
          "distance of BC: " + String(a) + "<br>"
          "distance of CA: " + String(b) + "<br>"
          "InCenter:" + String(eigen[0][0], 2) + ", " + String(eigen[1][0], 2) + "<br>";
      request->send(200, "text/html", result);
    } else {
      request->send(500, "text/plain", "Error allocating memory");
    }
    
    free(A);
    free(B);
    free(C);
    free(eigen);
    free(D);
    free(E);
    free(F);
  });

  server.onNotFound(notFound);

  server.begin();
}

void loop() {
  // Handle other tasks in the loop if needed
}

