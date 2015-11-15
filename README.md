MICO-RE
===
MICO-RE 【miko-li】  
Minimal Implementation of Cloud Optical Retrieval

MICO-RE estimates cloud optical thickness (COT) and
cloud droplet effective radius (CDER) from reflectances in two wavelengths
by using look-up table (Nakajima and King, 1990) and Gauss-Newton method.

Cloud Optical Retrieval Code
---



Look-up table
---
Look-up table used in this program should be like following:

ex.)

| COT | CDER | R1   | R2   |
|:----|:-----|:-----|:-----|
| 1.0 | 3.0  | 0.20 | 0.10 |
| 1.0 | 5.0  | 0.20 | 0.13 |
| 1.0 | 7.0  | 0.21 | 0.15 |
| 1.0 | 9.0  | 0.23 | 0.14 |
| ... | ...  | ...  | ...  |
| ... | ...  | ...  | ...  |
| 2.0 | 3.0  | 0.26 | 0.14 |
| 2.0 | 5.0  | 0.26 | 0.18 |
| 2.0 | 7.0  | 0.27 | 0.19 |
| 2.0 | 9.0  | 0.29 | 0.18 |
| ... | ...  | ...  | ...  |
| ... | ...  | ...  | ...  |


