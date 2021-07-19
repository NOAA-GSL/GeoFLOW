//==================================================================================
// Module       : gphysics.h
// Date         : 6/15/20 (DLR)
// Description  : Basic physics constants
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GPHYSICS_H)
#define _GPHYSICS_H

//#include "configure.hpp"

  #define CPD     1004.0    // sp. heat const press, dry air J/kg-K
  #define CPV     1885.0    // sp. heat const press, liq. water J/Kg-K
  #define CVD      717.0    // sp. heat const vol, dry air J/Kg-K
  #define CVV     1424.0    // sp. heat const vol, water vapor J/Kg-K
  #define CVI     2108.0    // sp. heat const vol, ice J/Kg-K
  #define CVL     4186.0    // sp. heat const press, liq. water J/Kg-K
  #define GG         9.81   // accel due to gravity m/s^2
  #define LV0        2.50e6 // latent heat of vaporization, reference J/kg
  #define RD       287.058  // gas const for dry air J/kg-K
  #define RV       461.0    // gas const for water vapor J/kg-K
  #define TKREF    273.15   // ref temp K
  #define KBOLTZ   1.38064852e-23 // Boltznann constat, kg m^2/s^2-K
#endif // !defined(_GPHYSICS_H)

