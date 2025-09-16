#ifndef CITY_H
#define CITY_H


//contiene le coordinate della city
struct City {
    double x;
    double y;
    
    City(double x = 0, double y = 0) : x(x), y(y) {}
};

#endif
