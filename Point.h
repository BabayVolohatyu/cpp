#ifndef POINT_H
#define POINT_H

#include "Vector3.h"

class Point {
private:
    int _id;
    static inline int idCounter = 0;
    Vector3 _position;

public:
    Point()
        : _position(0, 0, 0) {
        _id = idCounter++;
    }

    Point(const Vector3 &position)
        : _position(position) {
        _id = idCounter++;
    }

    Point(double x, double y, double z) {
        _position = Vector3(x, y, z);
        _id = idCounter++;
    }

    Vector3 position() const { return _position; }

    double x() const { return _position.x(); }

    double y() const { return _position.y(); }

    double z() const { return _position.z(); }

    void position(const Vector3 &position) { _position = position; }

    void x(double x) { _position.x(x); }

    void y(double y) { _position.y(y); }

    void z(double z) { _position.z(z); }
};

#endif
