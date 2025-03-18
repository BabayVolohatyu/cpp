#ifndef VECTOR3_H
#define VECTOR3_H

class Vector3 {
private:
    double _x, _y, _z;

public:
    Vector3(): _x(0.0), _y(0.0), _z(0.0) {
    }

    Vector3(double x, double y, double z): _x(x), _y(y), _z(z) {
    }

    Vector3 operator+(const Vector3 &other) const {
        return Vector3(_x + other._x, _y + other._y, _z + other._z);
    }

    Vector3 operator-(const Vector3 &other) const {
        return Vector3(_x - other._x, _y - other._y, _z - other._z);
    }

    Vector3 operator*(double scalar) const {
        return Vector3(_x * scalar, _y * scalar, _z * scalar);
    }

    Vector3 operator/(double scalar) const {
        return Vector3(_x / scalar, _y / scalar, _z / scalar);
    }

    void normalize() {
        double mag = magnitude();
        _x /= mag;
        _y /= mag;
        _z /= mag;
    }

    double magnitude() const {
        return std::sqrt(_x * _x + _y * _y + _z * _z);
    }

    double x() const { return _x; }

    double y() const { return _y; }

    double z() const { return _z; }

    void x(double x) { _x = x; }

    void y(double y) { _y = y; }

    void z(double z) { _z = z; }
};
#endif
