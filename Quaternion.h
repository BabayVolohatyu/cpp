#ifndef QUATERNION_H
#define QUATERNION_H

#include <iostream>
#include <cmath>

#include "Vector3.h"

class Quaternion {
private:
    double _w, _x, _y, _z;

public:
    Quaternion()
        : _w(1.0), _x(0.0), _y(0.0), _z(0.0) {
    }

    Quaternion(double w, const Vector3 &rotation)
        : Quaternion(1, rotation.x(), rotation.y(), rotation.z()) {
    }

    Quaternion(const Vector3 &rotation)
        : Quaternion(rotation.x(), rotation.y(), rotation.z()) {
    }

    Quaternion(double w, double x, double y, double z)
        : _w(w), _x(x), _y(y), _z(z) {
    }

    Quaternion(double xRotation, double yRotation, double zRotation) {
        double phi = degreesToRadians(xRotation);
        double theta = degreesToRadians(yRotation);
        double psi = degreesToRadians(zRotation);
        _w = std::cos(phi / 2) * std::cos(theta / 2) * std::cos(psi / 2) +
             std::sin(phi / 2) * std::sin(theta / 2) * std::sin(psi / 2);

        _x = std::sin(phi / 2) * std::cos(theta / 2) * std::cos(psi / 2)
             - std::cos(phi / 2) * std::sin(theta / 2) * std::sin(psi / 2);


        _y = std::cos(phi / 2) * std::sin(theta / 2) * std::cos(psi / 2)
             + std::sin(phi / 2) * std::cos(theta / 2) * std::sin(psi / 2);

        _z = std::cos(phi / 2) * std::cos(theta / 2) * std::sin(psi / 2)
             - std::sin(phi / 2) * std::sin(theta / 2) * std::cos(psi / 2);
    }

    ~Quaternion() = default;

    static double degreesToRadians(double degrees) {
        return degrees * M_PI / 180;
    }

    void normalize() {
        double mag = magnitude();
        _w /= mag;
        _x /= mag;
        _y /= mag;
        _z /= mag;
    }

    Quaternion conjugate() const {
        return Quaternion(_w, -_x, -_y, -_z);
    }

    double magnitude() const {
        return std::sqrt(_w * _w + _x * _x + _y * _y + _z * _z);
    }

    Quaternion operator*(const Quaternion &other) const noexcept {
        const double newW = _w * other._w - _x * other._x - _y * other._y - _z * other._z;
        const double newX = _w * other._x + _x * other._w + _y * other._z - _z * other._y;
        const double newY = _w * other._y - _x * other._z + _y * other._w + _z * other._x;
        const double newZ = _w * other._z + _x * other._y - _y * other._x + _z * other._w;
        return Quaternion(newW, newX, newY, newZ);
    }
    Vector3 operator*(const Vector3 &v) const {
        Quaternion q_v(0, v.x(), v.y(), v.z());
        Quaternion q_conj = conjugate();
        Quaternion q_rotated = *this * q_v * q_conj;
        return Vector3(q_rotated.x(), q_rotated.y(), q_rotated.z());
    }
    Quaternion &operator*=(const Quaternion &other) noexcept {
        const double oldW = _w;
        const double oldX = _x;
        const double oldY = _y;
        const double oldZ = _z;

        _w = oldW * other._w - oldX * other._x - oldY * other._y - oldZ * other._z;
        _x = oldW * other._x + oldX * other._w + oldY * other._z - oldZ * other._y;
        _y = oldW * other._y - oldX * other._z + oldY * other._w + oldZ * other._x;
        _z = oldW * other._z + oldX * other._y - oldY * other._x + oldZ * other._w;

        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const Quaternion &quaternion) {
        os << quaternion._w << ' ' << quaternion._x << ' ' << quaternion._y << ' ' << quaternion._z;
        return os;
    }

    friend std::istream &operator>>(std::istream &is, Quaternion &quaternion) = delete;

    double w() { return _w; }
    double x() { return _x; }
    double y() { return _y; }
    double z() { return _z; }

    void w(double w) { _w = w; }
    void x(double x) { _x = x; }
    void y(double y) { _y = y; }
    void z(double z) { _z = z; }
};

#endif
