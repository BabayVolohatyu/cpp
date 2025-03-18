#ifndef QUATERNION_H
#define QUATERNION_H

#include <iostream>
#include <cmath>

#include "Vector3.h"
#include "Matrix.h"

class Quaternion {
private:
    double _w, _x, _y, _z;

    static void approximate(double &w, double &x, double &y, double &z) {
        if (std::abs(w) <= 1e-10) w = 0;
        if (std::abs(x) <= 1e-10) x = 0;
        if (std::abs(y) <= 1e-10) y = 0;
        if (std::abs(z) <= 1e-10) z = 0;
    }

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

        approximate(_w, _x, _y, _z);
    }

    ~Quaternion() = default;

    static Matrix<double> toMatrix(const Quaternion& quaternion){
        Matrix<double> result;
        result.resize(4, 4);

        // Quaternion to rotation matrix conversion
        result[0][0] = 1 - 2 * (quaternion.y() * quaternion.y() + quaternion.z() * quaternion.z());
        result[0][1] = 2 * (quaternion.x() * quaternion.y() - quaternion.z() * quaternion.w());
        result[0][2] = 2 * (quaternion.x() * quaternion.z() + quaternion.y() * quaternion.w());
        result[1][0] = 2 * (quaternion.x() * quaternion.y() + quaternion.z() * quaternion.w());
        result[1][1] = 1 - 2 * (quaternion.x() * quaternion.x() + quaternion.z() * quaternion.z());
        result[1][2] = 2 * (quaternion.y() * quaternion.z() - quaternion.x() * quaternion.w());
        result[2][0] = 2 * (quaternion.x() * quaternion.z() - quaternion.y() * quaternion.w());
        result[2][1] = 2 * (quaternion.y() * quaternion.z() + quaternion.x() * quaternion.w());
        result[2][2] = 1 - 2 * (quaternion.x() * quaternion.x() + quaternion.y() * quaternion.y());

        // Setting the last row and column for the 4x4 matrix (homogeneous coordinates)
        result[3][0] = 0;
        result[3][1] = 0;
        result[3][2] = 0;
        result[3][3] = 1;

        return result;
    }

    static double degreesToRadians(double degrees) {
        return degrees * M_PI / 180;
    }

    void normalize() {
        double mag = magnitude();
        _w /= mag;
        _x /= mag;
        _y /= mag;
        _z /= mag;

        approximate(_w, _x, _y, _z);
    }

    Quaternion conjugate() const {
        return Quaternion(_w, -_x, -_y, -_z);
    }

    double magnitude() const {
        return std::sqrt(_w * _w + _x * _x + _y * _y + _z * _z);
    }

    Quaternion operator*(const Quaternion &other) const noexcept {
        double newW = _w * other._w - _x * other._x - _y * other._y - _z * other._z;
        double newX = _w * other._x + _x * other._w + _y * other._z - _z * other._y;
        double newY = _w * other._y - _x * other._z + _y * other._w + _z * other._x;
        double newZ = _w * other._z + _x * other._y - _y * other._x + _z * other._w;

        approximate(newW, newX, newY, newZ);

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

        approximate(_w, _x, _y, _z);

        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const Quaternion &quaternion) {
        os << quaternion._w << ' ' << quaternion._x << ' ' << quaternion._y << ' ' << quaternion._z;
        return os;
    }

    friend std::istream &operator>>(std::istream &is, Quaternion &quaternion) = delete;

    double w() const { return _w; }

    double x() const { return _x; }

    double y() const { return _y; }

    double z() const { return _z; }

    void w(double w) { _w = w; }

    void x(double x) { _x = x; }

    void y(double y) { _y = y; }

    void z(double z) { _z = z; }
};

#endif
