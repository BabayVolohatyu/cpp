#ifndef OBJECT_H
#define OBJECT_H

#include <iostream>

#include "Vector3.h"
#include "Point.h"
#include "Quaternion.h"

class Object {
private:
    Vector3 _position;
    Quaternion _rotation;
    std::vector<Point> _vertices;

public:
    Object()
        : _position(0, 0, 0), _rotation(1, 0, 0, 0) {
    }

    Object(Vector3 &position, Quaternion &rotation)
        : _position(position), _rotation(rotation) {
    }

    void addVertex(Point &vertex) {
        _vertices.push_back(vertex);
    }

    void addVertex(Vector3 &vertex) {
        _vertices.push_back(vertex);
    }

    void addVertex(double x, double y, double z) {
        _vertices.push_back(Point(x, y, z));
    }

    void rotate(const Quaternion &rotationQuaternion) {
        _rotation = rotationQuaternion * _rotation;
        for (Point &p: _vertices) {
            Vector3 localVertex(p.x() - _position.x(), p.y() - _position.y(), p.z() - _position.z());
            Vector3 rotatedVertex = _rotation * localVertex;
            p.position(rotatedVertex + _position);
        }
    }

    std::vector<Point> vertices() const { return _vertices; }

    Vector3 position() const { return _position; }

    Quaternion rotation() const { return _rotation; }

    void position(const Vector3 &position) {
        Vector3 translation = position - _position;
        for (Point &p: _vertices) {
            p.x(p.x() + translation.x());
            p.y(p.y() + translation.y());
            p.z(p.z() + translation.z());
        }
        _position = position;
    }

    void position(double x, double y, double z) {
        Vector3 newPos(x, y, z);
        Vector3 translation = newPos - _position;
        for (Point &p: _vertices) {
            p.x(p.x() + translation.x());
            p.y(p.y() + translation.y());
            p.z(p.z() + translation.z());
        }
        _position = newPos;
    }

    void rotation(const Quaternion &rotation) {
        _rotation = Quaternion(1, 0, 0, 0);
        rotate(rotation);
    }

    void vertices(const std::vector<Point> &vertices) { _vertices = vertices; }
};
#endif
