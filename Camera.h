#ifndef CAMERA_H
#define CAMERA_H

#include "Matrix.h"
#include "Vector3.h"

class Camera {
private:
    Vector3 _position; // Camera position in world space
    Quaternion _rotation; // Camera orientation
    double _fov; // Field of view (in degrees)
    double _aspectRatio; // Aspect ratio (width / height)
    double _nearClip; // Near clipping plane
    double _farClip; // Far clipping plane

public:
    const Vector3 forward{0, 0, -1};

    const Vector3 up{0, 1, 0};

    const Vector3 right{1, 0, 0};

    Camera(const Vector3 &position = {0, 0, 0},
           const Quaternion &rotation = {0, 0, 0},
           double fov = 80,
           double aspectRatio = 16.0 / 9.0,
           double nearClip = 1,
           double farClip = 100)
        : _position{position},
          _rotation{rotation},
          _fov{fov},
          _aspectRatio{aspectRatio},
          _nearClip{nearClip},
          _farClip{farClip} {
    }

    //function accepts two parameters for its variants that work accordingly to camera rotation
    void translatePosition(const Vector3 &direction, double amount) {
        _position = _position + (direction * amount);
    }

    void moveForward(double amount) { translatePosition(_rotation * forward, amount); }

    void moveRight(double amount) { translatePosition(_rotation * right, amount); }

    void moveUp(double amount) { translatePosition(_rotation * up, amount); }

    void rotate(const Vector3 &eulerAngles) {
        Quaternion deltaRotation(eulerAngles); // Uses your Euler constructor
        _rotation = deltaRotation * _rotation; // Apply new rotation
    }

    Vector3 position() const { return _position; }

    Quaternion rotation() const { return _rotation; }

    double fov() const { return _fov; }

    double aspectRatio() const { return _aspectRatio; }

    double nearClip() const { return _nearClip; }

    double farClip() const { return _farClip; }

    void position(const Vector3 &position) { _position = position; }

    void rotation(const Quaternion &rotation) { _rotation = rotation; }

    void fov(double fov) { _fov = fov; }

    void aspectRatio(double aspectRatio) { _aspectRatio = aspectRatio; }

    void nearClip(double nearClip) { _nearClip = nearClip; }

    void farClip(double farClip) { _farClip = farClip; }
};

#include "Quaternion.h"

#endif //CAMERA_H
