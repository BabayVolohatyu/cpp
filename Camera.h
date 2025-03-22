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

    Matrix<double> translate(const Vector3 &position) const{
        Matrix<double> result(4,4);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i == j) result[i][j] = 1;
            }
        }
        result[0][3] = position.x();
        result[1][3] = position.y();
        result[2][3] = position.z();
        return result;
    }

public:
    const Vector3 forward{0, 0, -1};

    const Vector3 up{0, 1, 0};

    const Vector3 right{1, 0, 0};

    Camera(const Vector3 &position = {0, 0, 0},
           const Quaternion &rotation = {0, 0, 0},
           double fov = 80,
           double aspectRatio = 16 / 9,
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

    Matrix<double> viewMatrix() const{
      Matrix<double> translationMatrix = translate(-_position);
      Matrix<double> rotationMatrix = Quaternion::toMatrix(_rotation).transpose();
      return translationMatrix * rotationMatrix;
    }

    Matrix<double> projectionMatrix() const {
        double fovRadians = _fov * (M_PI / 180);
        double tanHalfFov = tan(fovRadians / 2);
        double range = _nearClip - _farClip;

        Matrix<double> proj(4,4,0);
        proj[0][0] = 1 / (tanHalfFov * _aspectRatio);
        proj[1][1] = 1 / tanHalfFov;
        proj[2][2] = (-_nearClip - _farClip) / range;
        proj[2][3] = 2 * _farClip * _nearClip / range;
        proj[3][2] = 1;
        proj[3][3] = 0;

        return proj;
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
