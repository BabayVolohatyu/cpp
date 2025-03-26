#ifndef CAMERA_H
#define CAMERA_H

#include "Matrix.h"
#include "Quaternion.h"
#include "Vector3.h"
#include "Object.h"


class Camera {
private:
    Vector3 _position; // Camera position in world space
    Quaternion _rotation; // Camera orientation
    double _fov; // Field of view (in degrees)
    double _aspectRatio; // Aspect ratio (width / height)
    double _nearClip; // Near clipping plane
    double _farClip; // Far clipping plane

    static Matrix<double> translate(const Vector3 &position) {
        Matrix<double> result(4, 4);
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

    static void display(int screenX, int screenY) {
        std::cout << "\033[" << screenY << ";" << screenX << "H"; //i

        std::cout << "#" << std::flush;
        std::cout << "\033[0;0H";
    }

    void drawVertices(const Object &object, int screenWidth, int screenHeight) const {
        Matrix<double> modelMatrix = object.getModelMatrix();
        Matrix<double> viewMatrix = this->viewMatrix();
        Matrix<double> projectionMatrix = this->projectionMatrix();
        for (const Point &point: object.vertices()) {
            Matrix<double> worldPosVector(4, 1);
            worldPosVector[0][0] = point.x();
            worldPosVector[1][0] = point.y();
            worldPosVector[2][0] = point.z();
            worldPosVector[3][0] = 1;
            worldPosVector = modelMatrix * worldPosVector;
            Matrix<double> viewPosVector = viewMatrix * worldPosVector;
            Matrix<double> clipPosVector = projectionMatrix * viewPosVector;
            if (clipPosVector[3][0] != 0) {
                clipPosVector[0][0] /= clipPosVector[3][0];
                clipPosVector[1][0] /= clipPosVector[3][0];
                clipPosVector[2][0] /= clipPosVector[3][0];
            }
            // Convert to screen space (viewport transformation)
            double screenX = (clipPosVector[0][0] * 0.5 + 0.5) * screenWidth;
            double screenY = (clipPosVector[1][0] * 0.5 + 0.5) * screenHeight;//without 1.0 -... because inverse y-axis

            this->display(screenX, screenY);
        }
    }

public:
    const Vector3 worldForward{0, 0, -1};

    const Vector3 worldUp{0, 1, 0};

    const Vector3 worldRight{1, 0, 0};
    //add dependency between angles and those three vectors
    Vector3 forward{0, 0, -1};

    Vector3 up{0, 1, 0};

    Vector3 right{1, 0, 0};

    Camera(const Vector3 &position = {0, 0, 0},
           const Quaternion &rotation = {0, 0, 0},
           double fov = 90,
           double aspectRatio = 16 / 9,
           double nearClip = 0.1,
           double farClip = 1000)
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

    void lookAt(const Point &point) {
        Vector3 desiredForward = {
            point.x() - this->_position.x(),
            point.y() - this->_position.y(),
            point.z() - this->_position.z()
        };
        desiredForward.normalize();
        Vector3 currentForward = this->rotation() * this->forward;

        double cosTheta = currentForward.dot(desiredForward);

        if (cosTheta > 0.9999) return;

        if (cosTheta < -0.9999) {
            Quaternion turn = Quaternion::fromAxis(worldUp, M_PI);
            turn.normalize();
            Quaternion newRotation = turn * this->rotation();
            newRotation.normalize();
            this->_rotation = newRotation;
        }

        Vector3 rotAxis = currentForward.cross(desiredForward);
        rotAxis.normalize();

        double angle = acos(cosTheta);

        Quaternion rotQuat = Quaternion::fromAxis(rotAxis, angle);
        rotQuat.normalize();
        Quaternion newRotation = rotQuat * this->rotation();
        newRotation.normalize();
        this->_rotation = newRotation;
    }

    void rotate(const Vector3 &eulerAngles) {
        Quaternion deltaRotation(eulerAngles); // Uses your Euler constructor
        _rotation = deltaRotation * _rotation; // Apply new rotation
    }

    Matrix<double> viewMatrix() const {
        Matrix<double> translationMatrix = translate(-_position);
        Matrix<double> rotationMatrix = Quaternion::toMatrix(_rotation).transpose();
        return translationMatrix * rotationMatrix;
    }

    Matrix<double> projectionMatrix() const {
        double fovRadians = _fov * (M_PI / 180);
        double tanHalfFov = tan(fovRadians / 2);
        double range = _nearClip - _farClip;

        Matrix<double> proj(4, 4, 0);
        proj[0][0] = 1 / (tanHalfFov * _aspectRatio);
        proj[1][1] = 1 / tanHalfFov;
        proj[2][2] = (-_nearClip - _farClip) / range;
        proj[2][3] = 2 * _farClip * _nearClip / range;
        proj[3][2] = 1;
        proj[3][3] = 0;

        return proj;
    }

    //draws only vertices for now
    void drawObject(const Object &object, int screenWidth, int screenHeight) const {
        drawVertices(object, screenWidth, screenHeight);
    }

    Vector3 position() const { return _position; }

    Quaternion rotation() const { return _rotation; }

    double fov() const { return _fov; }

    double aspectRatio() const { return _aspectRatio; }

    double nearClip() const { return _nearClip; }

    double farClip() const { return _farClip; }

    void position(const Vector3 &position) { _position = position; }

    void rotation(const Quaternion &rotation) { _rotation = rotation; }

    void rotation(Quaternion &&rotation) { _rotation = std::move(rotation); }

    void fov(double fov) { _fov = fov; }

    void aspectRatio(double aspectRatio) { _aspectRatio = aspectRatio; }

    void nearClip(double nearClip) { _nearClip = nearClip; }

    void farClip(double farClip) { _farClip = farClip; }
};

#endif //CAMERA_H
