//
// Created by austin on 2/26/16.
//

#include "InputHandler.h"
#include <iostream>
#include <algorithm>

InputHandler::InputHandler() {

}

double InputHandler::x() const {
    return _mouseState.x;
}

double InputHandler::y() const {
    return _mouseState.y;
}

double InputHandler::delX() const {
    return _mouseState.delX;
}

double InputHandler::delY() const {
    return _mouseState.delY;
}

double InputHandler::delWheel() const {
    return _mouseState.delWheel;
}

bool InputHandler::leftDown() const {
    return _mouseState.leftDown;
}

bool InputHandler::wheelDown() const {
    return _mouseState.wheelDown;
}

bool InputHandler::rightDown() const {
    return _mouseState.rightDown;
}

bool InputHandler::key(int key) const {
    return _keyboard.find(key) != _keyboard.end();
}

void InputHandler::x(double val, bool events) {
    std::swap(_mouseState.x, val);

    if (events) {
        _mouseState.delX = _mouseState.x - val;

        mouseMoved();

        emit(_mouseState);
        _mouseState.leftDragInit = false;
        _mouseState.wheelDragInit = false;
        _mouseState.rightDragInit = false;
    }
}

void InputHandler::y(double val, bool events) {
    std::swap(_mouseState.y, val);

    if (events) {
        _mouseState.delY = _mouseState.y - val;

        mouseMoved();

        emit(_mouseState);
        _mouseState.leftDragInit = false;
        _mouseState.wheelDragInit = false;
        _mouseState.rightDragInit = false;
    }
}

void InputHandler::mouseMoved() {
    if (_mouseState.leftDown) {
        if (!_mouseState.leftDragInit) {
            _mouseState.leftDragInit = true;
            _mouseState.startLeftX = _mouseState.x;
            _mouseState.startLeftY = _mouseState.y;
        }
        _mouseState.leftDragging = true;
    }
    if (_mouseState.wheelDown) {
        if (!_mouseState.wheelDragInit) {
            _mouseState.wheelDragInit = true;
            _mouseState.startWheelX = _mouseState.x;
            _mouseState.startWheelY = _mouseState.y;
        }
        _mouseState.wheelDragging = true;
    }
    if (_mouseState.rightDown) {
        if (!_mouseState.rightDragInit) {
            _mouseState.rightDragInit = true;
            _mouseState.startRightX = _mouseState.x;
            _mouseState.startRightY = _mouseState.y;
        }
        _mouseState.rightDragging = true;
    }
}

void InputHandler::pos(double x, double y, bool events) {
    std::swap(_mouseState.x, x);
    std::swap(_mouseState.y, y);

    if (events) {
        _mouseState.delX = _mouseState.x - x;
        _mouseState.delY = _mouseState.y - y;

        mouseMoved();

        emit(_mouseState);
        _mouseState.leftDragInit = false;
        _mouseState.wheelDragInit = false;
        _mouseState.rightDragInit = false;
    }
}

void InputHandler::delX(double val, bool events) {
    _mouseState.delX = val;
    if (events) {
        emit(_mouseState);
    }
}

void InputHandler::delY(double val, bool events) {
    _mouseState.delY = val;
    if (events) {
        emit(_mouseState);
    }
}


void InputHandler::delWheel(double val, bool events) {
    _mouseState.delWheel = val;
    if (events) {
        emit(_mouseState);
    }
    _mouseState.delWheel = 0;
}

void InputHandler::leftDown(bool val, bool events) {
    _mouseState.leftDown = val;
    if (events) {
        if (!val) {
            _mouseState.leftDragging = false;
            _mouseState.leftDragFinish = true;
        }
        emit(_mouseState);
        _mouseState.leftDragFinish = false;
    }
}

void InputHandler::wheelDown(bool val, bool events) {
    _mouseState.wheelDown = val;
    if (events) {
        if (!val) {
            _mouseState.wheelDragging = false;
            _mouseState.wheelDragFinish = true;
        }
        emit(_mouseState);
        _mouseState.wheelDragFinish = false;
    }
}

void InputHandler::rightDown(bool val, bool events) {
    _mouseState.rightDown = val;
    if (events) {
        if (!val) {
            _mouseState.rightDragging = false;
            _mouseState.rightDragFinish = true;
        }
        emit(_mouseState);
        _mouseState.rightDragFinish = false;
    }
}

void InputHandler::key(int key, bool down, bool events) {
    if (down) {
        _keyboard.insert(key);
    } else {
        _keyboard.erase(_keyboard.find((key)));
    }
}

void InputHandler::emit(MouseState &event) {
    for (int i = 0; i < mouseSubscribers.size(); i++) {
        mouseSubscribers.at(i)(event);
    }
}

void InputHandler::registerMouseListener(InputHandler::MouseListener listener) {
    mouseSubscribers.push_back(listener);
}

//void InputHandler::deregisterMouseListener(MouseListener listener) {
//    mouseSubscribers.erase(std::remove(mouseSubscribers.begin(), mouseSubscribers.end(), listener), mouseSubscribers.end());
//}

void InputHandler::windowResized(int w, int h) {
    for (int i = 0; i < windowSubscribers.size(); i++) {
        windowSubscribers.at(i)(w, h);
    }
}

void InputHandler::registerWindowListener(InputHandler::WindowListener listener) {
    windowSubscribers.push_back(listener);
}

//void InputHandler::deregisterWindowListener(InputHandler::WindowListener listener) {
//    windowSubscribers.erase(std::remove(windowSubscribers.begin(), windowSubscribers.end(), listener), windowSubscribers.end());
//}
