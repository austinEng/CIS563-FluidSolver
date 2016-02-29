//
// Created by austin on 2/26/16.
//

#ifndef FLUIDSOLVER_INPUTHANDLER_H
#define FLUIDSOLVER_INPUTHANDLER_H

#include <set>
#include <vector>
#include <functional>

class InputHandler {

public:
    static InputHandler& getInputHandler() {
        static InputHandler inputHandler;
        return inputHandler;
    }

    double x() const;
    double y() const;
    double delX() const;
    double delY() const;
    double delWheel() const;
    bool leftDown() const;
    bool wheelDown() const;
    bool rightDown() const;
    bool key(int key) const;

    void x(double val, bool events = true);
    void y(double val, bool events = true);
    void pos(double x, double y, bool events = true);
    void delX(double val, bool events = true);
    void delY(double val, bool events = true);
    void delWheel(double val, bool events = true);
    void leftDown(bool val, bool events = true);
    void wheelDown(bool val, bool events = true);
    void rightDown(bool val, bool events = true);
    void key(int key, bool down, bool events = true);

    struct MouseState {
        double x;
        double y;
        double delX;
        double delY;
        double startLeftX;
        double startLeftY;
        double startWheelX;
        double startWheelY;
        double startRightX;
        double startRightY;
        double delWheel;
        bool leftDown;
        bool wheelDown;
        bool rightDown;
        bool leftDragInit;
        bool wheelDragInit;
        bool rightDragInit;
        bool leftDragging;
        bool wheelDragging;
        bool rightDragging;
        bool leftDragFinish;
        bool wheelDragFinish;
        bool rightDragFinish;
    };

    //typedef void(*MouseListener)(MouseState&);
    //typedef void (*WindowListener)(int w, int h);
    typedef std::function<void(int, int)> WindowListener;
    typedef std::function<void(MouseState&)> MouseListener;
    void registerMouseListener(MouseListener listener);
    //void deregisterMouseListener(MouseListener listener);
    void registerWindowListener(WindowListener listener);
    //void deregisterWindowListener(WindowListener listener);
    void windowResized(int w, int h);

private:
    InputHandler();
    InputHandler(InputHandler const&) {} // prevent copies
    void operator=(InputHandler const&) {} // prevent assignments

    std::set<int> _keyboard;

    MouseState _mouseState;

    void mouseMoved();

    std::vector<MouseListener> mouseSubscribers;
    std::vector<WindowListener> windowSubscribers;
    void emit(MouseState &event);
};


#endif //FLUIDSOLVER_INPUTHANDLER_H
