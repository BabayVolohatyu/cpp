//
// Created by victooor on 04.03.25.
//

#ifndef GNOMETERMINALWINDOW_H
#define GNOMETERMINALWINDOW_H



#include <iostream>
#include <sstream>
#include <cstdlib>
#include <csignal>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

class GNOMETerminalWindow {
private:
    pid_t pid; // Process ID of the terminal

public:
    GNOMETerminalWindow() : pid(-1) {}

    bool open(int width, int height, double zoom) {
        pid = fork();
        if (pid == 0) {
            // Child process: Execute gnome-terminal with specified dimensions and zoom
            std::string command = "gnome-terminal --geometry=" + std::to_string(width) + "x" + std::to_string(height) +
                                  " --zoom " + std::to_string(zoom);
            execlp("sh", "sh", "-c", command.c_str(), nullptr);
            exit(EXIT_FAILURE); // If execlp fails
        } else if (pid > 0) {
            std::cout << "Opened gnome-terminal with PID: " << pid << std::endl;
            return true;
        } else {
            std::cerr << "Failed to fork process" << std::endl;
            return false;
        }
    }

    void close() {
        if (pid > 0) {
            kill(pid, SIGTERM); // Send termination signal
            waitpid(pid, nullptr, 0); // Wait for the process to exit
            std::cout << "Closed gnome-terminal with PID: " << pid << std::endl;
            pid = -1;
        } else {
            std::cerr << "No active terminal to close" << std::endl;
        }
    }
    void outputText(const std::string& text) {
        if (pid > 0) {
            std::stringstream ss;
            ss << "echo \"" << text << "\" > /proc/" << pid << "/fd/0";
            std::string command = ss.str();
            system(command.c_str());
        } else {
            std::cerr << "No active terminal to output text to" << std::endl;
        }
    }
    ~GNOMETerminalWindow() {
        close(); // Ensure cleanup on destruction
    }
};



#endif //GNOMETERMINALWINDOW_H
