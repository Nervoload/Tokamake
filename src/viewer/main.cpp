#include <iostream>
#include <string>

#include "tokamak/viewer/viewer_app.hpp"

int main(int argc, char** argv) {
    tokamak::viewer::ViewerCliOptions options;
    std::string error;
    if (!tokamak::viewer::ParseViewerCliArgs(argc, argv, &options, &error)) {
        tokamak::viewer::PrintViewerUsage(argv[0]);
        if (!error.empty() && error != "help") {
            std::cerr << error << "\n";
            return 1;
        }
        return 0;
    }

    tokamak::viewer::ViewerApp app(options);
    return app.Run();
}
