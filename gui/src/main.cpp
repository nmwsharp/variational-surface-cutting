#include<viewer.h>

#include <gui.h>



int main(int argc, char* argv[]) {

    std::cout << std::endl << std::endl << "=== Grand Central Gui ===" << std::endl << std::endl;

    Gui g;

    // If the program was started with an argument, autoload the mesh
    if(argc > 1) {
        g.loadMesh(argv[1]);
    }   

    // Start the main program loop; run until it exits
    g.viewer->startMainLoop();

    exit(EXIT_SUCCESS);

}

