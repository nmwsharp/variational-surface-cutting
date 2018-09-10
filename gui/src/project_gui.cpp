#include <project_gui.h>

ProjectGui::ProjectGui(Gui* _mainGui, std::string _name)
   : name(_name), mainGui(_mainGui), window(nullptr), isVisible(false)
{}

ProjectGui::~ProjectGui()
{
   safeDelete( window );
}

void ProjectGui::toggleVisible(void) {

   if( !isVisible )
   {
      if( show() )
      {
         isVisible = true;
      }
   }
   else
   {
      hide();
      isVisible = false;
   }

}

void ProjectGui::hide(void) {

   if(window != nullptr)
   {
      window->dispose();
      mainGui->replaceMainSceneMesh(nullptr); // use default mesh visualization
      window = nullptr;
   }

}

bool ProjectGui::getVisibility(void) const {
   return isVisible;
}

