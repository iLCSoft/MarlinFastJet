#pragma once

#include <fastjet/contrib/AxesDefinition.hh>
#include <fastjet/JetDefinition.hh>

///------------------------------------------------------------------------
/// \class VLC_Axes
/// \brief Axes from exclusive VLC
///
/// Axes from VLC algorithm with E_scheme recombination.
///------------------------------------------------------------------------
class VLC_Axes : public fastjet::contrib::ExclusiveJetAxes {
public:
   /// Constructor
  VLC_Axes(fastjet::JetDefinition* jet_def)
    : fastjet::contrib::ExclusiveJetAxes(*jet_def)
  {
    //_plugin(plugin){
    
    setNPass(NO_REFINING);
   }

   /// Short description
   virtual std::string short_description() const {
      return "VLC";
   };
   
   /// Long description
   virtual std::string description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "VLC Axes";
      return stream.str();
   };
   
   /// For copying purposes
   virtual VLC_Axes* create() const {return new VLC_Axes(*this);}

  // ~VLC_Axes(){
  //   delete _plugin;
  // }

// private:
//   Plugin* _plugin;

};

