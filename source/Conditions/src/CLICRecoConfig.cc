#include "CLICRecoConfig.h"

#include <iomanip>


CLICRecoConfig aCLICRecoConfig;

CLICRecoConfig::CLICRecoConfig() : Processor("CLICRecoConfig") {

  // modify processor description
  _description = "CLICRecoConfig allows one to configure the processors to run via command line" ;

  std::stringstream trackingOptionDescription;
  trackingOptionDescription << "Which option to use for tracking: ";
  toString( trackingOptionDescription, m_trackingPossibleOptions );

  registerProcessorParameter("Tracking",
                             trackingOptionDescription.str(),
                             m_trackingOption,
                             m_trackingOption);


  std::stringstream overlayOptionDescription;
  overlayOptionDescription << "Which option to use for overlay: ";
  toString( overlayOptionDescription, m_overlayPossibleOptions );
  overlayOptionDescription <<". Then use, e.g., Config.Overlay3TeV in the condition";
  registerProcessorParameter("Overlay",
                             overlayOptionDescription.str(),
                             m_overlayChoice,
                             m_overlayChoice);

}


void CLICRecoConfig::init() {

  // Print the initial parameters
  printParameters() ;

  checkOptions( m_trackingOption, m_trackingPossibleOptions );

  if( m_trackingOption == "Truth" ) {

    m_truthTracking = true;
    m_conformalTracking = false;

  } else if( m_trackingOption == "Conformal" ) {

    m_truthTracking = false;
    m_conformalTracking = true;

  }



  // Treat option for overlay
  checkOptions( m_overlayChoice, m_overlayPossibleOptions );

  //fill map with overlay options
  for ( auto& overlayEnergy: m_overlayPossibleOptions ) {
    if( overlayEnergy == "False" ){ continue; }
    m_overlay["Overlay"+ overlayEnergy] = (overlayEnergy == m_overlayChoice);
  }

  for ( auto& overlayEnergy: m_overlay ) {
    streamlog_out(MESSAGE) << std::setw(15) << overlayEnergy.first << " "
                           << std::boolalpha << overlayEnergy.second  << std::endl;
  }

}

void CLICRecoConfig::processRunHeader( LCRunHeader* ){}

void CLICRecoConfig::processEvent( LCEvent* ) {
  modifyEvent( nullptr );
}


void CLICRecoConfig::modifyEvent( LCEvent* ) {
  streamlog_out(MESSAGE) << "Running Config" << std::endl;

  setReturnValue( "TruthTracking", m_truthTracking );
  setReturnValue( "ConformalTracking", m_conformalTracking );

  for (auto const& over : m_overlay) {
    setReturnValue( over.first, over.second );
  }

}

void CLICRecoConfig::check( LCEvent* ){}


void CLICRecoConfig::end(){}

void CLICRecoConfig::checkOptions( std::string const& optionValue, std::vector<std::string> const& options ) const {

  bool validOption = false;
  for(auto const& opt : options) {
    if( optionValue == opt) {
      validOption = true;
      break;
    }
  }
  if( not validOption ){
    std::stringstream errorMessage;
    errorMessage << "Invalid option value '" << optionValue << "'. Use one of [";
    toString( errorMessage, options );
    errorMessage << "]!";
    throw std::runtime_error( errorMessage.str() );
  }

}
