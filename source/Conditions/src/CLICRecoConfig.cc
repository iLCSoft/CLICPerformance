#include "CLICRecoConfig.h"

#include <iomanip>


CLICRecoConfig aCLICRecoConfig;

CLICRecoConfig::CLICRecoConfig() : Processor("CLICRecoConfig") {

  // modify processor description
  _description = "CLICRecoConfig allows one to configure the processors to run via command line" ;

  for (auto& entry : m_allOptions) {
    std::string const& name = std::get<0>(entry);
    std::string& value      = std::get<1>(entry);
    Choices& choices        = std::get<2>(entry);
    std::stringstream optionDescription, choicesDescription;
    optionDescription << "Which option to use for " << name << ": ";
    toString(optionDescription, choices);
    optionDescription <<". Then use, e.g., Config." << name << value << " in the condition";
    registerProcessorParameter(name, optionDescription.str(), value, value);

    choicesDescription << "Possible values and conditions for option " << name;
    registerProcessorParameter(name+"Choices", choicesDescription.str(), choices, choices);
  }

}


void CLICRecoConfig::init() {

  // Print the initial parameters
  printParameters() ;

  for (auto const& entry : m_allOptions) {
    std::string const& name   = std::get<0>(entry);
    std::string const& choice = std::get<1>(entry);
    Choices const& choices    = std::get<2>(entry);
    checkOptions(choice, choices);
    for (auto& option: choices ) {
      m_options[name+option] = (choice == option);
    }
  }

  //print all options and values
  for (auto& option: m_options) {
    std::stringstream output;
    output << std::left << std::setw(35) << option.first << " "
           << std::boolalpha << option.second;
    streamlog_out(MESSAGE) << output.str() << std::endl;
  }

}

void CLICRecoConfig::processRunHeader( LCRunHeader* ){}

void CLICRecoConfig::processEvent( LCEvent* ) {
  modifyEvent( nullptr );
}


void CLICRecoConfig::modifyEvent( LCEvent* ) {
  streamlog_out(DEBUG9) << "Running Config" << std::endl;

  for (auto const& option : m_options) {
    setReturnValue(option.first, option.second);
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
